//! NumPy `.npy` files: written as little-endian float32, read as float32 or float64.
//! Arrays are C-ordered (last index fastest).

use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

use anyhow::{Context, Result, bail, ensure};

pub fn write(path: &Path, shape: &[usize], data: &[f32]) -> Result<()> {
    let header = header(shape, data.len())?;
    let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    let mut w = BufWriter::new(file);
    w.write_all(&header)?;
    for chunk in data.chunks(1 << 14) {
        let bytes: Vec<u8> = chunk.iter().flat_map(|v| v.to_le_bytes()).collect();
        w.write_all(&bytes)?;
    }
    w.flush()?;
    Ok(())
}

/// The bytes of an `.npy` file, for sending rather than writing.
pub fn to_bytes(shape: &[usize], data: &[f32]) -> Result<Vec<u8>> {
    let mut bytes = header(shape, data.len())?;
    bytes.reserve(4 * data.len());
    bytes.extend(data.iter().flat_map(|v| v.to_le_bytes()));
    Ok(bytes)
}

/// Magic, version, header length and the header dictionary, padded so the data starts on a
/// 64-byte boundary.
fn header(shape: &[usize], len: usize) -> Result<Vec<u8>> {
    ensure!(
        shape.iter().product::<usize>() == len,
        "shape {shape:?} doesn't fit {len} values"
    );
    let dims: Vec<String> = shape.iter().map(usize::to_string).collect();
    let shape = if dims.len() == 1 {
        format!("({},)", dims[0])
    } else {
        format!("({})", dims.join(", "))
    };
    let mut dict = format!("{{'descr': '<f4', 'fortran_order': False, 'shape': {shape}, }}");
    let used = 10 + dict.len() + 1;
    dict.push_str(&" ".repeat((64 - used % 64) % 64));
    dict.push('\n');
    let mut bytes = b"\x93NUMPY\x01\x00".to_vec();
    bytes.extend((dict.len() as u16).to_le_bytes());
    bytes.extend(dict.as_bytes());
    Ok(bytes)
}

/// Reads an array of float32 or float64 values, as float32.
pub fn read(path: &Path) -> Result<(Vec<usize>, Vec<f32>)> {
    let bytes = fs::read(path)?;
    let h = parse_header(&bytes)?;
    let count: usize = h.shape.iter().product();
    let data = &bytes[h.data_start..];
    let width = if h.wide { 8 } else { 4 };
    ensure!(
        data.len() == width * count,
        "expected {count} float{} values",
        8 * width
    );
    Ok((h.shape, decode(data, h.wide)))
}

/// An array read a piece at a time, as float32: for scanning files too large to hold.
pub struct Stream {
    pub shape: Vec<usize>,
    reader: BufReader<File>,
    wide: bool,
    left: usize,
}

/// Opens an array for reading in pieces.
pub fn stream(path: &Path) -> Result<Stream> {
    let mut reader = BufReader::new(File::open(path)?);
    let mut prefix = vec![0u8; 10];
    reader.read_exact(&mut prefix).context("not an .npy file")?;
    ensure!(prefix.starts_with(b"\x93NUMPY"), "not an .npy file");
    let len = match prefix[6] {
        1 => u16::from_le_bytes([prefix[8], prefix[9]]) as usize,
        2 | 3 => {
            let mut more = [0u8; 2];
            reader.read_exact(&mut more).context("truncated header")?;
            prefix.extend(more);
            u32::from_le_bytes(prefix[8..12].try_into()?) as usize
        }
        v => bail!("unsupported .npy version {v}"),
    };
    let start = prefix.len();
    prefix.resize(start + len, 0);
    reader
        .read_exact(&mut prefix[start..])
        .context("truncated header")?;
    let h = parse_header(&prefix)?;
    let left = h.shape.iter().product();
    Ok(Stream {
        shape: h.shape,
        reader,
        wide: h.wide,
        left,
    })
}

impl Stream {
    /// The next `max` values or fewer, into `out`; false once there are none left.
    pub fn next_chunk(&mut self, out: &mut Vec<f32>, max: usize) -> Result<bool> {
        out.clear();
        let n = max.min(self.left);
        if n == 0 {
            return Ok(false);
        }
        let mut bytes = vec![0u8; n * if self.wide { 8 } else { 4 }];
        self.reader
            .read_exact(&mut bytes)
            .context("the file ends before its data does")?;
        out.extend(decode(&bytes, self.wide));
        self.left -= n;
        Ok(true)
    }
}

/// An array's shape and element type, and where its data starts.
struct Header {
    shape: Vec<usize>,
    /// float64 rather than float32.
    wide: bool,
    data_start: usize,
}

/// Parses the header at the start of `bytes`, which must hold all of it.
fn parse_header(bytes: &[u8]) -> Result<Header> {
    ensure!(
        bytes.len() >= 10 && bytes.starts_with(b"\x93NUMPY"),
        "not an .npy file"
    );
    let (len, start) = match bytes[6] {
        1 => (u16::from_le_bytes([bytes[8], bytes[9]]) as usize, 10),
        2 | 3 => {
            ensure!(bytes.len() >= 12, "truncated header");
            (u32::from_le_bytes(bytes[8..12].try_into()?) as usize, 12)
        }
        v => bail!("unsupported .npy version {v}"),
    };
    ensure!(bytes.len() >= start + len, "truncated header");
    let header = std::str::from_utf8(&bytes[start..start + len])?;
    let wide = match quoted(after(header, "descr")?)? {
        "<f4" => false,
        "<f8" => true,
        other => bail!("unsupported dtype {other}; use little-endian float32 or float64"),
    };
    ensure!(
        after(header, "fortran_order")?.starts_with("False"),
        "Fortran-ordered arrays aren't supported; save with numpy.ascontiguousarray"
    );
    let shape_text = after(header, "shape")?;
    let close = shape_text.find(')').context("malformed shape")?;
    let shape: Vec<usize> = shape_text[1..close]
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(str::parse)
        .collect::<Result<_, _>>()
        .context("malformed shape")?;
    Ok(Header {
        shape,
        wide,
        data_start: start + len,
    })
}

/// Little-endian float32 or float64 values, as float32.
fn decode(data: &[u8], wide: bool) -> Vec<f32> {
    if wide {
        data.as_chunks::<8>()
            .0
            .iter()
            .map(|&b| f64::from_le_bytes(b) as f32)
            .collect()
    } else {
        data.as_chunks::<4>()
            .0
            .iter()
            .map(|&b| f32::from_le_bytes(b))
            .collect()
    }
}

/// The text after `'key':` in a header dictionary.
fn after<'a>(header: &'a str, key: &str) -> Result<&'a str> {
    let tag = format!("'{key}':");
    let at = header
        .find(&tag)
        .with_context(|| format!("no {key} in the header"))?;
    Ok(header[at + tag.len()..].trim_start())
}

/// The contents of a leading single-quoted string.
fn quoted(text: &str) -> Result<&str> {
    let rest = text
        .strip_prefix('\'')
        .context("expected a quoted string")?;
    let end = rest.find('\'').context("unterminated string")?;
    Ok(&rest[..end])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips() {
        let dir = std::env::temp_dir().join(format!("vibe-rations-npy-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        for shape in [vec![7], vec![2, 3], vec![2, 3, 4]] {
            let data: Vec<f32> = (0..shape.iter().product::<usize>())
                .map(|v| v as f32 * 0.5)
                .collect();
            let path = dir.join("a.npy");
            write(&path, &shape, &data).unwrap();
            let bytes = fs::read(&path).unwrap();
            let header_len = u16::from_le_bytes([bytes[8], bytes[9]]) as usize;
            assert_eq!((10 + header_len) % 64, 0);
            assert_eq!(to_bytes(&shape, &data).unwrap(), bytes);
            assert_eq!(read(&path).unwrap(), (shape.clone(), data.clone()));
            // In pieces, the same values.
            let mut s = stream(&path).unwrap();
            assert_eq!(s.shape, shape);
            let (mut all, mut chunk) = (Vec::new(), Vec::new());
            while s.next_chunk(&mut chunk, 5).unwrap() {
                all.extend_from_slice(&chunk);
            }
            assert_eq!(all, data);
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn reads_float64() {
        let dir = std::env::temp_dir().join(format!("vibe-rations-npy64-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let header = "{'descr': '<f8', 'fortran_order': False, 'shape': (2,), }";
        let mut bytes = b"\x93NUMPY\x01\x00".to_vec();
        bytes.extend((header.len() as u16).to_le_bytes());
        bytes.extend(header.as_bytes());
        bytes.extend(1.5f64.to_le_bytes());
        bytes.extend((-2.0f64).to_le_bytes());
        let path = dir.join("b.npy");
        fs::write(&path, bytes).unwrap();
        assert_eq!(read(&path).unwrap(), (vec![2], vec![1.5, -2.0]));
        let mut chunk = Vec::new();
        assert!(stream(&path).unwrap().next_chunk(&mut chunk, 10).unwrap());
        assert_eq!(chunk, vec![1.5, -2.0]);
        fs::remove_dir_all(dir).unwrap();
    }
}
