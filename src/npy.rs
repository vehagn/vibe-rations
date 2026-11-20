//! NumPy `.npy` files: written as little-endian float32, read as float32 or float64.
//! Arrays are C-ordered (last index fastest).

use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result, bail, ensure};

pub fn write(path: &Path, shape: &[usize], data: &[f32]) -> Result<()> {
    ensure!(
        shape.iter().product::<usize>() == data.len(),
        "shape {shape:?} doesn't fit {} values",
        data.len()
    );
    let dims: Vec<String> = shape.iter().map(usize::to_string).collect();
    let shape = if dims.len() == 1 {
        format!("({},)", dims[0])
    } else {
        format!("({})", dims.join(", "))
    };
    let mut header = format!("{{'descr': '<f4', 'fortran_order': False, 'shape': {shape}, }}");
    // Magic, version and length take 10 bytes; pad so the data starts on a 64-byte boundary.
    let used = 10 + header.len() + 1;
    header.push_str(&" ".repeat((64 - used % 64) % 64));
    header.push('\n');

    let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    let mut w = BufWriter::new(file);
    w.write_all(b"\x93NUMPY\x01\x00")?;
    w.write_all(&(header.len() as u16).to_le_bytes())?;
    w.write_all(header.as_bytes())?;
    for chunk in data.chunks(1 << 14) {
        let bytes: Vec<u8> = chunk.iter().flat_map(|v| v.to_le_bytes()).collect();
        w.write_all(&bytes)?;
    }
    w.flush()?;
    Ok(())
}

/// Reads a C-ordered float32 or float64 array: its shape and values.
pub fn read(path: &Path) -> Result<(Vec<usize>, Vec<f32>)> {
    let bytes = fs::read(path)?;
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
    let descr = quoted(after(header, "descr")?)?;
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
    let count: usize = shape.iter().product();
    let data = &bytes[start + len..];
    let values: Vec<f32> = match descr {
        "<f4" => {
            ensure!(data.len() == 4 * count, "expected {count} float32 values");
            data.as_chunks::<4>()
                .0
                .iter()
                .map(|&b| f32::from_le_bytes(b))
                .collect()
        }
        "<f8" => {
            ensure!(data.len() == 8 * count, "expected {count} float64 values");
            data.as_chunks::<8>()
                .0
                .iter()
                .map(|&b| f64::from_le_bytes(b) as f32)
                .collect()
        }
        other => bail!("unsupported dtype {other}; use little-endian float32 or float64"),
    };
    Ok((shape, values))
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
            assert_eq!(read(&path).unwrap(), (shape, data));
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
        fs::remove_dir_all(dir).unwrap();
    }
}
