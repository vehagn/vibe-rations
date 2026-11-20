//! The computational grid: the interior nodes the user sees, padded by a halo.

/// Nodes beyond the interior on every side. The fourth-order stencil reaches two nodes past
/// the one it updates; the halo stays zero, or holds the free surface's mirror images.
pub const HALO: usize = 2;

/// A regular grid. Node (i, j, k) sits at (i·dx, j·dy, k·dz). z points down; with a free
/// surface, z = 0 is the surface.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Grid {
    /// Interior nodes along x, y and z.
    pub n: [usize; 3],
    /// Spacing along x, y and z, in metres.
    pub h: [f64; 3],
}

impl Grid {
    pub fn new(n: [usize; 3], h: [f64; 3]) -> Self {
        Self { n, h }
    }

    /// Nodes along each axis of a padded array.
    pub fn padded(&self) -> [usize; 3] {
        self.n.map(|n| n + 2 * HALO)
    }

    /// Distance in memory between neighbouring nodes along each axis of a padded array.
    pub fn strides(&self) -> [usize; 3] {
        let p = self.padded();
        [1, p[0], p[0] * p[1]]
    }

    /// Values in one z-plane of a padded array.
    pub fn plane_len(&self) -> usize {
        let p = self.padded();
        p[0] * p[1]
    }

    /// Values in a padded array.
    pub fn padded_len(&self) -> usize {
        self.padded().iter().product()
    }

    /// Values in an interior-only array.
    pub fn interior_len(&self) -> usize {
        self.n.iter().product()
    }

    /// Index of interior node (i, j, k) in a padded array, x fastest.
    #[inline]
    pub fn idx(&self, i: usize, j: usize, k: usize) -> usize {
        let p = self.padded();
        ((k + HALO) * p[1] + j + HALO) * p[0] + i + HALO
    }

    /// Index of node (i, j, k) in an interior-only array, x fastest.
    #[inline]
    pub fn interior_idx(&self, i: usize, j: usize, k: usize) -> usize {
        (k * self.n[1] + j) * self.n[0] + i
    }

    pub fn cell_volume(&self) -> f64 {
        self.h.iter().product()
    }

    /// Distance from the first to the last node along each axis.
    pub fn extent(&self) -> [f64; 3] {
        std::array::from_fn(|a| (self.n[a] - 1) as f64 * self.h[a])
    }

    pub fn position(&self, node: [usize; 3]) -> [f64; 3] {
        std::array::from_fn(|a| node[a] as f64 * self.h[a])
    }

    /// The node nearest a position, or None if the position is off the grid.
    pub fn nearest_node(&self, position: [f64; 3]) -> Option<[usize; 3]> {
        let mut node = [0; 3];
        for a in 0..3 {
            let r = (position[a] / self.h[a]).round();
            if !(0.0..=(self.n[a] - 1) as f64).contains(&r) {
                return None;
            }
            node[a] = r as usize;
        }
        Some(node)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn indexes_skip_the_halo() {
        let g = Grid::new([3, 4, 5], [1.0; 3]);
        assert_eq!(g.padded(), [7, 8, 9]);
        assert_eq!(g.idx(0, 0, 0), (2 * 8 + 2) * 7 + 2);
        assert_eq!(g.idx(1, 0, 0) - g.idx(0, 0, 0), g.strides()[0]);
        assert_eq!(g.idx(0, 1, 0) - g.idx(0, 0, 0), g.strides()[1]);
        assert_eq!(g.idx(0, 0, 1) - g.idx(0, 0, 0), g.strides()[2]);
        assert_eq!(g.interior_idx(2, 3, 4), g.interior_len() - 1);
    }

    #[test]
    fn snaps_positions_to_nodes() {
        let g = Grid::new([11, 11, 6], [10.0, 10.0, 5.0]);
        assert_eq!(g.nearest_node([14.9, 15.1, 12.4]), Some([1, 2, 2]));
        assert_eq!(g.nearest_node([100.0, 0.0, 25.0]), Some([10, 0, 5]));
        assert_eq!(g.nearest_node([106.0, 0.0, 0.0]), None);
        assert_eq!(g.nearest_node([0.0, -6.0, 0.0]), None);
    }
}
