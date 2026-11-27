//! The computational grid: the interior nodes the user sees, padded by a halo.
//!
//! A grid one node across y is a 2D grid, a plane of x and z: the fields are taken to be the
//! same at every y (plane strain), so nothing is differentiated along y, which needs no halo
//! there, and what the solver computes on it is exactly what its 3D scheme computes for fields
//! that don't vary along y. Its spacing along y is 1 m, so a source on it is a line source of
//! that strength per metre and a cell's volume is its area times 1 m.

/// A regular grid. Node (i, j, k) sits at (i·dx, j·dy, k·dz). z points down; with a free
/// surface, z = 0 is the surface.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Grid {
    /// Interior nodes along x, y and z; one along y for a 2D grid.
    pub n: [usize; 3],
    /// Spacing along x, y and z, in metres.
    pub h: [f64; 3],
    /// Nodes beyond the interior on every side of a padded array, but none across a 2D grid
    /// ([`Grid::halos`]). A stencil reaches its half-length past the node it updates; the halo
    /// stays zero, or holds the free surface's mirror images.
    pub halo: usize,
}

impl Grid {
    /// A grid with a two-node halo, enough for the fourth-order stencil. One node across y
    /// makes it a 2D grid, whose spacing along y is 1 m whatever `h` says.
    pub fn new(n: [usize; 3], h: [f64; 3]) -> Self {
        let h = if n[1] == 1 { [h[0], 1.0, h[2]] } else { h };
        Self { n, h, halo: 2 }
    }

    /// The same grid with a halo `halo` nodes wide (at least two).
    pub fn with_halo(self, halo: usize) -> Self {
        Self {
            halo: halo.max(2),
            ..self
        }
    }

    /// Whether it is a 2D grid: one node across y, the same at every y.
    #[inline(always)]
    pub fn is_2d(&self) -> bool {
        self.n[1] == 1
    }

    /// Whether anything varies along `axis`: every axis but y on a 2D grid.
    #[inline(always)]
    pub fn spans(&self, axis: usize) -> bool {
        axis != 1 || !self.is_2d()
    }

    /// How far a stencil `l` nodes long reaches along each axis: nowhere across a 2D grid.
    pub fn reach(&self, l: usize) -> [usize; 3] {
        std::array::from_fn(|a| if self.spans(a) { l } else { 0 })
    }

    /// The halo's width along each axis: [`Grid::halo`], and none across a 2D grid.
    pub fn halos(&self) -> [usize; 3] {
        self.reach(self.halo)
    }

    /// Nodes along each axis of a padded array.
    pub fn padded(&self) -> [usize; 3] {
        let halos = self.halos();
        std::array::from_fn(|a| self.n[a] + 2 * halos[a])
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

    /// The shape of an interior-only array with its slowest axis first, (nz, ny, nx), as
    /// `.npy` files hold it.
    pub fn shape(&self) -> [usize; 3] {
        let [nx, ny, nz] = self.n;
        [nz, ny, nx]
    }

    /// Index of interior node (i, j, k) in a padded array, x fastest.
    #[inline]
    pub fn idx(&self, i: usize, j: usize, k: usize) -> usize {
        let (p, [hx, hy, hz]) = (self.padded(), self.halos());
        ((k + hz) * p[1] + j + hy) * p[0] + i + hx
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

    /// A position's grid coordinates, position / spacing along each axis (a node's are whole
    /// numbers), or None if it is off the grid. Coordinates within 1e-9 of a whole number are
    /// that number, so a position on a node is read exactly as the node.
    pub fn coordinates(&self, position: [f64; 3]) -> Option<[f64; 3]> {
        let mut point = [0.0; 3];
        for a in 0..3 {
            let mut x = position[a] / self.h[a];
            if (x - x.round()).abs() < 1e-9 {
                x = x.round();
            }
            if !(0.0..=(self.n[a] - 1) as f64).contains(&x) {
                return None;
            }
            point[a] = x;
        }
        Some(point)
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

/// Something that stands for a point in grid coordinates: a node, or any point
/// ([`Grid::coordinates`]).
pub trait AsPoint {
    fn point(&self) -> [f64; 3];
}

impl AsPoint for [usize; 3] {
    fn point(&self) -> [f64; 3] {
        self.map(|n| n as f64)
    }
}

impl AsPoint for [f64; 3] {
    fn point(&self) -> [f64; 3] {
        *self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_wider_halo_pads_more() {
        let g = Grid::new([3, 4, 5], [1.0; 3]).with_halo(8);
        assert_eq!(g.padded(), [19, 20, 21]);
        assert_eq!(g.idx(0, 0, 0), (8 * 20 + 8) * 19 + 8);
        assert_eq!(Grid::new([3, 4, 5], [1.0; 3]).with_halo(1).halo, 2);
    }

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
    fn a_2d_grid_has_no_halo_across_y() {
        let g = Grid::new([5, 1, 6], [10.0, 10.0, 5.0]).with_halo(4);
        assert!(g.is_2d() && !Grid::new([5, 2, 6], [1.0; 3]).is_2d());
        assert_eq!(g.h, [10.0, 1.0, 5.0]);
        assert_eq!(g.cell_volume(), 50.0);
        assert_eq!(g.halos(), [4, 0, 4]);
        assert_eq!(g.padded(), [13, 1, 14]);
        assert_eq!(g.idx(0, 0, 0), 4 * 13 + 4);
        assert_eq!(g.idx(0, 0, 1) - g.idx(0, 0, 0), g.strides()[2]);
        assert_eq!(g.strides()[2], g.plane_len());
        assert_eq!(g.interior_idx(4, 0, 5), g.interior_len() - 1);
        assert_eq!(g.coordinates([12.5, 0.0, 25.0]), Some([1.25, 0.0, 5.0]));
        assert_eq!(g.coordinates([12.5, 0.5, 25.0]), None);
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
