use crate::{
    core::selection::select_nodes_in_plane_direction,
    io::{
        obj::{read_obj, write_obj},
        stl::{read_stl, write_stl, StlFormat},
    },
};
use crate::error::Result;

/// Represents a standard 3D mesh.
/// A `Mesh` consists of nodes representing points in 3D space,
/// and cells which are triangles connecting these points.
///
/// # Example
/// ```
/// # use some_module::Mesh;
/// let nodes = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
/// let cells = vec![[0, 1, 2]];
/// let mesh = Mesh::new(nodes, cells);
/// ```
#[derive(Default, Debug, Clone)]
pub struct Mesh {
    pub nodes: Vec<[f64; 3]>,
    pub cells: Vec<[usize; 3]>,
}

impl Mesh {
    /// Constructs a new `Mesh`.
    ///
    /// # Arguments
    /// * `nodes`: Nodes of the mesh.
    /// * `cells`: Faces of the mesh.
    pub fn new(nodes: Vec<[f64; 3]>, cells: Vec<[usize; 3]>) -> Self {
        Mesh { nodes, cells }
    }

    /// Reads a mesh from an ASCII STL file using `read_stl`.
    /// See `read_stl`
    pub fn from_stl(filename: &str) -> Result<Self> {
        let (nodes, cells) = read_stl(filename)?;
        Ok(Mesh::new(nodes, cells))
    }

    /// Reads a mesh from an OBJ file using `read_obj`.
    /// See `read_obj`
    pub fn from_obj(filename: &str) -> Result<Self> {
        let (nodes, cells) = read_obj(filename)?;
        Ok(Mesh::new(nodes, cells))
    }

    /// Writes the mesh to an STL file.
    /// See `write_stl`
    pub fn write_stl(
        &self,
        filename: Option<&str>,
        format: Option<StlFormat>,
    ) -> Result<()> {
        write_stl(&self.nodes, &self.cells, filename, format)
    }

    /// Writes the mesh to an OBJ file.
    /// See `write_obj`
    pub fn write_obj(&self, filename: Option<&str>) -> Result<()> {
        write_obj(&self.nodes, &self.cells, filename)
    }

    /// Get nodes that make up cell triangles.
    /// See `get_mesh_triangles`
    pub fn triangles(&self) -> Vec<[[f64; 3]; 3]> {
        get_mesh_triangles(&self.nodes, &self.cells)
    }

    /// Get cell normals.
    /// See `get_mesh_cell_normals`
    pub fn cell_normals(&self) -> Vec<[f64; 3]> {
        get_mesh_cell_normals(&self.nodes, &self.cells)
    }

    /// Get nodes that are interpolated onto a slicing plane.
    /// See `find_node_intersections_with_plane`
    pub fn slice(&self, origin: [f64; 3], normal: [f64; 3]) -> Option<Vec<[f64; 3]>> {
        find_mesh_intersections_with_plane(&self.nodes, &self.cells, origin, normal)
    }

    /// Get mesh in direction of plane.
    /// See `clip_mesh_from_plane`
    pub fn clip(&self, origin: [f64; 3], normal: [f64; 3]) -> Mesh {
        let (new_nodes, new_cells) =
            clip_mesh_from_plane(&self.nodes, &self.cells, origin, normal);
        Mesh::new(new_nodes, new_cells)
    }
}

/// Returns the coordinates of the vertices for a single cell.
///
/// # Arguments
/// * `nodes`: Nodes of the mesh.
/// * `cell`: Cell id.
///
/// # Returns
/// A triangle in the form `[[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]]`.
fn get_triangle(nodes: &[[f64; 3]], cell: &[usize; 3]) -> [[f64; 3]; 3] {
    [nodes[cell[0]], nodes[cell[1]], nodes[cell[2]]]
}

/// Returns the coordinates of the vertices for each cell in the mesh.
///
/// # Arguments
/// * `nodes`: Nodes of the mesh.
/// * `cells`: Cells of the mesh.
///
/// # Returns
/// A `Vec<[[f64; 3]; 3]>` where each element represents a triangle.
pub fn get_mesh_triangles(
    nodes: &[[f64; 3]],
    cells: &[[usize; 3]],
) -> Vec<[[f64; 3]; 3]> {
    cells.iter().map(|cell| get_triangle(nodes, cell)).collect()
}

/// Computes and returns the normals for each cell in the mesh.
/// The normal for each cell is computed using the right-hand rule.
///
/// # Arguments
/// * `nodes`: Nodes of the mesh.
/// * `cells`: Cells of the mesh.
///
/// # Returns
/// A `Vec<[f64; 3]>` where each element is a normal corresponding to a cell.
pub fn get_mesh_cell_normals(nodes: &[[f64; 3]], cells: &[[usize; 3]]) -> Vec<[f64; 3]> {
    cells
        .iter()
        .map(|cell| {
            let a = nodes[cell[0]];
            let b = nodes[cell[1]];
            let c = nodes[cell[2]];

            // Calculate vectors ab and ac
            let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];

            // Calculate the cross product of ab and ac to get the normal
            let normal = [
                ab[1] * ac[2] - ab[2] * ac[1],
                ab[2] * ac[0] - ab[0] * ac[2],
                ab[0] * ac[1] - ab[1] * ac[0],
            ];

            // Normalize the normal vector
            let length =
                (normal[0].powi(2) + normal[1].powi(2) + normal[2].powi(2)).sqrt();
            [normal[0] / length, normal[1] / length, normal[2] / length]
        })
        .collect()
}

/// Finds the intersection points of a mesh with a slicing plane.
///
/// # Arguments
/// * `nodes`: Nodes of the mesh.
/// * `cells`: Cells of the mesh.
/// * `origin`: A 3D point `[x, y, z]` representing a point through which the
///             slice plane passes.
/// * `normal`: A 3D vector `[x, y, z]` representing the normal to the slice plane.
///
/// # Returns
/// Returns a `Vec<[f64; 3]>` containing the intersection points of the mesh with
/// the slice plane.
/// Each intersection point is represented as a 3D point `[x, y, z]`.
///
/// # Example
/// ```
/// # Example usage
/// let nodes = vec![
///     [0.0, 0.0, 0.0],
///     [1.0, 0.0, 0.0],
///     [0.0, 1.0, 0.0],
///     [0.0, 0.0, 1.0],
/// ];
/// let cells = vec![[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]];
/// let mesh = Mesh { nodes, cells };
///
/// let origin = [0.5, 0.5, 0.5];
/// let normal = [0.0, 0.0, 1.0];
///
/// let result = find_intersections_with_slice_plane(mesh, origin, normal);
/// ```
pub fn find_mesh_intersections_with_plane(
    nodes: &[[f64; 3]],
    cells: &[[usize; 3]],
    origin: [f64; 3],
    normal: [f64; 3],
) -> Option<Vec<[f64; 3]>> {
    let d = origin
        .iter()
        .zip(normal.iter())
        .map(|(o, n)| o * n)
        .sum::<f64>();
    let triangles = get_mesh_triangles(nodes, cells);

    // This is an arbitrary size for initial preallocation; you might have a better estimate
    let mut intersections: Vec<[f64; 3]> = Vec::with_capacity(2 * triangles.len());

    for [a, b, c] in triangles {
        let mut edge_intersections: Vec<[f64; 3]> = Vec::with_capacity(2);

        for &[point1, point2] in &[[a, b], [b, c], [c, a]] {
            let [x1, y1, z1] = point1;
            let [x2, y2, z2] = point2;

            let t = (d - normal
                .iter()
                .zip(point1.iter())
                .map(|(n, p)| n * p)
                .sum::<f64>())
                / normal
                    .iter()
                    .zip(point2.iter())
                    .zip(point1.iter())
                    .map(|((n, p2), p1)| n * (p2 - p1))
                    .sum::<f64>();

            if (0.0..=1.0).contains(&t) {
                let x = x1 + (x2 - x1) * t;
                let y = y1 + (y2 - y1) * t;
                let z = z1 + (z2 - z1) * t;
                edge_intersections.push([x, y, z]);
            }
        }

        if edge_intersections.len() == 2 {
            intersections.push(edge_intersections[0]);
            intersections.push(edge_intersections[1]);
        }
    }

    if intersections.is_empty() {
        None
    } else {
        Some(intersections)
    }
}

/// Clips the mesh with a given slicing plane, keeping only the elements and nodes that
/// lie in the direction of the plane's normal.
///
/// # Arguments
/// * `origin`: A 3D point `[x, y, z]` representing a point through which the
///             slice plane passes.
/// * `normal`: A 3D vector `[x, y, z]` representing the normal to the slice plane.
///
/// # Returns
/// Returns a new `Mesh` object containing only the nodes and cells (triangles)
/// that are on the "positive" side of the slicing plane, i.e., in the direction
/// of the plane's normal.
pub fn clip_mesh_from_plane(
    nodes: &[[f64; 3]],
    cells: &[[usize; 3]],
    origin: [f64; 3],
    normal: [f64; 3],
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    // Get the indices of nodes that lie in the direction of the plane's normal
    let selected_node_indices = select_nodes_in_plane_direction(nodes, origin, normal);

    // Create a new list of nodes using the selected indices
    let new_nodes = selected_node_indices
        .iter()
        .map(|&i| nodes[i])
        .collect::<Vec<[f64; 3]>>();

    // Create a new list of cells (elements) using only the selected nodes and reindexing them
    let new_cells = cells
        .iter()
        .filter_map(|cell| {
            if cell.iter().all(|&i| selected_node_indices.contains(&i)) {
                Some([
                    selected_node_indices
                        .iter()
                        .position(|&x| x == cell[0])
                        .unwrap(),
                    selected_node_indices
                        .iter()
                        .position(|&x| x == cell[1])
                        .unwrap(),
                    selected_node_indices
                        .iter()
                        .position(|&x| x == cell[2])
                        .unwrap(),
                ])
            } else {
                None
            }
        })
        .collect::<Vec<[usize; 3]>>();

    // Create a new Mesh using the selected nodes and cells
    (new_nodes, new_cells)
}

#[cfg(test)]
mod tests {
    use super::*;
    const EPSILON: f64 = 1e-9;

    #[test]
    fn test_get_triangles() {
        let nodes = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let cells = vec![[0, 1, 2], [0, 2, 3]];
        let triangles = get_mesh_triangles(&nodes, &cells);

        let expected_triangles = vec![
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        ];

        assert_eq!(triangles, expected_triangles)
    }

    #[test]
    fn test_get_face_normals() {
        let nodes = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let cells = vec![[0, 1, 2], [0, 2, 3]];
        let normals = get_mesh_cell_normals(&nodes, &cells);

        let expected_normals = vec![[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]];

        assert_eq!(normals, expected_normals)
    }

    #[test]
    fn test_slice_mesh_with_plane() {
        // A simple tetrahedron for slicing
        let nodes = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, 0.5, 1.0],
        ];
        let cells = vec![[0, 1, 2], [0, 1, 3], [1, 2, 3], [0, 2, 3]];
        let mesh = Mesh::new(nodes, cells);

        // A plane that cuts horizontally through the tetrahedron
        let origin = [0.0, 0.0, 0.5];
        let normal = [0.0, 0.0, 1.0];

        let intersections = mesh
            .slice(origin, normal)
            .expect("Should find intersections");

        // The slice of a tetrahedron is a triangle, which is 3 line segments (6 points).
        assert_eq!(
            intersections.len(),
            6,
            "Slice should produce 3 line segments"
        );
        // All intersection points should have a z-coordinate of 0.5
        for point in intersections {
            assert!((point[2] - 0.5).abs() < EPSILON);
        }
    }

    #[test]
    fn test_slice_mesh_no_intersection() {
        let nodes = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let cells = vec![[0, 1, 2]];
        let mesh = Mesh::new(nodes, cells);

        // A plane far away from the triangle
        let origin = [0.0, 0.0, 10.0];
        let normal = [0.0, 0.0, 1.0];

        let result = mesh.slice(origin, normal);
        assert!(result.is_none(), "Should not find any intersections");
    }

    #[test]
    fn test_clip_mesh_from_plane() {
        // A 2x1x1 block centered at the origin
        let nodes = vec![
            [-1., -0.5, -0.5],
            [-1., -0.5, 0.5],
            [-1., 0.5, -0.5],
            [-1., 0.5, 0.5],
            [1., -0.5, -0.5],
            [1., -0.5, 0.5],
            [1., 0.5, -0.5],
            [1., 0.5, 0.5],
        ];
        let cells = vec![
            [0, 1, 3],
            [0, 3, 2],
            [4, 7, 5],
            [4, 6, 7],
            [0, 4, 5],
            [0, 5, 1],
            [2, 3, 7],
            [2, 7, 6],
            [0, 6, 4],
            [0, 2, 6],
            [1, 5, 7],
            [1, 7, 3],
        ];
        let mesh = Mesh::new(nodes, cells);

        // A plane that cuts the block in half at x=0
        let origin = [0.0, 0.0, 0.0];
        let normal = [1.0, 0.0, 0.0]; // Keep the +X side

        let new_mesh = mesh.clip(origin, normal);

        // Should keep the 4 nodes on the +X side
        assert_eq!(new_mesh.nodes.len(), 4);
        // Those 4 nodes should form a single quad face (2 triangles)
        assert_eq!(new_mesh.cells.len(), 2);

        // Verify all kept nodes have a positive or zero x-coordinate
        for node in new_mesh.nodes {
            assert!(node[0] >= -EPSILON);
        }
    }
}
