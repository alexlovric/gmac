use crate::core::{
    clusters::generate_block_cluster,
    transformation::{rotate_nodes, rotate_node, build_rotation_matrix},
    mesh::Mesh,
    selection::select_nodes_in_box,
};

/// `DesignBlock` represents a geometric block in 3D space.
///
/// # Arguments
/// * `nodes`: A vector of nodes (or points) inside the block, each represented by a `[f64; 3]` array.
/// * `length`: The lengths of the block in the x, y, and z directions as a `[f64; 3]` array.
/// * `centre`: The coordinates of the center of the block as a `[f64; 3]` array.
/// * `theta`: The angles (in radians) to rotate the block around the x, y, and z axes as a `[f64; 3]` array.
/// * `resolution`: The resolution of the block in each dimension, represented as a `[usize; 3]` array.
/// * `corner_nodes`: The 8 corner in x, y, and z directions as a `[[f64; 3]; 8]` array.
/// * `scaling_factors`: The scaling factors in x, y, and z directions as a `[f64; 3]` array.
/// * `local_coordinate_system`: A 4x3 array defining the local coordinate system of the block.
#[derive(Debug, Clone)]
pub struct DesignBlock {
    pub nodes: Vec<[f64; 3]>,
    pub length: [f64; 3],
    pub centre: [f64; 3],
    pub theta: [f64; 3],
    pub resolution: [usize; 3],
    pub corner_nodes: [[f64; 3]; 8],
    pub scaling_factors: [f64; 3],
    pub local_coordinate_system: [[f64; 3]; 4],
}

impl DesignBlock {
    /// Creates a new `DesignBlock` instance.
    ///
    /// # Arguments
    /// * `length`: The lengths of the block in the x, y, and z directions as a `[f64; 3]` array.
    /// * `centre`: The coordinates of the center of the block as a `[f64; 3]` array.
    /// * `theta`: The angles (in radians) to rotate the block around the x, y, and z axes as a `[f64; 3]` array.
    /// * `resolution`: The resolution of the block in each dimension, represented as a `[usize; 3]` array.
    ///
    /// # Returns
    /// A new instance of `DesignBlock`.
    pub fn new(
        length: [f64; 3],
        centre: [f64; 3],
        theta: [f64; 3],
        resolution: [usize; 3],
    ) -> Self {
        let nodes = generate_block_cluster(length, centre, theta, resolution);

        // Evaluate scaling factors
        let scaling_factors = [1. / length[0], 1. / length[1], 1. / length[2]];
        // let mut min_values = [f64::INFINITY; 3];
        // let mut max_values = [f64::NEG_INFINITY; 3];

        // for point in nodes.iter() {
        //     for (i, &val) in point.iter().enumerate() {
        //         if val < min_values[i] {
        //             min_values[i] = val;
        //         }
        //         if val > max_values[i] {
        //             max_values[i] = val;
        //         }
        //     }
        // }

        // let scaling_factors = [
        //     1. / (max_values[0] - min_values[0]),
        //     1. / (max_values[1] - min_values[1]),
        //     1. / (max_values[2] - min_values[2]),
        // ];

        // Evaluate corner points
        let [cx, cy, cz] = centre;
        let [lx, ly, lz] = length;

        let half_lx = lx / 2.0;
        let half_ly = ly / 2.0;
        let half_lz = lz / 2.0;

        let mut corner_nodes = [
            [cx - half_lx, cy - half_ly, cz - half_lz], // Corner 0
            [cx + half_lx, cy - half_ly, cz - half_lz], // Corner 1
            [cx - half_lx, cy + half_ly, cz - half_lz], // Corner 2
            [cx + half_lx, cy + half_ly, cz - half_lz], // Corner 3
            [cx - half_lx, cy - half_ly, cz + half_lz], // Corner 4
            [cx + half_lx, cy - half_ly, cz + half_lz], // Corner 5
            [cx - half_lx, cy + half_ly, cz + half_lz], // Corner 6
            [cx + half_lx, cy + half_ly, cz + half_lz], // Corner 7
        ];

        rotate_nodes(&mut corner_nodes, &theta, &centre);

        let local_coordinate_system = [
            corner_nodes[0],
            corner_nodes[1],
            corner_nodes[2],
            corner_nodes[4],
        ];

        let updated_coordinate_system: [[f64; 3]; 4] = {
            let mut result = [[0.0; 3]; 4];
            for (i, &node) in local_coordinate_system.iter().enumerate() {
                result[i] = [
                    node[0] - corner_nodes[0][0],
                    node[1] - corner_nodes[0][1],
                    node[2] - corner_nodes[0][2],
                ];
            }
            result
        };

        DesignBlock {
            nodes,
            length,
            centre,
            theta,
            resolution,
            corner_nodes,
            scaling_factors,
            local_coordinate_system: updated_coordinate_system,
        }
    }

    /// Selects free deformable control nodes by checking intersection with a mesh body.
    ///
    /// # Arguments
    /// * `target_mesh`: Mesh to check intersections with.
    /// * `fixed_layers`: How many layers to exclude from the intersection, for instance quadratic=2, linear=1.
    ///
    /// # Returns
    /// Node ids corresponding to free deformable `DesignBlock` nodes.
    pub fn select_free_design_nodes(
        &self,
        target_mesh: &Mesh,
        fixed_layers: Option<usize>,
    ) -> Result<Vec<usize>, String> {
        let fixed_layers = fixed_layers.unwrap_or(2);

        if self.resolution.iter().any(|&res| res < fixed_layers) {
            return Err(String::from(
                "Block resultion must be at least the size of the fixed layers!",
            ));
        }

        // Get sides
        let corner = self.corner_nodes;

        let sides = [
            [corner[0], corner[1], corner[3], corner[2]], // Side 0 (Front)
            [corner[4], corner[5], corner[7], corner[6]], // Side 1 (Back)
            [corner[0], corner[1], corner[5], corner[4]], // Side 2 (Bottom)
            [corner[2], corner[3], corner[7], corner[6]], // Side 3 (Top)
            [corner[0], corner[2], corner[6], corner[4]], // Side 4 (Left)
            [corner[1], corner[3], corner[7], corner[5]], // Side 5 (Right)
        ];

        let mut fixed_sides: Vec<usize> = Vec::new();
        for (i, side) in sides.iter().enumerate() {
            let a = [
                side[1][0] - side[0][0],
                side[1][1] - side[0][1],
                side[1][2] - side[0][2],
            ];
            let b = [
                side[2][0] - side[0][0],
                side[2][1] - side[0][1],
                side[2][2] - side[0][2],
            ];

            let normal = [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ];

            let mut centre = [0.0; 3];
            for point in side.iter() {
                centre[0] += point[0];
                centre[1] += point[1];
                centre[2] += point[2];
            }
            centre[0] /= 4.0;
            centre[1] /= 4.0;
            centre[2] /= 4.0;

            // Intersection found
            if target_mesh.slice(centre, normal).is_some() {
                fixed_sides.push(i)
            }
        }

        // Extract valid block nodes
        let mut new_length = self.length;
        let mut new_centre = self.centre;

        let step = [
            self.length[0] / (self.resolution[0] as f64),
            self.length[1] / (self.resolution[1] as f64),
            self.length[2] / (self.resolution[2] as f64),
        ];

        let correction_factor = 0.01;
        let correction = [
            step[0] * (fixed_layers - 1) as f64 + correction_factor * step[0],
            step[1] * (fixed_layers - 1) as f64 + correction_factor * step[1],
            step[2] * (fixed_layers - 1) as f64 + correction_factor * step[2],
        ];

        // Left check
        if fixed_sides.contains(&4) {
            new_length[0] -= correction[0];
            new_centre[0] += 0.5 * correction[0];
        }
        // Right check
        if fixed_sides.contains(&5) {
            new_length[0] -= correction[0];
            new_centre[0] -= 0.5 * correction[0];
        }

        // Bottom check
        if fixed_sides.contains(&2) {
            new_length[1] -= correction[1];
            new_centre[1] += 0.5 * correction[1];
        }
        // Top check
        if fixed_sides.contains(&3) {
            new_length[1] -= correction[1];
            new_centre[1] -= 0.5 * correction[1];
        }
        // Front check
        if fixed_sides.contains(&0) {
            new_length[2] -= correction[2];
            new_centre[2] += 0.5 * correction[2];
        }
        // Back check
        if fixed_sides.contains(&1) {
            new_length[2] -= correction[2];
            new_centre[2] -= 0.5 * correction[2];
        }

        new_length[0] += 0.1 * correction_factor * step[0];
        new_length[1] += 0.1 * correction_factor * step[1];
        new_length[2] += 0.1 * correction_factor * step[2];

        let rotation_matrix = build_rotation_matrix(&self.theta);
        rotate_node(&mut new_centre, &rotation_matrix, &self.centre);

        let free_node_ids =
            select_nodes_in_box(&self.nodes, new_length, new_centre, self.theta);

        Ok(free_node_ids)
    }

    /// Selects the node ids corresponding to the target mesh that fit inside the design block.
    ///
    /// # Arguments
    /// * `target_mesh_nodes`: Mesh nodes to check.
    ///
    /// # Returns
    /// Node ids corresponding to free deformable `Mesh` nodes.
    pub fn select_mesh_nodes_inside(
        &self,
        target_mesh_nodes: &[[f64; 3]],
    ) -> Result<Vec<usize>, String> {
        Ok(select_nodes_in_box(
            target_mesh_nodes,
            self.length,
            self.centre,
            self.theta,
        ))
    }
}

/// Creates a local coordinate system in 3D space given lengths, centre, and rotations.
/// The local coordinate system is composed of four points:
/// - Origin node (centre)
/// - X-axis node
/// - Y-axis node
/// - Z-axis node
///
/// These nodes represent the coordinate system after applying the specified
/// lengths and rotations around the centre.
///
/// # Arguments
/// * `length`: An array `[x_length, y_length, z_length]` that specifies the length of each axis.
/// * `centre`: An array `[x, y, z]` specifying the centre point of the orientation frame.
/// * `theta`: An array `[theta_x, theta_y, theta_z]` specifying rotations in radians for each axis.
///
/// # Returns
/// Returns a `Vec<[f64; 3]>` containing four nodes, each represented as `[x, y, z]`.
///
/// # Examples
/// ```
/// let frame = evaluate_lcs_from_components([1.0, 1.0, 1.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
/// ```
pub fn evaluate_lcs_from_components(
    length: [f64; 3],
    centre: [f64; 3],
    theta: [f64; 3],
) -> [[f64; 3]; 4] {
    let origin = [
        0.5 * length[0] - centre[0],
        0.5 * length[1] - centre[1],
        0.5 * length[2] - centre[2],
    ];
    let mut local_coordinate_system = [
        [centre[0], centre[1], centre[2]],
        [centre[0] + length[0], centre[1], centre[2]],
        [centre[0], centre[1] + length[1], centre[2]],
        [centre[0], centre[1], centre[2] + length[2]],
    ];

    rotate_nodes(&mut local_coordinate_system, &theta, &origin);

    local_coordinate_system
}
