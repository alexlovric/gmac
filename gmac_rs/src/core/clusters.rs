use crate::core::transformation::{build_rotation_matrix, rotate_node};

/// Creates a new 3D block of nodes.
///
/// # Arguments
/// * `length`: A `[f64; 3]` specifying the length of the grid in the x, y, and z directions.
/// * `centre`: A `[f64; 3]` specifying the coordinates of the centre of the grid.
/// * `theta`: A `[f64; 3]` specifying the angles of rotation of the grid.
/// * `resolution`: A `[usize; 3]` specifying the number of divisions in the x, y, and z directions.
///
/// # Returns
/// A new `Vec<[f64; 3]>` instance.
pub fn generate_block_cluster(
    length: [f64; 3],
    centre: [f64; 3],
    theta: [f64; 3],
    resolution: [usize; 3],
) -> Vec<[f64; 3]> {
    for &s in &length {
        if s <= 0.0 {
            panic!("Invalid length: dimensions must be positive");
        }
    }
    for &r in &resolution {
        if r == 0 {
            panic!("Invalid resolution: dimensions must be non-zero");
        }
    }

    let [rx, ry, rz] = resolution;

    let step = [
        length[0] / (rx as f64),
        length[1] / (ry as f64),
        length[2] / (rz as f64),
    ];

    let start = [
        centre[0] - length[0] * 0.5,
        centre[1] - length[1] * 0.5,
        centre[2] - length[2] * 0.5,
    ];

    let rotation_matrix = build_rotation_matrix(&theta);

    let mut nodes = vec![[0.0; 3]; (rx + 1) * (ry + 1) * (rz + 1)];

    let mut node_count = 0;
    for i in 0..=rx {
        for j in 0..=ry {
            for k in 0..=rz {
                nodes[node_count] = [
                    start[0] + i as f64 * step[0],
                    start[1] + j as f64 * step[1],
                    start[2] + k as f64 * step[2],
                ];

                if theta.iter().any(|angle| *angle != 0.0) {
                    rotate_node(&mut nodes[node_count], &rotation_matrix, &centre);
                }

                node_count += 1;
            }
        }
    }

    nodes
}
