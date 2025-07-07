use crate::core::transformation::{build_rotation_matrix, rotate_node};
use crate::error::{Error, Result};

/// Creates a new 3D block of nodes.
///
/// # Arguments
/// * `length`: A `[f64; 3]` specifying the length of the grid in the x, y, and z directions.
/// * `centre`: A `[f64; 3]` specifying the coordinates of the centre of the grid.
/// * `theta`: A `[f64; 3]` specifying the angles of rotation of the grid.
/// * `resolution`: A `[usize; 3]` specifying the number of divisions in the x, y, and z directions.
///
/// # Returns
/// A new `Result<Vec<[f64; 3]>>` instance.
pub fn generate_block_cluster(
    length: [f64; 3],
    centre: [f64; 3],
    theta: [f64; 3],
    resolution: [usize; 3],
) -> Result<Vec<[f64; 3]>> {
    for &s in &length {
        if s <= 0.0 {
            return Err(Error::MeshGeneration(
                "Invalid length: dimensions must be positive.".to_string(),
            ));
        }
    }
    for &r in &resolution {
        if r == 0 {
            return Err(Error::MeshGeneration(
                "Invalid resolution: dimensions must be non-zero.".to_string(),
            ));
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

    Ok(nodes)
}

/// Creates a new solid sphere of nodes.
///
/// Generates a grid of points within a bounding box and keeps only those
/// that fall within the sphere's radius.
///
/// # Arguments
/// * `radius`: The radius of the sphere.
/// * `centre`: A `[f64; 3]` specifying the coordinates of the centre of the sphere.
/// * `resolution`: An integer specifying the number of divisions along the bounding box axes.
///
/// # Returns
/// A new `Result<Vec<[f64; 3]>>` instance.
pub fn generate_sphere_cluster(
    radius: f64,
    centre: [f64; 3],
    resolution: usize,
) -> Result<Vec<[f64; 3]>> {
    if radius <= 0.0 {
        return Err(Error::MeshGeneration(
            "Invalid radius: must be positive.".to_string(),
        ));
    }
    if resolution == 0 {
        return Err(Error::MeshGeneration(
            "Invalid resolution: must be non-zero.".to_string(),
        ));
    }

    let mut nodes = Vec::new();
    let length = radius * 2.0;
    let step = length / (resolution as f64);
    let start = [centre[0] - radius, centre[1] - radius, centre[2] - radius];

    for i in 0..=resolution {
        for j in 0..=resolution {
            for k in 0..=resolution {
                let point = [
                    start[0] + i as f64 * step,
                    start[1] + j as f64 * step,
                    start[2] + k as f64 * step,
                ];

                // Check if the point is inside the sphere
                let dist_sq = (point[0] - centre[0]).powi(2)
                    + (point[1] - centre[1]).powi(2)
                    + (point[2] - centre[2]).powi(2);

                if dist_sq <= radius.powi(2) {
                    nodes.push(point);
                }
            }
        }
    }

    Ok(nodes)
}

#[cfg(test)]
mod tests {
    use super::*;
    const EPSILON: f64 = 1e-9;

    #[test]
    fn test_generate_block_cluster() {
        let length = [2.0, 2.0, 2.0];
        let centre = [1.0, 1.0, 1.0];
        let theta = [0.0, 0.0, 0.0];
        let resolution = [2, 2, 2];
        let nodes = generate_block_cluster(length, centre, theta, resolution).unwrap();

        // Expected nodes = (rx+1) * (ry+1) * (rz+1)
        assert_eq!(nodes.len(), 3 * 3 * 3);

        // Check corner nodes
        // Min corner should be at centre - length/2 = [0,0,0]
        assert!((nodes[0][0] - 0.0).abs() < EPSILON);
        assert!((nodes[0][1] - 0.0).abs() < EPSILON);
        assert!((nodes[0][2] - 0.0).abs() < EPSILON);

        // Max corner should be at centre + length/2 = [2,2,2]
        let last_node = nodes.last().unwrap();
        assert!((last_node[0] - 2.0).abs() < EPSILON);
        assert!((last_node[1] - 2.0).abs() < EPSILON);
        assert!((last_node[2] - 2.0).abs() < EPSILON);

        // Invalid length
        let result =
            generate_block_cluster([0.0, 1.0, 1.0], [0.0; 3], [0.0; 3], [1, 1, 1]);
        assert!(result.is_err());
        assert!(matches!(result, Err(Error::MeshGeneration(_))));

        // Invalid resolution
        let result =
            generate_block_cluster([1.0, 1.0, 1.0], [0.0; 3], [0.0; 3], [1, 0, 1]);
        assert!(result.is_err());
        assert!(matches!(result, Err(Error::MeshGeneration(_))));
    }

    #[test]
    fn test_generate_sphere_cluster() {
        let radius = 1.0;
        let centre = [0.0, 0.0, 0.0];
        let resolution = 10;
        let nodes = generate_sphere_cluster(radius, centre, resolution).unwrap();

        // The number of nodes should be > 0 but < total points in bounding box
        let bounding_box_nodes = (resolution + 1).pow(3);
        assert!(!nodes.is_empty());
        assert!(nodes.len() < bounding_box_nodes);

        // Check that all points are within the radius
        for node in nodes {
            let dist_sq = node[0].powi(2) + node[1].powi(2) + node[2].powi(2);
            assert!(dist_sq <= radius.powi(2) + EPSILON);
        }

        // Invalid radius
        let result = generate_sphere_cluster(-1.0, [0.0; 3], 10);
        assert!(result.is_err());
        assert!(matches!(result, Err(Error::MeshGeneration(_))));

        // Invalid resolution
        let result = generate_sphere_cluster(1.0, [0.0; 3], 0);
        assert!(result.is_err());
        assert!(matches!(result, Err(Error::MeshGeneration(_))));
    }
}
