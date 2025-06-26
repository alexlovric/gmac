use crate::{
    core::utilities::{squared_euclidean_distance_3d, euclidean_distance_3d},
    core::transformation::{build_rotation_matrix, rotate_node},
};

/// Selects the indices of nodes closest to a point.
///
/// # Arguments
/// * `nodes`: The list of 3D points to select from.
/// * `point`: The reference point.
///
/// # Returns
/// A vector of indices corresponding to the nodes closest to the point.
pub fn select_nodes_closest_to_point(nodes: &[[f64; 3]], point: [f64; 3]) -> Vec<usize> {
    let mut min_distance = f64::INFINITY;
    let mut closest_indices: Vec<usize> = Vec::new();

    for (i, node) in nodes.iter().enumerate() {
        let distance = squared_euclidean_distance_3d(*node, point);

        if (distance - min_distance).abs() <= f64::EPSILON {
            closest_indices.push(i);
        } else if distance < min_distance {
            closest_indices.clear();
            closest_indices.push(i);
            min_distance = distance;
        }
    }
    closest_indices
}

/// Selects the indices of nodes that pass through a line provided by two points.
///
/// # Arguments
/// * `nodes`: The list of 3D points to select from.
/// * `point_a`: The first endpoint of the line.
/// * `point_b`: The second endpoint of the line.
///
/// # Returns
/// A vector of indices corresponding to the nodes passing through the line.
pub fn select_nodes_on_line(
    nodes: &[[f64; 3]],
    point_a: [f64; 3],
    point_b: [f64; 3],
) -> Vec<usize> {
    let ab = euclidean_distance_3d(point_a, point_b);
    (0..nodes.len())
        .filter(|&id| {
            let node = nodes[id];
            let ap = euclidean_distance_3d(point_a, node);
            let bp = euclidean_distance_3d(point_b, node);
            (ap + bp - ab).abs() <= f64::EPSILON && ap <= ab && bp <= ab
        })
        .collect()
}

/// Selects the indices of nodes in a sphere with a radius and origin as provided.
///
/// # Arguments
/// * `nodes`: The list of 3D points to select from.
/// * `radius`: The radius of the sphere.
/// * `centre`: The centre of the sphere.
///
/// # Returns
/// A vector of indices corresponding to the nodes in the sphere.
pub fn select_nodes_in_sphere(
    nodes: &[[f64; 3]],
    radius: f64,
    centre: [f64; 3],
) -> Vec<usize> {
    let radius_squared = radius * radius;
    (0..nodes.len())
        .filter(|&id| {
            let node = nodes[id];
            squared_euclidean_distance_3d(node, centre) <= radius_squared
        })
        .collect()
}

/// Selects the indices of the nodes in a box with dimensions
/// and origin as provided.
///
/// # Arguments
/// * `nodes`: The list of 3D points to select from.
/// * `length`: The dimensions of the box [length, width, height].
/// * `centre`: The centre of the box.
///
/// # Returns
/// A vector of indices corresponding to the nodes in the box.
pub fn select_nodes_in_box(
    nodes: &[[f64; 3]],
    length: [f64; 3],
    centre: [f64; 3],
    theta: [f64; 3],
) -> Vec<usize> {
    let reversed_theta = [-theta[0], -theta[1], -theta[2]];
    let rotation_matrix = build_rotation_matrix(&reversed_theta);

    (0..nodes.len())
        .filter(|&id| {
            let mut node = nodes[id];

            rotate_node(&mut node, &rotation_matrix, &centre);

            (node[0] >= centre[0] - length[0] * 0.5
                && node[0] <= centre[0] + length[0] * 0.5)
                && (node[1] >= centre[1] - length[1] * 0.5
                    && node[1] <= centre[1] + length[1] * 0.5)
                && (node[2] >= centre[2] - length[2] * 0.5
                    && node[2] <= centre[2] + length[2] * 0.5)
        })
        .collect()
}

/// Selects the node indices in the mesh closest to the specified plane.
///
/// # Arguments
/// * `nodes`: The list of 3D points to select from.
/// * `origin`: An array of x, y, z coordinates representing a point on the plane.
/// * `normal`: An array representing the normal vector to the plane.
///
/// # Returns
/// A vector of references to nodes that are closest to the plane.
pub fn select_nodes_closest_to_plane(
    nodes: &[[f64; 3]],
    origin: [f64; 3],
    normal: [f64; 3],
) -> Vec<usize> {
    // Calculate constant term 'd' for the plane equation
    let d = origin
        .iter()
        .zip(normal.iter())
        .map(|(o, n)| o * n)
        .sum::<f64>();

    let mut closest_node_ids: Vec<usize> = Vec::new();
    let mut min_distance = f64::MAX;

    for (id, node) in nodes.iter().enumerate() {
        // Evaluate the absolute distance of the node to the plane using the plane equation
        let distance = (normal
            .iter()
            .zip(node.iter())
            .map(|(n, x)| n * x)
            .sum::<f64>()
            - d)
            .abs();

        // Check if this distance is smaller than the minimum distance
        if distance < min_distance {
            min_distance = distance;
            closest_node_ids.clear();
            closest_node_ids.push(id);
        } else if (distance - min_distance).abs() < f64::EPSILON {
            // Handle the case where multiple nodes have the same minimum distance to the plane
            closest_node_ids.push(id);
        }
    }
    closest_node_ids
}

/// Selects the indices of nodes that lie in the direction of the plane's normal.
///
/// # Arguments
/// * `nodes`: A slice of 3D points `[x, y, z]` representing the nodes of a mesh.
/// * `origin`: A 3D point `[x, y, z]` representing a point through which the slice plane passes.
/// * `normal`: A 3D vector `[x, y, z]` representing the normal to the slice plane.
///
/// # Returns
/// The indices of the nodes that are in the direction of the plane's normal.
pub fn select_nodes_in_plane_direction(
    nodes: &[[f64; 3]],
    origin: [f64; 3],
    normal: [f64; 3],
) -> Vec<usize> {
    (0..nodes.len())
        .filter(|&id| is_node_on_positive_side(nodes[id], origin, normal))
        .collect()
}

/// Determines if a given node lies in the direction of a plane's normal.
///
/// # Arguments
/// * `point`: A 3D point `[x, y, z]` representing the node.
/// * `origin`: A 3D point `[x, y, z]` representing a point through which the slice plane passes.
/// * `normal`: A 3D vector `[x, y, z]` representing the normal to the slice plane.
///
/// # Returns
/// Returns `true` if the node is on the positive side of the plane (i.e., in the direction
/// of the plane's normal). Otherwise, returns `false`.
pub fn is_node_on_positive_side(
    node: [f64; 3],
    origin: [f64; 3],
    normal: [f64; 3],
) -> bool {
    let dot_product = normal
        .iter()
        .zip(node.iter())
        .map(|(n, p)| n * p)
        .sum::<f64>()
        - origin
            .iter()
            .zip(normal.iter())
            .map(|(o, n)| o * n)
            .sum::<f64>();

    dot_product >= 0.0
}

#[cfg(test)]
mod tests {
    use crate::core::{primitives::generate_box, mesh::Mesh};
    use super::*;

    fn new_block() -> Mesh {
        generate_box([1., 1., 1.], [0.5, 0.5, 0.5], [0.0, 0.0, 0.0], [2, 2, 1])
    }

    #[test]
    fn test_select_node_closest_to() {
        let block = new_block();
        let closest_nodes = select_nodes_closest_to_point(&block.nodes, [0.4, 0.4, 0.4]);
        assert_eq!(closest_nodes[0], 8);
    }

    #[test]
    fn test_select_nodes_on_line() {
        let block = new_block();
        let selected_nodes =
            select_nodes_on_line(&block.nodes, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);

        assert_eq!(selected_nodes.len(), 3);
        assert_eq!(selected_nodes[0], 0);
        assert_eq!(selected_nodes[1], 6);
        assert_eq!(selected_nodes[2], 12);
    }

    #[test]
    fn test_select_nodes_in_sphere() {
        let block = new_block();
        let selected_nodes = select_nodes_in_sphere(&block.nodes, 0.5, [0.0, 0.0, 0.0]);

        assert_eq!(selected_nodes.len(), 3);
        assert_eq!(selected_nodes[0], 0);
        assert_eq!(selected_nodes[1], 2);
        assert_eq!(selected_nodes[2], 6);
    }

    #[test]
    fn test_select_nodes_in_box() {
        let block = new_block();
        let selected_nodes = select_nodes_in_box(
            &block.nodes,
            [1.0, 0.5, 0.5],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        );

        assert_eq!(selected_nodes.len(), 2);
        assert_eq!(selected_nodes[0], 0);
        assert_eq!(selected_nodes[1], 6);
    }
}
