use std::f64::consts::PI;

use crate::core::algebra::{mat3_mat3_mul, mat3_vec3_mul, mat4_vec4_mul};

/// Transforms a list of 3D points using a 4x4 transformation matrix.
///
/// # Arguments
/// * `nodes`: A mutable vector of 3D nodes represented as arrays of f64.
/// * `transformation_matrix`: A 4x4 array representing the transformation matrix.
/// * `origin`: A reference to an `[f64; 3]` specifying the origin of rotation.
pub fn transform_nodes(
    nodes: &mut [[f64; 3]],
    transformation_matrix: &[[f64; 4]; 4],
    origin: &[f64; 3],
) {
    for node in nodes.iter_mut() {
        transform_node(node, transformation_matrix, origin)
    }
}

/// Creates a new set of nodes by applying a transformation to a specified subset.
///
/// This function is non-mutating; it returns a new vector with the transformed points.
///
/// # Arguments
/// * `original_nodes`: A slice of the initial node positions.
/// * `node_ids`: A slice of indices specifying which nodes to transform.
/// * `transform_matrix`: The 4x4 affine transformation matrix to apply.
/// * `origin`: The point around which the transformation (e.g., rotation) should occur.
///
/// # Returns
/// A new `Vec<[f64; 3]>` containing the full set of nodes with the specified subset transformed.
pub fn transform_selected_nodes(
    original_nodes: &[[f64; 3]],
    node_ids: &[usize],
    transform_matrix: &[[f64; 4]; 4],
    origin: &[f64; 3],
) -> Vec<[f64; 3]> {
    let mut transformed_nodes = original_nodes.to_vec();

    for &id in node_ids {
        if id < transformed_nodes.len() {
            transform_node(&mut transformed_nodes[id], transform_matrix, origin);
        }
    }

    transformed_nodes
}

/// Transforms a node using a 4x4 transformation matrix.
///
/// # Arguments
/// * `node`: A mutable 3D node.
/// * `transformation_matrix`: A 4x4 array representing the transformation matrix.
/// * `origin`: A reference to an `[f64; 3]` specifying the origin of rotation.
pub fn transform_node(
    node: &mut [f64; 3],
    transformation_matrix: &[[f64; 4]; 4],
    origin: &[f64; 3],
) {
    let new_node = [
        node[0] - origin[0],
        node[1] - origin[1],
        node[2] - origin[2],
        1.0,
    ];

    let product = mat4_vec4_mul(transformation_matrix, &new_node);

    node[0] = product[0] + origin[0];
    node[1] = product[1] + origin[1];
    node[2] = product[2] + origin[2];
}

/// Creates a 4x4 transformation matrix from translation, rotation, and scaling vectors.
///
/// # Arguments
/// * `translation`: An array representing the translation.
/// * `rotation`: An array representing the rotation angles.
/// * `scaling`: An array representing the scaling factors.
///
/// # Returns
/// * Returns a 4x4 array representing the transformation matrix.
pub fn build_transformation_matrix(
    translation: [f64; 3],
    rotation: [f64; 3],
    scaling: [f64; 3],
) -> [[f64; 4]; 4] {
    let mut tmat: [[f64; 4]; 4] = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ];

    // Apply translation
    tmat[0][3] = translation[0];
    tmat[1][3] = translation[1];
    tmat[2][3] = translation[2];

    // Apply rotation
    let rot_mat = build_rotation_matrix(&rotation);

    for i in 0..3 {
        for j in 0..3 {
            tmat[i][j] = rot_mat[i][j];
        }
    }

    // Apply scaling
    for row in tmat.iter_mut().take(3) {
        row[0] *= scaling[0];
        row[1] *= scaling[1];
        row[2] *= scaling[2];
    }

    tmat
}

/// Translates a list of nodes by the given position.
///
/// # Arguments
/// * `nodes`: A mutable vector of 3D nodes represented as arrays of f64.
/// * `translation`: An array of f64 that represents the translation.
pub fn translate_nodes(nodes: &mut [[f64; 3]], translation: &[f64; 3]) {
    for node in nodes.iter_mut() {
        for j in 0..3 {
            node[j] += translation[j];
        }
    }
}

/// Rotate a set of nodes around an origin by given angles.
///
/// # Arguments
/// * `nodes`: A mutable reference to a `Vec<[f64; 3]>` containing nodes to rotate.
/// * `theta`: A reference to an `[f64; 3]` specifying rotation angles for each axis in radians.
/// * `origin`: A reference to an `[f64; 3]` specifying the origin of rotation.
pub fn rotate_nodes(nodes: &mut [[f64; 3]], theta: &[f64; 3], origin: &[f64; 3]) {
    let rotation_matrix = build_rotation_matrix(theta);

    nodes
        .iter_mut()
        .for_each(|node| rotate_node(node, &rotation_matrix, origin));
}

/// Rotate a node in 3D space using a given rotation matrix and an origin point.
///
/// # Arguments
/// * `node`: Mutable reference to the node to be rotated.
///             This node will be modified in-place.
/// * `rotation_matrix`: Reference to a 3x3 rotation matrix.
/// * `origin`: Reference to the origin point around which the node will be rotated.
///
/// # Examples
/// ```
/// let mut node = [1.0, 0.0, 0.0];
/// let rotation_matrix = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
/// let origin = [0.0, 0.0, 0.0];
///
/// rotate_node(&mut node, &rotation_matrix, &origin);
/// ```
pub fn rotate_node(
    node: &mut [f64; 3],
    rotation_matrix: &[[f64; 3]; 3],
    origin: &[f64; 3],
) {
    node[0] -= origin[0];
    node[1] -= origin[1];
    node[2] -= origin[2];

    let rotated_node = mat3_vec3_mul(rotation_matrix, node);

    node[0] = rotated_node[0] + origin[0];
    node[1] = rotated_node[1] + origin[1];
    node[2] = rotated_node[2] + origin[2];
}

/// Computes the rotation matrix from Euler angles.
///
/// # Arguments
/// * `theta`: A list of Euler angles [alpha, beta, gamma] in degrees.
///
/// # Returns
/// * A 3x3 rotation matrix represented as [[f64; 3]; 3].
///
/// # Example
/// ```
/// let angles = [45.0, 30.0, 60.0];
/// let r_matrix = build_rotation_matrix(angles);
/// ```
pub fn build_rotation_matrix(theta: &[f64; 3]) -> [[f64; 3]; 3] {
    let (alpha, beta, gamma) = (
        theta[0] * PI / 180.0,
        theta[1] * PI / 180.0,
        theta[2] * PI / 180.0,
    );

    let mut rotation_matrix = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    if gamma != 0.0 {
        let (c, s) = (gamma.cos(), gamma.sin());
        let r = [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]];
        rotation_matrix = mat3_mat3_mul(&r, &rotation_matrix);
    }
    if beta != 0.0 {
        let (c, s) = (beta.cos(), beta.sin());
        let r = [[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]];
        rotation_matrix = mat3_mat3_mul(&r, &rotation_matrix);
    }
    if alpha != 0.0 {
        let (c, s) = (alpha.cos(), alpha.sin());
        let r = [[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]];
        rotation_matrix = mat3_mat3_mul(&r, &rotation_matrix);
    }
    rotation_matrix
}

/// Scales a list of nodes by the given scaling factor, relative to a given origin.
///
/// # Arguments
/// * `nodes`: A mutable vector of 3D nodes represented as arrays of f64.
/// * `scaling`: An array of f64 that represents the scaling factors for each dimension.
/// * `origin`: The point relative to which the scaling should be performed.
pub fn scale_nodes(nodes: &mut [[f64; 3]], scaling: &[f64; 3], origin: &[f64; 3]) {
    for node in nodes.iter_mut() {
        scale_node(node, scaling, origin);
    }
}

/// Scales a node by the given scaling factor, relative to a given origin.
///
/// # Arguments
/// * `node`: A mutable 3D node.
/// * `scaling`: An array of f64 that represents the scaling factors for each dimension.
/// * `origin`: The point relative to which the scaling should be performed.
pub fn scale_node(node: &mut [f64; 3], scaling: &[f64; 3], origin: &[f64; 3]) {
    for j in 0..3 {
        node[j] -= origin[j];
        node[j] *= scaling[j];
        node[j] += origin[j];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_transform_nodes() {
        let transformation_matrix =
            build_transformation_matrix([1.3, 15., 2.5], [25., 25., 25.], [2., 5., 2.]);

        let mut nodes = vec![[1.4, 1.2, 3.5], [34., 23., 53.]];

        transform_nodes(&mut nodes, &transformation_matrix, &[0., 0., 0.]);

        let expected_nodes = vec![
            [4.260097156389118, 18.32001817946411, 11.047239560979023],
            [57.90475899451524, 97.23229418127102, 140.77057189748564],
        ];

        let mut equal = true;
        'outer: for i in 0..2 {
            for j in 0..3 {
                if (nodes[i][j] - expected_nodes[i][j]).abs() > 1e-5 {
                    equal = false;
                    break 'outer;
                }
            }
        }

        assert!(equal);
    }

    #[test]
    fn test_rotation_matrix() {
        let theta: [f64; 3] = [45.0, 45.0, 45.0];
        let rotation_matrix = build_rotation_matrix(&theta);

        let expected_matrix = [
            [0.5, -0.5, 0.70710678],
            [0.85355339, 0.14644661, -0.5],
            [0.14644661, 0.85355339, 0.5],
        ];

        let mut equal = true;
        'outer: for i in 0..3 {
            for j in 0..3 {
                if (rotation_matrix[i][j] - expected_matrix[i][j]).abs() > 1e-5 {
                    equal = false;
                    break 'outer;
                }
            }
        }

        assert!(equal);
    }

    #[test]
    fn test_rotate_node() {
        let mut node: [f64; 3] = [1.0, 0.0, 0.0];
        let origin: [f64; 3] = [0.0, 0.0, 0.0];
        let rotation_matrix = [
            [0.5, -0.5, 0.70710678],
            [0.85355339, 0.14644661, -0.5],
            [0.14644661, 0.85355339, 0.5],
        ];

        rotate_node(&mut node, &rotation_matrix, &origin);
        let expected_node: [f64; 3] = [0.5, 0.85355339, 0.14644661];

        assert!((node[0] - expected_node[0]).abs() < 1e-5);
        assert!((node[1] - expected_node[1]).abs() < 1e-5);
        assert!((node[2] - expected_node[2]).abs() < 1e-5);
    }

    #[test]
    fn test_rotate_nodes() {
        let mut nodes = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let theta = [0.0, 0.0, 90.0];
        let origin = [0.0, 0.0, 0.0];

        rotate_nodes(&mut nodes, &theta, &origin);

        let expected_nodes = vec![[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]];

        let mut equal = true;
        'outer: for i in 0..nodes.len() {
            for j in 0..3 {
                if (nodes[i][j] - expected_nodes[i][j]).abs() > 1e-5 {
                    equal = false;
                    break 'outer;
                }
            }
        }

        assert!(equal);
    }

    #[test]
    fn test_translate_nodes() {
        let mut nodes = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]];
        let translation = [1.0, 1.0, 1.0];
        translate_nodes(&mut nodes, &translation);
        assert_eq!(nodes, vec![[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]]);
    }

    #[test]
    fn test_scale_nodes() {
        let mut nodes = vec![[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]];
        let origin = [1.0, 1.0, 1.0];
        let scaling = [2.0, 2.0, 2.0];
        scale_nodes(&mut nodes, &scaling, &origin);
        assert_eq!(nodes, vec![[1.0, 1.0, 1.0], [3.0, 3.0, 3.0]]);
    }
}
