#[cfg(feature = "rayon")]
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::error::Result;

/// Apply an affine transformation to a set of 3D points.
///
/// # Arguments
/// * `points`: A vector of points in 3D space, represented as arrays `[x, y, z]`.
/// * `affine_weights`: A 4x4 array representing the affine transformation weights.
///
/// # Returns
/// * `Result<Vec<[f64; 3]>>`: A Result containing either:
/// A vector of transformed points (`Ok`)
/// An error if something goes wrong (`Err`)
pub fn apply_affine_transform(
    points: &[[f64; 3]],
    affine_weights: &[[f64; 4]; 4],
) -> Result<Vec<[f64; 3]>> {
    let padded_points = points
        .iter()
        .map(|&[x1, x2, x3]| vec![x1, x2, x3, 1.0])
        .collect::<Vec<Vec<f64>>>();

    let mut result = Vec::new();
    for point in padded_points {
        let mut transformed = [0.0; 4];
        for i in 0..4 {
            for j in 0..4 {
                transformed[i] += point[j] * affine_weights[j][i];
            }
        }
        result.push([transformed[0], transformed[1], transformed[2]]);
    }

    Ok(result)
}

/// Apply a Bernstein transform to a set of 3D points.
///
/// # Arguments
/// * `points`: A vector of points in 3D space, represented as arrays of f64 numbers `[x, y, z]`.
/// * `deltas`: A vector of delta shifts for each point, also in 3D `[dx, dy, dz]`.
/// * `resolution`: An array specifying the resolution in each dimension `[res_x, res_y, res_z]`.
///
/// # Returns
///
/// * `Result<Vec<[f64; 3]>, String>`: A Result containing either:
/// A vector of transformed points in the same format as the input (`Ok`)
/// An error message if something goes wrong (`Err`)
pub fn apply_bernstein_transform(
    points: &[[f64; 3]],
    deltas: &[[f64; 3]],
    resolution: &[usize; 3],
) -> Result<Vec<[f64; 3]>> {
    let dimension = [resolution[0] + 1, resolution[1] + 1, resolution[2] + 1];

    // Pre-compute all binomial coefficients once (common logic)
    let coeffs_x: Vec<f64> = (0..dimension[0])
        .map(|i| binomial_coefficient(dimension[0] - 1, i))
        .collect();
    let coeffs_y: Vec<f64> = (0..dimension[1])
        .map(|j| binomial_coefficient(dimension[1] - 1, j))
        .collect();
    let coeffs_z: Vec<f64> = (0..dimension[2])
        .map(|k| binomial_coefficient(dimension[2] - 1, k))
        .collect();

    // Process all points using the appropriate iterator
    #[cfg(feature = "rayon")]
    let transformed_points: Vec<[f64; 3]> = points
        .par_iter()
        .map(|point| {
            transform_single_point(
                point, deltas, &dimension, &coeffs_x, &coeffs_y, &coeffs_z,
            )
        })
        .collect();

    #[cfg(not(feature = "rayon"))]
    let transformed_points: Vec<[f64; 3]> = points
        .iter()
        .map(|point| {
            transform_single_point(
                point, deltas, &dimension, &coeffs_x, &coeffs_y, &coeffs_z,
            )
        })
        .collect();

    Ok(transformed_points)
}

/// Helper function containing the core logic to transform a single point.
fn transform_single_point(
    point: &[f64; 3],
    deltas: &[[f64; 3]],
    dimension: &[usize; 3],
    coeffs_x: &[f64],
    coeffs_y: &[f64],
    coeffs_z: &[f64],
) -> [f64; 3] {
    // Pre-compute 1D Bernstein basis values for this point
    let bernstein_x: Vec<f64> = (0..dimension[0])
        .map(|i| {
            let p = point[0];
            coeffs_x[i] * (1.0 - p).powi((dimension[0] - 1 - i) as i32) * p.powi(i as i32)
        })
        .collect();

    let bernstein_y: Vec<f64> = (0..dimension[1])
        .map(|j| {
            let p = point[1];
            coeffs_y[j] * (1.0 - p).powi((dimension[1] - 1 - j) as i32) * p.powi(j as i32)
        })
        .collect();

    let bernstein_z: Vec<f64> = (0..dimension[2])
        .map(|k| {
            let p = point[2];
            coeffs_z[k] * (1.0 - p).powi((dimension[2] - 1 - k) as i32) * p.powi(k as i32)
        })
        .collect();

    // Perform the summation using the pre-computed values
    let mut aux_shift = [0.0; 3];
    for i in 0..dimension[0] {
        for j in 0..dimension[1] {
            for k in 0..dimension[2] {
                let bernstein_prod = bernstein_x[i] * bernstein_y[j] * bernstein_z[k];
                let delta_id = i * dimension[1] * dimension[2] + j * dimension[2] + k;
                let delta = deltas[delta_id];

                aux_shift[0] += bernstein_prod * delta[0];
                aux_shift[1] += bernstein_prod * delta[1];
                aux_shift[2] += bernstein_prod * delta[2];
            }
        }
    }

    // Add the final shift to the original point
    [
        point[0] + aux_shift[0],
        point[1] + aux_shift[1],
        point[2] + aux_shift[2],
    ]
}

/// Compute the binomial coefficient "n choose k".
///
/// # Arguments
/// * `n`: The total number of items.
/// * `k`: The number of items to choose.
///
/// # Returns
/// * `f64`: The computed binomial coefficient.
fn binomial_coefficient(n: usize, k: usize) -> f64 {
    let mut coeff = 1.0;
    for i in 0..k {
        coeff *= (n - i) as f64 / (k - i) as f64;
    }
    coeff
}

#[cfg(test)]
mod tests {
    use super::*;
    const EPSILON: f64 = 1e-9;

    #[test]
    fn test_affine_identity() {
        let points = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let identity_matrix = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let transformed = apply_affine_transform(&points, &identity_matrix).unwrap();
        assert_eq!(
            points, transformed,
            "Identity matrix should not change points"
        );
    }

    #[test]
    fn test_affine_translation() {
        let points = vec![[1.0, 2.0, 3.0]];
        let translation_matrix = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [10.0, -5.0, 2.0, 1.0], // Translate by (10, -5, 2)
        ];
        let transformed = apply_affine_transform(&points, &translation_matrix).unwrap();
        let expected = vec![[11.0, -3.0, 5.0]];

        for i in 0..3 {
            assert!(
                (transformed[0][i] - expected[0][i]).abs() < EPSILON,
                "Translation failed at index {}",
                i
            );
        }
    }

    #[test]
    fn test_bernstein_zero_deltas() {
        let points = vec![[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]];
        let resolution = [2, 2, 2];
        let num_deltas = (resolution[0] + 1) * (resolution[1] + 1) * (resolution[2] + 1);
        let deltas = vec![[0.0, 0.0, 0.0]; num_deltas];

        let transformed =
            apply_bernstein_transform(&points, &deltas, &resolution).unwrap();
        assert_eq!(
            points, transformed,
            "Zero deltas should result in no change"
        );
    }

    #[test]
    fn test_bernstein_simple_linear() {
        let resolution = [1, 1, 1];
        let point = vec![[0.5, 0.5, 0.5]];
        let num_deltas = 2 * 2 * 2;
        let deltas = vec![[1.0, 2.0, 3.0]; num_deltas];

        let transformed =
            apply_bernstein_transform(&point, &deltas, &resolution).unwrap();

        let expected = vec![[0.5 + 1.0, 0.5 + 2.0, 0.5 + 3.0]];

        for i in 0..3 {
            assert!(
                (transformed[0][i] - expected[0][i]).abs() < EPSILON,
                "Linear Bernstein failed at index {}",
                i
            );
        }
    }
}
