/// Apply an affine transformation to a set of 3D points.
///
/// # Arguments
/// * `points`: A vector of points in 3D space, represented as arrays `[x, y, z]`.
/// * `affine_weights`: A 4x4 array representing the affine transformation weights.
/// * `origin`: A 3D array `[x, y, z]` representing the origin of the transformation.
///
/// # Returns
/// * `Result<Vec<[f64; 3]>, String>`: A Result containing either:
/// A vector of transformed points (`Ok`)
/// An error message if something goes wrong (`Err`)
pub fn apply_affine_transform(
    points: &[[f64; 3]],
    affine_weights: &[[f64; 4]; 4],
) -> Result<Vec<[f64; 3]>, String> {
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
///
pub fn apply_bernstein_transform(
    points: &[[f64; 3]],
    deltas: &[[f64; 3]],
    resolution: &[usize; 3],
) -> Result<Vec<[f64; 3]>, String> {
    let dimension = [resolution[0] + 1, resolution[1] + 1, resolution[2] + 1];
    let mut transformed_points = vec![[0.0; 3]; points.len()];

    for point_idx in 0..points.len() {
        let point = points[point_idx];
        let mut aux_shift = [0.0; 3]; // This holds sum of shifts for x, y, z

        for i in 0..dimension[0] {
            // Compute bernstein for x dimension
            let aux1_x = (1.0 - point[0]).powi((dimension[0] - 1 - i) as i32);
            let aux2_x = point[0].powi(i as i32);
            let bernstein_x = binomial_coefficient(dimension[0] - 1, i) * aux1_x * aux2_x;

            for j in 0..dimension[1] {
                // Compute bernstein for y dimension
                let aux1_y = (1.0 - point[1]).powi((dimension[1] - 1 - j) as i32);
                let aux2_y = point[1].powi(j as i32);
                let bernstein_y =
                    binomial_coefficient(dimension[1] - 1, j) * aux1_y * aux2_y;

                for k in 0..dimension[2] {
                    // Compute bernstein for z dimension
                    let aux1_z = (1.0 - point[2]).powi((dimension[2] - 1 - k) as i32);
                    let aux2_z = point[2].powi(k as i32);
                    let bernstein_z =
                        binomial_coefficient(dimension[2] - 1, k) * aux1_z * aux2_z;

                    // Summation step: this is where we sum up the shifts for each dimension
                    let delta_id = i * dimension[1] * dimension[2] + j * dimension[2] + k;

                    for d in 0..3 {
                        aux_shift[d] +=
                            bernstein_x * bernstein_y * bernstein_z * deltas[delta_id][d];
                    }
                }
            }
        }

        // Add the shifts to the original point
        for d in 0..3 {
            transformed_points[point_idx][d] = point[d] + aux_shift[d];
        }
    }

    Ok(transformed_points)
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
