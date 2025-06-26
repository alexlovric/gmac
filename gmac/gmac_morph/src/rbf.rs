use ndarray::{Array2, Axis};
use ndarray_linalg::Solve;

pub struct RbfDeformer {
    /// Input control points.
    pub x: Vec<[f64; 3]>,
    /// Kernel function.
    pub kernel: fn(f64, f64) -> f64,
    /// Kernel bandwidth parameter.
    pub epsilon: f64,
    /// Interpolation coefficients.
    weights: Vec<[f64; 3]>,
    /// Columns removed from the design matrix, i.e. unchanging
    pub removed_columns: Vec<usize>,
    /// Standard deviation of the original control points.
    pub x_std: [f64; 3],
    /// Mean of the deformed control points.
    pub y_mean: [f64; 3],
    /// Standard deviation of the deformed control points.
    pub y_std: [f64; 3],
}

impl RbfDeformer {
    /// Instantiates a new `RbfDeformer` instance.
    ///
    /// # Arguments
    /// * `original_control_points`: A n*3 array containing the original node positions.
    /// * `deformed_control_points`: A n*3 array containing the deformed node positions.
    /// * `kernel`: An optional function that computes the kernel function value.
    ///             Will default to Gaussian kernel if `None` given.
    /// * `epsilon`: An optional bandwidth parameter for the kernel.
    ///              Defaults to 1. if `None` given.
    ///
    /// # Returns
    /// A new `RbfDeformer` instance.
    pub fn new(
        original_control_points: Vec<[f64; 3]>,
        deformed_control_points: Vec<[f64; 3]>,
        kernel: Option<fn(f64, f64) -> f64>,
        epsilon: Option<f64>,
    ) -> Result<Self, String> {
        let x = original_control_points;
        let y = deformed_control_points;
        assert_eq!(x.len(), y.len(), "x and y must have the same length");
        let n = x.len();

        // Set kernel and epsilon (default to Gaussian and 1.0)
        let kernel = kernel.unwrap_or(gaussian_kernel);
        let epsilon = epsilon.unwrap_or(1.0);

        // Compute standard deviations for x (per dimension)
        let mut x_std = [0.0; 3];
        for d in 0..3 {
            let mean = x.iter().map(|p| p[d]).sum::<f64>() / n as f64;
            let variance =
                x.iter().map(|p| (p[d] - mean).powi(2)).sum::<f64>() / n as f64;
            let std = variance.sqrt();
            x_std[d] = if std < 1e-8 { 1.0 } else { std };
        }

        // Compute means and stds for y, track removed columns
        let mut y_mean = [0.0; 3];
        let mut y_std = [1.0; 3];
        let mut removed_columns = Vec::new();
        let mut normalised_y = Array2::zeros((n, 3));

        for d in 0..3 {
            let mean = y.iter().map(|p| p[d]).sum::<f64>() / n as f64;
            let variance =
                y.iter().map(|p| (p[d] - mean).powi(2)).sum::<f64>() / n as f64;
            let std = variance.sqrt();

            y_mean[d] = mean;
            if std < 1e-8 {
                removed_columns.push(d);
            } else {
                y_std[d] = std;
                for i in 0..n {
                    normalised_y[(i, d)] = (y[i][d] - mean) / std;
                }
            }
        }

        // Build design matrix (n x n) using Euclidean distances
        let mut design_matrix = Array2::zeros((n, n));
        for i in 0..n {
            for j in 0..n {
                let mut dist_sq = 0.0;
                for d in 0..3 {
                    let diff = (x[i][d] - x[j][d]) / x_std[d];
                    dist_sq += diff * diff;
                }
                design_matrix[(i, j)] = kernel(dist_sq.max(f64::EPSILON), epsilon);
            }
        }

        // Filter y columns and solve for weights
        let kept_columns: Vec<usize> =
            (0..3).filter(|d| !removed_columns.contains(d)).collect();
        let filtered_y = normalised_y.select(Axis(1), &kept_columns);
        let mut weights_matrix = Array2::zeros((n, kept_columns.len()));

        for (col_idx, col) in filtered_y.columns().into_iter().enumerate() {
            let w = design_matrix
                .solve(&col)
                .map_err(|e| format!("Linear solve failed: {}", e))?;
            weights_matrix.column_mut(col_idx).assign(&w);
        }

        // Reconstruct full weights (3D) with zeros for removed columns
        let mut weights = vec![[0.0; 3]; n];
        for (col_idx, &d) in kept_columns.iter().enumerate() {
            for i in 0..n {
                weights[i][d] = weights_matrix[(i, col_idx)];
            }
        }

        Ok(Self {
            x,
            kernel,
            epsilon,
            weights,
            removed_columns,
            x_std,
            y_mean,
            y_std,
        })
    }

    pub fn deform(&self, points: &[[f64; 3]]) -> Result<Vec<[f64; 3]>, String> {
        let mut result = vec![[0.0; 3]; points.len()];

        for (idx, point) in points.iter().enumerate() {
            let mut sum_contrib = [0.0; 3];

            // For each control point in training set
            for (i, train_point) in self.x.iter().enumerate() {
                // Compute squared Euclidean distance (normalized by x_std)
                let mut dist_sq = 0.0;
                for d in 0..3 {
                    let diff = (point[d] - train_point[d]) / self.x_std[d];
                    dist_sq += diff * diff;
                }
                let kernel_val = (self.kernel)(dist_sq.max(f64::EPSILON), self.epsilon);

                // Accumulate weighted kernel contributions
                for d in 0..3 {
                    if !self.removed_columns.contains(&d) {
                        sum_contrib[d] += self.weights[i][d] * kernel_val;
                    }
                }
            }

            // De-normalize results
            for d in 0..3 {
                result[idx][d] = if !self.removed_columns.contains(&d) {
                    sum_contrib[d] * self.y_std[d] + self.y_mean[d]
                } else {
                    self.y_mean[d]
                };
            }

            // Check for NaNs
            if result[idx].iter().any(|&v| v.is_nan()) {
                return Err("NaN values in output".to_string());
            }
        }

        Ok(result)
    }
}

/// kernels (All kernels accept squared Euclidean for the `distance`)
/// Guassian
pub fn gaussian_kernel(distance: f64, bandwidth: f64) -> f64 {
    (-0.5 * distance / bandwidth.powi(2)).exp()
}

/// Multiquadratic
pub fn multiquadric_kernel(distance: f64, bandwidth: f64) -> f64 {
    (distance + bandwidth.powi(2)).sqrt()
}

/// Inverse
pub fn inverse_multi_kernel(distance: f64, bandwidth: f64) -> f64 {
    1.0 / multiquadric_kernel(distance, bandwidth)
}

pub fn thin_plate_spline_kernel(distance: f64, _bandwidth: f64) -> f64 {
    if distance == 0.0 {
        0.0
    } else {
        distance.powi(2) * distance.ln()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_gaussian_kernel() {
        // Test known values
        assert_eq!(gaussian_kernel(0.0, 1.0), 1.0);
        assert_relative_eq!(gaussian_kernel(1.0, 1.0), (-0.5f64).exp(), epsilon = 1e-10);
    }

    #[test]
    fn test_thin_plate_spline_kernel() {
        assert_eq!(thin_plate_spline_kernel(0.0, 1.0), 0.0);
        assert_relative_eq!(thin_plate_spline_kernel(1.0, 1.0), 0.0, epsilon = 1e-10);
        assert_relative_eq!(thin_plate_spline_kernel(2.0, 1.0), 2.0_f64.ln() * 4.0, epsilon = 1e-10);
    }

    #[test]
    fn test_single_point() {
        let rbf = RbfDeformer::new(
            vec![[1.0, 2.0, 3.0]],
            vec![[2.0, 3.0, 4.0]],
            None,
            None,
        ).unwrap();
        
        // Should return exact deformation for training points
        let result = rbf.deform(&[[1.0, 2.0, 3.0]]).unwrap();
        assert_eq!(result[0], [2.0, 3.0, 4.0]);
    }

    #[test]
    fn test_identity_deformation() {
        let points = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let rbf = RbfDeformer::new(
            points.clone(),
            points.clone(),
            None,
            None,
        ).unwrap();
        
        // Should return exact same points
        let result = rbf.deform(&points).unwrap();
        for (res, pt) in result.iter().zip(points.iter()) {
            assert_relative_eq!(res[0], pt[0], epsilon = 1e-10);
            assert_relative_eq!(res[1], pt[1], epsilon = 1e-10);
            assert_relative_eq!(res[2], pt[2], epsilon = 1e-10);
        }
    }

    #[test]
    fn test_constant_deformation() {
        let original = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let deformed = vec![[10.0, 10.0, 10.0], [10.0, 10.0, 10.0]];
        let rbf = RbfDeformer::new(
            original,
            deformed,
            None,
            None,
        ).unwrap();
        
        // All points should map to [10.0, 10.0, 10.0]
        let result = rbf.deform(&[[2.0, 3.0, 4.0], [5.0, 6.0, 7.0]]).unwrap();
        assert_eq!(result, vec![[10.0, 10.0, 10.0], [10.0, 10.0, 10.0]]);
    }

    #[test]
    fn test_deform_standard() {
        let rbf = RbfDeformer::new(
            vec![[1.0, 2.0, 1.0], [3.0, 4.0, 2.0]],
            vec![[2.0, 3.0, 2.0], [4.0, 5.0, 3.0]],
            None,
            None,
        ).unwrap();

        let x_new = vec![[1.5, 2.6, 1.8]];
        let prediction = rbf.deform(&x_new).unwrap();

        // Compare the predicted result with the expected result
        assert_relative_eq!(
            prediction[0][0], 2.9073001606088247, epsilon = 1e-10
        );
        assert_relative_eq!(
            prediction[0][1], 3.9073001606088247, epsilon = 1e-10
        );
        assert_relative_eq!(
            prediction[0][2], 2.4536500803044126, epsilon = 1e-10
        );
    }

    #[test]
    fn test_different_kernels() {
        let points = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        
        // Test with each kernel type
        for kernel in &[gaussian_kernel, multiquadric_kernel, inverse_multi_kernel, thin_plate_spline_kernel] {
            let rbf = RbfDeformer::new(
                points.clone(),
                points.clone(),
                Some(*kernel),
                None,
            ).unwrap();
            
            let result = rbf.deform(&points).unwrap();
            for (res, pt) in result.iter().zip(points.iter()) {
                assert_relative_eq!(res[0], pt[0], epsilon = 1e-10);
                assert_relative_eq!(res[1], pt[1], epsilon = 1e-10);
                assert_relative_eq!(res[2], pt[2], epsilon = 1e-10);
            }
        }
    }

    #[test]
    #[should_panic(expected = "x and y must have the same length")]
    fn test_mismatched_lengths() {
        RbfDeformer::new(
            vec![[1.0, 2.0, 3.0]],
            vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
            None,
            None,
        ).unwrap();
    }
}
