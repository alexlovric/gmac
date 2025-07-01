use hologram::{
    kernels::{
        cubic_kernel, gaussian_kernel, inverse_multi_kernel, linear_kernel,
        multiquadric_kernel, thin_plate_spline_kernel,
    },
    rbf::Rbf,
    Interpolator,
};

/// A Radial Basis Function (RBF) deformer for 3D point transformations.
///
/// This struct implements a deformable model that can smoothly interpolate between
/// a set of control points in 3D space. It's particularly useful for mesh deformation,
/// shape morphing, and other spatial transformations.
///
/// # Fields
/// - `x_mean`: Mean of the input control points for normalization
/// - `x_std`: Standard deviation of the input control points for normalization
/// - `y_mean`: Mean of the output control points for denormalization
/// - `y_std`: Standard deviation of the output control points for denormalization
/// - `removed_columns`: Indices of dimensions with zero variance in the output
/// - `rbf`: The underlying RBF interpolator
pub struct RbfDeformer {
    x_mean: [f64; 3],
    x_std: [f64; 3],
    y_mean: [f64; 3],
    y_std: [f64; 3],
    removed_columns: Vec<usize>,
    rbf: Rbf<[f64; 3], [f64; 3]>,
}

impl RbfDeformer {
    /// Creates a new RbfDeformer instance.
    ///
    /// # Arguments
    /// * `x` - Input control points (n×3 array)
    /// * `y` - Corresponding output control points (n×3 array)
    /// * `kernel_name` - Name of the kernel function to use (optional, defaults to "gaussian"):
    ///   - "linear": Linear kernel
    ///   - "cubic": Cubic kernel
    ///   - "gaussian": Gaussian kernel (default)
    ///   - "multiquadric": Multiquadric kernel
    ///   - "inverse_multiquadratic": Inverse multiquadric kernel
    ///   - "thin_plate_spline": Thin plate spline kernel
    /// * `epsilon` - Bandwidth parameter for the kernel (optional, defaults to 1.0)
    ///
    /// # Returns
    /// A new `RbfDeformer` instance or an error string if creation fails.
    ///
    /// # Panics
    /// Panics if `x` and `y` have different lengths.
    ///
    /// # Example
    /// ```
    /// use gmac::morph::rbf::RbfDeformer;
    ///
    /// let x = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]];
    /// let y = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
    /// let deformer = RbfDeformer::new(x, y, Some("gaussian"), Some(1.0)).unwrap();
    /// ```
    pub fn new(
        x: Vec<[f64; 3]>,
        y: Vec<[f64; 3]>,
        kernel_name: Option<&str>,
        epsilon: Option<f64>,
    ) -> Result<Self, String> {
        assert_eq!(x.len(), y.len(), "x and y must have the same length");

        let epsilon = epsilon.unwrap_or(1.0);
        let kernel: fn(f64, f64) -> f64 = match kernel_name.unwrap_or("gaussian") {
            "linear" => linear_kernel,
            "cubic" => cubic_kernel,
            "gaussian" => gaussian_kernel,
            "multiquadric" => multiquadric_kernel,
            "inverse_multiquadratic" => inverse_multi_kernel,
            "thin_plate_spline" => thin_plate_spline_kernel,
            other => return Err(format!("Unsupported kernel: {other}")),
        };

        let n = x.len();

        // Compute x mean and std
        let mut x_mean = [0.0; 3];
        let mut x_std = [1.0; 3];
        for d in 0..3 {
            let mean = x.iter().map(|p| p[d]).sum::<f64>() / n as f64;
            let std =
                (x.iter().map(|p| (p[d] - mean).powi(2)).sum::<f64>() / n as f64).sqrt();
            x_mean[d] = mean;
            x_std[d] = if std < 1e-8 { 1.0 } else { std };
        }

        let normalized_x: Vec<[f64; 3]> = x
            .iter()
            .map(|p| {
                let mut np = [0.0; 3];
                for d in 0..3 {
                    np[d] = (p[d] - x_mean[d]) / x_std[d];
                }
                np
            })
            .collect();

        // Normalize y and detect constant columns
        let mut y_mean = [0.0; 3];
        let mut y_std = [1.0; 3];
        let mut removed_columns = Vec::new();

        for d in 0..3 {
            let mean = y.iter().map(|p| p[d]).sum::<f64>() / n as f64;
            let std =
                (y.iter().map(|p| (p[d] - mean).powi(2)).sum::<f64>() / n as f64).sqrt();
            y_mean[d] = mean;
            if std < 1e-8 {
                removed_columns.push(d);
            } else {
                y_std[d] = std;
            }
        }

        let normalized_y: Vec<[f64; 3]> = y
            .iter()
            .map(|p| {
                let mut np = [0.0; 3];
                for d in 0..3 {
                    if !removed_columns.contains(&d) {
                        np[d] = (p[d] - y_mean[d]) / y_std[d];
                    }
                }
                np
            })
            .collect();

        let rbf = Rbf::new(normalized_x, normalized_y, Some(kernel), Some(epsilon))
            .map_err(|e| format!("Failed to create RBF: {e}"))?;

        Ok(Self {
            x_mean,
            x_std,
            y_mean,
            y_std,
            removed_columns,
            rbf,
        })
    }

    /// Deforms input points using the learned RBF transformation.
    ///
    /// # Arguments
    /// * `points` - A slice of 3D points to transform
    ///
    /// # Returns
    /// A `Vec` of transformed points with the same length as the input, or an error string
    /// if the transformation fails.
    ///
    /// # Example
    /// ```
    /// # use gmac::morph::rbf::RbfDeformer;
    /// # let x = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]];
    /// # let y = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
    /// # let deformer = RbfDeformer::new(x, y, Some("gaussian"), Some(1.0)).unwrap();
    /// let points = [[0.5, 0.5, 0.5], [0.2, 0.8, 0.4]];
    /// let deformed = deformer.deform(&points).unwrap();
    /// assert_eq!(deformed.len(), 2);
    /// ```
    pub fn deform(&self, points: &[[f64; 3]]) -> Result<Vec<[f64; 3]>, String> {
        let normalized_input: Vec<[f64; 3]> = points
            .iter()
            .map(|p| {
                let mut np = [0.0; 3];
                for d in 0..3 {
                    np[d] = (p[d] - self.x_mean[d]) / self.x_std[d];
                }
                np
            })
            .collect();

        let normalized_output = self
            .rbf
            .predict(&normalized_input)
            .map_err(|e| format!("Prediction failed: {e}"))?;

        let mut result = vec![[0.0; 3]; points.len()];
        for (i, p) in normalized_output.iter().enumerate() {
            for d in 0..3 {
                result[i][d] = if self.removed_columns.contains(&d) {
                    self.y_mean[d]
                } else {
                    p[d] * self.y_std[d] + self.y_mean[d]
                };
            }
        }

        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_single_point() {
        let rbf =
            RbfDeformer::new(vec![[1.0, 2.0, 3.0]], vec![[2.0, 3.0, 4.0]], None, None)
                .unwrap();

        // Should return exact deformation for training points
        let result = rbf.deform(&[[1.0, 2.0, 3.0]]).unwrap();
        assert_eq!(result[0], [2.0, 3.0, 4.0]);
    }

    #[test]
    fn test_constant_deformation() {
        let original = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let deformed = vec![[10.0, 10.0, 10.0], [10.0, 10.0, 10.0]];
        let rbf = RbfDeformer::new(original, deformed, None, None).unwrap();

        // All points should map to [10.0, 10.0, 10.0]
        let result = rbf.deform(&[[2.0, 3.0, 4.0], [5.0, 6.0, 7.0]]).unwrap();
        assert_eq!(result, vec![[10.0, 10.0, 10.0], [10.0, 10.0, 10.0]]);
    }

    #[test]
    fn test_identity_deformation() {
        let points = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let rbf = RbfDeformer::new(points.clone(), points.clone(), None, None).unwrap();

        // Should return exact same points
        let result = rbf.deform(&points).unwrap();
        for (res, pt) in result.iter().zip(points.iter()) {
            assert_relative_eq!(res[0], pt[0], epsilon = 1e-10);
            assert_relative_eq!(res[1], pt[1], epsilon = 1e-10);
            assert_relative_eq!(res[2], pt[2], epsilon = 1e-10);
        }
    }

    #[test]
    fn test_deform_standard() {
        let rbf = RbfDeformer::new(
            vec![[1.0, 2.0, 1.0], [3.0, 4.0, 2.0]],
            vec![[2.0, 3.0, 2.0], [4.0, 5.0, 3.0]],
            None,
            None,
        )
        .unwrap();

        let x_new = vec![[1.5, 2.6, 1.8]];
        let prediction = rbf.deform(&x_new).unwrap();

        // Compare the predicted result with the expected result
        assert_relative_eq!(prediction[0][0], 2.9073001606088247, epsilon = 1e-10);
        assert_relative_eq!(prediction[0][1], 3.9073001606088247, epsilon = 1e-10);
        assert_relative_eq!(prediction[0][2], 2.4536500803044126, epsilon = 1e-10);
    }

    #[test]
    fn test_different_kernels() {
        let points = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];

        // Test with each kernel type
        for kernel in &[
            "gaussian",
            "multiquadric",
            "inverse_multiquadratic",
            "thin_plate_spline",
        ] {
            let rbf =
                RbfDeformer::new(points.clone(), points.clone(), Some(*kernel), None)
                    .unwrap();

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
        )
        .unwrap();
    }
}
