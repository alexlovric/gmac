use hologram::linear_algebra::lu_linear_solver;

use crate::error::Result;

/// Solves a system of linear equations using least squares.
///
/// # Arguments
/// * `mat`: The design matrix.
/// * `rhs`: A vector of n-dimensional points as the right-hand side.
///
/// # Returns
/// * `Ok(Vec<Vec<f64>>)`: A Result containing either:
/// A vector of transformed points (`Ok`)
/// An error if something goes wrong (`Err`)
pub fn least_squares_solver(mat: &[Vec<f64>], rhs: &[Vec<f64>]) -> Result<Vec<Vec<f64>>> {
    let mat_rows = mat.len();
    let mat_cols = mat[0].len();
    let rhs_dim = rhs[0].len();

    // Calculate A^T A
    let mut ata = vec![vec![0.0; mat_cols]; mat_cols];
    for i in 0..mat_cols {
        for j in 0..mat_cols {
            for k in 0..mat_rows {
                ata[i][j] += mat[k][i] * mat[k][j];
            }
        }
    }

    // Calculate A^T b for each dimension
    let mut atb = vec![vec![0.0; rhs_dim]; mat_cols];
    for i in 0..mat_cols {
        for d in 0..rhs_dim {
            for k in 0..mat_rows {
                atb[i][d] += mat[k][i] * rhs[k][d];
            }
        }
    }

    let x = lu_linear_solver(&ata, &atb)?;

    Ok(x)
}
