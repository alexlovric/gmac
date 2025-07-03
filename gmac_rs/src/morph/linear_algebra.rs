/// Solves a system of linear equations using least squares.
///
/// # Arguments
/// * `mat`: The design matrix.
/// * `rhs`: A vector of n-dimensional points as the right-hand side.
///
/// # Returns
/// * `Ok(Vec<[f64; 3]>)` if the system is successfully solved.
/// * `Err(String)` if the system cannot be solved, with an error message.
#[allow(clippy::needless_range_loop)]
pub fn least_squares_solver(
    mat: &[Vec<f64>],
    rhs: &[Vec<f64>],
) -> Result<Vec<Vec<f64>>, String> {
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

/// Solves a system of linear equations using LU decomposition.
///
/// # Arguments
/// * `mat`: The design matrix.
/// * `rhs`: A vector of n-dimensional points as the right-hand side.
///
/// # Returns
/// * `Ok(Vec<<Vec<f64>>>)` if the system is successfully solved.
/// * `Err(String)` if the system cannot be solved, with an error message.
pub fn lu_linear_solver(
    mat: &Vec<Vec<f64>>,
    rhs: &[Vec<f64>],
) -> Result<Vec<Vec<f64>>, String> {
    let mat_rows = mat.len();
    let rhs_dim = if let Some(first_row) = rhs.first() {
        first_row.len()
    } else {
        return Err(String::from("No rhs dimensions!"));
    };
    if mat_rows != rhs.len() {
        return Err(String::from(
            "Incompatible design matrix and right-hand side sizes!",
        ));
    }

    let (lu, p) = lu_decomposition(mat)?;

    // Initialize solution vectors for each dimension
    let mut x = vec![vec![0.0; rhs_dim]; mat_rows];

    for d in 0..rhs_dim {
        // Solve Ly = Pb for each dimension
        let mut y = vec![0.0; mat_rows];
        for i in 0..mat_rows {
            let mut sum = 0.0;
            for (j, yy) in y.iter().enumerate().take(i) {
                sum += lu[i][j] * yy;
            }
            y[i] = rhs[p[i]][d] - sum;
        }

        // Solve Ux = y for each dimension
        for i in (0..mat_rows).rev() {
            let mut sum = 0.0;
            for (j, xx) in x.iter().enumerate().take(mat_rows).skip(i + 1) {
                sum += lu[i][j] * xx[d];
            }
            x[i][d] = (y[i] - sum) / lu[i][i];
        }
    }

    Ok(x)
}

/// Performs LU decomposition with partial pivoting on a given square matrix.
/// The function decomposes a square matrix `mat` into its lower triangular
/// (L) and upper triangular (U) components, along with a permutation vector `p`
/// to handle row swaps due to partial pivoting.
///
/// # Parameters
/// * `mat: A reference to a square matrix.
///
/// # Returns
/// A tuple of `lu` and `p` as `Ok((lu, p))`
/// - `lu`` is a matrix where the diagonal and above contain elements of U,
///                       and below the diagonal contains elements of L.
/// - `p` is a permutation vector, representing the row swaps applied to `mat`.
///
/// Or an `Err(String)` if the decomposition fails, such as encountering a zero pivot.
///
/// # Example
/// ```
/// let mat = vec![
///     vec![1.0, 2.0, 3.0],
///     vec![4.0, 5.0, 6.0],
///     vec![7.0, 8.0, 9.0],
/// ];
/// let result = lu_decomposition(&mat);
/// ```
pub fn lu_decomposition(
    mat: &Vec<Vec<f64>>,
) -> Result<(Vec<Vec<f64>>, Vec<usize>), String> {
    let mut lu = mat.to_vec();
    let mut p = (0..mat.len()).collect::<Vec<usize>>();

    for k in 0..mat.len() - 1 {
        let mut max_row = k;
        for i in k + 1..mat.len() {
            if lu[i][k].abs() > lu[max_row][k].abs() {
                max_row = i;
            }
        }
        if lu[max_row][k] == 0.0 {
            return Err(String::from("Zero obtained in LU[max_row][k]!"));
        }
        if max_row != k {
            lu.swap(k, max_row);
            p.swap(k, max_row);
        }
        for i in k + 1..mat.len() {
            let factor = lu[i][k] / lu[k][k];
            lu[i][k] = factor;
            for j in k + 1..mat.len() {
                lu[i][j] -= factor * lu[k][j];
            }
        }
    }

    Ok((lu, p))
}
