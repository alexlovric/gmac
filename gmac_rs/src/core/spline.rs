use std::iter::repeat_with;

/// Generate a knot vector.
///
/// # Arguments
/// * `n`: Number of control points minus 1.
/// * `k`: Order of the spline.
///
/// # Returns
/// A vector of f64 representing the knot vector.
fn knot_vector(n: usize, k: usize) -> Vec<f64> {
    let mut t = vec![0.0; n + k + 1];
    for (ii, elem) in t.iter_mut().enumerate().skip(k).take(n + k + 1 - k) {
        *elem = if ii >= k && ii <= n {
            (ii as f64) - (k as f64) + 1.0
        } else {
            (n as f64) - (k as f64) + 2.0
        };
    }
    t
}
/// Calculate the B-spline basis function.
///
/// # Arguments
/// * `u`: The point where the basis function is evaluated.
/// * `t`: The knot vector.
/// * `i`: The index of the basis function.
/// * `k`: The order of the spline.
///
/// # Returns
/// A float representing the value of the basis function at `u`.
fn b_spline(u: f64, t: &[f64], i: usize, k: usize) -> f64 {
    if k == 1 {
        return if u >= t[i] && u < t[i + 1] { 1.0 } else { 0.0 };
    }
    let aa = t[i + k - 1] - t[i];
    let bb = t[i + k] - t[i + 1];
    let cc = if aa != 0.0 {
        (u - t[i]) * b_spline(u, t, i, k - 1) / aa
    } else {
        0.0
    };
    let dd = if bb != 0.0 {
        (t[i + k] - u) * b_spline(u, t, i + 1, k - 1) / bb
    } else {
        0.0
    };
    cc + dd
}

/// Calculate the B-spline basis function.
///
/// # Arguments
/// * `u`: The point where the basis function is evaluated.
/// * `t`: The knot vector.
/// * `i`: The index of the basis function.
/// * `k`: The order of the spline.
///
/// # Returns
/// A float representing the value of the basis function at `u`.
pub fn generate_bspline(
    order: usize,
    ngen: usize,
    mut cpts: Vec<Vec<f64>>,
    weights: Option<Vec<f64>>,
    knots: Option<Vec<f64>>,
    clamped: bool,
    closed: bool,
) -> Vec<Vec<f64>> {
    if closed {
        let initial_elements = cpts[0..order - 1].to_vec();
        cpts.extend(initial_elements);
    }

    let n = cpts.len() - 1;
    let weights = weights.unwrap_or_else(|| repeat_with(|| 1.0).take(n + 1).collect());
    let knots = match knots {
        Some(k) => k,
        None => {
            if clamped {
                knot_vector(n, order)
            } else {
                let mut v = vec![];
                for i in 0..(n + order + 1) {
                    v.push(i as f64 / (n + order) as f64);
                }
                v
            }
        }
    };

    let knot1 = knots[order - 1];
    let knot2 = knots[n + 1];
    let frac = (knot2 - knot1) / ngen as f64;
    let mut u = vec![];

    let mut val = knot1;
    while val < knot2 {
        u.push(val);
        val += frac;
    }

    let mut pts = vec![vec![0.0; cpts[0].len()]; ngen];

    for jj in 0..ngen {
        for ii in 0..(n + 1) {
            let bs = weights[ii] * b_spline(u[jj], &knots, ii, order);
            for kk in 0..cpts[0].len() {
                pts[jj][kk] += bs * cpts[ii][kk];
            }
        }
    }

    if closed {
        pts.push(pts[0].clone());
    }

    pts
}
