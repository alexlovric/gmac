/// Computes the squared Euclidean distance between two array of equal length.
///
/// # Arguments
/// * `a`: A floating point array.
/// * `b`: A second floating point array.
///
/// # Returns
/// The squared Euclidean distance.
pub fn squared_euclidean_distance_3d(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

/// Computes the Euclidean distance between two array of equal length.
///
/// # Arguments
/// * `a`: A floating point array.
/// * `b`: A second floating point array.
///
/// # Returns
/// The Euclidean distance.
pub fn euclidean_distance_3d(a: [f64; 3], b: [f64; 3]) -> f64 {
    squared_euclidean_distance_3d(a, b).sqrt()
}
