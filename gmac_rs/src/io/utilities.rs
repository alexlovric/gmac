/// Converts a 32-bit floating-point number to its little-endian byte representation.
///
/// # Arguments
/// * `f`: A 32-bit floating-point number.
///
/// # Returns
/// Returns an array of 4 bytes representing the little-endian byte order of the input floating-point number.
///
/// # Examples
/// ```
/// let bytes = f32_to_bytes(42.0);
/// ```
pub fn f32_to_bytes(f: f32) -> [u8; 4] {
    f.to_le_bytes()
}

/// Converts a 32-bit integer to its little-endian byte representation.
///
/// # Arguments
/// * `i`: A 32-bit integer.
///
/// # Returns
/// Returns an array of 4 bytes representing the little-endian byte order of the input integer.
///
/// # Examples
/// ```
/// let bytes = i32_to_bytes(42);
/// ```
pub fn i32_to_bytes(i: i32) -> [u8; 4] {
    i.to_le_bytes()
}
