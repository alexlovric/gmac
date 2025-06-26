/// Character set for Base64 encoding.
const CHARSET: &[u8; 64] =
    b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

/// Padding character for Base64 encoding.
const PADDING: char = '=';

/// Encodes a slice of bytes into a Base64 string.
///
/// # Arguments
/// * `input`: A slice of bytes to encode.
///
/// # Returns
/// Returns a Base64-encoded `String`.
///
/// # Examples
/// ```
/// let text = "Hello world";
/// let encoded = base64_encode(text.as_bytes());
/// ```
pub fn encode(input: &[u8]) -> String {
    let mut output = String::new();
    let mut buffer: u32 = 0;
    let mut buffer_size = 0;
    for &byte in input {
        buffer = (buffer << 8) | byte as u32;
        buffer_size += 8;
        while buffer_size >= 6 {
            let index = (buffer >> (buffer_size - 6)) as usize & 0x3F;
            output.push(CHARSET[index] as char);
            buffer_size -= 6;
        }
    }

    if buffer_size > 0 {
        buffer <<= 6 - buffer_size;
        let index = buffer as usize & 0x3F;
        output.push(CHARSET[index] as char);
    }

    while output.len() % 4 != 0 {
        output.push(PADDING);
    }

    output
}

/// Decodes a Base64 string into a `Vec<u8>`.
///
/// # Arguments
/// * `input`: A Base64-encoded `&str`.
///
/// # Returns
/// Returns an `Option` containing a `Vec<u8>` if the decoding is successful,
/// or `None` if the input contains invalid characters.
///
/// # Examples
/// ```
/// if let Some(decoded_bytes) = base64_decode("SGVsbG8gd29ybGQ=") {
///     let decoded_string = String::from_utf8(decoded_bytes).unwrap();
/// }
/// ```
pub fn decode(input: &str) -> Option<Vec<u8>> {
    let mut output = Vec::new();
    let mut buffer = 0u32;
    let mut buffer_size = 0usize;

    for c in input.chars() {
        if c == PADDING {
            break;
        }

        if let Some(pos) = CHARSET.iter().position(|&ch| ch == c as u8) {
            buffer = (buffer << 6) | pos as u32;
            buffer_size += 6;

            if buffer_size >= 8 {
                let decoded_byte = ((buffer >> (buffer_size - 8)) & 0xFF) as u8;
                output.push(decoded_byte);
                buffer_size -= 8;
            }
        } else {
            return None; // Invalid character in input
        }
    }

    Some(output)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tests base64 encoding and decoding for an example text.
    #[test]
    fn test_base64_encode_decode() {
        let text = "Hello, world!";
        let encoded = encode(text.as_bytes());
        let decoded = decode(&encoded).expect("Failed to decode");

        assert_eq!(text.as_bytes(), decoded.as_slice());
    }

    /// Tests base64 decoding for an invalid input.
    #[test]
    fn test_base64_decode_invalid_input() {
        let encoded = "SGVs@sbG8="; // Invalid character '@'
        let decoded = decode(encoded);

        assert!(decoded.is_none());
    }
}
