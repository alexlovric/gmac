use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Result, Seek, SeekFrom, Write};

/// Enumeration of supported STL file formats.
///
/// - `Ascii`: Human-readable ASCII STL format.
/// - `Binary`: Compact binary STL format, default choice for efficient storage and faster parsing.
pub enum StlFormat {
    Ascii,
    Binary,
}

/// Writes a mesh to STL file. Defaults to binary format.
///
/// # Arguments
/// * `nodes` - List of 3D vertex positions.
/// * `cells` - List of triangles (indices into nodes).
/// * `filename` - Optional file path. Defaults to `"mesh.stl"`.
/// * `format` - Optional STL format. Defaults to `StlFormat::Binary`.
///
/// # Returns
/// `Ok(())` on success, or an IO error.
pub fn write_stl(
    nodes: &[[f64; 3]],
    cells: &[[usize; 3]],
    filename: Option<&str>,
    format: Option<StlFormat>,
) -> Result<()> {
    match format.unwrap_or(StlFormat::Binary) {
        StlFormat::Ascii => write_stl_ascii(nodes, cells, filename),
        StlFormat::Binary => write_stl_binary(nodes, cells, filename),
    }
}

/// Writes a triangle mesh to an ASCII STL file.
///
/// The STL format represents a 3D surface mesh as a series of triangles.
/// Each triangle is described with a normal vector and three vertices.
///
/// # Arguments
/// * `nodes` - A list of 3D vertex positions.
/// * `cells` - A list of triangles, each defined by indices into the `nodes` array.
/// * `filename` - Optional file path. Defaults to `"mesh.stl"` if `None`.
///
/// # Returns
/// Returns `Ok(())` on success, or an `std::io::Error` on failure.
pub fn write_stl_ascii(
    nodes: &[[f64; 3]],
    cells: &[[usize; 3]],
    filename: Option<&str>,
) -> Result<()> {
    let mut file = File::create(filename.unwrap_or("mesh.stl"))?;
    writeln!(file, "solid exported_grid")?;

    for cell in cells {
        let (normal, [p1, p2, p3]) = compute_facet_data(nodes, cell);

        writeln!(
            file,
            "  facet normal {} {} {}",
            normal[0], normal[1], normal[2]
        )?;
        writeln!(file, "    outer loop")?;
        writeln!(file, "      vertex {} {} {}", p1[0], p1[1], p1[2])?;
        writeln!(file, "      vertex {} {} {}", p2[0], p2[1], p2[2])?;
        writeln!(file, "      vertex {} {} {}", p3[0], p3[1], p3[2])?;
        writeln!(file, "    endloop")?;
        writeln!(file, "  endfacet")?;
    }

    writeln!(file, "endsolid exported_grid")?;
    Ok(())
}

/// Writes a triangle mesh to a binary STL file.
///
/// Binary STL is a compact format used for 3D geometry exchange.
/// Each triangle includes a normal and three vertices stored as 32-bit floats.
/// An 80-byte header and a 4-byte triangle count are written before the triangle data.
///
/// # Arguments
/// * `nodes` - A list of 3D vertex positions.
/// * `cells` - A list of triangles, each defined by indices into the `nodes` array.
/// * `filename` - Optional file path. Defaults to `"mesh.stl"` if `None`.
///
/// # Returns
/// Returns `Ok(())` on success, or an `std::io::Error` on failure.
pub fn write_stl_binary(
    nodes: &[[f64; 3]],
    cells: &[[usize; 3]],
    filename: Option<&str>,
) -> Result<()> {
    // This closure defines the work to be done for a single triangle.
    let process_cell = |cell: &[usize; 3]| -> [u8; 50] {
        let (normal, [p1, p2, p3]) = compute_facet_data(nodes, cell);
        let mut buffer = [0u8; 50];
        let mut cursor = 0;

        // Write normal and vertices into the in-memory buffer
        for f in normal
            .iter()
            .chain(p1.iter())
            .chain(p2.iter())
            .chain(p3.iter())
        {
            let bytes = (*f as f32).to_le_bytes();
            buffer[cursor..cursor + 4].copy_from_slice(&bytes);
            cursor += 4;
        }
        buffer
    };

    // Compute all triangle byte data
    let all_triangle_data: Vec<[u8; 50]> = cells.iter().map(process_cell).collect();
    // #[cfg(feature = "rayon")]
    // let all_triangle_data: Vec<[u8; 50]> = cells.par_iter().map(process_cell).collect();

    // Write the final data to the file sequentially
    let mut file = File::create(filename.unwrap_or("mesh.stl"))?;

    // 80-byte header
    let mut header = [0u8; 80];
    header[..30].copy_from_slice(b"Binary STL generated with GMAC");
    file.write_all(&header)?;

    // Number of triangles
    let num_triangles = cells.len() as u32;
    file.write_all(&num_triangles.to_le_bytes())?;

    // Write all pre-computed triangle data
    let flat_buffer: Vec<u8> = all_triangle_data.into_iter().flatten().collect();
    file.write_all(&flat_buffer)?;

    Ok(())
}

/// Computes the facet normal and vertex coordinates for a triangle cell.
///
/// Given a triangle defined by three indices into the `nodes` array,
/// this function calculates the outward-facing normal using the right-hand rule,
/// and returns the normal along with the three vertex coordinates.
///
/// # Arguments
/// * `nodes` - A list of 3D vertex positions.
/// * `cell` - A triangle defined by indices into the `nodes` array (must have length 3).
///
/// # Returns
/// A tuple containing:
/// - `[f64; 3]`: The unit normal vector of the triangle.
/// - `[[f64; 3]; 3]`: The three vertex coordinates of the triangle.
fn compute_facet_data(
    nodes: &[[f64; 3]],
    cell: &[usize; 3],
) -> ([f64; 3], [[f64; 3]; 3]) {
    let p1 = nodes[cell[0]];
    let p2 = nodes[cell[1]];
    let p3 = nodes[cell[2]];

    // Vectors for normal
    let u = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
    let v = [p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]];

    let mut normal = [
        (u[1] * v[2]) - (u[2] * v[1]),
        (u[2] * v[0]) - (u[0] * v[2]),
        (u[0] * v[1]) - (u[1] * v[0]),
    ];

    // Normalize
    let norm =
        (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
    if norm > f64::EPSILON {
        normal = [normal[0] / norm, normal[1] / norm, normal[2] / norm];
    }

    (normal, [p1, p2, p3])
}

/// Reads an STL file (ASCII or binary) and extracts its mesh.
///
/// This function automatically detects the file format (ASCII or binary)
/// based on the initial bytes and calls the appropriate parser.
///
/// # Arguments
/// * `filename` - Path to the STL file.
///
/// # Returns
/// Returns a tuple `(nodes, cells)` where:
/// - `nodes`: Unique list of 3D points (`Vec<[f64; 3]>`)
/// - `cells`: Triangles as indices into `nodes` (`Vec<[usize; 3]>`)
///
/// # Errors
/// Returns an error if the file cannot be opened, read, or parsed.
pub fn read_stl(filename: &str) -> Result<(Vec<[f64; 3]>, Vec<[usize; 3]>)> {
    let file = File::open(filename)?;
    let mut reader = BufReader::new(file);

    // Peek first 512 bytes to detect ASCII or binary
    let mut peek_buf = [0u8; 512];
    let bytes_read = reader.read(&mut peek_buf)?;

    let header_str = std::str::from_utf8(&peek_buf[..bytes_read]).unwrap_or("");

    let is_ascii = header_str.trim_start().starts_with("solid")
        && header_str.contains("facet normal");

    // Reset reader to start for proper parsing
    reader.seek(SeekFrom::Start(0))?;

    if is_ascii {
        read_stl_ascii_from_buf(reader)
    } else {
        // For binary, skip the 80-byte header before reading triangles
        reader.seek(SeekFrom::Start(80))?;

        read_stl_binary_from_buf(reader)
    }
}

/// Parses an ASCII STL file and extracts the mesh data.
///
/// This parser reads each line looking for `vertex` entries,
/// deduplicates vertices using a hash map, and constructs a list of triangles.
///
/// # Arguments
/// * `filename` - Path to the ASCII STL file.
///
/// # Returns
/// Returns a tuple `(nodes, cells)` where:
/// - `nodes`: Unique list of 3D points (`Vec<[f64; 3]>`)
/// - `cells`: Triangles as indices into `nodes` (`Vec<[usize; 3]>`)
///
/// # Errors
/// Returns an error if the file cannot be opened or parsed correctly.
pub fn read_stl_ascii_from_buf<R: BufRead>(
    reader: R,
) -> Result<(Vec<[f64; 3]>, Vec<[usize; 3]>)> {
    let mut node_map: HashMap<(u64, u64, u64), usize> = HashMap::new();
    let mut nodes: Vec<[f64; 3]> = Vec::new();
    let mut cells: Vec<[usize; 3]> = Vec::new();

    let mut current_triangle = [0usize; 3];
    let mut vertex_index = 0;

    for line in reader.lines() {
        let line = line?.trim().to_string();
        if line.starts_with("vertex") {
            let parts: Vec<f64> = line
                .split_whitespace()
                .skip(1)
                .filter_map(|s| s.parse::<f64>().ok())
                .collect();
            if parts.len() == 3 {
                let point = [parts[0], parts[1], parts[2]];
                let key = (point[0].to_bits(), point[1].to_bits(), point[2].to_bits());
                let idx = *node_map.entry(key).or_insert_with(|| {
                    let new_idx = nodes.len();
                    nodes.push(point);
                    new_idx
                });
                current_triangle[vertex_index] = idx;
                vertex_index += 1;
                if vertex_index == 3 {
                    cells.push(current_triangle);
                    vertex_index = 0;
                }
            }
        }
    }

    Ok((nodes, cells))
}

/// Parses a binary STL file and extracts the mesh data.
///
/// This parser reads fixed-size 50-byte chunks for each triangle,
/// decodes 3 vertices and deduplicates them, building the mesh.
///
/// # Arguments
/// * `filename` - Path to the binary STL file.
///
/// # Returns
/// Returns a tuple `(nodes, cells)` where:
/// - `nodes`: Unique list of 3D points (`Vec<[f64; 3]>`)
/// - `cells`: Triangles as indices into `nodes` (`Vec<[usize; 3]>`)
///
/// # Errors
/// Returns an error if the file cannot be opened or parsed.
pub fn read_stl_binary_from_buf<R: Read>(
    mut reader: R,
) -> Result<(Vec<[f64; 3]>, Vec<[usize; 3]>)> {
    // Read number of triangles (4 bytes)
    let mut count_buf = [0u8; 4];
    reader.read_exact(&mut count_buf)?;
    let num_triangles = u32::from_le_bytes(count_buf) as usize;

    // Sanity check on number of triangles
    if num_triangles > 10_000_000 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("Too many triangles: {}", num_triangles),
        ));
    }

    let mut nodes: Vec<[f64; 3]> = Vec::with_capacity(num_triangles * 3);
    let mut cells: Vec<[usize; 3]> = Vec::with_capacity(num_triangles);
    let mut node_map: HashMap<(u64, u64, u64), usize> = HashMap::new();

    let mut tri_buf = [0u8; 50];

    for _ in 0..num_triangles {
        reader.read_exact(&mut tri_buf)?;

        let mut triangle = [0usize; 3];
        // Skip first 12 bytes (normal vector), then read vertices
        for (i, chunk) in tri_buf[12..48].chunks_exact(12).enumerate() {
            let x = f32::from_le_bytes(chunk[0..4].try_into().unwrap()) as f64;
            let y = f32::from_le_bytes(chunk[4..8].try_into().unwrap()) as f64;
            let z = f32::from_le_bytes(chunk[8..12].try_into().unwrap()) as f64;

            let key = (x.to_bits(), y.to_bits(), z.to_bits());
            let index = *node_map.entry(key).or_insert_with(|| {
                let idx = nodes.len();
                nodes.push([x, y, z]);
                idx
            });

            triangle[i] = index;
        }

        cells.push(triangle);
    }

    Ok((nodes, cells))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::Path;

    // Sample simple mesh: one triangle
    fn sample_mesh() -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
        (
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        )
    }

    #[test]
    fn test_write_stl_ascii_creates_file() {
        let (nodes, cells) = sample_mesh();
        let filename = "test_ascii.stl";
        let result = write_stl(&nodes, &cells, Some(filename), Some(StlFormat::Ascii));
        assert!(result.is_ok());
        assert!(Path::new(filename).exists());

        let content = fs::read_to_string(filename).expect("Failed to read ASCII STL");
        assert!(content.contains("solid"));
        assert!(content.contains("facet normal"));
        fs::remove_file(filename).unwrap();
    }

    #[test]
    fn test_write_stl_binary_creates_file() {
        let (nodes, cells) = sample_mesh();
        let filename = "test_binary.stl";
        let result = write_stl(&nodes, &cells, Some(filename), Some(StlFormat::Binary));
        assert!(result.is_ok());
        assert!(Path::new(filename).exists());

        let metadata = fs::metadata(filename).expect("Failed to get metadata");
        // Binary STL must be larger than minimal ASCII STL length (header + data)
        assert!(metadata.len() > 80);

        fs::remove_file(filename).unwrap();
    }

    #[test]
    fn test_write_stl_defaults_to_binary() {
        let (nodes, cells) = sample_mesh();
        let filename = "test_default.stl";
        let result = write_stl(&nodes, &cells, Some(filename), None);
        assert!(result.is_ok());
        assert!(Path::new(filename).exists());

        // Check header for binary STL: file size > 80 (header) + 4 (count) bytes
        let metadata = fs::metadata(filename).unwrap();
        assert!(metadata.len() > 84);

        fs::remove_file(filename).unwrap();
    }
}
