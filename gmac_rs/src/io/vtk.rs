use std::fs::File;
use std::io::{Write, BufReader, BufRead};

use crate::io::utilities::{f32_to_bytes, i32_to_bytes};
use crate::io::base64::encode;
use crate::error::{Result};

/// Reads a 3D mesh from a VTU (VTK UnstructuredGrid) file.
///
/// # Arguments
/// * `filename`: The path of the file to read from.
///
/// # Returns
/// Returns a `Result` which contains a tuple (`Vec<[f64; 3]>`, `Vec<[usize; 3]>`)
/// if the file is successfully read, or contains an error otherwise.
#[allow(clippy::type_complexity)]
pub fn read_vtu(filename: &str) -> Result<(Vec<[f64; 3]>, Vec<[usize; 3]>)> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let mut nodes = Vec::new();
    let mut cells = Vec::new();

    let mut inside_points_data_array = false;
    let mut inside_cells_data_array = false;

    for line in reader.lines() {
        let line = line?;
        if line.contains("<DataArray") && line.contains("Name=\"Points\"") {
            inside_points_data_array = true;
        } else if inside_points_data_array && line.contains("</DataArray>") {
            inside_points_data_array = false;
        } else if line.contains("<DataArray") && line.contains("Name=\"connectivity\"") {
            inside_cells_data_array = true;
        } else if inside_cells_data_array && line.contains("</DataArray>") {
            inside_cells_data_array = false;
        } else if inside_points_data_array {
            let coords: Vec<f64> = line
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
            if coords.len() % 3 == 0 {
                for i in (0..coords.len()).step_by(3) {
                    nodes.push([coords[i], coords[i + 1], coords[i + 2]]);
                }
            }
        } else if inside_cells_data_array {
            let indices: Vec<usize> = line
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
            if indices.len() % 3 == 0 {
                for i in (0..indices.len()).step_by(3) {
                    cells.push([indices[i], indices[i + 1], indices[i + 2]]);
                }
            }
        }
    }

    Ok((nodes, cells))
}

/// Reads a 3D mesh from a VTP (VTK PolyData) file.
///
/// # Arguments
/// * `filename`: The path of the file to read from.
///
/// # Returns
/// Returns a `Result` which contains a `Vec<[f64; 3]>` if the file is successfully read,
/// or contains an error otherwise.
pub fn read_vtp(filename: &str) -> Result<Vec<[f64; 3]>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let mut nodes = Vec::new();
    let mut inside_data_array = false;

    for line in reader.lines() {
        let line = line?;
        if line.contains("<DataArray") && line.contains("Name=\"Points\"") {
            inside_data_array = true;
        } else if inside_data_array && line.contains("</DataArray>") {
            inside_data_array = false;
        } else if inside_data_array {
            let coords: Vec<f64> = line
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();

            if coords.len() % 3 == 0 {
                for i in (0..coords.len()).step_by(3) {
                    nodes.push([coords[i], coords[i + 1], coords[i + 2]]);
                }
            }
        }
    }

    Ok(nodes)
}

/// Writes the given 3D mesh to a VTU (VTK Unstructured Grid) file.
///
/// # Arguments
/// * `nodes`: A reference to a vector of coordinates.
/// * `cells`: A reference to a vector of cells.
/// * `filename`: An `Option` containing the path of the file to write to.
///
/// # Returns
/// Returns a `Result` which is `Ok` if the file is successfully written,
/// or contains an error otherwise.
///
/// # Example
/// ```ignore
/// let nodes = vec![vec![0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0], vec![0.0, 1.0, 0.0]];
/// let cells = vec![vec![0, 1, 2]];
/// write_vtu(&nodes, &cells, None).expect("Failed to write VTU file");
/// ```
pub fn write_vtu(
    nodes: &Vec<[f64; 3]>,
    cells: &Vec<[usize; 3]>,
    filename: Option<&str>,
) -> Result<()> {
    let mut file = File::create(filename.unwrap_or("mesh.vtu"))
        .map_err(|e| e.to_string())
        .unwrap();

    writeln!(
        file,
        r#"<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
  <UnstructuredGrid>"#
    )?;
    writeln!(
        file,
        "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">",
        nodes.len(),
        cells.len()
    )?;

    write_nodes(&mut file, nodes)?;
    write_cells(&mut file, cells)?;

    writeln!(file, "    </Piece>")?;
    writeln!(file, "  </UnstructuredGrid>")?;
    writeln!(file, "</VTKFile>")?;

    Ok(())
}

/// Writes the given 3D mesh to a VTP (VTK PolyData) file.
///
/// # Arguments
/// * `nodes`: A reference to a 3D coordinate nodes vector.
/// * `filename`: An `Option` containing the path of the file to write to.
///
/// # Returns
/// Returns a `Result` which is `Ok` if the file is successfully written,
/// or contains an error otherwise.
pub fn write_vtp(nodes: &[[f64; 3]], filename: Option<&str>) -> Result<()> {
    let mut file = File::create(filename.unwrap_or("grid.vtp"))
        .map_err(|e| e.to_string())
        .unwrap();

    writeln!(
        file,
        r#"<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">
  <PolyData>"#
    )?;

    writeln!(
        file,
        "    <Piece NumberOfPoints=\"{}\"  NumberOfVerts=\"{}\">",
        nodes.len(),
        nodes.len()
    )?;
    write_nodes(&mut file, nodes)?;

    writeln!(file, "      <Verts>")?;
    writeln!(
        file,
        "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
    )?;
    write!(file, "          ")?;
    for i in 0..nodes.len() {
        write!(file, "{} ", i)?; // VTK_TRIANGLE is 5
    }
    writeln!(file, "\n        </DataArray>")?;
    writeln!(
        file,
        "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
    )?;
    write!(file, "          ")?;
    for i in 0..nodes.len() {
        write!(file, "{} ", i + 1)?; // VTK_TRIANGLE is 5
    }
    writeln!(file, "\n        </DataArray>")?;
    writeln!(file, "      </Verts>")?;
    writeln!(file, "    </Piece>")?;
    writeln!(file, "  </PolyData>")?;
    writeln!(file, "</VTKFile>")?;

    Ok(())
}

pub fn write_vtp_binary(nodes: &Vec<[f64; 3]>, filename: Option<&str>) -> Result<()> {
    let file_name = filename.unwrap_or("output.vtp");
    let mut file = File::create(file_name)?;

    // Write VTP header
    writeln!(
        file,
        r#"<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">
  <PolyData>
    <Piece NumberOfPoints="{}" NumberOfVerts="{}">"#,
        nodes.len(),
        nodes.len()
    )?;

    // Write Points and Verts DataArray reference
    writeln!(file, "      <Points>")?;
    writeln!(file, "        <DataArray type='Float32' NumberOfComponents='3' format='appended' offset='0' />")?;
    writeln!(file, "      </Points>")?;

    writeln!(file, "      <Verts>")?;
    writeln!(file, "        <DataArray type='Int32' Name='connectivity' format='appended' offset='{}' />", 4 + 12 * nodes.len())?;
    writeln!(
        file,
        "        <DataArray type='Int32' Name='offsets' format='appended' offset='{}' />",
        4 + 12 * nodes.len() + 4 + 4 * nodes.len()
    )?;
    writeln!(file, "      </Verts>")?;

    writeln!(file, "    </Piece>\n  </PolyData>")?;
    writeln!(file, "  <AppendedData encoding=\"base64\">")?;
    write!(file, "_")?;

    // Prepare binary data for Points
    let mut points_data: Vec<u8> = Vec::new();
    points_data.extend_from_slice(&i32_to_bytes((12 * nodes.len()) as i32));
    for &[x, y, z] in nodes {
        points_data.extend_from_slice(&f32_to_bytes(x as f32));
        points_data.extend_from_slice(&f32_to_bytes(y as f32));
        points_data.extend_from_slice(&f32_to_bytes(z as f32));
    }

    // Prepare binary data for Connectivity and Offsets
    let mut connectivity_data: Vec<u8> = Vec::new();
    connectivity_data.extend_from_slice(&i32_to_bytes((4 * nodes.len()) as i32));
    for i in 0..nodes.len() {
        connectivity_data.extend_from_slice(&i32_to_bytes(i as i32));
    }

    let mut offsets_data: Vec<u8> = Vec::new();
    offsets_data.extend_from_slice(&i32_to_bytes((4 * nodes.len()) as i32));
    for i in 1..=nodes.len() {
        offsets_data.extend_from_slice(&i32_to_bytes(i as i32));
    }

    // Combine all binary data and base64 encode
    let mut all_data = Vec::new();
    all_data.extend_from_slice(&points_data);
    all_data.extend_from_slice(&connectivity_data);
    all_data.extend_from_slice(&offsets_data);
    let base64_encoded = encode(&all_data);

    // Write base64 encoded data to file
    writeln!(file, "{}", base64_encoded)?;

    // Close the XML tags
    writeln!(file, "  </AppendedData>\n</VTKFile>")?;

    Ok(())
}

/// Writes nodes to the current VTK file.
///
/// # Arguments
/// * `file`: Std fs file.
///
/// # Returns
/// Returns a `Result` which is `Ok` if the file is successfully written,
/// or contains an error otherwise.
fn write_nodes(file: &mut File, nodes: &[[f64; 3]]) -> std::io::Result<()> {
    writeln!(file, "      <Points>")?;
    writeln!(
        file,
        "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">"
    )?;
    for point in nodes {
        writeln!(file, "          {} {} {}", point[0], point[1], point[2])?;
    }
    writeln!(file, "        </DataArray>")?;
    writeln!(file, "      </Points>")
}

fn write_cells(file: &mut File, cells: &Vec<[usize; 3]>) -> std::io::Result<()> {
    writeln!(file, "      <Cells>")?;
    writeln!(
        file,
        "        <DataArray type=\"UInt32\" Name=\"connectivity\" format=\"ascii\">"
    )?;
    for cell in cells {
        writeln!(file, "          {} {} {}", cell[0], cell[1], cell[2])?;
    }
    writeln!(file, "        </DataArray>")?;
    writeln!(
        file,
        "        <DataArray type=\"UInt32\" Name=\"offsets\" format=\"ascii\">"
    )?;
    let mut offset = 0;
    for cell in cells {
        offset += cell.len();
        writeln!(file, "          {}", offset)?;
    }
    writeln!(file, "        </DataArray>")?;
    writeln!(
        file,
        "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"
    )?;
    write!(file, "          ")?;
    for _ in cells {
        write!(file, "5 ")?; // VTK_TRIANGLE is 5
    }
    writeln!(file, "\n        </DataArray>")?;
    writeln!(file, "      </Cells>")
}
