use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use crate::error::{Result};

#[cfg(feature = "rayon")]
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};

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
pub fn write_vtu(
    nodes: &[[f64; 3]],
    cells: &[[usize; 3]],
    filename: Option<&str>,
) -> Result<()> {
    // Pre-format all data into strings using the helpers
    let nodes_str = format_nodes_string(nodes);
    let (conn_str, offsets_str, types_str) = format_cells_strings(cells);

    // Write the pre-formatted strings to a buffered file
    let mut writer = BufWriter::new(File::create(filename.unwrap_or("mesh.vtu"))?);

    writeln!(writer, r#"<?xml version="1.0"?>"#)?;
    writeln!(
        writer,
        r#"<VTKFile type="UnstructuredGrid" version="1.1" byte_order="LittleEndian">"#
    )?;
    writeln!(writer, "  <UnstructuredGrid>")?;
    writeln!(
        writer,
        "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">",
        nodes.len(),
        cells.len()
    )?;

    // Write the Points block
    write_points_block(&mut writer, &nodes_str)?;

    // Write the Cells block
    writeln!(writer, "      <Cells>")?;
    write_cells_block(&mut writer, &conn_str, &offsets_str, &types_str)?;
    writeln!(writer, "      </Cells>")?;

    writeln!(writer, "    </Piece>")?;
    writeln!(writer, "  </UnstructuredGrid>")?;
    writeln!(writer, "</VTKFile>")?;

    Ok(())
}

/// Writes a point cloud to a VTP file, using parallel computation if the `rayon` feature is enabled.
/// This function stores a list of vertices, which is suitable for visualizing point clouds.
///
/// # Arguments
/// * `nodes` - A slice of 3D vertex positions.
/// * `filename` - Optional file path. Defaults to `"points.vtp"` if `None`.
///
/// # Returns
/// Returns `Ok(())` on success, or an `std::io::Error` on failure.
pub fn write_vtp(nodes: &[[f64; 3]], filename: Option<&str>) -> Result<()> {
    // Pre-format data using helpers
    let nodes_str = format_nodes_string(nodes);
    let (conn_str, offsets_str) = format_point_cloud_cells_strings(nodes.len());

    // Write the XML structure to a buffered file
    let mut writer = BufWriter::new(File::create(filename.unwrap_or("points.vtp"))?);

    writeln!(writer, r#"<?xml version="1.0"?>"#)?;
    writeln!(
        writer,
        r#"<VTKFile type="PolyData" version="1.1" byte_order="LittleEndian">"#
    )?;
    writeln!(writer, "  <PolyData>")?;
    writeln!(
        writer,
        "    <Piece NumberOfPoints=\"{}\" NumberOfVerts=\"{}\">",
        nodes.len(),
        nodes.len()
    )?;

    // Write data blocks
    write_points_block(&mut writer, &nodes_str)?;
    write_verts_block(&mut writer, &conn_str, &offsets_str)?;

    writeln!(writer, "    </Piece>")?;
    writeln!(writer, "  </PolyData>")?;
    writeln!(writer, "</VTKFile>")?;

    Ok(())
}

/// (Internal Helper) Formats the node data into a single string.
fn format_nodes_string(nodes: &[[f64; 3]]) -> String {
    let iter = {
        #[cfg(feature = "rayon")]
        {
            nodes.par_iter()
        }
        #[cfg(not(feature = "rayon"))]
        {
            nodes.iter()
        }
    };
    iter.map(|p| format!("          {} {} {}", p[0], p[1], p[2]))
        .collect::<Vec<_>>()
        .join("\n")
}

/// (Internal Helper) Formats cell data for a mesh (<Polys>).
fn format_cells_strings(cells: &[[usize; 3]]) -> (String, String, String) {
    let (conn_iter, offsets_iter, types_iter) = {
        #[cfg(feature = "rayon")]
        {
            (
                cells.par_iter(),
                (0..cells.len()).into_par_iter(),
                cells.par_iter(),
            )
        }
        #[cfg(not(feature = "rayon"))]
        {
            (cells.iter(), (0..cells.len()), cells.iter())
        }
    };

    let conn_str = conn_iter
        .map(|c| format!("{} {} {}", c[0], c[1], c[2]))
        .collect::<Vec<_>>()
        .join(" ");

    let offsets_str = offsets_iter
        .map(|i| ((i + 1) * 3).to_string())
        .collect::<Vec<_>>()
        .join(" ");

    let types_str = types_iter.map(|_| "5").collect::<Vec<_>>().join(" ");

    (conn_str, offsets_str, types_str)
}

/// (Internal Helper) Formats cell data for a point cloud (<Verts>).
fn format_point_cloud_cells_strings(num_nodes: usize) -> (String, String) {
    let (conn_iter, offsets_iter) = {
        #[cfg(feature = "rayon")]
        {
            (
                (0..num_nodes).into_par_iter(),
                (1..=num_nodes).into_par_iter(),
            )
        }
        #[cfg(not(feature = "rayon"))]
        {
            ((0..num_nodes), (1..=num_nodes))
        }
    };

    let conn_str = conn_iter
        .map(|i| i.to_string())
        .collect::<Vec<_>>()
        .join(" ");
    let offsets_str = offsets_iter
        .map(|i| i.to_string())
        .collect::<Vec<_>>()
        .join(" ");

    (conn_str, offsets_str)
}

/// (Internal Helper) Writes the complete <Points> block to the file.
fn write_points_block(
    writer: &mut BufWriter<File>,
    nodes_str: &str,
) -> std::io::Result<()> {
    writeln!(writer, "      <Points>")?;
    writeln!(writer, "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">")?;
    writeln!(writer, "{}", nodes_str)?;
    writeln!(writer, "        </DataArray>")?;
    writeln!(writer, "      </Points>")
}

/// (Internal Helper) Writes the complete <Verts> block for a point cloud.
fn write_verts_block(
    writer: &mut BufWriter<File>,
    conn_str: &str,
    offsets_str: &str,
) -> std::io::Result<()> {
    writeln!(writer, "      <Verts>")?;
    writeln!(
        writer,
        "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">"
    )?;
    writeln!(writer, "          {}", conn_str)?;
    writeln!(writer, "        </DataArray>")?;
    writeln!(
        writer,
        "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">"
    )?;
    writeln!(writer, "          {}", offsets_str)?;
    writeln!(writer, "        </DataArray>")?;
    writeln!(writer, "      </Verts>")
}

/// (Internal Helper) Writes the complete <Cells> or <Polys> block for a mesh.
/// Note: For VTP, the outer tag is <Polys>, for VTU it's <Cells>. This helper writes the inner content.
fn write_cells_block(
    writer: &mut BufWriter<File>,
    conn_str: &str,
    offsets_str: &str,
    types_str: &str,
) -> std::io::Result<()> {
    // This helper can be used for both VTU's <Cells> and VTP's <Polys>
    writeln!(
        writer,
        "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">"
    )?;
    writeln!(writer, "{}", conn_str)?;
    writeln!(writer, "        </DataArray>")?;
    writeln!(
        writer,
        "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">"
    )?;
    writeln!(writer, "          {}", offsets_str)?;
    writeln!(writer, "        </DataArray>")?;
    writeln!(
        writer,
        "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"
    )?;
    writeln!(writer, "          {}", types_str)?;
    writeln!(writer, "        </DataArray>")
}

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
