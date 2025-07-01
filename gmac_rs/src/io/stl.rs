use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Result, Write};

pub fn write_stl(
    nodes: &[[f64; 3]],
    cells: &[[usize; 3]],
    filename: Option<&str>,
) -> Result<()> {
    let mut file = File::create(filename.unwrap_or("mesh.stl"))?;

    writeln!(file, "solid exported_grid")?;

    for cell in cells.iter() {
        let p1 = nodes[cell[0]];
        let p2 = nodes[cell[1]];
        let p3 = nodes[cell[2]];

        // Calculate normal
        let u = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
        let v = [p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]];

        let n = [
            (u[1] * v[2]) - (u[2] * v[1]),
            (u[2] * v[0]) - (u[0] * v[2]),
            (u[0] * v[1]) - (u[1] * v[0]),
        ];

        let norm = (n[0].powi(2) + n[1].powi(2) + n[2].powi(2)).sqrt();

        writeln!(
            file,
            "  facet normal {} {} {}",
            n[0] / norm,
            n[1] / norm,
            n[2] / norm
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

/// Reads an ASCII STL file and extracts triangles as `Vec<[[f64; 3]; 3]>`.
///
/// # Arguments
/// * `filename` - Path to the STL file.
///
/// # Returns
/// A vector of triangles, where each triangle is represented as 3 points `[f64; 3]`.
pub fn read_stl_ascii(filename: &str) -> Result<(Vec<[f64; 3]>, Vec<[usize; 3]>)> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

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
