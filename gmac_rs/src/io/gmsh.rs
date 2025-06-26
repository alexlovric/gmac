use std::fs::File;
use std::io::{BufRead, BufReader, Error};

/// Reads a 3D mesh from a GMSH file.
///
/// # Arguments
/// * `filename`: The path of the file to read from.
///
/// # Returns
/// Returns a `Result` which contains a tuple (`Vec<[f64; 3]>`, `Vec<[usize; 3]>`)
/// if the file is successfully read, or contains an error otherwise.
#[allow(clippy::type_complexity)]
pub fn read_gmsh(filename: &str) -> Result<(Vec<[f64; 3]>, Vec<[usize; 3]>), Error> {
    let mut nodes = Vec::new();
    let mut elements = Vec::new();
    let mut phys_names = false;

    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let mut is_reading_nodes = false;
    let mut is_reading_elements = false;

    for line in reader.lines() {
        let line = line?;
        if line == "$Nodes" {
            is_reading_nodes = true;
            continue;
        } else if line == "$EndNodes" {
            is_reading_nodes = false;
        } else if line == "$Elements" {
            is_reading_elements = true;
            continue;
        } else if line == "$EndElements" {
            is_reading_elements = false;
        } else if line == "$PhysicalNames" {
            phys_names = true;
        }

        if is_reading_nodes {
            let values: Vec<f64> = line
                .split_whitespace()
                .skip(1) // skip index
                .map(|s| s.parse().unwrap())
                .collect();
            if values.len() >= 3 {
                nodes.push([values[0], values[1], values[2]]);
            }
        } else if is_reading_elements {
            let values: Vec<usize> = line
                .split_whitespace()
                .skip(1) // skip element type and other metadata
                .map(|s| s.parse().unwrap())
                .collect();
            if values.len() >= 3 {
                elements.push([values[0] - 1, values[1] - 1, values[2] - 1]);
            }
        }
    }

    if !phys_names {
        println!("Warning: No physical names found in the mesh.");
    }

    Ok((nodes, elements))
}
