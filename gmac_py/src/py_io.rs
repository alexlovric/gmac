use pyo3::prelude::*;

use gmac::io::{
    obj::{read_obj, write_obj},
    stl::{read_stl, write_stl, StlFormat},
    vtk::write_vtp,
};

/// Write Ascii stl file
#[pyfunction(name = "write_stl")]
pub fn py_write_stl(
    nodes: Vec<[f64; 3]>,
    cells: Vec<[usize; 3]>,
    filename: Option<&str>,
    format: Option<&str>,
) -> PyResult<()> {
    let format = Some(match format {
        Some("ascii") => StlFormat::Ascii,
        _ => StlFormat::Binary,
    });
    write_stl(&nodes, &cells, filename, format).unwrap();
    Ok(())
}

/// Read stl file
#[pyfunction(name = "read_stl")]
pub fn py_read_stl(filename: &str) -> PyResult<(Vec<[f64; 3]>, Vec<[usize; 3]>)> {
    let (nodes, cells) = read_stl(filename).unwrap();
    Ok((nodes, cells))
}

/// Write Ascii vtp file
#[pyfunction(name = "write_vtp")]
pub fn py_write_vtp(nodes: Vec<[f64; 3]>, filename: Option<&str>) -> PyResult<()> {
    write_vtp(&nodes, filename).unwrap();
    Ok(())
}

/// Write obj file
#[pyfunction(name = "write_obj")]
pub fn py_write_obj(
    nodes: Vec<[f64; 3]>,
    cells: Vec<[usize; 3]>,
    filename: Option<&str>,
) -> PyResult<()> {
    write_obj(&nodes, &cells, filename).unwrap();
    Ok(())
}

/// Read obj file
#[pyfunction(name = "read_obj")]
pub fn py_read_obj(filename: &str) -> PyResult<(Vec<[f64; 3]>, Vec<[usize; 3]>)> {
    let (nodes, cells) = read_obj(filename).unwrap();
    Ok((nodes, cells))
}
