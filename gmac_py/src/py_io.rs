use pyo3::prelude::*;

use gmac::io::{stl::write_stl, vtk::write_vtp};

/// Write Ascii stl file
#[pyfunction(name = "write_stl")]
pub fn py_write_stl(
    nodes: Vec<[f64; 3]>,
    cells: Vec<[usize; 3]>,
    filename: Option<&str>,
) -> PyResult<()> {
    write_stl(&nodes, &cells, filename).unwrap();
    Ok(())
}

/// Write Ascii vtp file
#[pyfunction(name = "write_vtp")]
pub fn py_write_vtp(nodes: Vec<[f64; 3]>, filename: Option<&str>) -> PyResult<()> {
    write_vtp(&nodes, filename).unwrap();
    Ok(())
}
