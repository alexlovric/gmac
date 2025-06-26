use pyo3::prelude::*;
use crate::py_mesh::PyMesh;

use gmac_core::{
    primitives::{generate_box, generate_naca_wing},
    clusters::generate_block_cluster,
};

/// Generate block primative
#[pyfunction(name = "generate_box")]
pub fn py_generate_box(
    length: [f64; 3],
    centre: [f64; 3],
    theta: [f64; 3],
    resolution: [usize; 3],
) -> PyResult<PyMesh> {
    Ok(PyMesh::from(generate_box(
        length, centre, theta, resolution,
    )))
}

/// Generate naca wing primative
#[pyfunction(name = "generate_naca_wing")]
pub fn py_generate_naca_wing(
    maximum_camber: f64,
    camber_distance: f64,
    maximum_thickness: f64,
    n_points: usize,
    wing_span: (f64, f64),
) -> PyResult<PyMesh> {
    Ok(PyMesh::from(generate_naca_wing(
        maximum_camber,
        camber_distance,
        maximum_thickness,
        n_points,
        wing_span,
    )))
}

/// Generate block node cluster
#[pyfunction(name = "generate_block_cluster")]
pub fn py_generate_block_cluster(
    length: [f64; 3],
    centre: [f64; 3],
    theta: [f64; 3],
    resolution: [usize; 3],
) -> PyResult<Vec<[f64; 3]>> {
    Ok(generate_block_cluster(length, centre, theta, resolution))
}
