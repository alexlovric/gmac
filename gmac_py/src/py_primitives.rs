use pyo3::prelude::*;
use crate::py_mesh::PyMesh;

use gmac::core::{
    clusters::{generate_block_cluster, generate_sphere_cluster},
    primitives::{
        generate_box, generate_capsule, generate_cone, generate_cylinder,
        generate_icosphere, generate_naca_wing, generate_torus, generate_uvsphere,
    },
};

/// Generate box primative
#[pyfunction(name = "generate_box")]
pub fn py_generate_box(
    length: [f64; 3],
    centre: [f64; 3],
    theta: [f64; 3],
    resolution: [usize; 3],
) -> PyResult<PyMesh> {
    Ok(PyMesh::from(
        generate_box(length, centre, theta, resolution).unwrap(),
    ))
}

/// Generate icosphere primative
#[pyfunction(name = "generate_icosphere")]
pub fn py_generate_icosphere(
    radius: f64,
    centre: [f64; 3],
    subdivisions: u32,
) -> PyResult<PyMesh> {
    Ok(PyMesh::from(
        generate_icosphere(radius, centre, subdivisions).unwrap(),
    ))
}

/// Generate uv sphere primative
#[pyfunction(name = "generate_uvsphere")]
pub fn py_generate_uvsphere(
    radius: f64,
    centre: [f64; 3],
    sectors: u32,
    stacks: u32,
) -> PyResult<PyMesh> {
    Ok(PyMesh::from(
        generate_uvsphere(radius, centre, sectors, stacks).unwrap(),
    ))
}

/// Generate cylinder primative
#[pyfunction(name = "generate_cylinder")]
pub fn py_generate_cylinder(
    radius: f64,
    height: f64,
    centre: [f64; 3],
    sectors: u32,
    stacks: u32,
) -> PyResult<PyMesh> {
    Ok(PyMesh::from(
        generate_cylinder(radius, height, centre, sectors, stacks).unwrap(),
    ))
}

/// Generate torus primative
#[pyfunction(name = "generate_torus")]
pub fn py_generate_torus(
    major_radius: f64,
    minor_radius: f64,
    centre: [f64; 3],
    major_segments: u32,
    minor_segments: u32,
) -> PyResult<PyMesh> {
    Ok(PyMesh::from(
        generate_torus(
            major_radius,
            minor_radius,
            centre,
            major_segments,
            minor_segments,
        )
        .unwrap(),
    ))
}

/// Generate cone primative
#[pyfunction(name = "generate_cone")]
pub fn py_generate_cone(
    radius: f64,
    height: f64,
    centre: [f64; 3],
    sectors: u32,
) -> PyResult<PyMesh> {
    Ok(PyMesh::from(
        generate_cone(radius, height, centre, sectors).unwrap(),
    ))
}

/// Generate capsule primative
#[pyfunction(name = "generate_capsule")]
pub fn py_generate_capsule(
    radius: f64,
    cylinder_height: f64,
    centre: [f64; 3],
    sectors: u32,
    stacks: u32,
) -> PyResult<PyMesh> {
    Ok(PyMesh::from(
        generate_capsule(radius, cylinder_height, centre, sectors, stacks).unwrap(),
    ))
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
    Ok(PyMesh::from(
        generate_naca_wing(
            maximum_camber,
            camber_distance,
            maximum_thickness,
            n_points,
            wing_span,
        )
        .unwrap(),
    ))
}

/// Generate block node cluster
#[pyfunction(name = "generate_block_cluster")]
pub fn py_generate_block_cluster(
    length: [f64; 3],
    centre: [f64; 3],
    theta: [f64; 3],
    resolution: [usize; 3],
) -> PyResult<Vec<[f64; 3]>> {
    Ok(generate_block_cluster(length, centre, theta, resolution).unwrap())
}

/// Generate sphere node cluster
#[pyfunction(name = "generate_sphere_cluster")]
pub fn py_generate_sphere_cluster(
    radius: f64,
    centre: [f64; 3],
    resolution: usize,
) -> PyResult<Vec<[f64; 3]>> {
    Ok(generate_sphere_cluster(radius, centre, resolution).unwrap())
}
