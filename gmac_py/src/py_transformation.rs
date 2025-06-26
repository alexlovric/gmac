use pyo3::prelude::*;

use gmac_core::transformation::{
    translate_nodes, rotate_nodes, scale_nodes, transform_nodes,
    build_transformation_matrix,
};

/// Translate nodes
#[pyfunction(name = "translate_nodes")]
pub fn py_translate_nodes(
    mut nodes: Vec<[f64; 3]>,
    translation: [f64; 3],
) -> PyResult<Vec<[f64; 3]>> {
    translate_nodes(&mut nodes, &translation);
    Ok(nodes)
}

/// Rotate nodes
#[pyfunction(name = "rotate_nodes")]
pub fn py_rotate_nodes(
    mut nodes: Vec<[f64; 3]>,
    rotation: [f64; 3],
    origin: [f64; 3],
) -> PyResult<Vec<[f64; 3]>> {
    rotate_nodes(&mut nodes, &rotation, &origin);
    Ok(nodes)
}

/// Scale nodes
#[pyfunction(name = "scale_nodes")]
pub fn py_scale_nodes(
    mut nodes: Vec<[f64; 3]>,
    scaling: [f64; 3],
    origin: [f64; 3],
) -> PyResult<Vec<[f64; 3]>> {
    scale_nodes(&mut nodes, &scaling, &origin);
    Ok(nodes)
}

/// Transform nodes
#[pyfunction(name = "transform_nodes")]
pub fn py_transform_nodes(
    mut nodes: Vec<[f64; 3]>,
    transformation_matrix: [[f64; 4]; 4],
    origin: [f64; 3],
) -> PyResult<Vec<[f64; 3]>> {
    transform_nodes(&mut nodes, &transformation_matrix, &origin);
    Ok(nodes)
}

/// Transformation matrix
#[pyfunction(name = "build_transformation_matrix")]
pub fn py_build_transformation_matrix(
    translation: [f64; 3],
    rotation: [f64; 3],
    scaling: [f64; 3],
) -> PyResult<[[f64; 4]; 4]> {
    Ok(build_transformation_matrix(translation, rotation, scaling))
}
