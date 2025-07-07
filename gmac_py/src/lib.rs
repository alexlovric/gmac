use pyo3::prelude::*;

pub mod py_mesh;
use py_mesh::PyMesh;

pub mod py_selection;
use py_selection::{
    py_select_nodes_closest_to_point, py_select_nodes_on_line, py_select_nodes_in_sphere,
    py_select_nodes_in_box, py_select_nodes_closest_to_plane,
    py_select_nodes_in_plane_direction,
};

pub mod py_transformation;
use py_transformation::{
    py_translate_nodes, py_rotate_nodes, py_scale_nodes, py_transform_nodes,
    py_build_transformation_matrix,
};

pub mod py_primitives;
use py_primitives::{
    py_generate_box, py_generate_capsule, py_generate_cone, py_generate_cylinder,
    py_generate_icosphere, py_generate_naca_wing, py_generate_torus,
    py_generate_uvsphere, py_generate_block_cluster, py_generate_sphere_cluster,
};

pub mod py_io;
use py_io::{py_write_stl, py_write_vtp, py_write_obj, py_read_stl, py_read_obj};

pub mod py_morph;
use py_morph::{PyRbfDeformer, PyDesignBlock, PyFreeFormDeformer};

#[pymodule]
fn gmac(_py: Python, m: &PyModule) -> PyResult<()> {
    // Mesh
    m.add_class::<PyMesh>()?;

    // Selection
    m.add_function(wrap_pyfunction!(py_select_nodes_closest_to_point, m)?)?;
    m.add_function(wrap_pyfunction!(py_select_nodes_on_line, m)?)?;
    m.add_function(wrap_pyfunction!(py_select_nodes_in_sphere, m)?)?;
    m.add_function(wrap_pyfunction!(py_select_nodes_in_box, m)?)?;
    m.add_function(wrap_pyfunction!(py_select_nodes_closest_to_plane, m)?)?;
    m.add_function(wrap_pyfunction!(py_select_nodes_in_plane_direction, m)?)?;

    // Transformation
    m.add_function(wrap_pyfunction!(py_translate_nodes, m)?)?;
    m.add_function(wrap_pyfunction!(py_rotate_nodes, m)?)?;
    m.add_function(wrap_pyfunction!(py_scale_nodes, m)?)?;
    m.add_function(wrap_pyfunction!(py_transform_nodes, m)?)?;
    m.add_function(wrap_pyfunction!(py_build_transformation_matrix, m)?)?;

    // Primitives
    m.add_function(wrap_pyfunction!(py_generate_box, m)?)?;
    m.add_function(wrap_pyfunction!(py_generate_capsule, m)?)?;
    m.add_function(wrap_pyfunction!(py_generate_cone, m)?)?;
    m.add_function(wrap_pyfunction!(py_generate_cylinder, m)?)?;
    m.add_function(wrap_pyfunction!(py_generate_icosphere, m)?)?;
    m.add_function(wrap_pyfunction!(py_generate_naca_wing, m)?)?;
    m.add_function(wrap_pyfunction!(py_generate_torus, m)?)?;
    m.add_function(wrap_pyfunction!(py_generate_uvsphere, m)?)?;
    m.add_function(wrap_pyfunction!(py_generate_block_cluster, m)?)?;
    m.add_function(wrap_pyfunction!(py_generate_sphere_cluster, m)?)?;

    // IO
    m.add_function(wrap_pyfunction!(py_write_stl, m)?)?;
    m.add_function(wrap_pyfunction!(py_write_vtp, m)?)?;
    m.add_function(wrap_pyfunction!(py_write_obj, m)?)?;
    m.add_function(wrap_pyfunction!(py_read_stl, m)?)?;
    m.add_function(wrap_pyfunction!(py_read_obj, m)?)?;

    // Deformers
    m.add_class::<PyRbfDeformer>()?;
    m.add_class::<PyDesignBlock>()?;
    m.add_class::<PyFreeFormDeformer>()?;

    Ok(())
}
