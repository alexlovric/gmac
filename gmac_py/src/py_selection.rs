use pyo3::prelude::*;

use gmac_core::selection::{
    select_nodes_closest_to_point, select_nodes_on_line, select_nodes_in_sphere,
    select_nodes_in_box, select_nodes_closest_to_plane, select_nodes_in_plane_direction,
};

/// Select node closest to point
#[pyfunction(name = "select_nodes_closest_to_point")]
pub fn py_select_nodes_closest_to_point(
    nodes: Vec<[f64; 3]>,
    point: [f64; 3],
) -> PyResult<Vec<usize>> {
    Ok(select_nodes_closest_to_point(&nodes, point))
}

/// Select nodes on line
#[pyfunction(name = "select_nodes_on_line")]
pub fn py_select_nodes_on_line(
    nodes: Vec<[f64; 3]>,
    point_a: [f64; 3],
    point_b: [f64; 3],
) -> PyResult<Vec<usize>> {
    Ok(select_nodes_on_line(&nodes, point_a, point_b))
}

/// Select nodes in sphere
#[pyfunction(name = "select_nodes_in_sphere")]
pub fn py_select_nodes_in_sphere(
    nodes: Vec<[f64; 3]>,
    radius: f64,
    centre: [f64; 3],
) -> PyResult<Vec<usize>> {
    Ok(select_nodes_in_sphere(&nodes, radius, centre))
}

/// Select nodes in sphere
#[pyfunction(name = "select_nodes_in_box")]
pub fn py_select_nodes_in_box(
    nodes: Vec<[f64; 3]>,
    length: [f64; 3],
    centre: [f64; 3],
    theta: [f64; 3],
) -> PyResult<Vec<usize>> {
    Ok(select_nodes_in_box(&nodes, length, centre, theta))
}

/// Select nodes closest to plane
#[pyfunction(name = "select_nodes_closest_to_plane")]
pub fn py_select_nodes_closest_to_plane(
    nodes: Vec<[f64; 3]>,
    origin: [f64; 3],
    normal: [f64; 3],
) -> PyResult<Vec<usize>> {
    Ok(select_nodes_closest_to_plane(&nodes, origin, normal))
}

/// Select nodes in plane direction
#[pyfunction(name = "select_nodes_in_plane_direction")]
pub fn py_select_nodes_in_plane_direction(
    nodes: Vec<[f64; 3]>,
    origin: [f64; 3],
    normal: [f64; 3],
) -> PyResult<Vec<usize>> {
    Ok(select_nodes_in_plane_direction(&nodes, origin, normal))
}
