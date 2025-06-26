#![allow(non_local_definitions)]
use pyo3::prelude::*;

use gmac_core::mesh::{
    get_mesh_triangles, get_mesh_cell_normals, Mesh, clip_mesh_from_plane,
    find_mesh_intersections_with_plane,
};

use gmac_io::{stl::write_stl, vtk::write_vtp};

/// Mesh
#[derive(Clone, Debug, Default)]
#[pyclass(name = "Mesh")]
pub struct PyMesh {
    #[pyo3(get, set)]
    nodes: Vec<[f64; 3]>,
    #[pyo3(get, set)]
    cells: Vec<[usize; 3]>,
}

#[pymethods]
impl PyMesh {
    #[new]
    pub fn new(nodes: Vec<[f64; 3]>, cells: Vec<[usize; 3]>) -> Self {
        PyMesh { nodes, cells }
    }

    pub fn triangles(&self) -> PyResult<Vec<[[f64; 3]; 3]>> {
        Ok(get_mesh_triangles(&self.nodes, &self.cells))
    }

    pub fn cell_normals(&self) -> PyResult<Vec<[f64; 3]>> {
        Ok(get_mesh_cell_normals(&self.nodes, &self.cells))
    }

    pub fn slice_from_plane(
        &self,
        origin: [f64; 3],
        normal: [f64; 3],
    ) -> PyResult<Vec<[f64; 3]>> {
        Ok(
            find_mesh_intersections_with_plane(&self.nodes, &self.cells, origin, normal)
                .unwrap(),
        )
    }

    pub fn clip_from_plane(
        &self,
        origin: [f64; 3],
        normal: [f64; 3],
    ) -> PyResult<PyMesh> {
        let (new_nodes, new_cells) =
            clip_mesh_from_plane(&self.nodes, &self.cells, origin, normal);
        Ok(PyMesh::new(new_nodes, new_cells))
    }

    pub fn write_stl(&self, filename: Option<&str>) -> PyResult<()> {
        write_stl(&self.nodes, &self.cells, filename).unwrap();
        Ok(())
    }

    pub fn write_vtp(&self, filename: Option<&str>) -> PyResult<()> {
        write_vtp(&self.nodes, filename).unwrap();
        Ok(())
    }

    pub fn __repr__(&self) -> String {
        format!(
            "┌{}┐\n│{: <48}│\n╞{}╡\n│{} {: <41}│\n│{} {: <41}│\n└{}┘",
            "─".repeat(48),
            "Mesh",
            "═".repeat(48),
            "Nodes:",
            self.nodes.len(),
            "Faces:",
            self.cells.len(),
            "─".repeat(48)
        )
    }
}

impl From<Mesh> for PyMesh {
    /// Converter to `PyMesh`.
    fn from(mesh: Mesh) -> Self {
        PyMesh::new(mesh.nodes, mesh.cells)
    }
}

impl From<PyMesh> for Mesh {
    /// Converter to `PyMesh`.
    fn from(mesh: PyMesh) -> Self {
        Mesh::new(mesh.nodes, mesh.cells)
    }
}

impl From<&PyAny> for PyMesh {
    /// Converter to `PyMesh`.
    fn from(mesh: &PyAny) -> Self {
        if let Ok(mesh) =
            mesh.extract::<PyMesh>()
        {
            mesh.clone()
        } else {
            panic!("Unknown mesh!")
        }
    }
}
