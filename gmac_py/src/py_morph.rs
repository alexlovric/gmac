#![allow(non_local_definitions)]
use pyo3::{prelude::*, exceptions::PyValueError};

use gmac::morph::{design_block::DesignBlock, ffd::FreeFormDeformer, rbf::RbfDeformer};

use crate::py_mesh::PyMesh;

#[derive(Clone, Debug)]
#[pyclass(name = "DesignBlock")]
pub struct PyDesignBlock {
    pub inner: DesignBlock,
}

#[pymethods]
impl PyDesignBlock {
    #[new]
    pub fn new(
        length: [f64; 3],
        centre: [f64; 3],
        theta: [f64; 3],
        resolution: [usize; 3],
    ) -> Self {
        PyDesignBlock {
            inner: DesignBlock::new(length, centre, theta, resolution),
        }
    }

    #[getter]
    pub fn nodes(&self) -> Vec<[f64; 3]> {
        self.inner.nodes.clone()
    }

    #[getter]
    pub fn length(&self) -> [f64; 3] {
        self.inner.length
    }

    #[getter]
    pub fn centre(&self) -> [f64; 3] {
        self.inner.centre
    }

    #[getter]
    pub fn theta(&self) -> [f64; 3] {
        self.inner.theta
    }

    #[getter]
    pub fn resolution(&self) -> [usize; 3] {
        self.inner.resolution
    }

    pub fn select_free_design_nodes(
        &self,
        target_mesh: &PyAny,
        fixed_layers: Option<usize>,
    ) -> PyResult<Vec<usize>> {
        let mesh = target_mesh.extract::<PyMesh>()?;
        match self
            .inner
            .select_free_design_nodes(&mesh.into(), fixed_layers)
        {
            Ok(result) => Ok(result),
            Err(err_str) => Err(PyValueError::new_err(err_str)),
        }
    }
}

impl From<PyDesignBlock> for DesignBlock {
    fn from(py_design_block: PyDesignBlock) -> Self {
        py_design_block.inner
    }
}

impl From<DesignBlock> for PyDesignBlock {
    fn from(design_block: DesignBlock) -> Self {
        PyDesignBlock {
            inner: design_block,
        }
    }
}

impl From<&PyAny> for PyDesignBlock {
    /// Converter to `PyDesignBlock`.
    fn from(py_design_block: &PyAny) -> Self {
        if let Ok(py_design_block) = py_design_block.extract::<PyDesignBlock>() {
            py_design_block.clone()
        } else {
            panic!("Unknown design block!")
        }
    }
}

#[pyclass(name = "FreeFormDeformer")]
pub struct PyFreeFormDeformer {
    pub ffd: FreeFormDeformer,
}

#[pymethods]
impl PyFreeFormDeformer {
    #[new]
    pub fn new(original_design_block: &PyAny) -> Self {
        let design_block = original_design_block.extract::<PyDesignBlock>().unwrap();
        PyFreeFormDeformer {
            ffd: FreeFormDeformer::new(design_block.inner),
        }
    }

    fn deform(
        &self,
        points: Vec<[f64; 3]>,
        deformed_design_nodes: Vec<[f64; 3]>,
    ) -> PyResult<Vec<[f64; 3]>> {
        match self.ffd.deform(&points, &deformed_design_nodes) {
            Ok(result) => Ok(result),
            Err(err_str) => Err(PyValueError::new_err(err_str)),
        }
    }
}

/// Rbf deformer.
#[pyclass(name = "RbfDeformer")]
pub struct PyRbfDeformer {
    pub inner: RbfDeformer,
}

#[pymethods]
impl PyRbfDeformer {
    #[new]
    pub fn new(
        original_control_points: &PyAny,
        deformed_control_points: &PyAny,
        kernel: Option<&str>,
        epsilon: Option<f64>,
    ) -> Self {
        PyRbfDeformer {
            inner: RbfDeformer::new(
                original_control_points.extract().unwrap(),
                deformed_control_points.extract().unwrap(),
                Some(kernel.unwrap_or("gaussian")),
                epsilon,
            )
            .unwrap(),
        }
    }

    fn deform(&self, points: Vec<[f64; 3]>) -> PyResult<Vec<[f64; 3]>> {
        match self.inner.deform(&points) {
            Ok(result) => Ok(result),
            Err(err_str) => Err(PyValueError::new_err(err_str)),
        }
    }
}
