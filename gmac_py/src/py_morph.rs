#![allow(non_local_definitions)]
use pyo3::{prelude::*, exceptions::PyValueError};

use gmac::morph::{
    design_block::DesignBlock,
    ffd::FreeFormDeformer,
    rbf::{
        gaussian_kernel, inverse_multi_kernel, multiquadric_kernel,
        thin_plate_spline_kernel, RbfDeformer,
    },
};

use crate::py_mesh::PyMesh;

/// Rbf deformer.
#[pyclass(name = "RbfDeformer")]
pub struct PyRbfDeformer {
    pub rbf: RbfDeformer,
    pub kernel_type: String,
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
        let kernel_type = kernel.unwrap_or("gaussian");
        let kernel = match kernel_type {
            "gaussian" => gaussian_kernel,
            "multiquadric" => multiquadric_kernel,
            "inverse_multiquadratic" => inverse_multi_kernel,
            "thin_plate_spline" => thin_plate_spline_kernel,
            _ => panic!("Kernel not implemented"),
        };

        PyRbfDeformer {
            rbf: RbfDeformer::new(
                original_control_points.extract().unwrap(),
                deformed_control_points.extract().unwrap(),
                Some(kernel),
                epsilon,
            )
            .unwrap(),
            kernel_type: String::from(kernel_type),
        }
    }

    fn deform(&self, points: Vec<[f64; 3]>) -> PyResult<Vec<[f64; 3]>> {
        match self.rbf.deform(&points) {
            Ok(result) => Ok(result),
            Err(err_str) => Err(PyValueError::new_err(err_str)),
        }
    }

    pub fn __repr__(&self) -> String {
        format!(
            "┌{}┐\n│{: <48}│\n╞{}╡\n│{} {: <40}│\n│{} {: <39}│\n└{}┘",
            "─".repeat(48),
            "Rbf Model",
            "═".repeat(48),
            "Kernel:",
            self.kernel_type,
            "Epsilon:",
            self.rbf.epsilon,
            "─".repeat(48)
        )
    }
}

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
        self.inner.length.clone()
    }

    #[getter]
    pub fn centre(&self) -> [f64; 3] {
        self.inner.centre.clone()
    }

    #[getter]
    pub fn theta(&self) -> [f64; 3] {
        self.inner.theta.clone()
    }

    #[getter]
    pub fn resolution(&self) -> [usize; 3] {
        self.inner.resolution.clone()
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
