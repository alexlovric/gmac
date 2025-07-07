use crate::{
    core::transformation::rotate_nodes,
    morph::{
        linear_algebra::least_squares_solver,
        transforms::{apply_affine_transform, apply_bernstein_transform},
        design_block::DesignBlock,
    },
    error::{Error, Result},
};

const REFERENCE_COORDINATE_SYSTEM: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];

/// Free-form deformer
pub struct FreeFormDeformer {
    /// Input control node block.
    pub original_design_block: DesignBlock,
    /// Affine weights pre transformation
    pub affine_weights_pre: [[f64; 4]; 4],
    /// Affine weights post transformation
    pub affine_weights_post: [[f64; 4]; 4],
}

impl FreeFormDeformer {
    /// Instantiates a new `FreeFormDeformer` instance.
    ///
    /// # Arguments
    /// * `original_design_block`: A `DesignBlock` containing the original node positions.
    ///
    /// # Returns
    /// A new `FreeFormDeformer` instance or an error string if creation fails.
    pub fn new(original_design_block: DesignBlock) -> Result<Self> {
        // Evaluate affine weights
        let affine_weights_pre = evaluate_affine_weights(
            &original_design_block.local_coordinate_system,
            &REFERENCE_COORDINATE_SYSTEM,
        )?;

        let affine_weights_post = evaluate_affine_weights(
            &REFERENCE_COORDINATE_SYSTEM,
            &original_design_block.local_coordinate_system,
        )?;

        Ok(FreeFormDeformer {
            original_design_block,
            affine_weights_pre,
            affine_weights_post,
        })
    }

    /// Deform points based on deltas of the control points.
    ///
    /// # Arguments
    /// * `points`: New input mesh nodes to predict the deformed positions for.
    /// * `deformed_design_nodes`: New nodal positions of design box.
    ///
    /// # Returns
    /// New points or an error string if deformation fails.
    pub fn deform(
        &self,
        points: &[[f64; 3]],
        deformed_design_nodes: &[[f64; 3]],
    ) -> Result<Vec<[f64; 3]>> {
        if self.original_design_block.nodes.len() != deformed_design_nodes.len() {
            return Err(Error::Deformation(
                "The number of original and deformed design points must match!"
                    .to_string(),
            ));
        }

        // Apply affine transformation to reference
        let corner_point = self.original_design_block.corner_nodes[0];
        let scaled_points = apply_affine_transform(
            &points
                .iter()
                .map(|pt| {
                    [
                        pt[0] - corner_point[0],
                        pt[1] - corner_point[1],
                        pt[2] - corner_point[2],
                    ]
                })
                .collect::<Vec<[f64; 3]>>(),
            &self.affine_weights_pre,
        )?;

        let point_ids_in_unit_cube = indices_inside_unit_cube(&scaled_points);

        let scaled_points_in_cube = point_ids_in_unit_cube
            .iter()
            .map(|&index| scaled_points[index])
            .collect::<Vec<[f64; 3]>>();

        // Evaluate deltas
        let deltas = evaluate_deltas(
            &self.original_design_block.nodes,
            deformed_design_nodes,
            self.original_design_block.theta,
            self.original_design_block.centre,
            self.original_design_block.scaling_factors,
        );

        // Apply Bernstein transform
        let shifted_points = apply_bernstein_transform(
            &scaled_points_in_cube,
            &deltas,
            &self.original_design_block.resolution,
        )?;

        // Reverse affine transform
        let reverse_scaled_points =
            apply_affine_transform(&shifted_points, &self.affine_weights_post)?;

        let mut new_points = points.to_vec();

        point_ids_in_unit_cube
            .iter()
            .zip(reverse_scaled_points)
            .for_each(|(&id, pt)| {
                new_points[id] = [
                    pt[0] + corner_point[0],
                    pt[1] + corner_point[1],
                    pt[2] + corner_point[2],
                ]
            });

        Ok(new_points)
    }
}

/// Evaluates the affine transformation weights between two coordinate systems.
///
/// # Arguments
/// * `source_coordinate_system`: A 4x3 array representing the source coordinate system.
///                               Each row is a point in the source space.
/// * `target_coordinate_system`: A 4x3 array representing the target coordinate system.
///                               Each row is a point in the target space.
///
/// # Returns
/// Returns a 4x4 array representing the affine transformation weights.
/// Otherwise, returns an error message as a `String`.
///
/// # Example
/// ```
/// let source = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]];
/// let target = [[1.5, 0.0, 0.0], [0.0, 1.5, 0.0], [0.0, 0.0, 1.5], [0.0, 0.0, 0.0]];
/// let result = evaluate_affine_weights(&source, &target);
/// // ... (handle the Result)
/// ```
fn evaluate_affine_weights(
    source_coordinate_system: &[[f64; 3]; 4],
    target_coordinate_system: &[[f64; 3]; 4],
) -> Result<[[f64; 4]; 4]> {
    if source_coordinate_system.len() != target_coordinate_system.len() {
        return Err(Error::Deformation(
            "source_coordinate_system and target_coordinate_system must be of the same size!".to_string(),
        ));
    }

    let padded_source_coordinate_system = source_coordinate_system
        .iter()
        .map(|&[x1, x2, x3]| vec![x1, x2, x3, 1.0])
        .collect::<Vec<Vec<f64>>>();

    let padded_target_coordinate_system = target_coordinate_system
        .iter()
        .map(|&[x1, x2, x3]| vec![x1, x2, x3, 1.0])
        .collect::<Vec<Vec<f64>>>();

    let a = least_squares_solver(
        &padded_source_coordinate_system,
        &padded_target_coordinate_system,
    )?;

    let mut arr = [[0.0; 4]; 4];
    for (i, row) in a.iter().enumerate() {
        for (j, &value) in row.iter().enumerate() {
            arr[i][j] = value;
        }
    }

    Ok(arr)
}

/// Find the indices of points that are inside the unit cube `[0, 1]^3`.
///
/// # Arguments
/// * `points`: A vector of points in 3D space, represented as arrays `[x, y, z]`.
///
/// # Returns
/// * `Vec<usize>`: A vector of indices of points that are inside the unit cube.
fn indices_inside_unit_cube(points: &[[f64; 3]]) -> Vec<usize> {
    let mut indx_points_inside = Vec::new();

    for (index, &point) in points.iter().enumerate() {
        if point.iter().all(|&coord| (0.0..=1.0).contains(&coord)) {
            indx_points_inside.push(index);
        }
    }

    indx_points_inside
}

/// Evaluate the deltas between original and deformed design nodes.
///
/// # Arguments
/// * `original_nodes`: A vector of original node coordinates. Each node is represented as [f64; 3].
/// * `deformed_nodes`: A vector of deformed node coordinates. Each node is represented as [f64; 3].
/// * `theta`: The rotation vector as [f64; 3].
/// * `centre`: The centre of rotation as [f64; 3].
/// * `scaling_factors`: The scaling factors as [f64; 3].
///
/// # Returns
/// A vector of deltas. Each delta is represented as [f64; 3].
pub fn evaluate_deltas(
    original_nodes: &[[f64; 3]],
    deformed_nodes: &[[f64; 3]],
    theta: [f64; 3],
    centre: [f64; 3],
    scaling_factors: [f64; 3],
) -> Vec<[f64; 3]> {
    let mut original = original_nodes.to_vec();
    let mut deformed = deformed_nodes.to_vec();
    let neg_theta = [-theta[0], -theta[1], -theta[2]];

    rotate_nodes(&mut original, &neg_theta, &centre);
    rotate_nodes(&mut deformed, &neg_theta, &centre);

    original
        .iter()
        .zip(deformed.iter())
        .map(|(&ori, &def)| {
            [
                scaling_factors[0] * (def[0] - ori[0]),
                scaling_factors[1] * (def[1] - ori[1]),
                scaling_factors[2] * (def[2] - ori[2]),
            ]
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_evaluate_deltas() {
        let original_nodes = &[[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]];

        let deformed_nodes = &[[2.0, 2.0, 2.0], [3.0, 3.0, 3.0]];

        let theta = [0.0, 0.0, 0.0];
        let centre = [0.0, 0.0, 0.0];
        let scaling_factors = [1.0, 1.0, 1.0];

        let deltas = evaluate_deltas(
            original_nodes,
            deformed_nodes,
            theta,
            centre,
            scaling_factors,
        );

        assert_eq!(deltas, vec![[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]);
    }
}
