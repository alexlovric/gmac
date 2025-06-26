use crate::{
    mesh::Mesh,
    transformation::{build_rotation_matrix, rotate_node},
};

/// Creates a new box mesh. A box here is a shell, in that it is only the walls of
/// a block.
///
/// # Arguments
/// * `length`: An array `[x, y, z]` that specifies the dimensions of the box.
/// * `centre`: An array `[x, y, z]` that specifies the origin of the box.
/// * `theta`: An array `[x, y, z]` that specifies the rotation of the box.
/// * `resolution`: An array `[x, y, z]` that specifies the number of nodes
///                 along each axis.
///
/// # Returns
/// Returns a `Result` containing the `Mesh` instance or an error message.
///
/// # Example
/// ```
/// let box = generate_box([1.0, 1.0, 1.0], [0.0, 0.0, 0.0], [2, 2, 2]);
/// ```    
pub fn generate_box(
    length: [f64; 3],
    centre: [f64; 3],
    theta: [f64; 3],
    resolution: [usize; 3],
) -> Mesh {
    // Validation checks
    if length.iter().any(|&ll| ll <= 0.0) {
        panic!("Invalid length: dimensions must be positive");
    }
    if resolution.iter().any(|&rr| rr == 0) {
        panic!("Invalid resolution: dimensions must be non-zero");
    }

    let [rx, ry, rz] = resolution;

    let step = [
        length[0] / (rx as f64),
        length[1] / (ry as f64),
        length[2] / (rz as f64),
    ];

    let start = [
        centre[0] - length[0] * 0.5,
        centre[1] - length[1] * 0.5,
        centre[2] - length[2] * 0.5,
    ];

    let rotation_matrix = build_rotation_matrix(&theta);

    let n_nodes = (rx + 1) * (ry + 1) * (rz + 1) - (rx - 1) * (ry - 1) * (rz - 1);
    let mut nodes = vec![[0.0; 3]; n_nodes];

    let mut sides = vec![
        Vec::with_capacity((ry + 1) * (rz + 1)), // Min x
        Vec::with_capacity((ry + 1) * (rz + 1)), // Max x
        Vec::with_capacity((rx + 1) * (rz + 1)), // Min y
        Vec::with_capacity((rx + 1) * (rz + 1)), // Max y
        Vec::with_capacity((rx + 1) * (ry + 1)), // Min z
        Vec::with_capacity((rx + 1) * (ry + 1)), // Max z
    ];

    let mut node_count = 0;

    for ii in 0..=rx {
        for jj in 0..=ry {
            for kk in 0..=rz {
                // Add nodes on the boundary
                if ii == 0 || ii == rx || jj == 0 || jj == ry || kk == 0 || kk == rz {
                    nodes[node_count] = [
                        ii as f64 * step[0] + start[0],
                        jj as f64 * step[1] + start[1],
                        kk as f64 * step[2] + start[2],
                    ];

                    if theta.iter().any(|angle| *angle != 0.0) {
                        rotate_node(&mut nodes[node_count], &rotation_matrix, &centre);
                    }

                    if ii == 0 {
                        sides[0].push(node_count);
                    }
                    if ii == rx {
                        sides[1].push(node_count);
                    }
                    if jj == 0 {
                        sides[2].push(node_count);
                    } else if jj == ry {
                        sides[3].push(node_count);
                    }
                    if kk == 0 {
                        sides[4].push(node_count);
                    } else if kk == rz {
                        sides[5].push(node_count);
                    }

                    node_count += 1;
                }
            }
        }
    }

    // Cells
    let n_cells = 4 * (rx * ry + ry * rz + rx * rz);
    let mut cells: Vec<[usize; 3]> = Vec::with_capacity(n_cells);

    // X - min/max
    for jj in 0..resolution[1] {
        for kk in 0..resolution[2] {
            let p1 = sides[0][jj * (rz + 1) + kk];
            let p2 = sides[0][jj * (rz + 1) + (kk + 1)];
            let p3 = sides[0][(jj + 1) * (rz + 1) + kk];
            let p4 = sides[0][(jj + 1) * (rz + 1) + (kk + 1)];

            cells.push([p1, p2, p3]);
            cells.push([p3, p2, p4]);

            let p1 = sides[1][jj * (rz + 1) + kk];
            let p2 = sides[1][jj * (rz + 1) + (kk + 1)];
            let p3 = sides[1][(jj + 1) * (rz + 1) + kk];
            let p4 = sides[1][(jj + 1) * (rz + 1) + (kk + 1)];

            cells.push([p1, p3, p4]);
            cells.push([p4, p2, p1]);
        }
    }

    // Y - min/max
    for jj in 0..resolution[0] {
        for kk in 0..resolution[2] {
            let p1 = sides[2][jj * (rz + 1) + kk];
            let p2 = sides[2][jj * (rz + 1) + (kk + 1)];
            let p3 = sides[2][(jj + 1) * (rz + 1) + kk];
            let p4 = sides[2][(jj + 1) * (rz + 1) + (kk + 1)];

            cells.push([p1, p3, p4]);
            cells.push([p4, p2, p1]);

            let p1 = sides[3][jj * (rz + 1) + kk];
            let p2 = sides[3][jj * (rz + 1) + (kk + 1)];
            let p3 = sides[3][(jj + 1) * (rz + 1) + kk];
            let p4 = sides[3][(jj + 1) * (rz + 1) + (kk + 1)];

            cells.push([p1, p3, p2]);
            cells.push([p3, p4, p2]);
        }
    }

    // Z - min/max
    for jj in 0..resolution[0] {
        for kk in 0..resolution[1] {
            let p1 = sides[4][jj * (ry + 1) + kk];
            let p2 = sides[4][jj * (ry + 1) + (kk + 1)];
            let p3 = sides[4][(jj + 1) * (ry + 1) + kk];
            let p4 = sides[4][(jj + 1) * (ry + 1) + (kk + 1)];

            cells.push([p1, p2, p3]);
            cells.push([p3, p2, p4]);

            let p1 = sides[5][jj * (ry + 1) + kk];
            let p2 = sides[5][jj * (ry + 1) + (kk + 1)];
            let p3 = sides[5][(jj + 1) * (ry + 1) + kk];
            let p4 = sides[5][(jj + 1) * (ry + 1) + (kk + 1)];

            cells.push([p1, p3, p4]);
            cells.push([p4, p2, p1]);
        }
    }

    Mesh::new(nodes, cells)
}

/// Creates a new NACA wing mesh.
///
/// # Arguments
/// * `maximum_camber`: The maximum camber of the airfoil, usually expressed as a
///                     percentage of the chord length.
/// * `camber_distance`: The distance from the leading edge to the location of
///                      maximum camber, usually as a percentage of chord length.
/// * `maximum_thickness`: The maximum thickness of the airfoil, usually as a
///                        percentage of chord length.
/// * `n_points`: The number of points to generate for each half of the airfoil
///               (upper and lower).
/// * `wing_span`: A tuple representing the z-coordinates for the beginning and
///                end of the wing span.
///
/// # Returns
/// A `Mesh` struct containing the generated nodes and cells for the NACA wing mesh.
///
pub fn generate_naca_wing(
    maximum_camber: f64,
    camber_distance: f64,
    maximum_thickness: f64,
    n_points: usize,
    wing_span: (f64, f64),
) -> Mesh {
    let nodes_2d = generate_naca_profile(
        maximum_camber,
        camber_distance,
        maximum_thickness,
        n_points,
        false,
    );

    let nodes = nodes_2d
        .iter()
        .flat_map(|&[x, y, _]| vec![[x, y, wing_span.0], [x, y, wing_span.1]])
        .collect::<Vec<[f64; 3]>>();

    let mut cells = Vec::with_capacity(4 * (n_points - 1));

    for i in 0..(n_points) {
        let i0 = 2 * i;
        let i1 = i0 + 1;
        let i2 = i0 + 2;
        let i3 = i0 + 3;

        // Upper surface
        cells.push([i0, i2, i1]); // First triangle
        cells.push([i3, i1, i2]); // Second triangle
    }

    for i in (n_points + 1)..2 * (n_points) {
        let i0 = 2 * i;
        let i1 = i0 + 1;
        let i2 = i0 + 2;
        let i3 = i0 + 3;

        // Upper surface
        cells.push([i0, i2, i1]); // First triangle
        cells.push([i3, i1, i2]); // Second triangle
    }

    // Trailing edge
    cells.push([2 * (n_points + 1), 0, 2 * (n_points + 1) + 1]);
    cells.push([1, 2 * (n_points + 1) + 1, 0]);

    Mesh::new(nodes, cells)
}

fn generate_naca_profile(
    maximum_camber: f64,
    camber_distance: f64,
    maximum_thickness: f64,
    n_points: usize,
    finite_trailing_edge: bool,
) -> Vec<[f64; 3]> {
    let mm = 0.01 * maximum_camber;
    let pp = 0.10 * camber_distance;
    let tt = 0.01 * maximum_thickness;

    let a0 = 0.2969;
    let a1 = -0.1260;
    let a2 = -0.3516;
    let a3 = 0.2843;
    let a4 = if finite_trailing_edge {
        -0.1015
    } else {
        -0.1036
    };

    let mut xu = Vec::with_capacity(n_points + 1);
    let mut yu = Vec::with_capacity(n_points + 1);
    let mut xl = Vec::new();
    let mut yl = Vec::new();

    for i in (0..=n_points).rev() {
        let x = i as f64 / n_points as f64;
        let yt = 5.0
            * tt
            * (a0 * x.sqrt() + a1 * x + a2 * x.powi(2) + a3 * x.powi(3) + a4 * x.powi(4));

        let (yc, dyc_dx) = if pp == 0.0 {
            (0.0, 0.0) // Replace with whatever is appropriate when pp is zero
        } else if x <= pp {
            (
                mm / pp.powi(2) * x * (2.0 * pp - x),
                mm / pp.powi(2) * (2.0 * pp - 2.0 * x),
            )
        } else {
            (
                mm / (1.0 - pp).powi(2) * (1.0 - 2.0 * pp + x) * (1.0 - x),
                mm / (1.0 - pp).powi(2) * (2.0 * pp - 2.0 * x),
            )
        };

        let theta = dyc_dx.atan();

        xu.push(x - yt * theta.sin());
        yu.push(yc + yt * theta.cos());

        xl.push(x + yt * theta.sin());
        yl.push(yc - yt * theta.cos());
    }

    let combined: Vec<[f64; 3]> = xu
        .iter()
        .chain(xl.iter().skip(1))
        .zip(yu.iter().chain(yl.iter().skip(1)))
        .map(|(&x, &y)| [x - 0.25, y, 0.0])
        .collect();

    combined
}

#[cfg(test)]
mod tests {
    use super::generate_box;

    #[test]
    fn test_generate_box() {
        let uniform_box =
            generate_box([1.0, 1.0, 1.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1, 1, 1]);

        assert_eq!(uniform_box.nodes.len(), 8);

        for node in &uniform_box.nodes {
            assert_eq!(node.len(), 3);
        }

        assert_eq!(uniform_box.cells.len(), 12);

        for cell in &uniform_box.cells {
            assert_eq!(cell.len(), 3);
        }
    }
}
