"""
Simple Block Deformation using Radial Basis Functions (RBF)

This example demonstrates how to use the GMAC library to perform mesh deformation
using Radial Basis Functions (RBF). RBF is a powerful technique for creating smooth
deformations by interpolating the movement of control points throughout the mesh.
"""

import gmac as gm
import numpy as np


# Create a simple geometry or import one
mesh = gm.generate_torus(1.0, 0.3, [0.0, 0.0, 0.0], 64, 64)

# Save the original geometry as an STL file for visualization
mesh.write_stl("original_torus.stl", "binary")

# Create a set of control points that will be used to drive the deformation
# These points form a lattice around the geometry
original_control_points = np.array(
    gm.generate_block_cluster(
        length=[2.8, 1.0, 2.8],
        centre=[0.0, 0.0, 0.0],
        theta=[0.0, 0.0, 0.0],
        resolution=[3, 3, 3],
    )
)

# Save the original control points for visualization
gm.io.write_vtp(original_control_points, "original_control_points.vtp")

# Select control points that lie in a plane defined by an origin point and normal vector
# This will identify which control points we want to move to create the deformation
target_control_point_ids = gm.select_nodes_in_plane_direction(
    nodes=original_control_points,
    origin=[0.0, 0.0, 0.0],
    normal=[1.0, 0.0, 0.0],
)

# Save the selected target points for visualization
gm.io.write_vtp(
    nodes=original_control_points[target_control_point_ids],
    filename="target_points.vtp",
)

# Apply a transformation to the selected control points
deformed_control_points = original_control_points.copy()
deformed_control_points[target_control_point_ids] = gm.transform_nodes(
    nodes=deformed_control_points[target_control_point_ids],
    transformation_matrix=gm.build_transformation_matrix(
        translation=[0.0, 0.0, 0.0],
        rotation=[45.0, 20.0, 0.0],
        scaling=[1.0, 1.0, 1.0],
    ),
    origin=[0.0, 0.0, 0.0],
)

# Save the deformed control points for visualization
gm.io.write_vtp(nodes=deformed_control_points, filename="deformed_control_points.vtp")

# Create an RBF deformer with the control point configurations
rbf = gm.morph.RbfDeformer(
    original_control_points=original_control_points,
    deformed_control_points=deformed_control_points,
    kernel="thin_plate_spline",
    epsilon=1.0,
)

# Apply the RBF deformation to the original box geometry
mesh.nodes = rbf.deform(points=mesh.nodes)

# Save the final deformed geometry as an STL file
mesh.write_stl("deformed_torus.stl", "binary")
