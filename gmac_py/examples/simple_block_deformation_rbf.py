"""
Simple Block Deformation using Radial Basis Functions (RBF)

This example demonstrates how to use the GMAC library to perform mesh deformation
using Radial Basis Functions (RBF). RBF is a powerful technique for creating smooth
deformations by interpolating the movement of control points throughout the mesh.
"""

import gmac
import gmac.morph as morph
import gmac.io as io
import numpy as np

# Create a simple box geometry with specified dimensions, center, orientation, and resolution
box = gmac.generate_box(
    length=[1.0, 1.0, 1.0],  # Dimensions (length, width, height)
    centre=[0.0, 0.0, 0.0],  # Center coordinates
    theta=[0.0, 0.0, 0.0],   # Rotation angles (degrees)
    resolution=[5, 5, 5],     # Number of divisions in each direction
)

# Save the original box geometry as an STL file for visualization
io.write_stl(nodes=box.nodes, cells=box.cells, filename="original_box.stl")

# Create a set of control points that will be used to drive the deformation
# These points form a lattice around the geometry
original_control_points = np.array(
    gmac.generate_block_cluster(
        length=[1.2, 1.2, 1.2],  # Slightly larger than the original box
        centre=[0.0, 0.0, 0.0],  # Centered at the origin
        theta=[0.0, 0.0, 0.0],   # No rotation
        resolution=[2, 2, 2],     # 2x2x2 grid of control points (8 points total)
    )
)

# Save the original control points for visualization
io.write_vtp(original_control_points, "original_control_points.vtp")

# Select control points that lie in a plane defined by an origin point and normal vector
# This will identify which control points we want to move to create the deformation
target_control_point_ids = gmac.select_nodes_in_plane_direction(
    nodes=original_control_points,  # The control points to select from
    origin=[0.3, 0.0, 0.0],        # A point in the plane
    normal=[1.0, 0.0, 0.0]         # Normal vector defining the plane direction (x-axis in this case)
)

# Save the selected target points for visualization
io.write_vtp(
    nodes=original_control_points[target_control_point_ids],
    filename="target_points.vtp",
)

# Apply a transformation to the selected control points
deformed_control_points = original_control_points.copy() # Create a copy of the original control points to modify
deformed_control_points[target_control_point_ids] = gmac.transform_nodes(
    nodes=deformed_control_points[target_control_point_ids],  # Only transform selected points
    transformation_matrix=gmac.build_transformation_matrix(
        translation=[1.0, 0.0, 0.0],   # Move 1 unit in x-direction
        rotation=[45.0, 0.0, 0.0],     # Rotate 45 degrees around x-axis
        scaling=[1.0, 0.75, 0.75],     # Scale y and z dimensions to 75%
    ),
    origin=[0.0, 0.0, 0.0],  # Origin point for the transformation
)

# Save the deformed control points for visualization
io.write_vtp(nodes=deformed_control_points, filename="deformed_control_points.vtp")

# Create an RBF deformer with the control point configurations
rbf = morph.RbfDeformer(
    original_control_points=original_control_points,  # Original positions of control points
    deformed_control_points=deformed_control_points,  # Target positions of control points
    kernel="gaussian",  # Type of RBF kernel to use
    epsilon=1.0,                 # Shape parameter for the RBF kernel
)

# Apply the RBF deformation to the original box geometry
box.nodes = rbf.deform(points=box.nodes)

# Save the final deformed geometry as an STL file
io.write_stl(nodes=box.nodes, cells=box.cells, filename="deformed_box.stl")