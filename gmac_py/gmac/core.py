# This file defines the gmac.io submodule.

from .gmac import (
    # Mesh
    Mesh,

    # Selection
    select_nodes_closest_to_point,
    select_nodes_closest_to_plane,
    select_nodes_in_plane_direction,
    select_nodes_on_line,
    select_nodes_in_sphere,
    select_nodes_in_box,

    # Transformations
    translate_nodes,
    rotate_nodes,
    scale_nodes,
    transform_nodes,
    transform_selected_nodes,
    build_transformation_matrix,
  
    # Primitives
    generate_box,
    generate_capsule,
    generate_cone,
    generate_cylinder,
    generate_icosphere,
    generate_naca_wing,
    generate_torus,
    generate_uvsphere,
    generate_block_cluster,
    generate_sphere_cluster,
)

# Other python specific functions
def info():
    print("This is the GMAC core submodule.")