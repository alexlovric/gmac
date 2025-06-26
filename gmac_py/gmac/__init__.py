from .gmac import (
    Mesh,
    DesignBlock,
    FreeFormDeformer,
    RbfDeformer,
    select_nodes_closest_to_point,
    select_nodes_closest_to_plane,
    select_nodes_in_plane_direction,
    select_nodes_on_line,
    select_nodes_in_sphere,
    select_nodes_in_box,
    translate_nodes,
    rotate_nodes,
    scale_nodes,
    transform_nodes,
    build_transformation_matrix,
    generate_box,
    generate_naca_wing,
    generate_block_cluster,
)


# def mesh_from_trimesh(trimesh: object) -> Mesh:
#     """
#     Get mesh from trimesh object.

#     Parameters
#     ----------
#     trimesh : object

#     Returns
#     -------
#     Mesh
#         Gmac mesh.
#     """
#     return Mesh(trimesh.vertices, trimesh.faces)

# gmac.__all__.append("mesh_from_trimesh")

__doc__ = gmac.__doc__
if hasattr(gmac, "__all__"):
    __all__ = gmac.__all__
