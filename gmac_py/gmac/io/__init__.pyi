def write_vtp(
    nodes: list[list[float]],
    filename: str | None,
) -> None:
    """
    Write a ASCII VTP poly file for nodes.

    Parameters
    ----------
    nodes : list[list[float]]
        The list of 3D points to write.
    filename : str | None
        Optional filename. Defaults to grid.vtp
    """

def write_stl(
    nodes: list[list[float]],
    cells: list[list[int]],
    filename: str | None,
) -> None:
    """
    Write a ASCII stl file from nodes and cells.

    Parameters
    ----------
    nodes : list[list[float]]
        The list of 3D points to write.
    cells : list[list[int]]
        The list of cell indices of connecting nodes.
    filename : str | None
        Optional filename. Defaults to mesh.stl
    """
