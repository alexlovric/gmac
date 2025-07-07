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
    format: str = "binary",
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
    format : str
        Format of the stl file. Defaults to "binary".
    """

def read_stl(filename: str) -> tuple[list[list[float]], list[list[int]]]:
    """
    Read a stl file and return nodes and cells.

    Parameters
    ----------
    filename : str
        The filename of the stl file.

    Returns
    -------
    tuple[list[list[float]], list[list[int]]]
        The nodes and cells of the stl file.
    """

def read_obj(filename: str) -> tuple[list[list[float]], list[list[int]]]:
    """
    Read a obj file and return nodes and cells.

    Parameters
    ----------
    filename : str
        The filename of the obj file.

    Returns
    -------
    tuple[list[list[float]], list[list[int]]]
        The nodes and cells of the obj file.
    """
