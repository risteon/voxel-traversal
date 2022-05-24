import numpy as np
import pytraversal


def test_create_grid_without_args():
    _ = pytraversal.Grid3D()


def test_single_traversal():
    # define grid
    grid = pytraversal.Grid3D([0., 0., -50.], [4., 2., -15.], [4, 2, 1])
    # define expectation
    expected = np.asarray([[1, 0, 0], [2, 0, 0], [2, 1, 0], [3, 1, 0]], np.int64)

    traversed = grid.traverse([1.5, 0.8, -25.0], [6.0, 1.7, -25.0])
    np.testing.assert_array_equal(expected, traversed)

    traversed = grid.traverse([6.0, 1.7, -25.0], [1.5, 0.8, -25.0])
    np.testing.assert_array_equal(expected[::-1], traversed)
