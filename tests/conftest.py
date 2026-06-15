import pytest
import dionysus as d


@pytest.fixture
def prime():
    return 5


@pytest.fixture(params=["matrix_v", "matrix_u"])
def reduction_method(request):
    return request.param


@pytest.fixture
def triangle_cells():
    return [
        d.Simplex([0], 0.0),
        d.Simplex([1], 1.0),
        d.Simplex([2], 2.0),
        d.Simplex([0, 1], 3.0),
        d.Simplex([0, 2], 4.0),
        d.Simplex([1, 2], 5.0),
        d.Simplex([0, 1, 2], 6.0),
    ]
