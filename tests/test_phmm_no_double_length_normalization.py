import inspect

from strpg.phmm import MotifAwarePHMM, _forward_array


def test_phmm_no_double_length_normalization():
    source = inspect.getsource(_forward_array)
    python_source = inspect.getsource(
        MotifAwarePHMM._forward_log_likelihood_python
    )
    assert "- math.log(float(m))" not in source
    assert "end_log = -math.log" not in python_source
