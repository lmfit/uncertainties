import timeit

import pytest

from uncertainties import ufloat


def ufloat_sum_benchmark(num):
    sum(ufloat(1, 1) for _ in range(num)).std_dev


def time_ufloat_sum_benchmark(num):
    reps = int(100000 / num)
    t = timeit.timeit(lambda: ufloat_sum_benchmark(num), number=reps)
    return t / reps


def test_algorithm_complexity():
    result_dict = {}
    n_list = (10, 100, 1000, 10000, 100000)
    for n in n_list:
        result_dict[n] = time_ufloat_sum_benchmark(n)
    n0 = n_list[0]
    t0 = result_dict[n0]
    for n, t in result_dict.items():
        assert 1 / 4 < (t / n) / (t0 / n0) < 4
    assert result_dict[100000] < 0.75 * 4


@pytest.mark.benchmark
def test_speed():
    time_ufloat_sum_benchmark(100000)
