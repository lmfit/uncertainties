from math import log10
import time
import timeit

import pytest

from uncertainties import ufloat


def ufloat_sum_benchmark(num):
    sum(ufloat(1, 0.1) for _ in range(num)).std_dev


def time_ufloat_sum_benchmark(num):
    # T ~ N * 10 us, so if we do 1001`000 / N repetitions the test takes ~ 1 s
    reps = int(100000 / num)
    t = timeit.timeit(
        lambda: ufloat_sum_benchmark(num),
        number=reps,
        timer=time.process_time,
    )
    return t / reps


def test_complexity():
    """
    Test that the execution time is linear in problem size
    """
    result_dict = {}
    n_list = (10, 100, 1000, 10000, 100000)
    for n in n_list:
        result_dict[n] = time_ufloat_sum_benchmark(n)
    n0 = n_list[0]
    t0 = result_dict[n0]
    for n, t in result_dict.items():
        if n == n0:
            continue
        # Check that the plot of t vs n is linear on a log scale to within 10%
        assert 0.9 * log10(n / n0) < log10(t / t0) < 1.1 * log10(n / n0)


@pytest.mark.benchmark
def test_speed():
    time_ufloat_sum_benchmark(100000)
