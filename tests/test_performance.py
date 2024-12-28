from math import log10
import time
import timeit

import pytest

from uncertainties import ufloat


def repeated_summation(num):
    """
    generate and sum many floats together, then calculate the standard deviation of the
    output. Under the lazy expansion algorithm, the uncertainty remains non-expanded
    until a request is made to calculate the standard deviation.
    """
    result = sum(ufloat(1, 0.1) for _ in range(num)).std_dev
    return result


def test_repeated_summation_complexity():
    """
    Test that the execution time is linear in summation length
    """
    approx_execution_time_per_n = 10e-6  # 10 us
    target_test_duration = 1  # 1 s

    n_list = [10, 100, 1000, 10000, 100000]
    t_list = []
    for n in n_list:
        """
        Choose the number of repetitions so that the test takes target_test_duration
        assuming the timing of a single run is approximately
        N * approx_execution_time_per_n
        """
        # Choose the number of repetitions so that the test
        single_rep_duration = n * approx_execution_time_per_n
        num_reps = int(target_test_duration / single_rep_duration)

        t_tot = timeit.timeit(
            lambda: repeated_summation(n),
            number=num_reps,
            timer=time.process_time,
        )
        t_single = t_tot / num_reps
        t_list.append(t_single)
    n0 = n_list[0]
    t0 = t_list[0]
    for n, t in zip(n_list[1:], t_list[1:]):
        # Check that the plot of t vs n is linear on a log scale to within 10%
        # See PR 275
        assert 0.9 * log10(n / n0) < log10(t / t0) < 1.1 * log10(n / n0)


@pytest.mark.parametrize("num", (10, 100, 1000, 10000, 100000))
@pytest.mark.benchmark
def test_repeated_summation_speed(num):
    repeated_summation(num)
