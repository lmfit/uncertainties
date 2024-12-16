import platform
import timeit

import psutil

from uncertainties import ufloat


def print_system_info():
    sys_info = {
        "platform": platform.system(),
        "platform-release": platform.release(),
        "platform-version": platform.version(),
        "architecture": platform.machine(),
        "processor": platform.processor(),
        "ram": str(round(psutil.virtual_memory().total / (1024.0**3))) + " GB",
    }

    for key, value in sys_info.items():
        print(f"{key:17}: {value}")


def ufloat_sum_benchmark(num):
    str(sum(ufloat(1, 1) for _ in range(num)))


def test_algorithm_performance():
    print_system_info()
    result_dict = {}
    n_list = (10, 100, 1000, 10000, 100000)
    for n in n_list:
        print(f"### {n=} ###")
        reps = int(100000 / n)
        t = timeit.timeit(lambda: ufloat_sum_benchmark(n), number=reps)
        t_avg = t / reps
        print(f"    Test duration: {t:.2f} s, Repetitions: {reps}")
        print(f"    Average execution time: {t_avg:.4f} s")
        result_dict[n] = t_avg
    n0 = n_list[0]
    t0 = result_dict[n0]
    for n, t in result_dict.items():
        assert 1 / 4 < (t / n) / (t0 / n0) < 4
    assert result_dict[100000] < 0.75 * 4
