import json
from pathlib import Path
import random


unary_functions = ["__neg__", "__pos__"]
binary_functions = [
    "__add__",
    "__radd__",
    "__sub__",
    "__rsub__",
    "__mul__",
    "__rmul__",
    "__truediv__",  # Will fail if second value is 0 but that is unlikely
    "__rtruediv__",  # Will fail if first value is 0 but that is unlikely
    "__pow__",  # We will coerce the base to be positive
    "__rpow__",  # We will coerce the base to be positive
]

ufloat_method_cases_dict = {}

for func_name in unary_functions:
    cases_list = []
    for _ in range(10):
        nominal_value = random.uniform(-100, 100)
        std_dev = random.uniform(0, 100)
        ufloat_tuple = (nominal_value, std_dev)
        case = (ufloat_tuple,)
        cases_list.append(case)
    ufloat_method_cases_dict[func_name] = cases_list

for func_name in binary_functions:
    cases_list = []
    for _ in range(10):
        nominal_value_0 = random.uniform(-100, 100)
        if func_name == "__pow__":
            nominal_value_0 = abs(nominal_value_0)
        std_dev_0 = random.uniform(0, 100)
        ufloat_tuple_0 = (nominal_value_0, std_dev_0)

        nominal_value_1 = random.uniform(-100, 100)
        if func_name == "__rpow__":
            nominal_value_1 = abs(nominal_value_1)
        std_dev_1 = random.uniform(0, 100)
        ufloat_tuple_1 = (nominal_value_1, std_dev_1)

        case = (ufloat_tuple_0, ufloat_tuple_1)
        cases_list.append(case)
    ufloat_method_cases_dict[func_name] = cases_list

ufloat_method_cases_json_path = Path(Path.cwd(), "ufloat_method_cases.json")

if __name__ == "__main__":
    with open(ufloat_method_cases_json_path, "w") as f:
        json.dump(ufloat_method_cases_dict, f, indent=4)
