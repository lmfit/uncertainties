import json
from pathlib import Path
import random


real_single_input_funcs = [
    "asinh",
    "atan",
    "cos",
    "cosh",
    "degrees",
    "erf",
    "erfc",
    "exp",
    "expm1",
    "gamma",
    "lgamma",
    "radians",
    "sin",
    "sinh",
    "tan",
    "tanh",
]
positive_single_input_funcs = [
    "log",
    "log1p",
    "log10",
    "sqrt",
]
minus_one_to_plus_one_single_input_funcs = [
    "acos",
    "asin",
    "atanh",
]
greater_than_one_single_input_funcs = ["acosh"]
single_input_funcs = (
    real_single_input_funcs
    + positive_single_input_funcs
    + minus_one_to_plus_one_single_input_funcs
    + greater_than_one_single_input_funcs
)

real_double_input_funcs = (
    "atan2",
    "pow",  # We will coerce the base to be positive
)
positive_double_input_funcs = (
    "hypot",
    "log",
)
double_input_funcs = real_double_input_funcs + positive_double_input_funcs

umath_function_cases_dict = {}

for func_name in single_input_funcs:
    cases_list = []
    for _ in range(10):
        if func_name in real_single_input_funcs:
            min_val = -100
            max_val = +100
        elif func_name in positive_single_input_funcs:
            min_val = 0
            max_val = +100
        elif func_name in minus_one_to_plus_one_single_input_funcs:
            min_val = -1
            max_val = +1
        elif func_name in greater_than_one_single_input_funcs:
            min_val = +1
            max_val = +100
        else:
            raise ValueError
        nominal_value = random.uniform(min_val, max_val)
        std_dev = random.uniform(0, 100)
        ufloat_tuple = (nominal_value, std_dev)
        case = (ufloat_tuple,)
        cases_list.append(case)
    umath_function_cases_dict[func_name] = cases_list

for func_name in double_input_funcs:
    cases_list = []
    for _ in range(10):
        if func_name in real_double_input_funcs:
            min_val = -100
        elif func_name in positive_double_input_funcs:
            min_val = 0
        else:
            raise ValueError
        max_val = +100

        nominal_value_0 = random.uniform(min_val, max_val)
        if func_name == "pow":
            nominal_value_0 = abs(nominal_value_0)
        std_dev_0 = random.uniform(0, 100)
        ufloat_tuple_0 = (nominal_value_0, std_dev_0)

        nominal_value_1 = random.uniform(min_val, max_val)
        std_dev_1 = random.uniform(0, 100)
        ufloat_tuple_1 = (nominal_value_1, std_dev_1)

        case = (ufloat_tuple_0, ufloat_tuple_1)
        cases_list.append(case)
    umath_function_cases_dict[func_name] = cases_list

umath_function_cases_json_path = Path(Path.cwd(), "umath_function_cases.json")

if __name__ == "__main__":
    with open(umath_function_cases_json_path, "w") as f:
        json.dump(umath_function_cases_dict, f, indent=4)
