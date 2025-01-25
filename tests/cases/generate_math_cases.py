import json
import random


valid_inputs_dict = {}


def main():
    num_reps = 10

    single_inputs_dict = {
        "real": [random.uniform(-100, 100) for _ in range(num_reps)],
        "positive": [random.uniform(0, 100) for _ in range(num_reps)],
        "minus_one_to_plus_one": [random.uniform(-1, +1) for _ in range(num_reps)],
        "greater_than_one": [random.uniform(+1, 100) for _ in range(num_reps)],
    }
    with open("single_inputs.json", "w+") as f:
        json.dump(single_inputs_dict, f, indent=True)

    double_inputs_dict = {
        "real": [
            [random.uniform(-100, 100), random.uniform(-100, +100)]
            for _ in range(num_reps)
        ],
        "positive": [
            [random.uniform(0, 100), random.uniform(0, +100)] for _ in range(num_reps)
        ],
    }
    with open("double_inputs.json", "w+") as f:
        json.dump(double_inputs_dict, f, indent=True)


if __name__ == "__main__":
    main()
