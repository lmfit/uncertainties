import pytest
from unittest.mock import patch


orig_import = __import__


def no_numpy_import(*args, **kwargs):
    print(args[0])
    if args and args[0].startswith("numpy"):
        raise ImportError
    elif "name" in kwargs and kwargs["name"].startswith("numpy"):
        raise ImportError
    else:
        return orig_import(*args, **kwargs)


with patch("builtins.__import__", no_numpy_import):
    from uncertainties import (
        correlated_values,
        correlated_values_norm,
        correlation_matrix,
        ufloat,
    )


def test_no_numpy():
    nom_values = [1, 2, 3]
    std_devs = [0.1, 0.2, 0.3]
    cov = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ]

    with pytest.raises(
        NotImplementedError,
        match="not able to import numpy",
    ):
        _ = correlated_values(nom_values, cov)

    with pytest.raises(
        NotImplementedError,
        match="not able to import numpy",
    ):
        _ = correlated_values_norm(
            list(zip(nom_values, std_devs)),
            cov,
        )

    x = ufloat(1, 0.1)
    y = ufloat(2, 0.2)
    z = ufloat(3, 0.3)

    with pytest.raises(
        NotImplementedError,
        match="not able to import numpy",
    ):
        _ = correlation_matrix([x, y, z])
