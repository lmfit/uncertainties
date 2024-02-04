from . base import AffineScalarFuncBase, LinearCombination
from . formatting import UFloatFormatting
from . ops import AffineScalarFuncOps
from . numpy import UFloatNumpy

class AffineScalarFunc(
        AffineScalarFuncOps, UFloatFormatting, UFloatNumpy,
        ):
    pass
AffineScalarFunc._add_arithmetic_ops()
AffineScalarFunc._add_comparative_ops()
AffineScalarFunc._add_numpy_ufuncs()
AffineScalarFunc._add_numpy_comparative_ufuncs()