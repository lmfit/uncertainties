from . base import AffineScalarFuncBase, LinearCombination
from . formatting import UFloatFormatting
from . ops import AffineScalarFuncOps

class AffineScalarFunc(
        AffineScalarFuncOps, UFloatFormatting
        ):
    pass

AffineScalarFunc._add_arithmetic_ops()
AffineScalarFunc._add_comparative_ops()