from __future__ import annotations

from collections import defaultdict
from functools import cached_property
from dataclasses import dataclass, field
from math import sqrt
from numbers import Real
from typing import Dict, Optional, Tuple, TypeVar, Union
import uuid


@dataclass(frozen=True)
class UAtom:
    uuid: uuid.UUID = field(init=False, default_factory=uuid.uuid4)
    tag: Optional[str] = None

    def __str__(self):
        uuid_abbrev = f"{str(self.uuid)[0:2]}..{str(self.uuid)[-3:-1]}"
        if self.tag is not None:
            label = f"{self.tag}, {uuid_abbrev}"
        else:
            label = uuid_abbrev
        return f"{self.__class__.__name__}({label})"


Self = TypeVar("Self", bound="UCombo")  # TODO: typing.Self introduced in Python 3.11


# TODO: Right now UCombo lacks __slots__. Python 3.10 allows slot=True input argument to
#   dataclass. Until then the easiest way to get __slots__ back would be to not use a
#   dataclass here.
@dataclass(frozen=True)
class UCombo:
    ucombo_tuple: Tuple[Tuple[Union[UAtom, UCombo], float], ...]

    # TODO: Using cached_property instead of lru_cache. This misses the opportunity to
    #   cache across separate instances.
    @cached_property
    def expanded_dict(self: Self) -> Dict[UAtom, float]:
        expanded_dict: Dict[UAtom, float] = defaultdict(float)

        for term, term_weight in self:
            if isinstance(term, UAtom):
                expanded_dict[term] += term_weight
            else:
                expanded_term = term.expanded_dict
                for atom, atom_weight in expanded_term.items():
                    expanded_dict[atom] += term_weight * atom_weight

        # pruned_expanded_dict = defaultdict(float)
        # pruned_expanded_dict.update({
        #     atom: float(weight) for atom, weight in expanded_dict.items() if weight != 0
        # })
        # # pruned_expanded_dict = defaultdict(float)

        return expanded_dict

    @property
    def expanded(self: Self) -> Dict[UAtom, float]:
        return self.expanded_dict

    @cached_property
    def std_dev(self: Self) -> float:
        return sqrt(sum(weight**2 for weight in self.expanded_dict.values()))

    def __add__(self: Self, other) -> Self:
        if not isinstance(other, UCombo):
            return NotImplemented
        return UCombo(((self, 1.0), (other, 1.0)))

    def __radd__(self: Self, other):
        return self.__add__(other)

    def __mul__(self: Self, scalar: Real):
        if not isinstance(scalar, Real):
            return NotImplemented
        return UCombo(((self, float(scalar)),))

    def __rmul__(self: Self, scalar: Real):
        return self.__mul__(scalar)

    def __iter__(self: Self):
        return iter(self.ucombo_tuple)

    def __str__(self: Self) -> str:
        ret_str = ""
        first = True
        for term, weight in self:
            if not first:
                ret_str += " + "
            else:
                first = False

            if isinstance(term, UAtom):
                ret_str += f"{weight}×{str(term)}"
            else:
                ret_str += f"{weight}×({str(term)})"
        return ret_str
