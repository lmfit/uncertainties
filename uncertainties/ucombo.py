from __future__ import annotations

from collections import defaultdict
from math import sqrt
from numbers import Real
from typing import Tuple, TypeVar, Union
import uuid


class UAtom:
    __slots__ = ["uuid", "tag", "hash"]

    def __init__(self, tag: str = None):
        self.tag = tag
        self.uuid: uuid.UUID = uuid.uuid4()
        self.hash = hash(self.uuid)  # memoize the hash

    def __eq__(self, other):
        return self.hash == other.hash

    def __hash__(self):
        return self.hash

    def __str__(self):
        uuid_abbrev = f"{str(self.uuid)[0:2]}..{str(self.uuid)[-3:-1]}"
        if self.tag is not None:
            label = f"{self.tag}, {uuid_abbrev}"
        else:
            label = uuid_abbrev
        return f"{self.__class__.__name__}({label})"


Self = TypeVar("Self", bound="UCombo")  # TODO: typing.Self introduced in Python 3.11


class UCombo:
    __slots__ = ["ucombo_tuple", "_std_dev", "_expanded_dict"]

    def __init__(self, ucombo_tuple: Tuple[Tuple[Union[UAtom, UCombo], float], ...]):
        self.ucombo_tuple = ucombo_tuple
        self._std_dev = None
        self._expanded_dict = None

    @property
    def expanded_dict(self: Self) -> dict[UAtom, float]:
        if self._expanded_dict is None:
            term_list = list(self.ucombo_tuple)
            self._expanded_dict = defaultdict(float)
            while term_list:
                term, weight = term_list.pop()
                if isinstance(term, UAtom):
                    self._expanded_dict[term] += weight
                elif term.expanded_dict is not None:
                    for sub_term, sub_weight in term.expanded_dict.items():
                        self._expanded_dict[sub_term] += weight * sub_weight
                else:
                    for sub_term, sub_weight in term.ucombo_tuple:
                        term_list.append((sub_term, weight * sub_weight))
        return self._expanded_dict

    @property
    def expanded(self: Self) -> dict[UAtom, float]:
        return self.expanded_dict

    @property
    def std_dev(self: Self) -> float:
        if self._std_dev is None:
            self._std_dev = sqrt(
                sum(weight**2 for weight in self.expanded_dict.values())
            )
        return self._std_dev

    def covariance(self: Self, other: UCombo):
        # TODO: pull out to module function and cache
        self_uatoms = set(self.expanded_dict.keys())
        other_uatoms = set(other.expanded_dict.keys())
        shared_uatoms = self_uatoms.intersection(other_uatoms)
        covariance = 0
        for uatom in shared_uatoms:
            covariance += self.expanded_dict[uatom] * other.expanded_dict[uatom]
        return covariance

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

    def __bool__(self):
        return bool(self.ucombo_tuple)

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
