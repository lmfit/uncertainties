from __future__ import annotations

from collections import defaultdict
from math import sqrt
from typing import Optional, Tuple, Union
import uuid


class UAtom:
    __slots__ = ["uuid", "tag", "_hash"]

    def __init__(self: UAtom, tag: Optional[str] = None):
        self.tag = tag
        self.uuid: uuid.UUID = uuid.uuid4()
        self._hash = hash(self.uuid)

    def __eq__(self: UAtom, other: UAtom) -> bool:
        return hash(self) == hash(other)

    def __hash__(self: UAtom) -> int:
        return self._hash

    def __lt__(self: UAtom, other: UAtom) -> bool:
        return self.uuid < other.uuid

    def __str__(self: UAtom) -> str:
        uuid_abbrev = f"{str(self.uuid)[0:2]}..{str(self.uuid)[-2:]}"
        if self.tag is not None:
            label = f"{self.tag}, {uuid_abbrev}"
        else:
            label = uuid_abbrev
        return f"{self.__class__.__name__}({label})"


class UCombo:
    __slots__ = ["ucombo_tuple", "is_expanded", "_std_dev", "_expanded_dict", "_hash"]

    def __init__(
        self: UCombo,
        ucombo_tuple: Tuple[Tuple[Union[UAtom, UCombo], float], ...],
    ):
        self.ucombo_tuple = ucombo_tuple
        self.is_expanded = False
        self._std_dev = None
        self._expanded_dict = None
        self._hash = None

    @property
    def expanded(self: UCombo) -> dict[UAtom, float]:
        if self._expanded_dict is None:
            term_list = list(self.ucombo_tuple)
            self._expanded_dict = defaultdict(float)
            while term_list:
                term, weight = term_list.pop()
                if isinstance(term, UAtom):
                    self._expanded_dict[term] += weight
                elif term.is_expanded:
                    for sub_term, sub_weight in term.expanded.items():
                        self._expanded_dict[sub_term] += weight * sub_weight
                else:
                    for sub_term, sub_weight in term.ucombo_tuple:
                        term_list.append((sub_term, weight * sub_weight))
            self.is_expanded = True
        return self._expanded_dict

    @property
    def std_dev(self: UCombo) -> float:
        if self._std_dev is None:
            self._std_dev = sqrt(sum(weight**2 for weight in self.expanded.values()))
        return self._std_dev

    def covariance(self: UCombo, other: UCombo):
        # TODO: pull out to module function and cache
        self_uatoms = set(self.expanded.keys())
        other_uatoms = set(other.expanded.keys())
        shared_uatoms = self_uatoms.intersection(other_uatoms)
        covariance = 0
        for uatom in shared_uatoms:
            covariance += self.expanded[uatom] * other.expanded[uatom]
        return covariance

    def __add__(self: UCombo, other) -> UCombo:
        if not isinstance(other, UCombo):
            return NotImplemented
        if not other:
            return self
        if not self:
            return other
        return UCombo(((self, 1.0), (other, 1.0)))

    def __radd__(self: UCombo, other: UCombo) -> UCombo:
        return self.__add__(other)

    def __mul__(self: UCombo, scalar: Union[float, int]) -> UCombo:
        if not isinstance(scalar, (float, int)):
            return NotImplemented
        if scalar == 0 or not self:
            return UCombo(())
        return UCombo(((self, float(scalar)),))

    def __rmul__(self: UCombo, scalar: Union[float, int]) -> UCombo:
        return self.__mul__(scalar)

    def __iter__(self: UCombo):
        return iter(self.ucombo_tuple)

    def __bool__(self: UCombo) -> bool:
        return bool(self.ucombo_tuple)

    def __str__(self: UCombo) -> str:
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

    def __hash__(self: UCombo) -> int:
        if self._hash is None:
            self._hash = hash(
                tuple(sorted((key, value) for key, value in self.expanded.items()))
            )
        return self._hash

    def __eq__(self: UCombo, other: UCombo) -> bool:
        return hash(self) == hash(other)
