# Abstractions of numeric types to finite unions of closed intervals
# and abstractions of Boolean to tristate (namely for comparisons of intervals)

import numbers
from itertools import product
from enum import Enum


# Tristate implements an abstract domain for booleans:
# TRUE, FALSE, or UNKNOWN
#
# Note that there is no notion of overloading 'and', 'or', and 'not' in Python.
# They can still be used soundly, but a ValueError is raised
# if bool(UNKNOWN) is invoked.
# So we implement '&', '|', '-' as well, for convenience --
# but note these operators do not behave "correctly" on actual bools
class Tristate(Enum):
    TRUE = True
    FALSE = False
    UNKNOWN = None

    def __bool__(self):
        if self == Tristate.TRUE:
            return True
        elif self == Tristate.FALSE:
            return False
        else:
            raise ValueError("Can't cast " + str(self) + " to bool")

    def __and__(self, other):
        if isinstance(other, bool) or other is None:
            return self.__and__(Tristate(other))
        elif isinstance(other, Tristate):
            try:
                return self and other
            except ValueError: # self is Tristate.UNKNOWN
                return Tristate.FALSE if other == Tristate.FALSE else Tristate.UNKNOWN
        else:
            return NotImplemented
    def __rand__(self, other):
        return self.__and__(other)

    def __or__(self, other):
        if isinstance(other, bool) or other is None:
            return self.__or__(Tristate(other))
        elif isinstance(other, Tristate):
            try:
                return self or other
            except ValueError: # self is Tristate.UNKNOWN
                return Tristate.TRUE if other == Tristate.TRUE else Tristate.UNKNOWN
        else:
            return NotImplemented
    def __ror__(self, other):
        return self.__or__(other)

    def __neg__(self):
        assert isinstance(self, Tristate)
        if self == Tristate.TRUE:
            return Tristate.FALSE
        elif self == Tristate.FALSE:
            return Tristate.TRUE
        else:
            return Tristate.UNKNOWN


class Interval(object):

    def __init__(self, *pairs):
        assert len(pairs) > 0
        for p in pairs:
            assert p[0] <= p[1]
        self._pairs = Interval._normalize(list(pairs))

    @property
    def pairs(self):
        return self._pairs
    # Mutating the pairs should call _normalize
    @pairs.setter
    def pairs(self, val):
        if not (isinstance(val, list) and len(val) > 0):
            raise ValueError
        self._pairs = Interval._normalize(val.copy())

    # Properties for lower/upper bounds
    @property
    def upper(self):
        return self._pairs[-1][1]
    @property
    def lower(self):
        return self._pairs[0][0]

    # Collapse a list of intervals into a minimal normal form
    # Note: mutates the input (and returns it too)
    def _normalize(pairs):
        pairs.sort()
        i = 0
        while i < len(pairs):
            while i + 1 < len(pairs) and pairs[i+1][0] <= pairs[i][1]:
                pairs[i] = (pairs[i][0], max(pairs[i][1], pairs[i+1][1]))
                del pairs[i+1]
            i += 1
        return pairs

    def __str__(self):
        return str(self.pairs)

    def __contains__(self, item):
        if not isinstance(item, numbers.Number):
            raise ValueError()
        for pair in self.pairs:
            if pair[0] <= item and item <= pair[1]:
                return True
        return False

    # Arithmetic emulation
    
    def __add__(self, other):
        if isinstance(other, Interval):
            ret = []
            for i,j in product(self.pairs, other.pairs):
                ret.append((i[0]+j[0], i[1]+j[1]))
            return Interval(*ret)
        elif isinstance(other, numbers.Number):
            return self.__add__(Interval((other,other)))
        else:
            return NotImplemented
    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self.__add__(-other)
    def __rsub__(self, other):
        temp = -self
        return temp.__add__(other)

    def __mul__(self, other):
        if isinstance(other, Interval):
            ret = []
            for i,j in product(self.pairs, other.pairs):
                vals = [i[0]*j[0], i[0]*j[1], i[1]*j[0], i[1]*j[1]]
                lb = min(vals)
                ub = max(vals)
                ret.append((lb,ub))
            return Interval(*ret)
        elif isinstance(other, numbers.Number):
            return self.__mul__(Interval((other,other)))
        else:
            return NotImplemented
    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, Interval):
            if 0 in other:
                raise ZeroDivisionError() # TODO make this a warning a support numeric extensions to infinity
            ret = []
            for i,j in product(self.pairs, other.pairs):
                vals = [i[0]/j[0], i[0]/j[1], i[1]/j[0], i[1]/j[1]]
                lb = min(vals)
                ub = max(vals)
                ret.append((lb,ub))
            return Interval(*ret)
        elif isinstance(other, numbers.Number):
            return self.__truediv__(Interval((other,other)))
        else:
            return NotImplemented
    def __rtruediv__(self, other):
        # isinstance(other, Interval) is always False
        if isinstance(other, numbers.Number):
            return Interval((other,other)).__truediv__(self)
        else:
            return NotImplemented

    def __neg__(self):
        inv = [(-p[1], -p[0]) for p in self.pairs]
        return Interval(*inv)

    # Comparisons -- these all return Tristate objects, not bool
    # and are comparing possibilities for the concrete values;
    # NOT partially ordering the abstract values

    def __lt__(self, other):
        if isinstance(other, Interval):
            if self.pairs[-1][1] < other.pairs[0][0]:
                return Tristate.TRUE
            elif other.pairs[-1][1] <= self.pairs[0][0]:
                return Tristate.FALSE
            else:
                return Tristate.UNKNOWN
        elif isinstance(other, numbers.Number):
            return self.__lt__(Interval((other,other)))
        else:
            return NotImplemented

    def __le__(self, other):
        if isinstance(other, Interval):
            if self.pairs[-1][1] <= other.pairs[0][0]:
                return Tristate.TRUE
            elif other.pairs[-1][1] < self.pairs[0][0]:
                return Tristate.FALSE
            else:
                return Tristate.UNKNOWN
        elif isinstance(other, numbers.Number):
            return self.__le__(Interval((other,other)))
        else:
            return NotImplemented

    # Note: this checks if the Interval objects have the same concrete value
    # To check if two Interval objects have the same abstract value,
    # check self.pairs == other.pairs
    def __eq__(self, other):
        if isinstance(other, Interval):
            # If they are both singletons, equality check
            if len(self.pairs) == 1 and len(other.pairs) == 1:
                if self.pairs[0][0] == self.pairs[0][1] and other.pairs[0][0] == other.pairs[0][1]:
                    return Tristate(self.pairs[0][0] == other.pairs[0][0])
            # If there is overlap, then unknown
            for i,j in product(self.pairs, other.pairs):
                if i[1] < j[0] or j[1] < i[0]:
                    continue
                else:
                    return Tristate.UNKNOWN
            # If there is no overlap, then false
            return Tristate.FALSE
        elif isinstance(other, numbers.Number):
            return self.__eq__(Interval((other,other)))
        else:
            return NotImplemented

    def __ne__(self, other):
        inv = self.__eq__(other)
        if inv == NotImplemented:
            return NotImplemented
        return Tristate.__neg__(inv)

    def __gt__(self, other):
        if isinstance(other, Interval):
            return other.__lt__(self)
        elif isinstance(other, numbers.Number):
            return self.__gt__(Interval((other,other)))
        else:
            return NotImplemented

    def __ge__(self, other):
        if isinstance(other, Interval):
            return other.__le__(self)
        elif isinstance(other, numbers.Number):
            return self.__ge__(Interval((other,other)))
        else:
            return NotImplemented


