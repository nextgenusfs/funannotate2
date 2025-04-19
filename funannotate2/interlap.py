from itertools import groupby
from operator import itemgetter

__all__ = ["InterLap"]

__version__ = "0.2.6"

try:
    int_types = (int, int)
except NameError:
    int_types = (int,)


def binsearch_left_start(intervals, x, lo, hi):
    while lo < hi:
        mid = (lo + hi) // 2
        f = intervals[mid]
        if f[0] < x:
            lo = mid + 1
        else:
            hi = mid
    return lo


# like python's bisect_right find the _highest_ index where the value x
# could be inserted to maintain order in the list intervals


def binsearch_right_end(intervals, x, lo, hi):
    while lo < hi:
        mid = (lo + hi) // 2
        f = intervals[mid]
        if x < f[0]:
            hi = mid
        else:
            lo = mid + 1
    return lo


class InterLap(object):
    """Create an Interlap object. (See module docstring)."""

    def __init__(self, ranges=()):
        """The ranges are tuples of (start, end, *whatever)."""
        self._iset = sorted(ranges)
        self._maxlen = max(r[1] - r[0] for r in (ranges or [[0, 0]]))

    def add(self, ranges):
        r"""Add a single (or many) [start, end, \*] item to the tree."""
        if len(ranges) and isinstance(ranges[0], int_types):
            ranges = [ranges]
        iset = self._iset
        self._maxlen = max(self._maxlen, max(r[1] - r[0] + 1 for r in ranges))

        if len(ranges) > 30 or len(iset) < len(ranges):
            iset.extend(ranges)
            iset.sort()
        else:
            for o in ranges:
                iset.insert(binsearch_left_start(iset, o[0], 0, len(iset)), o)

    update = add

    def __len__(self):
        """Return number of intervals."""
        return len(self._iset)

    def find(self, other):
        """Return an interable of elements that overlap other in the tree."""
        iset = self._iset
        left_idx = binsearch_left_start(iset, other[0] - self._maxlen, 0, len(iset))
        right_idx = binsearch_right_end(iset, other[1], 0, len(iset))
        iopts = iset[left_idx:right_idx]
        iiter = (s for s in iopts if s[0] <= other[1] and s[1] >= other[0])
        for o in iiter:
            yield o

    def closest(self, other):
        iset = self._iset
        left_idx = binsearch_left_start(iset, other[0] - self._maxlen, 0, len(iset))
        right_idx = binsearch_right_end(iset, other[1], 0, len(iset))
        left_idx, right_idx = max(0, left_idx - 1), min(len(iset), right_idx + 2)

        while right_idx < len(iset) and iset[right_idx - 1][0] == iset[right_idx][0]:
            right_idx += 1

        while left_idx > 1 and iset[left_idx][1] == iset[left_idx + 1][1]:
            left_idx -= 1
        iopts = iset[left_idx:right_idx]
        ovls = [s for s in iopts if s[0] <= other[1] and s[1] >= other[0]]
        if ovls:
            for o in ovls:
                yield o
        else:
            iopts = sorted([(min(abs(i[0] - other[1]), abs(other[0] - i[1])), i) for i in iopts])
            for _, g in groupby(iopts, itemgetter(0)):
                # only yield the closest intervals
                for _, ival in g:
                    yield ival
                break

    def __contains__(self, other):
        """Indicate whether `other` overlaps any elements in the tree."""
        iset = self._iset
        left_idx = binsearch_left_start(iset, other[0] - self._maxlen, 0, len(iset))
        # since often the found interval will overlap, we short cut that
        # case of speed.
        max_search = 8
        if left_idx == len(iset):
            return False
        for left in iset[left_idx : left_idx + max_search]:
            if left[1] >= other[0] and left[0] <= other[1]:
                return True
            if left[0] > other[1]:
                return False

        right_idx = binsearch_right_end(iset, other[1], 0, len(iset))
        return any(
            s[0] <= other[1] and s[1] >= other[0] for s in iset[left_idx + max_search : right_idx]
        )

    def __iter__(self):
        return iter(self._iset)


def overlaps(s1, e1, s2, e2):
    """
    >>> overlaps(2, 4, 3, 5)
    True
    >>> overlaps(2, 4, 4, 5)
    False
    >>> overlaps(2, 200, 3, 5)
    True
    >>> overlaps(3, 5, 2, 200)
    True
    >>> overlaps (3, 5, 5, 6)
    False
    >>> overlaps (5, 6, 3, 5)
    False
    >>> overlaps(2, 4, 2, 4)
    True
    """
    return not (e1 <= s2 or s1 >= e2)


def reduce(args):
    """
    >>> reduce([(2, 4), (4, 9)])
    [(2, 4), (4, 9)]

    >>> reduce([(2, 6), (4, 10)])
    [(2, 10)]
    """
    if len(args) < 2:
        return args
    args.sort()
    ret = [args[0]]
    for next_i, (_, e) in enumerate(args, start=1):
        if next_i == len(args):
            ret[-1] = ret[-1][0], max(ret[-1][1], e)
            break

        ns, ne = args[next_i]
        if e > ns or ret[-1][1] > ns:
            ret[-1] = ret[-1][0], max(e, ne, ret[-1][1])
        else:
            ret.append((ns, ne))
    return ret


class Interval(object):
    """
    >>> i = Interval([(2, 10), (8, 20), (30, 40)])
    >>> i
    Interval([(2, 20), (30, 40)])

    >>> i.add([(20, 22)])
    >>> i
    Interval([(2, 20), (20, 22), (30, 40)])

    >>> i.add([(10, 31)])
    >>> i
    Interval([(2, 40)])

    >>> i.add([(55, 65), (65, 75), (75, 85), (85, 95), (95, 100)])
    >>> i
    Interval([(2, 40), (55, 65), (65, 75), (75, 85), (85, 95), (95, 100)])

    >>> i.add([(1, 95)])
    >>> i
    Interval([(1, 95), (95, 100)])

    >>> i.add(Interval([(90, 100)]))
    >>> i
    Interval([(1, 100)])

    >>> Interval()
    Interval([])

    """

    __slots__ = ("_vals", "_fixed")

    def __init__(self, args=None):
        self._vals = []
        if args is None:
            return
        assert isinstance(args, list)
        if len(args) > 0:
            assert isinstance(args[0], tuple), args
            assert isinstance(args[0][0], int)
            self._vals = reduce(args)

    def _as_tuples(self, args):
        vals = []
        if isinstance(args, Interval):
            vals.extend(args._vals)
        else:
            for a in args:
                if isinstance(a, Interval):
                    vals.extend(a._vals)
                else:
                    vals.append(a)
        return vals

    def add(self, args):
        self._vals = reduce(self._vals + self._as_tuples(args))

    def __repr__(self):
        return "%s(%r)" % (self.__class__.__name__, self._vals)


if __name__ == "__main__":
    import time

    t0 = time.time()
    import doctest

    print((doctest.testmod(verbose=0, optionflags=doctest.REPORT_ONLY_FIRST_FAILURE)))
    print((time.time() - t0))
