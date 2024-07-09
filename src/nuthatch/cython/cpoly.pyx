import cython
def c_mon_mul(self, other):
    if self.packed and other.packed:
        s: cython.int = self.weight
        o: cython.int = other.weight
        return other.__class__(s + o)
    new_list = []
    self_index: cython.int = 0
    other_index: cython.int = 0
    self_len: cython.int = len(self.degrees)
    other_len: cython.int = len(other.degrees)
    while self_index < self_len and other_index < other_len:
        if self.degrees[self_index] < other.degrees[other_index]:
            new_list.extend(self.degrees[self_index : self_index + 2])
            self_index += 2
        elif self.degrees[self_index] == other.degrees[other_index]:
            x = self.degrees[self_index + 1] + other.degrees[other_index + 1]
            if x != 0:
                new_list.append(self.degrees[self_index])
                new_list.append(x)
            self_index += 2
            other_index += 2
        else:
            new_list.extend(other.degrees[other_index : other_index + 2])
            other_index += 2

    # Now that we've reached the end of at least one list, we add
    # everything from the other list.
    while self_index < self_len:
        new_list.extend(self.degrees[self_index : self_index + 2])
        self_index += 2
    while other_index < other_len:
        new_list.extend(other.degrees[other_index : other_index + 2])
        other_index += 2


def c_poly_mul(self, other):
    new_dict = {}
    for m, c in self.monomial_dictionary.items():
        for n, d in other.monomial_dictionary.items():
            k = c_mon_mul(m, n)
            e = c * d
            if k in new_dict:
                new_dict[k] += e
            else:
                new_dict[k] = e
    return other.__class__(other.base_ring, new_dict)



