"""
C-MAT
Cython implementation of matrix multiplication methods. 
Compiled in C in the file cmat.c
"""
import cython
def c_generic_mul(se, oe, sr, on):
    new_entries = []
    for i in range(sr):
        new_entries.append([])
        for j in range(on):
            new_entry = se[i][0] * oe[0][j]
            for k in range(1, self.ncols):
                new_entry += se[i][k] * oe[k][j]
            new_entries[i].append(new_entry)
    return new_entries

