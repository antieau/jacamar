"""
MATRICES

A module for dense matrices. We create a single abstract matrix class which
holds a `data` attribute. We assume that this points to a class on which all
matrix operations can be performed and we overload matrices in that way. We
provide a generic implementation of these operations as `_MatrixGeneric`. Special
examples are provided by `FLINT`.

AUTHORS:
- Benjamin Antieau (2024): initial version.
- Luca Zerega (2024): major revisions.
"""

import flint
import numpy as np
import ray
from nuthatch.rings.integers import ZZ, ZZ_py
from nuthatch.rings.reals import RR, RR_py
from nuthatch.rings.complexes import CC
from nuthatch.rings.rationals import QQ
from nuthatch.rings.polynomials import PolynomialRing


class _MatrixGenericData:
    """
    This class must be called with all arguments. It is assumed, but not
    checked, that the shape of `entries` conforms with this. In the case that
    `nrows==0` or `ncols==0`, one should have `entries=[]`.
    """

    def __init__(self, *, base_ring, nrows, ncols, entries):
        self.base_ring = base_ring
        self.nrows = nrows
        self.ncols = ncols
        self.entries = entries

    def det(self):
        """Alias for `determinant` method."""
        return self.determinant()

    def determinant(self):
        if self.nrows != self.nrows:
            raise ValueError("Matrix must be square.")
        if self.nrows == 0:
            return self.base_ring(0)
        if self.nrows == 1:
            return self.entries[0][0]
        entries = self.entries
        def minor(array,i,j):
            c = array
            c = c[:i] + c[i+1:]
            for k in range(0,len(c)):
                c[k] = c[k][:j]+c[k][j+1:]
            return c
        def det(array,n):
            if n == 1 :
                return array[0][0]
            if n == 2 :
                return array[0][0]*array[1][1] - array[0][1]*array[1][0]
            det = PolynomialRing(base_ring=RR, ngens=3, prefix="x", packed=False)
            det = det.zero
            for i in range(0,n):
                m = minor(array,0,i)
                det = det + RR((-1)**i) * array[0][i] * det(m,n-1)
            return det
        return det(entries, self.nrows)

    def transpose(self):
        """Returns copy of the transposed matrix data."""
        entries = self.entries
        new_entries = []
        for i in range(self.ncols):
            new_entries.append([])
            for j in range(self.nrows):
                new_entries[-1].append(entries[j][i])

        return _MatrixGenericData(
            base_ring=self.base_ring,
            nrows=self.ncols,
            ncols=self.nrows,
            entries=new_entries
        )
    def size(self):
        """Returns size of a matrix as a tuple (rows, cols)."""
        return (self.nrows, self.ncols)

    def __add__(self, other):
        if self.nrows != other.nrows or self.ncols != other.ncols:
            raise ValueError(
                f"Cannot add matrix of size {self.nrows}x{self.ncols} to matrix of size {other.nrows}x{other.ncols}."
            )

        if self.nrows == 0 or other.nrows == 0:
            return other.__class__(
                base_ring=other.base_ring,
                nrows=other.nrows,
                ncols=other.ncols,
                entries=[],
            )
        new_entries = [self.entries[i].copy() for i in range(self.nrows)]
        for i in range(self.nrows):
            for j in range(self.ncols):
                new_entries[i][j] += other.entries[i][j]

        return other.__class__(
            base_ring=other.base_ring,
            nrows=other.nrows,
            ncols=other.ncols,
            entries=new_entries,
        )

    def __mul__(self, other):
        if str(type(other)) == "<class 'nuthatch.rings.reals.RealNumber'>" or str(type(other)) == "<class 'nuthatch.rings.complexes.ComplexNumber'>" or str(type(other)) == "<class 'nuthatch.rings.integers.Integer'>" or str(type(other)) == "<class 'nuthatch.rings.rationals.Rational'>":
            new_entries = []
            for i in self.entries:
                new_entries.append([])
                for j in i:
                    new_entries[-1].append(other * j)
                    
            return self.__class__(
                base_ring=self.base_ring,
                nrows=self.nrows, 
                ncols=self.ncols,
                entries=new_entries
            )
        if self.ncols != other.nrows:
            raise ValueError(
                f"Cannot multiply matrix of size {self.nrows}x{self.ncols} with matrix of size {other.nrows}x{other.ncols}."
            )

        if self.nrows == 0 or other.nrows == 0:
            
            return other.__class__(
                base_ring=other.base_ring,
                nrows=self.nrows,
                ncols=other.ncols,
                entries=[],
            )
        # KAUERS - MOOSBAUER
        if self.ncols % 5 == 0 and self.nrows % 5 == 0 and other.ncols % 5 == 0 and other.nrows % 5 == 0 and False ==  True:
            def kauers_moosbauer(a, b):
                la = a.nrows
                lb = b.ncols

                # Partition A and B into 25 blocks
                a11 = a[int(0*la/5):int(1*la/5), int(0*la/5):int(1*la/5)]
                b11 = b[int(0*lb/5):int(1*lb/5), int(0*lb/5):int(1*lb/5)]
                a12 = a[int(0*la/5):int(1*la/5), int(1*la/5):int(2*la/5)]
                b12 = b[int(0*lb/5):int(1*lb/5), int(1*lb/5):int(2*lb/5)]
                a13 = a[int(0*la/5):int(1*la/5), int(2*la/5):int(3*la/5)]
                b13 = b[int(0*lb/5):int(1*lb/5), int(2*lb/5):int(3*lb/5)]
                a14 = a[int(0*la/5):int(1*la/5), int(3*la/5):int(4*la/5)]
                b14 = b[int(0*lb/5):int(1*lb/5), int(3*lb/5):int(4*lb/5)]
                a15 = a[int(0*la/5):int(1*la/5), int(4*la/5):int(5*la/5)]
                b15 = b[int(0*lb/5):int(1*lb/5), int(4*lb/5):int(5*lb/5)]
                a21 = a[int(1*la/5):int(2*la/5), int(0*la/5):int(1*la/5)]
                b21 = b[int(1*lb/5):int(2*lb/5), int(0*lb/5):int(1*lb/5)]
                a22 = a[int(1*la/5):int(2*la/5), int(1*la/5):int(2*la/5)]
                b22 = b[int(1*lb/5):int(2*lb/5), int(1*lb/5):int(2*lb/5)]
                a23 = a[int(1*la/5):int(2*la/5), int(2*la/5):int(3*la/5)]
                b23 = b[int(1*lb/5):int(2*lb/5), int(2*lb/5):int(3*lb/5)]
                a24 = a[int(1*la/5):int(2*la/5), int(3*la/5):int(4*la/5)]
                b24 = b[int(1*lb/5):int(2*lb/5), int(3*lb/5):int(4*lb/5)]
                a25 = a[int(1*la/5):int(2*la/5), int(4*la/5):int(5*la/5)]
                b25 = b[int(1*lb/5):int(2*lb/5), int(4*lb/5):int(5*lb/5)]
                a31 = a[int(2*la/5):int(3*la/5), int(0*la/5):int(1*la/5)]
                b31 = b[int(2*lb/5):int(3*lb/5), int(0*lb/5):int(1*lb/5)]
                a32 = a[int(2*la/5):int(3*la/5), int(1*la/5):int(2*la/5)]
                b32 = b[int(2*lb/5):int(3*lb/5), int(1*lb/5):int(2*lb/5)]
                a33 = a[int(2*la/5):int(3*la/5), int(2*la/5):int(3*la/5)]
                b33 = b[int(2*lb/5):int(3*lb/5), int(2*lb/5):int(3*lb/5)]
                a34 = a[int(2*la/5):int(3*la/5), int(3*la/5):int(4*la/5)]
                b34 = b[int(2*lb/5):int(3*lb/5), int(3*lb/5):int(4*lb/5)]
                a35 = a[int(2*la/5):int(3*la/5), int(4*la/5):int(5*la/5)]
                b35 = b[int(2*lb/5):int(3*lb/5), int(4*lb/5):int(5*lb/5)]
                a41 = a[int(3*la/5):int(4*la/5), int(0*la/5):int(1*la/5)]
                b41 = b[int(3*lb/5):int(4*lb/5), int(0*lb/5):int(1*lb/5)]
                a42 = a[int(3*la/5):int(4*la/5), int(1*la/5):int(2*la/5)]
                b42 = b[int(3*lb/5):int(4*lb/5), int(1*lb/5):int(2*lb/5)]
                a43 = a[int(3*la/5):int(4*la/5), int(2*la/5):int(3*la/5)]
                b43 = b[int(3*lb/5):int(4*lb/5), int(2*lb/5):int(3*lb/5)]
                a44 = a[int(3*la/5):int(4*la/5), int(3*la/5):int(4*la/5)]
                b44 = b[int(3*lb/5):int(4*lb/5), int(3*lb/5):int(4*lb/5)]
                a45 = a[int(3*la/5):int(4*la/5), int(4*la/5):int(5*la/5)]
                b45 = b[int(3*lb/5):int(4*lb/5), int(4*lb/5):int(5*lb/5)]
                a51 = a[int(4*la/5):int(5*la/5), int(0*la/5):int(1*la/5)]
                b51 = b[int(4*lb/5):int(5*lb/5), int(0*lb/5):int(1*lb/5)]
                a52 = a[int(4*la/5):int(5*la/5), int(1*la/5):int(2*la/5)]
                b52 = b[int(4*lb/5):int(5*lb/5), int(1*lb/5):int(2*lb/5)]
                a53 = a[int(4*la/5):int(5*la/5), int(2*la/5):int(3*la/5)]
                b53 = b[int(4*lb/5):int(5*lb/5), int(2*lb/5):int(3*lb/5)]
                a54 = a[int(4*la/5):int(5*la/5), int(3*la/5):int(4*la/5)]
                b54 = b[int(4*lb/5):int(5*lb/5), int(3*lb/5):int(4*lb/5)]
                a55 = a[int(4*la/5):int(5*la/5), int(4*la/5):int(5*la/5)]
                b55 = b[int(4*lb/5):int(5*lb/5), int(4*lb/5):int(5*lb/5)]

                if la == 5 or lb == 5:
                    # Base case (5x5)
                    k1 =(a51*ZZ(-1)+a52*ZZ(2)+a54*ZZ(2)-a55*ZZ(2))*(b11*ZZ(3)+b21+b41)
                    k2 =(a53*ZZ(-1))*(b12*ZZ(2)-b14*ZZ(2)+b15*ZZ(4)-b32+b34-b35*ZZ(2)-b52+b54-b55*ZZ(2))
                    k3 =(a51-a52*ZZ(2)-a54+a55*ZZ(2))*(b11*ZZ(2)+b41)
                    k4 =(a53*ZZ(-1)+a55)*(b21+b51)
                    k6=(a41-a42-a44*ZZ(2)+a45*ZZ(2))*(b11*ZZ(-2)+b12*ZZ(2)-b21+b22)
                    k7=a44*(b12*ZZ(3)+b22+b42)
                    k8=(a44*ZZ(-1)-a51+a52*ZZ(2)+a54*ZZ(2)-a55*ZZ(2))*(b11*ZZ(2)+b12+b22+b41)
                    k9=(a44*ZZ(-1)+a54)*(b11*ZZ(-2)+b12*ZZ(2)-b41+b42)
                    k10=(a41*ZZ(-1)+a42*ZZ(2)+a44-a45*ZZ(2))*(b12+b22)
                    k11=(a45*ZZ(-1))*(b11*ZZ(2)-b12*ZZ(2)-b31+b32-b51+b52)
                    k12=(a43*ZZ(-1)+a45)*(b11*ZZ(-1)+b12-b21+b22-b31+b32)
                    k13=(a41-a42*ZZ(2)+a43-a44*ZZ(2)+a45-a51+a52*ZZ(2)+a54*ZZ(2)-a55*ZZ(2))*(b11*ZZ(-1)+b12-b21+b22)
                    k14=(a31*ZZ(-1)+a32*ZZ(2)+a34*ZZ(2)-a35*ZZ(2))*(b14*ZZ(3)-b15*ZZ(6)+b24-b25*ZZ(2)+b44-b45*ZZ(2))
                    k15=(a31-a32*ZZ(2)-a34*ZZ(2)+a35*ZZ(2)-a41+a42*ZZ(2)+a44*ZZ(2)-a45*ZZ(2))*(b11*ZZ(-1)+b14-b15*ZZ(2)-b21+b24-b25*ZZ(2)-b31+b32+b34-b35*ZZ(2))
                    k16=(a31*ZZ(-1)+a32*ZZ(2)+a34*ZZ(2)-a35*ZZ(2)-a44)*(b12+b14*ZZ(2)-b15*ZZ(4)+b22+b44-b45*ZZ(2))
                    k17=(a31*ZZ(-1)+a32*ZZ(2)+a34*ZZ(2)-a35*ZZ(2)+a41-a42*ZZ(2)+a43-a44*ZZ(2)+a45)*(b11*ZZ(-1)+b12-b21+b22-b31+b32+b34-b35*ZZ(2))
                    k18=(a31*ZZ(-1)+a32+a34*ZZ(2)-a35*ZZ(2)+a51-a52-a54*ZZ(2)+a55*ZZ(2))*(b22+b52)
                    k19=(a31*ZZ(-1)+a32+a34*ZZ(2)-a35*ZZ(2)+a41-a42-a44*ZZ(2)+a45*ZZ(2))*(b21-b22+b24-b25*ZZ(2)+b31-b32+b34-b35*ZZ(2)+b51-b52+b54-b55*ZZ(2))
                    k20=(a33-a53)*(b12*ZZ(-2)+b32+b52)
                    k21=(a33*ZZ(-1)+a43+a53)*b32
                    k22=(a34*ZZ(-1)+a54)*(b12-b14+b15-b42+b44-b45+b52-b54+b55)
                    k23=(a34*ZZ(-1)+a44)*(b11*ZZ(-2)+b12*ZZ(2)-b41+b42)
                    k24=(a34*ZZ(-1)+a44)*(b12*ZZ(-2)+b14*ZZ(2)-b15*ZZ(4)-b42+b44-b45*ZZ(2))
                    k25=(a31-a32*ZZ(2)-a34+a35*ZZ(2))*(b14*ZZ(2)-b15*ZZ(4)+b44-b45*ZZ(2))
                    k26=a35*(b24-b25*ZZ(2)+b34-b35*ZZ(2)+b54-b55*ZZ(2))
                    k27=(a35-a45-a53)*(b21+b31-b32+b51)
                    k28=(a31*ZZ(-1)+a32+a34*ZZ(2)-a35)*(b14*ZZ(2)-b15*ZZ(4)+b24-b25*ZZ(2))
                    k29=(a31-a32-a34*ZZ(2)+a35-a41+a42+a44*ZZ(2)-a45*ZZ(2))*(b11*ZZ(-2)+b12*ZZ(2)+b24-b25*ZZ(2)+b31-b32+b34-b35*ZZ(2)+b51-b52+b54-b55*ZZ(2))
                    k30=(a31-a32-a34*ZZ(2)+a35+a45-a51+a52+a54*ZZ(2)-a55)*(b11*ZZ(2)-b12*ZZ(2)+b21+b52)
                    k31=(a31*ZZ(-1)+a32+a34*ZZ(2)-a35+a41-a42-a44*ZZ(2)+a45)*(b11*ZZ(-2)+b12*ZZ(2)-b14*ZZ(2)+b15*ZZ(4)+b31-b32+b34-b35*ZZ(2)+b51-b52+b54-b55*ZZ(2))
                    k32=(a31-a32-a34*ZZ(2)+a35-a41+a42+a44*ZZ(2)-a45+a53)*(b12*ZZ(2)-b14*ZZ(2)+b15*ZZ(4)+b21+b31-b32+b34-b35*ZZ(2)+b51-b52+b54-b55*ZZ(2))
                    k33=(a31*ZZ(-1)+a32+a34*ZZ(2)-a35+a41-a42-a44*ZZ(2)+a45+a51-a52-a54*ZZ(2)+a55)*(b11*ZZ(2)+b21)
                    k34=(a33*ZZ(-1)+a35)*(b34-b35*ZZ(2))
                    k35=(a31*ZZ(-1)+a32*ZZ(2)-a33+a34*ZZ(2)-a35+a41-a42*ZZ(2)+a43-a44*ZZ(2)+a45)*(b31*ZZ(-1)+b32+b34-b35*ZZ(2))
                    k36=(a31*ZZ(-1)+a32-a33+a34*ZZ(2)-a35+a51-a52+a53-a54*ZZ(2)+a55)*(b12*ZZ(-2)+b52)
                    k37=(a34*ZZ(-1)-a35+a54+a55)*(b12*ZZ(-2)+b14*ZZ(2)-b15*ZZ(2)+b32-b34+b35+b52-b54+b55)
                    k38=a23*(b21-b22-b23+b25+b31-b32-b33+b35)
                    k39=(a22*ZZ(-1)+a23)*(b11*ZZ(-1)+b13+b21-b23-b41+b43)
                    k40=(a22-a23+a24*ZZ(2)+a25*ZZ(2)+a41-a42*ZZ(2)-a44*ZZ(4)+a45*ZZ(2))*(b11-b12-b13+b15+b21-b22-b23+b25+b31-b32*ZZ(2)-b33+b35)
                    k41=(a22-a23-a42+a43)*(b13*ZZ(-2)+b33+b53)
                    k42=(a22-a23+a24*ZZ(2)+a25*ZZ(2)+a41-a42*ZZ(2)+a43-a44*ZZ(4)+a45)*(b12*ZZ(-2)-b13*ZZ(2)+b14*ZZ(2)-b15*ZZ(2)-b21*ZZ(2)+b22*ZZ(2)+b23*ZZ(2)-b25*ZZ(2)-b31+b32*ZZ(3)+b33*ZZ(2)-b34+b41-b42-b43+b45+b52+b53-b54+b55)
                    k43=(a24*ZZ(-1)+a44)*(b15-b45+b55)
                    k44=(a25*ZZ(-1)*(b11*ZZ(-1)+b12*ZZ(2)+b13-b15*ZZ(2)+b41-b42*ZZ(2)-b43+b45*ZZ(2)-b51+b52*ZZ(2)+b53-b55*ZZ(2)))
                    k45=(a21-a22*ZZ(2)-a24*ZZ(4)-a25*ZZ(5)-a41+a42*ZZ(2)+a44*ZZ(4)-a45*ZZ(2))*(b11-b12*ZZ(2)-b13+b15+b21-b22*ZZ(2)-b23+b25+b31-b32*ZZ(2)-b33+b35)
                    k46=(a21*ZZ(-1)+a22+a24*ZZ(2)+a25*ZZ(3))*(b11*ZZ(2)-b12*ZZ(4)-b13*ZZ(2)+b15*ZZ(2)+b21-b22*ZZ(2)-b23+b25+b31-b32*ZZ(2)-b33+b35)
                    k47=(a24+a25-a44)*(b15-b22*ZZ(2)+b42-b45+b55)
                    k48=(a23*ZZ(-1)+a24+a25+a43-a44-a45)*(b12*ZZ(3)-b14*ZZ(3)+b15*ZZ(3)+b22-b24+b25+b32-b34+b35)
                    k49=(a22-a23+a24+a25)*(b11*ZZ(-2)+b13*ZZ(2)-b41+b43)
                    k50=(a22-a23+a24+a25+a41-a42*ZZ(2)+a43-a44*ZZ(3)+a45)*(b12*ZZ(-4)-b13*ZZ(4)+b14*ZZ(4)-b15*ZZ(4)-b21*ZZ(2)+b22+b23+b24-b25*ZZ(3)+b32+b33-b34+b35+b41-b42-b43+b45+b52+b53-b54+b55)
                    k51=(a22*ZZ(-1)+a23-a24-a25+a42-a43+a44+a45)*(b23+b33+b53)
                    k52=(a11-a12*ZZ(2)+a14*ZZ(2)+a15*ZZ(2)-a41+a42*ZZ(2)-a44*ZZ(2)-a45*ZZ(2)-a51+a52*ZZ(2)-a54*ZZ(2)-a55*ZZ(2))*(b12*ZZ(3)+b13*ZZ(3)-b14*ZZ(3)+b15*ZZ(4)+b22+b23-b24+b25+b32+b33-b34+b35-b43-b45+b53+b55)
                    k53=(a11-a12*ZZ(2)+a14*ZZ(2)+a15*ZZ(2)-a21+a22*ZZ(2)-a24*ZZ(2)-a25*ZZ(2)-a51+a52*ZZ(2)-a54*ZZ(2)-a55*ZZ(2))*b13
                    k54=(a11*ZZ(-1)+a12*ZZ(2)-a14*ZZ(2)-a15*ZZ(2)+a21-a22*ZZ(2)+a24*ZZ(2)+a25*ZZ(2)+a31-a32*ZZ(2)+a34*ZZ(2)+a35*ZZ(2))*(b12*ZZ(3)+b13*ZZ(3)-b14*ZZ(3)+b15*ZZ(4)+b22+b23-b24+b25-b42-b43+b44-b45+b52+b53-b54+b55)
                    k55=(a11-a12*ZZ(2)+a14*ZZ(2)+a15*ZZ(2)-a21+a22*ZZ(2)-a24*ZZ(2)-a25*ZZ(2)-a31+a32*ZZ(2)-a34*ZZ(3)-a35*ZZ(2)+a54)*(b12-b14+b15-b42-b43+b44-b45+b52+b53-b54+b55)
                    k56=(a11-a12*ZZ(2)+a14*ZZ(2)+a15*ZZ(2)+a24-a41+a42*ZZ(2)-a44*ZZ(3)-a45*ZZ(2)-a51+a52*ZZ(2)-a54*ZZ(2)-a55*ZZ(2))*(b15-b43-b45+b53+b55)
                    k57=(a11-a12*ZZ(2)+a14*ZZ(2)+a15*ZZ(2)-a23+a24+a25-a41+a42*ZZ(2)+a43-a44*ZZ(3)-a45*ZZ(3)-a51+a52*ZZ(2)-a54*ZZ(2)-a55*ZZ(2))*(b12*ZZ(3)+b13*ZZ(2)-b14*ZZ(3)+b15*ZZ(3)+b22+b23-b24+b25+b32+b33-b34+b35)
                    k58=(a11*ZZ(-1)+a12*ZZ(2)-a14*ZZ(2)-a15*ZZ(2)+a23-a24-a25+a31-a32*ZZ(2)+a34*ZZ(3)+a35*ZZ(2)+a41-a42*ZZ(2)-a43+a44*ZZ(3)+a45*ZZ(3)-a54)*(b12-b14+b15)
                    k59=(a11*ZZ(-1)+a12*ZZ(2)-a14*ZZ(2)-a15*ZZ(2)+a23-a24-a25+a33-a34-a35+a41-a42*ZZ(2)-a43+a44*ZZ(3)+a45*ZZ(3)+a51-a52*ZZ(2)-a53+a54*ZZ(3)+a55*ZZ(3))*(b12*ZZ(2)+b13*ZZ(2)-b14*ZZ(2)+b15*ZZ(2)+b22+b23-b24+b25+b32+b33-b34+b35)
                    k60=(a12-a22+a31-a32*ZZ(2)-a34*ZZ(2)+a35*ZZ(2)-a41+a42+a44*ZZ(2)-a45*ZZ(2))*(b12*ZZ(-5)+b14*ZZ(4)-b15*ZZ(3)+b24-b25*ZZ(2)-b42+b45+b52-b54+b55)
                    k61=(a12*ZZ(-1)+a23-a24-a25+a42-a43+a44+a45+a52)*(b13*ZZ(2)+b23+b33)
                    k5=(a12-a23+a24+a25-a32-a42+a43-a44-a45)*(b12*ZZ(2)-b14*ZZ(2)+b15*ZZ(2)+b22-b24+b25)
                    k62=(a12-a23+a24+a25-a33+a34+a35-a42+a43-a44-a45-a52+a53-a54-a55)*(b12*ZZ(2)+b13*ZZ(2)-b14*ZZ(2)+b15*ZZ(2)+b22+b23-b24+b25+b33)
                    k63=(a11-a12-a14*ZZ(2)-a15*ZZ(2)-a21+a22+a24*ZZ(2)+a25*ZZ(2)-a41+a42+a44*ZZ(2)+a45*ZZ(2))*(b12*ZZ(-1)+b15+b42-b45-b52+b55)
                    k64=(a11-a12-a14*ZZ(2)-a15*ZZ(2)+a23-a41+a42+a44*ZZ(2)+a45*ZZ(2))*(b11-b12-b13+b15)
                    k65=(a11*ZZ(-1)+a12+a14*ZZ(2)+a15*ZZ(2)+a21-a22-a24*ZZ(2)-a25*ZZ(3)+a41-a42-a44*ZZ(2)-a45*ZZ(2))*(b11-b12*ZZ(2)-b13+b15*ZZ(2)+b42-b45-b52+b55)
                    k66=(a13-a43)*(b15*ZZ(2)+b25+b35)
                    k67=(a13*ZZ(-1)-a22+a23*ZZ(2)-a24*ZZ(2)-a25*ZZ(2)-a41+a42*ZZ(2)-a43+a44*ZZ(4))*(b12*ZZ(-3)-b13*ZZ(4)+b14*ZZ(2)-b15-b21*ZZ(2)+b22+b23-b25+b32+b33-b34+b35+b41-b42-b43+b45+b52+b53-b54+b55)
                    k68=(a13-a23-a43)*(b11*ZZ(-1)+b12*ZZ(2)+b13-b15*ZZ(2)-b21+b22*ZZ(2)+b23-b25*ZZ(2)-b31+b32*ZZ(2)+b33-b35*ZZ(2))
                    k69=(a13*ZZ(-1)+a23+a43)*(b12*ZZ(-1)+b15-b22+b25-b32+b35)
                    k70=(a13-a23-a45)*(b12*ZZ(-4)+b14*ZZ(2)-b34+b35*ZZ(2)-b42+b45+b52-b54+b55)
                    k71=(a13-a21+a22+a24*ZZ(2)+a25*ZZ(3)-a43)*(b11-b12*ZZ(2)-b13+b21-b22*ZZ(2)-b23+b25+b31-b32*ZZ(2)-b33+b35)
                    k72=(a13-a23+a24+a25-a44-a45)*(b13*ZZ(-4)-b21*ZZ(2)+b23+b33+b41-b43+b53)
                    k73=(a12+a13-a23*ZZ(2)+a24*ZZ(2)+a25*ZZ(2)-a42+a43-a44*ZZ(2)-a45*ZZ(2)+a51-a52*ZZ(2)-a54*ZZ(2)+a55*ZZ(2))*(b11*ZZ(3)-b13*ZZ(3)+b21+b41-b43+b53)
                    k74=(a12*ZZ(-1)+a13+a42-a43+a52)*(b15*ZZ(2)+b25)
                    k75=(a12-a13-a22+a23-a42+a43)*(b12*ZZ(-1)+b15-b22+b25)
                    k76=(a14-a44-a54)*(b42*ZZ(-1)+b44-b45+b52-b54+b55)
                    k77=(a14-a24-a54)*(b43*ZZ(-1)+b53)
                    k78=(a11*ZZ(-1)+a12*ZZ(2)-a14*ZZ(3)-a15*ZZ(2)+a41-a42*ZZ(2)+a44*ZZ(3)+a45*ZZ(2)+a51-a52*ZZ(2)+a54*ZZ(3)+a55*ZZ(2))*(b43*ZZ(-1)-b45+b53+b55)
                    k79=(a11*ZZ(-1)+a12*ZZ(2)-a14*ZZ(3)-a15*ZZ(2)+a21-a22*ZZ(2)+a24*ZZ(3)+a25*ZZ(2)+a31-a32*ZZ(2)+a34*ZZ(3)+a35*ZZ(2))*(b42*ZZ(-1)-b43+b44-b45+b52+b53-b54+b55)
                    k80=(a15-a45)*(b41*ZZ(-1)+b42+b43-b45+b51-b52-b53+b55)
                    k81=(a11-a12-a14*ZZ(2)-a15*ZZ(3)-a21+a22+a24*ZZ(2)+a25*ZZ(3)-a41+a42+a44*ZZ(2)+a45*ZZ(3))*(b42-b45-b52+b55)
                    k82=(a14+a15-a44-a45)*(b11-b12-b13+b15-b21+b22+b23-b25+b41-b42-b43+b45)
                    k83=(a14+a15-a44-a45-a54)*(b12*ZZ(3)-b15*ZZ(3)+b22+b42-b45+b55)
                    k84=(a14*ZZ(-1)-a15+a44+a45+a54+a55)*(b13*ZZ(-2)+b33+b53)
                    k85=(a14+a15-a44-a45-a54-a55)*(b12*ZZ(-2)+b14*ZZ(2)-b15*ZZ(2)+b32-b34+b35+b52-b54+b55)
                    k86=(a14+a15-a44-a45+a53-a54-a55)*(b12*ZZ(2)-b14*ZZ(2)+b15*ZZ(4)-b32+b34-b35-b52+b54-b55*ZZ(2))
                    k87=(a14+a15+a22-a23-a44-a45)*(b11*ZZ(-1)-b12+b13+b15+b21-b22-b23+b25-b41+b43)
                    k88=(a14+a15+a22-a23-a34-a35-a42+a43-a44-a45)*(b12*ZZ(2)-b14*ZZ(2)+b15*ZZ(2)+b22+b23-b24+b25+b33+b53)
                    k89=(a14*ZZ(-1)-a15+a24+a25+a44+a45)*(b12*ZZ(-2)+b15*ZZ(2)-b42+b45)
                    k90=(a14*ZZ(-1)-a15+a24+a25+a34+a35)*(b22+b23-b24+b25+b32+b33-b34+b35+b52+b53-b54+b55)
                    k91=(a11-a12*ZZ(2)+a14*ZZ(3)+a15*ZZ(3)-a21+a22*ZZ(3)-a23-a24*ZZ(2)-a25*ZZ(2)-a31+a32*ZZ(3)-a33-a34*ZZ(2)-a35*ZZ(2)-a42+a43-a44-a45-a52+a53-a54-a55)*(b12*ZZ(2)+b13*ZZ(2)-b14*ZZ(2)+b15*ZZ(2)+b22+b23-b24+b25)
                    k92=(a12*ZZ(-1)-a13*ZZ(2)+a14+a15+a23*ZZ(3)-a24*ZZ(3)-a25*ZZ(3)+a42-a43+a44*ZZ(2)+a45*ZZ(2)-a51+a52*ZZ(2)+a54-a55*ZZ(2))*(b11*ZZ(2)-b13*ZZ(4)+b41-b43+b53)
                    k93=(a13-a14-a15-a23+a24+a25-a33+a34+a35)*(b32-b34+b35)
                    k94=(a12*ZZ(-1)+a13-a14-a15+a42-a43+a44+a45)*(b11-b12-b13+b15+b21-b22-b23+b25)
                    k95=(a12-a13-a14-a15-a42+a43+a44+a45+a51-a52*ZZ(2)-a54+a55*ZZ(2))*(b12+b15+b22)
                    k96=(a12-a13+a14+a15-a42+a43-a44-a45-a52+a53-a54-a55)*b33
                    k97=(a12*ZZ(-1)+a13-a14-a15+a22-a23+a24+a25-a31+a32*ZZ(2)+a34*ZZ(3)-a35*ZZ(2)+a41-a42-a44*ZZ(2)+a45*ZZ(2))*(b12*ZZ(-4)+b14*ZZ(2)-b42+b45+b52-b54+b55)


                    c11=k1*ZZ(2) +k3 -k8+k10+k11+k12+k13+k20+k21+k30+k33+k36+k38*ZZ(2)-k43+k46+k48+k52-k56-k57+k64*ZZ(2)+k65-k66-k68*ZZ(2)+k69*ZZ(2)-k71+k72+k73-k77+k80+k81+k92
                    c12=k1 -k8+k10+k12+k13+k38*ZZ(2)+k40+k42-k43+k44+k46+k48+k52-k53-k56-k57+k63*ZZ(2)+k64*ZZ(2)+k65*ZZ(2)-k66+k67-k68+k69+k70-k71+k72
                    c13=k1 +k2 -k8+k11+k13+k16-k17+k20+k23+k24+k25+k27+k30+k31+k32+k33-k34+k35+k36
                    c14=k1 -k8+k10+k11+k12+k13+k20+k21+k30+k33+k36
                    c15=k1 +k2 +k3 +k4 +k31+k32+k33
                    c21=k3 +k6+k7*ZZ(2)+k8+k9+k10+k20+k21+k30+k33+k36+k38+k43+k46+k56+k64+k65-k68+k69-k71+k77+k78+k81+k83+k95
                    c22=k7+k10+k38+k40+k43+k45+k46+k47+k63+k64+k65
                    c23=k2 +k6+k7+k11+k16+k19+k20+k24+k25+k26+k27+k29+k30+k31+k32+k33+k36
                    c24=k6+k7+k10+k20+k21+k30+k33+k36
                    c25=k2 +k3 +k6+k7+k8+k9+k11+k18+k19+k26+k27+k29+k30+k31+k32+k33
                    c31=k1*ZZ(2) +k3 -k8+k10+k12+k13+k39*ZZ(2)+k40+k41+k42-k43+k48+k49+k51+k52-k56-k57+k67+k69+k70+k72*ZZ(2)+k73+k75-k77+k87+k92+k94
                    c32=k1-k8+k10+k12+k13+k39+k40+k42-k43+k48+k49+k52-k53-k56-k57+k67+k69+k70+k72
                    c33=k1 +k3 -k22+k39+k41+k53+k54*ZZ(2)+k55+k5-k62+k72+k73+k75+k87-k88+k91*ZZ(2)+k92+k94+k96
                    c34=k1 -k8+k10+k12+k13+k39+k40+k41+k42+k49+k51+k67+k69+k70+k72
                    c35=k1 +k3 +k39+k61+k72+k73+k75+k84+k87+k92+k94+k96
                    c41=k3*ZZ(2) +k6*ZZ(2)+k7*ZZ(4)+k8*ZZ(2)+k9*ZZ(2)+k10*ZZ(3)+k11+k12+k14*ZZ(2)+k15-k16+k17+k20*ZZ(2)+k21*ZZ(2)+k22+k25+k28+k29+k30*ZZ(2)+k31+k33*ZZ(2)+k36*ZZ(2)+k37*ZZ(2)+k38+k39*ZZ(4)+k43+k46+k48+k50*ZZ(3)+k51*ZZ(2)+k56*ZZ(2)+k58+k60+k64+k65+k66+k67*ZZ(3)-k68+k69*ZZ(3)+k70*ZZ(2)-k71+k72*ZZ(3)+k75*ZZ(2)+k76+k77*ZZ(2)+k78*ZZ(2)+k81+k82*ZZ(4)+k83*ZZ(2)+k85*ZZ(2)+k87*ZZ(4)+k88*ZZ(2)+k89*ZZ(2)+k90*ZZ(2)+k95*ZZ(2)+k97
                    c42=k7*ZZ(2)+k10*ZZ(3)+k12+k14+k15-k16+k17+k38+k39*ZZ(2)+k40*ZZ(2)+k43*ZZ(2)+k45*ZZ(2)+k46+k47*ZZ(2)+k48+k50+k55+k58+k63+k64+k65+k66+k67+k68+k69+k70+k71+k72+k76+k77+k79+k82*ZZ(2)+k87*ZZ(2)+k89*ZZ(2)
                    c43=k2*ZZ(2) +k6*ZZ(2)+k7*ZZ(2)+k11*ZZ(2)+k14+k16*ZZ(2)+k19*ZZ(2)+k20*ZZ(2)+k24*ZZ(2)+k25*ZZ(3)+k26+k27*ZZ(2)+k28+k29*ZZ(2)+k30*ZZ(2)+k31*ZZ(2)+k32*ZZ(2)+k33*ZZ(2)+k34+k36*ZZ(2)+k37*ZZ(2)+k39*ZZ(2)+k50*ZZ(2)+k60*ZZ(2)+k5*ZZ(2)+k67*ZZ(2)+k72*ZZ(2)+k75*ZZ(2)+k82*ZZ(2)+k85*ZZ(2)+k87*ZZ(2)+k93*ZZ(2)+k97*ZZ(2)
                    c44=k6*ZZ(2)+k7*ZZ(2)+k10*ZZ(3)+k11+k12+k14+k15-k16+k17+k20*ZZ(2)+k21*ZZ(2)+k28+k29+k30*ZZ(2)+k31+k33*ZZ(2)+k36*ZZ(2)+k37*ZZ(2)+k39*ZZ(2)+k50*ZZ(2)+k51*ZZ(2)+k67*ZZ(2)+k69*ZZ(2)+k70*ZZ(2)+k72*ZZ(2)+k82*ZZ(2)+k85*ZZ(2)+k87*ZZ(2)+k88*ZZ(2)+k89*ZZ(2)+k90*ZZ(2)
                    c45=k2*ZZ(2) +k3*ZZ(2) +k6+k7*ZZ(2)+k8*ZZ(2)+k9*ZZ(2)+k11+k14+k18+k19+k22+k25+k26+k27+k29+k30+k31+k32+k33+k39*ZZ(2)+k50+k57+k58+k59+k60+k61+k62+k67+k72+k74+k75*ZZ(2)+k82*ZZ(2)+k83+k85*ZZ(2)+k86+k87*ZZ(2)+k93+k95+k97
                    c51=k3 +k6+k7*ZZ(2)+k8+k9+k10+k11+k20+k21+k28+k29+k30+k31+k33+k36+k37+k39*ZZ(2)+k43+k50+k51+k56+k66+k67+k69+k70+k72+k75+k77+k78+k82*ZZ(2)+k83+k85+k87*ZZ(2)+k88+k89+k90+k95
                    c52=k7+k10+k39+k40+k43+k45+k47+k66+k68+k71+k82+k87+k89
                    c53=k2 +k6+k7+k11+k16+k19+k20+k24+k25+k27+k28+k29+k30+k31+k32+k33+k34+k36+k37+k39+k50+k60+k5+k67+k72+k75+k82+k85+k87+k93+k97
                    c54=k6+k7+k10+k11+k20+k21+k28+k29+k30+k31+k33+k36+k37+k39+k50+k51+k67+k69+k70+k72+k82+k85+k87+k88+k89+k90
                    c55=k2 +k3 +k7+k8+k9+k39+k74+k75+k82+k83+k85+k86+k87+k95

                    # Create product matrix C entries
                    c = [[c11.entries[0][0], c12.entries[0][0], c13.entries[0][0], c14.entries[0][0], c15.entries[0][0]],
                        [c21.entries[0][0], c22.entries[0][0], c23.entries[0][0], c24.entries[0][0], c25.entries[0][0]],
                        [c31.entries[0][0], c32.entries[0][0], c33.entries[0][0], c34.entries[0][0], c35.entries[0][0]],
                        [c41.entries[0][0], c42.entries[0][0], c43.entries[0][0], c44.entries[0][0], c45.entries[0][0]],
                        [c51.entries[0][0], c52.entries[0][0], c53.entries[0][0], c54.entries[0][0], c55.entries[0][0]]]
                    
                    # Return C as a generic matrix
                    return _MatrixGenericData(
                    base_ring=a.base_ring,
                    nrows=len(c),
                    ncols=len(c[0]),
                    entries=c
                    ) 
                

                k1=kauers_moosbauer(a51*ZZ(-1)+a52*ZZ(2)+a54*ZZ(2)-a55*ZZ(2), b11*ZZ(3)+b21+b41)
                k2=kauers_moosbauer(a53*ZZ(-1), b12*ZZ(2)-b14*ZZ(2)+b15*ZZ(4)-b32+b34-b35*ZZ(2)-b52+b54-b55*ZZ(2))
                k3=kauers_moosbauer(a51-a52*ZZ(2)-a54+a55*ZZ(2), b11*ZZ(2)+b41)
                k4=kauers_moosbauer(a53*ZZ(-1)+a55, b21+b51)
                k6=kauers_moosbauer(a41-a42-a44*ZZ(2)+a45*ZZ(2), b11*ZZ(-2)+b12*ZZ(2)-b21+b22)
                k7=kauers_moosbauer(a44, b12*ZZ(3)+b22+b42)
                k8=kauers_moosbauer(a44*ZZ(-1)-a51+a52*ZZ(2)+a54*ZZ(2)-a55*ZZ(2), b11*ZZ(2)+b12+b22+b41)
                k9=kauers_moosbauer(a44*ZZ(-1)+a54, b11*ZZ(-2)+b12*ZZ(2)-b41+b42)
                k10=kauers_moosbauer(a41*ZZ(-1)+a42*ZZ(2)+a44-a45*ZZ(2), b12+b22)
                k11=kauers_moosbauer(a45*ZZ(-1), b11*ZZ(2)-b12*ZZ(2)-b31+b32-b51+b52)
                k12=kauers_moosbauer(a43*ZZ(-1)+a45, b11*ZZ(-1)+b12-b21+b22-b31+b32)
                k13=kauers_moosbauer(a41-a42*ZZ(2)+a43-a44*ZZ(2)+a45-a51+a52*ZZ(2)+a54*ZZ(2)-a55*ZZ(2), b11*ZZ(-1)+b12-b21+b22)
                k14=kauers_moosbauer(a31*ZZ(-1)+a32*ZZ(2)+a34*ZZ(2)-a35*ZZ(2), b14*ZZ(3)-b15*ZZ(6)+b24-b25*ZZ(2)+b44-b45*ZZ(2))
                k15=kauers_moosbauer(a31-a32*ZZ(2)-a34*ZZ(2)+a35*ZZ(2)-a41+a42*ZZ(2)+a44*ZZ(2)-a45*ZZ(2), b11*ZZ(-1)+b14-b15*ZZ(2)-b21+b24-b25*ZZ(2)-b31+b32+b34-b35*ZZ(2))
                k16=kauers_moosbauer(a31*ZZ(-1)+a32*ZZ(2)+a34*ZZ(2)-a35*ZZ(2)-a44, b12+b14*ZZ(2)-b15*ZZ(4)+b22+b44-b45*ZZ(2))
                k17=kauers_moosbauer(a31*ZZ(-1)+a32*ZZ(2)+a34*ZZ(2)-a35*ZZ(2)+a41-a42*ZZ(2)+a43-a44*ZZ(2)+a45, b11*ZZ(-1)+b12-b21+b22-b31+b32+b34-b35*ZZ(2))
                k18=kauers_moosbauer(a31*ZZ(-1)+a32+a34*ZZ(2)-a35*ZZ(2)+a51-a52-a54*ZZ(2)+a55*ZZ(2), b22+b52)
                k19=kauers_moosbauer(a31*ZZ(-1)+a32+a34*ZZ(2)-a35*ZZ(2)+a41-a42-a44*ZZ(2)+a45*ZZ(2), b21-b22+b24-b25*ZZ(2)+b31-b32+b34-b35*ZZ(2)+b51-b52+b54-b55*ZZ(2))
                k20=kauers_moosbauer(a33-a53, b12*ZZ(-2)+b32+b52)
                k21=kauers_moosbauer(a33*ZZ(-1)+a43+a53, b32)
                k22=kauers_moosbauer(a34*ZZ(-1)+a54, b12-b14+b15-b42+b44-b45+b52-b54+b55)
                k23=kauers_moosbauer(a34*ZZ(-1)+a44, b11*ZZ(-2)+b12*ZZ(2)-b41+b42)
                k24=kauers_moosbauer(a34*ZZ(-1)+a44, b12*ZZ(-2)+b14*ZZ(2)-b15*ZZ(4)-b42+b44-b45*ZZ(2))
                k25=kauers_moosbauer(a31-a32*ZZ(2)-a34+a35*ZZ(2), b14*ZZ(2)-b15*ZZ(4)+b44-b45*ZZ(2))
                k26=kauers_moosbauer(a35, b24-b25*ZZ(2)+b34-b35*ZZ(2)+b54-b55*ZZ(2))
                k27=kauers_moosbauer(a35-a45-a53, b21+b31-b32+b51)
                k28=kauers_moosbauer(a31*ZZ(-1)+a32+a34*ZZ(2)-a35, b14*ZZ(2)-b15*ZZ(4)+b24-b25*ZZ(2))
                k29=kauers_moosbauer(a31-a32-a34*ZZ(2)+a35-a41+a42+a44*ZZ(2)-a45*ZZ(2), b11*ZZ(-2)+b12*ZZ(2)+b24-b25*ZZ(2)+b31-b32+b34-b35*ZZ(2)+b51-b52+b54-b55*ZZ(2))
                k30=kauers_moosbauer(a31-a32-a34*ZZ(2)+a35+a45-a51+a52+a54*ZZ(2)-a55, b11*ZZ(2)-b12*ZZ(2)+b21+b52)
                k31=kauers_moosbauer(a31*ZZ(-1)+a32+a34*ZZ(2)-a35+a41-a42-a44*ZZ(2)+a45, b11*ZZ(-2)+b12*ZZ(2)-b14*ZZ(2)+b15*ZZ(4)+b31-b32+b34-b35*ZZ(2)+b51-b52+b54-b55*ZZ(2))
                k32=kauers_moosbauer(a31-a32-a34*ZZ(2)+a35-a41+a42+a44*ZZ(2)-a45+a53, b12*ZZ(2)-b14*ZZ(2)+b15*ZZ(4)+b21+b31-b32+b34-b35*ZZ(2)+b51-b52+b54-b55*ZZ(2))
                k33=kauers_moosbauer(a31*ZZ(-1)+a32+a34*ZZ(2)-a35+a41-a42-a44*ZZ(2)+a45+a51-a52-a54*ZZ(2)+a55, b11*ZZ(2)+b21)
                k34=kauers_moosbauer(a33*ZZ(-1)+a35, b34-b35*ZZ(2))
                k35=kauers_moosbauer(a31*ZZ(-1)+a32*ZZ(2)-a33+a34*ZZ(2)-a35+a41-a42*ZZ(2)+a43-a44*ZZ(2)+a45, b31*ZZ(-1)+b32+b34-b35*ZZ(2))
                k36=kauers_moosbauer(a31*ZZ(-1)+a32-a33+a34*ZZ(2)-a35+a51-a52+a53-a54*ZZ(2)+a55, b12*ZZ(-2)+b52)
                k37=kauers_moosbauer(a34*ZZ(-1)-a35+a54+a55, b12*ZZ(-2)+b14*ZZ(2)-b15*ZZ(2)+b32-b34+b35+b52-b54+b55)
                k38=kauers_moosbauer(a23, b21-b22-b23+b25+b31-b32-b33+b35)
                k39=kauers_moosbauer(a22*ZZ(-1)+a23, b11*ZZ(-1)+b13+b21-b23-b41+b43)
                k40=kauers_moosbauer(a22-a23+a24*ZZ(2)+a25*ZZ(2)+a41-a42*ZZ(2)-a44*ZZ(4)+a45*ZZ(2), b11-b12-b13+b15+b21-b22-b23+b25+b31-b32*ZZ(2)-b33+b35)
                k41=kauers_moosbauer(a22-a23-a42+a43, b13*ZZ(-2)+b33+b53)
                k42=kauers_moosbauer(a22-a23+a24*ZZ(2)+a25*ZZ(2)+a41-a42*ZZ(2)+a43-a44*ZZ(4)+a45, b12*ZZ(-2)-b13*ZZ(2)+b14*ZZ(2)-b15*ZZ(2)-b21*ZZ(2)+b22*ZZ(2)+b23*ZZ(2)-b25*ZZ(2)-b31+b32*ZZ(3)+b33*ZZ(2)-b34+b41-b42-b43+b45+b52+b53-b54+b55)
                k43=kauers_moosbauer(a24*ZZ(-1)+a44, b15-b45+b55)
                k44=kauers_moosbauer(a25*ZZ(-1), b11*ZZ(-1)+b12*ZZ(2)+b13-b15*ZZ(2)+b41-b42*ZZ(2)-b43+b45*ZZ(2)-b51+b52*ZZ(2)+b53-b55*ZZ(2))
                k45=kauers_moosbauer(a21-a22*ZZ(2)-a24*ZZ(4)-a25*ZZ(5)-a41+a42*ZZ(2)+a44*ZZ(4)-a45*ZZ(2), b11-b12*ZZ(2)-b13+b15+b21-b22*ZZ(2)-b23+b25+b31-b32*ZZ(2)-b33+b35)
                k46=kauers_moosbauer(a21*ZZ(-1)+a22+a24*ZZ(2)+a25*ZZ(3), b11*ZZ(2)-b12*ZZ(4)-b13*ZZ(2)+b15*ZZ(2)+b21-b22*ZZ(2)-b23+b25+b31-b32*ZZ(2)-b33+b35)
                k47=kauers_moosbauer(a24+a25-a44, b15-b22*ZZ(2)+b42-b45+b55)
                k48=kauers_moosbauer(a23*ZZ(-1)+a24+a25+a43-a44-a45, b12*ZZ(3)-b14*ZZ(3)+b15*ZZ(3)+b22-b24+b25+b32-b34+b35)
                k49=kauers_moosbauer(a22-a23+a24+a25, b11*ZZ(-2)+b13*ZZ(2)-b41+b43)
                k50=kauers_moosbauer(a22-a23+a24+a25+a41-a42*ZZ(2)+a43-a44*ZZ(3)+a45, b12*ZZ(-4)-b13*ZZ(4)+b14*ZZ(4)-b15*ZZ(4)-b21*ZZ(2)+b22+b23+b24-b25*ZZ(3)+b32+b33-b34+b35+b41-b42-b43+b45+b52+b53-b54+b55)
                k51=kauers_moosbauer(a22*ZZ(-1)+a23-a24-a25+a42-a43+a44+a45, b23+b33+b53)
                k52=kauers_moosbauer(a11-a12*ZZ(2)+a14*ZZ(2)+a15*ZZ(2)-a41+a42*ZZ(2)-a44*ZZ(2)-a45*ZZ(2)-a51+a52*ZZ(2)-a54*ZZ(2)-a55*ZZ(2), b12*ZZ(3)+b13*ZZ(3)-b14*ZZ(3)+b15*ZZ(4)+b22+b23-b24+b25+b32+b33-b34+b35-b43-b45+b53+b55)
                k53=kauers_moosbauer(a11-a12*ZZ(2)+a14*ZZ(2)+a15*ZZ(2)-a21+a22*ZZ(2)-a24*ZZ(2)-a25*ZZ(2)-a51+a52*ZZ(2)-a54*ZZ(2)-a55*ZZ(2), b13)
                k54=kauers_moosbauer(a11*ZZ(-1)+a12*ZZ(2)-a14*ZZ(2)-a15*ZZ(2)+a21-a22*ZZ(2)+a24*ZZ(2)+a25*ZZ(2)+a31-a32*ZZ(2)+a34*ZZ(2)+a35*ZZ(2), b12*ZZ(3)+b13*ZZ(3)-b14*ZZ(3)+b15*ZZ(4)+b22+b23-b24+b25-b42-b43+b44-b45+b52+b53-b54+b55)
                k55=kauers_moosbauer(a11-a12*ZZ(2)+a14*ZZ(2)+a15*ZZ(2)-a21+a22*ZZ(2)-a24*ZZ(2)-a25*ZZ(2)-a31+a32*ZZ(2)-a34*ZZ(3)-a35*ZZ(2)+a54, b12-b14+b15-b42-b43+b44-b45+b52+b53-b54+b55)
                k56=kauers_moosbauer(a11-a12*ZZ(2)+a14*ZZ(2)+a15*ZZ(2)+a24-a41+a42*ZZ(2)-a44*ZZ(3)-a45*ZZ(2)-a51+a52*ZZ(2)-a54*ZZ(2)-a55*ZZ(2), b15-b43-b45+b53+b55)
                k57=kauers_moosbauer(a11-a12*ZZ(2)+a14*ZZ(2)+a15*ZZ(2)-a23+a24+a25-a41+a42*ZZ(2)+a43-a44*ZZ(3)-a45*ZZ(3)-a51+a52*ZZ(2)-a54*ZZ(2)-a55*ZZ(2), b12*ZZ(3)+b13*ZZ(2)-b14*ZZ(3)+b15*ZZ(3)+b22+b23-b24+b25+b32+b33-b34+b35)
                k58=kauers_moosbauer(a11*ZZ(-1)+a12*ZZ(2)-a14*ZZ(2)-a15*ZZ(2)+a23-a24-a25+a31-a32*ZZ(2)+a34*ZZ(3)+a35*ZZ(2)+a41-a42*ZZ(2)-a43+a44*ZZ(3)+a45*ZZ(3)-a54, b12-b14+b15)
                k59=kauers_moosbauer(a11*ZZ(-1)+a12*ZZ(2)-a14*ZZ(2)-a15*ZZ(2)+a23-a24-a25+a33-a34-a35+a41-a42*ZZ(2)-a43+a44*ZZ(3)+a45*ZZ(3)+a51-a52*ZZ(2)-a53+a54*ZZ(3)+a55*ZZ(3), b12*ZZ(2)+b13*ZZ(2)-b14*ZZ(2)+b15*ZZ(2)+b22+b23-b24+b25+b32+b33-b34+b35)
                k60=kauers_moosbauer(a12-a22+a31-a32*ZZ(2)-a34*ZZ(2)+a35*ZZ(2)-a41+a42+a44*ZZ(2)-a45*ZZ(2), b12*ZZ(-5)+b14*ZZ(4)-b15*ZZ(3)+b24-b25*ZZ(2)-b42+b45+b52-b54+b55)
                k61=kauers_moosbauer(a12*ZZ(-1)+a23-a24-a25+a42-a43+a44+a45+a52, b13*ZZ(2)+b23+b33)
                k5=kauers_moosbauer(a12-a23+a24+a25-a32-a42+a43-a44-a45, b12*ZZ(2)-b14*ZZ(2)+b15*ZZ(2)+b22-b24+b25)
                k62=kauers_moosbauer(a12-a23+a24+a25-a33+a34+a35-a42+a43-a44-a45-a52+a53-a54-a55, b12*ZZ(2)+b13*ZZ(2)-b14*ZZ(2)+b15*ZZ(2)+b22+b23-b24+b25+b33)
                k63=kauers_moosbauer(a11-a12-a14*ZZ(2)-a15*ZZ(2)-a21+a22+a24*ZZ(2)+a25*ZZ(2)-a41+a42+a44*ZZ(2)+a45*ZZ(2), b12*ZZ(-1)+b15+b42-b45-b52+b55)
                k64=kauers_moosbauer(a11-a12-a14*ZZ(2)-a15*ZZ(2)+a23-a41+a42+a44*ZZ(2)+a45*ZZ(2), b11-b12-b13+b15)
                k65=kauers_moosbauer(a11*ZZ(-1)+a12+a14*ZZ(2)+a15*ZZ(2)+a21-a22-a24*ZZ(2)-a25*ZZ(3)+a41-a42-a44*ZZ(2)-a45*ZZ(2), b11-b12*ZZ(2)-b13+b15*ZZ(2)+b42-b45-b52+b55)
                k66=kauers_moosbauer(a13-a43, b15*ZZ(2)+b25+b35)
                k67=kauers_moosbauer(a13*ZZ(-1)-a22+a23*ZZ(2)-a24*ZZ(2)-a25*ZZ(2)-a41+a42*ZZ(2)-a43+a44*ZZ(4), b12*ZZ(-3)-b13*ZZ(4)+b14*ZZ(2)-b15-b21*ZZ(2)+b22+b23-b25+b32+b33-b34+b35+b41-b42-b43+b45+b52+b53-b54+b55)
                k68=kauers_moosbauer(a13-a23-a43, b11*ZZ(-1)+b12*ZZ(2)+b13-b15*ZZ(2)-b21+b22*ZZ(2)+b23-b25*ZZ(2)-b31+b32*ZZ(2)+b33-b35*ZZ(2))
                k69=kauers_moosbauer(a13*ZZ(-1)+a23+a43, b12*ZZ(-1)+b15-b22+b25-b32+b35)
                k70=kauers_moosbauer(a13-a23-a45, b12*ZZ(-4)+b14*ZZ(2)-b34+b35*ZZ(2)-b42+b45+b52-b54+b55)
                k71=kauers_moosbauer(a13-a21+a22+a24*ZZ(2)+a25*ZZ(3)-a43, b11-b12*ZZ(2)-b13+b21-b22*ZZ(2)-b23+b25+b31-b32*ZZ(2)-b33+b35)
                k72=kauers_moosbauer(a13-a23+a24+a25-a44-a45, b13*ZZ(-4)-b21*ZZ(2)+b23+b33+b41-b43+b53)
                k73=kauers_moosbauer(a12+a13-a23*ZZ(2)+a24*ZZ(2)+a25*ZZ(2)-a42+a43-a44*ZZ(2)-a45*ZZ(2)+a51-a52*ZZ(2)-a54*ZZ(2)+a55*ZZ(2), b11*ZZ(3)-b13*ZZ(3)+b21+b41-b43+b53)
                k74=kauers_moosbauer(a12*ZZ(-1)+a13+a42-a43+a52, b15*ZZ(2)+b25)
                k75=kauers_moosbauer(a12-a13-a22+a23-a42+a43, b12*ZZ(-1)+b15-b22+b25)
                k76=kauers_moosbauer(a14-a44-a54, b42*ZZ(-1)+b44-b45+b52-b54+b55)
                k77=kauers_moosbauer(a14-a24-a54, b43*ZZ(-1)+b53)
                k78=kauers_moosbauer(a11*ZZ(-1)+a12*ZZ(2)-a14*ZZ(3)-a15*ZZ(2)+a41-a42*ZZ(2)+a44*ZZ(3)+a45*ZZ(2)+a51-a52*ZZ(2)+a54*ZZ(3)+a55*ZZ(2), b43*ZZ(-1)-b45+b53+b55)
                k79=kauers_moosbauer(a11*ZZ(-1)+a12*ZZ(2)-a14*ZZ(3)-a15*ZZ(2)+a21-a22*ZZ(2)+a24*ZZ(3)+a25*ZZ(2)+a31-a32*ZZ(2)+a34*ZZ(3)+a35*ZZ(2), b42*ZZ(-1)-b43+b44-b45+b52+b53-b54+b55)
                k80=kauers_moosbauer(a15-a45, b41*ZZ(-1)+b42+b43-b45+b51-b52-b53+b55)
                k81=kauers_moosbauer(a11-a12-a14*ZZ(2)-a15*ZZ(3)-a21+a22+a24*ZZ(2)+a25*ZZ(3)-a41+a42+a44*ZZ(2)+a45*ZZ(3), b42-b45-b52+b55)
                k82=kauers_moosbauer(a14+a15-a44-a45, b11-b12-b13+b15-b21+b22+b23-b25+b41-b42-b43+b45)
                k83=kauers_moosbauer(a14+a15-a44-a45-a54, b12*ZZ(3)-b15*ZZ(3)+b22+b42-b45+b55)
                k84=kauers_moosbauer(a14*ZZ(-1)-a15+a44+a45+a54+a55, b13*ZZ(-2)+b33+b53)
                k85=kauers_moosbauer(a14+a15-a44-a45-a54-a55, b12*ZZ(-2)+b14*ZZ(2)-b15*ZZ(2)+b32-b34+b35+b52-b54+b55)
                k86=kauers_moosbauer(a14+a15-a44-a45+a53-a54-a55, b12*ZZ(2)-b14*ZZ(2)+b15*ZZ(4)-b32+b34-b35-b52+b54-b55*ZZ(2))
                k87=kauers_moosbauer(a14+a15+a22-a23-a44-a45, b11*ZZ(-1)-b12+b13+b15+b21-b22-b23+b25-b41+b43)
                k88=kauers_moosbauer(a14+a15+a22-a23-a34-a35-a42+a43-a44-a45, b12*ZZ(2)-b14*ZZ(2)+b15*ZZ(2)+b22+b23-b24+b25+b33+b53)
                k89=kauers_moosbauer(a14*ZZ(-1)-a15+a24+a25+a44+a45, b12*ZZ(-2)+b15*ZZ(2)-b42+b45)
                k90=kauers_moosbauer(a14*ZZ(-1)-a15+a24+a25+a34+a35, b22+b23-b24+b25+b32+b33-b34+b35+b52+b53-b54+b55)
                k91=kauers_moosbauer(a11-a12*ZZ(2)+a14*ZZ(3)+a15*ZZ(3)-a21+a22*ZZ(3)-a23-a24*ZZ(2)-a25*ZZ(2)-a31+a32*ZZ(3)-a33-a34*ZZ(2)-a35*ZZ(2)-a42+a43-a44-a45-a52+a53-a54-a55, b12*ZZ(2)+b13*ZZ(2)-b14*ZZ(2)+b15*ZZ(2)+b22+b23-b24+b25)
                k92=kauers_moosbauer(a12*ZZ(-1)-a13*ZZ(2)+a14+a15+a23*ZZ(3)-a24*ZZ(3)-a25*ZZ(3)+a42-a43+a44*ZZ(2)+a45*ZZ(2)-a51+a52*ZZ(2)+a54-a55*ZZ(2), b11*ZZ(2)-b13*ZZ(4)+b41-b43+b53)
                k93=kauers_moosbauer(a13-a14-a15-a23+a24+a25-a33+a34+a35, b32-b34+b35)
                k94=kauers_moosbauer(a12*ZZ(-1)+a13-a14-a15+a42-a43+a44+a45, b11-b12-b13+b15+b21-b22-b23+b25)
                k95=kauers_moosbauer(a12-a13-a14-a15-a42+a43+a44+a45+a51-a52*ZZ(2)-a54+a55*ZZ(2), b12+b15+b22)
                k96=kauers_moosbauer(a12-a13+a14+a15-a42+a43-a44-a45-a52+a53-a54-a55, b33)
                k97=kauers_moosbauer(a12*ZZ(-1)+a13-a14-a15+a22-a23+a24+a25-a31+a32*ZZ(2)+a34*ZZ(3)-a35*ZZ(2)+a41-a42-a44*ZZ(2)+a45*ZZ(2), b12*ZZ(-4)+b14*ZZ(2)-b42+b45+b52-b54+b55)

                c11=k1*ZZ(2) +k3 -k8+k10+k11+k12+k13+k20+k21+k30+k33+k36+k38*ZZ(2)-k43+k46+k48+k52-k56-k57+k64*ZZ(2)+k65-k66-k68*ZZ(2)+k69*ZZ(2)-k71+k72+k73-k77+k80+k81+k92
                c12=k1 -k8+k10+k12+k13+k38*ZZ(2)+k40+k42-k43+k44+k46+k48+k52-k53-k56-k57+k63*ZZ(2)+k64*ZZ(2)+k65*ZZ(2)-k66+k67-k68+k69+k70-k71+k72
                c13=k1 +k2 -k8+k11+k13+k16-k17+k20+k23+k24+k25+k27+k30+k31+k32+k33-k34+k35+k36
                c14=k1 -k8+k10+k11+k12+k13+k20+k21+k30+k33+k36
                c15=k1 +k2 +k3 +k4 +k31+k32+k33
                c21=k3 +k6+k7*ZZ(2)+k8+k9+k10+k20+k21+k30+k33+k36+k38+k43+k46+k56+k64+k65-k68+k69-k71+k77+k78+k81+k83+k95
                c22=k7+k10+k38+k40+k43+k45+k46+k47+k63+k64+k65
                c23=k2 +k6+k7+k11+k16+k19+k20+k24+k25+k26+k27+k29+k30+k31+k32+k33+k36
                c24=k6+k7+k10+k20+k21+k30+k33+k36
                c25=k2 +k3 +k6+k7+k8+k9+k11+k18+k19+k26+k27+k29+k30+k31+k32+k33
                c31=k1*ZZ(2) +k3 -k8+k10+k12+k13+k39*ZZ(2)+k40+k41+k42-k43+k48+k49+k51+k52-k56-k57+k67+k69+k70+k72*ZZ(2)+k73+k75-k77+k87+k92+k94
                c32=k1-k8+k10+k12+k13+k39+k40+k42-k43+k48+k49+k52-k53-k56-k57+k67+k69+k70+k72
                c33=k1 +k3 -k22+k39+k41+k53+k54*ZZ(2)+k55+k5-k62+k72+k73+k75+k87-k88+k91*ZZ(2)+k92+k94+k96
                c34=k1 -k8+k10+k12+k13+k39+k40+k41+k42+k49+k51+k67+k69+k70+k72
                c35=k1 +k3 +k39+k61+k72+k73+k75+k84+k87+k92+k94+k96
                c41=k3*ZZ(2) +k6*ZZ(2)+k7*ZZ(4)+k8*ZZ(2)+k9*ZZ(2)+k10*ZZ(3)+k11+k12+k14*ZZ(2)+k15-k16+k17+k20*ZZ(2)+k21*ZZ(2)+k22+k25+k28+k29+k30*ZZ(2)+k31+k33*ZZ(2)+k36*ZZ(2)+k37*ZZ(2)+k38+k39*ZZ(4)+k43+k46+k48+k50*ZZ(3)+k51*ZZ(2)+k56*ZZ(2)+k58+k60+k64+k65+k66+k67*ZZ(3)-k68+k69*ZZ(3)+k70*ZZ(2)-k71+k72*ZZ(3)+k75*ZZ(2)+k76+k77*ZZ(2)+k78*ZZ(2)+k81+k82*ZZ(4)+k83*ZZ(2)+k85*ZZ(2)+k87*ZZ(4)+k88*ZZ(2)+k89*ZZ(2)+k90*ZZ(2)+k95*ZZ(2)+k97
                c42=k7*ZZ(2)+k10*ZZ(3)+k12+k14+k15-k16+k17+k38+k39*ZZ(2)+k40*ZZ(2)+k43*ZZ(2)+k45*ZZ(2)+k46+k47*ZZ(2)+k48+k50+k55+k58+k63+k64+k65+k66+k67+k68+k69+k70+k71+k72+k76+k77+k79+k82*ZZ(2)+k87*ZZ(2)+k89*ZZ(2)
                c43=k2*ZZ(2) +k6*ZZ(2)+k7*ZZ(2)+k11*ZZ(2)+k14+k16*ZZ(2)+k19*ZZ(2)+k20*ZZ(2)+k24*ZZ(2)+k25*ZZ(3)+k26+k27*ZZ(2)+k28+k29*ZZ(2)+k30*ZZ(2)+k31*ZZ(2)+k32*ZZ(2)+k33*ZZ(2)+k34+k36*ZZ(2)+k37*ZZ(2)+k39*ZZ(2)+k50*ZZ(2)+k60*ZZ(2)+k5*ZZ(2)+k67*ZZ(2)+k72*ZZ(2)+k75*ZZ(2)+k82*ZZ(2)+k85*ZZ(2)+k87*ZZ(2)+k93*ZZ(2)+k97*ZZ(2)
                c44=k6*ZZ(2)+k7*ZZ(2)+k10*ZZ(3)+k11+k12+k14+k15-k16+k17+k20*ZZ(2)+k21*ZZ(2)+k28+k29+k30*ZZ(2)+k31+k33*ZZ(2)+k36*ZZ(2)+k37*ZZ(2)+k39*ZZ(2)+k50*ZZ(2)+k51*ZZ(2)+k67*ZZ(2)+k69*ZZ(2)+k70*ZZ(2)+k72*ZZ(2)+k82*ZZ(2)+k85*ZZ(2)+k87*ZZ(2)+k88*ZZ(2)+k89*ZZ(2)+k90*ZZ(2)
                c45=k2*ZZ(2) +k3*ZZ(2) +k6+k7*ZZ(2)+k8*ZZ(2)+k9*ZZ(2)+k11+k14+k18+k19+k22+k25+k26+k27+k29+k30+k31+k32+k33+k39*ZZ(2)+k50+k57+k58+k59+k60+k61+k62+k67+k72+k74+k75*ZZ(2)+k82*ZZ(2)+k83+k85*ZZ(2)+k86+k87*ZZ(2)+k93+k95+k97
                c51=k3 +k6+k7*ZZ(2)+k8+k9+k10+k11+k20+k21+k28+k29+k30+k31+k33+k36+k37+k39*ZZ(2)+k43+k50+k51+k56+k66+k67+k69+k70+k72+k75+k77+k78+k82*ZZ(2)+k83+k85+k87*ZZ(2)+k88+k89+k90+k95
                c52=k7+k10+k39+k40+k43+k45+k47+k66+k68+k71+k82+k87+k89
                c53=k2 +k6+k7+k11+k16+k19+k20+k24+k25+k27+k28+k29+k30+k31+k32+k33+k34+k36+k37+k39+k50+k60+k5+k67+k72+k75+k82+k85+k87+k93+k97
                c54=k6+k7+k10+k11+k20+k21+k28+k29+k30+k31+k33+k36+k37+k39+k50+k51+k67+k69+k70+k72+k82+k85+k87+k88+k89+k90
                c55=k2 +k3 +k7+k8+k9+k39+k74+k75+k82+k83+k85+k86+k87+k95

                c = [[c11.entries[0][0], c12.entries[0][0], c13.entries[0][0], c14.entries[0][0], c15.entries[0][0]],
                    [c21.entries[0][0], c22.entries[0][0], c23.entries[0][0], c24.entries[0][0], c25.entries[0][0]],
                    [c31.entries[0][0], c32.entries[0][0], c33.entries[0][0], c34.entries[0][0], c35.entries[0][0]],
                    [c41.entries[0][0], c42.entries[0][0], c43.entries[0][0], c44.entries[0][0], c45.entries[0][0]],
                    [c51.entries[0][0], c52.entries[0][0], c53.entries[0][0], c54.entries[0][0], c55.entries[0][0]]]

                return _MatrixGenericData(
                    base_ring=a.base_ring,
                    nrows=len(c),
                    ncols=len(c[0]),
                    entries=c
                )

            # Splitting matrix into easily parellizible tasks
            new_matrix = kauers_moosbauer(self, other)
            return new_matrix


        # STRASSEN
        elif self.ncols % 2 == 0 and self.nrows % 2 == 0 and other.ncols % 2 == 0 and other.nrows % 2 == 0 and self.size() == other.size():
            def strassen(A, B):
                n = len(A.entries)
                if n <= 2:  # Base case
                    mid = n // 2

                    A11 = A[:mid, :mid]
                    A12 = A[:mid, mid:]
                    A21 = A[mid:, :mid]
                    A22 = A[mid:, mid:]
                    B11 = B[:mid, :mid]
                    B12 = B[:mid, mid:]
                    B21 = B[mid:, :mid]
                    B22 = B[mid:, mid:]
                    P1 = A11 * (B12 - B22)
                    P2 = (A11 + A12) * (B22)
                    P3 = (A21 + A22) * B11
                    P4 = A22 * (B21 - B11)
                    P5 = (A11 + A22) * (B11 + B22)
                    P6 = (A12 - A22) * ( B21 + B22)
                    P7 = (A11 - A21) * (B11 + B12)

                    C11 = P5 + P4 - P2 + P6
                    C12 = P1 + P2
                    C21 = P3 + P4
                    C22 = P5 + P1 - P3 - P7
                    
                    C = [[C11.entries[0][0], C12.entries[0][0]], [C21.entries[0][0], C22.entries[0][0]]]
                    return B.__class__(
                        base_ring=B.base_ring,
                        nrows=A.nrows,
                        ncols=A.ncols,
                        entries=C,
                    )

                # Partitions
                mid = n // 2
                A11 = A[:mid, :mid]
                A12 = A[:mid, mid:]
                A21 = A[mid:, :mid]
                A22 = A[mid:, mid:]
                B11 = B[:mid, :mid]
                B12 = B[:mid, mid:]
                B21 = B[mid:, :mid]
                B22 = B[mid:, mid:]
                # Recursions
                P1 = strassen(A11, B12 - B22)
                P2 = strassen(A11 + A12, B22)
                P3 = strassen(A21 + A22, B11)
                P4 = strassen(A22, B21 - B11)
                P5 = strassen(A11 + A22, B11 + B22)
                P6 = strassen(A12 - A22, B21 + B22)
                P7 = strassen(A11 - A21, B11 + B12)

                # Combine results to form C
                C11 = (P5 + P4 - P2 + P6).entries
                C12 = (P1 + P2).entries
                C21 = (P3 + P4).entries
                C22 = (P5 + P1 - P3 - P7).entries

                C = []
                lc = len(C11)
                lc0 = len(C11[0])
                for i in range(lc):
                    C.append([])
                    for j in range(lc0):
                        C[-1].append(C11[i][j])
                    for j in range(lc0):
                        C[-1].append(C12[i][j])
                for i in range(lc):
                    C.append([])
                    for j in range(lc0):
                        C[-1].append(C21[i][j])
                    for j in range(lc0):
                        C[-1].append(C22[i][j])

                # Combine quadrants to form C
                C = _MatrixGenericData(
                    base_ring=A.base_ring,
                    nrows=len(C),
                    ncols=len(C[0]),
                    entries=C
                )
                return C
            new_matrix = strassen(self, other)
            return new_matrix

        # standard matrix multiplication
        else:
            new_entries = []
            for i in range(self.nrows):
                new_entries.append([])
                for j in range(other.ncols):
                    new_entry = self.entries[i][0] * other.entries[0][j]
                    for k in range(1, self.ncols):
                        new_entry += self.entries[i][k] * other.entries[k][j]
                    new_entries[i].append(new_entry)


        return other.__class__(
            base_ring=other.base_ring,
            nrows=self.nrows,
            ncols=other.ncols,
            entries=new_entries,
        )

    def __sub__(self, other):
        if self.nrows != other.nrows or self.ncols != other.ncols:
            raise ValueError(
                f"Cannot add matrix of size {self.nrows}x{self.ncols} to matrix of size {other.nrows}x{other.ncols}."
            )

        if self.nrows == 0 or other.nrows == 0:
            return other.__class__(
                base_ring=other.base_ring,
                nrows=other.nrows,
                ncols=other.ncols,
                entries=[],
            )

        new_entries = [self.entries[i].copy() for i in range(self.nrows)]
        for i in range(self.nrows):
            for j in range(self.ncols):
                new_entries[i][j] -= other.entries[i][j]

        return other.__class__(
            base_ring=other.base_ring,
            nrows=other.nrows,
            ncols=other.ncols,
            entries=new_entries,
        )

    def __str__(self):
        return str(self.entries)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        # TODO: when should two empty matrices be equal?
        if self.nrows != other.nrows or self.ncols != other.ncols:
            return False
        if self.nrows == 0 or self.ncols == 0:
            return self.base_ring == other.base_ring
        return self.entries == other.entries

    def __getitem__(self, args):
        r, c = args
        entries = self.entries
        new_entries = entries[r]
        new_data = []

        if not isinstance(new_entries[0], list):
            new_data = [new_entries[c]]
            if not isinstance(new_data[0], list):
                new_data = [new_data]
        else:
            for i in new_entries:
                if not isinstance(i[c], list):
                    new_data.append([i[c]])
                else:
                    new_data.append(i[c])

        ncols = len(new_data[0])
        nrows = len(new_data)
        return self.__class__(
            base_ring=self.base_ring,
            nrows=nrows,
            ncols=ncols,
            entries=new_data,
        )



class _MatrixPythonData:
    def __init__(self, *, base_ring, shape, entries):
        self.base_ring = base_ring
        self.shape = shape
        self.entries = entries



class Matrix:
    """
    Base class for matrices of numbers or polynomials.
    Numerical matrices are built off of flint mat objects
    while generic (polynomial) matrices use the 
    _GenericMatrixData format.
    """
    def __init__(
        self,
        *,
        base_ring,
        nrows=None,
        ncols=None,
        entries=None,
        data=None,
    ):
        self.base_ring = base_ring
        if self.base_ring in {ZZ, QQ, RR, CC, RR_py, ZZ_py}:
            self._is_generic = False
            self._is_python = False

            if self.base_ring in {RR_py, ZZ_py}:
                self._is_python = True
        else:
            self._is_generic = True
            self._is_python = False

        # If data is provided, there is a fast constructor.
        if data is not None:
            self.data = data
            if self._is_generic:
                self.nrows = data.nrows
                self.ncols = data.ncols
            else:
                if self._is_python:
                    self.nrows = data.shape[0]
                    self.ncols = data.shape[1]
                else:

                    self.nrows = data.nrows()
                    self.ncols = data.ncols()
        # Else, if the rows or columns are zero, construct an empty matrix of
        # the appropriate size.
        elif nrows == 0 or ncols == 0:
            if nrows is None or ncols is None:
                raise ValueError("Not enough information provided to construct Matrix.")
            self.nrows = nrows
            self.ncols = ncols

            if self.base_ring == ZZ:
                self.data = flint.fmpz_mat(self.nrows, self.ncols)
            elif self.base_ring == QQ:
                self.data = flint.fmpq_mat(self.nrows, self.ncols)
            elif self.base_ring == RR:
                self.data = flint.arb_mat(self.nrows, self.ncols)
            elif self.base_ring == CC:
                self.data = flint.acb_mat(self.nrows, self.ncols)
            elif self.base_ring == ZZ_py:
                self.data = np.ndarray([self.nrows, self.ncols])
            elif self.base_ring == RR_py:
                self.data = np.ndarray([self.nrows, self.ncols])


            else:
                self.data = _MatrixGenericData(
                    base_ring=self.base_ring,
                    nrows=self.nrows,
                    ncols=self.ncols,
                    entries=[],
                )
        else:
            if (nrows is None or ncols is None) and entries is None:
                raise ValueError("Not enough information provided to construct Matrix.")

            if entries is None:
                self.nrows = nrows
                self.ncols = ncols
                new_entries = self._zero_entries()
            else:
                if isinstance(entries, dict):
                    if nrows is None or ncols is None:
                        raise ValueError(
                            "Not enough information to determine the size of the desired matrix."
                        )
                    self.nrows = nrows
                    self.ncols = ncols
                    new_entries = self._zero_entries()
                    # COERCE the given entries into the base ring, or it's data
                    # class in the non-generic case.
                    for t, e in entries.items():
                        if self._is_generic:
                            new_entries[t[0]][t[1]] = base_ring(e)
                        else:
                            new_entries[t[0]][t[1]] = base_ring(e).data

                # Else, assume that it is a list of lists of base_ring elements.
                else:
                    self.nrows = len(entries)
                    if self.nrows == 0:
                        self.ncols = 0
                    else:
                        self.ncols = len(entries[0])
                    # Now, COERCE the given entries into the base_ring.
                    if self._is_generic:
                        new_entries = []
                        for i in range(self.nrows):
                            new_entries.append([])
                            for j in range(self.ncols):
                                new_entries[i].append(self.base_ring(entries[i][j]))
                    # Unless, the base_ring is special, in which case we use
                    # the underlying data to later construct a FLINT matrix.
                    else:
                        new_entries = []
                        for i in range(self.nrows):
                            new_entries.append([])
                            for j in range(self.ncols):
                                # TODO: fix the coercion here. If the input
                                # entries happen to already be NUTHATCH
                                # classes wrapping FLINT classes, then the
                                # following code does not work!
                                new_entries[i].append(
                                    self.base_ring(entries[i][j]).data
                                )

            if self.base_ring == ZZ:
                self.data = flint.fmpz_mat(new_entries)
            elif self.base_ring == QQ:
                self.data = flint.fmpq_mat(new_entries)
            elif self.base_ring == RR:
                self.data = flint.arb_mat(new_entries)
            elif self.base_ring == CC:
                self.data = flint.acb_mat(new_entries)
            elif self.base_ring == ZZ_py:
                self.data = np.array(new_entries)
            elif self.base_ring == RR_py:
                self.data = np.array(new_entries)

            else:
                self.data = _MatrixGenericData(
                    base_ring=self.base_ring,
                    nrows=self.nrows,
                    ncols=self.ncols,
                    entries=new_entries,
                )

    # Constructor helper function.
    def _zero_entries(self):
        entries = []
        for i in range(self.nrows):
            i_row = []
            for j in range(self.ncols):
                if self._is_generic:
                    i_row.append(self.base_ring(0))
                else:
                    i_row.append(self.base_ring(0).data)
            entries.append(i_row)
        return entries

    def det(self):
        """Alias for `determinenant` method."""
        return self.determinant()

    def determinant(self):
        """Determinant method. Uses flint if generic, if not, see _GenericMatrixData.determinant()"""
        if self.nrows != self.ncols:
            raise ValueError("Matrix must be square.")
        if not self._is_generic:
            if self._is_python:
                return self.base_ring(np.linalg.det(self.data))
            return self.base_ring(self.data.det())
        return self.data.det()
    
    def T(self):
        """Alias for the transpose() method."""
        return self.transpose()
    
    def transpose(self):
        """Returns a transposed copy of self."""
        return Matrix(base_ring=self.base_ring,
                      nrows=self.ncols,
                      ncols=self.nrows,
                      data=self.data.transpose())
    
    def size(self):
        """Returns size of a matrix as a tuple (rows, cols)."""
        return (self.nrows, self.ncols)

    def __add__(self, other):
        """Returns self + other with base ring that of other."""
        return other.__class__(
            base_ring=other.base_ring,
            data=self.data + other.data,
        )

    def __sub__(self, other):
        """Returns self - other with base ring that of other."""
        return other.__class__(
            base_ring=other.base_ring,
            data=self.data - other.data,
        )

    def __mul__(self, other):
        if str(type(other)) == "<class 'nuthatch.rings.reals.RealNumber'>" or str(type(other)) == "<class 'nuthatch.rings.complexes.ComplexNumber'>" or str(type(other)) == "<class 'nuthatch.rings.integers.Integer'>" or str(type(other)) == "<class 'nuthatch.rings.rationals.Rational'>":
            if str(self.base_ring) == "<class 'nuthatch.rings.polynomials.PolynomialRing'>":
                return self.__class__(
                    base_ring=self.base_ring,
                    nrows=self.nrows,
                    ncols=self.ncols,
                    data=self.data * other
                )

            return self.__class__(
                base_ring=self.base_ring,
                nrows=self.nrows,
                ncols=self.ncols,
                data=self.data * other.data
            )
        
        if self.ncols != other.nrows:
            raise ValueError(
                f"Cannot multiply matrix of size {self.nrows}x{self.ncols} with matrix of size {other.nrows}x{other.ncols}."
            )
        
        if self.nrows == 0 or other.nrows == 0 or other.ncols == 0:
            return other.__class__(
                base_ring=other.base_ring,
                nrows=self.nrows,
                ncols=other.ncols,
            )
        ss = self.size()
        so = other.size()
        # if ss == so and ss[0] % 2 == 0 and ss[1] % 2 == 0:
        #     mid = self.size()[0] // 2
        #     A11 = self[:mid, :mid].data
        #     A12 = self[:mid, mid:].data
        #     A21 = self[mid:, :mid].data
        #     A22 = self[mid:, mid:].data
        #     B11 = other[:mid, :mid].data
        #     B12 = other[:mid, mid:].data
        #     B21 = other[mid:, :mid].data
        #     B22 = other[mid:, mid:].data
        #     P1 = A11 * (B12 - B22)
        #     P2 = (A11 + A12) * (B22)
        #     P3 = (A21 + A22) * B11
        #     P4 = A22 * (B21 - B11)
        #     P5 = (A11 + A22) * (B11 + B22)
        #     P6 = (A12 - A22) * ( B21 + B22)
        #     P7 = (A11 - A21) * (B11 + B12)

        #     C11 = P5 + P4 - P2 + P6
        #     C12 = P1 + P2
        #     C21 = P3 + P4
        #     C22 = P5 + P1 - P3 - P7

        #     C11 = C11.tolist()
        #     C12 = C12.tolist()
        #     C21 = C21.tolist()
        #     C22 = C22.tolist()

        #     C1 = np.concatenate((C11, C12), axis=1)
        #     C2 = np.concatenate((C21, C22), axis=1)
        #     C = np.concatenate((C1, C2), axis=0)

        #     C = flint.arb_mat(C.tolist())
        #     return other.__class__(
        #         base_ring=self.base_ring,
        #         nrows=self.nrows,
        #         ncols=other.ncols,
        #         data=C
        #         )
        


        return other.__class__(
            base_ring=other.base_ring,
            data=self.data * other.data,
        )

    def __str__(self):
        return self.data.__str__()

    def __repr__(self):
        return self.data.__repr__()

    def __eq__(self, other):
        if self._is_python and other._is_python:
            return self.data.all() == other.data.all()
        return self.data == other.data
    

    def __getitem__(self, args):
        if type(args) == int or type(args) == slice:
            return ValueError (
                f'The matrix slice method takes 2 args [rows, columns], but 1 were given.'
            )

        elif len(args) == 2:
            r, c = args
            if not self._is_generic:
                r, c = args
                entries = self.data.tolist()
            else:
                entries = self.data.entries
            new_entries = entries[r]
            new_data = []

            if not isinstance(new_entries[0], list):
                new_data = [new_entries[c]]
                if not isinstance(new_data[0], list):
                    new_data = [new_data]
            else:
                for i in new_entries:
                    if not isinstance(i[c], list):
                        new_data.append([i[c]])
                    else:
                        new_data.append(i[c])

            ncols = len(new_data[0])
            nrows = len(new_data)
            if not self._is_generic:
                r, c = args
                entries = self.data.tolist()
                return self.__class__(
                base_ring=self.base_ring,
                nrows=nrows,
                ncols=ncols,
                entries=new_data,
                )

            return self.__class__(
                base_ring=self.base_ring,
                nrows=nrows,
                ncols=ncols,
                entries=new_data,
                data=_MatrixGenericData(
                    base_ring=self.base_ring,
                    nrows=nrows,
                    ncols=ncols,
                    entries=new_data,
                    )
                )

        return ValueError (
            f'The matrix slice method takes 2 args [rows, columns], but {len(args)} were given.'
        )


# Functions for matricies

def generate(value, nrows, ncols):
    """
    Generates a repeating matrix of size (rows, cols). Can be generic or not.
    """
    entries = []
    for i in range(nrows):
        entries.append([])
        for j in range(ncols):
            entries[-1].append(value)
    return Matrix(
                base_ring=value.ring,
                nrows=nrows,
                ncols=ncols,
                entries=entries,
                data=_MatrixGenericData(
                    base_ring=value.ring,
                    nrows=nrows,
                    ncols=ncols,
                    entries=entries,
                    )
                )

def random(base_ring, max_val, nrows, ncols, poly=0, ngens=1):
    """
    Generates a random matrix of size nrows x ncols using np.random.
    Values go from [0, max_val]. 
    Polynomials are of form n0*x0**r0 + n1*x1**r1 ....
    """
    entries = []

    if poly:
        p = PolynomialRing(base_ring=base_ring, ngens=ngens, prefix="x", packed=False)
        for i in range(nrows):
            entries.append([])
            for j in range(ncols):
                val = p(0)
                for k in range(ngens):
                    val = val + base_ring(np.random.random() * max_val) * p.gens[k] ** ZZ(int(np.random.random() * max_val))
                entries[-1].append(val)
                

        return Matrix(
                base_ring=PolynomialRing,
                nrows=nrows,
                ncols=ncols,
                entries=entries,
                data=_MatrixGenericData(
                    base_ring=PolynomialRing,
                    nrows=nrows,
                    ncols=ncols,
                    entries=entries,
                    )
                )
    if base_ring == RR:
        vals = np.random.random(nrows * ncols) * max_val
        for i in range(nrows):
            entries.append([])
            for j in range(ncols):
                entries[-1].append(base_ring(vals[(i + 1) * (j)]))

    elif base_ring == ZZ:
        vals = np.random.randint(max_val, size=nrows * ncols)
        for i in range(nrows):
            entries.append([])
            for j in range(ncols):
                entries[-1].append(base_ring(int(vals[(i + 1) * (j)])))
    
    elif base_ring == RR_py:
        vals = np.random.random(nrows * ncols) * max_val
        for i in range(nrows):
            entries.append([])
            for j in range(ncols):
                entries[-1].append(vals[(i + 1) * (j)])
        
        return Matrix(
                    base_ring=base_ring,
                    entries=entries,
                    )
    elif base_ring == ZZ_py:
        vals = np.random.randint(max_val, size=nrows * ncols)
        for i in range(nrows):
            entries.append([])
            for j in range(ncols):
                entries[-1].append(int(vals[(i + 1) * (j)]))
    
        return Matrix(
                    base_ring=base_ring,
                    entries=entries,
                    )

    return Matrix(
                base_ring=base_ring,
                nrows=nrows,
                ncols=ncols,
                entries=entries,
                )
