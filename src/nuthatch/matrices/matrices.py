"""
MATRICES

A module for dense matrices. We create a single abstract matrix class which
holds a `data` attribute. We assume that this points to a class on which all
matrix operations can be performed and we overload matrices in that way. We
provide a generic implementation of these operations as `_MatrixGeneric`. Special
examples are provided by `FLINT`.

AUTHORS:
- Benjamin Antieau (2024): initial version.
"""

import flint
import numpy as np
from nuthatch.rings.integers import ZZ
from nuthatch.rings.reals import RR, RR_py
from nuthatch.rings.complexes import CC

from nuthatch.rings.rationals import QQ


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
            raise ValueError("matrix must be square")
        if self.nrows == 0:
            return self.base_ring(0)
        if self.nrows == 1:
            return self.entries[0][0]
        # TODO: implement cofactor definition of determinant for starters.

        return NotImplemented

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
        if self.ncols % 5 == 0 and self.nrows % 5 == 0:
            print('kmb')
            def kauers_moosbauer(a, b):
                la = a.nrows
                lb = b.ncols
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

                k1 =(-a51+2*a52+2*a54-2*a55)*(3*b11+b21+b41)
                k2 =(-1*a53)*(2*b12-2*b14+4*b15-b32+b34-2*b35-b52+b54-2*b55)
                k3 =(a51-2*a52-a54+2*a55)*(2*b11+b41)
                k4 =(-1*a53+a55)*(b21+b51)
                k6=(a41-a42-2*a44+2*a45)*(-2*b11+2*b12-b21+b22)
                k7=a44*(3*b12+b22+b42)
                k8=(-1*a44-a51+2*a52+2*a54-2*a55)*(2*b11+b12+b22+b41)
                k9=(-1*a44+a54)*(-2*b11+2*b12-b41+b42)
                k10=(-1*a41+2*a42+a44-2*a45)*(b12+b22)
                k11=(-1*a45)*(2*b11-2*b12-b31+b32-b51+b52)
                k12=(-1*a43+a45)*(-1*b11+b12-b21+b22-b31+b32)
                k13=(a41-2*a42+a43-2*a44+a45-a51+2*a52+2*a54-2*a55)*(-1*b11+b12-b21+b22)
                k14=(-1*a31+2*a32+2*a34-2*a35)*(3*b14-6*b15+b24-2*b25+b44-2*b45)
                k15=(a31-2*a32-2*a34+2*a35-a41+2*a42+2*a44-2*a45)*(-1*b11+b14-2*b15-b21+b24-2*b25-b31+b32+b34-2*b35)
                k16=(-1*a31+2*a32+2*a34-2*a35-a44)*(b12+2*b14-4*b15+b22+b44-2*b45)
                k17=(-1*a31+2*a32+2*a34-2*a35+a41-2*a42+a43-2*a44+a45)*(-1*b11+b12-b21+b22-b31+b32+b34-2*b35)
                k18=(-1*a31+a32+2*a34-2*a35+a51-a52-2*a54+2*a55)*(b22+b52)
                k19=(-1*a31+a32+2*a34-2*a35+a41-a42-2*a44+2*a45)*(b21-b22+b24-2*b25+b31-b32+b34-2*b35+b51-b52+b54-2*b55)
                k20=(a33-a53)*(-2*b12+b32+b52)
                k21=(-1*a33+a43+a53)*b32
                k22=(-1*a34+a54)*(b12-b14+b15-b42+b44-b45+b52-b54+b55)
                k23=(-1*a34+a44)*(-2*b11+2*b12-b41+b42)
                k24=(-a34+a44)*(-2*b12+2*b14-4*b15-b42+b44-2*b45)
                k25=(a31-2*a32-a34+2*a35)*(2*b14-4*b15+b44-2*b45)
                k26=a35*(b24-2*b25+b34-2*b35+b54-2*b55)
                k27=(a35-a45-a53)*(b21+b31-b32+b51)
                k28=(-a31+a32+2*a34-a35)*(2*b14-4*b15+b24-2*b25)
                k29=(a31-a32-2*a34+a35-a41+a42+2*a44-2*a45)*(-2*b11+2*b12+b24-2*b25+b31-b32+b34-2*b35+b51-b52+b54-2*b55)
                k30=(a31-a32-2*a34+a35+a45-a51+a52+2*a54-a55)*(2*b11-2*b12+b21+b52)
                k31=(-a31+a32+2*a34-a35+a41-a42-2*a44+a45)*(-2*b11+2*b12-2*b14+4*b15+b31-b32+b34-2*b35+b51-b52+b54-2*b55)
                k32=(a31-a32-2*a34+a35-a41+a42+2*a44-a45+a53)*(2*b12-2*b14+4*b15+b21+b31-b32+b34-2*b35+b51-b52+b54-2*b55)
                k33=(-a31+a32+2*a34-a35+a41-a42-2*a44+a45+a51-a52-2*a54+a55)*(2*b11+b21)
                k34=(-a33+a35)*(b34-2*b35)
                k35=(-a31+2*a32-a33+2*a34-a35+a41-2*a42+a43-2*a44+a45)*(-b31+b32+b34-2*b35)
                k36=(-a31+a32-a33+2*a34-a35+a51-a52+a53-2*a54+a55)*(-2*b12+b52)
                k37=(-a34-a35+a54+a55)*(-2*b12+2*b14-2*b15+b32-b34+b35+b52-b54+b55)
                k38=a23*(b21-b22-b23+b25+b31-b32-b33+b35)
                k39=(-a22+a23)*(-b11+b13+b21-b23-b41+b43)
                k40=(a22-a23+2*a24+2*a25+a41-2*a42-4*a44+2*a45)*(b11-b12-b13+b15+b21-b22-b23+b25+b31-2*b32-b33+b35)
                k41=(a22-a23-a42+a43)*(-2*b13+b33+b53)
                k42=(a22-a23+2*a24+2*a25+a41-2*a42+a43-4*a44+a45)*(-2*b12-2*b13+2*b14-2*b15-2*b21+2*b22+2*b23-2*b25-b31+3*b32+2*b33-b34+b41-b42-b43+b45+b52+b53-b54+b55)
                k43=(-a24+a44)*(b15-b45+b55)
                k44=-(a25*(-b11+2*b12+b13-2*b15+b41-2*b42-b43+2*b45-b51+2*b52+b53-2*b55))
                k45=(a21-2*a22-4*a24-5*a25-a41+2*a42+4*a44-2*a45)*(b11-2*b12-b13+b15+b21-2*b22-b23+b25+b31-2*b32-b33+b35)
                k46=(-a21+a22+2*a24+3*a25)*(2*b11-4*b12-2*b13+2*b15+b21-2*b22-b23+b25+b31-2*b32-b33+b35)
                k47=(a24+a25-a44)*(b15-2*b22+b42-b45+b55)
                k48=(-a23+a24+a25+a43-a44-a45)*(3*b12-3*b14+3*b15+b22-b24+b25+b32-b34+b35)
                k49=(a22-a23+a24+a25)*(-2*b11+2*b13-b41+b43)
                k50=(a22-a23+a24+a25+a41-2*a42+a43-3*a44+a45)*(-4*b12-4*b13+4*b14-4*b15-2*b21+b22+b23+b24-3*b25+b32+b33-b34+b35+b41-b42-b43+b45+b52+b53-b54+b55)
                k51=(-a22+a23-a24-a25+a42-a43+a44+a45)*(b23+b33+b53)
                k52=(a11-2*a12+2*a14+2*a15-a41+2*a42-2*a44-2*a45-a51+2*a52-2*a54-2*a55)*(3*b12+3*b13-3*b14+4*b15+b22+b23-b24+b25+b32+b33-b34+b35-b43-b45+b53+b55)
                k53=(a11-2*a12+2*a14+2*a15-a21+2*a22-2*a24-2*a25-a51+2*a52-2*a54-2*a55)*b13
                k54=(-a11+2*a12-2*a14-2*a15+a21-2*a22+2*a24+2*a25+a31-2*a32+2*a34+2*a35)*(3*b12+3*b13-3*b14+3*b15+b22+b23-b24+b25-b42-b43+b44-b45+b52+b53-b54+b55)
                k55=(a11-2*a12+2*a14+2*a15-a21+2*a22-2*a24-2*a25-a31+2*a32-3*a34-2*a35+a54)*(b12-b14+b15-b42-b43+b44-b45+b52+b53-b54+b55)
                k56=(a11-2*a12+2*a14+2*a15+a24-a41+2*a42-3*a44-2*a45-a51+2*a52-2*a54-2*a55)*(b15-b43-b45+b53+b55)
                k57=(a11-2*a12+2*a14+2*a15-a23+a24+a25-a41+2*a42+a43-3*a44-3*a45-a51+2*a52-2*a54-2*a55)*(3*b12+2*b13-3*b14+3*b15+b22+b23-b24+b25+b32+b33-b34+b35)
                k58=(-a11+2*a12-2*a14-2*a15+a23-a24-a25+a31-2*a32+3*a34+2*a35+a41-2*a42-a43+3*a44+3*a45-a54)*(b12-b14+b15)
                k59=(-a11+2*a12-2*a14-2*a15+a23-a24-a25+a33-a34-a35+a41-2*a42-a43+3*a44+3*a45+a51-2*a52-a53+3*a54+3*a55)*(2*b12+2*b13-2*b14+2*b15+b22+b23-b24+b25+b32+b33-b34+b35)
                k60=(a12-a22+a31-2*a32-2*a34+2*a35-a41+a42+2*a44-2*a45)*(-5*b12+4*b14-3*b15+b24-2*b25-b42+b45+b52-b54+b55)
                k61=(-a12+a23-a24-a25+a42-a43+a44+a45+a52)*(2*b13+b23+b33)
                k5=(a12-a23+a24+a25-a32-a42+a43-a44-a45)*(2*b12-2*b14+2*b15+b22-b24+b25)
                k62=(a12-a23+a24+a25-a33+a34+a35-a42+a43-a44-a45-a52+a53-a54-a55)*(2*b12+2*b13-2*b14+2*b15+b22+b23-b24+b25+b33)
                k63=(a11-a12-2*a14-2*a15-a21+a22+2*a24+2*a25-a41+a42+2*a44+2*a45)*(-b12+b15+b42-b45-b52+b55)
                k64=(a11-a12-2*a14-2*a15+a23-a41+a42+2*a44+2*a45)*(b11-b12-b13+b15)
                k65=(-a11+a12+2*a14+2*a15+a21-a22-2*a24-3*a25+a41-a42-2*a44-2*a45)*(b11-2*b12-b13+2*b15+b42-b45-b52+b55)
                k66=(a13-a43)*(2*b15+b25+b35)
                k67=(-a13-a22+2*a23-2*a24-2*a25-a41+2*a42-a43+4*a44)*(-3*b12-4*b13+2*b14-b15-2*b21+b22+b23-b25+b32+b33-b34+b35+b41-b42-b43+b45+b52+b53-b54+b55)
                k68=(a13-a23-a43)*(-b11+2*b12+b13-2*b15-b21+2*b22+b23-2*b25-b31+2*b32+b33-2*b35)
                k69=(-a13+a23+a43)*(-b12+b15-b22+b25-b32+b35)
                k70=(a13-a23-a45)*(-4*b12+2*b14-b34+2*b35-b42+b45+b52-b54+b55)
                k71=(a13-a21+a22+2*a24+3*a25-a43)*(b11-2*b12-b13+b21-2*b22-b23+b25+b31-2*b32-b33+b35)
                k72=(a13-a23+a24+a25-a44-a45)*(-4*b13-2*b21+b23+b33+b41-b43+b53)
                k73=(a12+a13-2*a23+2*a24+2*a25-a42+a43-2*a44-2*a45+a51-2*a52-2*a54+2*a55)*(3*b11-3*b13+b21+b41-b43+b53)
                k74=(-a12+a13+a42-a43+a52)*(2*b15+b25)
                k75=(a12-a13-a22+a23-a42+a43)*(-b12+b15-b22+b25)
                k76=(a14-a44-a54)*(-b42+b44-b45+b52-b54+b55)
                k77=(a14-a24-a54)*(-b43+b53)
                k78=(-a11+2*a12-3*a14-2*a15+a41-2*a42+3*a44+2*a45+a51-2*a52+3*a54+2*a55)*(-b43-b45+b53+b55)
                k79=(-a11+2*a12-3*a14-2*a15+a21-2*a22+3*a24+2*a25+a31-2*a32+3*a34+2*a35)*(-b42-b43+b44-b45+b52+b53-b54+b55)
                k80=(a15-a45)*(-b41+b42+b43-b45+b51-b52-b53+b55)
                k81=(a11-a12-2*a14-3*a15-a21+a22+2*a24+3*a25-a41+a42+2*a44+3*a45)*(b42-b45-b52+b55)
                k82=(a14+a15-a44-a45)*(b11-b12-b13+b15-b21+b22+b23-b25+b41-b42-b43+b45)
                k83=(a14+a15-a44-a45-a54)*(3*b12-3*b15+b22+b42-b45+b55)
                k84=(-a14-a15+a44+a45+a54+a55)*(-2*b13+b33+b53)
                k85=(a14+a15-a44-a45-a54-a55)*(-2*b12+2*b14-2*b15+b32-b34+b35+b52-b54+b55)
                k86=(a14+a15-a44-a45+a53-a54-a55)*(2*b12-2*b14+4*b15-b32+b34-b35-b52+b54-2*b55)
                k87=(a14+a15+a22-a23-a44-a45)*(-b11-b12+b13+b15+b21-b22-b23+b25-b41+b43)
                k88=(a14+a15+a22-a23-a34-a35-a42+a43-a44-a45)*(2*b12-2*b14+2*b15+b22+b23-b24+b25+b33+b53)
                k89=(-a14-a15+a24+a25+a44+a45)*(-2*b12+2*b15-b42+b45)
                k90=(-a14-a15+a24+a25+a34+a35)*(b22+b23-b24+b25+b32+b33-b34+b35+b52+b53-b54+b55)
                k91=(a11-2*a12+3*a14+3*a15-a21+3*a22-a23-2*a24-2*a25-a31+3*a32-a33-2*a34-2*a35-a42+a43-a44-a45-a52+a53-a54-a55)*(2*b12+2*b13-2*b14+2*b15+b22+b23-b24+b25)
                k92=(-a12-2*a13+a14+a15+3*a23-3*a24-3*a25+a42-a43+2*a44+2*a45-a51+2*a52+a54-2*a55)*(2*b11-4*b13+b41-b43+b53)
                k93=(a13-a14-a15-a23+a24+a25-a33+a34+a35)*(b32-b34+b35)
                k94=(-a12+a13-a14-a15+a42-a43+a44+a45)*(b11-b12-b13+b15+b21-b22-b23+b25)
                k95=(a12-a13-a14-a15-a42+a43+a44+a45+a51-2*a52-a54+2*a55)*(b12+b15+b22)
                k96=(a12-a13+a14+a15-a42+a43-a44-a45-a52+a53-a54-a55)*b33
                k97=(-a12+a13-a14-a15+a22-a23+a24+a25-a31+2*a32+3*a34-2*a35+a41-a42-2*a44+2*a45)*(-4*b12+2*b14-b42+b45+b52-b54+b55)


                c11=2*k1 +k3 -k8+k10+k11+k12+k13+k20+k21+k30+k33+k36+2*k38-k43+k46+k48+k52-k56-k57+2*k64+k65-k66-2*k68+2*k69-k71+k72+k73-k77+k80+k81+k92
                c12=k1 -k8+k10+k12+k13+2*k38+k40+k42-k43+k44+k46+k48+k52-k53-k56-k57+2*k63+2*k64+2*k65-k66+k67-k68+k69+k70-k71+k72
                c13=k1 +k2 -k8+k11+k13+k16-k17+k20+k23+k24+k25+k27+k30+k31+k32+k33-k34+k35+k36
                c14=k1 -k8+k10+k11+k12+k13+k20+k21+k30+k33+k36
                c15=k1 +k2 +k3 +k4 +k31+k32+k33
                c21=k3 +k6+2*k7+k8+k9+k10+k20+k21+k30+k33+k36+k38+k43+k46+k56+k64+k65-k68+k69-k71+k77+k78+k81+k83+k95
                c22=k7+k10+k38+k40+k43+k45+k46+k47+k63+k64+k65
                c23=k2 +k6+k7+k11+k16+k19+k20+k24+k25+k26+k27+k29+k30+k31+k32+k33+k36
                c24=k6+k7+k10+k20+k21+k30+k33+k36
                c25=k2 +k3 +k6+k7+k8+k9+k11+k18+k19+k26+k27+k29+k30+k31+k32+k33
                c31=2*k1 +k3 -k8+k10+k12+k13+2*k39+k40+k41+k42-k43+k48+k49+k51+k52-k56-k57+k67+k69+k70+2*k72+k73+k75-k77+k87+k92+k94
                c32=k1-k8+k10+k12+k13+k39+k40+k42-k43+k48+k49+k52-k53-k56-k57+k67+k69+k70+k72
                c33=k1 +k3 -k22+k39+k41+k53+k54+k55+k5-k62+k72+k73+k75+k87-k88+k91+k92+k94+k96
                c34=k1 -k8+k10+k12+k13+k39+k40+k41+k42+k49+k51+k67+k69+k70+k72
                c35=k1 +k3 +k39+k61+k72+k73+k75+k84+k87+k92+k94+k96
                c41=2*k3 +2*k6+4*k7+2*k8+2*k9+3*k10+k11+k12+2*k14+k15-k16+k17+2*k20+2*k21+k22+k25+k28+k29+2*k30+k31+2*k33+2*k36+2*k37+k38+4*k39+2*k43+k46+k48+3*k50+2*k51+2*k56+k58+k60+k64+k65+k66+3*k67-k68+3*k69+2*k70-k71+3*k72+2*k75+k76+2*k77+2*k78+k81+4*k82+2*k83+2*k85+4*k87+2*k88+2*k89+2*k90+2*k95+k97
                c42=2*k7+3*k10+k12+k14+k15-k16+k17+k38+2*k39+2*k40+2*k43+2*k45+k46+2*k47+k48+k50+k55+k58+k63+k64+k65+k66+k67+k68+k69+k70+k71+k72+k76+k77+k79+2*k82+2*k87+2*k89
                c43=2*k2 +2*k6+2*k7+2*k11+k14+2*k16+2*k19+2*k20+2*k24+3*k25+k26+2*k27+k28+2*k29+2*k30+2*k31+2*k32+2*k33+k34+2*k36+2*k37+2*k39+2*k50+2*k60+2*k5+2*k67+2*k72+2*k75+2*k82+2*k85+2*k87+2*k93+2*k97
                c44=2*k6+2*k7+3*k10+k11+k12+k14+k15-k16+k17+2*k20+2*k21+k28+k29+2*k30+k31+2*k33+2*k36+2*k37+2*k39+2*k50+2*k51+2*k67+2*k69+2*k70+2*k72+2*k82+2*k85+2*k87+2*k88+2*k89+2*k90
                c45=2*k2 +2*k3 +k6+2*k7+2*k8+2*k9+k11+k14+k18+k19+k22+k25+k26+k27+k29+k30+k31+k32+k33+2*k39+k50+k57+k58+k59+k60+k61+k62+k67+k72+k74+2*k75+2*k82+k83+2*k85+k86+2*k87+k93+k95+k97
                c51=k3 +k6+2*k7+k8+k9+k10+k11+k20+k21+k28+k29+k30+k31+k33+k36+k37+2*k39+k43+k50+k51+k56+k66+k67+k69+k70+k72+k75+k77+k78+2*k82+k83+k85+2*k87+k88+k89+k90+k95
                c52=k7+k10+k39+k40+k43+k45+k47+k66+k68+k71+k82+k87+k89
                c53=k2 +k6+k7+k11+k16+k19+k20+k24+k25+k27+k28+k29+k30+k31+k32+k33+k34+k36+k37+k39+k50+k60+k5+k67+k72+k75+k82+k85+k87+k93+k97
                c54=k6+k7+k10+k11+k20+k21+k28+k29+k30+k31+k33+k36+k37+k39+k50+k51+k67+k69+k70+k72+k82+k85+k87+k88+k89+k90
                c55=k2 +k3 +k7+k8+k9+k39+k74+k75+k82+k83+k85+k86+k87+k95
                c = [[c11, c12, c13, c14, c15], [c21, c22, c23, c24, c25], [c31, c32, c33, c34, c35], [c41, c42, c43, c44, c45], [c51, c52, c53, c54, c55]]
                
                return c
            
            new_entries = kauers_moosbauer(self, other)

        # STRASSEN
        elif self.ncols % 2 == 0 and self.nrows % 2 == 0:
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
                    
                    C = [[C11, C12], [C21, C22]]
                    print(C)
                    return C
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
                C11 = P5 + P4 - P2 + P6
                C12 = P1 + P2
                C21 = P3 + P4
                C22 = P5 + P1 - P3 - P7

                # Combine quadrants to form C
                C = np.vstack((np.hstack((C11, C12)), np.hstack((C21, C22))))
                return C
            new_entries = strassen(self, other)
            print(new_entries)
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
        return self.__class__(
            base_ring=self.base_ring,
            nrows=self.nrows,
            ncols=self.ncols,
            entries=self.entries[r][c],
        )


        
class Matrix:
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
        if self.base_ring in {ZZ, QQ, RR, CC}:
            self._is_generic = False
        else:
            self._is_generic = True

        # If data is provided, there is a fast constructor.
        if data is not None:
            self.data = data
            if self._is_generic:
                self.nrows = data.nrows
                self.ncols = data.ncols
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
        if self.nrows != self.ncols:
            raise ValueError("matrix must be square")
        return self.base_ring(self.data.det())
    
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
        if True:
            return self.__class__(
                base_ring=self.base_ring,
                nrows=self.nrows, 
                ncols=self.ncols,
                data=self.data * flint.arb()
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

        return other.__class__(
            base_ring=other.base_ring,
            data=self.data * other.data,
        )

    def __str__(self):
        return self.data.__str__()

    def __repr__(self):
        return self.data.__repr__()

    def __eq__(self, other):
        return self.data == other.data

    def __getitem__(self, args):
        if len(args) == 2:
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


    # def concat(self, other, axis=0):
    #     """Concatinates matrices along an axis."""
    #     if axis != 0 and axis != 1:
    #         return ValueError (
    #             f'Axis must be 0 (vertical), or 1 (horizontal), a value of {axis} was given.'
    #         )
        
    #     if axis == 0 and self.data.ncols != other.data.ncols:
    #         return ValueError (
    #             'For vertical concatination, the number of columns on both matrices must match.'
    #         )
        
    #     if axis == 1 and self.data.nrows != other.data.nrows:
    #         return ValueError (
    #             'For horizontal concatination, the number of rows on both matrices must match.'
    #         )

    #     s = self.data.entries
    #     o = other.data.entries
    #     if axis == 0:
    #         for row in range(ncols(s)):
    #             s.append(o[row])

    #     else:
    #         len_s = s.nrows
    #         len_o = o.ncols
    #         for i in range(len_s):
    #             for j in range(len_o):
    #                 s[i].append(o[i][j])
                    
    #     nrows = s.nrows
    #     ncols = s.ncols

    #     return self.__class__(
    #         base_ring=self.base_ring,
    #         nrows=nrows,
    #         ncols=ncols,
    #         entries=s,
    #         data=_MatrixGenericData(
    #             base_ring=self.base_ring,
    #             nrows=nrows,
    #             ncols=ncols,
    #             entries=s,
    #             )
    #         )

        

"""Functions for matricies."""