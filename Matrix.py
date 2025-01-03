import copy as cp
from Expressions import *
from IPython.display import display, Latex

Ex_MNotSq = Exception("The matrix isn't square")
Ex_MDim = Exception("Dimensions don't match")

class Matrix:
    def __init__(self, data:list):
        self.data = data
        self.n = len(data)
        self.p = len(data[0])

    def __getitem__(self, idx:(int, int)):
        i, j = idx
        return self.data[i][j]

    def __setitem__(self, idx:(int, int), val):
        i, j = idx
        self.data[i][j] = val

    def __neg__(self):
        return self*-1

    def __iadd__(self, other):
        if type(other) != Matrix:
            return 
        if self.n != other.n or self.p != other.p:
            raise Ex_MDim
        
        for i in range(self.n):
            for j in range(self.p):
                self[i,j] += other[i,j]

        return self

    def __add__(self, other):
        copy = cp.deepcopy(self)
        copy += other
        return copy

    def __sub__(self, other):
        return self + -other

    def __mul__(self, other):
        if isNumber(other):
            copy = cp.deepcopy(self)
            for i in range(copy.n):
                for j in range(copy.p):
                    copy[i,j] *= other
            return copy

        elif type(other) == Matrix:
            if self.p != other.n:
                raise Ex_MDim

            res = Matrix.zero(self.n, other.p)
            for i in range(self.n):
                for j in range(other.p):
                    for k in range(self.p):
                        res[i, j] += self[i,k]*other[k,j]
            return res

    def __rmul__(self, other):
        if isNumber(other):
            return self*other

    def __pow__(self, n:int):
        if self.n != self.p:
            raise Ex_MNotSq
        
        if n == 0:
            return Matrix.id(self.n)
        else:
            return self**(n-1) * self

    def T(self):
        data = self.data
        res = Matrix.zero(self.p, self.n)
        for i in range(self.p):
            for j in range(self.n):
                res[i,j] = self[j, i]
        return res

    def zero(n:int, p:int = None):
        p = n if p == None else p
        mtx = [[Scalar(0)] * p for _ in range(n)]
        return Matrix(mtx)

    def id(n:int, p:int = None):
        p = n if p == None else p
        idm = Matrix.zero(n, p)
        for i in range(min(n, p)):
            idm[i, i] = Scalar(1)
        return idm

    def __str__(self):
        tex = r"\pmatrix{"
        for i in range(self.n):
            for j in range(self.p):
                tex += str(self[i, j]) + " & "
            tex = tex[:-2] + r" \\ "
        tex = tex[:-3] + r"}"
        return tex

    def __repr__(self):
        latex(str(self))
        return ""
        
def Diag(*eigens):
    res = Matrix.id(len(eigens))
    for i in range(res.n):
        res[i, i] *= eigens[i]
    return res

def tr(mtx):
    if mtx.n != mtx.p:
        raise Ex_MNotSq
    
    res = Scalar(0)
    for i in range(mtx.n):
        res += mtx[i, i]
    return res

def displayMul(A, B):
    latex(str(A) + str(B) + " = " + str(A*B))
