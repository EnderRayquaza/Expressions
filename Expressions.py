import copy as cp
from abc import ABC
from IPython.display import display, Latex

def round_(val):
    N = 3 #Change here to increase "accuracy"
    return remFloat(round(val, N))

def latex(text:str):
    text = "$" + text + "$"
    display(Latex(text))

def remFloat(val):
    return int(val) if int(val) == val else val

class Tex(ABC):
    def __repr__(self):
        latex(str(self))
        return ""

class Scalar(Tex):
    def __init__(self, val):
        if type(val) in (int, float):
            self.real = remFloat(val)
            self.img = 0
        elif type(val) == complex:
            self.real = remFloat(val.real)
            self.img = remFloat(val.imag)
        elif type(val) == Scalar:
            self.real = val.real
            self.img = val.img
        else:
            raise TypeError(type(val))

    def __complex__(self):
        return self.real + self.img*1j

    def __neg__(self):
        return Scalar(-complex(self))
    
    def __add__(self, other):
        if isNumber(other):
            return Scalar(complex(self) + other)
        elif isExpr(other):
            return other+self

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + -other

    def __mul__(self, other):
        if isNumber(other):
            return Scalar(complex(self) * other)
        elif isExpr(other):
            return other*self

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if isNumber(other):
            return Scalar(complex(self) / other)

    def __rtruediv__(self, other):
        if isNumber(other):
            return Scalar(other / complex(self))
        elif isExpr(other):
            return other*(1/self)

    def __pow__(self, other):
        if isNumber(other):
            return Scalar(complex(self)**other)

    def __rpow__(self, other):
        if isNumber(other):
            return Scalar(other**complex(self))
        elif isExpr(other):
            return other**self

    def __abs__(self):
        return Scalar(abs(complex(self)))

    def toStr(self, hide=False, exp=False):
        real, img = "", ""
        if self.real == self.img == 0:
            return "0"
            
        if self.real == 0:
            img = str(round_(self.img)) if abs(round_(self.img)) != 1 else str(round_(self.img))[:-1] #The 1 is removed to avoid "1x" or "-1x"
            if exp:
                img += "1" if self.img < 0 else ""
            return img + "i"
        elif self.img == 0:
            if hide:
                real = str(round_(self.real)) if abs(round_(self.real)) != 1 else str(round_(self.real))[:-1] #The 1 is removed to avoid "1x" or "-1x"
            else:
                real = str(round_(self.real))
            return real
        else:
            real = str(round_(self.real))
            img = str(abs(round_(self.img))) if abs(round_(self.img)) != 1 else str(abs(round_(self.img)))[:-1] #The 1 is removed to avoid "1x" or "-1x"
            op = "+" if self.img > 0 else "-"
            img += "i"
            if hide:
                return "(" + real + op + img + ")"
            else:
                return real + op + img
                
    def __str__(self):
        return self.toStr()

    def __eq__(self, other):
        return complex(self) == other

def isNumber(obj):
    return type(obj) in (int, float, complex) or type(obj) == Scalar

def checkScalar(obj):
    if type(obj) in (int, float, complex):
        return Scalar(obj)
    else:
        return obj

class Var(Tex):    
    def __init__(self, name:str):
        self.name = name

    def __deepcopy__(self, memo):
        return self

    def __str__(self):
        return self.name

NoneVar = Var("∅")


class Monom(Tex):
    def __init__(self, var:Var, scalar:Scalar = Scalar(1), exponent:Scalar = Scalar(1)):
        self.var = var
        self.scalar = scalar
        self.exponent = exponent

    def isScalar(self):
        return self.var == NoneVar

    def toScalar(self):
        if not self.isScalar():
            raise Exception("This monom isn't a scalar")
        else:
            return self.scalar

    def __neg__(self):
        return Monom(self.var, self.scalar*-1, self.exponent)

    def __add__(self, other):
        other = checkScalar(other)
        
        if type(other) != Monom: #Polynom manage all additions
            return Polynom(self, other)
        
        if prop(self, other): #Though if the monoms are proportional, the sum can still be a monom
            if self.scalar+other.scalar == 0:
                return Scalar(0)
            else:
                return Monom(self.var, self.scalar+other.scalar, self.exponent)
        else:
            return Polynom(self, other)

    def __radd__(self, other):
        return self+other

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        other = checkScalar(other)
        
        if type(other) == Scalar:
            if other == 0:
                return Scalar(0)
            else:
                return Monom(self.var, self.scalar*other, self.exponent)
        
        if type(other) == Monom:
            if other.isScalar():
                return self*other.toScalar()
                
            if self.var == other.var :
                if self.exponent+other.exponent == 0: #The var cancels so it's only a scalar now
                    return Scalar(self.scalar*other.scalar)
                else:
                    return Monom(self.var, self.scalar*other.scalar, self.exponent+other.exponent)
            else: #Multipling different vars make a multinom
                return Multinom(self, other)
        if isExpr(other):
            return other*self
        else:
            raise TypeError(type(other))

    def __rmul__(self, other):
        return self*other

    def __truediv__(self, other):
        return self * other**(-1)

    def __rtruediv__(self, other):
        return self**(-1) * other

    def __pow__(self, ex:Scalar):
        return Monom(self.var, self.scalar**ex, self.exponent*ex)

    def __str__(self):
        sc = Scalar(self.scalar).toStr(True)
        ex = r"^{" + Scalar(self.exponent).toStr(True, True) + r"}"
        return sc + self.var.name + ex

    def unit(self):
        return Monom(self.var, 1, self.exponent)

ZeroMonom = Monom(NoneVar, 0, 1)
UnitMonom = Monom(NoneVar, 1, 1)

class Multinom(Tex):
    def __init__(self, *monoms):
        self.monoms = []
        self.scalar = Scalar(1)
        for monom in monoms:
            self *= monom

    def __neg__(self):
        return self*-1

    def __add__(self, other):
        if type(other) in (Monom, Multinom) and prop(self, other): #Roughly the same as monoms
            copy = cp.copy(self)
            copy.scalar += other.scalar
            if(copy.scalar == 0):
                return Scalar(0)
            else:
                return copy
        else:
            return Polynom(self, other)

    def __radd__(self, other):
        return self+other

    def __sub__(self, other):
        return self + (-other)

    def __imul__(self, other):
        if isNumber(other):
            self.scalar *= Scalar(other)
            if self.scalar == 0:
                return Scalar(0)
            else:
                return self

        if type(other) == Monom:
            self.scalar *= other.scalar
            
            if self.scalar == 0: #Special cases
                return Scalar(0)
            if other.var == NoneVar:
                print("NoneVar used !")
                return self
                
            for i in range(len(self.monoms)):
                if other.var == self.monoms[i].var:
                    self.monoms[i] *= other.unit()
                    if isNumber(self.monoms[i]):
                        self.scalar *= self.monoms[i]
                        del self.monoms[i]
                    return self
            self.monoms.append(other.unit())

        elif type(other) == Multinom:
            for monom in other.monoms:
                self *= monom
                
        return self

    def __mul__(self, other):
        copy = cp.deepcopy(self)
        copy *= other                
        return copy

    def __rmul__(self, other):
        return self*other

    def __truediv__(self, other):
        return self * other**(-1)

    def __rtruediv__(self, other):
        return self**(-1) * other

    def __pow__(self, ex:Scalar):
        res = cp.copy(self)
        res.scalar **= ex
        for i in range(len(res.monoms)):
            res.monoms[i] **= ex
            if isNumber(res.monoms[i]):
                res.scalar *= res.monoms[i]
                del res.monoms[i]
        return res

    def __str__(self):
        res = self.scalar.toStr(True)
        for monom in self.monoms:
            res += str(monom) + " "
        return res

    def isMonom(self):
        return len(self.monoms) == 1

    def getMonom(self):
        if not self.isMonom:
            raise Exception("The multinom isn't a monom : " + str(self.monoms))
        return Monom(self.monoms[0].var, self.scalar, self.monoms[0].exponent)

class Polynom(Tex):
    def __init__(self, *multinoms):
        self.scalar = Scalar(0) #Additive scalar !
        self.multinoms = []
        for multinom in multinoms:
            self += multinom

    def __neg__(self):
        return self*-1

    def __iadd__(self, other):
        if isNumber(other):
            self.scalar += Scalar(other)
            return self

        if type(other) in [Monom, Multinom]:
            if type(other) == Multinom and other.isMonom():
                self += other.getMonom()
                return self
            else:
                for multinom in self.multinoms:
                    if type(multinom) in [Monom, Multinom] and prop(multinom, other):
                        temp = multinom+other
                        if type(temp) in (Monom, Multinom) and temp.scalar == 0 or isNumber(temp) and temp == 0:
                            self.multinoms.remove(multinom)
                            if self.multinoms == []:
                                return Scalar(0)
                        else:
                            multinom += other
                        return self
                
                self.multinoms.append(other if type(other) == Multinom else Multinom(other))

        elif type(other) == Polynom:
            for multinom in other.multinoms:
                self += multinom

        return self

    def __add__(self, other):
        copy = cp.deepcopy(self)
        copy += other
        return copy

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + (-other)

    def __imul__(self, other):
        if isNumber(other) or type(other) in (Monom, Multinom):
            for multinom in self.multinoms:
                temp = multinom*other
                if type(temp) in (Monom, Multinom) and temp.scalar == 0 or isNumber(temp) and temp == 0:
                    self.multinoms.remove(multinom)
                    if self.multinoms == []:
                        return Scalar(0)
                else:
                    multinom *= other

        elif type(other) == Polynom :
            for multinom in other.multinoms:
                self *= multinom
                
        return self

    def __mul__(self, other):
        copy = cp.deepcopy(self)
        copy *= other
        return copy

    def __rmul__(self, other):
        return self*other

    def __str__(self):
        res = self.scalar.toStr() if self.scalar != 0 else ""
        for multinom in self.multinoms:
            if multinom.scalar.real < 0:
                res += " - " + str(-multinom)
            else:
                res += " + " + str(multinom)
        if res[:3] == " + ":
            res = res[3:]
        return res

def isExpr(obj):
    return type(obj) in (Monom, Multinom, Polynom)

def prop(a, b):
    if type(a) not in (Monom, Multinom) or type(b) not in (Monom, Multinom):
        print("Invalid arg for 'prop' : ", type(a), " or ", type(b), ", False returned")
        return False

    if type(a) == type(b) == Monom:
        return a.var == b.var and a.exponent == b.exponent

    if type(a) == Monom:
        return prop(b, a)
    
    if type(b) == Monom:
        return len(a.monoms) == 1 and prop(a.monoms[0], b)

    else:
        if len(a.monoms) != len(b.monoms):
            return False
        aMonomsS = sorted(a.monoms, key = lambda monom : monom.var.name)
        bMonomsS = sorted(b.monoms, key = lambda monom : monom.var.name)

        for i in range(len(aMonomsS)):
            if not prop(aMonomsS[i], bMonomsS[i]):
                return False
                
        return True