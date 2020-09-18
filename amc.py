class Fraction:
    def __init__(self, n, d):
        self.num = int(n / self.gcd(abs(n), abs(d)))
        self.denom = int(d / self.gcd(abs(n), abs(d)))
        if self.denom < 0:
            self.denom = abs(self.denom)
            self.num = -1*self.num

    def gcd(self, a,b):
        while b > 0:
            a, b = b, a % b
        return a

    def Add(self, other):
        return Fraction(self.num*other.denom + self.denom*other.num, self.denom*other.denom)

    def Sub(self, other):
        return Fraction(self.num*other.denom - self.denom*other.num, self.denom*other.denom)

    def Mul(self, other):
        return Fraction(self.num*other.num, self.denom*other.denom)

    def Div(self, other):
        return Fraction(self.num*other.denom, self.denom*other.num)


def num_of_transients(m):
    for r in range(len(m)):
        for c in range(len(m[r])):
            if m[r][c] != 0: break
        else: return r

def decompose(m):
    t = num_of_transients(m)
    Q = []
    for r in range(t):
        qRow = []
        for c in range(t):
            qRow.append(m[r][c])
        Q.append(qRow)

    R = []
    for r in range(t):
        rRow = []
        for c in range(t, len(m[r])):
            rRow.append(m[r][c])
        R.append(rRow)
    return Q, R

def identity(t):
    m = []
    for i in range(t):
        r = []
        for j in range(t):
            r.append(int(i == j))
        m.append(r)
    return m

def swap(m, i, j):
    n = []
    s = len(m)

    if i == j: return m

    for r in range(s):
        nRow = []
        tmpRow = m[r]
        if r == i:
            tmpRow = m[j]
        if r == j:
            tmpRow = m[i]
        for c in range(s):
            tmpEl = tmpRow[c]
            if c == i:
                tmpEl = tmpRow[j]
            if c == j:
                tmpEl = tmpRow[i]
            nRow.append(tmpEl)
        n.append(nRow)
    return n

def sort(m):
    size = len(m)
    zero_row = -1
    for r in range(size):
        sum = 0
        for c in range(size):
            sum += m[r][c]
        if sum == 0:
            zero_row = r
        if sum != 0 and zero_row > -1:
            n = swap(m, r, zero_row)
            return sort(n)
    return m

def normalize(m):
    n = []
    for r in range(len(m)):
        sum = 0
        cols = len(m[r])
        for c in range(cols):
            sum += m[r][c]

        nRow = []

        if sum == 0:
            nRow = m[r]
        else:
            for c in range(cols):
                nRow.append(Fraction(m[r][c], sum))
        n.append(nRow)
    return n

def subtract(i, q):
    s = []
    for r in range(len(i)):
        sRow = []
        for c in range(len(i[r])):
            sRow.append(Fraction.Sub(Fraction(i[r][c], 1), q[r][c]))
        s.append(sRow)
    return s

def multiply(a, b):
    m = []
    rows = len(a)
    cols = len(b[0])
    iters = len(a[0])

    for r in range(rows):
        row = []
        for c in range(cols):
            sum = Fraction(0,1)
            for i in range(iters):
                sum = Fraction.Add(sum,Fraction.Mul(a[r][i],b[i][c]))
            row.append(sum)
        m.append(row)
    return m

def transposeMatrix(m):
    t = []
    for r in range(len(m)):
        row = []
        for c in range(len(m[r])):
            if c == r:
                row.append(m[r][c])
            else:
                row.append(m[c][r])
        t.append(row)
    return t

def getMatrixMinor(m,i,j):
    return [row[:j] + row[j+1:] for row in (m[:i]+m[i+1:])]

def getMatrixDeternminant(m):
    if len(m) == 2:
        return Fraction.Sub(Fraction.Mul(m[0][0],m[1][1]),Fraction.Mul(m[0][1],m[1][0]))

    d = Fraction(0,1)
    for c in range(len(m)):
        d = Fraction.Add(d,Fraction.Mul(Fraction(m[0][c].num,m[0][c].denom*((-1)**c)),getMatrixDeternminant(getMatrixMinor(m,0,c))))
    return d

def getMatrixInverse(m):
    d = getMatrixDeternminant(m)

    if len(m) == 2:
        return [[Fraction.Div(m[1][1],d), Fraction.Div(Fraction.Mul(Fraction(-1,1),m[0][1]),d)], 
                [Fraction.Div(Fraction.Mul(Fraction(-1,1),m[1][0]),d), Fraction.Div(m[0][0],d)]]

    cofactors = []
    for r in range(len(m)):
        cofactorRow = []
        for c in range(len(m)):
            minor = getMatrixMinor(m,r,c)
            cofactorRow.append(Fraction.Mul(Fraction((-1)**(r+c), 1), getMatrixDeternminant(minor)))
        cofactors.append(cofactorRow)
    cofactors = transposeMatrix(cofactors)
    for r in range(len(cofactors)):
        for c in range(len(cofactors)):
            cofactors[r][c] = Fraction.Div(cofactors[r][c],d)
    return cofactors

def convert_to_lcd(probs):
    ret = []

    lcm = reduce(lambda x, y: least_common_multiple(x, y), [f.denom for f in probs])
    for f in probs:
        if f.denom != lcm:
            ret.append(lcm / f.denom * f.num)
        else:
            ret.append(f.num)

    ret.append(lcm)
    return ret

def least_common_multiple(a, b):
    greater = a if a > b else b

    while True:
        if greater % a == 0 and greater % b == 0:
            lcm = greater
            break
        greater += 1

    return lcm

def solution(m):
    if len(m) == 1: return [1,1]
    if len(m) == 2:
        return [0, 1] if m[0][1] == 0 else [1,1]
    m = sort(m)
    n = normalize(m)
    (q, r) = decompose(n)
    i = identity(len(q))
    s = subtract(i, q)
    v = getMatrixInverse(s)
    b = multiply(v, r)
    return convert_to_lcd(b[0])


m = [[0, 1, 0, 0, 0, 1], [4, 0, 0, 3, 2, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]
print(solution(m))

m = [[0,1], [1,2]]
print(solution(m))

m = [[2,0], [0,0]]
print(solution(m))

m = [[1,2,3], [0,0,0], [3,2,1]]
print(solution(m))

m = [[1,3,2],
        [3,1,2],
        [0,0,0]]
print(solution(m))

m = [[1,1,1,1],
        [1,4,2,3],
        [0,0,0,0],
        [0,0,0,0]]
print(solution(m))

m = [[0,1,0,0,0,1],
        [4,0,0,3,2,0],
        [0,0,0,0,0,0],
        [0,0,0,0,0,0],
        [0,0,0,0,0,0],
        [0,0,0,0,0,0]]
print(solution(m))

m = [[1,1,1,1],
        [0,0,0,0],
        [0,0,0,0],
        [1,2,3,4]]
print(solution(m))

m = [[0, 2, 1, 0, 0], [0, 0, 0, 3, 4], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]
print(solution(m))

m = [[0, 1, 0, 0, 0, 1], [4, 0, 0, 3, 2, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]
print(solution(m))

m = [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
print(solution(m))

m = [[1, 1, 1, 0, 1, 0, 1, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 1, 1, 1, 0, 1, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 1, 0, 1, 1, 1, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 1, 0, 1, 0, 1, 1, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 1, 0, 1, 0, 1, 0, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
print(solution(m))
