import numpy as np
from collections import defaultdict


class Tableau:
    # init works
    def __init__(self, A, unit_cost, b):  # A->matrix b->vector unit_cost->vector
        self.unit_cost = unit_cost
        self.A = A
        self.b = b
        self.basis_mapping = {} # will be initialized in reduce function
        #This won't work. Storing the identity columns as a map does not give us any idea about their arrangement
        #But we need to know the ordering to undertstand which column to remove from basis while swapping.
        self.c, self.d = A.shape  # c constraints d variables
        self.basis = {}
        self.Matrix = np.zeros((self.c + 2, self.d + 2))
        self.Matrix[1 : self.c + 1, 1 : self.d + 1] = self.A
        self.Matrix[self.c + 1, 1 : self.d + 1] = unit_cost  # unit cost -> row vector
        self.Matrix[1 : self.c + 1, self.d + 1] = np.transpose(b)
        self.Matrix[self.c + 1, self.d + 1] = 0
        print(self.Matrix)
        self.reduce()

    def pivot(self, p, q):  # pivot about (p,q)
        self.Matrix[p] /= self.Matrix[p, q]
        for row in range(1, self.c + 2):
            if row != p:
                self.Matrix[row] -= self.Matrix[p] * self.Matrix[row][q]
                print(self.Matrix)
        print("Pivot done brdr")
        self.printMat()

    def find_q(self):  # find column in non-basis to swap
        min_ind = 0
        min_val = -1
        for i in range(1, self.d + 1):
            if min_val == -1:
                min_val = self.Matrix[self.c + 1][i]
                min_ind = i
            if self.Matrix[self.c + 1][i] < min_val and (
                self.basis.get(i) == None or self.basis.get(i) == False
            ):
                min_val = self.Matrix[self.c + 1][i]
                min_ind = i
        if min_ind == 0 or min_val >= 0:
            return [False, -1]
        return [True, min_ind]

    def find_p(self, q):  # find column to swap with q
        min_ind = 0
        min_val = -1
        for i in range(1, self.c + 1):
            if self.Matrix[i][q] > 0:
                if min_val == -1:
                    min_val = self.Matrix[i][self.d + 1] / self.Matrix[i][q]
                    min_ind = i
                if self.Matrix[i][self.d + 1] / self.Matrix[i][q] < min_val:
                    min_ind = i
                    min_val = self.Matrix[i][self.d + 1] / self.Matrix[i][q]
        if min_ind == 0:
            return [False, -1]
        return [True, min_ind]

    def reduce(self):  # find an initial bfs for the problem, //artifical problem
        pass

    def printMat(self):
        V = self.Matrix[1:, 0:]
        V = V[0:, 1:]
        print(V)

    def solve(self):
        count = 0
        while 1:
            q = self.find_q()
            if q[0] == False:
                print("Optimum")
                return
            p = self.find_p(q[1])
            if p[0] == False:
                print("Unbounded")
                return
            print(p)
            print(q)
            self.pivot(p[1], q[1])
            print("self.basis_mapping:", self.basis_mapping)
            print("p:", p)
            r = self.basis_mapping.get(p[1])
            print(r)
            self.basis[r] = False
            self.basis[q[1]] = True
        self.printMat()
        pass


def transform_to_standard_lp():  # min ,all vars > 0, slack vars
    pass


A=np.array([[1,0,1,0,0],
            [0,1,0,1,0],
            [1,1,0,0,1]])
b=np.array([4,6,8])
print(b.shape)
unit_cost=np.array([-2,-5,0,0,0])
tab1=Tableau(A,unit_cost,b)
tab1.basis_mapping[1] = 3
tab1.basis_mapping[2] = 4
tab1.basis_mapping[3] = 5
tab1.basis[3] = True
tab1.basis[4] = True
tab1.basis[5] = True
tab1.solve()
# # tab1.pivot(2,2)
# # tab1.pivot(3,1)
# print(tab1.find_q())
# print(tab1.find_p(tab1.find_q()))
# np.set_printoptions(suppress=True)
# A = np.array([[1, 0, -2 / 7, 1 / 7], [0, 1, 1 / 14, -2 / 7]])
# b = np.array([18 / 7, 6 / 7])
# unit_cost = np.array([2, 3, 0, 0])
# tab2 = Tableau(A, unit_cost, b)
# tab2.basis_mapping[1] = True
# tab2.basis_mapping[2] = True
# tab2.printMat()
# print(tab2.find_q())
# print(tab2.find_p(tab2.find_q()[1]))
# if __name__=='__main__':
#     pass
