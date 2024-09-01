import numpy as np
from collections import defaultdict


class Tableau:
    # init works
    def __init__(self,A,b,unit_cost):  # A->matrix b->vector unit_cost->vector
        self.unit_cost = unit_cost
        self.A = A
        self.b = np.transpose(b)
        self.basis_ordering = {} # will be initialized in reduce function
        #This won't work. Storing the identity columns as a map does not give us any idea about their arrangement
        #But we need to know the ordering to undertstand which column to remove from basis while swapping.
        self.c, self.d = A.shape  # c constraints d variables
        self.is_basis = {}
        self.Matrix = np.zeros((self.c + 2, self.d + 2))
        self.Matrix[1 : self.c + 1, 1 : self.d + 1] = self.A
        self.Matrix[self.c + 1, 1 : self.d + 1] = unit_cost  # unit cost -> row vector
        self.Matrix[1 : self.c + 1, self.d + 1] = b
        self.Matrix[self.c + 1, self.d + 1] = 0
        # print(self.Matrix)
        self.find_initial_bfs()

    def pivot(self, p, q):  # pivot about (p,q)
        self.Matrix[p] /= self.Matrix[p, q]
        for row in range(1, self.c + 2):
            if row != p:
                self.Matrix[row] -= self.Matrix[p] * self.Matrix[row][q]
                # print(self.Matrix)
        print("Pivot done brdr")
        self.printMat()

    def find_q(self):  # find column in non-basis to swap
        min_ind = 0
        min_val = None
        for i in range(1, self.d + 1):
            print('considering column ', i)
            print(self.Matrix[self.c+1][i])
            if min_val is None:
                min_val = self.Matrix[self.c + 1][i]
                min_ind = i
            if self.Matrix[self.c + 1][i] < min_val and (
                self.is_basis.get(i) == None or self.is_basis.get(i) == False
            ):
                min_val = self.Matrix[self.c + 1][i]
                min_ind = i
        if min_ind == 0 or min_val >= 0:
            print('returning false, q not found:')
            print(f'min_val={min_val} min_ind ={min_ind}')
            return [False, -1]
        return [True, min_ind]

    def find_p(self, q):  # find column to swap with q
        min_ind = 0
        min_val = None
        print('while finding p self.d is ', self.d)
        for i in range(1, self.c + 1):
            if self.Matrix[i][q] > 0:
                if min_val==None or self.Matrix[i][self.d + 1] / self.Matrix[i][q] < min_val: ########+2
                    min_ind = i
                    min_val = self.Matrix[i][self.d + 1] / self.Matrix[i][q]
                    print(f'updated min_val {min_val} min_ind {min_ind}')
                    print(f'divided {self.Matrix[i][self.d+1]} with {self.Matrix[i][q]}')
        print('min val',min_val)
        print('min_ind',min_ind)
        if min_ind == 0:
            return [False, -1]
        return [True, min_ind]

    def find_initial_bfs(self):  # find an initial bfs for the problem, //artifical problem
        artificial_A=self.A
        Ic=np.identity(A.shape[0])
        print('art shape: ',  artificial_A.shape)
        print('bshape', b.shape)
        btmp = b
        print('b.shape is', b.shape)
        btmp = np.reshape(btmp, (b.shape[0],1))
        
        artificial_A=np.hstack((artificial_A,Ic,btmp)) #make sure 1 based indexing is followed
        lower_row=np.hstack((np.zeros(self.d),np.ones(self.c),np.zeros(1)))
        artificial_A=np.vstack((artificial_A,lower_row))
        self.artificialMatrix=np.zeros((self.c+2,self.d+2+self.c))
        self.artificialMatrix[1:,1:]=artificial_A #deal with the last row and define basis
        
        for col_no in range(self.d+1,self.d+1+self.c):
            self.is_basis[col_no]=True
            self.basis_ordering[col_no-self.d]=col_no
            self.artificialMatrix[self.c+1]-=self.artificialMatrix[col_no-self.d] #dealing with the last row making reduced cost corresponding to the basis vectors = 0
        self.original_c=self.c
        self.original_d=self.d
        self.d=self.c+self.d
        self.c=self.c
        self.original_matrix=self.artificialMatrix #wastage of space
        self.Matrix=self.artificialMatrix
        print('self.c', self.c)
        print('self.d', self.d)
        self.printMat()

        self.solve()
        print('Initial bfs found brdr')
        self.printMat()
        print(self.basis_ordering)
        self.ready_to_solve()
        

    def ready_to_solve(self): #retrieve the original matrix from the artificial one
        print("ye lo sizes")
        print(self.Matrix[0:,:self.original_d+2].shape)
        print(self.Matrix[:,self.d+1:self.d+2].shape) 
        tmp=np.hstack((self.Matrix[0:,:self.original_d+1],self.Matrix[:,self.d+1:self.d+2]))
        self.printMat(tmp)
        self.Matrix=tmp
        self.Matrix[-1,1:-1]=unit_cost
        self.Matrix[-1,-1]=0
        self.printMat(self.Matrix)
        for key,val in self.basis_ordering.items():
           print('key, val is -------------', key, val)
           last=self.Matrix[-1,val]
           self.Matrix[-1]-=self.Matrix[key]*last
        print('bhai party dede')
        self.c=self.original_c
        self.d=self.original_d
        self.printMat()


    def printMat(self, matr=None):
        if matr is None:
            V = self.Matrix[1:, 0:]
            V = V[0:, 1:]
            print(V)
        else:
            V = matr[1:, 0:]
            V = V[0:, 1:]
            print(V)


    def solve(self):
        while True:
            q = self.find_q()
            if q[0] == False:
                print("Optimum")
                return
            p = self.find_p(q[1])
            if p[0] == False:
                print("Unbounded")
                return
            print('p is',p)
            print('q is',q)
            self.pivot(p[1], q[1])
            print("self.basis_ordering:", self.basis_ordering)
            leaving_column = self.basis_ordering.get(p[1])
            print('leaving column: ',leaving_column)
            self.is_basis[leaving_column] = False
            self.is_basis[q[1]] = True
            self.basis_ordering[p[1]]=q[1]


def transform_to_standard_lp():  # min ,all vars > 0, slack vars
    #todo
    pass

#test1 yay this works
A=np.array([[1,0,1,0,0],
            [0,1,0,1,0],
            [1,1,0,0,1]])
b=np.array([4,6,8])
# print(b.shape)
unit_cost=np.array([-2,-5,0,0,0])
tab1=Tableau(A=A,unit_cost=unit_cost,b=b)
# tab1.basis_ordering[1] = 3
# tab1.basis_ordering[2] = 4
# tab1.basis_ordering[3] = 5
# tab1.basis[3] = True
# tab1.basis[4] = True
# tab1.basis[5] = True
tab1.solve()

#test2
# A=np.array([[2,1,1,0],
#             [1,4,0,1]])
# b=np.array([3,4])
# unit_cost=np.array([-7,-6,0,0])
# tab2=Tableau(A=A,b=b,unit_cost=unit_cost)
# # tab2.basis_ordering[1]=3
# # tab2.basis_ordering[2]=4 #1 based indexing
# # tab2.basis={1:False,2:False,3:True,4:True}
# tab2.solve()

#test3
# A=np.array([[-2/3,0,1,1/3],
#             [5/3,1,0,-1/3]])
# b=np.array([4/3,8/3])
# unit_cost=np.array([-3,-5,0,0])
# tab3=Tableau(A=A,b=b,unit_cost=unit_cost)
# tab3.basis_ordering={1:3,2:2}
# tab3.basis={1:False,2:True,3:True,4:False}
# tab3.printMat()
# tab3.solve()
# tab3.printMat()
# print(tab3.basis_ordering)

#test4
# A=np.array([[1,1,1,0],
#             [2,1,0,1]])
# b=np.array([12,16])

# unit_cost=np.array([-40,-30,0,0])
# tab4=Tableau(A=A,b=b,unit_cost=unit_cost)

# tab4.solve()

# # tab1.pivot(2,2)
# # tab1.pivot(3,1)
# print(tab1.find_q())
# print(tab1.find_p(tab1.find_q()))
# np.set_printoptions(suppress=True)
# A = np.array([[1, 0, -2 / 7, 1 / 7], [0, 1, 1 / 14, -2 / 7]])
# b = np.array([18 / 7, 6 / 7])
# unit_cost = np.array([2, 3, 0, 0])
# tab2 = Tableau(A, unit_cost, b)
# tab2.basis_ordering[1] = True
# tab2.basis_ordering[2] = True
# tab2.printMat()
# print(tab2.find_q())
# print(tab2.find_p(tab2.find_q()[1]))
# if __name__=='__main__':
#     pass


#test5
# A=np.array([[4,2,-1,0],
#             [1,4,0,-1]])
# b=np.array([12,6])
# unit_cost=np.array([2,3,0,0])
# tab5=Tableau(A=A,unit_cost=unit_cost,b=b)
# tab5.solve()



#doesn't work
#test3
# A=np.array([[-2/3,0,1,1/3],
#             [5/3,1,0,-1/3]])
# b=np.array([4/3,8/3])
# unit_cost=np.array([-3,-5,0,0])
# tab3=Tableau(A=A,b=b,unit_cost=unit_cost)
# # tab3.basis_ordering={1:3,2:2}
# # tab3.basis={1:False,2:True,3:True,4:False}
# # tab3.printMat()
# tab3.solve()


# A = np.array([[1, 1, 1, 0], [5, 3, 0, -1]])
# b = np.array([4, 8])
# unit_cost = [-3, -5, 0, 0]
# tabnu = Tableau(A=A, b=b, unit_cost=unit_cost)
# tabnu.solve()
# print(tabnu.basis_ordering)

#works 
# A = np.array([[4, 2, -1, 0], [1, 4, 0, -1]])
# b = np.array([12, 6])
# unit_cost = np.array([2,3, 0, 0])
# thetab = Tableau(A=A, b=b, unit_cost=unit_cost)
# thetab.solve()


#AtoZ math.com
# A=np.array([[3,5,2,1,0,0],
#             [4,4,4,0,1,0],
#             [2,4,5,0,0,1]])
# b=np.array([60,72,100])
# unit_cost=np.array([-5,-10,-8,0,0,0])
# newtab=Tableau(A=A,unit_cost=unit_cost,b=b)
# newtab.solve()
# newtab.printMat()

A=np.array([[1,0,1,0],
            [0,1,0,1]])
unit_cost=np.array([1,-1,0,0])
b=np.array([1,1])
newtab=Tableau(unit_cost=unit_cost,A=A,b=b)
newtab.solve()
newtab.printMat()