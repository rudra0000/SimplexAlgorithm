import numpy as np
from collections import defaultdict
class Tableau:
    #init works
    def __init__(self,A,unit_cost,b): #A->matrix b->vector unit_cost->vector 
        self.unit_cost=unit_cost
        self.A=A
        self.b=b
        self.basis_mapping=defaultdict() #will be initialized in reduce function
        self.c,self.d=A.shape #c constraints d variables
        self.Matrix=np.zeros((self.c+2,self.d+2))
        self.Matrix[1:self.c+1,1:self.d+1]=self.A
        self.Matrix[self.c+1,1:self.d+1]=unit_cost #unit cost -> row vector
        self.Matrix[1:self.c+1,self.d+1]=np.transpose(b)
        self.Matrix[self.c+1,self.d+1]=0
        print(self.Matrix)
        self.reduce()

    def pivot(self,p,q): # pivot about (p,q)
        self.Matrix[p]/=self.Matrix[p,q]
        for row in range(1,self.c+2):
            if row!=p:
                self.Matrix[row]-=self.Matrix[p]*self.Matrix[row][q]
                print(self.Matrix)
        print("Pivot done brdr")
        print(self.Matrix)

    def find_q(self): #find column in non-basis to swap
        pass
    def find_p(self): #find column to swap with q
        pass
    def reduce(self):# find an initial bfs for the problem, //artifical problem
        pass

def transform_to_standard_lp():# min ,all vars > 0, slack vars
    pass

A=np.array([[1,0,1,0,0],
            [0,1,0,1,0],
            [1,1,0,0,1]])
b=np.array([4,6,8])
print(b.shape)
unit_cost=np.array([-2,-5,0,0,0])
tab1=Tableau(A,unit_cost,b)
tab1.pivot(2,2)
# if __name__=='__main__':
#     pass