import numpy as np
from collections import defaultdict
import csv 
  
# Open file  
with open('inp.csv') as file_obj: 
      
    # Create reader object by passing the file  
    # object to reader method 
    reader_obj = csv.reader(file_obj) 
    reader_obj=list(reader_obj)
      
    # Iterate over each row in the csv  
    # file using reader object 
    # for idx, row in enumerate(reader_obj): 
    #     print(row)
        



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
        # self.Matrix[-1][-1] = -34

    def pivot(self, p, q):  # pivot about (p,q)
        EPSILON = 1e-8  # Tolerance for small values
        self.Matrix[p] /= self.Matrix[p, q]
        for row in range(1, self.c + 2):
            if row != p:
                self.Matrix[row] -= self.Matrix[p] * self.Matrix[row][q]
                self.Matrix[row][abs(self.Matrix[row]) < EPSILON] = 0
                # print(self.Matrix)
        print("Pivot done brdr")
        self.printMat()

    def find_q(self):  # find column in non-basis to swap
        min_ind = 0
        min_val = None
        for i in range(1, self.d + 1):
            # print('considering column ', i)
            # print(self.Matrix[self.c+1][i])
            if min_val is None:
                min_val = self.Matrix[self.c + 1][i]
                min_ind = i
            if self.Matrix[self.c + 1][i] < min_val and (
                self.is_basis.get(i) == None or self.is_basis.get(i) == False
            ):
                min_val = self.Matrix[self.c + 1][i]
                min_ind = i
        if min_ind == 0 or min_val >= 0:
            # print('returning false, q not found:')
            # print(f'min_val={min_val} min_ind ={min_ind}')
            return [False, -1]
        return [True, min_ind]

    def find_p(self, q):  # find column to swap with q
        min_ind = 0
        min_val = None
        # print('while finding p self.d is ', self.d)
        for i in range(1, self.c + 1):
            if self.Matrix[i][q] > 0:
                if min_val==None or self.Matrix[i][self.d + 1] / self.Matrix[i][q] < min_val: ########+2
                    min_ind = i
                    min_val = self.Matrix[i][self.d + 1] / self.Matrix[i][q]
                    # print(f'updated min_val {min_val} min_ind {min_ind}')
                    # print(f'divided {self.Matrix[i][self.d+1]} with {self.Matrix[i][q]}')
        # print('min val',min_val)
        # print('min_ind',min_ind)
        if min_ind == 0:
            return [False, -1]
        return [True, min_ind]

    def find_initial_bfs(self):  # find an initial bfs for the problem, //artifical problem
        artificial_A=self.A
        Ic=np.identity(self.A.shape[0])
        # print('art shape: ',  artificial_A.shape)
        # print('bshape', self.b.shape)
        btmp = self.b
        # print('b.shape is', self.b.shape)
        btmp = np.reshape(btmp, (self.b.shape[0],1))
        
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
        # print('self.c', self.c)
        # print('self.d', self.d)
        self.printMat()

        self.solve()
        print('Initial bfs found brdr')
        self.printMat()
        for key,val in self.basis_ordering.items():
            if (val>self.d):
                print("Infeasible")
                exit(0)
        print(self.basis_ordering)
        self.ready_to_solve()
        

    def ready_to_solve(self): #retrieve the original matrix from the artificial one
        print('artificial problem solved--------------------------###################')
        # print("ye lo sizes")
        # print(self.Matrix[0:,:self.original_d+2].shape)
        # print(self.Matrix[:,self.d+1:self.d+2].shape) 
        tmp=np.hstack((self.Matrix[0:,:self.original_d+1],self.Matrix[:,self.d+1:self.d+2]))
        # self.printMat(tmp)
        self.Matrix=tmp
        self.Matrix[-1,1:-1]=self.unit_cost
        self.Matrix[-1,-1]=0
        # self.printMat(self.Matrix)
        for key,val in self.basis_ordering.items():
           print('key, val is -------------', key, val)
           last=self.Matrix[-1,val]
           self.Matrix[-1]-=self.Matrix[key]*last
        print('bhai party dede')
        self.c=self.original_c
        self.d=self.original_d
        # self.printMat()


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
                self.printMat()
                return
            p = self.find_p(q[1])
            if p[0] == False:
                print(p, q)
                self.printMat()
                print("Unbounded")
                exit(0)
            print('p is',p)
            print('q is',q)
            self.pivot(p[1], q[1])
            print("self.basis_ordering:", self.basis_ordering)
            leaving_column = self.basis_ordering.get(p[1])
            print('leaving column: ',leaving_column)
            self.is_basis[leaving_column] = False
            self.is_basis[q[1]] = True
            self.basis_ordering[p[1]]=q[1]
            # self.Matrix = np.round(self.Matrix, 1)
    
    def displayans(self, slack_start):
        final_ans = []
        answer = [0]*(self.d+1)
        for key, value in self.basis_ordering.items():
            answer[value] = self.Matrix[key][-1]
        print('The final answer is:', answer)
        print(f'self.basis rows {self.basis_ordering}')
        print(slack_start)
        for i in range(1, slack_start, 2):
            print(f'variable {(i+1)//2} is ', answer[i]-answer[i+1])
            final_ans.append(answer[i]-answer[i+1])
        return final_ans


def transform_to_standard_lp(type, d, lessthan_list, greaterthan_list, eq_list, unit_cost):  # min ,all vars > 0, slack vars
    #todo
    print('at the start', eq_list)
    #change the objective if needeed based on the type
    eq_list  = eq_list[:]
    lessthan_list = lessthan_list[:]
    greaterthan_list = greaterthan_list[:]
    if type=="MAX":
        unit_cost=[-x for x in unit_cost]

    # remove stuff from the greater than list and put it in the less than list
    for eqn in greaterthan_list:
        eqn=[-elem for elem in eqn]
        lessthan_list.append(eqn)
    # remove all redundant inequations using the O(n^2d) trick 

    
    # replace each variable by u - v
    # change the rows of A
    for i in range(0, len(lessthan_list)):
        list  = lessthan_list[i]
        last_el = list[-1]
        list = list[:-1]
        neglist = [-1*x for x in list]
        interleaved = [x for pair in zip(list, neglist) for x in pair]
        lessthan_list[i] = interleaved
        lessthan_list[i].append(last_el) 
    #u-v for lesshthanlist added

    #change unit costs
    negunit_cost = [-1*x for x in unit_cost]
    interleaved = [x for pair in zip(unit_cost, negunit_cost) for x in pair]
    unit_cost = interleaved
    #added u-v for unit cost
    print(f'here by magic {eq_list}')
    #  now add the slack variables
    sz=len(lessthan_list)
    slack_start = 2*d+1
    for index,value in enumerate(lessthan_list): #adding slack variables to less than lists
        to_add=[0]*len(lessthan_list)
        to_add[index]=1
        lessthan_list[index]=lessthan_list[index][:-1]+to_add+lessthan_list[index][-1:]
        print(f'here by human clay magic {eq_list}')

    
    print(f'b4 u-v eq_list is', eq_list)
    for i in range(0, len(eq_list)):
        list  = eq_list[i]
        last_el = list[-1]
        list = list[:-1]
        neglist = [-1*x for x in list]
        interleaved = [x for pair in zip(list, neglist) for x in pair]
        eq_list[i] = interleaved
        eq_list[i].append(last_el) 
    #u-v fo requal list
    print(f'after u-v eq_list is', eq_list)
    
    print(len(lessthan_list[0]) - len(eq_list[0]), '0s added manually')
    eq_list += [0] * (len(lessthan_list[0]) - len(eq_list[0]))
    print(f'kombucha {eq_list}')

    
    eq_list=lessthan_list #all equations 
    print(sz)
    print('less than list[0] is', lessthan_list[0])
    print('lenght of ith row of less than list is:', len(lessthan_list[0]))
    print(unit_cost)
    unit_cost+=[0]*sz
    print(f'eq_list {eq_list}')
    A=np.array([row[:-1] for row in eq_list])
    b=np.array([row[-1] for row in eq_list])
    unit_cost=np.array(unit_cost)
    return A,b,unit_cost,slack_start


    
#function to take input for a single test case
#will be called multiple times inside while loop
def take_input_of_one_test_case(test_case,equal_list,greater_than_list,less_than_list):
    #read the first two lists
    first_list=test_case[0]
    second_list=test_case[1]
    type=first_list[0]
    number_of_actual_variables=int(first_list[1])

    #auxillary function to parse a single row
    def parse_input_row(row):
        idx=0
        print(f'row {row}')
        v=[0]*number_of_actual_variables
        while idx<len(row):
            # print(row[idx])
            splitted=row[idx].split('@') #splitted is a list
            if len(splitted)>1:
                x,j=int(splitted[0]),int(splitted[1])
                if j<=idx:
                    print(f'j {j} idx {idx}')
                    print("Error while parsing")
                    exit(0)
                else: 
                    idx=j
                v[idx-1]=x
            else:
                idx+=1
                v[idx-1]=int(splitted[0])
        return v #returns a list
    
    unit_cost=parse_input_row(second_list)
    curr_row=2
    while curr_row<len(test_case):
        print(f'test_case[curr_row][1] {test_case[curr_row][1]}')
        last_element=int(test_case[curr_row][0])
        parsed_row=parse_input_row(test_case[curr_row][2:])
        parsed_row.append(last_element)
        print(f'parsedd row {parsed_row}') #thik hai bhai
        if test_case[curr_row][1]=='=':
           print(f'parsed row {parsed_row} went to equal list')
           equal_list.append(parsed_row) 
        elif test_case[curr_row][1]=='>=':
            print(f'parsed row {parsed_row} went to great list')
            greater_than_list.append(parsed_row)
        else:
            print(f'parsed row {parsed_row} went to less list')
            less_than_list.append(parsed_row)
        curr_row+=1
    
    print(f'type {type}')
    print(f'd {number_of_actual_variables}')
    print(f'lessthan_list {less_than_list}')
    print(f'greaterthan_list {greater_than_list}')
    print(f'eq_list {equal_list}')
    print(f'unit_cost {unit_cost}')
    return type,number_of_actual_variables, less_than_list, greater_than_list, equal_list, unit_cost

def take_input(reader_obj):
    current_index=0
    while current_index<len(reader_obj):
        print(f'current index {current_index}')
        equal_list=[]
        greater_than_list=[]
        less_than_list=[]
        start=current_index
        current_index+=1
        while current_index<len(reader_obj) and 'M' not in reader_obj[current_index][0]:
            current_index+=1
        end=current_index
        print(reader_obj)
        type,number_of_actual_variables, less_than_list, greater_than_list, equal_list, unit_cost=take_input_of_one_test_case(reader_obj[start:end],equal_list,greater_than_list,less_than_list)
        print("we will do this correctly")
        A,b,unit_cost,slack_start=transform_to_standard_lp(type,number_of_actual_variables,less_than_list,greater_than_list,equal_list,unit_cost)
        print("md salah")
        tableau=Tableau(A=A,b=b,unit_cost=unit_cost)
        tableau.solve()
        tableau.printMat()
        tableau.displayans(slack_start)
        print("Input taken successfully")
        print(f'current index sar {current_index}')

    


print(reader_obj)
take_input(reader_obj=reader_obj)

    




#test1 yay this works
# A=np.array([[1,0,1,0,0],
#             [0,1,0,1,0],
#             [1,1,0,0,1]])
# b=np.array([4,6,8])
# # print(b.shape)
# unit_cost=np.array([-2,-5,0,0,0])
# tab1=Tableau(A=A,unit_cost=unit_cost,b=b)
# # tab1.basis_ordering[1] = 3
# # tab1.basis_ordering[2] = 4
# # tab1.basis_ordering[3] = 5
# # tab1.basis[3] = True
# # tab1.basis[4] = True
# # tab1.basis[5] = True
# tab1.solve()

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

# A=np.array([[1,0,1,0],
#             [0,1,0,1]])
# unit_cost=np.array([1,-1,0,0])
# b=np.array([1,1])
# newtab=Tableau(unit_cost=unit_cost,A=A,b=b)
# newtab.solve()
# newtab.printMat()


#test from https://www.uobabylon.edu.iq/eprints/publication_3_29932_132.pdf
# A = np.array([
#                 [1, -1, 1, -1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0]
#               , [1, -1, 1, -1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0]
#               , [4, -4, 2, -2, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0]
#               , [2, -2, 4, -4, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0]
#               , [1, -1, 1, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
#               ]
#               )
# b = np.array([6, 7, 8, 8, 5])
# unit_cost = np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1])

# tabket = Tableau(A=A, b=b, unit_cost=unit_cost)
# print('after init')
# tabket.printMat()
# tabket.solve()
# print('after soling')
# tabket.printMat()
# print('adfsa adf4r3 dwkfjdsfj')
# print(tabket.basis_ordering)

# A = np.array([
#                 [4, 2, 1, 0]
#               , [2, 4, 0, 1]
#               ]
#             )
# b = np.array([8, 8])
# unit_cost = np.array([-6, -5, 0, 0])
# tabket = Tableau(A=A, b=b, unit_cost=unit_cost)
# print('after init')
# tabket.printMat()
# tabket.solve()
# print('after soling')
# tabket.printMat()
# A,b,unit_cost,slack_start=transform_to_standard_lp('MIN',2, [[1,1,6],[1,1,7],[4,2,8],[2,4,8]],[[-1,-1,-5]],[],[-6,-5])
# # A, b, unit_cost, slack_start = transform_to_standard_lp('MAX', 3, [[2, 1, 3, 10], [3, 4, 2, 12], [1, 2, 4, 8]], [[-1, -2, -3, -5]], [], [-8, -6, -7])
# print(A)
# print(b)
# print(unit_cost)
# cool_tab=Tableau(A=A,unit_cost=unit_cost,b=b)
# cool_tab.solve()
# cool_tab.printMat()
# print(cool_tab.basis_ordering)

# print('###############################################$$')
# cool_tab.displayans(slack_start=slack_start)

# print('###############################################$$')