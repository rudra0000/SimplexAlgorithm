import numpy as np
from collections import defaultdict
import csv,math
  
# Open file  
with open('inp.csv') as file_obj: 
      
    # Create reader object by passing the file  
    # object to reader method 
    reader_obj = csv.reader(file_obj) 
    reader_obj=list(reader_obj)
        



class Tableau:
    # init works
    def __init__(self,A,b,unit_cost):  # A->matrix b->vector unit_cost->vector
        self.unit_cost = unit_cost
        self.A = A
        self.b = np.transpose(b)
        self.is_valid=True
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
        if not self.find_initial_bfs():
            self.is_valid=False
        # self.Matrix[-1][-1] = -34

    def pivot(self, p, q):  # pivot about (p,q)
        EPSILON = 1e-8  # Tolerance for small values
        self.Matrix[p] /= self.Matrix[p, q]
        for row in range(1, self.c + 2):
            if row != p:
                self.Matrix[row] -= self.Matrix[p] * self.Matrix[row][q]
                self.Matrix[row][abs(self.Matrix[row]) < EPSILON] = 0

    def find_q(self):  # find column in non-basis to swap
        min_ind = 0
        min_val = None
        for i in range(1, self.d + 1):
            if min_val is None:
                min_val = self.Matrix[self.c + 1][i]
                min_ind = i
            if self.Matrix[self.c + 1][i] < min_val and (
                self.is_basis.get(i) == None or self.is_basis.get(i) == False
            ):
                min_val = self.Matrix[self.c + 1][i]
                min_ind = i
        if min_ind == 0 or min_val >= 0:
            return [False, -1]
        return [True, min_ind]

    def find_p(self, q):  # find column to swap with q
        min_ind = 0
        min_val = None
        for i in range(1, self.c + 1):
            if self.Matrix[i][q] > 0:
                if min_val==None or self.Matrix[i][self.d + 1] / self.Matrix[i][q] < min_val: ########+2
                    min_ind = i
                    min_val = self.Matrix[i][self.d + 1] / self.Matrix[i][q]
        if min_ind == 0:
            return [False, -1]
        return [True, min_ind]

    def find_initial_bfs(self):  # find an initial bfs for the problem, //artifical problem
        artificial_A=self.A
        Ic=np.identity(self.A.shape[0])
        btmp = self.b
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
        self.printMat()

        self.solve()
        for key,val in self.basis_ordering.items():
            if (val>self.d):
                print("Infeasible")
                return False
        self.ready_to_solve()
        

    def ready_to_solve(self): #retrieve the original matrix from the artificial one
        print('artificial problem solved--------------------------###################')
        tmp=np.hstack((self.Matrix[0:,:self.original_d+1],self.Matrix[:,self.d+1:self.d+2]))
        self.Matrix=tmp
        self.Matrix[-1,1:-1]=self.unit_cost
        self.Matrix[-1,-1]=0
        for key,val in self.basis_ordering.items():
           last=self.Matrix[-1,val]
           self.Matrix[-1]-=self.Matrix[key]*last
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
        if not self.is_valid:
            return None
        while True:
            q = self.find_q()
            if q[0] == False:
                print("Optimum")
                self.printMat()
                return
            p = self.find_p(q[1])
            if p[0] == False:
                self.printMat()
                print("Unbounded")
                return None
            self.pivot(p[1], q[1])
            leaving_column = self.basis_ordering.get(p[1])
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

def remove_redundant_equalities_rref(eq_list):
    if len(eq_list)==0 :
        return []
    cpy=eq_list[:]
    def pivot(row,col):
        eq_list[row]=[x/eq_list[row][col] for x in eq_list[row]]
        for r in range(0,len(eq_list)):
            if r!=row:
                eq_list[r]=[x-y*eq_list[r][col] for x,y in zip(eq_list[r],eq_list[row])]
    
    is_available={i:True for i in range(0,len(eq_list))}
    for col in range(len(eq_list[0])):
        #loop to find the first non zero row
        for row in range(len(eq_list)):
            if is_available[row] and eq_list[row][col]!=0:
                pivot(row,col)
                is_available[row]=False
    filtered_list = [lst2 for lst, lst2 in zip(eq_list,cpy) if any(lst)]
    return filtered_list



def transform_to_standard_lp(type_of_lp, d, lessthan_list, greaterthan_list, eq_list, unit_cost):  # min ,all vars > 0, slack vars
    #todo
    #removing redundancies in the equals list
    # equality_set=set() #not doing as of now 
    # for eqn in eq_list:
    #     gcd=0
    #     for elem in eqn:
    #         gcd=math.gcd(elem,gcd)
    #     for elem in eqn:
    #         elem/=gcd
    #     print(f'eqn is {eqn}')
    #     equality_set.add(tuple(eqn))
        
    # eq_list=[[x for x in tup] for tup in equality_set]
    eq_list=remove_redundant_equalities_rref(eq_list)
    if len(eq_list)>d:
        print("Overconstrained ,infeasible")
        return False,None,None,None
    #change the objective if needeed based on the type_of_lp
    eq_list  = eq_list[:]
    lessthan_list = lessthan_list[:]
    greaterthan_list = greaterthan_list[:]
    if type_of_lp=="MAX":
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
    #  now add the slack variables
    sz=len(lessthan_list)
    slack_start = 2*d+1
    for index,value in enumerate(lessthan_list): #adding slack variables to less than lists
        to_add=[0]*len(lessthan_list)
        to_add[index]=1
        lessthan_list[index]=lessthan_list[index][:-1]+to_add+lessthan_list[index][-1:]
    
    for i in range(0, len(eq_list)):
        list1  = eq_list[i]
        last_el = list1[-1]
        list1 = list1[:-1]
        neglist = [-1*x for x in list1]
        interleaved = [x for pair in zip(list1, neglist) for x in pair]
        eq_list[i] = interleaved
        eq_list[i].append(last_el) 
    #u-v fo requal list
    
    # print(len(lessthan_list[0]) - len(eq_list[0]), '0s added manually')
    
    # if(len(lessthan_list) != 0):

    # print('that loser idff is', len(lessthan_list[0]) - len(eq_list[0]))

    # print(len(lessthan_list[0]))
    for i in range(len(eq_list)):
        if len(lessthan_list) == 0:
            break
        tmp_list = eq_list[i][:-1] + [0] * (len(lessthan_list[0]) - len(eq_list[i])) + [eq_list[i][-1]]
        eq_list[i] = tmp_list
    


    for list1 in lessthan_list:
        eq_list.append(list1)
    unit_cost+=[0]*sz
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
    type_of_lp=first_list[0]
    number_of_actual_variables=int(first_list[1])

    #auxillary function to parse a single row
    def parse_input_row(row):
        idx=0
        v=[0]*number_of_actual_variables
        while idx<len(row):
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
        last_element=int(test_case[curr_row][0])
        parsed_row=parse_input_row(test_case[curr_row][2:])
        parsed_row.append(last_element)
        if test_case[curr_row][1]=='=':
           equal_list.append(parsed_row) 
        elif test_case[curr_row][1]=='>=':
            greater_than_list.append(parsed_row)
        else:
            less_than_list.append(parsed_row)
        curr_row+=1
    
    # print(f'type_of_lp {type_of_lp}')
    # print(f'd {number_of_actual_variables}')
    # print(f'lessthan_list {less_than_list}')
    # print(f'greaterthan_list {greater_than_list}')
    # print(f'eq_list {equal_list}')
    # print(f'unit_cost {unit_cost}')
    return type_of_lp,number_of_actual_variables, less_than_list, greater_than_list, equal_list, unit_cost

def take_input(reader_obj):
    current_index=0
    test_case=1
    while current_index<len(reader_obj):
        equal_list=[]
        greater_than_list=[]
        less_than_list=[]
        start=current_index
        current_index+=1
        while current_index<len(reader_obj) and 'M' not in reader_obj[current_index][0]:
            current_index+=1
        end=current_index
        type_of_lp,number_of_actual_variables, less_than_list, greater_than_list, equal_list, unit_cost=take_input_of_one_test_case(reader_obj[start:end],equal_list,greater_than_list,less_than_list)
        A,b,unit_cost,slack_start=transform_to_standard_lp(type_of_lp,number_of_actual_variables,less_than_list,greater_than_list,equal_list,unit_cost)
        if type(A)==bool and not A:
            test_case+=1
            continue
        print(f'Test Case {test_case} #####################################################################################################################################################')
        tableau=Tableau(A=A,b=b,unit_cost=unit_cost)
        test_case+=1
        tableau.solve()
        tableau.printMat()
        tableau.displayans(slack_start)

    


take_input(reader_obj=reader_obj)

    
