"""
    Program by:
        
        THUMATI UJJIEVE CS16BTECH11039
        P RAMKISHAN     CS16BTECH11029
        S SHANMUKHA RAO CS16BTECH11034
"""

import numpy as np
import logging
print("Please ignore the conditions like x_1,x_2,x_3,... >=0")
print("see test_case.txt file to see how the A and B matrix are represented")
print("Please don't include the above said constraints in the vector A and B but you can include constraints like 3x_1 <= 2 in A and B matrix")

print("Enter Matrix A")
m = int(input("Enter the number of rows in A:"))
n = int(input("Enter the number of columns in A:"))
A = []
for i in range(m):
    temp = []
    for j in range(n):
        temp.append(int(input("enter element A_{}{}:".format(i,j))))
    A.append(temp)
print("Enter for B")
b = []
for i in range(m):
    b.append(int(input("enter element B_{}:".format(i))))
print("Enter for C")
c = []
for i in range(n):
    c.append(int(input("enter element C_{}:".format(i))))

A = np.array(A)
b = np.array(b)
c = np.array(c)

"""

    The program implements simplex algorithm according to notes taught.

        AX <=B  and  maximize CX

    Note ** -> The input is assumed to exempt the conditions x_1,x_2,x_3,... >=0 so the vector A and B is assumed to contain non trivial conditions
                The trivial conditions are added automatically


    Pseudo code :
    * choose x to be the initial feasible solution
    * find A' and A'' given A'X = b' and A''X < b''   
    * find neighbour vectors given by the columns of the negative of inverse of A'
    * Check which neighbours will give greater cost if  C(x_i - x) > 0 and choose the neighbour that gives the max cost
    * Then find the neighbour point given by x' = x + tv_i   v_i is the vector containing the direction of the selected neighbour
    * to find the value of t 
             t = min over s (  (b_s - (A_s)x)/ ((A_s)v_i) )  s ranges over the rows of A'' and the corresponding b values
    * set x = x' and repeat the steps until not neighbour that gives greater cost is found


    To solve degeneracy blands rule is used
    
    If all the values in the vector B are positive then [0,0,0,....] is chosen as initial feasible solution
    if any value in the vector B is negative two phase method is used which uses the same simplex algorithm more details given below    


"""

cond = input("Do you want to pring DEBUG statements 0 for False and 1 for True:")

if int(cond):
    logging.basicConfig(level=0)


def SolveSimplex(A,ext_point,b,c,phase1 = False):
    """ 
        Finding A' and A''. A' contains the rows of A for which A'X_0 = b' and the other rows are present in A''. 
    
    """
    iter = 1
    old_A_dash = np.zeros(shape = (A.shape[1],A.shape[1]))
    while iter < 100:
        logging.debug("iteration no.: {}".format(iter))
        logging.debug("current extreme point: {}".format(ext_point))

        ## Finding A' and A'' and the corresponding rows in b given by b' and b''.
        new_A_ind = np.matmul(A,ext_point.reshape(-1,1)) - b.reshape(-1,1) 
        new_A_ind = np.around(new_A_ind,decimals=5)
        inds = ((new_A_ind == 0).reshape(-1,))
        A_dash = A[inds,:]
        # print(A_dash)
        if A_dash.shape[1] != A_dash.shape[0]:
            ## Extra basis, leave one basis
            if iter !=1:
                for a in old_A_dash.tolist():
                    # print(a)
                    # print(A_dash.tolist())
                    if a in A_dash.tolist():
                        A_dash = np.delete(A_dash,A_dash.tolist().index(a),axis = 0)
                        break
            else:
                A_dash = np.delete(A_dash,0,axis = 0)
        # print(A_dash)
        old_A_dash = A_dash
        A_ddash = np.array([x for x,i in zip(A,range(A.shape[0])) if not inds[i] ])
        b_dash = (b.reshape(-1,1)[inds,:]).reshape(-1,)
        b_ddash = np.array([ b[x] for x in range(len(b)) if not inds[x]])
        ## Finding the direction of the neighbours by calculating the inverse of the negative of A' matrix.
        ## The columns of dir_v gives the direction of neighbours
        dir_v = -np.linalg.inv(A_dash)
        # logging.debug("columns of the below matrix represent the direction of neighbours of {}:".format(ext_point))
        # logging.debug("{}".format(dir_v))

        # ## Finding which neighbours will increase the objective function
        cost = np.matmul(c,dir_v)
        ## Choosing the neighbour which gives the most cost
        indi = np.argmax(cost)
        max_val = np.max(cost)
        indi = [ i for i in range(len(cost)) if cost[i] == max_val ]
        logging.debug("number of neighbours that gives the greater cost : {}".format(len(indi)))
        selc_indi = indi[0]
        # print(cost)
        if cost[selc_indi] <=0:
            if cost[selc_indi] == 0 and not phase1:
                logging.debug("*****NOTE: Another feasible solution of equivalent cost exists")
            if not phase1:
                logging.debug("breaking optimum reached")
            break
        dir_neigh = dir_v[:,selc_indi] 
        logging.debug("selected dir: {}".format(dir_neigh))
        
        
        """ Unbounded checking  """
        
        """ 
            (A'')V_i  V_i is the direction vector of selected neighbour is being calculated
             and only the positive rows are being selected 
        """
        asvi = np.matmul(A_ddash,dir_neigh.reshape(-1,1))
        cond = (asvi>0).reshape(-1,)
        req_asvi = asvi[cond,:].reshape(-1,1)
        if len(req_asvi) == 0:
            #print("Unbounded")
            return -1
        reqb_ddash = (b_ddash.reshape(-1,1))[cond,:]
        reqa_ddash = A_ddash[cond,:]
        
        ## finding (b_s - ((A_s)x_0)) s ranges over all the rows of A which belong to A'' and for which  (A_s)V_i > 0
        aa = reqb_ddash - np.matmul(reqa_ddash,ext_point.reshape(-1,1)) 
        ## Finding the value of t to find the neighbour point x_0 <-  x_0 + t(V_i)  | V_i is the vector of the neigbour direction
        t = aa/req_asvi
        # print(t)
        t[t<0] = 0
        t = np.min(t)
        if t == 0:
            """ This condition will not occur given assumptions """
            logging.debug("Cannot move further ERROR")
            return -2

        logging.debug("Moving {} units in that direction".format(t))
        ## Getting the next extreme point which is the neighbour of the current extreme point
        ext_point = ext_point + (t * dir_neigh.reshape(-1,))
        # print("new_point",ext_point)
        logging.debug("---------X----------")
        iter+=1
    return ext_point

## Used to add constraints x1>=0 x2>=0 x3>=0 .... 
neg_I = -1 * np.identity(A.shape[1])
A = np.vstack((A,neg_I))
temp = np.zeros(shape = (A.shape[1],))
b = np.concatenate((b,temp))
ext_point = np.zeros(shape=(A.shape[1],))

"""
If any value in b vector is negative then [0,0,...] is not a feasible solution. In order to find the feasible solution
we use two phase method in the first phase the initial feasible solution is found and in the second phase using the obtained solution 
we again run the simplex algorithm with original constraints. This gives the correct solution.

For the first phase to work:

we change the constraints.

New objective function :->   maximize  -x_0   where x_0 is the artificial variable/dimension introduced

New constraits :->    AX - X_0 <= B  where X_0 is the vector of shape (number of rows in A,1) so X_0 = [[x_0],[x_0],...] 
                      and  x_0 >= 0

    and the initial feasible solution is where x_1,x_2,... is set to zeros and x_0 is set to the absolute value of min(B)

"""

logging.debug("Finding initial feasible solution")

""" When [0,0,0...] is not initial feasible solution """
if np.min(b) < 0:
    """ Changing the A,B,C matrix to the new objective function and new constraints """
    c_temp = np.zeros(shape = (len(c)+1,))
    c_temp[-1] = -1
    dims= A.shape[1]
    temp = -1 * np.ones(shape=(A.shape[0] - dims,1),dtype = int)
    tempA = np.concatenate((A[:-dims,:],temp),axis = 1)
    tempAA = np.concatenate((A[A.shape[0]-dims:,:],np.zeros(shape=(A.shape[0] - dims,1))),axis = 1)
    tempAAA = np.zeros(shape=(A.shape[1] +1,))
    tempAAA[-1] = -1
    A_temp = np.vstack((tempA,tempAA,tempAAA))
    b_temp = np.concatenate((b,np.array([0])))
    ext_point = np.zeros(shape = (A.shape[1] + 1,))
    ext_point[-1] = abs(np.min(b))

    logging.debug("### PHASE 1 to find the initial feasible solution")
    ext_point = SolveSimplex(A_temp,ext_point,b_temp,c_temp,True)

    if ext_point[-1] !=0:
        print("Infeasible Solution")
    else:
        print("The initial feasible solution occurs at : {}".format(ext_point))
        logging.debug("### Phase 1 completed")
        logging.debug("### Phase 2 \n--------------X-----------")
        solution = SolveSimplex(A,ext_point[:-1],b,c)
        if type(solution) == int and solution == -1:
            print ("Unbounded solution")
        else:
            print("optimum solution",solution,"optimal_value:",np.matmul(solution,c.reshape(-1,1)))

else:
    """ When [0,0,0...] is the feasible solution """

    print("The initial feasible solution occurs at : {}".format(ext_point))
    solution = SolveSimplex(A,ext_point,b,c)
    if type(solution) == int and solution == -1:
        print("Unbounded solution")
    else:
        print("Optimum solution occurs at:",solution,"optimal_value:",np.matmul(solution,c.reshape(-1,1)))
