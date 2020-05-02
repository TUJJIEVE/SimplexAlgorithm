# The program is done by:
    * THUMATI UJJIEVE CS16BTECH11039
    * P RAMKISHAN     CS16BTECH11029
    * S SHANMUKHA RAO CS16BTECH11034

* use python 3+

* Please read this readme correctly to understand the functioning and assumptions made by the program.

* testcase.txt is provided. Please enter the values given in the txt file for checking for various cases if you want to use them.

* Please insert the correct input as asked by the program. Do not enter wrong inputs program will crash.

* The code is commented refer to it for any doubts.

* The code implements the algorithm as discussed in class

* The input for matrix A and B should ignore these condition:
    x_1>=0,x_2>=0,x_3>=0 .....
    These conditions are automatically inserted in the A and B matrix

* degeneracy is solved using bland's rule

* The degenerate.simplex.py file handles the degenerate examples and also outputs wether the solution
  is unbounded or infeasible.

* The simplex.py file handles non degenerate examples and also outputs wether the solution is unbounded
or infeasible

* Do not enter degenerate question in the simplex.py file the program will crash.

* The program implements simplex algorithm according to notes taught.

    AX <=B  and  maximize CX

## Note ** -> The input is assumed to exempt the conditions x_1,x_2,x_3,... >=0 so the vector A and B is assumed to contain non trivial conditions.
##            The trivial conditions are added automatically


## Pseudo code :
* choose x to be the initial feasible solution
* find A' and A'' given A'X = b' and A''X < b''   
* find neighbour vectors given by the columns of the negative of inverse of A'
* Check which neighbours will give greater cost if  C(x_i - x) > 0 and choose the neighbour that gives the max cost
* Then find the neighbour point given by x' = x + tv_i   v_i is the vector containing the direction of the selected neighbour
* to find the value of t 
            t = min over s (  (b_s - (A_s)x)/ ((A_s)v_i) )  s ranges over the rows of A'' and the corresponding b values
* set x = x' and repeat the steps until not neighbour that gives greater cost is found
* To solve degeneracy blands rule is used

Note :
* If all the values in the vector B are positive then [0,0,0,....] is chosen as initial feasible solution
* if any value in the vector B is negative two phase method is used which uses the same simplex algorithm more details given below 


## Details to tackle if [0,0,0...] is not the initial feasible solution:

    If any value in b vector is negative then [0,0,...] is not a feasible solution. In order to find the feasible solution
    we use two phase method in the first phase the initial feasible solution is found using the same simplex algorithm and in the second phase using the obtained solution
    of the previous phase as initial feasible solution we again run the simplex algorithm with original constraints. This gives the correct solution.

    For the first phase to work:

    we change the constraints.

    New objective function :->   maximize  -x_0   where x_0 is the artificial variable/dimension introduced

    New constraits :->    AX - X_0 <= B  where X_0 is the vector of shape (number of rows in A,1) so X_0 = [[x_0],[x_0],...] 
                        and  x_0 >= 0

        and the initial feasible solution for this phase is where x_1,x_2,... is set to zeros and x_0 is set to the absolute value of min(B)
