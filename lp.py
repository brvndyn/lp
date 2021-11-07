## Brandyn Butchard  ##
## V00906030         ##
## CSC 445 LP Solver ##

import sys
import numpy as np
import math


global epsilon
epsilon = 0.00001


def parse():

    linenum = 0
    c_row = []
    constraints = []

    for line in sys.stdin:
        line = line.split()
        if line:
            if linenum == 0:
                c_row = line
            else:
                constraints.append(line)
            linenum = linenum+1

    return c_row, constraints



def aux(a,b,original_c,l):
    for i in range(len(b)):
        if b[i] == l:
            leaving_index = i

    leaving = []
    for i in range(len(a[leaving_index])):
        leaving.append(-a[leaving_index][i])
    c=[]



    for i in range(len(leaving)):
        #leaving[i] = -1*leaving[i]
        c.append(leaving[i])



    for i in range(len(a)):
        b[i] = b[i] - l
        for j in range(len(a[i])):
            a[i][j] = -a[i][j]
            a[i][j] = a[i][j] - leaving[j]

    for i in range(len(a)):
        for j in range(len(a[i])):
            a[i][j] = -a[i][j]


    obj_val = l
    b[leaving_index] = -l
    a[leaving_index] = leaving

    c.append(-1)
    for i in range(len(a)):
        a[i].append(-1)



    nonbasic_variables = np.array([[0]*len(c)])     # keeps track of position of original nb variables (x1, x2, ..., xn)
    row = nonbasic_variables
    for i in range(len(b)):
        nonbasic_variables = np.vstack([nonbasic_variables, row])


    for i in range(len(nonbasic_variables[0])):
        nonbasic_variables[0][i] = i+1

    nonbasic_variables[leaving_index+1][0] = -1     # Adding omega to the nonbasic variables tracker


    # Use the simplex function to solve the auxilliary problem
    a,b,c,obj_val,nonbasic = simplex(c,a,b,obj_val,nonbasic_variables)


    if obj_val > 0+epsilon or obj_val < 0-epsilon:
        print("infeasible")
        quit()

    omega = -1
    for i in range(len(nonbasic_variables[0])):
        if nonbasic_variables[0][i] == -1:
            omega = i                               # Finding omega to be removed




    # If omega is in the basis, pivot it out
    if omega == -1:
        for i in range(len(nonbasic_variables)):
            if nonbasic_variables[i][0] == -1:
                leaving_linenum = i-1


        leaving = a[leaving_linenum]

        # Pivot omega on the first non zero coefficient in it's row
        for i in range(len(leaving)):
            if leaving[i] != 0:
                omega = i
                break

        temp = nonbasic_variables[leaving_linenum+1][0]
        nonbasic_variables[leaving_linenum+1][0] = nonbasic_variables[0][omega]
        nonbasic_variables[0][omega] = temp


        objective_increase = c[omega]
        c[omega] = c[omega]/leaving[omega]
        final_c_omega = -c[omega]
        divisor = leaving[omega]
        b_leaving = b[leaving_linenum]/divisor


        # Pivot
        for i in range(len(leaving)):
            leaving[i] = leaving[i]/divisor

        for i in range(len(c)):
            c[i] = c[i] - leaving[i]*objective_increase

        linenum = 0
        for line in a:
            b[linenum] = b[linenum] - b_leaving*line[omega]
            a_omega = a[linenum][omega]
            for i in range(len(line)):
                if linenum != leaving_linenum:
                    a[linenum][i] = a[linenum][i] - a_omega*leaving[i]

            if linenum == leaving_linenum:
                a[linenum][omega] = a_omega/divisor
            else:
                a[linenum][omega] = -a_omega/divisor
            linenum = linenum+1


        b[leaving_linenum] = b_leaving
        obj_val = obj_val + objective_increase*b[leaving_linenum]
        c[omega] = final_c_omega







    c[omega] = 0
    for linenum in range(len(a)):
        a[linenum][omega] = 0



    for i in range(len(nonbasic_variables)):
        for j in range(len(nonbasic_variables[i])):
            if nonbasic_variables[i][j] == -1 or nonbasic_variables[i][j] > len(original_c):
                nonbasic_variables[i][j] = 0



    # Restore the c vector and objective value
    for i in range(len(nonbasic_variables[0])):
        if nonbasic_variables[0][i]:
            v = nonbasic_variables[0][i]-1
            c[i] = original_c[v]

    for i in range(len(c)):
        for linenum in range(1, len(nonbasic_variables)):
            if nonbasic_variables[linenum][0]:
                v = nonbasic_variables[linenum][0]-1
                #print(original_c[v]," * ", a[linenum-1][i])
                c[i] += original_c[v]*-a[linenum-1][i]
                #print("c",i," = ",c[i])


    for linenum in range(1,len(nonbasic_variables)):
        if nonbasic_variables[linenum][0]:
            v = nonbasic_variables[linenum][0]-1
            obj_val += original_c[v]*b[linenum-1]





    return a,b,c,obj_val,nonbasic_variables




def simplex(c,a,b,obj_val,nonbasic_variables):

    while True:
        #blands rule
        entering = -1
        for i in range(0,len(c)):
            if c[i] > 0+epsilon:  # some floating point errors are solved by checking above an epsilon value, I defined it globally as 0.0000001
                entering = i
                break

        # if there's no more entering variabls, we're at the optimal
        if entering == -1:
            # Return everything so it can be used for auxilliary problems as well
            return a,b,c,obj_val,nonbasic_variables

        #print("entering: ", entering)


        linenum = 0
        tightest_bound = 10000000000000000000000000000000000000000
        leaving = None
        unbounded = True
        for line in a:
            # Check upper bounds
            if line[entering] > 0:
                ratio = b[linenum]/line[entering]
                unbounded = False
            else:
                ratio = 100000000000000000000000000000000000000000 # ignores 0 coefficient


            if ratio < tightest_bound:
                tightest_bound = ratio
                leaving = line
                leaving_linenum = linenum
            linenum = linenum+1

        if unbounded == True:
            return(a,b,c,"unbounded",nonbasic_variables)

        # Swap the variable locations in the nonbasic_variables list
        temp = nonbasic_variables[leaving_linenum+1][0]
        nonbasic_variables[leaving_linenum+1][0] = nonbasic_variables[0][entering]
        nonbasic_variables[0][entering] = temp


        # Save some values for later use
        objective_increase = c[entering]
        c[entering] = c[entering]/leaving[entering]
        final_c_entering = -c[entering]
        divisor = leaving[entering]
        b_leaving = b[leaving_linenum]/divisor


        # Pivot
        for i in range(len(leaving)):
            leaving[i] = leaving[i]/divisor

        for i in range(len(c)):
            c[i] = c[i] - leaving[i]*objective_increase

        linenum = 0
        for line in a:
            b[linenum] = b[linenum] - b_leaving*line[entering]
            a_entering = a[linenum][entering]
            for i in range(len(line)):
                if linenum != leaving_linenum:
                    a[linenum][i] = a[linenum][i] - a_entering*leaving[i]

            if linenum == leaving_linenum:
                a[linenum][entering] = a_entering/divisor
            else:
                a[linenum][entering] = -a_entering/divisor
            linenum = linenum+1


        b[leaving_linenum] = b_leaving
        obj_val = obj_val + objective_increase*b[leaving_linenum]
        c[entering] = final_c_entering






def main():
    c,constraints = parse()


    # Prepare input data
    a = []
    b = []
    for i in constraints:
        a.append(i[0:len(i)-1])
        b.append(float(i[len(i)-1]))


    for i in range(len(c)):
        c[i] = float(c[i])

    linenum = 0
    for i in a:
        for j in range(len(i)):
            a[linenum][j] = float(a[linenum][j])
        linenum = linenum+1



    nonbasic_variables = np.array([[0]*len(c)]) # keeps track of position of original nb variables (x1, x2, ..., xn)
    row = nonbasic_variables
    for i in range(len(b)):
        nonbasic_variables = np.vstack([nonbasic_variables, row])


    for i in range(len(nonbasic_variables[0])):
        nonbasic_variables[0][i] = i+1

    obj_val = 0

    nonbasic_variables = np.array([[0]*len(c)]) # keeps track of position of original nb variables (x1, x2, ..., xn)
    row = nonbasic_variables
    for i in range(len(b)):
        nonbasic_variables = np.vstack([nonbasic_variables, row])


    for i in range(len(nonbasic_variables[0])):
        nonbasic_variables[0][i] = i+1

    lowest_value = 0
    for i in b:
        if i < lowest_value-epsilon:
            lowest_value = i
    if lowest_value != 0:
        a,b,c,obj_val,nonbasic_variables = aux(a,b,c,lowest_value)




    # Let it rip
    a,b,c,obj_val,nonbasic = simplex(c,a,b,obj_val,nonbasic_variables)


    if obj_val == "unbounded":
        print("unbounded")
        quit()

    # Prepare output + output

    optimal_assignments = [0]
    for i in range(len(nonbasic)):
        for j in nonbasic[i]:
            if j:
                optimal_assignments.append(0)
    #optimal_assignments = [0]*(len(c)+1)
    for i in range((len(b))):
        v = nonbasic[i+1][0]
        if v:
            optimal_assignments[v] = b[i]


    print("optimal")
    if obj_val == int(obj_val):
        print(int(obj_val))
    else:
        print("%.7g" % obj_val)

    optimal_assignments = optimal_assignments[1:]
    for assignment in optimal_assignments:
        if assignment == int(assignment):
            print(int(assignment),end =" ")
        else:
            print("%.7g" % assignment,end =" ")

if __name__ == '__main__':
    main()
