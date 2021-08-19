import numpy as np
import scipy as sci
import matplotlib.pyplot as plt

# This function builds a linear system Y = Ax + Bu for whatever timespan T_span

def buildlinsys (A, B, X0, U, t_span):

    n = np.size(A,1)
    m = np.size(t_span,0)
    Y = np.zeros((n,m-1))
    Y[:, 0] = X0
    for i in range(0, m - 2):

        Y[:, i + 1] = A @ Y[:, i] + B @ U[:, i]

    return[Y]


#Plots comparisons between two sets of data assuming each row corresponds to the same state
#wrapping to Pi (or 2 Pi) shoul be done prior to this function

def comparestates (X1, X2, t_span, title_main, title1, title2, figsize, dpi ):
    n = len(X1)
    m = len(X2[1])
    # figsize = (22, 25)
    # dpi = 600

    plt.figure(figsize=figsize, dpi=dpi)

    for i in range(1, n + 1):
        plt.subplot(n, 1, i)
        plt.title(title_main + " State " + str(i))
        label1 = title1 + " State " + str(i)
        label2 = title2 + " State " + str(i)
        plt.plot(t_span[0:m], X1[i - 1, 0:m], 'r-', label=label1)
        plt.plot(t_span[0:m], X2[i - 1, 0:m], 'b--', label=label2)
        plt.xlabel(r"$t$ [s]")
        plt.ylabel(" State " + str(i))
        plt.legend(loc='best')

def showstate (x, t_span, title, figsize, dpi):
    n = len(x)
    m = len(x[1])
    plt.figure(figsize=figsize, dpi=dpi)
    for i in range(1, n + 1):
        plt.subplot(n, 1, i)
        plt.title(title + " State " + str(i))
        plt.plot(t_span[0:m], x[i - 1, 0:m], 'r-', label=r'State'+str(i))
        plt.xlabel(r"$t$ [s]")
        plt.ylabel(" State " + str(i))
        plt.legend(loc='best')


#Takes in a matrix where the X axis is the time step and an Y Axis is the individual states

# Here we will make our 'lifting' function. In this example I will be using a pool data function from
#  ["Discovering Governing Equations from Data: Sparse Identification of Nonlinear Dynamical Systems" by S. L. Brunton,
#  J. L. Proctor, and J. N. Kutz](10.1073/pnas.1517384113)
#
# # Note for later
# Double check the higher order sizing aspect, worlframalpha might be wrong (or just put a cap at order 5?)
# Also look back into the Fouier xpac, it looks right and tests right at order 1 not too sure after..
############################################################################################################################
# yin - input n x m matrix
# nVars - Number of independent variables should multiple columns be some prior combination of original states

def pool(yin, nVars, porder, sorder):
    if porder == 2:
        col = nVars + (nVars*(nVars+1)/2)
    elif porder ==3:
        col = nVars + (nVars*(nVars+1)/2)+(nVars*(nVars+1)*(nVars+2)/(2*3))
    elif porder >= 4:
        col = nVars + (nVars*(nVars+1)/2)+(nVars*(nVars+1)*(nVars+2)/(2*3))+11 #might need to double check this algo l8r but simple test says A.OK
    else:
        col = nVars

    if sorder > 0:
        col = col + sorder * 2

    yout = np.zeros((len(yin),int(col)))
    ind = 0
    for i in (0,nVars-1):#    poly order 1
        #print(ind)
        yout[:,ind] = yin[:,i]
        ind = ind+1
        #print('i is {}'.format(i))
        if ind >= (col-1):
            print('For loops failed, breaking out of this')
            break
    if porder >= 2: # poly order 2
        for i in (0,nVars-1):
            for j in (i,nVars-1):
                #print('i is {}'.format(i))
                #print('j is {}'.format(j))
                #print(ind)
                yout[:,ind] = yin[:,i]*yin[:,j]
                ind = ind+1
                #print(ind+10)
                if ind >= (col-1):
                    print('For loops failed, breaking out of this')
                    break
    if porder >=3: # poly order 3
        for i in (0,nVars-1):
            for j in (i,nVars-1):
                for k in (j,nVars-1):
                    yout[:,ind] = yin[:,i] * yin[:,j] * yin[:,k]
                    ind = ind+1
                    if ind >= (col-1):
                        print('For loops failed, breaking out of this')
                        break
    if porder>=4: # poly order 4
        for i in (0,nVars-1):
            for j in (i,nVars-1):
                for k in (j,nVars-1):
                    for q in (k,nVars-1):
                        yout[:,ind] = yin[:,i]*yin[:,j]*yin[:,k]*yin[:,q]
                        ind = ind+1
                        if ind >= (col-1):
                            print('For loops failed, breaking out of this')
                            break
    if sorder >= 1:
        for k in (0,sorder-1):
            yout = np.matrix[yout, np.sin(k*yin), np.cos(k*yin)]
    return yout