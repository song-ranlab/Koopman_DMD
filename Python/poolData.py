import numpy as np
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