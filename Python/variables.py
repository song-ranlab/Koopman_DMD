
import numpy as np
import scipy
import scipy.linalg as linalg
#import control
#from control.matlab import *
dt = 0.01 # timestep size
l = 2.0 # arm length
m = 1.0 # pendulum mass
g = 9.81 # gravity acceleration
b = 0 #Damping Factor
duration = 20.0 # duration of simulation in seconds
thresh = 0.0 # DMD Threshold: 0 for no truncation, <0 sets the eigen value minimum threshold, >0 sets the rank as long as thresh is lower than max rank
x1_init = np.pi/4.0
x2_init = np.pi/4.0
x0 = [x1_init, x2_init]
#x0 = np.matrix([[x1_init],
#                [x2_init]]) # Initial conditions for pendulum simulation; theta, thetadot
x1_fin = 0.0
x2_fin = 0.0
#xf = np.matrix([[x1_fin],
#              [x2_fin]])
xf = [x1_fin, x2_fin]
steps = int(duration * 1.0 / dt)
t_span = np.linspace(0.0, duration, steps)

Qmod = 1.0
R = 10.0  # Process Cost
Qstate =([1.0, 0.1])
Q = np.diag(Qstate) # Control Cost

N = 20 # prediction horizon
Nu = 20 # Control Horizon
Ru = 1.0 # Control Cost

nvari = 2
order = 2
sine = 1

# linear system
A = np.matrix([[0.,1.],
               [-l/g,0.]])
B = np.matrix([[0.],
               [1./(m*l**2)]])


def dLQR(A,B,Q,R):
    #Solve the continuous time lqr controller.

    #dx/dt = A x + B u

    #cost = integral x.T*Q*x + u.T*R*u
    
    #ref Bertsekas, p.151

    #first, try to solve the ricatti equation
    P = np.matrix(scipy.linalg.solve_discrete_are(A, B, Q, R))
    # compute the LQR gain
    K = np.matrix(scipy.linalg.inv(B.T@P@B+R)@(B.T@P@A))

    eigVals, eigVecs = linalg.eig(A-B@K)
    return K, P, eigVals
K,X,E = dLQR(A,B,Q,R)
print(K)

K = [1.6902, 3.6785]