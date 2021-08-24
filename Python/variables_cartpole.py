import numpy as np
import scipy
import scipy.linalg as linalg

# import control
# from control.matlab import *
dt = 0.01  # timestep size
l = 2.0  # arm length
m = 1.0  # pendulum mass
M = 5.0
g = 9.81  # gravity acceleration
b = 1 #above or below axis
d= 0 # Damping Factor
duration = 20.0  # duration of simulation in seconds
thresh = 0.0  # DMD Threshold: 0 for no truncation, <0 sets the eigen value minimum threshold, >0 sets the rank as long as thresh is lower than max rank
x1_init = -1.0
x2_init = 0.0
x3_init = np.pi + 0.00001
x4_init = 0.0
x0 = [x1_init, x2_init, x3_init, x4_init]
# x0 = np.matrix([[x1_init],
#                [x2_init]]) # Initial conditions for pendulum simulation; theta, thetadot
x1_fin = 1.0
x2_fin = 0.0
x3_fin = np.pi
x4_fin = 0.0
# xf = np.matrix([[x1_fin],
#              [x2_fin]])
xf = [x1_fin, x2_fin, x3_fin, x4_fin]
steps = int(duration * 1.0 / dt)
t_span = np.linspace(0.0, duration, steps)

Qmod = 1.0
R = 10.0  # Process Cost
Qstate = ([1.0, 10, 1.0, 1.0])
Q = np.diag(Qstate)  # Control Cost

N = 20  # prediction horizon
Nu = 20  # Control Horizon
Ru = 1.0  # Control Cost

nvari = 4
order = 2
sine = 0

# linear system
A = np.matrix([[0, 1, 0, 0],
               [0, -d/M, b*m*g/M, 0],
               [0,0,0,1],
               [0, -b*d/(M*l),-b*(m+M)*g/(M*l), 0]])
B = np.matrix([[0],
               [1/M],
               [0],
               [b/(M*l)]])

"""
def LQR(A,B,Q,R):
    #Solve the continuous time lqr controller.

    #dx/dt = A x + B u

    #cost = integral x.T*Q*x + u.T*R*u

    #ref Bertsekas, p.151

    #first, try to solve the ricatti equation
    X = np.matrix(linalg.solve_continuous_are(A, B, Q, R))

    #compute the LQR gain
    K = np.matrix(linalg.inv(R)*(B.T*X))

    eigVals, eigVecs = linalg.eig(A-B*K)

    return K, X, eigVals
K,X,E = LQR(A,B,Q,R)
print(K)
"""

K = np.array([-10.0, -40.7584, 512.9406, 199.6643],ndmin=2)
#K = np.array([-.3162, -2.1850, 140.814, 57.35],ndmin=2)