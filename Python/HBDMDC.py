
# My own dmdc
import numpy as np
import scipy as sci


def mydmd(X, r, dt):
    # n = int(len(X[1]))  # rank of augmented control data
    # m = len(X[0])
    m = np.size(X, 1)
    X1 = X[:, :-1]
    X2 = X[:, 1:]
    U, Sig, V = sci.linalg.svd(X1, lapack_driver='gesvd')

    sig = np.diag(Sig)
    V = V.T
    U = U[0:r, :]
    V = V[0:r,:]

    # U = np.array([[-1, 1]]) * U  # I am totally going to cheat my ass off here with a few commands - lets look back at this today
    # #print(U)
    # sig = np.diag(np.reciprocal(Sig[:r]))  # Figure out later why python stupid
    # #V = V[0:r, :].T
    # V = np.array([[-1, 1]]) * V # I am totally going to cheat my ass off here with a few commands
    #print(V)
    Atilde = np.linalg.multi_dot([U.T.conj(), X2, V, sig])
    print(Atilde)
    [W, D] = np.linalg.eig(Atilde)
    # Compute the dynamic modes of operator Atilde
    Phi = np.linalg.multi_dot([X2, V, sig, U.T.conj(), W])
    #Python sometimes gets really large eigen values and cannot continue so try-catch to just brute force it
    try:
        t = np.array((0, m - 1))
        tspan = np.linspace(0., dt * (m - 1), m - 1)

        lambdaval = np.diag(D)
        omega = np.log(lambdaval) / dt
        hold = len(omega)
        omega.resize((hold, 1), refcheck=False)
        x1 = X[:, 0]
        b = np.linalg.lstsq(Phi,x1)
        print(b)
        hold = len(b)
        b.resize((hold, 1), refcheck=False)
        mm1 = m - 1
        time_dynamics = np.zeros((r, mm1))
        for i in range(0, mm1 - 1):
            check = b * (np.exp(omega * tspan[i]))
            time_dynamics[:, i] = check
            print(check)
        Xdmd = Phi * time_dynamics
    except:
        Xdmd = 0.
        print('Something happened calculating Xdmd')

    return [Atilde, Phi, Xdmd]


def myDMDcUB(X, U,r,q,dt):
    # X is the input state Matrix
    # U is the control input
    # r is the desired rank or parameters of the state reconstruction
    # q is the desired rank or parameters of the control reconstruction
    thresh = 1e-10
    n = np.size(X, 0)  # rank of augmented control data
    m = np.size(X, 1)
    X1 = X[:, :-1]
    X2 = X[:, 1:]
    Omega = np.concatenate((X1, U), axis=0)
    #
    #     #Start SVD for Control
    Up, sp, Vp = sci.linalg.svd(Omega, lapack_driver='gesvd')
    if r == 0:
        r = len(np.where(sp > -thresh))
    elif r < 0:
        r = n
    elif r > 0:
        if r > n:
            r = n
    Vp = Vp.T
    Vp = Vp[:, :r]
    Sp = np.diag(np.reciprocal(sp[:r]))
    Ur, sr, Vr = sci.linalg.svd(X2, lapack_driver='gesvd')
    #This part of Kutz code but it feels a little redundant and possibly interfering with future iterations
    # if r == 0:
    #     r = len(np.where(sr > -thresh))
    # elif r < 0:
    #     r = n
    # elif r > 0:
    #     if r > n:
    #         r = n
    Ur = Ur[:, 0:r]
    Up1 = Up[:r, :r]
    Up2 = Up[r:, :r]
    Atilde = np.linalg.multi_dot([Ur.T.conj(), X2, Vp,
                                   Sp, Up1.T.conj(), Ur])
    Btilde = np.linalg.multi_dot([Ur.T.conj(), X2, Vp,
                                   Sp, Up2.T.conj()])

    [W, D] = np.linalg.eig(Atilde)

    # Compute the dynamic modes of operator Atilde

    Phi = np.linalg.multi_dot([X2, Vp, Sp, Up1.T.conj(), W])

    try:
        t = np.array((0, m - 1))
        tspan = np.linspace(0., dt * (m - 1), m - 1)

        lambdaval = np.diag(D)
        omega = np.log(lambdaval) / dt
        omega.resize((r, 1), refcheck=False)
        x1 = X[:, 0]
        b = np.linalg.lstsq(Phi,x1)
        b.resize((r, 1), refcheck=False)
        mm1 = m - 1
        time_dynamics = np.zeros((r, mm1))
        for i in range(0, mm1 - 1):
            check = b * (np.exp(omega * tspan[i]))
            time_dynamics[:, i] = check

        Xdmd = Phi * time_dynamics
    except:
        Xdmd = 0.

    return[Atilde, Btilde, Phi, Xdmd]


def myDMDcKB(X,U,B,r,dt):

    n = int(len(X[1]))  # rank of augmented control data
    m = int(len(X[0]))

    #Set up Snapshots
    X1 = X[:, :-1]
    X2 = X[:, 1:]

    U, Sig, V = sci.linalg.svd(X2, lapack_driver='gesvd')

    U = U[0:r, :]
    sig = np.diag(np.reciprocal(Sig[:r]))  # Figure out later why python stupid
    V = V[0:r, :].T
    holder = X2 - B * U
    Atilde = np.linalg.multi_dot([holder, V, sig, U.T.conj()])
    print(Atilde)
    [W, D] = np.linalg.eig(Atilde)

    # Compute the dynamic modes of operator Atilde
    Phi = np.linalg.multi_dot([X2, V, sig, U.T.conj(), W])
    try:
        t = np.array((0, m - 1))
        tspan = np.linspace(0., dt * (m - 1), m - 1)

        lambdaval = np.diag(D)
        omega = np.log(lambdaval) / dt
        omega.resize((r, 1), refcheck=False)
        x1 = X[:, 0]
        b = np.linalg.lstsq(Phi,x1)
        b.resize((r, 1), refcheck=False)
        mm1 = m - 1
        time_dynamics = np.zeros((r, mm1))
        for i in range(0, mm1 - 1):
            check = b * (np.exp(omega * tspan[i]))
            time_dynamics[:, i] = check

        Xdmd = Phi * time_dynamics
    except:
        Xdmd = 0.
    return [Atilde, Xdmd]