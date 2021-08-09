
# My own dmdc
import numpy as np
import scipy as sci


def myDMDcUB(X, U,r,q,dt):
    # X is the input state Matrix
    # U is the control input
    # r is the desired rank or parameters of the state reconstruction
    # q is the desired rank or parameters of the control reconstruction
    thresh = 1e-10
    n = int(len(X[1]))  # rank of augmented control data
    m = int(len(X[0]))
    Omega = np.zeros((n + q, m - 1))

    #Set up Snapshots
    X1 = X[:, :-1]
    X2 = X[:, 1:]
    Omega = np.concatenate((X1, U), axis=0)

    #Start SVD for Control
    Up, sp, Vp = sci.linalg.svd(Omega, lapack_driver='gesvd')

    if r == 0 :
        r = len(np.where(sp > -thresh))
    elif r < 0 :
        r = n
    elif r > 0:
        if r > n:
            r = n

    Vp = Vp.T
    Vp = Vp[:, :r]
    Sp = np.diag(np.reciprocal(sp[:r]))


    #SVD for State Evolution
    Ur, sr, Vr = sci.linalg.svd(X2, lapack_driver='gesvd')
    if r == 0 :
        r = len(np.where(sr > -thresh))
    elif r < 0 :
        r = n
    elif r > 0:
        if r > n:
            r = n
    Ur = Ur[:,0:r]
    Up1 = Up[:r, :r]
    Up2 = Up[r:, :r]
    #Up2 = Up[r+1:r+q, :]
    Atilde = np.linalg.multi_dot([Ur.T.conj(), X2, Vp,
                                  Sp, Up1.T.conj(), Ur])
    Btilde = np.linalg.multi_dot([Ur.T.conj(), X2, Vp,
                                  Sp, Up2.T.conj()])

    [W, D] = np.linalg.eig(Atilde)

    #Phi = Xp * Vtil * inv(Sigtil) * U_1.T *Uhat * W;
    Phi = np.linalg.multi_dot([X2,Vp, Sp, Up1.T.conj(), W])
    t = np.array(0, m-1)
    tspan =


    return[Atilde, Btilde, Phi, Xdmd]


def myDMDcKB(X,U,,B,r,dt):

    n = int(len(X[1]))  # rank of augmented control data
    m = int(len(X[0]))

    #Set up Snapshots
    X1 = X[:, :-1]
    X2 = X[:, 1:]

    U, Sig, V = sci.linalg.svd(X2, lapack_driver='gesvd')

    U = U[0:r, :].T
    sig = np.diag(np.reciprocal(Sig[:r]))  # Figure out later why python stupid
    V = V[0:r, :].T
    holder = X2 - B * U
    A_DMDc = np.linalg.multi_dot([holder, V, sig, U.T.conj()])
    print(A_DMDc)

    return [Atilde, Xdmd]