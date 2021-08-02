
# My own dmdc


def myDMDcKB(X, U,A_rank,B_rank):

    Omega = np.zeros((3, m - 1))
    for i in range(0, m - 1):
        # print(i)
        Omega[0, i] = X[0, i]
        Omega[1, i] = X[1, i]
        Omega[2, i] = u[0, i]

    Omega2 = np.vstack([X[:, :m - 1], u])
    Up, sp, Vp = sci.linalg.svd(Omega, lapack_driver='gesvd')
    Vp = Vp.T
    Vp = Vp[:, :r]
    Sp = np.diag(np.reciprocal(sp[:r]))
    Up1 = Up[:r, :r]
    Up2 = Up[r:, :r]
    Ur, sr, Vr = sci.linalg.svd(X2, lapack_driver='gesvd')

    Atilde = np.linalg.multi_dot([Ur.T.conj(), X2, Vp,
                                  Sp, Up1.T.conj(), Ur])
    Btilde = np.linalg.multi_dot([Ur.T.conj(), X2, Vp,
                                  Sp, Up2.T.conj()])

    s1 = np.zeros((1, m))
    s2 = np.zeros((1, m))
    s1[0, 0] = 4
    s2[0, 0] = 7
    print('This is Atilde', Atilde)
    print('this is Btilde:', Btilde)

    for i in range(0, 4):
        useme = np.array([[Xplot1[0, i]], [Xplot1[1, i]]])
        checkme = np.matmul(Atilde, useme) - Btilde * u
        s1[0, i + 1] = checkme[0, 0]
        s2[0, i + 1] = checkme[1, 0]

    return [x_1, x_2+ctrl]


def myDMDcUB(t,xk):
    x_1 = xk[1]
    x_2 = -b/m * xk[1] - g/l * np.sin(xk[0])
    ctrl = 1/(m*l**2) * u(xk)
    #print(ctrl)
    return [x_1, x_2+ctrl]