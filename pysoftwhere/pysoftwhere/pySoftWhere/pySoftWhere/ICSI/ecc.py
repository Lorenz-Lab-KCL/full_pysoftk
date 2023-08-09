
###check this gives the same as mda

def inertia_tensor(pos,masses):
    tens = np.zeros((3, 3), dtype=np.float64)
    # xx
    tens[0][0] = (masses * (pos[:, 1] ** 2 + pos[:, 2] ** 2)).sum()
    # xy & yx
    tens[0][1] = tens[1][0] = - (masses * pos[:, 0] * pos[:, 1]).sum()
    # xz & zx
    tens[0][2] = tens[2][0] = - (masses * pos[:, 0] * pos[:, 2]).sum()
    # yy
    tens[1][1] = (masses * (pos[:, 0] ** 2 + pos[:, 2] ** 2)).sum()
    # yz + zy
    tens[1][2] = tens[2][1] = - (masses * pos[:, 1] * pos[:, 2]).sum()
    # zz
    tens[2][2] = (masses * (pos[:, 0] ** 2 + pos[:, 1] ** 2)).sum()

    return tens

def eccentricity(pos,masses):
    tens = inertia_tensor(pos,masses)

    return 1-3*min(np.linalg.eigvals(tens))/np.sum(np.linalg.eigvals(tens))
