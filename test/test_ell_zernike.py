from aotools import functions
import matplotlib.pyplot as plt

if __name__ == '__main__':
    a = 1
    b = 0.8

    rmax = 7
    npix = 256
    l = 35

    ell_zern = functions.ZernikeEllipticalaperture(rmax, npix, a, b, l)

    Ell = ell_zern.CalculateEllipticalZernike()
    plt.imshow(Ell[2])
    plt.show()

    phi = ell_zern.EllZernikeMap()
    plt.imshow(phi)
    plt.show()
