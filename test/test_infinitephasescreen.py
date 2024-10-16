from aotools.turbulence import infinitephasescreen
import numpy as np

def testVKInitScreen():

    scrn = infinitephasescreen.PhaseScreenVonKarman(128, 4./64, 0.2, 50, n_columns=4)

def testVKAddRow():

    scrn = infinitephasescreen.PhaseScreenVonKarman(128, 4./64, 0.2, 50, n_columns=4)
    scrn.add_row()

def testLeastSquaresSolver():
    airmass = 1.0 / np.cos(30.0 / 180. * np.pi)
    r0 = 0.9759 * 0.5 / (0.7 * 4.848) * airmass ** (-3. / 5.)
    r0wavelength = r0 * (500 / 500.0) ** (6. / 5.)
    screen = infinitephasescreen.PhaseScreenVonKarman(120, 1./120, r0wavelength, 50, 1, n_columns=2)
    mean_scrn = np.sqrt(np.mean(screen.scrn**2))
    for i in range(500):
        screen.add_row()
    assert(0.5 <= mean_scrn/np.sqrt(np.mean(screen.scrn**2)) <= 1.5)


# Test of Kolmogoroc screen
def testKInitScreen():

    scrn = infinitephasescreen.PhaseScreenKolmogorov(128, 4./64, 0.2, 50, stencil_length_factor=4)

def testKAddRow():

    scrn = infinitephasescreen.PhaseScreenKolmogorov(128, 4./64, 0.2, 50, stencil_length_factor=4)
    scrn.add_row()

if __name__ == "__main__":

    from matplotlib import pyplot

    screen = infinitephasescreen.PhaseScreenVonKarman(64, 8./32, 0.2, 40, 2)

    pyplot.ion()
    pyplot.imshow(screen.stencil)

    pyplot.figure()
    pyplot.imshow(screen.scrn)
    pyplot.colorbar()
    for i in range(100):
        screen.add_row()

        pyplot.clf()
        pyplot.imshow(screen.scrn)
        pyplot.colorbar()
        pyplot.draw()
        pyplot.pause(0.01)
