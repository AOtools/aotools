from distutils.core import setup
import versioneer


with open("README.rst", "r") as readme:
    long_description = readme.read()

setup(
    name='aotools',
    author_email='andrewpaulreeves@gmail.com, matthew.townson@durham.ac.uk',
    url='https://github.com/aotools/aotools',
    packages=['aotools',
              'aotools.astronomy',
              'aotools.functions',
              'aotools.image_processing',
              'aotools.turbulence',
              'aotools.wfs',
              ],
    description='A set of useful functions for Adaptive Optics in Python',
    long_description=long_description,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    install_requires=[
        'numpy',
        'scipy',
        'astropy',
        'astroquery',
    ],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
)
