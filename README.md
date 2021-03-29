This code is a port of the [LSST](https://github.com/lsst) codebase that deals
with finding and masking cosmic ray events on images. This is a python/C++ based project
which is built with the python setuptools and uses cmake in the background.

## Contents

- [Requirements](#requirements)
- [Installing](#installing)

## Requirements
The following packages are required in order to compile the C++ code:

- [Boost](https://www.boost.org/) (version 1.36.0 or newer)
- OpenMP (any flavor that uses the OpenMP standard API should work)
- [Gnu Scientific Library](https://www.gnu.org/software/gsl/) 
- [Cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/)
- [Minuit2](https://root.cern.ch/doc/master/Minuit2Page.html)
- [AST](https://github.com/lsst/starlink_ast) (from the LSST git fork)
- [FFTW](http://www.fftw.org/)
- [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page)
- Python (version 3.7 or newer)
- [pybind11](https://github.com/pybind/pybind11)
- [ndarray](https://github.com/ndarray/ndarray)

The following Python packages are required to run the python side of the code:

- numpy
- astropy
- setuptools

## Installing

- Clone this repository
- Run the install script

```bash
git clone https://github.com/DarkEnergySurvey/cosmicRays.git
cd cosmicRays
python setup.py build
python setup.py install
```
You can specify an install location with
```bash
python setup.py install --prefix=<path to install>
```
