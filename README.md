# STIGA

**Author:** Alen Kushova  
**Year:** 2025  
**License:**  GNU GPL license (v3).

---

## Description
This project contains MATLAB code developed for scientific research purposes.  
The main goal is [...].

The code is provided for anyone who finds it useful, with proper attribution to the author.

---

## Dependencies
This project uses the following open-source libraries:

- [**GeoPDEs**](https://github.com/rafavzqz/geopdes) – research and teaching of Isogeometric Analysis, written in Octave and fully compatible with Matlab.  
- [**nurbs**](https://sourceforge.net/p/octave/nurbs/ci/default/tree/) – library for NURBS curves and surfaces

Make sure to add these libraries to the MATLAB path before running the code, i.e. add this path in the initialize script:

```matlab
addpath(genpath('path/to/geopdes'))
addpath(genpath('path/to/nurbs'))
