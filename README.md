# STIGA
Space-Time IsoGeonetric Analysis of evolutionary PDEs.

## Description
This project contains MATLAB code developed mainly for scientific research purposes.  

The main goal is [...].

---

## Dependencies
This project uses the following open-source libraries:

- [**GeoPDEs**](https://github.com/rafavzqz/geopdes) – research and teaching of Isogeometric Analysis, written in Octave and fully compatible with Matlab.  
- [**nurbs**](https://sourceforge.net/p/octave/nurbs/ci/default/tree/) – library for NURBS curves and surfaces

Make sure to add these libraries to the MATLAB path before running the code, i.e. add this path in the initialize script:

```matlab
addpath(genpath('path/to/geopdes'))
addpath(genpath('path/to/nurbs'))
```

---

## License 
Because the code depends on GPL-licensed libraries, the entire project is distributed under the GPLv3 terms.

You are free to use, study, modify, and redistribute this software, provided that any distributed version (modified or not) is also released under the GPLv3. 

Terms and conditions are reported in the LICENSE file but can also be found here [**GNU GENERAL PUBLIC LICENSE**](https://www.gnu.org/licenses/gpl-3.0-standalone.html).
