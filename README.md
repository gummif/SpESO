SpESO
============

Spectrum Estimation by Sparse Optimization (SpESO). MATLAB code for energy spectrum estimation of turbulence by compressive sampling (compressed sensing). The core of the method is Quasi-Oracle Multilevel Orthogonal Matching Pursuit (QOMOMP), which is a multilevel variation of OMP using a priori information about the signal.


About
----------------------

Written by Gudmundur F. Adalsteinsson (2013-2014). This code is intended as a supplement to a journal article and a thesis to promote reproducible research in computational science.

Copyright (c) 2014 Gudmundur Adalsteinsson (MIT License). See file LICENSE for license and warranty details.

Installation and requirements
----------------------

To use, it is simplest to change the current directory to the folder SpESO. The main files can then be executed directly (a path command temporarily adds the functions folder to the current path). The execution time might be about an hour for SpESO. Once complete, the plotting scripts can be used to visualize the saved results (which are in .mat files).

Requires a wavelet package (MATLAB's Wavelet Toolbox, WaveLab) but can also be used with the supplied CDF 9/7 transform. Package type is set with the DWTTYPE parameter.


Reproducible Research
----------------------

The classical branches of science are deductive and empirical science. In the last few decades, a third branch has emerged: computational science.

"Scientific Computation is emerging as absolutely central to the scientific method. Unfortunately,
it is error-prone and currently immature: traditional scientific publication is
incapable of finding and rooting out errors in scientific computation; this must be recognized
as a crisis. Reproducible computational research, in which the full computational environment
that produces a result is published along with the article, is an important recent
development, and a necessary response to this crisis." [15 Years of Reproducible Research in
Computational Harmonic Analysis](http://statweb.stanford.edu/~donoho/Reports/2008/15YrsReproResch-20080426.pdf)
