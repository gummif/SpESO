SpESO
============

Spectrum Estimation by Sparse Optimization (SpESO). MATLAB code for energy spectrum estimation of turbulence by compressive sampling (compressed sensing). The core of the method is Quasi-Oracle Multilevel Orthogonal Matching Pursuit (QOMOMP), which is a multilevel variation of OMP using a priori information about the signal.

Requires a wavelet package (MATLAB's Wavelet Toolbox, WaveLab) but can also be used with the supplied CDF 9/7 transform. Package type is set with the DWTTYPE parameter.


About
----------------------

Written by Gudmundur F. Adalsteinsson (2013-2014). This code is intended as a supplement to a journal article and a thesis to promote reproducible research in computational science, the third branch of science, the other two being deductive and empirical science.


Reproducible Research
----------------------

"Scientific Computation is emerging as absolutely central to the scientific method. Unfortunately,
it is error-prone and currently immature: traditional scientific publication is
incapable of finding and rooting out errors in scientific computation; this must be recognized
as a crisis. Reproducible computational research, in which the full computational environment
that produces a result is published along with the article, is an important recent
development, and a necessary response to this crisis." [15 Years of Reproducible Research in
Computational Harmonic Analysis](http://statweb.stanford.edu/~donoho/Reports/2008/15YrsReproResch-20080426.pdf)
