# Spectroscopy-time-dynamics-fit

![github3](https://user-images.githubusercontent.com/92934177/236469577-a7172534-9bcc-41a2-a4c2-8b6710f5044e.png)

Tool used for fitting functions on time dynamics obtained by measurements.

Singular value decomposition (SVD) analysis is automatically done at start.
Singular values and vectors (spectra and time dynamics) can be explored by choosing the number of SVD values.
SVD analysis is of great help when considering doing global analysis.

Function for fitting: convolution of exponential function (decay or rising dynamics) with a gaussian (representing instrument response function).
Additionally, a gaussian and first two derivatives are summed with mentioned convolution for the description of coherent artifact near time zero.

Residuals of fitting and subtracted coherent artifact can be explored.

Details are given in MANUAL.
  -Choosing the 'wavelength to fit' plots a time-cut at chosen wavelength.
  -Changing 'Nr. of exponentials' and initial coefficients plots the initial guess for a fit.
  -'Fit Curve' fits the curve and in pop-up windows gives the fit results and goodness of fit.
