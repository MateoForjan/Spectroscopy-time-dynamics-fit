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


----------------------------------------------------------------------------------------------
Copyright (c) 2023 Mateo Forjan

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
