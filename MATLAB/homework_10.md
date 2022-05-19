# Homework 10: Multidimensional Least-Squares Kronecker Factorization  (MLS-KronF)

## Problem 1 - For randomly chosen matrices compute the implementation of MLS-KronF for a random number of matrices solving the following problem and compare the original matrices with the ones estimated by the algorithm. What can you conclude? 

<p align="center">
  <img src="https://latex.codecogs.com/svg.image?\inline&space;\LARGE&space;&space;\left(\hat{\boldsymbol{A}}^{(1)},&space;\cdots,&space;\hat{\boldsymbol{A}}^{(N)}\right)&space;=&space;\underset{\boldsymbol{A}^{(1)},&space;\cdots,&space;\boldsymbol{A}^{(N)}}{\text{min}}&space;\left|\left|&space;\boldsymbol{X}&space;-&space;\boldsymbol{A}^{(1)}&space;\otimes&space;\cdots&space;\otimes&space;\boldsymbol{A}^{(N)}&space;\right|\right|^2_{\text{F}}" title="https://latex.codecogs.com/svg.image?\inline \LARGE \left(\hat{\boldsymbol{A}}^{(1)}, \cdots, \hat{\boldsymbol{A}}^{(N)}\right) = \underset{\boldsymbol{A}^{(1)}, \cdots, \boldsymbol{A}^{(N)}}{\text{min}} \left|\left| \boldsymbol{X} - \boldsymbol{A}^{(1)} \otimes \cdots \otimes \boldsymbol{A}^{(N)} \right|\right|^2_{\text{F}}" />
</p>

### Solution:

| NMSE(**X**,**Xhat**)  | NMSE(**A1**,**A1hat**) | NMSE(**A2**,**A2hat**) | NMSE(**A3**,**A3hat**) |
| --------------------- | -------------------- | -------------------- | -------------------- |
|  -605.1941 | +11.9214 | +11.5548 | +6.0950 |

## Problem 2 - Now assuming 1000 Monte Carlo rounds analyze the impact of the reconstruction using the MLS-KronF under the presence of noisy signals in the SNR range of [0, 5, 10, 15, 20, 25, 30] dB for the scenario defined by the matrices of dimmensions: I1 = J1 = 2, I2 = J2 = 3 and I3 = J3 = 4. Plot the curve for NMSE vs. SNR considering the reconstructed matrix.

### Solution:

<p align="center">
  <img src="https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/Images/hw10a1.png" />
</p>

&nbsp;&nbsp;&nbsp;&nbsp; The code for this results can be acessed in this [link](https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/MATLAB/homework10.m).
