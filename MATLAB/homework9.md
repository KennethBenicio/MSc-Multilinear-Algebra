# Homework 9: Multidimensional Least-Squares Khatri-Rao Factorization (MLS-KRF)

## Problem 1 - For randomly chosen matrices compute the implementation of MLSKRF for a random number of matrices solving the following problem and compare the original matrices with the ones estimated by the algorithm. What can you conclude? 

<p align="center">
  <img src="https://latex.codecogs.com/svg.image?\inline&space;\LARGE&space;&space;\left(\hat{\boldsymbol{A}}^{(1)},&space;\cdots,&space;\hat{\boldsymbol{A}}^{(N)}\right)&space;=&space;\underset{\boldsymbol{A}^{(1)},&space;\cdots,&space;\boldsymbol{A}^{(N)}}{\text{min}}&space;\left|\left|&space;\boldsymbol{X}&space;-&space;\boldsymbol{A}^{(1)}&space;\diamond&space;\cdots&space;\diamond&space;\boldsymbol{A}^{(N)}&space;\right|\right|^2_{\text{F}}" title="https://latex.codecogs.com/svg.image?\inline \LARGE \left(\hat{\boldsymbol{A}}^{(1)}, \cdots, \hat{\boldsymbol{A}}^{(N)}\right) = \underset{\boldsymbol{A}^{(1)}, \cdots, \boldsymbol{A}^{(N)}}{\text{min}} \left|\left| \boldsymbol{X} - \boldsymbol{A}^{(1)} \diamond \cdots \diamond \boldsymbol{A}^{(N)} \right|\right|^2_{\text{F}}" />
</p>

### Solution:

| NMSE(**X**,**Xhat**)  | NMSE(**A1**,**A1hat**) | NMSE(**A2**,**A2hat**) | NMSE(**A3**,**A3hat**) |
| --------------------- | -------------------- | -------------------- | -------------------- |
|  -606.2255 | +3.1432 | +5.0165 | +4.7297 |

## Problem 2 - Now assuming 1000 Monte Carlo rounds analyze the impact of the reconstruction using the MLS-KRF under the presence of noisy signals in the SNR range of [0, 5, 10, 15, 20, 25, 30] dB for the scenario defined by the matrices of dimmensions: I1 = 2, I2 = 3, I3 = 4 and R = 4. Plot the curve for NMSE vs. SNR considering the reconstructed matrix.

### Solution:

<p align="center">
  <img src="https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/Images/hw9a1.png" />
</p>

&nbsp;&nbsp;&nbsp;&nbsp; The code for this results can be acessed in this [link](https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/MATLAB/homework9.m).
