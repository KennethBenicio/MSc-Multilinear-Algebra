# Homework 4: Least Squares Kronecker Product Factorization (LSKronF)

## Problem 1 - For randomly chosen matrices compute the implementation of LSKronF and compare the original matrices with the ones estimated by the algorithm. What can you conclude? 

### Solution:

| NMSE(**X**,**Xhat**)  | NMSE(**A**,**Ahat**) | NMSE(**B**,**Bhat**) |
| --------------------- | -------------------- | -------------------- |
|  -619.2196 | +13.5472 | +9.5922 |

## Problem 2 - Now assuming 1000 Monte Carlo rounds analyze the impact of the reconstruction under the presence of noisy signals for two different scenarios: (I,J,P,Q) = (2,4,3,5) and (I,J,P,Q) = (4,8,3,5). Considering the SNR range [0, 5, 10, 15, 20, 25, 30] dB plot the curves NMSE vs. SNR for the reconstruction of the original matrix.

### Solution:

<p align="center">
  <img src="https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/Images/hw4a1.png" />
</p>

<p align="center">
  <img src="https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/Images/hw4a2.png" />
</p>

&nbsp;&nbsp;&nbsp;&nbsp; The code for this results can be acessed in this [link](https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/MATLAB/homework4.m).
