# Homework 3: Least-Squares Khatri-Rao Factorization (LSKRF)

## Problem 1 - For randomly chosen matrices compute the implementation of LSKRF and compare the original matrices with the ones estimated by the algorithm. What can you conclude? 

### Solution:

| NMSE(**X**,**Xhat**)  | NMSE(**A**,**Ahat**) | NMSE(**B**,**Bhat**) |
| --------------------- | -------------------- | -------------------- |
|  -623.4093 | +11.5658 | +7.8479 |

## Problem 2 - Now assuming 1000 Monte Carlo rounds analyze the impact of the reconstruction under the presence of noisy signals for two different scenarios: (I,J) = (10,10) and (I,J) = (30,10) both with R = 4. Considering the SNR range [0, 5, 10, 15, 20, 25, 30] dB plot the curves NMSE vs. SNR for the reconstruction of the original matrix.

### Solution:

<p align="center">
  <img src="https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/Images/hw3a1.png" />
</p>

<p align="center">
  <img src="https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/Images/hw3a2.png" />
</p>

&nbsp;&nbsp;&nbsp;&nbsp; The code for this results can be acessed in this [link](https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/MATLAB/homework3.m).
