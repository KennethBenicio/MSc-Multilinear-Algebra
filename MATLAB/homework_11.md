# Homework 11: Alternating Least Squares (ALS) Algorithm

## Problem 1 - Implement the ALS algorithm for the third-order tensor provided in the files. Estimate the factor matrices by solving the following problem and compare the estimations with the original ones. Explain the results.

<p align="center">
  <img src="https://latex.codecogs.com/svg.image?\inline&space;\LARGE&space;\left(\hat{\boldsymbol{A}},&space;\hat{\boldsymbol{B}},&space;\hat{\boldsymbol{C}}\right)&space;=&space;\underset{\boldsymbol{A},&space;\boldsymbol{B},&space;\boldsymbol{C}}{\text{min}}&space;\left|&space;\left|&space;\mathcal{X}&space;-&space;\sum^{R}_{r&space;=&space;1}&space;\boldsymbol{a}_{r}&space;\circ&space;\boldsymbol{b}_{r}&space;\circ&space;\boldsymbol{c}_{r}&space;\right|&space;\right|^{2}_{\text{F}}" title="https://latex.codecogs.com/svg.image?\inline \LARGE \left(\hat{\boldsymbol{A}}, \hat{\boldsymbol{B}}, \hat{\boldsymbol{C}}\right) = \underset{\boldsymbol{A}, \boldsymbol{B}, \boldsymbol{C}}{\text{min}} \left| \left| \mathcal{X} - \sum^{R}_{r = 1} \boldsymbol{a}_{r} \circ \boldsymbol{b}_{r} \circ \boldsymbol{c}_{r} \right| \right|^{2}_{\text{F}}" />
</p>

### Solution:

The error of the i-th iteration is obtained with 

<p align="center">
  <img src="https://latex.codecogs.com/svg.image?\inline&space;\LARGE&space;e_{i}&space;=&space;\left|&space;\left|&space;[\mathcal{X}]_{(1)}&space;-&space;\hat{\boldsymbol{A}}&space;(&space;\hat{\boldsymbol{C}}&space;\diamond&space;\hat{\boldsymbol{B}}&space;)^{\text{T}}&space;\right|\right|_{\text{F}}" title="https://latex.codecogs.com/svg.image?\inline \LARGE e_{i} = \left| \left| [\mathcal{X}]_{(1)} - \hat{\boldsymbol{A}} ( \hat{\boldsymbol{C}} \diamond \hat{\boldsymbol{B}} )^{\text{T}} \right|\right|_{\text{F}}" />
</p>

<p align="center">
  <img src="https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/Images/hw11a1.png" />
</p>

## Problem 2 - Now assuming 1000 Monte Carlo rounds analyze the impact of the reconstruction using the ALS under the presence of noisy signals in the SNR range of [0, 5, 10, 15, 20, 25, 30] dB for the scenario defined by the matrices of dimmensions: (I, J, K, R) = (10, 4, 2, 3). Plot the curve for NMSE vs. SNR considering the reconstructed tensor and its factor matrices.

### Solution:

<p align="center">
  <img src="https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/Images/hw11a2.png" />
</p>

&nbsp;&nbsp;&nbsp;&nbsp; The code for this results can be acessed in this [link](https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/MATLAB/homework11.m).
