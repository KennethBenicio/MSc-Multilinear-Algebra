# Homework 2: Khatri-Rao Product Properties

## Problem 1 - Evaluate the performance of the following methods to compute the pseudoinverse of a given Khatri-Rao product.

- 1st Method:

<p align="center">
<img src="https://latex.codecogs.com/svg.image?\LARGE&space;\bg{white}\begin{align*}&space;(\boldsymbol{A}&space;\diamond&space;\boldsymbol{B})^{\dagger}&space;=&space;\text{pinv}(\boldsymbol{A}&space;\diamond&space;\boldsymbol{B})&space;\end{align*}" title="https://latex.codecogs.com/svg.image?\LARGE \bg{white}\begin{align*} (\boldsymbol{A} \diamond \boldsymbol{B})^{\dagger} = \text{pinv}(\boldsymbol{A} \diamond \boldsymbol{B}) \end{align*}" />
</p>

- 2nd Method:

<p align="center">
<img src="https://latex.codecogs.com/svg.image?\LARGE&space;\bg{white}\begin{align*}&space;(\boldsymbol{A}&space;\diamond&space;\boldsymbol{B})^{\dagger}&space;=&space;[(\boldsymbol{A}&space;\diamond&space;\boldsymbol{B})^{\text{T}}&space;(\boldsymbol{A}&space;\diamond&space;\boldsymbol{B})]^{-1}&space;(\boldsymbol{A}&space;\diamond&space;\boldsymbol{B})^{\text{T}}&space;\end{align*}" title="https://latex.codecogs.com/svg.image?\LARGE \bg{white}\begin{align*} (\boldsymbol{A} \diamond \boldsymbol{B})^{\dagger} = [(\boldsymbol{A} \diamond \boldsymbol{B})^{\text{T}} (\boldsymbol{A} \diamond \boldsymbol{B})]^{-1} (\boldsymbol{A} \diamond \boldsymbol{B})^{\text{T}} \end{align*}" />
</p>

- 3rd Method:

<p align="center">
<img src="https://latex.codecogs.com/svg.image?\LARGE&space;\bg{white}\begin{align*}&space;(\boldsymbol{A}&space;\diamond&space;\boldsymbol{B})^{\dagger}&space;=&space;[(\boldsymbol{A}^{\text{T}}\boldsymbol{A})(\boldsymbol{B}^{\text{T}}\boldsymbol{B})]^{-1}&space;(\boldsymbol{A}&space;\diamond&space;\boldsymbol{B})^{\text{T}}&space;\end{align*}" title="https://latex.codecogs.com/svg.image?\LARGE \bg{white}\begin{align*} (\boldsymbol{A} \diamond \boldsymbol{B})^{\dagger} = [(\boldsymbol{A}^{\text{T}}\boldsymbol{A})(\boldsymbol{B}^{\text{T}}\boldsymbol{B})]^{-1} (\boldsymbol{A} \diamond \boldsymbol{B})^{\text{T}} \end{align*}" />
</p>

### Solution:

<p align="center">
  <img src="https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/Images/hw2a1.png" />
</p>

<p align="center">
  <img src="https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/Images/hw2a2.png" />
</p>

## Problem 2 - Evaluate the performance of the Khatri-Rao as the number of products increases.

### Solution:

<p align="center">
  <img src="https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/Images/hw2a3.png" />
</p>

&nbsp;&nbsp;&nbsp;&nbsp; The code for this results can be acessed in this [link](https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/MATLAB/homework2.m).
