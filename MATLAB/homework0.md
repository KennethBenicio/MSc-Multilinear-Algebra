# Homework 0: Kronecker Product Properties

## Problem 1 - Evaluate the computational performance of the following Kronecker Product propertie by varying the number of columns or the number of products for randomly generated mattrices


### Solution:

&nbsp;&nbsp;&nbsp;&nbsp; In here I will briefly analyze the run time performance of the inverse operator while also using the Kronecker Product. In the first case the number of products is fixed while the number of columns is varying. In the second case we have a varying number of products for a fixed number of columns. In both cases is possible to see that is preferable to first invert the matrices before applying the Kronecker operator.  

<p align="center">
  <img src="https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/Images/hw0a1.png" />
</p>

<p align="center">
  <img src="https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/Images/hw0a2.png" />
</p>

## Problem 2 - Show algebraically that the following expression holds true

<p align="center">
  <img src="https://latex.codecogs.com/svg.image?\LARGE&space;\bg{white}\begin{align*}&space;\text{eig}(\boldsymbol{A}&space;\otimes&space;\boldsymbol{B})&space;=&space;\text{eig}(\boldsymbol{A})&space;\otimes&space;\text{eig}(\boldsymbol{B})&space;\end{align*}" title="https://latex.codecogs.com/svg.image?\LARGE \bg{white}\begin{align*} \text{eig}(\boldsymbol{A} \otimes \boldsymbol{B}) = \text{eig}(\boldsymbol{A}) \otimes \text{eig}(\boldsymbol{B}) \end{align*}" />
</p>

### Solution:

&nbsp;&nbsp;&nbsp;&nbsp; By using the eigenvalue decomposition (eig) of two matrices and apply the Kronecker Product to them it is possible to reach the intended result

<p align="center">
<img src="https://latex.codecogs.com/svg.image?\LARGE&space;\bg{white}\begin{align*}&space;\boldsymbol{A}&space;\otimes&space;\boldsymbol{B}&space;&=&space;(\boldsymbol{C}_{1}&space;\boldsymbol{\Lambda_{1}}&space;\boldsymbol{C}^{-1}_{1})&space;\otimes&space;(\boldsymbol{C}_{2}&space;\boldsymbol{\Lambda_{2}}&space;\boldsymbol{C}^{-1}_{2}),&space;\\\boldsymbol{A}&space;\otimes&space;\boldsymbol{B}&space;&=&space;(\boldsymbol{C}_{1}&space;\boldsymbol{\Lambda_{1}}&space;\otimes&space;\boldsymbol{C}_{2}&space;\boldsymbol{\Lambda_{2}})&space;(\boldsymbol{C}^{-1}_{2}&space;\otimes&space;\boldsymbol{C}^{-1}_{2}),&space;\\\boldsymbol{A}&space;\otimes&space;\boldsymbol{B}&space;&=&space;(\boldsymbol{C}_{1}&space;\otimes&space;\boldsymbol{C}_{2})&space;(\boldsymbol{\Lambda_{1}}&space;\otimes&space;\boldsymbol{\Lambda_{2}})&space;(\boldsymbol{C}^{-1}_{2}&space;\otimes&space;\boldsymbol{C}^{-1}_{2}),&space;\\\text{eig}(\boldsymbol{A}&space;\otimes&space;\boldsymbol{B})&space;&=&space;(\boldsymbol{\Lambda_{1}}&space;\otimes&space;\boldsymbol{\Lambda_{2}})&space;=&space;\text{eig}(\boldsymbol{A})&space;\otimes&space;\text{eig}(\boldsymbol{B})&space;\end{align*}" title="https://latex.codecogs.com/svg.image?\LARGE \bg{white}\begin{align*} \boldsymbol{A} \otimes \boldsymbol{B} &= (\boldsymbol{C}_{1} \boldsymbol{\Lambda_{1}} \boldsymbol{C}^{-1}_{1}) \otimes (\boldsymbol{C}_{2} \boldsymbol{\Lambda_{2}} \boldsymbol{C}^{-1}_{2}), \\\boldsymbol{A} \otimes \boldsymbol{B} &= (\boldsymbol{C}_{1} \boldsymbol{\Lambda_{1}} \otimes \boldsymbol{C}_{2} \boldsymbol{\Lambda_{2}}) (\boldsymbol{C}^{-1}_{2} \otimes \boldsymbol{C}^{-1}_{2}), \\\boldsymbol{A} \otimes \boldsymbol{B} &= (\boldsymbol{C}_{1} \otimes \boldsymbol{C}_{2}) (\boldsymbol{\Lambda_{1}} \otimes \boldsymbol{\Lambda_{2}}) (\boldsymbol{C}^{-1}_{2} \otimes \boldsymbol{C}^{-1}_{2}), \\\text{eig}(\boldsymbol{A} \otimes \boldsymbol{B}) &= (\boldsymbol{\Lambda_{1}} \otimes \boldsymbol{\Lambda_{2}}) = \text{eig}(\boldsymbol{A}) \otimes \text{eig}(\boldsymbol{B}) \end{align*}" />
</p>

&nbsp;&nbsp;&nbsp;&nbsp; The code for this results can be acessed in this [link](https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/MATLAB/homework0.m).
