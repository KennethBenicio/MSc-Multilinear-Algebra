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

# Homework 1: Hadamard, Kronecker and Khatri-Rao Products

## Problem 1 - Create algorithms to compute the Hadamard, Kronecker and Khatri-Rao products and compare the run time performance with the native functions available on MATLAB. 

### Solution:

<p align="center">
  <img src="https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/Images/hw1a1.png" />
</p>

<p align="center">
  <img src="https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/Images/hw1a2.png" />
</p>

<p align="center">
  <img src="https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/Images/hw1a3.png" />
</p>

&nbsp;&nbsp;&nbsp;&nbsp; The code for this results can be acessed in this [link](https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/MATLAB/homework1.m).

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

# Homework 5: Kronecker Product Single Value Decomposition (KPSVD)

## Problem 1 - Implement the generalization of the LSKronF operation, the KPSVD, according the following structure

<p align="center">
<img src="https://latex.codecogs.com/svg.image?\LARGE&space;\boldsymbol{X}&space;=&space;\sum^{r_{kp}}_{k&space;=&space;1}&space;\sigma_{k}&space;\boldsymbol{U}_{k}&space;\otimes&space;\boldsymbol{V}_{k}" title="https://latex.codecogs.com/svg.image?\LARGE \boldsymbol{X} = \sum^{r_{kp}}_{k = 1} \sigma_{k} \boldsymbol{U}_{k} \otimes \boldsymbol{V}_{k}" />
</p>

### Solution: 

&nbsp;&nbsp;&nbsp;&nbsp; The implementation of this operation can be acessed in this [link](https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/MATLAB/homework5.m).

## Problem 2 - Considering a random matrix compute its nearest rank-r approximation using the KPSVD prototype.

### Solution: 

&nbsp;&nbsp;&nbsp;&nbsp; The implementation of this operation can be acessed in this [link](https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/MATLAB/homework5.m). In the example is considered that the original matrix has a full rank equals to 9, then two approximations using the KPSVD are provided: One using the full-rank approximation and other using a r-rank approximation (r < 9). The NMSE between some approximations and its original matrix are provided in the sequence

| Full-rank  | 7-rank | 5-rank | 3-rank | 1-rank |
| --------------------- | -------------------- | -------------------- | -------------------- | -------------------- |
|  -604.8023 |  -51.3722 | -26.7589 | -16.9054 | -6.3068 |

# Homework 6: Unfolding, Folding, and n-mode Product

## Problem 1 - For a random order tensor implement the operation unfolding according to the following prototype 

<p align="center">
  <img src="https://latex.codecogs.com/svg.image?\LARGE&space;[\mathcal{X}]_{n}&space;=&space;\text{unfold}(\mathcal{X},&space;[I_{1}&space;\cdots&space;I_{N}&space;],&space;n)&space;\in&space;\mathbb{C}^{I_{n}&space;\times&space;I_{1}&space;\cdots&space;I_{n-1}&space;I_{n&plus;1}&space;\cdots&space;I_{N}}" title="https://latex.codecogs.com/svg.image?\LARGE [\mathcal{X}]_{n} = \text{unfold}(\mathcal{X}, [I_{1} \cdots I_{N} ], n) \in \mathbb{C}^{I_{n} \times I_{1} \cdots I_{n-1} I_{n+1} \cdots I_{N}}" />
</p>

### Solution:

&nbsp;&nbsp;&nbsp;&nbsp; The generic order unfolding was implemented but just for the sake of simplicity the next examples will be using a third-order tensor. Thus, I choose to define the following tensor and provide its file inside the folders of the project.

<p align="center">
<img src="https://latex.codecogs.com/svg.image?\LARGE&space;\bg{white}\begin{align*}\boldsymbol{X}_{1}&space;=&space;\begin{bmatrix}&space;1&space;&&space;4&space;&&space;7&space;&&space;10&space;\\&space;2&space;&&space;5&space;&&space;8&space;&&space;11&space;\\&space;3&space;&&space;6&space;&&space;9&space;&&space;12&space;\\\end{bmatrix},\boldsymbol{X}_{2}&space;=&space;\begin{bmatrix}&space;13&space;&&space;16&space;&&space;19&space;&&space;22&space;\\&space;14&space;&&space;17&space;&&space;20&space;&&space;23&space;\\&space;15&space;&&space;18&space;&&space;21&space;&&space;24&space;\\\end{bmatrix}&space;\\\end{align*}" title="https://latex.codecogs.com/svg.image?\LARGE \bg{white}\begin{align*}\boldsymbol{X}_{1} = \begin{bmatrix} 1 & 4 & 7 & 10 \\ 2 & 5 & 8 & 11 \\ 3 & 6 & 9 & 12 \\\end{bmatrix},\boldsymbol{X}_{2} = \begin{bmatrix} 13 & 16 & 19 & 22 \\ 14 & 17 & 20 & 23 \\ 15 & 18 & 21 & 24 \\\end{bmatrix} \\\end{align*}" />
</p>

&nbsp;&nbsp;&nbsp;&nbsp; From the tensor just above it is possible to obtain the following results for its unfoldings by calling the method unfolding(Tensor, Mode) inside the tensor class provided in the link [link](https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/MATLAB/homework6.m).

<p align="center">
<img src="https://latex.codecogs.com/svg.image?\LARGE&space;\bg{white}[\boldsymbol{X}]_{(1)}&space;=\begin{bmatrix}&space;1&space;&&space;4&space;&&space;7&space;&&space;10&space;&&space;13&space;&&space;16&space;&&space;19&space;&&space;22&space;\\&space;2&space;&&space;5&space;&&space;8&space;&&space;11&space;&&space;14&space;&&space;17&space;&&space;20&space;&&space;23&space;\\&space;3&space;&&space;6&space;&&space;9&space;&&space;12&space;&&space;15&space;&&space;18&space;&&space;21&space;&&space;24&space;\\\end{bmatrix}" title="https://latex.codecogs.com/svg.image?\LARGE \bg{white}[\boldsymbol{X}]_{(1)} =\begin{bmatrix} 1 & 4 & 7 & 10 & 13 & 16 & 19 & 22 \\ 2 & 5 & 8 & 11 & 14 & 17 & 20 & 23 \\ 3 & 6 & 9 & 12 & 15 & 18 & 21 & 24 \\\end{bmatrix}" />
</p>

<p align="center">
<img src="https://latex.codecogs.com/svg.image?\LARGE&space;\bg{white}[\boldsymbol{X}]_{(2)}&space;=\begin{bmatrix}&space;1&space;&space;&&space;2&space;&&space;3&space;&&space;13&space;&&space;14&space;&&space;15&space;\\&space;4&space;&space;&&space;5&space;&&space;6&space;&&space;16&space;&&space;17&space;&&space;18&space;\\&space;7&space;&space;&&space;8&space;&&space;9&space;&&space;19&space;&&space;20&space;&&space;21&space;\\&space;10&space;&&space;11&space;&&space;12&space;&&space;22&space;&&space;23&space;&&space;24&space;\\\end{bmatrix}" title="https://latex.codecogs.com/svg.image?\LARGE \bg{white}[\boldsymbol{X}]_{(2)} =\begin{bmatrix} 1 & 2 & 3 & 13 & 14 & 15 \\ 4 & 5 & 6 & 16 & 17 & 18 \\ 7 & 8 & 9 & 19 & 20 & 21 \\ 10 & 11 & 12 & 22 & 23 & 24 \\\end{bmatrix}" />
</p>

<p align="center">
<img src="https://latex.codecogs.com/svg.image?\LARGE&space;\bg{white}[\boldsymbol{X}]_{(3)}&space;=\begin{bmatrix}1&space;&&space;2&space;&&space;3&space;&&space;4&space;&&space;5&space;&&space;6&space;&&space;7&space;&&space;8&space;&&space;9&space;&&space;10&space;&&space;11&space;&&space;12&space;\\13&space;&&space;14&space;&&space;15&space;&&space;16&space;&&space;17&space;&&space;18&space;&&space;19&space;&&space;20&space;&&space;21&space;&&space;22&space;&&space;23&space;&&space;24&space;\\\end{bmatrix}" title="https://latex.codecogs.com/svg.image?\LARGE \bg{white}[\boldsymbol{X}]_{(3)} =\begin{bmatrix}1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & 10 & 11 & 12 \\13 & 14 & 15 & 16 & 17 & 18 & 19 & 20 & 21 & 22 & 23 & 24 \\\end{bmatrix}" />
</p>

## Problem 2 - In similar fashion implement the folding operation that will convert back the result from the unfolding operation to its original tensor according to the following prototype

<p align="center">
  <img src="https://latex.codecogs.com/svg.image?\LARGE&space;\mathcal{X}&space;=&space;\text{fold}([\mathcal{X}]_{n},&space;[I_{1}&space;\cdots&space;I_{N}&space;],&space;n)&space;\in&space;\mathbb{C}^{I_{1}&space;\times&space;\cdots&space;\times&space;I_{N}}" title="https://latex.codecogs.com/svg.image?\LARGE \mathcal{X} = \text{fold}([\mathcal{X}]_{n}, [I_{1} \cdots I_{N} ], n) \in \mathbb{C}^{I_{1} \times \cdots \times I_{N}}" />
</p>

### Solution:

&nbsp;&nbsp;&nbsp;&nbsp; In a similar fashion it is possible to reobtain the original tensor from each one of its unfolding by calling the method folding(Unfolding, Dimensions, Mode) inside the tensor class. From this will be possible to confirm that the tensor is correctly reconstruct. Inside this [link](https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/MATLAB/homework6.m) will be possible to find an example of using both the unfolding and folding methods.

## Problem 3 - For a given number of matrices implement the operation of multilinear tensor product according the following prototype

<p align="center">
    <img src="https://latex.codecogs.com/svg.image?\LARGE&space;\mathcal{Y}&space;=&space;\mathcal{X}&space;\times_{1}&space;\boldsymbol{U}_{1}&space;\times_{2}&space;\cdots&space;\times_{N}&space;\boldsymbol{U}_{N}" title="https://latex.codecogs.com/svg.image?\LARGE \mathcal{Y} = \mathcal{X} \times_{1} \boldsymbol{U}_{1} \times_{2} \cdots \times_{N} \boldsymbol{U}_{N}" />
</p>

### Solution:

&nbsp;&nbsp;&nbsp;&nbsp; Once more it will be possible to find an example in this [link](https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/MATLAB/homework6.m) about the usage of the method n_mod_prod(Tensor, Matrices, Modes) from the class tensor. Besides there is also a validation file available inside the folders of the project.

# Homework 7: High Order Singular Value Decomposition (HOSVD)

## Problem 1 - For a random order tensor implement the HOSVD operation according the following prototype

<p align="center">
    <img src="https://latex.codecogs.com/svg.image?\LARGE&space;[\mathcal{S},&space;\boldsymbol{U}^{(1)},&space;\cdots,&space;\boldsymbol{U}^{(N)}]&space;=&space;\text{HOSVD}({\mathcal{X}}),&space;\forall\mathcal{X}&space;\in&space;\mathbb{C}^{I_{n}&space;\times&space;I_{1}&space;\times&space;\cdots&space;\times&space;I_{N}}" title="https://latex.codecogs.com/svg.image?\LARGE [\mathcal{S}, \boldsymbol{U}^{(1)}, \cdots, \boldsymbol{U}^{(N)}] = \text{HOSVD}({\mathcal{X}}), \forall\mathcal{X} \in \mathbb{C}^{I_{n} \times I_{1} \times \cdots \times I_{N}}" />
</p>

### Solution: For the sake of simplicity it was provided in this [link](https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/MATLAB/homework7.m) an example of usage of the HOSVD of a third-order tensor where its file is available inside the folders of the project. However, the implemented function does work with any order tensor just fine. Just bellow there is the NMSE analysis of the previous mentioned tensor 

| NMSE(**tenX**,**tenX_hat**)  | NMSE(**tenS**,**tenS_hat**) | NMSE(**U1**,**U1_hat**) | NMSE(**U2**,**U2_hat**) | NMSE(**U3**,**U3_hat**) |
| --------------------- | -------------------- | -------------------- | -------------------- | -------------------- |
| -609.1300 | +2.0933 | +2.6551 | +1.9105 | +1.9316 |

## Problem 2 - Considering two random tesors define then their low multilinear rank approximation according the following structure and computes the NMSE error between the origin tensor and its approximation

<p align="center">
  <img src="https://latex.codecogs.com/svg.image?\LARGE&space;\begin{align*}&space;\mathcal{X}&space;&\in&space;\mathbb{C}^{8&space;\times&space;4&space;\times&space;10}&space;\to&space;\hat{\mathcal{X}}&space;\in&space;\mathbb{C}^{R_{1}&space;\times&space;R_{2}&space;\times&space;R_{3}},&space;\\&space;\mathcal{X}&space;&\in&space;\mathbb{C}^{5&space;\times&space;5&space;\times&space;5}&space;\to&space;\hat{\mathcal{Y}}&space;\in&space;\mathbb{C}^{P_{1}&space;\times&space;P_{2}&space;\times&space;P_{3}},&space;\end{align*}" title="https://latex.codecogs.com/svg.image?\LARGE \begin{align*} \mathcal{X} &\in \mathbb{C}^{8 \times 4 \times 10} \to \hat{\mathcal{X}} \in \mathbb{C}^{R_{1} \times R_{2} \times R_{3}}, \\ \mathcal{X} &\in \mathbb{C}^{5 \times 5 \times 5} \to \hat{\mathcal{Y}} \in \mathbb{C}^{P_{1} \times P_{2} \times P_{3}}, \end{align*}" />
</p>

### Solution:

&nbsp;&nbsp;&nbsp;&nbsp; The multilinear rank approximation example can be acessed in this [link](https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/MATLAB/homework7.m). Its advantage is maximized when it comes to process sparce tensors since the dimmensions can be greatly reduced without losing of too much relevant information by analysing the profile of its multiples unfoldings.

# Homework 8: Higher-Order Orthogonal Iteration (HOOI)

## Problem 1 -

### Solution:

## Problem 2 -

### Solution:

# Homework 9: Multidimensional Least-Squares Khatri-Rao Factorization (MLS-KRF)

## Problem 1 - For randomly chosen matrices compute the implementation of MLSKRF for a random number of matrices solving the following problem and compare the original matrices with the ones estimated by the algorithm. What can you conclude? 

<p align="center">
  <img src="https://latex.codecogs.com/svg.image?\inline&space;\large&space;(\hat{\boldsymbol{A}}^{(1)},&space;\cdots,&space;\hat{\boldsymbol{A}}^{(N)})&space;=&space;\underset{\boldsymbol{A}^{(1)},&space;\cdots,&space;\boldsymbol{A}^{(N)}}{\text{min}}&space;||&space;\boldsymbol{X}&space;-&space;\boldsymbol{A}^{(1)}&space;\diamond&space;\cdots&space;\diamond&space;\boldsymbol{A}^{(N)}&space;||^2_{\text{F}}" title="https://latex.codecogs.com/svg.image?\inline \large (\hat{\boldsymbol{A}}^{(1)}, \cdots, \hat{\boldsymbol{A}}^{(N)}) = \underset{\boldsymbol{A}^{(1)}, \cdots, \boldsymbol{A}^{(N)}}{\text{min}} || \boldsymbol{X} - \boldsymbol{A}^{(1)} \diamond \cdots \diamond \boldsymbol{A}^{(N)} ||^2_{\text{F}}" />
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

# Homework 10: Multidimensional Least-Squares Kronecker Factorization  (MLS-KronF)

## Problem 1 - For randomly chosen matrices compute the implementation of MLS-KronF for a random number of matrices solving the following problem and compare the original matrices with the ones estimated by the algorithm. What can you conclude? 

<p align="center">
  <img src="https://latex.codecogs.com/svg.image?\inline&space;\large&space;(\hat{\boldsymbol{A}}^{(1)},&space;\cdots,&space;\hat{\boldsymbol{A}}^{(N)})&space;=&space;\underset{\boldsymbol{A}^{(1)},&space;\cdots,&space;\boldsymbol{A}^{(N)}}{\text{min}}&space;||&space;\boldsymbol{X}&space;-&space;\boldsymbol{A}^{(1)}&space;\otimes&space;\cdots&space;\otimes&space;\boldsymbol{A}^{(N)}&space;||^2_{\text{F}}" title="https://latex.codecogs.com/svg.image?\inline \large (\hat{\boldsymbol{A}}^{(1)}, \cdots, \hat{\boldsymbol{A}}^{(N)}) = \underset{\boldsymbol{A}^{(1)}, \cdots, \boldsymbol{A}^{(N)}}{\text{min}} || \boldsymbol{X} - \boldsymbol{A}^{(1)} \otimes \cdots \otimes \boldsymbol{A}^{(N)} ||^2_{\text{F}}" />
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

# Homework 11: Alternating Least Squares (ALS) Algorithm

## Problem 1 -

### Solution:

<p align="center">
  <img src="https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/Images/hw11a1.png" />
</p>

## Problem 2 -

### Solution:

<p align="center">
  <img src="https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/Images/hw11a2.png" />
</p>

&nbsp;&nbsp;&nbsp;&nbsp; The code for this results can be acessed in this [link](https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/MATLAB/homework11.m).
