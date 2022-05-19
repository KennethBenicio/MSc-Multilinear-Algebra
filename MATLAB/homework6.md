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
