# Homework 7: High Order Singular Value Decomposition (HOSVD)

## Problem 1 - For a random order tensor implement the HOSVD operation according the following prototype

<p align="center">
    <img src="https://latex.codecogs.com/svg.image?\LARGE&space;[\mathcal{S},&space;\boldsymbol{U}^{(1)},&space;\cdots,&space;\boldsymbol{U}^{(N)}]&space;=&space;\text{HOSVD}({\mathcal{X}}),&space;\forall\mathcal{X}&space;\in&space;\mathbb{C}^{I_{n}&space;\times&space;I_{1}&space;\times&space;\cdots&space;\times&space;I_{N}}" title="https://latex.codecogs.com/svg.image?\LARGE [\mathcal{S}, \boldsymbol{U}^{(1)}, \cdots, \boldsymbol{U}^{(N)}] = \text{HOSVD}({\mathcal{X}}), \forall\mathcal{X} \in \mathbb{C}^{I_{n} \times I_{1} \times \cdots \times I_{N}}" />
</p>

### Solution: For the sake of simplicity it was provided in this [link](https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/MATLAB/homework7.m) an example of usage of the HOSVD of a third-order tensor where its file is available inside the folders of the project. However, the implemented function does work with any order tensor just fine. Just bellow there is the NMSE analysis of the previous mentioned tensor 

| NMSE(**tenX**,**tenX_hat**)  | NMSE(**tenS**,**tenS_hat**) | NMSE(**U1**,**U1_hat**) | NMSE(**U2**,**U2_hat**) | NMSE(**U3**,**U3_hat**) |
| --------------------- | -------------------- | -------------------- | -------------------- | -------------------- |
| -611.2162 | +7.7656 | +2.6667 | +2.0000 | +1.6063 |

## Problem 2 - Considering two random tesors define then their low multilinear rank approximation according the following structure and computes the NMSE error between the origin tensor and its approximation

<p align="center">
  <img src="https://latex.codecogs.com/svg.image?\LARGE&space;\begin{align*}&space;\mathcal{X}&space;&\in&space;\mathbb{C}^{8&space;\times&space;4&space;\times&space;10}&space;\to&space;\hat{\mathcal{X}}&space;\in&space;\mathbb{C}^{R_{1}&space;\times&space;R_{2}&space;\times&space;R_{3}},&space;\\&space;\mathcal{X}&space;&\in&space;\mathbb{C}^{5&space;\times&space;5&space;\times&space;5}&space;\to&space;\hat{\mathcal{Y}}&space;\in&space;\mathbb{C}^{P_{1}&space;\times&space;P_{2}&space;\times&space;P_{3}},&space;\end{align*}" title="https://latex.codecogs.com/svg.image?\LARGE \begin{align*} \mathcal{X} &\in \mathbb{C}^{8 \times 4 \times 10} \to \hat{\mathcal{X}} \in \mathbb{C}^{R_{1} \times R_{2} \times R_{3}}, \\ \mathcal{X} &\in \mathbb{C}^{5 \times 5 \times 5} \to \hat{\mathcal{Y}} \in \mathbb{C}^{P_{1} \times P_{2} \times P_{3}}, \end{align*}" />
</p>

### Solution:

&nbsp;&nbsp;&nbsp;&nbsp; The multilinear rank approximation example can be acessed in this [link](https://github.com/KennethBenicio/MSc-Multilinear-Algebra/blob/master/MATLAB/homework7.m). Its advantage is maximized when it comes to process sparce tensors since the dimmensions can be greatly reduced without losing of too much relevant information by analysing the profile of its multiples unfoldings.
