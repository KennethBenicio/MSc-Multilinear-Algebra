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
