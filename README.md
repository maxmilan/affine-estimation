# Estimating affine models parameters using closed-form likelihood approximation

## Formulas
![equation](http://latex.codecogs.com/gif.latex?%5CGamma%20_%7B0%7D%5Cleft%20%28%20%5Ctheta%20%5Cright%20%29%3D%5Cbegin%7Bbmatrix%7D-A%5Cleft%20%28%20%5Ctau%20_%7B1%7D%20%5Cright%20%29%20%5C%5C%20...%20%5C%5C%20-A%5Cleft%20%28%20%5Ctau%20_%7BN%7D%20%5Cright%20%29%20%5Cend%7Bbmatrix%7D)

![equation](http://latex.codecogs.com/gif.latex?%5CGamma%20%5Cleft%20%28%20%5Ctheta%20%5Cright%20%29%3D%5Cbegin%7Bbmatrix%7DB%5Cleft%20%28%20%5Ctau%20_%7B1%7D%20%5Cright%20%29%5E%7BT%7D%20%5C%5C%20...%20%5C%5C%20B%5Cleft%20%28%20%5Ctau%20_%7BN%7D%20%5Cright%20%29%5E%7BT%7D%20%5Cend%7Bbmatrix%7D)


## Parameters binding
| Code | Theory |
|------|---------------|
|θ<sub>0</sub>|b<sub>11</sub>|
|θ<sub>1</sub>|b<sub>21</sub>|
|θ<sub>2</sub>|b<sub>22</sub>|
|θ<sub>3</sub>|λ<sub>1</sub>|
|θ<sub>4</sub>|λ<sub>2</sub>|

## Pricing functions
Where δ<sub>1</sub>=0 and δ<sub>2</sub>=1
### A<sub>0</sub>(2)
![equation](http://latex.codecogs.com/gif.latex?B%20%5Cleft%20%28%20%5Ctau%20%5Cright%20%29%3D%5Cbegin%7Bbmatrix%7D-%5Cfrac%7B%5Cleft%20%28%20e%5E%7Bb_%7B11%7D%7D%20-%20e%5E%7Bb_%7B22%7D%7D%20%5Cright%20%29%5Ctau%20b_%7B21%7D%5E%7B2%7D%7D%7Bb_%7B11%7Db_%7B22%7D%5Cleft%20%28%20b_%7B11%7D-b_%7B22%7D%20%5Cright%20%29%7D%20%5C%5C%20-%5Cfrac%7B1-%5Ctau%20e%5E%7Bb_%7B22%7D%7D%7D%7Bb_%7B22%7D%7D%20%5Cend%7Bbmatrix%7D)

![equation](http://latex.codecogs.com/gif.latex?A%5Cleft%20%28%20%5Ctau%20%5Cright%20%29%3D%20%5Cbegin%7Bmatrix%7D%20%5Cfrac%7B%5Ctau%7D%7Bb_%7B11%7D%5E%7B2%7Db_%7B22%7D%5E%7B2%7D%28%20b_%7B11%7D%5E%7B2%7D-2b_%7B11%7Db_%7B22%7D&plus;b_%7B22%7D%5E%7B2%7D%20%29%7D%28%28%20%5Cfrac%7B1%7D%7B6%7D%20e%5E%7B2b_%7B11%7D%7D&plus;%5Cfrac%7B1%7D%7B6%7D%20e%5E%7B2b_%7B22%7D%7D-%5Cfrac%7B1%7D%7B3%7D%20e%5E%7Bb_%7B11%7D&plus;b_%7B22%7D%7D%29%5Ctau%5E%7B2%7D%20b_%7B21%7D%5E%7B4%7D&plus;%28%5Cfrac%7B1%7D%7B2%7De%5E%7Bb_%7B11%7D%7D-%5Cfrac%7B1%7D%7B2%7De%5E%7Bb_%7B22%7D%7D%29%5Ctau%20%5Clambda_%7B1%7D%20b_%7B11%7Db_%7B21%7D%5E%7B2%7Db_%7B22%7D%5E%7B2%7D&plus;%20%5C%5C&plus;b_%7B11%7D%5E%7B3%7Db_%7B22%7D%28%20%5Ctau%20e%5E%7Bb_%7B22%7D%7D-%5Cfrac%7B1%7D%7B3%7D%20%5Ctau%5E%7B2%7D%20e%5E%7B2b_%7B22%7D%7D&plus;b_%7B22%7D%5Clambda_%7B2%7D%20%28%202-%5Ctau%20e%5E%7Bb_%7B22%7D%7D%29-1%29%20&plus;%20%5C%5C&plus;b_%7B11%7D%5E%7B4%7D%28%20%5Cfrac%7B1%7D%7B2%7D-%20%5Cfrac%7B1%7D%7B2%7D%20%5Ctau%20e%5E%7Bb_%7B22%7D%7D&plus;%5Cfrac%7B1%7D%7B6%7D%5Ctau%5E%7B2%7D%20e%5E%7B2b_%7B22%7D%7D%20&plus;%20b_%7B22%7D%20%5Clambda_%7B2%7D%28%5Cfrac%7B1%7D%7B2%7D%20%5Ctau%20e%5E%7Bb_%7B22%7D%7D%20-1%29%20&plus;%20%5C%5C%20&plus;%20b_%7B11%7D%5E%7B2%7Db_%7B22%7D%28%28%5Cfrac%7B1%7D%7B2%7D-%5Cfrac%7B1%7D%7B2%7D%5Ctau%20e%5E%7Bb_%7B22%7D%7D%20&plus;%20%5Cfrac%7B1%7D%7B6%7D%20%5Ctau%5E%7B2%7D%20e%5E%7B2b_%7B22%7D%7D%29b_%7B22%7D%20&plus;%28%5Cfrac%7B1%7D%7B2%7De%5E%7Bb_%7B22%7D%7D-%5Cfrac%7B1%7D%7B2%7De%5E%7Bb_%7B11%7D%7D%29%5Ctau%20b_%7B21%7D%5E%7B2%7D%20%5Clambda_%7B1%7D%20&plus;%20%28%5Cfrac%7B1%7D%7B2%7D%20%5Ctau%20e%5E%7Bb_%7B22%7D%7D%20-%201%29%20b_%7B22%7D%5E%7B2%7D%20%5Clambda_%7B2%7D%29%20%29%20%5Cend%7Bmatrix%7D)
