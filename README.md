# Estimating affine models parameters using closed-form likelihood approximation

* [Formulas](#formulas)
* [Parameters binding](#parameters-binding)
* [Pricing functions](#pricing-functions)
* [Approximation coefficients](#approximation-coefficients)

## Formulas
![equation](http://latex.codecogs.com/gif.latex?%5CGamma%20_%7B0%7D%5Cleft%20%28%20%5Ctheta%20%5Cright%20%29%3D%5Cbegin%7Bbmatrix%7D-A%5Cleft%20%28%20%5Ctau%20_%7B1%7D%20%5Cright%20%29%20%5C%5C%20...%20%5C%5C%20-A%5Cleft%20%28%20%5Ctau%20_%7BN%7D%20%5Cright%20%29%20%5Cend%7Bbmatrix%7D)

![equation](http://latex.codecogs.com/gif.latex?%5CGamma%20%5Cleft%20%28%20%5Ctheta%20%5Cright%20%29%3D%5Cbegin%7Bbmatrix%7DB%5Cleft%20%28%20%5Ctau%20_%7B1%7D%20%5Cright%20%29%5E%7BT%7D%20%5C%5C%20...%20%5C%5C%20B%5Cleft%20%28%20%5Ctau%20_%7BN%7D%20%5Cright%20%29%5E%7BT%7D%20%5Cend%7Bbmatrix%7D)

Extracting state vector from yields:

![equation](http://latex.codecogs.com/gif.latex?X_%7Bt%7D%3D%5Cleft%20%5B%20%5CGamma%28%5Ctheta%29%5E%7BT%7D%20%5Cright%20%5D%5E%7B-1%7D%5Cleft%20%5B%20g_%7Bt%7D-%5CGamma_%7B0%7D%28%5Ctheta%29%20%5Cright%20%5D)

Observation errors:

![equation](http://latex.codecogs.com/gif.latex?%5Cbegin%7Bbmatrix%7D%20%5Cvarepsilon%28t%2C%20t&plus;%5Ctau_%7BN&plus;1%7D%29%20%5C%5C%20...%20%5C%5C%20%5Cvarepsilon%28t%2C%20t&plus;%5Ctau_%7BN&plus;H%7D%29%20%5Cend%7Bbmatrix%7D%20%3D%20%5Cbegin%7Bbmatrix%7D%20%5Cgamma_%7B0%7D%28%5Ctau_%7BN&plus;1%7D%3B%5Ctheta%29%20%5C%5C%20...%20%5C%5C%20%5Cgamma_%7B0%7D%28%5Ctau_%7BN&plus;H%7D%3B%5Ctheta%29%20%5Cend%7Bbmatrix%7D%20&plus;%20%5Cbegin%7Bbmatrix%7D%20%5Cgamma%28%5Ctau_%7BN&plus;1%7D%3B%5Ctheta%29%5E%7BT%7D%20%5C%5C%20...%20%5C%5C%20%5Cgamma%28%5Ctau_%7BN&plus;H%7D%3B%5Ctheta%29%5E%7BT%7D%20%5Cend%7Bbmatrix%7D%20%5Cbegin%7Bbmatrix%7D%20X_%7B1t%7D%20%5C%5C%20...%20%5C%5C%20X_%7BNt%7D%20%5Cend%7Bbmatrix%7D%20-%20%5Cbegin%7Bbmatrix%7D%20g_%7B1t%7D%20%5C%5C%20...%20%5C%5C%20g_%7BNt%7D%20%5Cend%7Bbmatrix%7D)

Log-likelihood function for yields, observed without error

![equation](http://latex.codecogs.com/gif.latex?l_%7Bn%7D%28%5Ctheta%29%3D%5Cfrac%7B1%7D%7Bn%7D%5Csum_%7Bi%3D1%7D%5E%7Bn%7Dl_%7Bg%7D%28t_%7Bi%7D-t_%7Bi-1%7D%2Cg_%7Bt_%7Bi%7D%7D%7Cg_%7Bt_%7Bi-1%7D%7D%3B%5Ctheta%29)

![equation](http://latex.codecogs.com/gif.latex?l_%7Bg%7D%28%5CDelta%2Cg%7Cg_%7B0%7D%3B%5Ctheta%29%3Dln%5Cleft%20%28%20det%5Cleft%20%7C%20%5Cleft%28%20%5CGamma%28%5Ctheta%29%5E%7BT%7D%20%5Cright%20%29%5E%7B-1%7D%20%5Cright%20%7C%5Cright%20%29&plus;l_%7BX%7D%5Cleft%28%5CDelta%2Cx%7Cx_%7B0%7D%3B%5Ctheta%5Cright%29)

![equation](http://latex.codecogs.com/gif.latex?l_%7BX%7D%5E%7B%28K%29%7D%5Cleft%28%5CDelta%2Cx%7Cx_%7B0%7D%3B%5Ctheta%5Cright%29%3D-%5Cfrac%7Bm%7D%7B2%7Dln%282%5Cpi%20%5CDelta%29-%5Cfrac%7B1%7D%7B2%7Dln%20%5Cleft%20%28det%20%5Cleft%20%7C%20%5Csigma%28x%3B%5Ctheta%29%5Csigma%28x%3B%5Ctheta%29%5E%7BT%7D%5Cright%20%7C%20%5Cright%20%29%20&plus;%20%5Cfrac%7BC_%7BX%7D%5E%7B-1%7D%28x%7Cx_%7B0%7D%3B%5Ctheta%29%7D%7B%5CDelta%7D&plus;%5Csum_%7Bk%3D0%7D%5E%7BK%7DC_%7BX%7D%5E%7B%28K%29%7D%28x%7Cx_%7B0%7D%3B%5Ctheta%29%5Cfrac%7B%5CDelta%5E%7Bk%7D%7D%7Bk%21%7D)

## Parameters binding
| Code | Theory |
|------|---------------|
|θ<sub>0</sub>|b<sub>11</sub>|
|θ<sub>1</sub>|b<sub>21</sub>|
|θ<sub>2</sub>|b<sub>22</sub>|
|θ<sub>3</sub>|λ<sub>1</sub>|
|θ<sub>4</sub>|λ<sub>2</sub>|

## Pricing functions
### A<sub>0</sub>(1)
![equation](http://latex.codecogs.com/gif.latex?B%5Cleft%20%28%5Ctau%20%5Cright%20%29%3D%5Cfrac%7B%5Cdelta_%7B0%7D%5Cleft%20%28e%5E%7B%5Ctau%20b_%7B11%7D%7D%20-%201%20%5Cright%20%29%7D%7Bb_%7B11%7D%7D)

![equation](http://latex.codecogs.com/gif.latex?A%5Cleft%20%28%5Ctau%20%5Cright%20%29%3D%5Cfrac%7B%5Cdelta_%7B0%7D%5Cleft%20%28%20%5Cleft%20%28%203-4e%5E%7B%5Ctau%20b_%7B11%7D%7D&plus;e%5E%7B2%5Ctau%20b_%7B11%7D%7D&plus;2%5Ctau%20b_%7B11%7D%5Cright%20%29%5Cdelta_%7B0%7D-4b_%7B11%7D%20%5Clambda_%7B1%7D%20%5Cleft%20%281-e%5E%7B%5Ctau%20b_%7B11%7D%7D&plus;%5Ctau%20b_%7B11%7D%20%5Cright%20%29%5Cright%20%29%7D%7B4b_%7B11%7D%5E3%7D)

### A<sub>0</sub>(2)
Where δ<sub>1</sub>=0 and δ<sub>2</sub>=1


![equation](http://latex.codecogs.com/gif.latex?B%20%5Cleft%20%28%20%5Ctau%20%5Cright%20%29%3D%5Cbegin%7Bbmatrix%7D-%5Cfrac%7B%5Cleft%20%28%20e%5E%7Bb_%7B11%7D%7D%20-%20e%5E%7Bb_%7B22%7D%7D%20%5Cright%20%29%5Ctau%20b_%7B21%7D%5E%7B2%7D%7D%7Bb_%7B11%7Db_%7B22%7D%5Cleft%20%28%20b_%7B11%7D-b_%7B22%7D%20%5Cright%20%29%7D%20%5C%5C%20-%5Cfrac%7B1-%5Ctau%20e%5E%7Bb_%7B22%7D%7D%7D%7Bb_%7B22%7D%7D%20%5Cend%7Bbmatrix%7D)

![equation](http://latex.codecogs.com/gif.latex?A%5Cleft%20%28%20%5Ctau%20%5Cright%20%29%3D%20%5Cbegin%7Bmatrix%7D%20%5Cfrac%7B%5Ctau%7D%7Bb_%7B11%7D%5E%7B2%7Db_%7B22%7D%5E%7B2%7D%28%20b_%7B11%7D%5E%7B2%7D-2b_%7B11%7Db_%7B22%7D&plus;b_%7B22%7D%5E%7B2%7D%20%29%7D%28%28%20%5Cfrac%7B1%7D%7B6%7D%20e%5E%7B2b_%7B11%7D%7D&plus;%5Cfrac%7B1%7D%7B6%7D%20e%5E%7B2b_%7B22%7D%7D-%5Cfrac%7B1%7D%7B3%7D%20e%5E%7Bb_%7B11%7D&plus;b_%7B22%7D%7D%29%5Ctau%5E%7B2%7D%20b_%7B21%7D%5E%7B4%7D&plus;%28%5Cfrac%7B1%7D%7B2%7De%5E%7Bb_%7B11%7D%7D-%5Cfrac%7B1%7D%7B2%7De%5E%7Bb_%7B22%7D%7D%29%5Ctau%20%5Clambda_%7B1%7D%20b_%7B11%7Db_%7B21%7D%5E%7B2%7Db_%7B22%7D%5E%7B2%7D&plus;%20%5C%5C&plus;b_%7B11%7D%5E%7B3%7Db_%7B22%7D%28%20%5Ctau%20e%5E%7Bb_%7B22%7D%7D-%5Cfrac%7B1%7D%7B3%7D%20%5Ctau%5E%7B2%7D%20e%5E%7B2b_%7B22%7D%7D&plus;b_%7B22%7D%5Clambda_%7B2%7D%20%28%202-%5Ctau%20e%5E%7Bb_%7B22%7D%7D%29-1%29%20&plus;%20%5C%5C&plus;b_%7B11%7D%5E%7B4%7D%28%20%5Cfrac%7B1%7D%7B2%7D-%20%5Cfrac%7B1%7D%7B2%7D%20%5Ctau%20e%5E%7Bb_%7B22%7D%7D&plus;%5Cfrac%7B1%7D%7B6%7D%5Ctau%5E%7B2%7D%20e%5E%7B2b_%7B22%7D%7D%20&plus;%20b_%7B22%7D%20%5Clambda_%7B2%7D%28%5Cfrac%7B1%7D%7B2%7D%20%5Ctau%20e%5E%7Bb_%7B22%7D%7D%20-1%29%20&plus;%20%5C%5C%20&plus;%20b_%7B11%7D%5E%7B2%7Db_%7B22%7D%28%28%5Cfrac%7B1%7D%7B2%7D-%5Cfrac%7B1%7D%7B2%7D%5Ctau%20e%5E%7Bb_%7B22%7D%7D%20&plus;%20%5Cfrac%7B1%7D%7B6%7D%20%5Ctau%5E%7B2%7D%20e%5E%7B2b_%7B22%7D%7D%29b_%7B22%7D%20&plus;%28%5Cfrac%7B1%7D%7B2%7De%5E%7Bb_%7B22%7D%7D-%5Cfrac%7B1%7D%7B2%7De%5E%7Bb_%7B11%7D%7D%29%5Ctau%20b_%7B21%7D%5E%7B2%7D%20%5Clambda_%7B1%7D%20&plus;%20%28%5Cfrac%7B1%7D%7B2%7D%20%5Ctau%20e%5E%7Bb_%7B22%7D%7D%20-%201%29%20b_%7B22%7D%5E%7B2%7D%20%5Clambda_%7B2%7D%29%20%29%20%5Cend%7Bmatrix%7D)

## Approximation coefficients
### A<sub>0</sub>(1)

![equation](http://latex.codecogs.com/gif.latex?p_%7Bx%7D%5E%7B%282%29%7D%5Cleft%20%28%5CDelta%2Cx%7Cx_%7B0%7D%3B%5Ctheta%20%5Cright%20%29%3Dp_%7Bx%7D%5E%7B%280%29%7D%5Cleft%20%28%5CDelta%2Cx%7Cx_%7B0%7D%3B%5Ctheta%20%5Cright%20%29%5Cleft%20%5B1&plus;c_%7B1%7D%28x%7Cx_%7B0%7D%3B%5Ctheta%29%5CDelta&plus;c_%7B2%7D%28x%7Cx_%7B0%7D%3B%5Ctheta%29%5Cfrac%7B%5CDelta%5E%7B2%7D%7D%7B2%7D%20%5Cright%20%5D)

![equation](http://latex.codecogs.com/gif.latex?p_%7Bx%7D%5E%7B%280%29%7D%5Cleft%20%28%5CDelta%2Cx%7Cx_%7B0%7D%3B%5Ctheta%20%5Cright%20%29%3D%5Cfrac%7B1%7D%7B%5Csqrt%7B2%5Cpi%5CDelta%7D%7Dexp%5Cleft%20%5B-%5Cfrac%7B%28x-x_%7B0%7D%29%5E%7B2%7D%7D%7B2%5CDelta%7D&plus;%5Cfrac%7Bb_%7B11%7Dx%5E%7B2%7D%7D%7B2%7D-%5Cfrac%7Bb_%7B11%7Dx_%7B0%7D%5E%7B2%7D%7D%7B2%7D%20%5Cright%20%5D)

![equation](http://latex.codecogs.com/gif.latex?c_%7B1%7D%28x%7Cx_%7B0%7D%3B%5Ctheta%29%3D%5Cfrac%7Bb_%7B11%7D%7D%7B6%7D%5Cleft%20%283&plus;b_%7B11%7Dx%5E2&plus;b_%7B11%7Dxx_%7B0%7D&plus;b_%7B11%7Dx_%7B0%7D%5E%7B2%7D%20%5Cright%20%29)

![equation](http://latex.codecogs.com/gif.latex?c_%7B2%7D%28x%7Cx_%7B0%7D%3B%5Ctheta%29%3D%5Cfrac%7Bb_%7B11%7D%5E%7B2%7D%7D%7B36%7D%5Cleft%20%283&plus;6b_%7B11%7Dx%5E2&plus;b_%7B11%7D%5E%7B2%7Dx%5E%7B4%7D&plus;2b_%7B11%7Dxx_%7B0%7D%5Cleft%20%283&plus;b_%7B11%7Dx%5E%7B2%7D%20%5Cright%20%29&plus;3b_%7B11%7Dx_%7B0%7D%5E%7B2%7D%5Cleft%20%282&plus;b_%7B11%7Dx%5E%7B2%7D%20%5Cright%20%29-2b_%7B11%7D%5E%7B2%7Dxx_%7B0%7D%5E%7B3%7D&plus;b_%7B11%7D%5E%7B2%7Dx_%7B0%7D%5E%7B4%7D%20%5Cright%20%29)
