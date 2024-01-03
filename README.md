# ES
C++ implementation of Evolution Strategies

## Requirements

- The Pagmo 2.19.0 C++ library (https://esa.github.io/pagmo2/). Note that the Eigen3 dependency should be installed and enabled in the Pagmo library installation ( see flag PAGMO_WITH_EIGEN3).

- Boost Library (https://www.boost.org/): it is a requirement of the Pagmo library and it is also used for unit testing.


## How to use

- Make the target ES in the CMakeLists.txt.
- Run the executable ES with the following arguments:

| parameter | description                                                                                                                                                                                    |
|-----------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| -i        | path to the file containing the input data (training data, or both training and test data)                                                                                                     |
| -m        | model name                                                                                                                                                                                     |
| -e        | ODE solver name                                                                                                                                                                                |
| -a        | algorithm name                                                                                                                                                                                 |
| -n        | maximum number of objective function evaluations                                                                                                                                               |
| -s        | seed for random number generation                                                                                                                                                              |
| -sd       | distribution of training and test sets from the input file. Example: -sd 0 34 35 49 means training set contains the data points in the interval [0, 34] and test set is the interval [34, 49]. |
| -ts       | if present, receives a file containing the data points to be used as the test set and the file passed to the -i parameter is used as the training set.                                         |

Note:the parameters -sd and -ts are mutually exclusive. Also, if none are available, uses the full training set as the test set.\
\
ODE solvers:
- lsoda:
- rk4: Runge-Kutta 4th order

Optimization algorithms:
- de: Differential Evolution
- sade: Self-Adaptive Differential Evolution
- cmaes: Covariance Matrix Adaptation Evolution Strategy
- 1+1: 1+1 Evolution Strategy
- es-i: mu+lambda Evolution Strategy with isotropic mutation
- es-ni mu+lambda Evolution Strategy with non-isotropic mutation

Example:

```
./ES -i GRN5.txt -m grn5 -e lsoda -a de -n 10000 -s 0 -sd 0 49 0 49 
```

This will run the executable ES to adjust the coefficients of the model grn5 using the training set in the file GRN5.txt and the test set in the same file. The ODE solver used is lsoda and the algorithm is differential evolution. The maximum number of objective function evaluations is 10000 and the seed for random number generation is 0. The training set contains the data points in the interval [0, 49] and the test set is the interval [0, 49].

## Output


The output of the program is a line of comma-separated values containing the following information about the best individual found by the algorithm:

- seed used
- value of the objective function (evaluated using the test set)
- execution time (in seconds)
- number of objective function evaluations
- the values of the coefficients of the model