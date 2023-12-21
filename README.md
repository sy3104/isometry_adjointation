# isometry_adjointation

This repository contains

- SDP code to obtain the optimal success probability or approximation error of isometry adjointation, isometry inversion, and universal error detection.

This code is accompanied to the following paper:

- Satoshi Yoshida, Akihito Soeda, and Mio Murao, In preparation.

## Requirement

The SDP code is written in Matlab and requires the following interpreter:

- [CVX](http://cvxr.com): a Matlab-based convex modeling framework

These codes also use functions of QETLAB ([QETLAB](https://qetlab.com): A MATLAB Toolbox for Quantum Entanglement), but all used functions are contained in the subfolder [QETLAB_used_functions](https://github.com/sy3104/isometry_adjointation/tree/main/QETLAB_used_functions).

It has been tested on MATLAB R2021b and CVX 2.2.

## Description

The main script of this repository is

- [run_optimal_isometry_adjointation.m](https://github.com/sy3104/isometry_adjointation/blob/main/run_optimal_isometry_adjointation.m): Code to obtain the optimal approximation error or success probability of transforming n calls of any isometry operations with input dimension din and output dimension dout into its inverse map, corresponding POVM, or adjoint map.

The script [run_optimal_untiary_inversion.m](https://github.com/sy3104/deterministic_exact_unitary_inversion/blob/main/run_optimal_untiary_inversion.m) makes use of the functions in the subfolder [QETLAB_used_functions](https://github.com/sy3104/deterministic_exact_unitary_inversion/tree/main/QETLAB_used_functions) and the following Matlab codes in the subfolder [sdp](https://github.com/sy3104/isometry_adjointation/tree/main/sdp):

- [parallel_isometry_adjointation](https://github.com/sy3104/isometry_adjointation/blob/main/sdp/parallel_isometry_adjointation.m): Primal SDP code to obtain the optimal approximation error or success probability of transforming parallel n calls of any isometry operations with input dimension din and output dimension dout into its inverse map, corresponding POVM, or adjoint map.
- [parallel_isometry_adjointation_dual](https://github.com/sy3104/isometry_adjointation/blob/main/sdp/parallel_isometry_adjointation_dual.m): Dual SDP code to obtain the optimal approximation error or success probability of transforming parallel n calls of any isometry operations with input dimension din and output dimension dout into its inverse map, corresponding POVM, or adjoint map.
- [sequential_isometry_adjointation](https://github.com/sy3104/isometry_adjointation/blob/main/sdp/sequential_isometry_adjointation.m): Primal SDP code to obtain the optimal approximation error or success probability of transforming sequential n calls of any isometry operations with input dimension din and output dimension dout into its inverse map, corresponding POVM, or adjoint map.
- [sequential_isometry_adjointation_dual](https://github.com/sy3104/isometry_adjointation/blob/main/sdp/sequential_isometry_adjointation_dual.m): Dual SDP code to obtain the optimal approximation error or success probability of transforming sequential n calls of any isometry operations with input dimension din and output dimension dout into its inverse map, corresponding POVM, or adjoint map.
- [indefinite_isometry_adjointation_dual](https://github.com/sy3104/isometry_adjointation/blob/main/sdp/indefinite_isometry_adjointation_dual.m): Dual SDP code to obtain the optimal approximation error or success probability of transforming n calls of any isometry operations with input dimension din and output dimension dout into its inverse map, corresponding POVM, or adjoint map using general protocols including indefinite causal order.

These Matlab codes utilize group-theoretic calculations done by a [Sagemath](https://www.sagemath.org) code shown in [young_diagrams.sage](https://github.com/sy3104/isometry_adjointation/blob/main/sdp/young_diagrams.sage), which is taken from [deterministic_exact_unitry_inversion](https://github.com/sy3104/deterministic_exact_unitary_inversion).


## License

This code is under the [MIT license](https://opensource.org/licenses/MIT).
