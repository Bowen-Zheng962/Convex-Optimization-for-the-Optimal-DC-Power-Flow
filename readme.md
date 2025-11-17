# Convex Optimization for the Optimal DC Power Flow

A convex optimization approach for solving the Direct Current Optimal Power Flow (DC-OPF) problem in DC grids.

## Project Overview

The DC-OPF problem aims to minimize total power loss in DC networks while satisfying power balance, voltage, and operational constraints. This project implements and benchmarks four convex optimization algorithms for solving this problem.

## Key Features

- Implementation of four convex optimization algorithms for DC-OPF
- Benchmarking against nonlinear reference solutions
- Analysis of convergence characteristics and computational performance
- Support for SOCP, SDP, SCA, and linear formulations

## Code

four convex/
  ├── build_dc10.m/    # dc system
  ├── opf_dc_lin.m/         # linear approximation algorithm
  ├── opf_dc_nonlinear.m/         # dc optimization problem fomulation
  ├── opf_dc_sca.m/         # sca algorithm
  ├── opf_dc_sdp.m/         # sdp algorithm
  ├── opf_dc_soc.m/         # soc algorithm
  └── run_all.m /               # experiment with 4 algorithms
  
  
## Quick Start

### Prerequisites
- MATLAB or Python with CVXPY
- Optimization solver (e.g., MOSEK, Gurobi, or open-source alternatives)

### Basic Usage
```matlab
% Run the main optimization script
run_all.m