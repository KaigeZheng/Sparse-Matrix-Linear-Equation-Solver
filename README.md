# Sparse-Matrix-Linear-Equation-Solver
A C-based solver for sparse linear systems of equations Ax=b, utilizing iterative methods such as Jacobi or Gauss-Seidel. Designed for efficient handling of large, sparse matrices, this solver aims to provide fast and memory-efficient solutions.

# Build

## build and run a demo

```shell
make clean && make
./test
```

## test the compute_det function

```shell
make clean && make compute_det
./compute_det
```

## test the solver function

```shell
make clean && make solver
./solver
```

## run the baseline

```shell
make clean && make baseline
./baseline
```

## run the optimized-CSR version

```shell
make clean && make best
./CSR-best > log
```

> Note: It's essential to put sparse matrix (`.mtx`) file in ./data directory.
