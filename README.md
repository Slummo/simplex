# zmax

A program to find the optimal solution and value of a linear programming problem in standard form.

## Support
The program supports:
1) Primal simplex method:
    - PhaseI to find feasible base
    - PhaseII to solve
2) Branch and bound

## How to define a model
Create a `.txt` file with these values:
1) **n**: Number of constraint
2) **m**: Number of variables
3) **c**: Vector of reduced costs (size `m`)
4) **A**: Matrix of constraints coefficients (size `n*m`)
5) **b**: RHS vector
6) **variables**: `m` values to describe the variables;
    - 0 = real
    - 1 = integer
    - 2 = binary

**Example** (bb1 model):
1) n = 3
2) m = 5
3) c = 0 -1 0 -1 0
4) A = -1 0.5 0 -0.5 0
0 2.5 -1 -0.5 0
0 1 0 1 -1
5) b = -3.5 4.5 -15
6) variables = 0 1 0 1 0

*Notes*:
- Newlines between parameters are valid
- c, A and b can be on the same line or not

## How to use

1) Compile
    ```bash
    make
    ```

2) Run
    ```bash
    make run ARGS="model_filename"
    ```
    (or optionally `make run`
    to read from `stdin`).

    **Example**: `make run ARGS="bb1"`

    *Note*: Test models are located in test_models directory, and every new model file **MUST** be located there.

## Collaborate - How to debug with gdb
1) Compile
    ```bash
    make debug
    ```

2) Run
    ```bash
    make drun ARGS="model_filename"
    ```
    (or optionally `make drun`
    to read from `stdin`).

## Dependencies
1) [GSL - GNU Scientific Library (Linear algebra)](https://www.gnu.org/software/gsl/)
2) [rc_t (Single-threaded reference-counting pointer)](https://github.com/Slummo/rc_t)