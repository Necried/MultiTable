# Artifact for A Multi-Table Approach to Floating-Point Function Approximation (CASCON 2025 Technical Paper Track)

This artifact supports the results and claims in the CASCON 2025 paper titled "A Multi-Table Approach to Floating-Point Function Approximation" (Short Technical Paper).


## Contents
- `src/`: Contains two directories `Single/` and `Double/`, for single-precision and double-precision `recip`, `sqrt` and `rsqrt`.
- `Makefile` that executes and reports all tests for the single and double-precision functions.

## Requirements
- Clang compiler, Version 10 or later.
- Make

## Installation
1. Clone or download this repository.
2. Ensure that Clang and Make is on your system.

## Execution

To run the full test suite, just run the Makefile in the parent directory:

```bash
make
```

This will execute Makefiles in `src/Single/*` and `src/Double/*` directories.
A log directory will be generated in the parent directory `./log`, with associated log files
for each precision and function.

During execution, a summary of each function will be printed to the terminal. For example,
single-precision reciprocal will generate the following:

```bash
===================================================
RECIP SINGLE PRECISION LOG REPORT:
Results for recip single precision recorded in:
        log/recip_sp.log
Ulp error interval:
         [-1.00000, 3.00000]
===================================================
```

## Verification

In the `log` directory, open the associated log file for a given function and precision.

#### Single

There are two parts in the log output:
1. Detailed output for each interval, including values that produced the highest ulp error.
2. Summary histogram outlining ulp errors.

#### Double

The log file has the same structure as single-precision, but does not have detailed output as we performed random testing instead of exhaustive testing.

#### Histogram

An example of the histogram output taken from `recip_sp.log`:

```
Histogram of exponent 0
Number of test 8388608
===neg to pos===     ===count===  =====absolute===== 
[-inf,-10)   0.00%           0
  [-10,-5)   0.00%           0
   [-5,-2)   0.00%           0
   [-2,-1)   0.00%           0
 [-1,-0.5)   0.00%           3
  (-0.5,0]  37.76%     3167424
   (0,0.5]   0.00%           0     [0,0.5]  37.76%
   (0.5,1]  51.02%     4280014     (0.5,1]  51.02%
     (1,2]  11.22%      941099       (1,2]  11.22%
     (2,5]   0.00%          68       (2,5]   0.00%
    (5,10]   0.00%           0      (5,10]   0.00%
  (10,inf]   0.00%           0    (10,inf]   0.00%
```

From these, we see that the highest ulp errors are in `(2, 5]`. 
The detailed output above in the log file states that this has maximum ulp error of 3,
which lies in this interval.