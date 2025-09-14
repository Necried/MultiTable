# Artifact for A Multi-Table Approach to Floating-Point Function Approximation (CASCON 2025 Technical Paper Track)

This artifact supports the results and claims in the CASCON 2025 paper titled "A Multi-Table Approach to Floating-Point Function Approximation" (Short Technical Paper).


## Contents
- `src/`: Contains two directories `Single/` and `Double/`, for single-precision and double-precision `recip`, `sqrt` and `rsqrt`.
- `Makefile` that executes and reports all tests for the single and double-precision functions.

## Requirements
- Dependencies:
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



