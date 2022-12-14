# Monte Carlo simulation of coupled, two-dimensional XY models

This repository contains the code developed for the paper G. Bighin et al. Phys. Rev. Lett. **123**, 100601 (2019). Please refer to [the paper](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.123.100601), also available on the [arXiv](https://arxiv.org/abs/1907.06253), for more details about the physics of the system and about the novel Berezinskii-Kosterlitz-Thouless paired phase we report.

# Requirements

A modern C compiler (GCC or Clang) with OpenMP support and the following libraries:

- [OpenMP](https://www.openmp.org)
- NCurses

The following libraries/codes are included in the repository:

- [libprogressbar](https://github.com/doches/progressbar) to display a nice progress bar.
- thr_pool, a simple thread pool, from the OpenSolaris documentation.

# Usage

Compile the code using the provided Makefile (`make`), eventually adapting it to point to the correct location of your C compiler and of your libraries. Run the code as `./xypcc`.
