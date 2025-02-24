# Winkler-et-al.-2025
We supply the code necessary to recreate the results from the paper.
The codes are distributed into two folders corresponding to the respectively used software.

### Numerical Continuation
The codes placed in this folder are written for the use of the software [Auto-07p](https://github.com/auto-07p/auto-07p) which is freely available.
The code provides results for steadily moving cells (c.f. Fig 2) and implements time independent ordinary differential equations only.

### Fininte Volume Simulations
To simulate the full systems of PDEs we use the [fipy package](https://www.ctcms.nist.gov/fipy/) for python.
For different physcial limits we provide two analogous impementations for the initial system as well as the rigid limit.

### Environment Setup
For all codes we use a respective conda environment. In order for all files to be executed properly, please install the respective environments through the provided `.yaml` files.
