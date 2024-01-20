# TOEPLITZ INERTIA


![GitHub License](https://img.shields.io/github/license/GTAC-ITEAM-UPV/toeplitz_intertia)
![GitHub language count](https://img.shields.io/github/languages/count/GTAC-ITEAM-UPV/toeplitz_intertia)
![GitHub top language](https://img.shields.io/github/languages/top/GTAC-ITEAM-UPV/toeplitz_intertia)

![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/GTAC-ITEAM-UPV/toeplitz_intertia)
![GitHub pull requests](https://img.shields.io/github/issues-pr-raw/GTAC-ITEAM-UPV/toeplitz_intertia)
![GitHub issues](https://img.shields.io/github/issues/GTAC-ITEAM-UPV/toeplitz_intertia)

![GitHub watchers](https://img.shields.io/github/watchers/GTAC-ITEAM-UPV/toeplitz_intertia)
![GitHub forks](https://img.shields.io/github/forks/GTAC-ITEAM-UPV/toeplitz_intertia)
![GitHub Repo stars](https://img.shields.io/github/stars/GTAC-ITEAM-UPV/toeplitz_intertia)


### Software for computing the inertia of a Hermitian Toeplitz-block matrix.

The code in this zip file is based on algorithm 2.3 in the paper 
by M. Ng and W. Trench, [“Numerical solution of the eigenvalue problem for hermitian toeplitz matrices”](https://doi.org/10.1137/0610010), 
SIAM Journal on Matrix Analysis and Applications, vol. 10, 08 1997




It includes two mex files, [`inercia_toep_mex.c`](./inercia_toep_mex.c) and [`inercia_toep_mex2.c`](./inercia_toep_mex2.c). These two files must be compiled requiring an appropriate  C compiler.
Once done, you can examine and run the example file [`example_toep_block3x3.m`](./example_toep_block3x3.m).

### References
[1] Ng, Michael & Trench, William. (1997). _Numerical Solution of the Eigenvalue Problem for Hermitian Toeplitz Matrices_. SIAM Journal on Matrix Analysis and Applications. Vol. 10. DOI: [10.1137/0610010](https://doi.org/10.1137/0610010)

[2]  Ferrer Contreras, Miguel; García Mollá, Víctor Manuel; Vidal Maciá, Antonio Manuel; de Diego Antón, María ; Gonzalez, Alberto. _Assessment of stability of distributed FxLMS active noise control systems_.
Signal Processing (ISSN 0165-1684). Elsevier. Vol. 210, pp: 1-13, 2023. DOI: [10.1016/j.sigpro.2023.109087](https://doi.org/10.1016/j.sigpro.2023.109087)

