# Multigrid Method Library for 3D Poisson Solver in Fortran

## Developed by Dr. Ji-Hoon Kang

This library is the result of the advanced 3D multigrid code developed by Dr. Ji-Hoon Kang at KISTI. 
You can contact Dr. Kang for more information at jhkang@kisti.re.kr.

## Modifications by Dr. Mingyu Yang

Further modifications and enhancements were made by Dr. Mingyu Yang from Yonsei University. 
For inquiries related to these modifications, please reach out to Dr. Yang at yang926@yonsei.ac.kr.

## Program Description

The main procedure is in `poisson.f90` in `src` folder. Please check the algorithm procedure start from `poisson.f90`.

## Installation Instructions

1. Please use nvfortran (NVIDIA HPC SDK).
2. Use command at the home folder as follows:
    ```
	make all
    ```

## How to Run

1. Go to the `run` folder and use command as follows:
    ```
	mpirun -np 8 ./a.out ./PARA_INPUT.inp
    ```
2. The core number and npx*npy*npz on `PARA_INPUT.inp` should be identical!

## Uninstallation

1. Use command at the home folder as follows:
    ```
	make clean
    ```
