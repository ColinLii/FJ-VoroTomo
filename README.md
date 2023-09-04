# Example Codes for the FJ-VoroTomo

> Autor: Zhengbo Li et al.,
> Department of Earth and Space Sciences, SUSTech, China

>  Related Paper: Multiple Voronoi partition improves multimodal dispersion imaging from ambient noise: a case study of Lasso dense array (https://doi.org/10.1029/2022JB026081)

## Introduction
This is an example code of the implementation of the FJ-VoroTomo method for 3-D Vs imaging of a seismic array. The step-by-step workflow is shown below:

1. Before using this code, you need to prepare the calculated Noise Cross-correlation Functions (NCFs) first. The CC-FJpy package (https://github.com/ColinLii/CC-FJpy) can be used for the NCFs' calculation.
2. For our codes, the NCFs in frequency domain were saved as a HDF5 file, which contains "StationPairs", "f", "ncfs" and "r" four matrices, which were calculated by the CC-FJpy. 
3. After the NCFs preparation, use the "VoroTomo.py" to generate the VoroTomo partitions and calculate the F-J spectrum for each partition. 
Then, you will get all the output HDF5 files containing F-J spectrum in the out dir. Detailed parameters please see the comments in the "VoroTomo.py".
4. You can use DisperNet (https://github.com/Dongsh/DisperNet) to obtain the dispersion curves of each F-J spectrum. Saving these dispersion curves in the dir "curves" as ASCII files with the same file name but with different suffixes of the F-J spectrum file (e.g., output/vor5_2.h5 -> curves/vor5_2.txt). The dispersion curve file contains three columns: frequency, phase velocity, and mode number.
5. Using the "generate_curves.py" to generate Dispersion curves for each target point from the previous dispersion curves.
6. Finally, use the DisbaTomo (https://github.com/pan3rock/DisbaTomo) to invert the Vs structure of each target point.

## Python Package Used

- Python3
- Numpy
- Scipy
- os
- h5py
- Pandas
- CC-FJpy
- DisbaTomo
- DisperNet

