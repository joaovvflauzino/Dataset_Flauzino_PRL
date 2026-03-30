# Dataset for "Quantifying Disorder in Data" — Physical Review Letters

This repository contains the datasets and code related to the models used in the following publication:

**João Vitor Vieira Flauzino, Thiago Lima Prado, Norbert Marwan, Jürgen Kurths, and Sergio Roberto Lopes.**  
*Quantifying Disorder in Data*.  
*Physical Review Letters*. DOI: [https://doi.org/10.1103/1y98-x33s](https://doi.org/10.1103/1y98-x33s)

## Repository Structure

### `1_Generating_Dataset_Main-Text/`
Contains Julia scripts that generate the time series for the six models studied in the **main text** of the letter.

### `2_Dataset_End-Matter/`
Contains `.txt` files with time series data for the models used in the **End Matter** of the letter, for comparative analysis between methods.

- Each subfolder corresponds to a specific model.
- Filenames ending in `-p1-X.txt` or `-p2-X.txt` indicate the datasets generated using parameter sets $p_1$ and $p_2$, respectively (as specified in the letter), where `X` is the index of a sample in a total of 50 time series per parameter set.

### `3_Algorithm_for_Quantifying_Disorder`

- **`1_QuantifyingDisorderInData.jl`**  
This is an implementation of the proposed method, written in Julia. The script is ready to use — simply replace the time series in the code to compute the corresponding degree of disorder. Note that for correct execution, it requires the file `Labels_disorder_N4.txt` to be in the same folder.

- **`2_ComputingClassesInBinaryMatrices.jl`**  
  This script implements the procedure described in the manuscript to compute the equiprobability classes of binary matrices.

- **`Labels_disorder_N4.txt`**  
  This file contains the class index assigned to each binary matrix, ordered according to its row position.
---

If you use this algorithm or dataset, please cite the original paper where this method was introduced (linked above).
