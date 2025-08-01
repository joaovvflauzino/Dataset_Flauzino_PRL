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
- Filenames ending in `-p1.txt` or `-p2.txt` indicate the datasets generated using parameter sets $p_1$ and $p_2$, respectively (as specified in the letter).

---

If you use this dataset, please cite the paper linked above.
