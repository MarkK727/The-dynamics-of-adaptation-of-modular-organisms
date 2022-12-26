# Project description

A model of evolution of modular organisms with two traits.
Codes are organized into two main parts: (i) numerical simulation and analytical approximation of phenotypic evolution (MATLAB) and (ii) re-analysis of empirical avGFP fitness landscape data from https://doi.org/10.1038/nature17995 (R).

(i) MATLAB simulation codes are organized into three folders, each storing codes necessary to reproduce three models we used to get results (Constant s model, epistasis model, and Nested FGM).

Within each folder there are four items: (1) the main script that runs Wright-Fisher simulation and analytical calculations, sub-folders with (2) a set of functions necessary to run the main script ('Functions'), (3) codes that reproduce main figures ('Figures'), and (4) all the intermediate data required for reproducing our results ('Data').
To reproduce the data, put together the main script and the functions in a common working directory. Execute the main script to run the simulation and save all the intermediate variables in the same format as those saved in the 'Data' subfolder.
To generate the figures, put together the figure script in subfolder 'Figures' and the MAT-files saved in the 'Data' subfolder. Execute the figure script to reproduce the main and supplementary figures.
All parameters are set to values used to create the main figures by default. 

(ii) R codes are organized into two parts. The script 'Neighbor_Search' updates the raw avGFP data with the calculated meta information on the 1-mutant neighbors of each sampled variants of native avGFP. The code automatically saves a data frame of meta information in .csv format. The 'Figures' script download the meta information and runs all necessary codes to reproduce figures. 
