# Project description

A model of evolution of modular organisms with two traits. 
Codes are organized into two main parts: (i) numerical simulation and analytical approximation of phenotypic evolution (MATLAB) and (ii) re-analysis of empirical avGFP fitness landscape data from https://doi.org/10.1038/nature17995 (R).

(i) MATLAB simulation codes are organized into three folders, each storing codes necessary to reproduce three models we used to get results (Constant s model, epistasis model, and Nested FGM).

Within each folder there are four items: (1) the main script that runs Wright-Fisher simulation and analytical calculations, sub-folders with (2) a set of functions necessary to run the main script ('Functions'), (3) codes that reproduce main figures ('Figures'), and (4) all the intermediate data required for reproducing our results ('Data'). 
To reproduce the results, put all items in the sub-folders in a common directory. 
All parameters are set to values used to create the main figures by default. 

(ii) R codes are organized into two major parts. 
