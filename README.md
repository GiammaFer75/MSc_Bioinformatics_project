# MSc Bioinformatics with Systems Biology
This repository contains the Python module Proteogenome and a folder with proteomics data. 
These data are provided in order to test Proteogenome.


## Run the simulation
In order to run the simulation, create a woring folder.

Dowload the in the folder Proteogenome.py and the contend of TestFiles folder. 

## 1. Prepare input data
**1.1 Upload proteomics data**

    # Define the set of columns names in the csv input file, that match the      # required data for Proteogenome
	csv_cols = ['protein.Accession','peptide.seq','peptide.modification',
            'peptide.MatchedProducts','peptide.MatchedProductsSumInten']

    proteome = pg.load_proteomic_data('proteome_test.csv',
                                   target_cols=csv_cols,reshape=True)

#### Folder TestFiles:
- **PoGo_peptides.bed**      :  This file contains the peptide map performed by PoGo
- **PoGo_peptides_PTM.bed**  :  This file contains the map of peptide PTMs performed by PoGo
- **proteome_test.csv**      :  This file contains a PLGS table with peptides data from MS/MS analysis