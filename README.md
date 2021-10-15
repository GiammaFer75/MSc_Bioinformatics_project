# MSc Bioinformatics with Systems Biology
This repository contains the Python module Proteogenome and a folder with proteomics data. 
These data are provided in order to test Proteogenome.


## Run the simulation
In order to run the simulation, create a woring folder.

Dowload the in the folder Proteogenome.py and the contend of TestFiles folder. 

## 1. Prepare input data
**1.1 Upload proteomics data**

    # Define the set of columns names in the csv input file, that match the required data for Proteogenome.

	csv_cols = ['protein.Accession','peptide.seq','peptide.modification',
            'peptide.MatchedProducts','peptide.MatchedProductsSumInten']

    proteome = pg.load_proteomic_data('proteome_test.csv',
                                   target_cols=csv_cols,reshape=True)

**2	Apply the PTMs to the peptide sequences**

    remove_ PTMs = ['None','Carbamidomethyl']

    proteome_PTM = pg.PoGo_input_df(proteome, PTMs_to_remove = remove_ PTMs)

**1.3	Generate the PoGo input file**

    pg.PoGo_input_file('PoGo_input_file.txt',proteome_PTM,'Experiment_test')


## 2. Generate protein heat map
**2.1 Create protein reference tables for level of expression and peptides**

    # Protein expression levels dictionary
    prot_expre = pg.get_proteins_expr(proteomePTM)

    # Protein-Peptides dictionary
    prot_peptides = pg.get_proteins(proteomePTM)

**2.2 Fetching protein genomic coordinates**

    # Set the output filenames for the raw protein records and the mismatched codes.
    p_1_filename  = 'protein_M1.txt'
    p_M_filename  = 'p_notfound.txt'




#### Folder TestFiles:
- **PoGo_peptides.bed**      :  This file contains the peptide map performed by PoGo
- **PoGo_peptides_PTM.bed**  :  This file contains the map of peptide PTMs performed by PoGo
- **proteome_test.csv**      :  This file contains a PLGS table with peptides data from MS/MS analysis