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

**1.2	Apply the PTMs to the peptide sequences**

    remove_ PTMs = ['None','Carbamidomethyl']

    proteome_PTM = pg.PoGo_input_df(proteome, PTMs_to_remove = remove_ PTMs)

**1.3	Generate the PoGo input file**

    pg.PoGo_input_file('PoGo_input_file.txt',proteome_PTM,'Experiment_test')


## 2. Generate protein heat map
**METHOD 1 - 1. Create protein reference tables for level of expression and peptides**

    # Protein expression levels dictionary
    prot_expre = pg.get_proteins_expr(proteomePTM)

    # Protein-Peptides dictionary
    prot_peptides = pg.get_proteins(proteomePTM)

**METHOD 1 - 2. Fetching protein genomic coordinates**

    # Set the output filenames for the raw protein records and the mismatched codes.
    p_1_filename  = 'protein_M1.txt'
    p_M_filename  = 'p_notfound.txt'
**METHOD 1 - 3. Download the protein records**
     
    prot_raw_records,mismatch_protein = pg.download_proteins_set(prot_expre,prot_peptides, 
                         db_type='UPKB', print_out=True,
                         proteins_output       = p_1_filename, 
                         error_proteins_output = p_M_filename)

Now it is possible to try to recover the protein codes that did not found a match in the previous step. There are two METHODs that could be applied alternatively. The task include a conversion of UniProt codes in Ensembl codes. Then the UniProtKB database will be quired again with the Ensembl codes. The funtion used for this conversion is 'UP2Ens_code_conv'. Setting the parameter db_type='IDmap' will be used the API service for Retrieve / ID map for the conversion. Instead, setting the same parameter to 'UniParc', will be used the 'UniParc' API.

**METHOD 2/3 - 1. Create protein reference tables for level of expression and peptides**
     
    M_expr,M_pept = pg.generate_prot_expr_pept_dict(mismatch_protein)


**METHOD 2 - *Convert the UniProt codes in Ensembl codes - Retrieve / ID map service***

    M_expr,M_pept = pg.UP2Ens_code_conv(M_expr, M_pept, db_type='IDmap', print_out=True)

**METHOD 3 - *Convert the UniProt codes in Ensembl codes - UniParc service***

	M_expr,M_pept = pg.UP2Ens_code_conv(M_expr, M_pept, db_type='UniParc', print_out=True)


**METHOD 2/3 - 2. Download the protein records.**

	p_recovered_filename = 'p_recovered.txt'

	rec_raw_prot = pg.download_proteins_set(M_expr, M_pept, db_type='Ens', print_out=True, proteins_output=p_recovered_filename)

**2.1 	Create the protein databases for the methods that have been used**

	# Generate the protein database for the protein records fetched trough the Method 1
	prot_tab_M1, exon_tab_M1 = pg.prot_genome_loc(prot_raw_records, ens=False)

	# Generate the protein database for the protein records fetched trough the Method 2 or Method 3
	M_prot_tab, M_exon_tab = pg.prot_genome_loc(rec_raw_prot, ens=True)  

**2.2	Merge the two protein databases** 

	 prot_tab, exon_tab = pg.merge_prot_exon_tabs(prot_tab_M1, exon_tab_M1, M_prot_tab , M_exon_tab)

**2.3	Generate the protein heat map**

    protein_hm_filename  = 'protein_heat_map.bed'

	pg.proteome_Hmap_BED(prot_tab,exon_tab, protein_hm_filename, log_transf='2', rev_col_gradient=False)

             



#### Folder TestFiles:
- **PoGo_peptides.bed**      :  This file contains the peptide map performed by PoGo
- **PoGo_peptides_PTM.bed**  :  This file contains the map of peptide PTMs performed by PoGo
- **proteome_test.csv**      :  This file contains a PLGS table with peptides data from MS/MS analysis