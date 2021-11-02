"""
Module Name : Proteogenome
Author      : Giammarco Ferrari
Version     : 1.0
 
"""


                                                #***********************#
                                                #    UTILITY SECTION    #
                                                #***********************#


def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Version : 1.0

    Call in a loop to create terminal progress bar
    
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)

    Referencies: https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()


def print_lst(input_list, limit):
    """
    Version : 1.0

    Print the list content limitely to the element indicated by the parameter 'limit'. 
    INPUT  : input_list     List      Items to print
             limit          Int       Number of items to print 
    
    OUTPUT :
    """
    for ind, i in enumerate(input_list):
        if ind < limit:
            print(i)


def print_dict(input_dict, limit, separator = '-'):
    """
    Version : 1.0

    Print touples  key / value  of a dictionary separated by the variable separator.
    INPUT  : input_dict     Dictionary  to print   
             limit          Int         Number of items to print
             separator      Str         The separation element printed between key / value of the dictionary 
    
    OUTPUT :
    """
    separator = ' ' + separator + ' '
    for ind, (key,values) in enumerate(input_dict.items()):
        if ind == limit : break
        print('{}{}{}'.format(key, separator, values))


def blockPrint():
    """
    Version 1.0

    This function disable the print output to the system console.

    INPUT  :
    
    OUTPUT :  old_std_out  The previous configuration of the stantard output. Could be useful to save it in case it is needed to restore the configuration.
    """
    import sys
    import os
    old_std_out = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    return old_std_out


def enablePrint(old_std_out):
    """
    Version : 1.0

    This function enable the print output to the system console using the standard output configuration provided.

    INPUT  : old_std_out  Standard output configuration  
    
    OUTPUT :
    """
    import sys
    sys.stdout = old_std_out



def load_proteomic_data(filename, target_cols=[], reshape=True): 
    """
    Version      : 2.0

    Name History : PLGS_load_data 
    
    This function read a .csv file and generates a pandas dataframe.
    This data frame could be subset and renamed basing on the rename_list :
                                        ['protein.Accession','peptide.seq','peptide.modification','peptide.PSMs','peptide.intensity']
    
    If the user leave reshape=True, he must provide the 5 column names in the csv table that contain the same types of data 
    reported in the rename_list. This data are mandatory in order to run the Proteogenome workflow.
     
    INPUT:  filename        String      Contains the path for the csv file.
            target_cols     List        Strings that contain the column names in the csv table, that must be renamed by this function.
            reshape         Boolean     Enable the data frame slicing.
    
    OUTPUT: proteome        DataFrame   A pandas dataframe that contains the proteome data (reshaped).
    """
    import pandas as pd

    rename_list = ['protein.Accession','peptide.seq','peptide.modification','peptide.PSMs','peptide.intensity']
    
    proteome = pd.read_csv(filename, encoding='latin1')

    if reshape:
        proteome = proteome[target_cols]   # Subset the columns of the input csv table
        proteome.columns = rename_list     # Rename the columns with the Proteogenome column names required



    return proteome



def load_raw_proteins(rough_p_filename):
    """
    Version : 1.0

    Name History: load_rough_proteins

    This function upload the content of a file in a list.
    The format of the input file is defined by the function download_proteins_set. 
    What is expected from this format is a protein record collection where each record is separated by this string ********** NEW RECORD ********** 
    
    INPUT  : rough_p_filename        String  The file path to upload.
    
    OUTPUT : rough_proteins_string   List    The list that contains the input file content.
    """
    rough_proteins = open(rough_p_filename, 'r')
    rough_proteins_string = rough_proteins.read()
    rough_proteins_string = rough_proteins_string.split('\n ********** NEW RECORD ********** \n')
    rough_proteins_string.remove('')

    return rough_proteins_string



def up_PoGo_bed(file_name, df_columns =['chrom','chromStart','chromEnd','peptide','score','strand','thickStart',
                                        'thickEnd','itemRgb','blockCount','blockSizes','blockStarts']):
    """
    Version: 1.0

    This function upload in a Pandas Dataframe a file that has a .bed format.
    Each file row is splitted by '\t' separator, then is removed the '\n' character attached to the last element of the split.
    Each row converted in list is stored in a list (generating a matrix).
    Than the matrix is converted in Pandas Dataframe using the columns titles received as parameter.
    What is expected in input is a BED file. 
    The usage of this function is to upload the peptides track generated by PoGo. 

    INPUT  :  file_name   String     Path name for the file to upload
              df_columns  List       List of strings that contains the Dataframe columns names 
    
    OUTPUT :  bed_df      Dataframe  The Pandas dataframe containing all the peptides data from the input file

    
    References: https://stackoverflow.com/questions/13784192/creating-an-empty-pandas-dataframe-then-filling-it
    """
    import pandas as pd
    
    bed_df = pd.DataFrame()
    bed_handler = open(file_name,'r')
    list_rows = []
    
    for row in bed_handler:
        row = row.split('\t')              # Create a list with all the columns of the bed format
        try: row.remove('')
        except: pass
        row[-1] = row[-1].replace('\n','') # The last colum has the \n delimiter. Remove it!!!!!
        list_rows.append(row)
    
    bed_df = pd.DataFrame(list_rows, columns=df_columns) 
    
    return bed_df



def make_tab_file(out_file_name, input_lst):
    """
    Version : 1.0

    This function create a tab separated file from a list of data.
    Input list must be organized as list of lists. 
    Example: [[...],
              [...],
                .
              [...]]
    The function retrieve the number of columns looking at the number of elements stored in the first element (a list by itself) stored in the input list.
    The usage of this function is the creation of BED file

    INPUT  :  out_file_name     String   The file path where the new file will be created
              input_lst         List     Contains lists with data to save in the output file
    
    OUTPUT : 
    """
    
    out_file_hand = open(out_file_name, 'w')
    number_of_columns = len(input_lst[0])

    for row in input_lst:
        out_row =''
        for col_ind, col_data in enumerate(row):
            out_row += col_data
            if (col_ind < number_of_columns): out_row += '\t'
        out_row += '\n'
        out_file_hand.write(out_row)
    out_file_hand.close()    



def file2list(file_name):
    """
    Version : 1.0 

    This function upload in a list the file rows.
    For each row it is necessary to remove the last character that usually is a ('/n'). 

    INPUT  :  file_name   String   The name of the file to upload.
    
    OUTPUT :  file_in_lst List     The output list with the content of the file.
    """
    file_in_lst = []
    file_hand = open(file_name,'r')

    for row in file_hand:
        file_in_lst.append(row[:-1])

    return file_in_lst


def list2file(file_name, input_list):
    """
    Version : 1.0

    This function writes the list content in a file.
    For each row it is necessary to remove the last character ('/n').
    The variable acc indicates that the main usage is to save lists of protein accessions. 

    INPUT  :    file_name   String   The name of the file where the list content will be write.
                input_list  List     The list to save in the file.
    
    OUTPUT :  
    """
    fh = open(file_name, 'w')
    for acc in input_list:
        fh.write(acc + '\n')
    fh.close()


def dict2file(input_dict, dict_filename, sep= ' - '):
    """
    Version : 1.0

    This function write the dictionary content in a file.
    The key/value dictionary touples are separated by the string sep

    INPUT  : input_dict     Dictionary      The dictionaary to save in the file.
             sep            String          The separator string between each key and value in the dictionary
             dict_filename  String          The name of the output file 
    
    OUTPUT :
    """
    fh = open(dict_filename, 'w')
    for k, value in input_dict.items():
        fh.write(str(k) + sep + str(value) + '\n')
    fh.close()



def file_dict2dict_vlst(filename, sep = ' - '):
    """
    Version : 1.0

    This function uploads in a dictionary the content of a file.
    The dictionary touples in the file contain a list in the value field.
    It is used in case of uploading a list of protein with their exons from a file.
    In the file the data must have a format like that:
                        Q6IMN6 - ['ENSE00001522921','ENSE00003688381',.......]

    INPUT  : filename    String      The name of the source file.
             sep         String      The separator string between each key and value in the dictionary.
    
    OUTPUT : out_dict    Dictionary  [Key]      The protein accession number.
                                     [Value]    The list of exon codes.

    """
    import re

    ens_pat = re.compile(r".*?\'(ENS.*?\d)\'")
    out_dict = {}
    input_fh = open(filename, 'r')
    
    for row in input_fh:
        row = row.split(sep)
        #print(row)
        UniP_cod = row[0]
        codes_str = row[1]
        codes_lst = ens_pat.findall(codes_str)
        out_dict[UniP_cod] = codes_lst

    return out_dict



def file_dict2dict(filename, sep = ' - '):
    """
    Version : 1.0

    This function uploads in a dictionary the content of a file.
    Compared to the function file_dict2dict_vlst, this time the value contains a string.

    INPUT  : filename    String      The name of the source file.
             sep         String      The separator string between each key and value in the dictionary.

    OUTPUT : out_dict    Dictionary  [Key]      
                                     [Value]    

    """
    out_dict = {}
    input_fh = open(filename, 'r')

    for row in input_fh:
        row = row.split(sep)
        out_dict[row[0]] = row[1].replace('\n','')  # remove the '\n' from the file row 

    return out_dict

def generate_prot_expr_pept_dict(rough_prots, ens=False):
    """
    Version : 1.0

    This function receives a list of protein records. Then generates two dictionary one for the protein expression 
    and another for the peptides that belong to each protein.
    What is expected from each record is a format like that:

                        Protein accession code § Protein expression § Peptides § Raw protein record 

    INPUT  : rough_prots     List           Contains the protein record downloaded from the UniProtKB
             ens             Boolean        If True indicates that in the first position of the rough record there is not a UniProt code (instead an Ensembl).
                                            This boolean is used to change the parsing pattern in the function and look for the accession code 
                                            in the protein record, instead of in the first group of the pat_header

    OUTPUT : prot_expr_dict  Dictionary     [Key]   :   String      Protein accession code.
                                            [Value] :   String      Protein expression level in string format.

             prot_pept_dict  Dictionary     [Key]   :   String      Protein accession code.
                                            [Value] :   List        String with the set of peptides seuences that belong to the specific protein.
    """
    import re, copy
    
    if isinstance(rough_prots, dict) == True: # If rough_prots is a dictionary, then convert the dictinary in a list of rows
        rough_prots_lst = []
        for k, v in rough_prots.items():
            rough_prots_lst.append(' § '.join([k,str(v)]))
        rough_prots = copy.deepcopy(rough_prots_lst)
    
    pat_header    = re.compile(r'(.*?)\s§\s(.*?)\s§\s\[(.*?)\]')   # Pattern for record that belong to UniProt accession
    pat_accession = re.compile(r'.*?accession:(.*?),')             # Pattern for record that come from uniprot but have been associated with Ensembl codes (usually after the conversion)

    prot_expr_dict = {}         # Dictionary of protein expression levels
    prot_pept_dict = {}         # Dictionary of protein with their peptides list

    for row in rough_prots:
        row = row.replace('\"','')
        empty_record = False
        prot_data = pat_header.match(row)
        if ens:
            try:
                ens_accession = pat_accession.match(row)   # If ens is true do not consider the accession code at the beginning of the row
                prot_acc  = ens_accession.group(1)         # Instead use the accession present in the field 'accession' in the protein record 
            except:
                empty_record = True            
        else:
            prot_acc  = prot_data.group(1)    
            
        if not empty_record:
            prot_expr = prot_data.group(2)                            # Protein expression value       
            prot_pept = prot_data.group(3)
            prot_pept = prot_pept.replace("\'","").replace(' ','')
            prot_pept = prot_pept.split(',')                          # List of peptides for the current protein
            prot_expr_dict[prot_acc] = prot_expr
            prot_pept_dict[prot_acc] = prot_pept
    
    return prot_expr_dict,prot_pept_dict


def merge_prot_exon_tabs(main_prot_tab,main_exon_tab,prot_to_merge,exon_to_merge):
    """
    This function merge the rows of 2 proteins and 2 exons tables.
    It is used when are performed Method 2 and Method 3 in Proteogenome workflow.
    The formats of the dictionary values in input are:
                                                main_prot_tab | prot_to_merge
                                                
    'Protein accession code' : ['chromosome', 'genomic start', 'genomic end', 'strand', [exon list], [peptide list], 'expression level', 'RGB code for expression level']
        
                                                
                                                main_exon_tab | exon_to_merge
                                                
                                   'Ensembl exon code': ['genomic start', 'genomic end']
                                   
    When a UniProt code from the table to merge is already in the main table, then the function puts a tag to indicate that this is another version of the protein.
    This happend when the UniProt codes are converted in Ensembl codes. Indeed, an Ensembl code cound refer to multiple UniProt codes.
    
    INPUT  : main_prot_tab      Dictionary      It is the collection of data for the main protein set. Usually obtained through the Method 1 in Proteogenome 
                                                [Key]     String               UniProt accession code
                                                [Value]   List of Strings      refer to the formats above
             
             main_exon_tab      Dictionary      It is the collection of the exons that contribute to create the proteins in main_prot_tab
                                                [Key]     String               Ensembl exon code
                                                [Value]   List of strings      
                 
             prot_to_merge      Dictionary      It is the collection of data for the protein set obtained through the Method 2 or 3 in Proteogenome 
                                                Data structure equal to  main_prot_tab
                                           
             exon_to_merge      Dictionary      It is the collection of the exons that contribute to create the proteins in prot_to_merge
                                                Data structure equal to  main_exon_tab
                                           
                                           
    OUTPUT : prot_tab_complete  Dictionary    The protein dictionary merged with the same data structure. 
             exon_tab_complete  Dictionary    The exon dictionary merged with the same data structure. 
    """
    import copy
    prot_tab_complete = copy.deepcopy(main_prot_tab)
    exon_tab_complete = {}
     
    exon_tab_complete= dict(list(main_exon_tab.items()) + list(exon_to_merge.items()))   # The exons are merged directly in one dictionary

    for accession,prot_data in prot_to_merge.items():   
        count_isof = 0
        if accession in main_prot_tab:
            new_acc = copy.deepcopy(accession)
            while new_acc in main_prot_tab:
                count_isof+=1
                original_acc = new_acc.split('-')[0]
                new_acc = original_acc + '-' + str(count_isof)
            prot_tab_complete[new_acc] = prot_data

        else:
            prot_tab_complete[accession] = prot_data
          
            
    return prot_tab_complete, exon_tab_complete

                                                #************************************#
                                                #   PLGS DATA MANIPULATION SECTION   #
                                                #************************************#






def get_proteins_expr(proteome, columns=['protein.Accession','peptide.intensity'], keys=False):
    """
    Version  :   2.0
    This function generate a dictionary of proteins codes with their expression level.
    It is supposed that the expression level might be calculated from the sum of the peptides intensities already present in the dataframe as Proteogenome workflow input.
    Because the function uses a groupby instruction on a df, it needs the names of columns to group and sum. As default there are the names of PLGS table
    
    INPUT  : proteome   Dataframe   The proteogenomic data as input for the Proteogenome workflow.
             columns    List        Strings with the column names for the groupby instruction in the dataframe

    OUTPUT : prot_exp   Dictionary  [Key]    :  String      Protein accession code
                                    [Values] :  String      Expression value


    """
    prot_peptide_intensitie = proteome[columns]  # Subset the df to the proteins - peptides - peptides intensities columns.
    prot_exp = prot_peptide_intensitie.groupby(columns[0])[columns[1]].sum().to_dict() # Group by proteins codes and sum the peptides intensities to have the proteins expressions
    return prot_exp



def get_proteins(proteome, keys=False, col_prot_pept=['protein.Accession', 'peptide.seq']):
    """
    Version : 1.0

    This function fetch all the proteins UniProt accession codes for the proteins identified in the MS analysis.
    Each accession code is reported only one time in as a dictionary key. For each key is reported a list of peptides that refer to the 
    specific protein.

    INPUT:  proteome                DataFrame       A pandas dataframe that contains the proteome data
            keys                    Bool            Flag to enable the output limited to the list of proteins accession 
    
    OUTPUT: prot_ident_grouped      Dictionary      [Key]  :  String    UniProt accession for specific protein
                                                    [Item] :  List      peptides identified in this protein
    """
    proteins_identified = proteome.loc[:,col_prot_pept]
    prot_ident_grouped = proteins_identified.groupby(col_prot_pept[0])[col_prot_pept[1]].apply(list).to_dict() # Transform the .groupby result in a dictionary of lists
    if keys == True:
        prot_ident_grouped = prot_ident_grouped.keys() # If the user wants only a list of proteins accession
    return prot_ident_grouped



def get_peptides(proteome):
    """
    Version : 1.0

    This function fetch all the peptides identified in the MS/MS analysis.
    Each accession code is reported only one time in as a dictionary key. For each key is reported a list of peptides that refer to the specific protein.

    INPUT:  proteome                DataFrame       A pandas dataframe that contains the proteome data
    
    OUTPUT: pep_ident_grouped       Dictionary      [Key]  :  String    The peptide sequence
                                                    [Item] :  List      of proteins where the peptides has been identified

    Referncies: https://stackoverflow.com/questions/29876184/groupby-results-to-dictionary-of-lists
    """
    peptides_identified = proteome.loc[:,['protein.Accession', 'peptide.seq']]
    pep_ident_grouped = peptides_identified.groupby('peptide.seq')['protein.Accession'].apply(list).to_dict() # Transform the .groupby result in a dictionary of lists

    return pep_ident_grouped





                                                #***************************#
                                                #   QUERY DBs API SECTION   #
                                                #***************************#

def UniParc_API(accession_code): # UniProt2Ensembl 
    """
    Version: 1.0

    This function access the UniParc API and try to convert one accession code in an Ensembl code.
    If the conversion fails, it returns a void string.

    INPUT  :  accession_code      String      The protein accession code in UniProt format.  

    OUTPUT :  responseBody        String      The rough record from uniparc where to search for the Ensembl conversion.

    Referencies: https://www.ebi.ac.uk/proteins/api/doc/#!/coordinates/getByDbXRef
    """
    import requests, sys

    requestURL = "https://www.ebi.ac.uk/proteins/api/uniparc?accession=" + accession_code

    r = requests.get(requestURL, headers={ "Accept" : "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    responseBody = r.text  
      
    return(responseBody)



def RetrIDmap_API(prot_accession, ens_conv): # crosslink_protdb
    """
    Version: 1.0

    This function access the Retrieve/ID mapping service, to convert the UniProt protein code provided.

    INPUT  : prot_accession     String      The protein accession code in UniProt format.
             ens_conv           String      The encoded string that indicates which type of conversion is required. 

    OUTPUT : conv_list          List        The list of possible Ensembl codes for the protein accession in input.
    """
    
    import urllib.parse
    import urllib.request
   
    url = 'https://www.uniprot.org/uploadlists/'
    params = {'from':'ACC+ID',
              'to': ens_conv,
              'format': 'tab',
              'query': prot_accession
             }
    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    j=False
    conv_list = []
    try:
        req = urllib.request.Request(url, data) # Try to fetch the Ensembl codes
        with urllib.request.urlopen(req) as f:
            response = f.read()
        req_out = response.decode('utf-8')
        
        isoforms = req_out.split('\n')[1:]  # If it is present more than one isoforms expressed as Ensembl code, save them in a list. Skip the columns titles 'From To' in the first position of the list
        
        for isof in isoforms:
            conversion = isof.split('\t')[1]
            conv_list.append(conversion)
            j=True       
    except:
        if j == False:
            conversion = 'No match'
            conv_list.append(conversion)
       
    return conv_list






def UP2Ens_code_conv(prot_expre_dict, prot_peptides_dict, error_list = False, db_type='UniParc', down_sleep = 0.5, out_file = '', print_out = False): #convert_ERROR_proteins
    """
    Version: 1.0

    This function converts UniProt accession codes in Ensembl codes using the UniParc API or the Retrieve/ID map.
    Along with the codes conversion, the function creates two dictionaries:
        prot_expre_dict    ------> {accession: proteins expression, .....................}
        prot_peptides_dict ------> {accession: list of peptides found in the MS/MS,......}
    The information contained in these dictionaries will be merged with the raw protein records by the function download_proteins_set.


    FUNCTIONS USED : UniParc_API - RetrIDmap_API

     
    INPUT  : prot_expre_dict        Dictionary  [Key]   :       String      UniProt accession code.
                                                [Value] :       String      The protein expression level.
             prot_peptides_dict     Dictionary  [Key]   :       String      UniProt accession code.
                                                [Value] :       String      The list of peptides that belong to the protein accession.
             error_list             Boolean     Enable the output of the UniProt codes that have generated an error. 
                                                For the error there are two possible causes: the code was invalid or does not exist
                                                any Ensembl conversion.
             db_type                String      The type of service that must be used to perform the conversion from UniProt to Ensembl.
             down_sleep             Integer     The pause suggested for the programmatic accession to UniProt. Pre set to 1 second.
             out_file               String      The file name for saving the list proteins NOT CONVERTED (EnsemblProteinERR).
             print_out              Boolean     Enable the console text messages that show the ongoing download. 
    
    OUTPUT : Ensembl_exp_lvl        Dictionary  [Key]    :     String     The Ensemble code from the conversion operation. 
                                                [Values] :     String     The level of expression of the protein.
             Ensembl_peptide        Dictionary  [Key]    :     String     The Ensemble codes from the conversion operation. 
                                                [Values] :     List       The list of peptides sequnces that belong to the protein accession.

             EnsemblProteinERR      List        Stores The UniProt codes that generated an error during the conversion. Two possible causes the code was invalid or
                                                does not exist the Ensembl conversion for those UniProt code.
             
    
    Referencies: https://www.ebi.ac.uk/proteins/api/doc/#!/coordinates/getByDbXRef
    """

    import re
    import time
   
    ensembl_pat = re.compile(r'.*?\"Ensembl\",\"id\":\"(.*?)\"') 
   
    Ensembl_exp_lvl = {}
    Ensembl_peptide = {}
    EnsemblProteinERR = [] 
    

    tot_UP_input = len(prot_expre_dict)
    ind = 0
    if db_type=='UniParc':  print_dbtype = 'UniParc'
    elif db_type== 'IDmap': print_dbtype = 'Retrieve / ID map'
    print('CONVERSION PROCEDURE FROM UNIPROT TO ENSEMBL - ', print_dbtype)
    print('{} - UniProt codes to convert in Ensembl codes'.format(len(prot_expre_dict)))
    print('------------------------------------------------------------------------ \n')


    start = time.time()
    for ind, (acc,exp_lvl) in enumerate(prot_expre_dict.items()):
    
    #for acces_code in UP_list:           # If enabled print the progress bar
        ind += 1       
        rough_code = ''                
        
        # ------------------ UNIPARC
        if db_type == 'UniParc':
            try:
                rough_code = UniParc_API(acc)              # Try to convert the UniProt accession in an Ensembl accession using UniParc API.
            except:
                pass
            try:
                Ensembl_conversion = ensembl_pat.match(rough_code).group(1)
                Ensembl_exp_lvl[Ensembl_conversion] = prot_expre_dict[acc]               # Update the dictionary of protein - expression
                Ensembl_peptide[Ensembl_conversion] = prot_peptides_dict[acc]            # Update the dictionary of protein - peptides
                print('{} - {} - {} '.format(ind, acc, Ensembl_conversion), '                               ' * 1000, end="\r")
            except:
                EnsemblProteinERR.append(acc)            # If the conversion does not exist save the accession code in the list of MISMATCHED codes.
                print('{} - {} -  Ensembl conversion not found'.format(ind, acc), '             ' * 1000, end="\r") 
        
        # ------------------ RETRIEVE / ID MAP        
        elif  db_type == 'IDmap':                        # Try to convert the UniProt accession in an Ensembl accession using Retrieve / ID map.
                try:
                    Ensembl_conversions = RetrIDmap_API(acc,'ENSEMBL_PRO_ID')
                except:
                    pass        
                if Ensembl_conversions[0] == 'No match':
                    EnsemblProteinERR.append(acc)
                    print('{} - {} -  Ensembl conversion not found'.format(ind, acc), '             ' * 1000, end="\r") 

                else:
                    print('{} - {} - {} '.format(ind, acc, Ensembl_conversions), '                               ' * 1000, end="\r")

                    for Ens_cod in Ensembl_conversions:                            # Because could be more than one isofomr save each of them with the expression and peptides from the original UniProt code
                        Ensembl_exp_lvl[Ens_cod] = prot_expre_dict[acc]            # Update the dictionary of protein - expression
                        Ensembl_peptide[Ens_cod] = prot_peptides_dict[acc]         # Update the dictionary of protein - peptides       


        time.sleep(down_sleep)    # Pause before the next download

    end = time.time()
    print('\n')
    print('Elapsed time - {}'.format((end-start)//60))

    if out_file:                                    # If requested, write the output file with the codes not converted
        file2list(out_file,EnsemblProteinERR)

    if not error_list: return Ensembl_exp_lvl, Ensembl_peptide
    else: return Ensembl_exp_lvl, Ensembl_peptide, EnsemblProteinERR



def fetch_gene_coordinates(p_code, db_type='UPKB'):   
    """
    Version : 1.5

    This function use a protein accession code to retrieve the protein annotation record from UniProtKB.
    The user can query the database using two accession code formats:
      - UniProtKB  ----->   'UPKB'
      - Ensembl    ----->   'Ens'
    There must be a correspondence between input code format and method of querying.
    For instance, query the database with an Ensembl code selecting db_type='UPKB' will generate an empty responseBody.
    
    INPUT:      p_code          String     Could contain alternatively the UniProt or Ensembl accession code 

    OUTPUT:     responseBody    String     The record fetched from UniProtKB

    Referencies:    https://www.ebi.ac.uk/proteins/api/doc/#!/coordinates/getByDbXRef
    """
    import requests, sys

    if (db_type == 'UPKB'):
        requestURL = "https://www.ebi.ac.uk/proteins/api/coordinates/" + p_code # Attach the UniProt accession to the URL query string
    elif (db_type == 'Ens'):
        requestURL = "https://www.ebi.ac.uk/proteins/api/coordinates?ensembl=" + p_code # Attach the Ensembl accession to the URL query string

    r = requests.get(requestURL, headers={ "Accept" : "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    responseBody = r.text
    
    return responseBody


def extract_PTMs_type(proteome_tab, col_target, PTMs_remove):
    """
    Version : 1.0

    This function extracts a set of PTMs type contained into a Pandas dataframe.
    Because the tabels format is not defined, the user must indicates in which colum are stored the encoding for the PTMs (col_target).

    Moreover, cleanes the final set of PTMs found from the not desired modifications.

    INPUT  :  proteome_tab   Dataframe   A df with all the MS/MS data (PLGS table). 
              col_target     List        The column name in the proteome_tab that contains the PTMs. In PLGS is equal to ['peptide.modification']
              PTMs_remove    List        Set of strings with the name of PTMs to remove from the final PTMs set. In PLGS is equal to ['None','Carbamidomethyl']
    
    OUTPUT :  PTM_types      List        Set of string with the name of PTMs present in the peptides stored in proteome_tab
    """

    unique_modif = proteome_tab[col_target] # Consider only the column of the dataframe where to perform the filtering operation
                                            # Apply the unique function but the result might be something like that [Oxidation M(6), Oxidation M(10), Oxidation M(1),  ...] because the PTM |amino acid(position)| part might affect the uniqueness of the PTM type (Oxidation in this example) 
    
    PTM_types = []
    for ind, PTM in unique_modif.iterrows():    # Iterate over the PTMs list. Because the PTMs could appear like that [Oxidation M(6), Oxidation M(10), Oxidation M(1),  ...]
                                                # we are interested only in the PTM type (in this case 'Oxidation').
        PTM_types.append(PTM[0].split()[0])     # Splitting the original PTM and taking [0] we extract always the type (  Example of split result ['Oxidation', 'M(10)']  )
    PTM_types = list(set(PTM_types))            # Apply the unique function for the list to have only one element for each PTM type.
    
     
    for mod_type in PTMs_remove:       # Loop over the list of PTMs to exclude
        try:
            PTM_types.remove(mod_type) # Remove the not_desired PTMs 
        except:
            pass
            
    return PTM_types


def apply_PTMs_to_pep(proteome_tab, PTM_col, PTMs_to_remove):
    """
    Version : 1.0

    This function updates the peptide sequences with the PTMs encoding (see PoGo software instructions).

    INPUT  :  proteome_tab              Dataframe     The table with proteomics data.
              PTM_col                   String        The name of the column in the proteomics table that contains the PTMs encoding.
              PTMs_to_remove            List          Strings with the type of PTMs that must not be applied to the peptides.   
    OUTPUT :  peptides_PTMs_updated     List          The peptides sequence updated with the PTMs encoding (format in PoGo instruction)

    References: PoGo - https://github.com/cschlaffner/PoGo 
    """

    def PTMs_type_set(modifications, PTMs_remove):

                                                        # Apply the unique function but the result might be something like that [Oxidation M(6), Oxidation M(10), Oxidation M(1),  ...] because the PTM |amino acid(position)| part might affect the uniqueness of the PTM type (Oxidation in this example) 
        PTM_types = []
        
        for ind, PTM in modifications.iterrows(): # Iterate over the PTMs list. Because the PTMs could appear like that Example: [Oxidation M(6), Oxidation M(10), Oxidation M(1),  ...]
                                                 # we are interested only in the PTM type (in the example: 'Oxidation').
            PTM_types.append(PTM[0].split()[0])  # Splitting the original PTM and taking [0] we extract always the type (  Example of split result ['Oxidation', 'M(10)']  )
        PTM_types = list(set(PTM_types))         # Apply the unique (set) function for the list in order to have only one element for each PTM type.
        

        for mod_type in PTMs_remove:   # Loop over the list of PTMs to exclude
            if (mod_type in PTM_types): PTM_types.remove(mod_type) # Remove the PTM not desired
        
        return PTM_types
    
    PTMs_code = {'phosphorylation' : '(phospho)'   , 'acetylation' : '(acetyl)', 'acetylation'      : '(acetyl)'    , 'amidation' : '(amidated)',                                  
                 'oxidation'       : '(oxidation)' , 'methylation' : '(methyl)', 'ubiquitinylation' : '(glygly; gg)', 'sulfation' : '(sulfo)'   ,                 
                 'palmitoylation'  : '(palmitoyl)' , 'formylation' : '(formyl)', 'deamidation'      : '(deamidated)'                 
                 }#Any other post-translational modification
    
    peptides_PTMs_updated = []
    peptides_PTMs_updated_2 = {}
    
    PTM_column = proteome_tab[PTM_col] # Consider the dataframe column that contains the PTM encoding 
    # Create a list of valid PTMs to insert in the peptides sequences
    PTMs_types = PTMs_type_set(PTM_column, PTMs_to_remove)
   
    for ind, db_row in proteome_tab.iterrows():             # Loop over the proteomic database reading each row
        #print('Analyse - ', db_row)        
        peptide_PTM = db_row['peptide.seq']                    # From the row take only the peptide SEQUENCE
        #print(peptide_PTM)    
        for apply_PTM in PTMs_types:                            # Loop over the PTMslist
          
           
            peptide_modification = db_row['peptide.modification']  # From the row take only the peptide MODIFICATION
            peptide_modification = peptide_modification.split()    # Split the peptide modification string in the two main components: Modification Type - Modification Position           
            modificatio_type = peptide_modification[0]             # Modification Type
                        
            try:
                modification_position = peptide_modification[1] # Modification Position. Only for those modification other than 'None'    
                if (modificatio_type.lower() == apply_PTM.lower()):                                   # Find the specific PTM allowed to be applyed
                    modification_position = int(modification_position.split('(')[1].replace(')',''))  # Clean the Modification position. Example: in general modification position might appear like that M(10) and we want only '10'
                    PTM_encoded = PTMs_code[modificatio_type.lower()]                                 # Translate the PTM name in the PoGo PTM name
                    peptide_PTM = peptide_PTM[:modification_position] + PTM_encoded + peptide_PTM[modification_position:] #Insert the PTM encoded in the peptide sequence
            except:            
                pass  # If you are on a peptide that reports 'None' in the modification field you have to skip the entire updating part
                

        peptides_PTMs_updated.append(peptide_PTM)      # Update the array of peptide sequences updated with their own PTM
        
    return peptides_PTMs_updated



def PoGo_input_df(input_df, PTMs_to_remove = ['None','Carbamidomethyl']):
    """
    Version : 2.0

    This function create a peptides dataframe with the peptide modifications embedded in the peptide sequences.
    
    
    INPU  : input_df                Dataframe
            PTMs_to_remove          List
    OUTPUT: proteome_PTM_updated    Dataframe    The dataframe that contain the peptides sequences updated with their PTM and two columns with the data
                                                 for the peptides quantitation analysis in PoGo   
    """
    PTMs = extract_PTMs_type(input_df, ['peptide.modification'], PTMs_to_remove)
    peptides_updated = apply_PTMs_to_pep(input_df, ['peptide.modification'], PTMs_to_remove)

    proteome_PTM_updated = input_df[['protein.Accession','peptide.PSMs','peptide.intensity']]

    proteome_PTM_updated.insert(1, 'peptide.seq', peptides_updated) # peptide.PTMupdated


    return proteome_PTM_updated



def PoGo_input_file(filename, input_df, experiment, exclude_list=[], FileAppend=False):
    """
    Version : 1.0
    This function create a .txt file in the file input format suitable for the PoGo software.
    The format of the file is described in the PoGo documentation: https://github.com/cschlaffner/PoGo

    INPUT : filename      String    :  The name of the file to save. Could be also the name of the file to update with more peptides 'FileAppend' flag = True
            input_df      Dataframe :  The dataframe that conatins the peptides data in the format suitable for PoGo (see the PoGo documentation).
            experiment    String    :  The tag of the experiment from which the peptides come from. 
            exclude_list  List      :  The list of columns to skip in making the file. In case the df has additional columns not included in the PoGo format
            FileAppend    Bool      :  Enable to append the peptides in input_df to the existing file (named 'filename').
    OUTPUT:   
    """
    if FileAppend:
        append_data = 'a'
    else:
        append_data = 'w'
            
    file_hand = open(filename, append_data)

    input_df = input_df[['peptide.seq','peptide.PSMs','peptide.intensity']]
    
    number_of_columns = len(input_df.columns)        # Fetch the dataframe number of columns 
    for ind1, row_df in input_df.iterrows():
        row = ''                                     # Create a row to write in the output file
        row += experiment + '\t'                     # Insert the experiment/sample tag   # Peptide sequence updated with his PTM
        for col_ind, col_data in enumerate(row_df):  # Loop over the data in the actual row
            if (col_ind in exclude_list) == False:   # If the actual column is not in the list of columns to exclude then put the data in the row to write in the file
                row += str(col_data) 
                if (col_ind < number_of_columns): row += '\t' # If it is not the last column of the actual row, then separate the nex data with a tabulator
        row += '\n'                                  # Close the row to write in the file with the new line character 
        file_hand.write(str(row))
    
    file_hand.close()



                                                #*****************************#
                                                #   VISUALIZATION   SECTION   #
                                                #*****************************#



def UP2Ens_code_conv(prot_expre_dict, prot_peptides_dict, error_list = False, db_type='UniParc', down_sleep = 0.5, out_file = '', print_out = False): 
    """
    Version : 2.0
    
    Name History : convert_ERROR_proteins

    This function converts UniProt accession codes in Ensembl codes using the UniParc API or the Retrieve/ID map.
    Along with the ensembl codes conversion, the function creates two dictionaries:
        1 ------> accession - proteins expression 
        2 ------> accession - list of peptides found in the MS/MS

    FUNCTIONS USED: UniParc_API - RetrIDmap_API

     
    INPUT  : prot_expre_dict        Dictionary  [Key]   : UniProt accession codes.
                                                [Value] : protein expression level.
             prot_peptides_dict     Dictionary  [Key]   : UniProt accession codes.
                                                [Value] : protein expression level.
             error_list             Boolean     Enable the output of the UniProt codes that have generated an error. 
             db_type                String      The type of service that must be used to perform the conversion from UniProt to Ensembl.
             down_sleep             Integer     The pause suggested for the programmatic accession to UniProt. Pre set to 1 second.
             out_file               String      The file name for saving the proteins NOT CONVERTED.
             print_out              Boolean     Enable the console text messages that show the ongoing download. 
    OUTPUT : EnsemblProtein         Dictionary  Stores the Ensemble codes from the conversion operation and the level of expression of the protein.
             EnsemblProteinERR      List        Stores The UniProt codes that generated an error during the conversion. Two possible causes the code was invalid or
                                                does not exist the Ensembl conversion for those UniProt code.
             
    
    Referencies: https://www.ebi.ac.uk/proteins/api/doc/#!/coordinates/getByDbXRef
    """

    import re
    import time
   
    ensembl_pat = re.compile(r'.*?\"Ensembl\",\"id\":\"(.*?)\"') 
   
    Ensembl_exp_lvl = {}
    Ensembl_peptide = {}
    EnsemblProteinERR = [] 
    

    tot_UP_input = len(prot_expre_dict)
    ind = 0
    if db_type=='UniParc':  print_dbtype = 'UniParc'
    elif db_type== 'IDmap': print_dbtype = 'Retrieve / ID map'
    print('CONVERSION PROCEDURE FROM UNIPROT TO ENSEMBL - ', print_dbtype)
    print('{} - UniProt codes to convert in Ensembl codes'.format(len(prot_expre_dict)))
    print('------------------------------------------------------------------------ \n')


    start = time.time()
    for ind, (acc,exp_lvl) in enumerate(prot_expre_dict.items()):
    
        ind += 1       
        rough_code = ''                
        # ------------------ UNIPARC
        if db_type == 'UniParc':
            try:
                rough_code = UniParc_API(acc)              # Try to convert the UniProt accession in an Ensembl accession using UniParc API.
            except:
                pass
            try:
                Ensembl_conversion = ensembl_pat.match(rough_code).group(1)
                Ensembl_exp_lvl[Ensembl_conversion] = prot_expre_dict[acc]               # Update the dictionary of 
                Ensembl_peptide[Ensembl_conversion] = prot_peptides_dict[acc]
                print('{} - {} - {} '.format(ind, acc, Ensembl_conversion), '                               ' * 1000, end="\r")
            except:
                EnsemblProteinERR.append(acc)            # If the conversion does not exist save the accession code in the list of MISMATCHED codes.
                print('{} - {} -  Ensembl conversion not found'.format(ind, acc), '             ' * 1000, end="\r") 
        
        # ------------------ RETRIEVE / ID MAP        
        elif  db_type == 'IDmap':                         # Try to convert the UniProt accession in an Ensembl accession using Retrieve / ID map.
                try:
                    Ensembl_conversions = RetrIDmap_API(acc,'ENSEMBL_PRO_ID')
                except:
                    pass        
                if Ensembl_conversions[0] == 'No match':
                    EnsemblProteinERR.append(acc)
                    print('{} - {} -  Ensembl conversion not found'.format(ind, acc), '             ' * 1000, end="\r") 

                else:
                    print('{} - {} - {} '.format(ind, acc, Ensembl_conversions), '                               ' * 1000, end="\r")

                    for Ens_cod in Ensembl_conversions:                            # Because could be more than one isofomr save each of them with the expression and peptides from the original UniProt code
                        Ensembl_exp_lvl[Ens_cod] = prot_expre_dict[acc] 
                        Ensembl_peptide[Ens_cod] = prot_peptides_dict[acc]                


        time.sleep(down_sleep)                      # Pause before the next download

    end = time.time()
    print('\n')
    print('Elapsed time - {}'.format((end-start)//60))

    if out_file:                                    # If requested, write the output file with the codes not converted
        file2list(out_file,EnsemblProteinERR)

    if not error_list: return Ensembl_exp_lvl, Ensembl_peptide
    else: return Ensembl_exp_lvl, Ensembl_peptide, EnsemblProteinERR

#*****************************************************************************************************************************
#*****************************************************************************************************************************
#*****************************************************************************************************************************
#*****************************************************************************************************************************



def download_proteins_set(prot_expre_dict, prot_peptides_dict, down_sleep = 0 , db_type='UPKB', print_out = False, proteins_output='', error_proteins_output=''):
    """
    Version 4.0

    This function takes as input a list of proteins accession and generate two files with proteins accessions records found and not-found
    in the database UniProtKB.
    The records are saved in the output file with additional information in this format:

    [protein ID]   §   [protein expression level]   §   [List of peptides that belog to the protein]   §   [rough record from UniProtKB (responseBody)]

    The database could be query through different codes format in a way similar to the cross-link section present in the UniProt webpages. 
    At the moment we allowing to query the database through these protein codes formats:
                                                - UniProtKB  ----->   UPKB
                                                - Ensembl    ----->   Ens
    For each code, the function access the database through 'fetch_gene_coordinates function' (using the API) and fetch the record for the protein accession code. 
    If the function finds the protein than save the data in the file-path proteins_output. Otherwise, the protein is not found in the UniProtKB then 
    the protein code is saved in the error file (error_proteins_output).
    The function can print the download progression bar and the actual protein code that is fetching. 
    The UniProt rules for programmatic access data request a delay of 1sec between two access.

    FUNCTIONS USED : fetch_gene_coordinates
    
    INPUT  : prot_expre_dict        Dictionary      [Key]   : UniProt accession codes.
                                                    [Value] : protein expression level.
             prot_peptides_dict     Dictionary      [Key]   : UniProt accession codes.
                                                    [Value] : protein expression level.
             down_sleep             Int             Time delay (in seconds) between two proteins records download events. 
             db_type                String          It is the code for the specific databe to query. Default is on UniProtKB access.
             print_out              Boolean         Enable the display of the current protein that the function is trying to fetch from the database.
             proteins_output        String          File path for the collection of proteins records downloaded.
             error_proteins_output  String          File path for the collection of proteins codes not-found in UniProt.
    OUTPUT : rough_records_lst      List
             prot_nf_dict           Dictionary      The dictionary of proteins not found and their expression levels
    """
    import requests, sys, time 
    
    stand_out = sys.stdout         #  Save the current standard output configuration

    rough_records_lst = []
    prot_nf_dict       = {}          
 
    if not print_out:  stand_out = blockPrint()         # Disable the step by step printing output by default (change to True to enable)

    print('DOWNLOAD PROTEIN RECORDS FROM UniProtKB')
    print('{} - proteins coordinates to download'.format(len(prot_expre_dict)))
    print('------------------------------------------------------------------------ \n')
    
        
    start = time.time()
    for pos, (prot_ID,prot_exp) in enumerate(prot_expre_dict.items()):
        
        peptides_list = prot_peptides_dict[prot_ID]
       
        print('{} - Fetching protein accession --- {}         '.format(pos+1, prot_ID), end="\r")
        
        try:                                                  # Try to retrieve the information for the current protein
            responseBody = fetch_gene_coordinates(prot_ID,db_type)

                                                              # If the protein is found append the expression level of the protein at the beginning of the protein record      
            responseBody = prot_ID + ' § ' + str(prot_exp) + ' § ' + str(peptides_list) + ' § ' + responseBody # Create the record structure 
            
            rough_records_lst.append(responseBody)  # append the rough protein record             
        except:
            prot_nf_dict[prot_ID] = str(prot_exp) + ' § ' + str(peptides_list)  # But if something goes wrong than put the protein code and its expression level into the dictionary of proteins not found           
        time.sleep(down_sleep)                  # Pause before the next database query         
    
                                    # If the user indicate a file where to save the records FOUND
    if proteins_output != '':
        proteins_txt = open(proteins_output,"w")                                        
        for row in rough_records_lst:
            proteins_txt.write(row)
            proteins_txt.write('\n ********** NEW RECORD ********** \n')
        proteins_txt.close()

                                    # If the user indicate a file where to save the records NOT FOUND        
    if error_proteins_output != '':
        dict2file(prot_nf_dict, error_proteins_output,' § ')
  
    end = time.time()
    print('\n\n')
    print('Elapsed time - {}'.format((end-start)//60))
    
    
    enablePrint(stand_out)
    if error_proteins_output != '':
        return rough_records_lst, prot_nf_dict
    else: 
        return rough_records_lst




def fetch_gene_coordinates_2(codes, db_type='UPKB'):   
    """
    Version: 2.0

    This function use a protein accession code to retrieve the protein annotation record from UniProtKB.
    The user can query the database using two accession code formats:
      - UniProtKB  ----->   'UPKB'
      - Ensembl    ----->   'Ens'
    There must be a correspondence between input code format and method of querying.
    For instance, query the database with an Ensembl code selecting db_type='UPKB' will generate an empty responseBody.
    
    INPUT:      p_code          String     Could contain alternatively the UniProt or Ensembl accession code 

    OUTPUT:     responseBody    String     The record fetched from UniProtKB

    Referencies:    https://www.ebi.ac.uk/proteins/api/doc/#!/coordinates/getByDbXRef
    """
    import requests, sys

    if (db_type == 'UPKB'):
        requestURL = "https://www.ebi.ac.uk/proteins/api/coordinates?offset=0&size=100&accession=" + codes # Attach the UniProt accession to the URL query string
    elif (db_type == 'Ens'):
        requestURL = "https://www.ebi.ac.uk/proteins/api/coordinates?offset=0&size=100&ensembl=" + codes # Attach the UniProt accession to the URL query string

    r = requests.get(requestURL, headers={ "Accept" : "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    responseBody = r.text
    
    return responseBody






def prot_coordinates_2(prot_code, prot_location_info):#UP_ID
    """
    Version : 2.0
   
    This function takes as input the protein record fetched from UniProtKB API.
    Then it extracts only the genomic coordinates for the protein and its exons.

    INPUT:      prot_code               String      Contains the UniProt accession code for the specific protein
                prot_location_info      String      Part of the whole UniProtKB record that contains the protein and gene coordinates.
                
    OUTPUT:     protein_record          List        Contains protein genomic coordinates. 
                exons_dict              Dictionary  [Key]   :  Ensembl exon code.
                                                    [Value] :  List with exon genomic [start, end] coordinates.
    """
    import re
    protein_record = []
    exons_dict = {}
    
# --------- GENE    
    #gene_pattern = re.compile(r'.*?chromosome.*start:(\d+),end:(\d+),')                           # Pattern for absolute gene start-end
    gen_pattern = re.compile(r'.*?chromosome:(\d+).*start:(\d+),end:(\d+),reverseStrand:(\w+)')   # Pattern for chromosome, absolute coordinate in chromosome and strand direction
    gen_start_end = gen_pattern.match(prot_location_info)
    #print('FROM prot_coordinates_2', gen_start_end)
    if gen_start_end == None:
        gen_pattern = re.compile(r'.*?chromosome:(\w),start:(\d+),end:(\d+),reverseStrand:(\w+)') # In proteins that map on sexual chromosomes the pattern for parsing coordinates is different
        gen_start_end = gen_pattern.match(prot_location_info)
                                                                                             
    protein_record = [gen_start_end.group(1), gen_start_end.group(2), gen_start_end.group(3), gen_start_end.group(4)] # ------> [chr, start, end, strand]
   
    
    # --------- EXONS   
    exon_pattern = re.compile(r'.*?\{(proteinLocation.*?),chromosome') # Extract the rough exons positions and absolute locations  
    rough_exons = exon_pattern.match(prot_location_info)
       
    
    exon_id_pattern          = re.compile(r'.*?id:(ENSE.*?)\}')                                              # Pattern for exon ID extraction
                    # NOT IN THIS FUNCTION pos_pept_pattern = re.compile(r'.*?proteinLocation:\{begin:\{position:(\d*),.*?position:(\d*)') # Pattern for in-protein positions extraction
    exon_id_notvalid_pattern = re.compile(r'.*?genomeLocation:\{position:\{position:.*?id:(ENSE.*?)\}')
    pos_geno_pattern         = re.compile(r'.*?genomeLocation:\{begin:\{position:(\d*),.*?position:(\d*),') # Pattern for exon positions extraction

    exons_ids        = exon_id_pattern.findall(rough_exons.group(1))    # exon ID
    exons_ids_remove = exon_id_notvalid_pattern.findall(rough_exons.group(1))
    if exons_ids_remove != []:
        for e_id in exons_ids_remove:
            exons_ids.remove(e_id)

    exons_genomic_locations=[]
    exons_genomic_locations = pos_geno_pattern.findall(rough_exons.group(1)) # exon absolute position IN GENOME
    
    for ind, exon_id in enumerate(exons_ids):   # Create the dictionary with exon id as index and a list with start-end genomic location
       exons_dict[exon_id] = [exons_genomic_locations[ind][0],exons_genomic_locations[ind][1]]       # { ENSE00001 : [234252, 436654] }
    
    #*******************************
    protein_record.append(exons_ids)  # Link between Protein table and exon table
    #*******************************
    
    return protein_record, exons_dict


def parse_prot_info_2(responseBody):
    """
    Version: 3.0
    This function takes in input the record resulting from fetch_gene_coordinates function and extract the protein's genomic coordinates and details. 
    The name responseBody is related to the .text method that returns the whole record of a protein as a unique string.
    Due to the fact that the raw protein records could present different formats, two parsing pattern are applied alternatively. 
    PATTERN 1 is the most used. In case it fails a second parsing attempt is performed applying PATTERN 2.

    INPUT:      responseBody        String      The raw protein record. 
    OUTPUT:     gene_det.group(1)   String      The gene accession. 
                gene_det.group(3)   String      The genomic coordinates and annotations of the protein.
                
    """    
    import re    
    UPKB = responseBody.replace('\"','') # Cancel the single quotes " to avoid overcomplicated regular expression
    
    # --- PATTERN 1
    #pattern = re.compile(r".*?\{accession:.*?,name.*?,protein.*?(gene.*?chromosome.*?,start:.*?,end:.*?,reverseStrand:\w+)")
    pattern = re.compile(r".*?\{?accession:.*?,name.*?,protein.*?(gene.*?chromosome.*?,start:.*?,end:.*?,reverseStrand:\w+)")
    gen_det = pattern.match(UPKB)
    
    # --- PATTERN 2
    if gen_det == None: # This pattern is applied if the first one failed  
        #pattern = re.compile(r".*?\{accession:.*?,name.*?,protein.*?gnCoordinate.*?(chromosome.*?,start:.*?,end:.*?,reverseStrand:\w+)")
        pattern = re.compile(r".*?\{?accession:.*?,name.*?,protein.*?gnCoordinate.*?(chromosome.*?,start:.*?,end:.*?,reverseStrand:\w+)")
        gen_det = pattern.match(UPKB) 

    # --- PATTERN 3
    try:
        if 'ENSE' not in gen_det.group(1):                                                             # Some records match the PATTERN 2 fetching wrong information 
            #pattern = re.compile(r".*?\{accession:.*?,name.*?,protein.*?gnCoordinate(.*?chromosome.*?,start:.*?,end:.*?,reverseStrand:\w+)")
            pattern = re.compile(r".*?\{?accession:.*?,name.*?,protein.*?gnCoordinate(.*?chromosome.*?,start:.*?,end:.*?,reverseStrand:\w+)")
            gen_det = pattern.match(UPKB) 
    except:
        pass

    return gen_det.group(1)



def prot_genome_loc(rough_proteins_list, ens=False, error_list=False, progbar=True):
    """
    Version: 5.0
    This function create a dictionary of proteins accession codes with their genomic coordinates.
    Moreover it manage possible errors generated by parsing the rough data from UniProtKB database.
    For each protein code it is possible to retrieve the exons structure setting the boolean variable exons_structure = True 

    FUNCTION USED : parse_prot_info - prot_coordinates

    INPUT  : rough_proteins_list   List        A list with UniProt accession codes.
             error_list            Bool        Enable the output of those accession codes that generate errors. Defaoult is disabled.
             progbar               Bool        Enable the progress bar printing on the terminal.
    OUTPUT : protein_dictionary    Dictionary  Dictionary of proteins accession codes with their genome coordinates.
                                               [Key]    String   The protein accession code.
                                               [Value]  List     An integer with the number of the chromosome where the protein's gene has been found. 
                                                                 Two integers as protein genomic coordinates.
                                                                 A Boolean with value = True - if the peptide comes from reverse translation. False in the other case. 
             exons_dictionary      Dictionary  [Key]    String   The exon Ensembl code.
                                               [Value]  List     The start-end genomic coordinates 

             prot_mismatch         Dictionary  [Key]    String   The protein accession code.
                                               [Value]  String   The type of error for the accession code.  
    """
    import re 

    explode_prot_record_pat = re.compile(r'(.*?)\s§\s(.*?)\s§\s(.*?)\s§\s(.*)')
    extract_peptid_list_pat   = re.compile(r'\'(.*?)\'')
    
    pat_accession = re.compile(r'.*?accession:(.*?),')

    protein_dictionary = {}
    exons_dictionary   = {}
    prot_mismatch      = {}
    
    tot_prot = len(rough_proteins_list)
    bar_increase = 0
   
    tar = ''

    for ind_rough_pro, rough_pro_record in enumerate(rough_proteins_list):
        if progbar:
            bar_increase += 1
            printProgressBar(bar_increase, tot_prot)            # if enable, print the progress bar

        rough_pro_record = rough_pro_record.replace('\"','')
        
        expl_ro_pro_rec = explode_prot_record_pat.match(rough_pro_record)               # Explode the rough record in -----> protein ID - intensity - peptides - rough record.
       
        if ens:
            try:
                acc_match = pat_accession.match(rough_pro_record)
                protein_ID = acc_match.group(1)
            except:
                pass
                    
        else:
            protein_ID  = expl_ro_pro_rec.group(1)                                      # Extract the protein ID 
        
        pro_expre_lvl   = expl_ro_pro_rec.group(2)                                      # Extract the protein expression lvl
        #print('2 -------- ',expl_ro_pro_rec.group(2))
        peptide_lst     = extract_peptid_list_pat.findall(expl_ro_pro_rec.group(3))     # Extract the list of peptides for this protein. Refer to the record structure for expl_ro_pro_rec.        
        #print('3 -------- ',expl_ro_pro_rec.group(3))
        rough_pro       = expl_ro_pro_rec.group(4)                                      # Extract the original rough protein record from the UniProtKB.
        #print('4 -------- ',expl_ro_pro_rec.group(4)[0:150])
        prot_details =[]
                                    
                                                                #   ---  PROTEIN RECORD  ---   #
        try:
            prot_rough_coo= parse_prot_info_2(rough_pro)
            prot_genomic_rec, append_exons = prot_coordinates_2(protein_ID, prot_rough_coo)     # chr - start - end - strand - EXONS
            prot_genomic_rec.append(peptide_lst)                                                #                                    - PEPTIDES 
            prot_genomic_rec.append(pro_expre_lvl)                                              #                                               - expression level
            protein_dictionary[protein_ID] = prot_genomic_rec
                                                                  #   ---  EXON DICTIONARY  ---   #
            exons_dictionary = dict(list(exons_dictionary.items()) + list(append_exons.items()))
        except: pass

    if error_list:                                                          # Because the record is not parsed, we lack of the protein accession. So we save the error through the index.
        return protein_dictionary, exons_dictionary, prot_mismatch
    else:
        return protein_dictionary, exons_dictionary



def filter_PoGo_BED(PoGo_bed_df, protein_tab, exon_tab, PG_BED_filtered_filename, col_prot_pept=['protein.Accession', 'peptide.seq'], progbar = False):
    """
    Version: 4.0

    This function generates a bed file format that contains only peptides mapped by PoGo that map also in the PLGS proteins coordinates.
    Looping over PoGo dataframe every peptide coordinates span is compared against the dictionary of PLGS proteins coordinates.
    If the peptide coordinates are in one protein coordinates then the peptide is stored in the bed output file.
    
    INPUT  :  PoGo_bed_df                   Dataframe   The Pandas dataframe containing all the peptides data from the PoGo .bed file 
              protein_tab                   Dictionary  Dictionary of proteins accession codes with their genome coordinates.
                                                        [Key]    String   The protein accession code.
                                                        [Value]  List     An integer with the number of the chromosome where the protein's gene has been found. 
                                                                 Two integers as protein genomic coordinates.
                                                                 A Boolean with value = True - if the peptide comes from reverse translation. False in the other case. 
              exons_tab                     Dictionary  [Key]    String   The exon Ensembl code.
                                                        [Value]  List     The start-end genomic coordinates 
              PG_BED_filtered_filename      String      Name of the file where will be stored the peptide map
             
                                               
    OUTPUT : prots           List        Proteins accession related to the peptides that match with the proteins genomic coordinates. 
             pep_filtered    List        Contains lists, where each of them is formatted like a row of a bed file.  
    """
    import copy 

    prot_gen_coo               = {}
    protein_peptide_dictionary = {}
    pep_filtered = []
    prots = []
    tot_pep = len(PoGo_bed_df)    
    bar_increase = 0  
    prot_pep_link_complete=[]
                                            
    for prot_ID, prot_record_ek in protein_tab.items():
        prot_gen_coo[prot_ID] = [prot_record_ek[0], prot_record_ek[1], prot_record_ek[2], prot_record_ek[3]]   # prot_gen_coo               - The dictionary of protein genomic coordinates                                                  
        protein_peptide_dictionary[prot_ID] = copy.deepcopy(prot_record_ek[5])                                 # protein_peptide_dictionary - The dictionary of proteins and peptides
    
    # ------ PEPTIDES
    for pep_row in PoGo_bed_df.iterrows():        # Loop over the peptides in dataframe
        

        pep_coordinates = list(pep_row[1][0:3])   # Extract only the peptide genomic coordinates
        
        
        pep_coordinates[0] = pep_coordinates[0].replace('chr','')    # Clean the chromosome number
        pep_s_e = [int(str_num) for str_num in pep_coordinates[1:3]] # Convert a list of string coordinates in integer coordinates for the comparison. Example: from this ['11', '64879630', '64854336'] to that [64879630, 64854336]
        
        if progbar:
            bar_increase += 1
            printProgressBar(bar_increase, tot_pep)
        
        pep_start = min(pep_s_e)       # Extract peptide COORDINATES and order them in a positive starnd
        pep_end   = max(pep_s_e)
        pep_strand = pep_row[1][5]     # Extract peptide STRAND
        
        
        # ------ PROTEINS
        count = 0
        for protein, prot_record in prot_gen_coo.items():             # Loop over the proteins dictionary  
            count+=1
            prot_coordinates = prot_record[1:3]
            

            prot_s_e = [int(str_num) for str_num in prot_coordinates] #prot_coordinates[1:3] # Convert a list of string coordinates in integer coordinates for the comparison. Example: from this ['11', '64879630', '64854336'] to that [64879630, 64854336]
            prot_start = min(prot_s_e)          # Extract proteins COORDINATES
            prot_end   = max(prot_s_e)
            prot_strand = prot_record[3]        # Extract protein STRAND
            if (prot_strand == 'true'):         # Convert the strand codification from UniProtKB into BED strand codification
                prot_strand = '-'
            else:
                prot_strand ='+'
                                                                                                                                                         
            if (prot_start <= pep_start) & (prot_end >= pep_end) & (prot_record[0] == pep_coordinates[0]) & (prot_strand == pep_strand): #  Check if the peptide coordinates are included in the protein coordinates               
                actual_prot_pep_link = [protein, pep_row[1][3], pep_coordinates[1], pep_coordinates[2]] # Create the connection between protein and its peptide. Example: ['A0A024R571', 'IPTAR', '951219', '951234']
                
                if pep_row[1][3] in protein_peptide_dictionary[protein]:

                    if actual_prot_pep_link not in prot_pep_link_complete:  # Check if the combination of protein and peptide with its coordinates has been already filtered
                        protein_peptide_dictionary[protein].remove(pep_row[1][3])
                        pep_filtered.append(list(pep_row[1]))               # If not store the peptide row in a list because this combination of protein and peptide has not yet been seen
                        prots.append(protein)
                        prot_pep_link_complete.append(actual_prot_pep_link) # Put the combination protein-peptide-coordinates into the list of combinations already seen

    return prots, pep_filtered        





def proteome_Hmap_BED(protein_tab, exon_tab, prot_hm_filename, log_transf='', 
                      col_gradient=['black','blue','green','orange','red'], rev_col_gradient=False): #['black','purple','mediumslateblue','blue','deepskyblue','mediumseagreen','green','limegreen','yellow','orange','red']
    """
    Version : 3.0

    This function generate a .bed file with the protein genomic locations and their level of expression.
    Into the track it is highlighted the protein exons structures. 
    The colour and intensity of the proteins is relaated to the heatmap RGB scaled in refernce to the protein expression range.

    INPUT  :    protein_tab         Dictionary  Dictionary of proteins accession codes with their genome coordinates.
                                                [Key]    String   The protein accession code.
                                                [Value]  List     An integer with the number of the chromosome where the protein's gene has been found. 
                                                                  Two integers as protein genomic coordinates.
                                                                  A Boolean with value = True - if the peptide comes from reverse translation. False in the other case. 
                exons_tab           Dictionary  [Key]    String   The exon Ensembl code.
                                                [Value]  List     The start-end genomic coordinates 
                prot_hm_filename    String          The name of the bed file as output of the function.
    OUTPUT :
    Refernces: For the List of named colors in mathplotlib - https://matplotlib.org/stable/gallery/color/named_colors.html
    """
    import math 
    
    def generate_color_gradient(color_lst, revers_gradient=False):
        """
        INPUT  :  color_lst     List    List of strings with the colour names that will be compose the colour gradient
                                        Example: ['black','gray','blue','green','yellow','red']
        OUTPUT :
        """

        def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
            c1=np.array(mpl.colors.to_rgb(c1))
            c2=np.array(mpl.colors.to_rgb(c2))
            return mpl.colors.to_rgb((1-mix)*c1 + mix*c2)
            

        def rgb_block(color1,color2, n_data_point):
            col = []
            for data_point in range(n_data_point): 
                col.append(colorFader(color1,color2,data_point/n))
            return col

        n=80
        color_blocks_lst = []
        for ind, color in enumerate(color_lst):
            color_blocks_lst.append([color, color_lst[ind+1]])
            if ind == (len(color_lst)-2): break
        gradient =[]
        for block in color_blocks_lst:
            gradient += rgb_block(block[0],block[1],n)
        if revers_gradient: gradient.reverse()
        return gradient

    def sep_k_val_dict(dictionary):
        """
        """
        keys_lst =[]
        values_lst =[]
        for keys, values in dictionary.items():
            keys_lst.append(keys)
            values_lst.append(int(values))
        return keys_lst, values_lst
    
    def exprlev_resc_RGB(values,RGB_scale):
        """
        """
        old_max = min(values)
        old_min = max(values)
        new_max = len(RGB_scale)
        new_min = 0

        rescaled_values = []

        for old_value in values:
            try:
                NewValue = int((((old_value - old_min) * (new_max - new_min)) / (old_max - old_min)) + new_min)
            except:
                print('{} - {} - {}'.format(old_max, old_min, new_min))
                a=input()
            rescaled_values.append(NewValue)

        return rescaled_values

    def merge_genloc_RGB(prot_tab, prot_exp_RGB, error_lst=False):
        """
        Version: 2.0
        Name History: merge_genloc_protexp
        
        This subfunction updates the rows in the protein table, with the RGB codes that reflects the level of expression of each protein.
        
        INPUT :  prot_tab          Dictionary   It is the protein table from the prot_genome_loc function.
                                                [Key]    String   Protein ID
                                                [Value]  List     [[protein genomic coordinates], [Exon list], [Peptidelist], level of expression]
                 prot_exp_RGB      Dictionary   Contain protein codes with their level of expression.
        OUTUT :  hm_dict           Dictionary   
                                                 [Key]    String   Protein ID
                                                 [Value]  List     Protein genomic coordinates + the RGB code for the protein expression level
        """
        hm_dict = {}
        errorID = []
        count = 0
        for prot_ID, exp_level_RGB in prot_exp_RGB.items():
            if prot_ID in prot_tab: 
                count += 1
                genloc_data = prot_tab[prot_ID]#prot_tab[prot_ID][0:3] + [exp_level_RGB]# prot_tab[prot_ID][6]  # From the whole record of the protein consider only the genomic location and expression level.
                
                if len(genloc_data) == 8:      # If the genomic data vector is long 6 elements, then it means that the color has already been assigned previously 
                    del genloc_data[7]               # Cancel the previus color  
                    genloc_data.insert(7, exp_level_RGB) # Update with the new color

                else:
                    genloc_data.append(exp_level_RGB)
                    
                
                prot_tab[prot_ID] = genloc_data
            else:
                errorID.append(prot_ID)
        
        if error_lst: return prot_tab, errorID
        else: return prot_tab
        

    def prot_genloc2BED(prot_tab, exon_tab):#,prot_exp_RGB
        """
        This subfunction create a dicionary with the data suitable for the creation of a proteins BED file 
        """
        
        BED_rough_data = []
        for prot_ID, protein_record in prot_tab.items():

            exon_lst = protein_record[4]  # Extract the list of exons for the current protein
            
            exon_gen_loc_dict = {}        # Extract from the exon table the genomic coordinates for all the exons codes related to the current protein
            for exon_id in exon_lst:
                exon_gen_loc_dict[exon_id] = exon_tab[exon_id]
                
        
    

            itemRgb = protein_record[7]        # Extract the protein expression level in RGB format
            
            coord1 = int(protein_record[1])    # Extract the protein genomic coordinates
            coord2 = int(protein_record[2])
            if coord1 < coord2: prot_chromStart = coord1   # Order the protein genomic coordinates in order to have always a positive strand
            else: prot_chromStart = coord2
            
            blockList = []
            first_exon = True
            for exon_ID, exon_details in exon_gen_loc_dict.items():  # From the exon dictionary take the coordiantes list of each exon
                
                exon_to_insert = list(map(int,exon_details))  # Convert exons coordinates from strings to integers. Example--> 'ENSE00003750855': ['64854857', '64854336'] Result in --->[64854857, 64854336] 
                                
                if first_exon:
                    blockList.append(exon_to_insert) # If it is the first exon to be considered then append directly into the block list 
                    first_exon = False
                else:
                    inserted = False 
                    index_insert = 0                      # Compare the new exon with other exons coordinates already in blockList
                    for exon_ind, exon_block in enumerate(blockList):
                        if min(exon_to_insert) < min(exon_block):     # If the start coordinate of the new exon is lower than the start coordinate of the current exon in blockList 
                            #if prot_ID==xxx: print('{} < {}'.format(min(exon_to_insert),min(exon_block)))
                            index_insert = exon_ind                   # Save the position where insert the new exon
                            inserted = True
                            break
                    if inserted: blockList.insert(exon_ind,exon_to_insert) # If the new exon is to insert, then use the index saved before to insert  
                    if not inserted: blockList.append(exon_to_insert)      # If the for loop finished and no one exon is located forward the new exon, then append the new exon
            
                if prot_chromStart > min(exon_to_insert): prot_chromStart = min(exon_to_insert)

            blockCount = len(blockList)

            blockSize_lst = []
            blockStarts_lst = []
            for block in blockList:
                blockSize_lst.append(abs(block[0]-block[1])) # Calculate the block size through the coordinates subtraction in abs
                
                blockStarts_lst.append(abs(prot_chromStart - min(block)))

            prot_totalspan_details=protein_record[0:4] 
            
            BED_blocks_details = {'prot_ID': prot_ID, 'prot_gen_det': prot_totalspan_details, 'itemRgb': itemRgb, 'blockCount': blockCount, 
                                  'blockSize': blockSize_lst, 'blockStarts': blockStarts_lst}  
        

            BED_rough_data.append(BED_blocks_details)
        return BED_rough_data


    def prot_heatmap_BED(BED_rows_lst, BED_filename):
        """
        This function receive as input a list of data suitable to fill the BED file format columns and creates a .bed file.
        """
       
        BED_file_rows =[]
        for BED_data in BED_rows_lst:
            BED_row = []
            chrom = 'chr' + BED_data['prot_gen_det'][0]
            BED_row.append(chrom)
            prot_coord = list(map(int,BED_data['prot_gen_det'][1:3]))
            chromStart = min(prot_coord)
    
            chromEnd = max(prot_coord)
    
            strand = BED_data['prot_gen_det'][3]
            if chromStart < chromEnd:
                BED_row.append(chromStart)
                BED_row.append(chromEnd)
            else:
                BED_row.append(chromEnd)
                BED_row.append(chromStart)
               
            
            tickStart = chromStart
            tickEnd = chromEnd

            name = BED_data['prot_ID']                # ---- NAME
            BED_row.append(name)
            
            
            BED_row.append('1000')

            if strand == 'true': BED_row.append('-')
            else:               BED_row.append('+')

            BED_row.append(tickStart)
            BED_row.append(tickEnd)

            itemRgb = BED_data['itemRgb']
            BED_row.append(itemRgb)


            BED_row.append(BED_data['blockCount'])
            blockSize = ''
            for bs in BED_data['blockSize']:
                blockSize += str(bs) + ','
            BED_row.append(blockSize)

            blockStarts_lst= BED_data['blockStarts']   # ---- BLOCKSTARTS
            blockStarts_lst.sort()                     # In BED file the sequence of starts should be sorted also if in reverse strand. IGV will assign the correct stat based on the strand indication.
            blockStarts = ''
            for bs in blockStarts_lst:
                blockStarts += str(bs) + ','
            BED_row.append(blockStarts)

            BED_row = list(map(str,BED_row))

            BED_file_rows.append(BED_row)
            
            
        make_tab_file(BED_filename,BED_file_rows)



    def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
        """
        #    https://stackoverflow.com/questions/25668828/how-to-create-colour-gradient-in-python
        #    https://stackoverflow.com/questions/876853/generating-color-ranges-in-python

        """
        c1=np.array(mpl.colors.to_rgb(c1))
        c2=np.array(mpl.colors.to_rgb(c2))
        return mpl.colors.to_rgb((1-mix)*c1 + mix*c2)
    
                            #****************************   MAIN FUNCTION   ***************************#
    
        
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    
    RGB_tuples=[]
    
    # Create two lists one for the proteins codes and one for the protein expression levels
    pro_ID_lst   = []
    pro_expr_lst = []
    for prot_ID, prot_details in protein_tab.items():  
        pro_ID_lst.append(prot_ID)                    # Slit protein ID and expression level in two lists  
        
        if log_transf == '':
            pro_expr_lst.append(prot_details[6])
        else:            
            pro_expr_lst.append(math.log(int(log_transf), int(prot_details[6]))) #math.log(int(prot_details[prot_ID][4]))# The expression level in the protein details is in fourth position. Example: 'A0A024R571': ['11', '64879630', '64854336', 'true', '2307']
            
    #*****************
    # CREATE RGB SCALE for the heatmap
    #*****************
          
    # CALCULATE the EXPRESSION LEVEL RANGE 
    pro_expr_lst = list(map(float,pro_expr_lst))
    max_expression =max(pro_expr_lst)   
    min_expression=min(pro_expr_lst)
    
    expression_range= int(max_expression - min_expression)
    
    #******************************************************************************************************
    RGB_tuples = generate_color_gradient(col_gradient,rev_col_gradient)
    
    fig, ax = plt.subplots(figsize=(8, 5))
    for x, colorRGB in enumerate(RGB_tuples):
        ax.axvline(x, color=colorRGB, linewidth=4) 
        # REINDEX the protein expressions value based on the RGB values (previously fitted TO THE EXPRESSIONRANGE) 
    prot_expressions_rescaled = exprlev_resc_RGB(pro_expr_lst,RGB_tuples)

    
    plt.hist(prot_expressions_rescaled, density=True, bins=30)  # density=False would make counts
    plt.show()
    
    # VECTORIZE the RGB codes based on the expression values index previously set
    prot_expressions_RGB=[]

    for RGB_pos in prot_expressions_rescaled:
        try:
            RGB_toup = RGB_tuples[RGB_pos-1]
            RGB_code = str(RGB_toup[0]) + ',' + str(RGB_toup[1]) + ',' + str(RGB_toup[1])
            prot_expressions_RGB.append(RGB_code)
        except:
            pass
            
    prot_exp_IDRGB = dict(zip(pro_ID_lst,prot_expressions_RGB))
    
    
    # MERGE the proteins genomic details with the respective RGB value to represent graphically the protein expression level
    protein_tab = merge_genloc_RGB(protein_tab, prot_exp_IDRGB)
    
    
    # Due to the complex structure of the proteins_genomic_coord it is better to create a list of list with the rows for the BED file
    prot_gen_loc_BED_rows=prot_genloc2BED(protein_tab, exon_tab)#heatmap_dict

    # Create the BED file
    prot_heatmap_BED(prot_gen_loc_BED_rows,prot_hm_filename)

