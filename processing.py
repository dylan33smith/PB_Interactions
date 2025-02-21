import os  # Used for file management (removing files, path operations)
import pandas as pd  # Used for handling tabular data (reading/writing CSV files)
import subprocess  # Used to run PHANOTATE as an external process
from tqdm import tqdm  # Used to display a progress bar for tracking processing
from os import listdir  # Used to list files in a directory
from Bio import SeqIO  # Used to read and parse FASTA files
from Bio.Seq import Seq  # Used for reverse complementing sequences

def gene_to_protein(gene_sequence):
  """
  converts gene sequence to protein sequence
  """
  dna = gene_sequence.upper()
  dna_seq_obj = Seq(dna) # convert string to Biopython Seq object
  protein_sequence = str(dna_seq_obj.translate(to_stop=True)) # convert to protein sequence
  return protein_sequence

def phanotate_processing(phage_fasta_path, phanotate_path):
    """
    Processes a single phage genome using PHANOTATE and saves the gene predictions to a CSV.

    INPUTS:
    - input_fasta (str): Path to the single FASTA file containing the phage genome.
        - path in relation to current working directory
    - phanotate_path (str): Path to the PHANOTATE executable/script.

    OUTPUT:
    - CSV file with columns ['phage_ID', 'gene_ID', 'gene_sequence'].
    """

    # extract phage name from file (remove directory and .fasta)
    phage_name = os.path.basename(phage_fasta_path).replace('.fasta', '')
    
    # run phanotate on input fasta file
    ## shell command string to call phanotate with the input fasta file
    phanotate_shell_command = f"{phanotate_path} {phage_fasta_path}"
    ## running subprocess that executes phanotate
    process = subprocess.Popen(phanotate_shell_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    ## standard output (gene predictions from phanotate)
    ## wats for process to finish execution and then returns the output
    stdout, _ = process.communicate()
    if process.returncode != 0:
        raise RuntimeError(f"Error running PHANOTATE: {stdout.decode()} on {phage_name}")

    # process phanotate output (skip headers)
    # stdout is a byte string
    std_splits = stdout.split(sep=b'\n')[2:]
    
    temp_tsv_path = 'temp_phanotate_output.tsv'
    with open(temp_tsv_path, 'wb') as temp_tab:
        for split in std_splits:
            temp_tab.write(split.replace(b',', b'') + b'\n')
    
    """with open(temp_tsv_path, 'w') as temp_tab:
        for split in std_splits:
            split = split.replace(b',',b'') # remove commas for pandas compatability
            temp_tab.write(split.decode('utf-8') + '\n')"""

    orfs = pd.read_csv(temp_tsv_path, sep='\t', lineterminator='\n', index_col=False)
    
    sequence = str(SeqIO.read(phage_fasta_path, 'fasta').seq)

    name_list = []
    gene_list = []
    protein_list = []
    gene_ids = []
    count = 1

    for j, strand in enumerate(orfs['FRAME']):
        start = orfs['#START'][j]
        stop = orfs['STOP'][j]
        
        if strand == '+':
            gene = sequence[start-1:stop]
        else:
            sequence_part = sequence[stop-1:start]
            gene = str(Seq(sequence_part).reverse_complement())

        protein = gene_to_protein(gene)

        name_list.append(phage_name)
        gene_list.append(gene)
        protein_list.append(protein)
        gene_ids.append(f"{phage_name}_gp{count}")
        count += 1

    # remove temp tsv file
    os.remove(temp_tsv_path)

    genebase = pd.DataFrame(list(zip(name_list, gene_ids, gene_list, protein_list)), columns=['phage_ID', 'gene_ID', 'gene_sequence', 'protein_sequence'])
    return genebase



