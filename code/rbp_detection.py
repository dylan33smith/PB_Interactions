import os  # Used for file management (removing files, path operations)
import pandas as pd  # Used for handling tabular data (reading/writing CSV files)
import subprocess  # Used to run PHANOTATE as an external process
from tqdm import tqdm  # Used to display a progress bar for tracking processing
from os import listdir  # Used to list files in a directory
from Bio import SeqIO  # Used to read and parse FASTA files
from Bio.Seq import Seq  # Used for reverse complementing sequences
from Bio.SearchIO import HmmerIO
import glob

############################only needs to be run once ##############################
def hmmpress(hmm_path, pfam_file):
    """
    Prepares an HMM profiles database (pfam_file) for efficient querying by HMMERs hmmscan.
    - This is a one-time setup.
    - hmmpress compressses and indexes the .hmm file, generating auxilary files (.h3m, .h3i) needed for fast searches.

    Inputs:
        hmm_path: path to HMMER software
        pfam_file: Path to the HMM profiles 
    """
    cd_str = "cd " + hmm_path 
    press_str = 'hmmpress ' + pfam_file
    command = cd_str + ';' + press_str
    press_process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, _ = press_process.communicate()
    if press_process.returncode != 0:
        raise RuntimeError(f"Error running hmmpress: {stdout.decode()} ")


def single_hmm_scan(hmm_path, pfam_path, fasta_file):
    """
    Does a hmmscan for a given FASTA file of one (or multiple) sequences,
    against a given profile database. Assuming an already pressed profiles
    database (see function above).
    
    INPUT: all paths to the hmm, profiles_db, fasta and results file given as strings.
            results_file should be a .txt file
    OUPUT: ...
    """

    cd_str = 'cd ' + hmm_path
    cd_process = subprocess.Popen(cd_str, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    cd_out, cd_str = cd_process.communicate()

    # scan the sequences
    scan_str = 'hmmscan ' + pfam_path + ' ' + fasta_file + ' > hmmscan_out.txt'
    scan_process = subprocess.Popen(scan_str, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    scan_out, scan_err = scan_process.communicate()

    # process output
    results = open('hmmscan_out.txt')
    scan_results = HmmerIO.Hmmer3TextParser(results)
    os.remove('hmmscan_out.txt')
    return scan_out, scan_err, scan_results

def process_protein_fastas(dir_path, hmm_path, pfam_path):
    """
    Runs HMM scan on all protein FASTA files in the specified directory structure.
    
    Args:
        dir_path (str): Base directory path
        hmm_path (str): Path to HMMER software
        pfam_path (str): Path to HMM profiles database
    
    Returns:
        dict: Dictionary containing results for each file
    """
    # Get all protein fasta files in the data directory
    fasta_files = glob.glob(os.path.join(dir_path, 'data', 'annotated', '*_data', '*_proteins.fasta'))
    
    results_dict = {}
    
    # Process each fasta file
    for fasta_file in tqdm(fasta_files, desc="Processing protein files"):
        file_name = os.path.basename(fasta_file)
        results_dict[file_name] = []
        
        try:
            scan_out, scan_err, scan_results = single_hmm_scan(hmm_path, pfam_path, fasta_file)
            hmmer_results = list(scan_results)
            
            # Store results for each query sequence
            for query in hmmer_results:
                query_results = {
                    'query_id': query.id,
                    'num_hits': len(query),
                    'hits': []
                }
                
                for hit in query.hits:
                    query_results['hits'].append({
                        'hit_id': hit.id,
                        'evalue': hit.evalue
                    })
                    
                results_dict[file_name].append(query_results)
                
        except Exception as e:
            print(f"Error processing {file_name}: {str(e)}")
            
    return results_dict

def print_hmm_results(results_dict):
    """
    Prints HMM scan results in a readable format.
    
    Args:
        results_dict (dict): Results dictionary from process_protein_fastas
    """
    for file_name, results in results_dict.items():
        print(f"\nProcessing {file_name}...")
        
        for query_result in results:
            print(f"\nQuery: {query_result['query_id']}, Hits: {query_result['num_hits']}")
            for hit in query_result['hits']:
                print(f"  Hit ID: {hit['hit_id']}, E-value: {hit['evalue']}")
        
        print("-" * 80)

if __name__ == "__main__":
    # Example usage
    dir_path = '/home/dylan33smith/projects/Yuzhen/PB_interactions'
    hmm_path = '/home/dylan33smith/src/hmmer-3.4/src'
    pfam_path = '/home/dylan33smith/projects/Yuzhen/PB_interactions/hmm_files/RBPdetect_phageRBPs.hmm'
    
    results = process_protein_fastas(dir_path, hmm_path, pfam_path)
    print_hmm_results(results)
