import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from code.rbp_detection import process_protein_fastas, print_hmm_results

# Define paths
dir_path = '/home/dylan33smith/projects/Yuzhen/PB_interactions'
hmm_path = '/home/dylan33smith/src/hmmer-3.4/src'
pfam_path = '/home/dylan33smith/projects/Yuzhen/PB_interactions/hmm_files/RBPdetect_phageRBPs.hmm'

# Run the analysis
results = process_protein_fastas(dir_path, hmm_path, pfam_path)
print_hmm_results(results)
