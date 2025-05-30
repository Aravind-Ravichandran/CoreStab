import subprocess
import sys
import re
import os
import shutil

# Mapping of three-letter to one-letter amino acid codes
three_to_one = {
    'ALA': 'A', 'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'MET': 'M', 'PHE': 'F', 
    'TRP': 'W', 'PRO': 'P', 'GLY': 'G', 'SER': 'S', 'THR': 'T', 'CYS': 'C', 
    'TYR': 'Y', 'ASN': 'N', 'GLN': 'Q', 'ASP': 'D', 'GLU': 'E', 'LYS': 'K', 
    'ARG': 'R', 'HIS': 'H'
}

# Hydrophobic residues
hydrophobic_residues = ['A', 'V', 'L', 'I', 'F', 'W', 'M']

def clean_pdb(pdb_file):
    cleaned_lines = []
    chain_counter = 0
    current_chain = None
    chain_mapping = {}

    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith('ATOM'):
                chain_id = line[21]
                if chain_id == ' ':
                    if current_chain is None:
                        current_chain = chr(ord('A') + chain_counter)
                        chain_counter += 1
                    line = line[:21] + current_chain + line[22:]
                else:
                    if chain_id not in chain_mapping:
                        chain_mapping[chain_id] = chr(ord('A') + chain_counter)
                        chain_counter += 1
                    current_chain = chain_mapping[chain_id]
                    line = line[:21] + current_chain + line[22:]
                cleaned_lines.append(line)
            elif line.startswith('HETATM'):
                res_name = line[17:20].strip()
                if res_name == 'HOH':
                    continue
                chain_id = line[21]
                if chain_id == ' ':
                    if current_chain is None:
                        current_chain = chr(ord('A') + chain_counter)
                        chain_counter += 1
                    line = line[:21] + current_chain + line[22:]
                else:
                    if chain_id not in chain_mapping:
                        chain_mapping[chain_id] = chr(ord('A') + chain_counter)
                        chain_counter += 1
                    current_chain = chain_mapping[chain_id]
                    line = line[:21] + current_chain + line[22:]
                cleaned_lines.append(line)

    cleaned_pdb_file = pdb_file.replace('.pdb', '_cleaned.pdb')
    with open(cleaned_pdb_file, 'w') as file:
        file.writelines(cleaned_lines)
    
    return cleaned_pdb_file

def run_freesasa(pdb_file):
    result = subprocess.run(['freesasa', pdb_file, '--format=seq'], capture_output=True, text=True)
    return result.stdout

def parse_freesasa_output(output):
    core_residues = []
    for line in output.split('\n'):
        match = re.match(r'SEQ\s*(\w?)\s*(\d+)\s*(\w+)\s*:\s*([\d\.]+)', line)
        if match:
            chain, res_num, res_name, sasa = match.groups()
            chain = chain if chain else 'A'  
            if float(sasa) < 1.0:
                core_residues.append((chain, res_num, res_name))
    return core_residues

def generate_mutation_list(core_residues):
    mutations = []
    mutation_dict = {}

    for chain, res_num, res_name in core_residues:
        original_res = three_to_one.get(res_name, 'X')
        for target_res in hydrophobic_residues:
            if target_res != original_res:
                mutation = f'{original_res}{chain}{res_num}{target_res};'
                mutations.append(mutation)
                mutation_dict.setdefault(f'{chain}-{res_name}{res_num}', []).append(target_res)

    with open('individual_list.txt', 'w') as f:
        f.write('\n'.join(mutations))

    return mutations, mutation_dict

def run_foldx(pdb_file):
    with open('foldx_log.txt', 'w') as log_file:
        result = subprocess.run(
            [
                'foldx',
                '--command=BuildModel',
                '--pdb=' + pdb_file,
                '--mutant-file=individual_list.txt',
                '--numberOfRuns=5'  # Add this option
            ],
            stdout=log_file,
            stderr=log_file
        )
    
    # Check if FoldX ran successfully
    if result.returncode != 0:
        raise RuntimeError(
            f"FoldX failed with return code {result.returncode}. "
            "Check 'foldx_log.txt' for details."
        )

def parse_foldx_output(pdb_file, mutations):
    base_name = os.path.splitext(os.path.basename(pdb_file))[0]  # Get the base name (e.g., 'ub')
    output_file = f'Average_{base_name}.fxout'  # Construct the expected output file name

    # If the file doesn't exist, try the cleaned PDB file name
    if not os.path.exists(output_file):
        cleaned_base_name = base_name.replace('_cleaned', '')  # Remove '_cleaned' if present
        output_file = f'Average_{cleaned_base_name}.fxout'

    # If the file still doesn't exist, raise an error
    if not os.path.exists(output_file):
        raise FileNotFoundError(
            f"FoldX output file not found. Tried: 'Average_{base_name}.fxout' and 'Average_{cleaned_base_name}.fxout'. "
            "Ensure FoldX completed successfully and generated the expected output."
        )

    negative_ddg_mutations = []

    with open(output_file, 'r') as file:
        lines = file.readlines()
    
    # Extract mutation mapping from FoldX output
    mutation_map = {}
    for i, line in enumerate(lines):
        if line.startswith(f'{base_name}_'):
            parts = line.strip().split()
            if len(parts) > 2:  # Ensure there are enough columns
                try:
                    ddg = float(parts[2])  # Total energy is in the 3rd column
                    mutation_map[i] = ddg
                except ValueError:
                    continue
    
    # Match mutations with correct FoldX output lines
    mutation_index = 0
    for i in sorted(mutation_map.keys()):
        if mutation_index < len(mutations):
            ddg = mutation_map[i]
            if ddg < 0:
                negative_ddg_mutations.append((mutations[mutation_index], ddg))
            mutation_index += 1
    
    return negative_ddg_mutations

def create_run_directory(pdb_file):
    base_name = os.path.splitext(os.path.basename(pdb_file))[0]
    dir_name = base_name
    counter = 1

    # Find a unique directory name
    while os.path.exists(dir_name):
        dir_name = f"{base_name}_{counter}"
        counter += 1

    # Create the directory
    try:
        os.makedirs(dir_name, exist_ok=True)
    except Exception as e:
        raise RuntimeError(f"Failed to create directory '{dir_name}': {e}")

    return dir_name

def main(pdb_file):
    # Create a unique directory for this run
    run_dir = create_run_directory(pdb_file)
    original_dir = os.getcwd()
    os.chdir(run_dir)

    # Copy the PDB file to the run directory
    shutil.copy(os.path.join(original_dir, pdb_file), '.')

    # Clean the PDB file
    cleaned_pdb_file = clean_pdb(pdb_file)

    # Run the rest of the pipeline
    sasa_output = run_freesasa(cleaned_pdb_file)
    core_residues = parse_freesasa_output(sasa_output)
    mutations, mutation_dict = generate_mutation_list(core_residues)
    run_foldx(cleaned_pdb_file)  # Updated function with --numberOfRuns=5
    negative_ddg_mutations = parse_foldx_output(cleaned_pdb_file, mutations)

    print("Residues selected for mutation from PDB:")
    for chain, res_num, res_name in core_residues:
        print(f"{chain}-{res_name}{res_num}")

    print("Mutants created:")
    for key, value in mutation_dict.items():
        print(f"{key} --> {''.join(value)}")

    print("Mutations with negative ??G values:")
    for mutation, ddg in sorted(negative_ddg_mutations, key=lambda x: x[1]):
        print(f'{mutation}: {ddg}')

    # Return to the original directory
    os.chdir(original_dir)

if __name__ == '__main__':
    main(sys.argv[1])
