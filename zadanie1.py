#!/usr/bin/python
# -*- coding: latin-1 -*-
import requests
from Bio import PDB
import numpy as np
import matplotlib.pyplot as plt

def download_pdb(pdb_id, save_path):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    with open(save_path, 'wb') as pdb_file:
        pdb_file.write(response.content)

def calculate_distance(atom1, atom2):
    return np.linalg.norm(atom1.coord - atom2.coord)

def create_contact_map(structure, threshold=8.0):
    model = structure[0]
    chain = model['A']
    
    contact_map = np.zeros((len(chain), len(chain)), dtype=int)

    # tworzymy list� atom�w CA dla ka�dej reszty
    ca_atoms = [residue['CA'] for residue in chain if 'CA' in residue]

    #por�wnujemy ka�dy atom z ka�dym
    #unikamy przypadk�w gdzie i >= j, aby unikn�� powtarzania por�wnan�
    for i, ca1 in enumerate(ca_atoms):
        for j, ca2 in enumerate(ca_atoms):
            if i >= j:
                continue

            # Obliczmy odleg�o�ci mi�dzy dwoma atomami CA
            distance = calculate_distance(ca1, ca2)

            if distance <= threshold:
                contact_map[i, j] = 1
                contact_map[j, i] = 1
    print(contact_map)
    return contact_map

def plot_contact_map(contact_map):
    plt.imshow(contact_map, cmap='binary', interpolation='none')
    plt.title('Contact Map')
    plt.xlabel('Residue Index')
    plt.ylabel('Residue Index')
    plt.show()

def main():
    pdb_id = '1HHB'
    pdb_file_path = f'{pdb_id}.pdb'

    # Download PDB file
    download_pdb(pdb_id, pdb_file_path)

    # Load 1HHB structure
    pdb_parser = PDB.PDBParser(QUIET=True)
    structure = pdb_parser.get_structure(pdb_id, pdb_file_path)

    # Create contact map
    contact_map = create_contact_map(structure)

    # Plot contact map
    plot_contact_map(contact_map)

if __name__ == "__main__":
    main()