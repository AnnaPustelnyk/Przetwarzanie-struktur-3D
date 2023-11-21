#!/usr/bin/python
# -*- coding: latin-1 -*-
from Bio import PDB
import numpy as np
import requests
import csv

def download_pdb(pdb_id, save_path):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    with open(save_path, 'wb') as pdb_file:
        pdb_file.write(response.content)

def calculate_torsion_angles(residue):
    try:
        if 'O4\'' in residue and 'C1\'' in residue and 'N9' in residue and 'C4' in residue:
            chi_purines = calculate_angle(residue['O4\''].coord, residue['C1\''].coord, residue['N9'].coord, residue['C4'].coord)
            chi_purines += 180.0
        else:
            chi_purines = None

        if 'O4\'' in residue and 'C1\'' in residue and 'N1' in residue and 'C2' in residue:
            chi_pyrimidines = calculate_angle(residue['O4\''].coord, residue['C1\''].coord, residue['N1'].coord, residue['C2'].coord)
        else:
            chi_pyrimidines = None

        return chi_purines, chi_pyrimidines

    except KeyError:
        return None, None

def calculate_angle(a, b, c, d):
    #Obliczmy katy pomiedzy trzema punktami wykorzystyjac wektory
    vector_ab = b - a
    vector_cb = c - b
    vector_dc = d - c
    #Obliczamy projekcje na plaszczyznie prostopadla do wektora b
    Pa = vector_ab - np.dot(vector_ab, vector_cb) / np.dot(vector_cb, vector_cb) * vector_cb
    Pc = vector_dc - np.dot(vector_dc, vector_cb) / np.dot(vector_cb, vector_cb) * vector_cb

    #Obliczamy kat torsyjny(np.cross - iloczyn wektorowy, np.dot - iloczyn skalarny)
    angle_radians = np.arctan2(np.linalg.norm(np.cross(Pa, Pc)), np.dot(Pa, Pc))
    angle_degrees = np.degrees(angle_radians)

    return angle_degrees

def calculate_backbone_angles(residue):
    result_list = []
    try:
        alpha = calculate_angle(residue['O3\''].coord, residue['P'].coord, residue['O5\''].coord, residue['C5\''].coord)
        result_list.append(alpha)

        beta = calculate_angle(residue['P'].coord, residue['O5\''].coord, residue['C5\''].coord, residue['C4\''].coord)
        result_list.append(beta)

        gamma = calculate_angle(residue['O5\''].coord, residue['C5\''].coord, residue['C4\''].coord, residue['C3\''].coord)
        result_list.append(gamma)

        delta = calculate_angle(residue['C5\''].coord, residue['C4\''].coord, residue['C3\''].coord, residue['O3\''].coord)
        result_list.append(delta)

        epsilon = calculate_angle(residue['C4\''].coord, residue['C3\''].coord, residue['O3\''].coord, residue['P'].coord)
        result_list.append(epsilon)

        zeta = calculate_angle(residue['C3\''].coord, residue['O3\''].coord, residue['P'].coord, residue['O5\''].coord)
        result_list.append(zeta)

        chi_purines, chi_pyrimidines = calculate_torsion_angles(residue)
        result_list.append(chi_purines)
        result_list.append(chi_pyrimidines)

        return result_list
    except KeyError:
        return None

def save_to_csv(file_path, data, header):
    with open(file_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        writer.writerows(data)

pdb_id = '1EHZ'
pdb_file_path = f'1EHZ.pdb'

# Pobieranie pliku
download_pdb(pdb_id, pdb_file_path)

# Utworzenie PDB parser
parser = PDB.PDBParser(QUIET=True)

# Wczytanie struktury z pliku PDB
structure = parser.get_structure(pdb_id, pdb_file_path)

#Dane dla ka¿dej reszty
residue_data = []

for model in structure:
    for chain in model:
        for residue in chain:
            angles = calculate_backbone_angles(residue)
            if angles is not None:
                residue_data.append([residue.resname] + list(angles))

csv_file_path = f'{pdb_id}_angles.csv'
header = ['Residue'] + ['Alpha', 'Beta', 'Gamma', 'Delta', 'Epsilon', 'Zeta', 'Chi puryny', 'Chi pirymidyna']

save_to_csv(csv_file_path, residue_data, header)
print(f"Angles saved to {csv_file_path}")