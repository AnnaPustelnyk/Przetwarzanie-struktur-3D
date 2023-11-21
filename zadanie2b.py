import numpy as np
import matplotlib.pyplot as plt
from Bio import PDB
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
import requests
import Bio.PDB
import math
from RamachanDraw import fetch, phi_psi, plot
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def download_pdb(pdb_id, save_path):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    with open(save_path, 'wb') as pdb_file:
        pdb_file.write(response.content)
#konwertacja z radianow ta stopnie
def degrees(rad_angle) :
    if rad_angle is None :
        return None
    angle = rad_angle * 180 / math.pi
    while angle > 180 :
        angle = angle - 360
    while angle < -180 :
        angle = angle + 360
    return angle

def ramachandran(residue, next_residue) :
    if residue.resname.upper()=="GLY" :
        return "Glycine"
    elif residue.resname.upper()=="PRO" :
        return "Proline"
    elif next_residue is not None \
    and next_residue.resname.upper()=="PRO" :
        #jesli nastepna reszta to Prolina
        return "Pre-Pro"
    else :
        return "General"    

pdb_code = "1mbo"
pdb_file_path = f'MBO.pdb'

download_pdb(pdb_code, pdb_file_path)

import Bio.PDB
structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_file_path)
print ("Done")

phi_list = []
psi_list = []
secondary_structure_list = []

print ("About to save angles to file...")
output_file = open("%s_biopython.tsv" % pdb_code,"w")
for model in structure :
    for chain in model :
        #polipeptydy na podstawie lancucha aminokwasow
        polypeptides = Bio.PDB.CaPPBuilder().build_peptides(chain)
        for poly_index, poly in enumerate(polypeptides) :
            phi_psi = poly.get_phi_psi_list() #lista katow phi i psi dla reszt w danym polipeptydzie
            for res_index, residue in enumerate(poly) :
                phi, psi = phi_psi[res_index]
                if phi and psi :
                    phi_list.append(degrees(phi))
                    psi_list.append(degrees(psi))
                    secondary_structure_list.append(ramachandran(residue, poly[res_index+1]))
                    #Don't write output when missing an angle
                    output_file.write("%s:Chain%s:%s%i\t%f\t%f\t%s\n" \
                        % (pdb_code, str(chain.id), residue.resname,
                           residue.id[1], degrees(phi), degrees(psi),
                           ramachandran(residue, poly[res_index+1])))
output_file.close()
print ("Done")

df = pd.DataFrame({'phi': phi_list, 'psi': psi_list, 'secondary_structure': secondary_structure_list})
# Plot using Seaborn
plt.figure(figsize=(10, 8))

sns.scatterplot(x='phi', y='psi', hue='secondary_structure', data=df, palette='viridis')
plt.title('Ramachandran Plot')
plt.xlabel('Phi Angle (degrees)')
plt.ylabel('Psi Angle (degrees)')
plt.xlim(-180, 180)
plt.ylim(-180, 180)
plt.xticks(np.arange(-180, 181, 30))
plt.yticks(np.arange(-180, 181, 30))
plt.legend(title='Secondary Structure')
plt.grid(True)
plt.show()