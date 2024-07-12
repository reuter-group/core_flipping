from rdkit import Chem
from rdkit.Chem import AllChem
import sys

mol2 = Chem.MolFromMol2File(sys.argv[1], removeHs=False)
molecule = Chem.MolToSmiles(mol2)

molecule_smiles=Chem.MolFromSmiles(molecule)
chiralty=Chem.FindMolChiralCenters(molecule_smiles)
centers=len(chiralty)

for i in range(centers):
         current_center=chiralty[i][1]
         print(current_center)
         center_to_change=chiralty[i][0]
         print(center_to_change)


for a in molecule_smiles.GetAtoms():
     print(a.GetChiralTag())
