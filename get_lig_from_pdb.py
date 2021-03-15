import sys
from prody import *
from pdbe import pyPDBeREST
import json
from rdkit import Chem
from rdkit.Chem import AllChem
from io import StringIO
from collections import defaultdict
from rdkit.Chem import rdFMCS
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolAlign
IPythonConsole.ipython_3d = True
import py3Dmol
from openff.toolkit.topology import Molecule
from openff.toolkit.topology import Topology
from openff.toolkit.utils import RDKitToolkitWrapper
from openff.toolkit.typing.engines.smirnoff import ForceField
from simtk import openmm
from simtk import unit
from rdkit.Chem import rdMolAlign
import numpy as np
from simtk import openmm, unit
from datetime import timedelta
from timeit import time

def stopwatch(method):
    def timed(*args, **kw):
        ts = time.perf_counter()
        result = method(*args, **kw)
        te = time.perf_counter()
        duration = timedelta(seconds=te - ts)
        print(f"{method.__name__}: {duration}")
        return result
    return timed

@stopwatch
def drawit(m,p=None,confId=-1):
        mb = Chem.MolToMolBlock(m,confId=confId)
        if p is None:
            p = py3Dmol.view(width=400,height=400)
        p.removeAllModels()
        p.addModel(mb,'sdf')
        p.setStyle({'stick':{}})
        p.setBackgroundColor('0xeeeeee')
        p.zoomTo()
        return p.show()

@stopwatch
def get_bound_ligand_from_pdb(pdb_id):
    # get information associated with molecules from PDBe
    p = pyPDBeREST()
    molecule_info = p.PDB.getMolecules(pdbid=pdb_id)
    
    #get the pdb in a way similar to Pat's post
    ag = parsePDB(pdb_id)
    
    #load the json into python from PDBe
    compound_dict = json.loads(molecule_info)
    
    #use list comprehension to deal with the compounds in the json
    compounds = [compound_dict[pdb_id][l] for l in range(len(compound_dict[pdb_id]))]

    # create a dictionary with key pairs and use them to get ligand name and ligand_code for extraction
    # with prody. Loop of the list above and get the info
    Dict = {}
    Dict['molecule_name'] = {}
    for entries in compounds:    
        if entries['molecule_type'] == 'bound':
            ligand_name = entries['molecule_name'][0]
            ligand_code = entries['chem_comp_ids'][0]
            Dict['molecule_name'][ligand_name] = {'chem_comp_ids':ligand_code}
    ## use list comprehension to unwrap the compound names
    compound_names = [k for k in Dict['molecule_name']]
    
    ## create a new dictionary for RDKIT mol objects we might find.
    extracted_ligands = {}
    extracted_ligands['ligand_code'] = {}
    
    ##Now, lets loop over the compound_names list and use count/enumerate to unwrap 
    # ligand codes one by one and then pass them into our function
    for count,compound in enumerate(compound_names):
        ligand_code = Dict['molecule_name'][compound_names[count]]['chem_comp_ids']
        
        #get the info associated with the bound ligand and return the SMILES
        ligand_data = fetchPDBLigand(ligand_code)
        template = AllChem.MolFromSmiles(ligand_data['OpenEye_OEToolkits_SMILES_CANONICAL'])

        #now we need the prody command to pull the ligand
        arg = 'resname_{ligand_code}'.format(ligand_code = ligand_code)

        #use getattr to pass in arg the allow us to select the ligand in the pdb
        ligand = getattr(ag, arg)

        #use stringIO to collect the pdb information for the ligand
        output = StringIO()
        writePDBStream(output, ligand)
        pdb_string = output.getvalue()

        #make an RDKit mol
        rd_mol = AllChem.MolFromPDBBlock(pdb_string)
        
        #now we use our logic from before, but this time we put the rdkit mol objects
        #into a dictionary
        
        if len(AllChem.GetMolFrags(rd_mol)) == 1:
            mol_w_bond_orders = AllChem.AssignBondOrdersFromTemplate(template, rd_mol)
            mol = AllChem.AddHs(mol_w_bond_orders, addCoords=True)
            extracted_ligands['ligand_code'][ligand_code] = {'rdkmol_obj':mol}
        else:
            split_structures = Chem.GetMolFrags(rd_mol,asMols=True)
            if template.GetNumAtoms() == split_structures[0].GetNumAtoms():
                mol_w_bond_orders = AllChem.AssignBondOrdersFromTemplate(template, split_structures[0])
                mol = AllChem.AddHs(mol_w_bond_orders, addCoords=True)
                extracted_ligands['ligand_code'][ligand_code] = {'rdkmol_obj':mol}
            else:
                raise ValueError()
                print(f'template and identified mol object from {pdb_id} do not have the same number of atoms')

    return extracted_ligands

pdb_id = '7ckz'
extracted_ligands = get_bound_ligand_from_pdb(pdb_id)
rdkit_mol = [r for r in extracted_ligands['ligand_code']]
for count,mol in enumerate(rdkit_mol):
    mol_obj = extracted_ligands['ligand_code'][rdkit_mol[count]]['rdkmol_obj']
    print(rdkit_mol[count])
    print(mol_obj.GetNumAtoms())
    print(Chem.MolToMolBlock(mol_obj))

# @stopwatch
# def get_bound_ligand_from_pdb(pdb_id):
#     # get information associated with molecules from PDBe
#     p = pyPDBeREST()
#     molecule_info = p.PDB.getMolecules(pdbid=pdb_id)
    
#     #get the pdb in a way similar to Pat's post
#     ag = parsePDB(pdb_id)
    
#     #load the json into python from PDBe
#     compound_dict = json.loads(molecule_info)
    
#     #use list comprehension to deal with the compounds in the json
#     compounds = [compound_dict[pdb_id][l] for l in range(len(compound_dict[pdb_id]))]
    
#     #run a for loop and find the bound ligand
#     for entries in compounds:
#         if entries['molecule_type'] == 'bound':
#             ligand_code = entries['chem_comp_ids'][0]
#             break
#         else:
#             print('Not your ligand')
            
#     #get the info associated with the bound ligand and return the SMILES
#     ligand_data = fetchPDBLigand(ligand_code)
#     template = AllChem.MolFromSmiles(ligand_data['OpenEye_OEToolkits_SMILES_CANONICAL'])
    
#     #now we need the prody command to pull the ligand
#     arg = 'resname_{ligand_code}'.format(ligand_code = ligand_code)
    
#     #use getattr to pass in arg the allow us to select the ligand in the pdb
#     ligand = getattr(ag, arg)
    
#     #use stringIO to collect the pdb information for the ligand
#     output = StringIO()
#     writePDBStream(output, ligand)
#     pdb_string = output.getvalue()
    
#     #make an RDKit mol
#     rd_mol = AllChem.MolFromPDBBlock(pdb_string)

#     #now the weird part.. I wanted to start with 3pbl, and ETQ always came out as a dimer. I bet
#     #this can be made simplier if new more about how prody worked, but this was an ok solution
#     #for now. Basically we are making some logic up.. if we get mol objects as dimers, split them
#     #and make sure our smiles matches up. 
    
#     if len(AllChem.GetMolFrags(rd_mol)) == 1:
#         mol_w_bond_orders = AllChem.AssignBondOrdersFromTemplate(template, rd_mol)
#         mol = AllChem.AddHs(mol_w_bond_orders, addCoords=True)
#         return mol
#     else:
#         split_structures = Chem.GetMolFrags(rd_mol,asMols=True)
#         if template.GetNumAtoms() == split_structures[0].GetNumAtoms():
#             mol_w_bond_orders = AllChem.AssignBondOrdersFromTemplate(template, rs[0])
#             mol = AllChem.AddHs(mol_w_bond_orders, addCoords=True)
#             return mol
#         else:
#             raise ValueError()
#             print(f'template and identified mol object from {pdb_id} do not have the same number of atoms')

