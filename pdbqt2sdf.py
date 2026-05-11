from rdkit import Chem
from meeko import ResidueTemplate
import argparse
from meeko import PDBQTMolecule
from meeko import RDKitMolCreate

def cmd_lineparser():
    parser = argparse.ArgumentParser(
        description='Export docked ligand to SDF, and receptor to PDB',
    )
    parser.add_argument('-f', '--pdbqt_files',  nargs = "+",metavar='seperated pdbqt files', default='D:/ComData_PFAS/1erea-refER/out_pdbqt/CHEMBL132868_1erea_out_model1.pdbqt',
                        help="list of pdbqt files for the ligand")
    parser.add_argument(
        '-s',
        '--smile',default='CC12CCC3C(CCc4cc(O)ccc43)C1CCC2O',
        metavar='smile for the ligand',
        help="smile is used as template",
    )
    return parser.parse_args()

def fix_2943_75_1(mol_from_pdb):
    Oxygenindex = []
    for atom in mol_from_pdb.GetAtoms():
        if atom.GetSymbol() == 'O':
            Oxygenindex.append(atom.GetIdx())  
    emol = Chem.EditableMol(mol_from_pdb)
    emol.RemoveBond(Oxygenindex[0],Oxygenindex[1])
    emol.RemoveBond(Oxygenindex[0],Oxygenindex[2])
    emol.RemoveBond(Oxygenindex[1],Oxygenindex[2])
    new_mol_from_pdb = emol.GetMol()  
    return new_mol_from_pdb

def pdbqt2pdb(input_filename):

    with open(input_filename, "r") as f:
        pdbqt_lines = f.readlines()
    col_id = 77
    pdbqt_block = ""
#
    for line in pdbqt_lines: 
        if len(line) > col_id and line[col_id] == "A":
            # Replace A by C
            line = line[:col_id] + "C" + line[col_id + 1:]
        if line[col_id:col_id+2] in ['Cl','Br','Si']:
            line = line[:col_id-1] + line[col_id:]
        pdbqt_block += line
    try:
        mol_from_pdb = Chem.MolFromPDBBlock(pdbqt_block, removeHs = False)
        if "2943-75-1" in str(input_filename):
            mol_from_pdb = fix_2943_75_1(mol_from_pdb)
    except:
        mol_from_pdb = Chem.MolFromPDBBlock(pdbqt_block, sanitize=False,removeHs = False)
        if "2943-75-1" in str(input_filename):
            mol_from_pdb = fix_2943_75_1(mol_from_pdb)
    return mol_from_pdb

def smile2template(lig_smi):
    lig_rt = ResidueTemplate(lig_smi)
    return lig_rt

def mapindex(lig_rt,mol_from_pdb):
    mapping = lig_rt.match(input_mol = mol_from_pdb)
    zero_based = mapping[1]
    index_map = [val for k, v in zero_based.items() for val in (k + 1, v + 1)]
    h_parent = []
    zero_based_inv = {v: k for k,v in zero_based.items()}
    for atom in mol_from_pdb.GetAtoms():
        if atom.GetAtomicNum()==1:
            print(atom.GetSymbol())
            parent_atom = atom.GetNeighbors()[0]
            parent_index_in_smi = zero_based_inv[parent_atom.GetIdx()]
            h_parent.extend([parent_index_in_smi+1, atom.GetIdx()+1])    
    return index_map,h_parent

def pdbqt2sdf_fast(input_filename):
    input_filename = str(input_filename)
    pdbqt_mol = PDBQTMolecule.from_file(input_filename, skip_typing=True)
    rdkitmol_list = RDKitMolCreate.from_pdbqt_mol(
        pdbqt_mol,
        only_cluster_leads=False,
        keep_flexres=False,
    )
    sdf_fn =  input_filename.replace('pdbqt','sdf') 
    with Chem.SDWriter(sdf_fn) as writer:
        for mol in rdkitmol_list:
            writer.write(mol)

def pdbqt2sdf(lig_smi,lig_rt,input_filename):
    try:
        mol_from_pdb = pdbqt2pdb(input_filename)
        #print('sucess read as pdb')
    except Exception as e:
        print(e)
    index_map,h_parent = mapindex(lig_rt,mol_from_pdb)
    try:
        pdbqt_mol = PDBQTMolecule.from_file(input_filename, skip_typing=True)
        print('sucess read as pdbqt')
    except Exception as e:
        print(e)
    pdbqt_mol._pose_data['smiles'] = {0: lig_smi}
    pdbqt_mol._pose_data['smiles_index_map'] = {0: index_map}
    pdbqt_mol._pose_data['smiles_h_parent'] = {0: h_parent}  
    try:  
        rdkit_mol = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
        sdf_fn =  input_filename.replace('pdbqt','sdf') 
        with Chem.SDWriter(sdf_fn) as writer:
            for mol in rdkit_mol:
                writer.write(mol)
    except:
        sdf_fn =  input_filename.replace('pdbqt','sdf')
        print("RDKitMolCreate fails")
        template = Chem.MolFromSmiles(lig_smi)
        mol_from_pdb_nh = Chem.RemoveHs(mol_from_pdb)
        mol_from_pdb = Chem.AllChem.AssignBondOrdersFromTemplate(template,mol_from_pdb_nh)
        mol_from_pdbh = Chem.AddHs(mol_from_pdb, addCoords=True)
        writer = Chem.SDWriter(sdf_fn)
        writer.write(mol_from_pdbh)


    
def main():
    args = cmd_lineparser()
    input_filename_list = args.pdbqt_files
    input_filename_list = input_filename_list[0].split(' ')
    #print(len(input_filename_list))
    lig_smi = args.smile
    try:
        lig_rt = smile2template(lig_smi)
    except:
        print('{} has error in smile template'.format(lig_smi))
        raise
    for input_filename in input_filename_list:
        try:
            pdbqt2sdf(lig_smi,lig_rt,input_filename)
        except Exception as e:
            print('{} has error in pdbqt'.format(input_filename))
            print(e)
            break

if __name__ == "__main__":
    main()