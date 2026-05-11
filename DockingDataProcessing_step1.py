from pathlib import Path
import pandas as pd
import numpy as np
import argparse
import os
import pdbqt2sdf
import MDAnalysis as mda
import prolif as plf
import oddt
from oddt.fingerprints import PLEC,tanimoto,SPLIF,similarity_SPLIF
from rdkit import Chem
from rdkit import DataStructs

#ligand_folder: folder contains all ligand pdbqt file used as input
#protein_folder: folder contains all protein pdbqt and pdb(with explicit hydrogen and bond information) file used as input
#result_path: folder store output results including extracted IFP info.
#rawresult_path: folder contains all conf files and output pdbqt 

ligand_folder = Path(f'{ligand_folder}')
protein_folder = Path(f'{protein_folder}')
resultsummary_folder = Path(f'{result_path}')
rawresult_folder = Path(f'{rawresult_path}')

antprolist = ['1ERRA', '1ERRB', '1R5KA', '1R5KB', '1R5KC', '2IOKA', '2IOKB', '2OUZA', '3DT3A', '3DT3B', '3ERTA', '3UUCA', '3UUCB', '3UUCC', '3UUCD', '1SJ0A', '1XP1A', '1XP6A', '1XP9A', '1XPCA', '1YIMA', '1YINA', '2IOGA']
fpdf_path = {'prolif':protein_folder/'ifpfile/ref_ifp_prolif.csv',
             'plec':protein_folder/'ifpfile/ref_ifp_plec_16384.csv',
             'splif':protein_folder/'ifpfile/ref_ifp_splif_4096.pkl'}

def findallpose(resultfolder,pro,dataset):
    resultpath = resultfolder/'{}-{}_deltaG.csv'.format(pro,dataset)
    pathtarget = rawresult_folder/'{}-{}/log'.format(pro,dataset)
    agonistres = sorted(pathtarget.glob("*_{}_log.txt".format(pro)))    
    if len(agonistres) == 0:
        print("result for {}-{} is missing".format(pro,dataset))
        return None,None
    IDlist=[]
    posenumber = []
    resultdeltaG = []
    for path in agonistres:
        ligandcas = path.name.split('.')[0].rsplit(f'_{pro}')[0]
        starttorecord = 0
        with open(path,'r') as fp:
            lines = fp.readlines()
            for row in lines:
                if row.startswith('   1'):
                    starttorecord = 1
                if starttorecord == 1:
                    IDlist.append(ligandcas)
                    test_list = row.split(' ')
                    result = [i for i in test_list if i]
                    posenumber.append(result[0]) 
                    resultdeltaG.append(result[1])

        fp.close()   
    IDlist = pd.DataFrame(IDlist,columns=['CASRN'])
    posenumber = pd.DataFrame(posenumber,columns=['PoseID_{}'.format(pro)])
    resultdeltaG = pd.DataFrame(resultdeltaG,columns=['deltaG_{}'.format(pro)])
    adtscoring = pd.concat([IDlist,posenumber,resultdeltaG],axis=1)  
    adtscoring['deltaG_{}'.format(pro)] = adtscoring['deltaG_{}'.format(pro)].astype(float)
    adtscoring['PoseID_{}'.format(pro)] = adtscoring['PoseID_{}'.format(pro)].astype(int)

    adtscoring.to_csv(resultpath,index=False)

    return adtscoring

def pdbqtall2ind(path,folderpath):
    with open(path, "r") as f1:
        pdbqt_lines = f1.readlines()
        newfilename = path.name.split('.')[0] + '_model'+str(1) + '.pdbqt'
        tempfile = []
        newfilepath = folderpath / newfilename
        for line in pdbqt_lines:
            if line.startswith('MODEL'):
                PoseID = line.split(" ")[1].strip()
                if int(PoseID) > 1:
                    with open(newfilepath, 'w') as f2:
                        f2.writelines(tempfile)
                    f2.close()
                newfilename = path.name.split('.')[0] + '_model'+PoseID+ '.pdbqt'
                newfilepath = folderpath / newfilename
                tempfile = []
            tempfile.append(line)
    f1.close()

def splitpdbqt(pro,dataset):
    pair = '{}-{}'.format(pro,dataset)
    folderpath_out = rawresult_folder /pair
    folderpath_in = folderpath_out/'out_pdbqt'
    folderpath = folderpath_in/'out_pdbqt'
    ligandpath = sorted(folderpath.glob(r"*_out.pdbqt".format(folderpath)))

    for path in ligandpath: #range(len(ligandprotein))
        pdbqtall2ind(path,folderpath_in)

def getsmile_ligandlist(ligandgroup):
    templatepath = ligand_folder / f'{ligandgroup}.csv'
    ligandtemplate = pd.read_csv(templatepath) 
    #Ligand ID must have column name as CASRN, and smiles column must have a column name as SMILES
    smilelist = ligandtemplate['SMILES'].to_list()
    ligandlist = ligandtemplate['CASRN'].to_list()

    ligandlist = [x.split(',')[0] if type(x) == str else x for x in ligandlist ]
    smilelist = [x.split(',')[0] for x in smilelist]
    return ligandlist,smilelist

def pdbqt2sdf_func(pro,dataset):
    pair = '{}-{}'.format(pro,dataset)
    folderpath_out = rawresult_folder /pair
    pdbqtpath = folderpath_out / pair /'out_pdbqt'
    sdfpath = folderpath_out / pair / 'out_sdf'
    pdbpath = folderpath_out / pair / 'out_pdb'
    if not os.path.exists(sdfpath):
        os.mkdir(sdfpath)
    if not os.path.exists(pdbpath):
        os.mkdir(pdbpath)

    ligandlist,smilelist = getsmile_ligandlist(dataset)

    for ligand, smile in zip(ligandlist,smilelist):
        ligand = str(ligand)
        temp = sorted(pdbqtpath.glob("{}*out_model*.pdbqt".format(ligand)))
        if len(temp) == 0:
            print("{} data is missing".format(ligand))
            continue
        input_filename_list = temp
        lig_smi = smile
        method = "slow"
        #pdbqt2sdf must run on virtual env for meeko
        with open(input_filename_list[0],"r") as file:
            lines = file.readlines()
            for line in lines:
                if "REMARK SMILES" in line:
                    method = "fast"
                    break
        file.close()
        if method == "fast":
            try:
                for input_filename in input_filename_list:                                    
                    pdbqt2sdf.pdbqt2sdf_fast(input_filename)
            except:
                try:
                    lig_rt = pdbqt2sdf.smile2template(lig_smi)
                except:
                    print('{} has error in smile template'.format(lig_smi))
                for input_filename in input_filename_list:
                    input_filename = str(input_filename)
                    try:
                        pdbqt2sdf.pdbqt2sdf(lig_smi,lig_rt,input_filename)
                    except Exception as e:
                        print('{} has error in pdbqt'.format(input_filename))
                        print(e)
                        break                
        else:
            try:
                lig_rt = pdbqt2sdf.smile2template(lig_smi)
            except:
                print('{} has error in smile template'.format(lig_smi))
            for input_filename in input_filename_list:
                input_filename = str(input_filename)
                try:
                    pdbqt2sdf.pdbqt2sdf(lig_smi,lig_rt,input_filename)
                except Exception as e:
                    print('{} has error in pdbqt'.format(input_filename))
                    print(e)
                    break

def combinesdf(sdf_files,ligand,pair):
    #sdf_files = sorted(pathlib.Path(folderpath).glob("{}*_out_model*.sdf".format(ligandname)))
    pairname_temp = ligand + '_' + pair.split('-')[0]
    combinedfile = "{}_out_model.sdf".format(pairname_temp)
    combinesdf_folder = sdf_files[0].parent/'out_sdf'
    filecount = len(sdf_files)
    if filecount > 20:
        filecount = 20
    #filenamepattern = pairname_temp+'_out_model'
    if not os.path.exists(combinesdf_folder):
        os.mkdir(combinesdf_folder)
    sdf_combine = combinesdf_folder / combinedfile
    with open(sdf_combine, 'w') as outfile:
        for index in range(filecount):
            filename2find1 = pairname_temp
            filename2find2 = '_out_model'+str(index+1)+'.sdf'
            for sdf in sdf_files:
                if filename2find1 in sdf.name and filename2find2 in sdf.name:
                    file = sdf                    
                    with open(file,'r') as infile:
                        outfile.write(sdf.name+'\n')
                        infile.readline() # skip first line
                        for line in infile:
                            outfile.write(line)
                    infile.close()
    outfile.close()
    sdf_combine = str(sdf_combine)
    return sdf_combine

def read_dockedprotein2prolif(proteinname):
    protein_eqpdb_path =protein_folder/"{}.pdb".format(proteinname)
    upro = mda.Universe(protein_eqpdb_path,guess_bonds=False)
    protein_mol = plf.Molecule.from_mda(upro)   
    return protein_mol

def read_dockfile2prolif(pair,ligand,smile):
    pdbqtpath = rawresult_folder / pair /'out_pdbqt'
    sdfpath = rawresult_folder / pair / 'out_sdf'
    if not os.path.exists(pdbqtpath):
        return None
    template = Chem.MolFromSmiles(smile)
    sdf_files = sorted(sdfpath.glob("{}*out_model*.sdf".format(ligand)))
    sdf_files = [i for i in sdf_files if i.name.split('_')[0] == ligand]
    if len(sdf_files) == 0:
        #seems fail to convert into sdf with Meeko, so we will use pdbqt files plus smiles
        pdbqt_files = sorted(pdbqtpath.glob("{}*out_model*.pdbqt".format(ligand)))
        pose_iterable = plf.pdbqt_supplier(pdbqt_files, template)
    else:
        sdf_combine = combinesdf(sdf_files,ligand,pair)
        #print(sdf_combine)
        pose_iterable = plf.sdf_supplier(sdf_combine)
    return pose_iterable

def tcagonist_rdkit(x,bitvectors_ref,fpdf_prolif_name):
    tclist = []
    tcname = []
    refbitv = []

    for i in range(len(bitvectors_ref)):
        if fpdf_prolif_name[i] not in antprolist:
            refbitv.append(bitvectors_ref[i])           
            tcname.append(fpdf_prolif_name[i])
    tclist = DataStructs.BulkTanimotoSimilarity(x, refbitv)
    maxtc = max(tclist)
    maxname = tcname[np.argmax(tclist)]
    return maxname,maxtc

def tcantagonist_rdkit(x,bitvectors_ref,fpdf_prolif_name):
    tclist = []
    tcname = []
    refindex = []
    for i in range(len(bitvectors_ref)):
        if fpdf_prolif_name[i] in antprolist:
            refindex.append(bitvectors_ref[i])
            tcname.append(fpdf_prolif_name[i])
    tclist = DataStructs.BulkTanimotoSimilarity(x, refindex)
    maxtc = max(tclist)
    maxname = tcname[np.argmax(tclist)]
    return maxname,maxtc

def calculate_interaction_fingerprint_ADT(pair,ligand,smile,protein_mol,threshold=4.5):
    pose_iterable = read_dockfile2prolif(pair,ligand,smile)
    if pose_iterable is None:
        print("fails to read {}".format(ligand))
        return None
    try:
        fp = plf.Fingerprint(parameters={"Hydrophobic": {"distance": threshold},"VdWContact":{"vdwradii":{"C":threshold}}})
        fp.run_from_iterable(pose_iterable, protein_mol,n_jobs=1)
    except Exception as e:
        print(e)
    return fp

def renamecolumn(x):
    residue = x[1].split('.')[0]
    reaction = x[2]
    return residue+reaction

def findbest_prolif(fp,pdbid):
    pdbid = pdbid.split('_')[0].upper()
    df_lig = fp.to_dataframe(index_col="Pose")
    df_lig.columns = df_lig.columns.to_flat_index()
    df_lig = df_lig.rename(renamecolumn,axis=1)   
    bitvectors_sample = plf.to_bitvectors(df_lig) 
    ligdict = {}
    for index in range(len(df_lig)):
        moldict = {index+1:bitvectors_sample[index].ToList()}
        ligdict.update(moldict)
    fpdflig = pd.DataFrame.from_dict(ligdict,orient='index',columns=df_lig.columns)
    fpdflig = fpdflig.reset_index(drop=False)
    fpdf_prolif = pd.read_csv(fpdf_path['prolif'])
    summary = pd.concat([fpdf_prolif,fpdflig],axis=0).fillna(int(0)).reset_index(drop=True)
    sample = summary.loc[len(fpdf_prolif):,:]
    reflist = summary.loc[:len(fpdf_prolif)-1,:]
    fpdf_prolif_name = reflist.pop('index').to_list()
    poselist = sample.pop('index').to_list()
    bitvectors_sample = plf.to_bitvectors(sample)
    bitvectors_ref = plf.to_bitvectors(reflist)
    matchpdb = []
    besttc = []
    if pdbid in antprolist:
        for pose in  range(len(bitvectors_sample)):
            maxname,maxtc = tcantagonist_rdkit(bitvectors_sample[pose],bitvectors_ref,fpdf_prolif_name)
            matchpdb.append(maxname)
            besttc.append(maxtc)          
    else:
        for pose in  range(len(bitvectors_sample)):
            maxname,maxtc = tcagonist_rdkit(bitvectors_sample[pose],bitvectors_ref,fpdf_prolif_name)
            matchpdb.append(maxname)
            besttc.append(maxtc)             
    result = pd.DataFrame({"poseID":poselist,"MatchPDB":matchpdb,"bestTc":besttc})
    bestpose = result['poseID'][result['bestTc'].argmax()]
    matchpdb = result['MatchPDB'][result['bestTc'].argmax()]
    return bestpose,matchpdb,result['bestTc'].max(),result

def _get_records(ifp, all_metadata):
    records = []
    for (lig_resid, prot_resid), int_data in ifp.items():
        for int_name, metadata_tuple in int_data.items():
            entry = {
                "ligand": str(lig_resid),
                "protein": str(prot_resid),
                "interaction": int_name,
            }
            if all_metadata:
                for metadata in metadata_tuple:
                    records.append(
                        {
                            **entry,
                            "atoms": metadata["parent_indices"]["ligand"],
                            "distance": metadata.get("distance", 0),
                        },
                    )
            else:
                # extract interaction with shortest distance
                metadata = min(
                    metadata_tuple,
                    key=lambda m: m.get("distance", np.nan),
                )
                entry["atoms"] = metadata["parent_indices"]["ligand"]
                entry["distance"] = metadata.get("distance", 0)
                records.append(entry)
    return records

def extract_res_distance(fp,dockeddata = False):
    interactionlist = plf.Fingerprint.list_available()
    temp_data = []
    temp_df = pd.DataFrame()
    for ifp in fp.ifp.values():
        temp_data.extend(_get_records(ifp, all_metadata=False))
    temp_df = pd.DataFrame(temp_data)
    if dockeddata:
        posenum = []
        count=0
        for pose in range(len(fp.ifp)):
            for i in interactionlist:
                count += sum(1 for interaction  in fp.ifp[pose].values() if i in interaction )
            for i in range(count):
                posenum.append(pose+1)
            count=0        
        posenum = pd.DataFrame(posenum,columns=['pose'])
        temp_df = pd.concat([posenum,temp_df],axis=1)
    return temp_df

def extract_fp_df(temp_df,interactiontype,dockeddata = False):
    interact_resin = temp_df[temp_df['interaction']==interactiontype]['protein'].values
    interact_dist = temp_df[temp_df['interaction']==interactiontype]['distance'].values
    interact_resin = pd.DataFrame(interact_resin,columns=['resin'])
    interact_dist = pd.DataFrame(interact_dist,columns=['distance'])
    interactdf = pd.concat([interact_resin,interact_dist],axis=1)
    interactdf['distance'] = interactdf['distance'].astype(float)
    interactdf.loc[:,'resid'] = interactdf['resin'].str.extract(r'(\d+)')
    interactdf.dropna(subset=['resin','distance'],how='all',inplace=True)
    interactdf['resid'] = interactdf['resid'].astype(int)
    if dockeddata:
        interact_pose = temp_df[temp_df['interaction']==interactiontype]['pose'].values
        interact_pose = pd.DataFrame(interact_pose,columns=['pose'])
        interact_pose['pose'] = interact_pose['pose'].astype(int)
        interactdf = pd.concat([interactdf,interact_pose],axis=1)
    return interactdf

def save_fp_df_adv(interactdf,interactiontype,dataset,ligand,pair,progene):
    resultpath = resultsummary_folder/progene/dataset
    interesttype = ['vdw','Hpho','pi','Hbond','CloseContact']
    interactdf.loc[:,'CASN'] = ligand
    interactdf.loc[:,'protein_target'] = pair.split('-')[0]
    if interactiontype == 'VdWContact':
        interactdf.to_csv(resultpath/'{}_{}.csv'.format(pair,interesttype[0]),mode='a',index = False)
    elif interactiontype == 'Hydrophobic':
        interactdf.to_csv(resultpath/'{}_{}.csv'.format(pair,interesttype[1]),mode='a',index = False)
    elif interactiontype == 'PiStacking':
        interactdf.to_csv(resultpath/'{}_{}.csv'.format(pair,interesttype[2]),mode='a',index = False)                    
    elif interactiontype in ['HBAcceptor','HBDonor']:
        interactdf.to_csv(resultpath/'{}_{}.csv'.format(pair,interesttype[3]),mode='a',index = False)
    elif interactiontype == 'CloseContact':
        interactdf.to_csv(resultpath/'{}_{}.csv'.format(pair,'CloseContact'),mode='a',index = False)

def clean_df_csv(resultpath,filename):
    interesttype = ['vdw','Hpho','pi','Hbond','CloseContact','prolifpick','prolifpickplus','plecpick','splifpick']
    csvsearch = [resultpath/'{}_{}.csv'.format(filename,type) for type in interesttype]
    for file in csvsearch:
        try:
            df = pd.read_csv(file,low_memory=False)
            df.drop_duplicates(inplace=True,ignore_index=True)
            columnlist = df.columns.tolist()
            df = df[df[columnlist[0]]!=columnlist[0]]
            df.to_csv(file,index = False)
        except:
            continue

def batch_fp_extract_dock(pro,dataset,resultpath,progene,processdata=True):
    ligandlist,smilelist = getsmile_ligandlist(dataset)
    interesttype = ['vdw','Hpho','pi','Hbond']
    protein_mol = read_dockedprotein2prolif(pro)
    pair = '{}-{}'.format(pro,dataset)
    if not os.path.exists(resultpath):
        os.mkdir(resultpath)
    for type in interesttype:
        csvpath = resultpath/'{}_{}.csv'.format(pair,type)
        if not os.path.exists(csvpath):
            with open(csvpath,'w') as df:
                pass
            df.close()

    #ligandsublist = check_fp_df_adv(ligandgroup,pair)
    prolif_pick_pose = []
    prolif_match_pdb = []
    prolif_best_tc = []
    prolif_ligand = []
    prolif_pickdf = pd.DataFrame()
    detailifp = pd.DataFrame()

    for ligand,smile in zip(ligandlist,smilelist):
        ligand = str(ligand)
        try:
            fp = calculate_interaction_fingerprint_ADT(pair,ligand,smile,protein_mol,threshold=4.5)
            test=fp.ifp
            if processdata:
                bestpose,matchpdb,besttc,tempdf = findbest_prolif(fp,pro)
                tempdf['CASRN'] = ligand
                prolif_pickdf = pd.concat([prolif_pickdf,tempdf],axis=0,ignore_index=True)
                prolif_pick_pose.append(bestpose)
                prolif_match_pdb.append(matchpdb)
                prolif_best_tc.append(besttc)
                prolif_ligand.append(ligand)

                temp_df = pd.DataFrame()
                temp_df = extract_res_distance(fp,dockeddata = True)
                intertypes = temp_df['interaction'].unique()

                for interactiontype in intertypes:
                    interactdf = extract_fp_df(temp_df,interactiontype,dockeddata = True)
                    save_fp_df_adv(interactdf,interactiontype,dataset,ligand,pair,progene)
            else:
                return fp                                        
        except:
            print("{} fails".format(ligand))
            continue

            #print(temp_df)              

    csvpathpro = resultpath/'{}_fullprolififp.csv'.format(pair)        #print(temp_df)
    detailifp.to_csv(csvpathpro,index=False)

    prolif_path = resultpath/'{}_prolifpick.csv'.format(pair)
    prolif_df = pd.DataFrame({"CASN":prolif_ligand,"Prolif_pose":prolif_pick_pose,
                            "Match_pdb":prolif_match_pdb,"Best_Tc":prolif_best_tc})
    prolif_p_path = resultpath/'{}_prolifpluspick.csv'.format(pair)

    prolif_df.to_csv(prolif_path,index=False)
    prolifdetail_path = resultpath/'{}_prolifdetail.csv'.format(pair)   
    prolif_pickdf.to_csv(prolifdetail_path,index=False)

    clean_df_csv(resultpath,pair) 
    return None      

def tcantagonist_oddt(x,ref_name,fpdf_plec_list):
    queryifp = np.array(x)
    tclist = []
    tcname = []
    for key,value in zip(ref_name,fpdf_plec_list):
        if key in antprolist:
            tclist.append(tanimoto(queryifp, value))
            tcname.append(key)
    maxtc = max(tclist)
    maxname = tcname[np.argmax(tclist)]
    return [maxname,maxtc]

def tcagonist_oddt(x,ref_name,fpdf_plec_list):
    queryifp = np.array(x)
    tclist = []
    tcname = []
    for key,value in zip(ref_name,fpdf_plec_list):
        if key not in antprolist:
            tclist.append(tanimoto(queryifp, value))
            tcname.append(key)
    maxtc = max(tclist)
    maxname = tcname[np.argmax(tclist)]
    return [maxname,maxtc]

def findbest_plec(moldata,receptor,matchfunc,ref_name,fpdf_ref_list):
    moldata['plec'] = moldata['mol'].map(lambda x: list(PLEC(x, protein=receptor, 
                                                    size=16384, 
                                                    depth_protein=5,
                                                    depth_ligand=1,
                                                    distance_cutoff=4.5,
                                                    sparse=False
                                                    )))
    moldata[['MatchPDB','bestTc']]  = moldata.apply(lambda x: matchfunc(x['plec'],ref_name,fpdf_ref_list),axis=1,result_type='expand')
    moldata['PoseID'] = moldata.index+1
    result = moldata[['PoseID','MatchPDB','bestTc']]
    bestpose = moldata['bestTc'].argmax()+1
    matchpdb = moldata['MatchPDB'][moldata['bestTc'].argmax()]  
    besttc = moldata['bestTc'].max()
    return bestpose,matchpdb,besttc,result

def findbest_oddt(pair,resultpath):

    proteinname = pair.split('-')[0]
    sdfpath = rawresult_folder / pair / 'out_sdf'/ 'out_sdf'
    pdbid = proteinname.split('_')[0].upper()
    proteinpath = protein_folder/"{}.pdb".format(proteinname)
    sdffilelist = sorted(sdfpath.glob('*_out_model.sdf'))
    receptor = next(oddt.toolkit.readfile('pdb', str(proteinpath)))


    ifpmethod = 'plec'
    fpdf_plec = pd.read_csv(fpdf_path[ifpmethod])
    ref_name_plec = fpdf_plec.pop('index')
    fpdf_plec.at[0,ifpmethod] = 0
    fpdf_plec[ifpmethod] = fpdf_plec[ifpmethod].astype(object)
    fpdf_plec[ifpmethod] = fpdf_plec.apply(lambda x: np.array(x[:-1]),axis=1)    
    fpdf_plec_list = fpdf_plec[ifpmethod].to_list()
    fpdf_plec_list = [np.array(x) for x in fpdf_plec_list]
    pick_pose_plec = []
    match_pdb_plec = []
    best_tc_plec = []
    plec_ligand_plec = []  
    detail_plec = pd.DataFrame()    

    if pdbid in antprolist:
        func_oddt = tcantagonist_oddt
    else:
        func_oddt = tcagonist_oddt

    for sdf in sdffilelist:
        ligandcasn = sdf.name.split('_')[0]
        ligandmols = list(oddt.toolkit.readfile('sdf', str(sdf)))
        moldata = pd.DataFrame({"mol":ligandmols})

        try:
            bestpose,matchpdb,besttc,tempdf = findbest_plec(moldata,receptor,func_oddt,ref_name_plec,fpdf_plec_list)
            tempdf['CASRN'] = ligandcasn
            detail_plec = pd.concat([detail_plec,tempdf],axis=0)
            pick_pose_plec.append(bestpose)
            match_pdb_plec.append(matchpdb)
            best_tc_plec.append(besttc)
            plec_ligand_plec.append(ligandcasn)
        except:
            print("{} fail in plec".format(ligandcasn))

    plec_path = resultpath/'{}_plecpick.csv'.format(pair)
    plec_df = pd.DataFrame({"CASN":plec_ligand_plec,"Prolif_pose":pick_pose_plec,
                            "Match_pdb":match_pdb_plec,"Best_Tc":best_tc_plec})
    plec_df.to_csv(plec_path,index=False)         

    plecdetail_path = resultpath/'{}_plec_detail.csv'.format(pair)
    detail_plec.to_csv(plecdetail_path,index=False)

def main():
    parser = argparse.ArgumentParser(
                    prog='docking data processing workflow-step1',
                    description='''
One batch of docking defines as one dataset docking to one protein conformation
The corresponding folder should be named as {pro}-{dataset}
Inside folder, for each ligand_protein conformation docking, it should include
1 log file named {lig}_{pro}_log.txt
1 pdbqt output file named {lig}_{pro}_out.pdbqt
''',
                    epilog='Text at the bottom of help')    
    parser.add_argument('--proteinlist',default='1erea')
    parser.add_argument('--datasetlist',default='EDSP')
    parser.add_argument('--progene',default='esr1')
    args = parser.parse_args()
    proteinlist = args.proteinlist
    datasetlist = args.datasetlist  
    progene = args.progene
    proteinlist = proteinlist.split(',')
    datasetlist = datasetlist.split(',')
     
    for dataset in datasetlist:
        resultfolder = resultsummary_folder/progene/'{}'.format(dataset)
        for pro in proteinlist:
            proligpair = pro+'-'+dataset
            pro1erea_all = findallpose(resultfolder,pro,dataset)
            datafolder = rawresult_folder/proligpair

            #split each pose to individual pdbqt file
            splitpdbqt(pro,dataset)

            #convert pdbqt to sdf file
            pdbqt2sdf_func(pro,dataset)

            #extract protein-ligand fingerprint
            #extract using ProLIF
            batch_fp_extract_dock(pro,dataset,resultfolder,progene,processdata=True)
            #extract using ODDT
            findbest_oddt(proligpair,resultfolder)

if __name__ == '__main__':
    main()