from pathlib import Path
import pandas as pd
import argparse
import os
#ligand_folder: folder contains all ligand pdbqt file used as input
#result_path: folder store output results including extracted IFP info.
resultsummary_folder = Path(f'{result_path}')
ligand_folder = Path(f'{ligand_folder}')

def main():
    parser = argparse.ArgumentParser(
                    prog='docking data processing workflow-step2',
                    description='''
After step 2, 
a combined file will be created for each scoring combination (scoring method for pose selection x ranking)
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
    activityColumn = input('please type in the column name for activity label')
    resultsummary_folder1 = resultsummary_folder/progene
    #assign a sample ID to each pose
    for dataset in datasetlist:
        fp_cache = {}
        pose_cache = {}
        ifpdetail_cache = {}
        pairlist = []
        for pro in proteinlist:

            pair = pro+'-'+dataset
            key = f"{pair}"
            delGpath = resultsummary_folder1/dataset/'{}_deltaG.csv'.format(pair)    
            fp_cache[key] = pd.read_csv(delGpath,dtype={'CASRN':str,'PoseID_{}'.format(pro):str})
            fp_cache[key]['sampleID'] = fp_cache[key]['CASRN']  + 'p' + fp_cache[key]['PoseID_{}'.format(pro)]
            fp_cache[key].drop_duplicates(subset=['sampleID'], inplace=True,keep='last',ignore_index=True) 
            fp_cache[key].to_csv(delGpath, index=False)
            if len(fp_cache[key]) > 0:
                for ifp in ['plec','prolif']:
                    pickposepath = resultsummary_folder1/dataset/'{}_{}pick.csv'.format(pair,ifp)
                    if ifp in ['plec']:
                        ifpdetailpath = resultsummary_folder1/dataset/'{}_{}_detail.csv'.format(pair,ifp)
                    else:
                        ifpdetailpath = resultsummary_folder1/dataset/'{}_{}detail.csv'.format(pair,ifp)
                    if os.path.exists(pickposepath):
                        pose_cache[key] = pd.read_csv(pickposepath,dtype={'CASN':str,'Prolif_pose':str})  
                        pose_cache[key]['sampleID'] = pose_cache[key]['CASN']  + 'p' + pose_cache[key]['Prolif_pose']  
                        pose_cache[key].drop_duplicates(subset=['CASN'], inplace=True,keep='last',ignore_index=True)
                        pose_cache[key].to_csv(pickposepath, index=False)
                    if os.path.exists(ifpdetailpath):
                        ifpdetail_cache[key] = pd.read_csv(ifpdetailpath,dtype={'CASRN':str,'PoseID':str,'poseID':str})
                        
                        try:
                            ifpdetail_cache[key]['sampleID'] = ifpdetail_cache[key]['CASRN']  + 'p' + ifpdetail_cache[key]['PoseID'] 
                            ifpdetail_cache[key].rename(columns={'bestTc':'Best_Tc'}, inplace=True)
                        except:
                            ifpdetail_cache[key]['sampleID'] = ifpdetail_cache[key]['CASRN']  + 'p' + ifpdetail_cache[key]['poseID'] 
                            ifpdetail_cache[key].rename(columns={'poseID':'PoseID','bestTc':'Best_Tc'}, inplace=True)
                        
                        
                        ifpdetail_cache[key].drop_duplicates(subset=['sampleID'], inplace=True,keep='last',ignore_index=True) 
                        
                        ifpdetail_cache[key].to_csv(ifpdetailpath, index=False)        

    for dataset in datasetlist:
        datadf_path = ligand_folder/dataset/f'{dataset}.csv'
        labeldf = pd.read_csv(datadf_path,dtype={'CASRN':str})
        for ifpmethod in ['prolif','delta','plec']:
            fp_cache = {}
            pose_cache = {}
            merge_cache = {}
            pairlist = []
            df_all = pd.DataFrame()
            df_all_combined = pd.DataFrame()
            if ifpmethod in ['plec','prolif']:
                for pro in proteinlist: 
                    pair = pro+'-'+dataset
                    delGpath = resultsummary_folder1/dataset/'{}_deltaG.csv'.format(pair)
                    pickposepath = resultsummary_folder1/dataset/'{}_{}pick.csv'.format(pair,ifpmethod)
                    if (os.path.isfile(delGpath)) & (os.path.isfile(pickposepath)):
                        
                        if pair not in pairlist:
                            pairlist.append(pair)
                            key = f"{pair}"
                            fp_cache[key] = pd.read_csv(delGpath,dtype={'CASRN':str,'PoseID_{}'.format(pro):str})
                            pose_cache[key] = pd.read_csv(pickposepath,dtype={'CASN':str,'Prolif_pose':str})
                            if len(fp_cache[key]) > 0:
                                fp_cache[key]['sampleID'] = fp_cache[key]['CASRN']  + 'p' + fp_cache[key]['PoseID_{}'.format(pro)]
                                pose_cache[key]['sampleID'] = pose_cache[key]['CASN']  + 'p' + pose_cache[key]['Prolif_pose']
                                merge_cache[key] = fp_cache[key].merge(pose_cache[key],how='inner',on = 'sampleID')
                                merge_cache[key] = merge_cache[key][['CASRN',"deltaG_{}".format(pro),'Best_Tc']]
                                merge_cache[key]['deltaGxBestTc_{}'.format(pro)] = merge_cache[key]["deltaG_{}".format(pro)] * merge_cache[key]["Best_Tc"]
                                merge_cache[key].reset_index(drop=True,inplace=True)
                                merge_cache[key].rename(columns={"Best_Tc":"Best_Tc_{}".format(pro)},inplace=True)
                    else:
                        continue
                    if len(df_all) == 0:
                        df_all = merge_cache[key].copy(deep=True)
                    else:
                        df_all = df_all.merge(merge_cache[key],how='inner',on='CASRN')
                

                subcollist = [activityColumn]
                subcollist.append('CASRN')
                df_all = df_all.merge(labeldf[subcollist],how='left',on='CASRN')

                df_all_combined = df_all.copy(deep=True)
            elif ifpmethod == 'delta':
                for pro in proteinlist: 
                    pair = pro+'-'+dataset
                    delGpath = resultsummary_folder1/dataset/'{}_deltaG.csv'.format(pair)
                    if os.path.isfile(delGpath):
                        pairlist.append(pair)
                        key = f"{pair}"
                        fp_cache[key] = pd.read_csv(delGpath,dtype={'CASRN':str})
                        if len(fp_cache[key]) > 0:
                            df1 = fp_cache[key][fp_cache[key]['PoseID_{}'.format(pro)]==1][['CASRN','deltaG_{}'.format(pro),'PoseID_{}'.format(pro)]]
                            #df1 = fp_cache[key][fp_cache[key]['PoseID']==1][['CASRN','deltaG']]
                            df1.rename(columns={'PoseID_{}'.format(pro):"PoseID"},inplace=True)
                            if len(df_all) == 0:
                                df_all = df1.copy(deep=True)
                            else:
                                df_all = df_all.merge(df1,how='outer',on=['CASRN','PoseID'])

                subcollist = [activityColumn]
                subcollist.append('CASRN')
                df_all = df_all.merge(labeldf[subcollist],how='left',on='CASRN')
                #print(len(df_all_combined))
                df_all_combined = df_all.copy(deep=True)
            #print(len(df_all_combined))
                df_all_combined = df_all_combined.drop_duplicates().reset_index(drop=True)
            resultoutpath = resultsummary_folder1/dataset/f'{dataset}_combined{ifpmethod}G_Tc.csv'
            if not df_all_combined.isnull().values.any():
                df_all_combined.to_csv(resultoutpath,index=False)
            else:
                print('There are NaN values in the combined dataframe.')

    for dataset in datasetlist:
        datadf_path = ligand_folder/dataset/f'{dataset}.csv'
        labeldf = pd.read_csv(datadf_path,dtype={'CASRN':str})
        for ifpmethod in ['prolif','plec']:
            fp_cache = {}
            pose_cache = {}
            merge_cache = {}
            picked_cache = {}
            pairlist = []
            df_all = pd.DataFrame()
            df_all_combined = pd.DataFrame()

            for pro in proteinlist: 
                pair = pro+'-'+dataset
                delGpath = resultsummary_folder1/dataset/'{}_deltaG.csv'.format(pair)
                if ifpmethod in ['plec']:
                    ifpdetailpath = resultsummary_folder1/dataset/'{}_{}_detail.csv'.format(pair,ifpmethod)
                else:
                    ifpdetailpath = resultsummary_folder1/dataset/'{}_{}detail.csv'.format(pair,ifpmethod)
                pickposepath = resultsummary_folder1/dataset/'{}_{}pick.csv'.format(pair,ifpmethod)
                if (os.path.isfile(delGpath)) & (os.path.isfile(ifpdetailpath)):
                    
                    if pair not in pairlist:
                        pairlist.append(pair)
                        key = f"{pair}"
                        fp_cache[key] = pd.read_csv(delGpath,dtype={'CASRN':str,'PoseID_{}'.format(pro):str})
                        pose_cache[key] = pd.read_csv(ifpdetailpath,dtype={'CASRN':str,'PoseID':str})
                        if len(fp_cache[key]) > 0:
                            merge_cache[key] = fp_cache[key].merge(pose_cache[key],how='inner',on = 'sampleID')                       
                            merge_cache[key]['deltaGxBestTc_{}'.format(pro)] = merge_cache[key]["deltaG_{}".format(pro)] * merge_cache[key]["Best_Tc"]
                            merge_cache[key].rename(columns={"CASRN_x":"CASRN"},inplace=True)
                            merge_cache[key] = merge_cache[key][['CASRN','deltaGxBestTc_{}'.format(pro),"deltaG_{}".format(pro),"Best_Tc"]]
                            merge_cache[key].reset_index(drop=True,inplace=True)
                            picked_cache[key] = merge_cache[key].groupby('CASRN',as_index=False).agg(deltaGxBestTc=('deltaGxBestTc_{}'.format(pro),'min'),deltaGxBestTc_index=('deltaGxBestTc_{}'.format(pro),'idxmin'))
                            picked_cache[key].rename(columns={"deltaGxBestTc":'deltaGxBestTc_{}'.format(pro),},inplace=True)
                            picked_cache[key]["deltaG_{}".format(pro)] = merge_cache[key].loc[picked_cache[key]['deltaGxBestTc_index'], 'deltaG_{}'.format(pro)].values
                            picked_cache[key]["Best_Tc_{}".format(pro)] = merge_cache[key].loc[picked_cache[key]['deltaGxBestTc_index'], "Best_Tc"].values
                            picked_cache[key].drop(columns=['deltaGxBestTc_index'], inplace=True)
                            picked_cache[key].reset_index(drop=True,inplace=True) #pick the min deltaGxBestTc
                            if len(df_all) == 0:
                                df_all = picked_cache[key].copy(deep=True)
                            else:
                                df_all = df_all.merge(picked_cache[key],how='inner',on='CASRN')
                else:
                    continue
                if len(df_all) == 0:
                    df_all = picked_cache[key].copy(deep=True)
                else:
                    df_all = df_all.merge(picked_cache[key],how='inner',on='CASRN')
            

            subcollist = [activityColumn]
            subcollist.append('CASRN')
            df_all = df_all.merge(labeldf[subcollist],how='left',on='CASRN')

            df_all_combined = df_all.copy(deep=True)

            df_all_combined = df_all_combined.drop_duplicates().reset_index(drop=True)
            resultoutpath = resultsummary_folder1/dataset/f'{dataset}_combined{ifpmethod}GtTc.csv'
            if not df_all_combined.isnull().values.any():
                df_all_combined.to_csv(resultoutpath,index=False)
            else:
                print('There are NaN values in the combined dataframe.')


if __name__ == '__main__':
    main()