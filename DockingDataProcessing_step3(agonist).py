from pathlib import Path
import pandas as pd
import numpy as np
import argparse
import joblib
from pathlib import Path

from sklearn.metrics import confusion_matrix

from sklearn.base import BaseEstimator, TransformerMixin

resultsummary_folder = Path(f'{result_path}')

reslong = pd.read_csv('D:/ComData_PFAS/pdb_collect/3ert_1ere/res_chi1_group.csv')
reslist = reslong['resname'].unique()
interactionlist = ['Hpho','Hbond','vdw','pi']

class ParetoScaler(BaseEstimator, TransformerMixin):
    def fit(self, X, y=None):
        # Compute column-wise mean and sqrt(std)
        self.mean_ = np.mean(X, axis=0)
        self.scale_ = np.sqrt(np.std(X, axis=0, ddof=0))
        return self
    
    def transform(self, X):
        # Apply Pareto scaling
        return ((X - self.mean_) / self.scale_)
    
    def inverse_transform(self, X_scaled):
        # Reconstruct original values
        return (X_scaled * self.scale_ + self.mean_)
    
def contains_number(dataset):
    return any(char.isdigit() for char in dataset)

def balanced_classification_rate(y_true, y_pred):
    """
    Compute Balanced Classification Rate (BCR) using:
        BCR = (SE + SP) * (1 - |SE - SP|) / 2

    Parameters
    ----------
    y_true : array-like
        True binary labels.
    y_pred : array-like
        Predicted binary labels.

    Returns
    -------
    float
        Balanced Classification Rate (BCR). Returns np.nan if undefined.
    """
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()

    se = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    sp = tn / (tn + fp) if (tn + fp) > 0 else np.nan

    if np.isnan(se) or np.isnan(sp):
        return np.nan

    return (se + sp) * (1 - abs(se - sp)) / 2

def extract_allpose(proteinlist,activitylabel,resultsummary_folder1,dataset='EDSP'):
    deltaGlabel = pd.read_csv(resultsummary_folder1/dataset/f'{dataset}_combineddeltaG_Tc.csv',dtype={'CASRN':str})
    ressummary = pd.DataFrame()
    for interaction in interactionlist:
        interactionsum = pd.DataFrame()
        res_pdbis_cache = reslong['pdbid'].unique()
        res_pdbis_cache = [x.lower() for x in res_pdbis_cache]
        res_count = pd.DataFrame()
        for pro in proteinlist:
            closecontact = pd.read_csv(resultsummary_folder1/dataset/f'{pro}_noteq-{dataset}_{interaction}.csv',dtype={'CASN':str})
            closecontact.loc[:,'resin'] = closecontact['resin'].apply(lambda x: x.split('.')[0])
            h524index = closecontact[closecontact['resid']==524].index
            closecontact.loc[h524index,'resin'] = 'HIE524' #correct H524 to HIE524 regardless the ionization state
            closecontact.rename(columns={'CASN':'CASRN'},inplace=True)       
            deltaGlabelsub = pd.read_csv(resultsummary_folder1/dataset/f'{pro}_noteq-{dataset}_deltaG.csv',dtype={'CASRN':str})
            deltaGlabelsub.rename(columns={f'PoseID_{pro}_noteq':'pose'},inplace=True)  
            combined = deltaGlabelsub[['CASRN','pose',f'deltaG_{pro}_noteq']].merge(closecontact,how='left',on=['CASRN','pose']).dropna()
            if activitylabel in deltaGlabel.columns:
                deltaGlabeltemp = deltaGlabel[[activitylabel,'CASRN']] 
                combined = combined.merge(deltaGlabeltemp,how='left',on='CASRN').dropna()
            else:
                combined[activitylabel] = -1
            combined = combined.rename(columns={f'deltaG_{pro}_noteq':'deltaG'})
            combined = combined.merge(deltaGlabelsub.groupby('CASRN',as_index=False).agg(totalpose=('pose','count')),how='left',on='CASRN').dropna()
            combined = combined.drop_duplicates(subset=['CASRN','resid'])
            interactionsum = pd.concat([interactionsum,combined],axis=0)
        interactionsum['interaction'] = interaction
        ressummary = pd.concat([ressummary,interactionsum],axis=0) 
    ressummary['res_interact'] = ressummary['resin'] + ressummary['interaction']
    groupedsum = ressummary.groupby(by=['CASRN','res_interact',activitylabel,'totalpose'],as_index=False).agg(count_target=('protein_target','count'))
    groupedsum['count_target'] = groupedsum['count_target']/(len(proteinlist) * groupedsum['totalpose'])*100
    oplsdaready = groupedsum.pivot_table(index=['CASRN',activitylabel],values='count_target',columns='res_interact').reset_index().fillna(0).merge(ressummary.groupby(by=['CASRN'],as_index=False).agg(mean_deltaG=('deltaG','mean'),min_deltaG=('deltaG','min')),how='left',on='CASRN')    
    return oplsdaready



def main():
    parser = argparse.ArgumentParser(
                    prog='docking data processing workflow-step3(agonist)',
                    description='''
After step 3(agonist), 
agonist prediction will be generated
''',
                    epilog='Text at the bottom of help')    

    parser.add_argument('--datasetlist',default='EDSP')
    parser.add_argument('--progene',default='esr1')
    args = parser.parse_args()

    datasetlist = args.datasetlist  
    progene = args.progene

    datasetlist = datasetlist.split(',')
    resultsummary_folder1 = resultsummary_folder/progene

    proteinlist = ['3uucb','1qkua']
    activityColumn = input('please type in the column name for activity label')

    bestcolumn_ref = ['deltaG_3uucb_noteq','deltaG_1qkua_noteq','CASRN','PoseID'] #,'activity',f'{predicttype}_class',f'{predicttype}_class'

    final_model_ago = joblib.load(f'agonist{('').join(proteinlist)}.joblib') #agonist{('').join(proteinlist)} 'SBSagonist_classallposeBCR', SBSagonist_class3uucb1qkuc1gwqb1l2iaAP80
    selected_features = list(final_model_ago.named_steps['sfs'].k_feature_idx_)
    sfs = final_model_ago.named_steps['sfs']
    viporder = sfs.feature_names
    selected_feature_names = [viporder[i] for i in selected_features]    

    ad_checker = joblib.load('hotellingT2.joblib')
    for dataset in datasetlist:
        #score-based method
        df_all_combined_ref =pd.read_csv(resultsummary_folder1/dataset/f'{dataset}_combineddeltaG_Tc.csv')
        scoredf = df_all_combined_ref[bestcolumn_ref]
        scoredf['agVina1-score'] = (scoredf['deltaG_3uucb_noteq'] - scoredf['deltaG_1qkua_noteq'])
        se75_threshold = 0.1794749999999999
        scoredf['agVina1-score_pred'] = scoredf['agVina1-score'].apply(lambda x: 1 if x > se75_threshold else 0)
        scoredf = scoredf[['agVina1-score','agVina1-score_pred','CASRN',activityColumn,'PoseID']]

        #ML-based method
        oplsdaready = extract_allpose(proteinlist,activityColumn,resultsummary_folder1,dataset)
        sampleid = oplsdaready.pop('CASRN')    
        spectra = oplsdaready.copy(deep=True)
        target = spectra.pop(activityColumn)
        target = target.values  
        missingcollist = []
        for col in viporder:
            if col not in spectra.columns:
                missingcollist.append(col)  
        new_columns = {col: 0 for col in missingcollist}    
        spectra = pd.concat([
            spectra, 
            pd.DataFrame(new_columns, index=spectra.index)
        ], axis=1)
        X = spectra[viporder]
        y_pred = final_model_ago.predict(X)
        X_pred_selected = X[selected_feature_names]
        X_prednp = X_pred_selected.to_numpy()
        ad_pred = ad_checker.predict(X_prednp)
        tempdf = pd.concat([sampleid.loc[X.index],pd.DataFrame({f'agVina1-ML_AD':ad_pred},index=X.index),pd.DataFrame({f'agVina1-ML':y_pred},index=X.index)],axis=1)
        
        ago_predict_full = tempdf.merge(scoredf,how='inner',on='CASRN')
        ago_predict_full['agVina1-hybrid'] = ago_predict_full[['agVina1-score_pred','agVina1-ML','agVina1-ML_AD']].apply(lambda x: x['agVina1-ML'] if x['agVina1-ML_AD'] else x['agVina1-score_pred'],axis=1)

        ago_predict_full.to_csv(resultsummary_folder1/dataset/f'{dataset}_agonist_predition.csv',index=False)


