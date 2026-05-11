from pathlib import Path
import pandas as pd
import argparse
from pathlib import Path
import pandas as pd


resultsummary_folder = Path(f'{result_path}')

reslong = pd.read_csv('res_chi1_group.csv')
reslist = reslong['resname'].unique()
interactionlist = ['Hpho','Hbond','vdw','pi']


def main():
    parser = argparse.ArgumentParser(
                    prog='docking data processing workflow-step3(antagonist)',
                    description='''
After step 3(antagonist), 
antagonist prediction will be generated
''',
                    epilog='Text at the bottom of help')    

    parser.add_argument('--datasetlist',default='EDSP')
    parser.add_argument('--progene',default='esr1')
    args = parser.parse_args()

    datasetlist = args.datasetlist  
    progene = args.progene

    datasetlist = datasetlist.split(',')
    resultsummary_folder1 = resultsummary_folder/progene

    activityColumn = input('please type in the column name for activity label')

    bestcolumn_ref = ['Best_Tc_1r5kc_noteq','Best_Tc_1gwqb_noteq','CASRN','PoseID'] #,'activity',f'{predicttype}_class',f'{predicttype}_class'

    for dataset in datasetlist:
        #score-based method
        df_all_combined_ref =pd.read_csv(resultsummary_folder1/dataset/f'{dataset}_combinedplecG_Tc.csv')
        scoredf = df_all_combined_ref[bestcolumn_ref]
        scoredf['anPLEC1-score'] = (scoredf['Best_Tc_1r5kc_noteq'] + scoredf['Best_Tc_1gwqb_noteq']) * 0.5
        se75_threshold = 0.0222061826179831
        scoredf['anPLEC1-score_pred'] = scoredf['anPLEC1-score'].apply(lambda x: 1 if x > se75_threshold else 0)
        scoredf = scoredf[['anPLEC1-score','anPLEC1-score_pred','CASRN',activityColumn]]

        scoredf.to_csv(resultsummary_folder1/dataset/f'{dataset}_antagonist_predition.csv',index=False)


