import pandas as pd
import numpy as np
from Modules.PlotDataAnalysis import PlotDataAnalysis
import os
wd=os.getcwd()
bar='\\' if os.name=='nt' else '/'
result_path=wd+bar+'Results'
save_path=wd+bar+'DataAnalysis'
catalog_path=wd+bar+'Metadata'

plotter=PlotDataAnalysis(save_path,catalog_path,result_path)

df1=pd.read_csv(save_path+bar+'Final_results_simulations.dat',delim_whitespace=True)
df2=pd.read_csv(save_path+bar+'Metada_simulations_variables.dat',delim_whitespace=True)

final_df=pd.DataFrame([])
for simulation in df1.Simulation.values:
    ind=df2[df2.isin([simulation])].dropna(how='all').index
    first_half_row=df2.iloc[ind,:]
    if first_half_row.empty:
        continue
    second_half_row=df1[df1['Simulation']==simulation]
    first_half_row.reset_index(drop=True,inplace=True)
    second_half_row.reset_index(drop=True,inplace=True)
    df=pd.concat([first_half_row,second_half_row],axis=1)
    final_df=final_df.append(df)

final_df.reset_index(drop=True,inplace=True)
##Add new columns    
final_df['deltaJ']=np.sqrt(final_df['deltaJx']**2+ final_df['deltaJy']**2+final_df['deltaJz']**2)
final_df['J']=np.sqrt(final_df['Jx']**2+final_df['Jy']**2+final_df['Jz']**2)

final_df['deltaM']=-final_df['deltaM']
final_df['deltaM_Roch']=abs(final_df['final_mass']-final_df['initial_mass'])
final_df['P']=np.sqrt(final_df['Px']**2+ final_df['Py']**2+final_df['Pz']**2)
final_df['PdM']=final_df['P']/final_df['M']
df['Erad']=1-df['M']/df['initial_mass']

final_df['Jin']=np.sqrt(final_df['Jx_adm']**2+final_df['Jy_adm']**2+final_df['Jz_adm']**2)
final_df['deltaJ_Roch']=abs(final_df['Jf_Roch']-final_df['Jin'])
final_df['Sf']=np.sqrt(final_df['Sx']**2+final_df['Sy']**2+final_df['Sz']**2)
final_df['S1']=np.sqrt(final_df['s1x']**2+final_df['s1y']**2+final_df['s1z']**2)
final_df['S2']=np.sqrt(final_df['s2x']**2+final_df['s2y']**2+final_df['s2z']**2)

final_df['chi']=final_df['J']/final_df['M']**2
final_df['chiz']=final_df['Jz']/final_df['M']**2

final_df['q']=final_df['m1']/final_df['m2']

final_df['LinS1']=final_df['Lin']*final_df['s1z']/(np.sqrt(final_df['s1x']**2+final_df['s1y']**2+final_df['s1z']**2)*abs(final_df['Lin']))
final_df['LinS2']=final_df['Lin']*final_df['s2z']/(np.sqrt(final_df['s2x']**2+final_df['s2y']**2+final_df['s2z']**2)*abs(final_df['Lin']))
final_df['S1S2']=(final_df['s1x']*final_df['s2x']+final_df['s1y']*final_df['s2y']+final_df['s1z']*final_df['s2z'])/(np.sqrt(final_df['s1x']**2+final_df['s1y']**2+final_df['s1z']**2)*np.sqrt(final_df['s2x']**2+final_df['s2y']**2+final_df['s2z']**2))
final_df['LinS1']=np.arccos(final_df['LinS1'])
final_df['LinS2']=np.arccos(final_df['LinS2'])
final_df['S1S2']=np.arccos(final_df['S1S2'])

final_df.to_csv(save_path+bar+'All_final_variables.dat',sep=' ', index=False)
    
