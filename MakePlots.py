# File should be placed on first level project folders
import pandas as pd
from Modules.PlotDataAnalysis import PlotDataAnalysis
import os
wd=os.getcwd()
bar='\\' if os.name=='nt' else '/'
result_path=wd+bar+'Results'
save_path=wd+bar+'DataAnalysis'
catalog_path=wd+bar+'Metadata'

plotter=PlotDataAnalysis(save_path,catalog_path,result_path)
df=pd.read_csv(save_path+bar+'All_final_variables.dat',sep=' ',engine='python')

df.reset_index(drop=True,inplace=True)

pars=[['Jin','J','J_{in}','J^f'],
      ['J','S','J^f','S^f'],
      ['Jf_Roch','J','J_{Roch}^f','J^f'],
      ['chi_Roch','chi','\chi_{Roch}','\chi'],
      ['chi_Roch','chiz','\chi_{Roch}','\chi_z']]
#pars_kwargs=[['Jin','J','J_{in}','J^f',{'types':[5]}]]
#pars_fit=[['deltaJ_Roch','deltaJ','\Delta J_{Roch}','\Delta J',{'xtext':0,'ytext':1.4}]]
for parameter in pars:
    plotter.scatter_plot(df,parameter[0],parameter[1],x_label=parameter[2],y_label=parameter[3])
#for parameter in pars_kwargs:
#    plotter.scatter_plot(df,parameter[0],parameter[1],x_label=parameter[2],y_label=parameter[3],**parameter[4])
#for parameter in pars_fit:
#    plotter.scatter_plot_with_fit(df,parameter[0],parameter[1],x_label=parameter[2],y_label=parameter[3],**parameter[4])