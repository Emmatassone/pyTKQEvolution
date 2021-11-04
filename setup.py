import pandas as pd
import os
from Modules.Evolution_equations import EvolutionEquations 

bar='\\' if os.name=='nt' else '/'
path_of_run_file=os.getcwd()
MetadataPath=path_of_run_file+bar+'Metadata'
StrainDataPath=path_of_run_file+bar+'Strain_data'

if not os.path.isdir('Results'):os.makedirs('Results')
results_folder=path_of_run_file+bar+'Results'

def load_line(df,word):
    return float(df[df['[metadata]']==word].iloc[0,1].replace('=',''))

def read_metadata(step,MetadataPath,ADM_Fail=False):
    df=pd.read_fwf(MetadataPath,sep='/t',engine='python')
    
    rel_time,Jx,Jy,Jz=list(map(load_line,[df,df,df,df],['relaxed-time','initial-ADM-angular-momentum-x','initial-ADM-angular-momentum-y','initial-ADM-angular-momentum-x']))
    if ADM_Fail==False:
        initial_mass=load_line(df,'initial-ADM-energy')
    else:
        initial_mass=load_line(df,'initial-total-mass')
    rel_time*=(1/step)
    return int(rel_time),[Jx,Jy,Jz],initial_mass

def run_calculation(StrainFile,MetadataFile):
    sigma=EvolutionEquations(StrainFile,StrainDataPath)
    
    try:
        t_relax,J0,M0=read_metadata(sigma.dt,MetadataPath+bar+MetadataFile)
    except ValueError:
        print('Simulation '+MetadataFile+' can not load some variable from metadata')
        t_relax,J0,M0=read_metadata(sigma.dt,MetadataPath+bar+MetadataFile,ADM_Fail=True)

    D0,P0=[0,0,0],[0,0,0]
    result_directory=results_folder+bar+sigma.SIM_NAME
    if not os.path.isdir(result_directory): os.makedirs(result_directory)
    os.chdir(result_directory)
    sigma.save_all_physics(J0,D0,P0,M0)
    