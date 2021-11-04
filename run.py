import multiprocessing as mp
import os
import time
from setup import run_calculation 

if __name__ == "__main__":
    start_time=time.time()
    bar='\\' if os.name=='nt' else '/'
    path_of_run_file=os.getcwd()
    StrainFiles=os.listdir(path_of_run_file+bar+'Strain_data')
    MetadataFiles=os.listdir(path_of_run_file+bar+'Metadata')
    procs=int(mp.cpu_count()/2)
    if len(StrainFiles)!=len(MetadataFiles): raise IndexError('Number of files in folder Metadata and Strain_data do not match')
    with mp.Pool(procs) as p:
        p.starmap(run_calculation,zip(StrainFiles,MetadataFiles))
        
    print("-----------------------------------------")
    print("Ran in all folders in ",time.time()-start_time,"seconds")