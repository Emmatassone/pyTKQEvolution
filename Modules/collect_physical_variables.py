import os
import pandas as pd
import numpy as np 

class CollectPhysicalVariables:
    def __init__(self,save_path,catalog_path=None,result_path=None):
        self.catalog_path=catalog_path
        self.save_path=save_path
        self.result_path=result_path
        self.bar='\\' if os.name=='nt' else '/'
        pass
    
    def IsCatalogPathInput(self):
        if not self.catalog_path: raise ValueError('Catalog path must be an absolute path direction')
    
    def IsResultsPathInput(self):
        if not self.result_path: raise ValueError('Results path must be an absolute path direction')
        
    def load_line(self,df,word):
        return float(df[df['[metadata]']==word].iloc[0,1].replace('=',''))

    def classify_spins(self):
        self.IsCatalogPathInput()
        
        metafiles=os.listdir(self.catalog_path)
        nonspinning,aligned,precessing=[],[],[]
        for file in metafiles:
            new_path=self.catalog_path+self.bar+file
            df=pd.read_fwf(new_path,sep='/t',engine='python')
            system_type=df[df['[metadata]']=='system-type'].iloc[0,1].replace('=','').strip()
                
            if system_type=='Nonspinning':
                nonspinning.append(file)
            elif system_type=='Aligned':
                aligned.append(file)
            elif system_type=='Precessing':
                precessing.append(file)
            else:
                print('Simulation ',file,'has not been classified')
                      
        with open(self.save_path+self.bar+'Nonspinning_simulations.dat','w') as f:
            for lines in nonspinning:
                f.write(lines+'\n')
        with open(self.save_path+self.bar+'Aligned_simulations.dat','w') as f:
            for lines in aligned:
                f.write(lines+'\n')
        with open(self.save_path+self.bar+'Precessing_simulations.dat','w') as f:
            for lines in precessing:
                f.write(lines+'\n')        
                
    def classify_masses(self):
        self.IsCatalogPathInput()
        
        metafiles=os.listdir(self.catalog_path)
        equal_masses,nonequal_masses=[],[]
        
        for file in metafiles:
            new_path=self.catalog_path+self.bar+file
            df=pd.read_fwf(new_path,sep='/t',engine='python')
                
            mass_1=self.load_line(df,'initial-mass1')
            mass_2=self.load_line(df,'initial-mass2')
            
            if mass_1==mass_2:
                equal_masses.append(file)
            else:
                nonequal_masses.append(file)
                
        with open(self.save_path+self.bar+'Equal_masses_simulations.dat','w') as f:
            for lines in equal_masses:
                f.write(lines+'\n')
        with open(self.save_path+self.bar+'Non_equal_masses_simulations.dat','w') as f:
            for lines in nonequal_masses:
                f.write(lines+'\n')
                
    def classify_spins_and_masses(self):
        self.classify_spins()
        self.classify_masses()
    
    def ReadMetadata(self,df,ADM_Fail=False):
  
        system_type=df[df['[metadata]']=='system-type'].iloc[0,1].replace('= ','')
        variables_to_extract=['relaxed-time','initial-ADM-angular-momentum-x','initial-ADM-angular-momentum-y','initial-ADM-angular-momentum-z','final-chi','final-mass','final-kick','peak-luminosity-ergs-per-sec','initial-orbital-angular-momentum','eccentricity','initial-separation','initial-mass1','initial-mass2']
        rel_time,Jxr,Jyr,Jzr,finalChi,finalMass,kick,peakL,Lin,e,r0,mass1,mass2=list(map(self.load_line,[df,df,df,df,df,df,df,df,df,df,df,df,df],variables_to_extract))
        
        if ADM_Fail==False:
            initial_mass=self.load_line(df,'initial-ADM-energy')
        else:
            initial_mass=self.load_line(df,'initial-total-mass')
        
        if(system_type=='Precessing'):
            chix1,chiy1,chiz1,chix2,chiy2,chiz2=list(map(self.load_line,[df,df,df,df,df,df],['initial-bh-chi1x','initial-bh-chi1y','initial-bh-chi1z','initial-bh-chi2x','initial-bh-chi2y','initial-bh-chi2z']))
            
        else:
            chix1,chiy1,chiz1=0,0,self.load_line(df,'initial-bh-chi1z')
            chix2,chiy2,chiz2=0,0,self.load_line(df,'initial-bh-chi2z')
            
        s1x,s1y,s1z=chix1*mass1**2,chiy1*mass1**2,chiz1*mass1**2
        s2x,s2y,s2z=chix2*mass2**2,chiy2*mass2**2,chiz2*mass2**2
            
        return rel_time,initial_mass,Jxr,Jyr,Jzr,finalChi,finalMass,kick,peakL,Lin,e,r0,mass1,mass2,s1x,s1y,s1z,s2x,s2y,s2z

    def metadata_collect_physics(self):
        self.IsCatalogPathInput()
        metafiles=os.listdir(self.catalog_path)
        d = {'EM-NS': [], 'EM-A': [],'EM-P': [], 'NEM-NS': [],'NEM-A': [], 'NEM-P': [],'Relaxed_time':[],'Jx_adm':[],'Jy_adm':[],'Jz_adm':[],'Jf_Roch':[],'initial_mass':[],'final_mass':[],'deltaJ_Roch':[],'deltaM_Roch':[],'Kick':[],'Lmax':[],'chi_Roch':[],'Lin':[],'e':[],'r0':[],'m1':[],'m2':[],'s1x':[],'s1y':[],'s1z':[],'s2x':[],'s2y':[],'s2z':[]}
        final_quantities=pd.DataFrame(data=d)
        
        files=['Nonspinning_simulations.dat','Aligned_simulations.dat','Precessing_simulations.dat','Equal_masses_simulations.dat','Non_equal_masses_simulations.dat']
        files=[self.save_path+self.bar+file for file in files]
        for file in files:
            if not os.path.isfile(file):
                self.classify_spins_and_masses()
                break
        NS_list,A_list,P_list,EM_list,NEM_list=[pd.read_csv(file,header=None) for file in files]
        
        for file in metafiles:
            new_path=self.catalog_path+self.bar+file
            df=pd.read_fwf(new_path,sep='/t',engine='python')
            try:
                rel_time,initial_mass,Jxr,Jyr,Jzr,finalChi,finalMass,kick,peakL,Lin,e,r0,mass1,mass2,s1x,s1y,s1z,s2x,s2y,s2z=self.ReadMetadata(df)
            except ValueError:
                try:
                    rel_time,initial_mass,Jxr,Jyr,Jzr,finalChi,finalMass,kick,peakL,Lin,e,r0,mass1,mass2,s1x,s1y,s1z,s2x,s2y,s2z=self.ReadMetadata(df,ADM_Fail=True)
                except ValueError:
                    print('Simulation '+file+' can not load some variable from metadata')
                    continue
            
            Jf_Roch=finalChi*finalMass**2
            variables={'Relaxed_time':rel_time,'Jx_adm':Jxr,'Jy_adm':Jyr,'Jz_adm':Jzr,'Jf_Roch':Jf_Roch,'initial_mass':initial_mass,'final_mass':finalMass,'deltaJ_Roch':Jf_Roch-np.sqrt(Jxr**2+Jyr**2+Jzr**2),'deltaM_Roch':finalMass-initial_mass,'Kick':kick,'Lmax':peakL,'chi_Roch':finalChi,'Lin':Lin,'e':e,'r0':r0,'m1':mass1,'m2':mass2,'s1x':s1x,'s1y':s1y,'s1z':s1z,'s2x':s2x,'s2y':s2y,'s2z':s2z}
            if file in EM_list.values.flatten():
                    if file in NS_list.values.flatten():
                        d=dict({'EM-NS':'s'+file[8:12]},**variables)
                        final_quantities=final_quantities.append(d,ignore_index=True)
                    elif file in A_list.values.flatten():
                        d=dict({'EM-A':'s'+file[8:12]},**variables)
                        final_quantities=final_quantities.append(d,ignore_index=True)
                    elif file in P_list.values.flatten():
                        d=dict({'EM-P':'s'+file[8:12]},**variables)
                        final_quantities=final_quantities.append(d,ignore_index=True)
                    else:
                        print(file,'is not classified')
                        pass
                    
            elif file in NEM_list.values.flatten():
                    if file in NS_list.values.flatten():
                        d=dict({'NEM-NS':'s'+file[8:12]},**variables)
                        final_quantities=final_quantities.append(d,ignore_index=True)
                    elif file in A_list.values.flatten():
                        d=dict({'NEM-A':'s'+file[8:12]},**variables)
                        final_quantities=final_quantities.append(d,ignore_index=True)
                    elif file in P_list.values.flatten():
                        d=dict({'NEM-P':'s'+file[8:12]},**variables)
                        final_quantities=final_quantities.append(d,ignore_index=True)
                    else:
                        print(file,'is not classified')
                        pass
        final_quantities.fillna('nan').to_csv(self.save_path+self.bar+'Metada_simulations_variables.dat',sep=' ', index=False)
        
    def collect_final_files(self):
        self.IsResultsPathInput()
        
        folders=os.listdir(self.result_path)
        finalValues=pd.DataFrame([])
        for folder in folders:
            new_path=self.result_path+self.bar+folder
            df=pd.read_csv(new_path+self.bar+folder+'_final_values.dat',skiprows=1,header=None,delim_whitespace=True)
            M=np.loadtxt(new_path+self.bar+folder+'_dotM.dat',usecols=[0])
            max_M=min(M)
            new_df=pd.DataFrame([df[1].values],columns=df[0].values)
            new_df['Simulation']='s'+folder[8:12]
            new_df['maxM']=max_M
            finalValues=finalValues.append(new_df)
        
        finalValues.fillna('nan').to_csv(self.save_path+self.bar+'Final_results_simulations.dat',sep=' ', index=False)
    