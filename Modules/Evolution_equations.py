import numpy as np
from Modules.constants import G,c,epsilon
from Modules.sigma_class import Sigma
from scipy.integrate import cumtrapz

class EvolutionEquations(Sigma):
    
    def __init__(self,SimFile,StrainDataPath):
        super().__init__(SimFile,StrainDataPath)
        sigma2,dsigma2,ddsigma2=super().get_sigma2_matrix()
        self.sr=sigma2.real
        self.si=sigma2.imag
        self.dsr=dsigma2.real
        self.dsi=dsigma2.imag
        self.ddsr=ddsigma2.real
        self.ddsi=ddsigma2.imag
        
        self.SIM_NAME=SimFile.replace('ExtrapStrain_','').replace('.h5','')
    
    @staticmethod
    def dot_D(P):
        term1=np.sqrt(2)*P/c
        return term1
    
    @staticmethod
    def dot_J(sr,si,dsr,dsi):
        term1=(c**2/(5*G))*np.einsum('ijk,jl,lk->i',epsilon,si,dsi)
        term2=(c**2/(5*G))*np.einsum('ijk,jl,lk->i',epsilon,sr,dsr)
        return term1+term2
        
    @staticmethod
    def dot_P(dsr,dsi):
        term1=-np.sqrt(2)*c**2/(15*G)*np.einsum('ijk,jl,lk->i',epsilon,dsi,dsr)
        return term1
    
    @staticmethod
    def dot_M(dsr,dsi):
        term1=-(c**2/(10*np.sqrt(2)*G))*np.einsum('ij,ij',dsi,dsi)
        term2=-(c**2/(10*np.sqrt(2)*G))*np.einsum('ij,ij',dsr,dsr)
        return term1+term2
    
    @staticmethod
    def R(D,M,P,sr):#chequear formula
        term1=(1/M)*D
        term2=(1/M)*4*np.sqrt(2)*np.einsum('i,ij->j',P,sr)/(5*c)
        return term1+term2
    
    @staticmethod
    def S(J,R,P):
        term1=J
        term2=np.einsum('ijk,j,k->i',epsilon,P,R)/c
        return term1+term2
       
    def get_dot_D(self,P):
        return np.array(list(map(self.dot_D,P)))
    
    def get_dot_J(self):
        return np.array(list(map(self.dot_J,self.sr,self.si,self.dsr,self.dsi)))
    
    def get_dot_P(self):
        return np.array(list(map(self.dot_P,self.dsr,self.dsi)))
    
    def get_dot_M(self):
        return np.array(list(map(self.dot_M,self.dsr,self.dsi)))
    
    def get_R(self,D,M,P):
        return np.array(list(map(self.R,D,M,P,self.sr)))
    
    def get_S(self,J,R,P):
        return np.array(list(map(self.S,J,R,P)))
    
    def get_D(self,P,initial_D=0):
        return cumtrapz(self.get_dot_D(P).T,self.t,initial=initial_D)
    
    def get_J(self,initial_J=0):
        return cumtrapz(self.get_dot_J().T,self.t,initial=initial_J)
    
    def get_P(self,initial_P=0):
        return cumtrapz(self.get_dot_P().T,self.t,initial=initial_P)
    
    def get_M(self,initial_M=0):
        return cumtrapz(self.get_dot_M(),self.t,initial=initial_M)
    
    def save_all_physics(self,J0,D0,P0,M0,final_quantities=True,print_all_tables=True):
        
        if any(np.array(list(map(len,[J0,D0,P0])))!=3): raise ValueError('Initial conditions must have length 3')
        t=self.t
        dotP=self.get_dot_P()
        dotJ=self.get_dot_J()   
        dotM=self.get_dot_M()   
        deltaP=cumtrapz(dotP.T,self.t,initial=0)
        J0,D0,P0=list(map(np.array,[J0,D0,P0]))
        P=np.add(deltaP.T,P0)
        dotD=self.get_dot_D(P)
        deltaD=cumtrapz(dotD.T,self.t,initial=0)
        deltaJ=cumtrapz(dotJ.T,self.t,initial=0)
        deltaM=cumtrapz(dotM,self.t,initial=0)
        
        M=np.add(deltaM,M0)
        J=np.add(deltaJ.T,J0)
        D=np.add(deltaD.T,D0)
                
        R=self.get_R(D,M,P)
        S=self.get_S(J,R,P)
        
        phy_quant=[t,dotD,dotP,dotM,dotJ,D,P,M,J,R,S,deltaD.T,deltaP.T,deltaM.T,deltaJ.T] #Traspose are made to before trasposed quantities (deltaD,deltaP,deltaJ)
        file_names=['t','dotD','dotP','dotM','dotJ','D','P','M','J','R','S','deltaD','deltaP','deltaM','deltaJ']
        
        if print_all_tables:
            for i in range(len(file_names)):
                with open(self.SIM_NAME+'_'+file_names[i]+'.dat','w') as file_to_save:
                        file_to_save.writelines('###Rochester simulation nº'+self.SIM_NAME+'\n')    
                        file_to_save.writelines('#Evolution of '+file_names[i]+'\n')
                        try:
                            if len(phy_quant[i][0])==3: head='1:x component  2:y component  3:z component'
                        except TypeError:
                            head=''
                        np.savetxt(file_to_save, phy_quant[i], header=head)
        axis=['x','y','z']    
        if final_quantities:
            with open(self.SIM_NAME+'_final_values.dat','w') as file:
                file.writelines('##Final values of the simulation nº'+self.SIM_NAME+'\n')
                i=0    
                for quant in phy_quant:
                    try:
                        j=0
                        for comp in quant[-1]:
                            file.writelines([file_names[i],axis[j],' ',str(comp),'\n'])
                            j+=1
                        i+=1
                    except TypeError:
                        file.writelines([file_names[i],' ',str(quant[-1]),'\n'])
                        i+=1
        
    
    