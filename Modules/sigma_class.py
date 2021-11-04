import numpy as np
from numpy import sqrt,pi,cos,sin
from scipy.interpolate import UnivariateSpline
import os
import h5py
bar='\\' if os.name=='nt' else '/'
class Sigma():
    def __init__(self,SimFile,StrainDataPath):
        self.SimPath=StrainDataPath+bar+SimFile
        self.h5 = h5py.File(self.SimPath,'r')
        self.t=np.array(self.h5['NRTimes'])
        self.dt=abs(self.t[1]-self.t[0])
        self.GetStrainComponents()
        
    def GetStrainComponents(self):
        files=list(self.h5.keys())
        ampls=list(map(lambda x: x.startswith('amp_l2'),files))
        ampls=[files[i] for i in range(len(ampls)) if ampls[i]]        
        phs=list(map(lambda x: x.startswith('phase_l2'),files))
        phs=[files[i] for i in range(len(phs)) if phs[i]]
        self.h2m1,self.h2m2,self.h20,self.h21,self.h22=[self.StrainInterpolation(ampls[i],phs[i]) for i in range(len(ampls))]    
    
    def StrainInterpolation(self,amp,ph,k=5,s=0):
        t=self.t
        xa,ya,xp,yp=list(map(np.array,[self.h5[amp]['X'],self.h5[amp]['Y'],self.h5[ph]['X'],self.h5[ph]['Y']]))

        A = UnivariateSpline(xa, ya,k=k,s=s)
        Phi = UnivariateSpline(xp, yp,k=k,s=s)
        hp=A(t)*cos(Phi(t))
        hc=A(t)*sin(Phi(t))
        return hp+1j*hc
    
    @staticmethod
    def GetSigmaFromh(dt,hp,hc):
        sr=-hp
        si=-hc
        dsr=np.gradient(sr,dt)
        dsi=np.gradient(si,dt)
        ddsr=np.gradient(dsr,dt)
        ddsi=np.gradient(dsi,dt)
        
        return sr+1j*si,dsr+1j*dsi,ddsr+1j*ddsi
    
    def s2xx(self):
        
        hxx= (-1/4)*sqrt(5/pi)*(self.h2m2+self.h22)+(1/6)*sqrt(15/(2*pi))*self.h20
        return self.GetSigmaFromh(self.dt, hxx.real,hxx.imag)
    
    def s2xz(self):
    
        hxz=(1/4)*sqrt(5/pi)*(self.h21-self.h2m1)
        return self.GetSigmaFromh(self.dt, hxz.real,hxz.imag)
    
    def s2xy(self):
    
        hxy=(-1j/4)*sqrt(5/pi)*(self.h22-self.h2m2)
        return self.GetSigmaFromh(self.dt, hxy.real,hxy.imag)
    
    def s2yy(self):
        
        hyy=(1/4)*sqrt(5/pi)*(self.h2m2+self.h22)+(1/6)*sqrt(15/(2*pi))*self.h20
        return self.GetSigmaFromh(self.dt, hyy.real,hyy.imag)
    
    def s2yz(self):
        
        hyz=(1j/4)*sqrt(5/pi)*(self.h2m1+self.h21)
        return self.GetSigmaFromh(self.dt, hyz.real,hyz.imag)
    
 
    def get_sigma2_matrix(self):
        sxx,dsxx,ddsxx=self.s2xx()
        sxy,dsxy,ddsxy=self.s2xy()
        sxz,dsxz,ddsxz=self.s2xz()
        syy,dsyy,ddsyy=self.s2yy()
        syz,dsyz,ddsyz=self.s2yz()

        sigma2=[[sxx,sxy,sxz], \
                [sxy,syy,syz], \
                [sxz,syz,-sxx-syy]]
        
        #dot_sigma_ij matrix
        dsigma2=[[dsxx,dsxy,dsxz], \
                [dsxy,dsyy,dsyz], \
                [dsxz,dsyz,-dsxx-dsyy]]
        
        #ddot_sigma_ij matrix
        ddsigma2=[[ddsxx,ddsxy,ddsxz], \
                [ddsxy,ddsyy,ddsyz], \
                [ddsxz,ddsyz,-ddsxx-ddsyy]]
        
        sigma2,dsigma2,ddsigma2=map(np.array,[sigma2,dsigma2,ddsigma2])
                    
        return np.moveaxis(sigma2,-1,0),np.moveaxis(dsigma2,-1,0),np.moveaxis(ddsigma2,-1,0)
    
    def save_all_sigmas(self):
        if not os.isdir(self.SimPath.replace('.h5','')): os.mkdir(self.SimPath.replace('.h5',''))
        for name in dir(self):
            if name.startswith('s2'):
                method = getattr(self, name)
                sigma,dot_sigma,ddot_sigma=method()
            
                sigmas_list=[sigma,dot_sigma,ddot_sigma]
                der_order=['_','_dot_','_ddot_']
                for i in range(len(sigmas_list)):
                
                    with open(self.SIM_NAME+der_order[i]+name+'.dat','w') as file_to_save:
                        DataOut = np.column_stack((self.t, sigmas_list[i].real, sigmas_list[i].imag))
                        file_to_save.writelines('#Time evolution of '+der_order[i]+name+'\n')
                        np.savetxt(file_to_save, DataOut, header ='1:time  2:real  3:imag')

    