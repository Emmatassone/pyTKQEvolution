import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from scipy.optimize import curve_fit
from Modules.collect_physical_variables import CollectPhysicalVariables

class PlotDataAnalysis(CollectPhysicalVariables):
    def __init__(self,save_path,catalog_path=None,result_path=None):
        super().__init__(save_path,catalog_path,result_path)
        self.save_path=save_path
        self.bar='\\' if os.name=='nt' else '/'
        if not os.path.isdir(save_path): os.makedirs(save_path)
        if not os.path.isfile(save_path+self.bar+'Metada_simulations_variables.dat'): super().metadata_collect_physics()
        if not os.path.isfile(save_path+self.bar+'Final_results_simulations.dat'): super().collect_final_files()
        
    def linear_function(self,x,a):
        return a*x
    
    def scatter_plot(self,df,x_name,y_name,x_label=0,y_label=0,markers=['s','>','x','d','+','.'],**kwargs):
        if not isinstance(df, pd.DataFrame): raise TypeError('first argument must be of type Pandas DataFrame')
        
        if 'filter' in kwargs: df=df[df[kwargs['filter'][0]].astype(float)<kwargs['filter'][1]]
        x_pos=df.columns.get_loc(x_name)
        y_pos=df.columns.get_loc(y_name)
        all_x=np.sort(np.array(df.loc[:,x_name].values))
        fig, ax = plt.subplots(1, 1)
        types=range(6) if not 'types' in kwargs else kwargs['types']
        for i in types:
            to_plot=df.iloc[:,[i,x_pos,y_pos]].dropna()
            x_axis=to_plot.loc[:,x_name].values
            y_axis=to_plot.loc[:,y_name].values
            ax.plot(x_axis,y_axis, markers[i],label=to_plot.columns[0])
        
        if ('reference_line' in kwargs) and (kwargs['reference_line']==True):
            ax.plot(all_x,all_x)
        loc=kwargs['legend_loc'] if 'legend_loc' in kwargs else 'lower right'
        fig.legend(numpoints=1,loc=loc,bbox_to_anchor=(0.08,.98), frameon=False)
        if not(x_label or y_label): x_label,y_label=x_name,y_name  
        plt.xlabel('$'+x_label+'$')
        plt.ylabel('$'+y_label+'$', rotation=0)
        plt.tight_layout()
        plt.savefig(self.save_path+self.bar+x_name+'_vs_'+y_name+'.png') 
        
    def scatter_plot_with_fit(self,df,x_name,y_name,x_label=0,y_label=0,markers=['s','>','x','d','+','.'],**kwargs):
        if not isinstance(df, pd.DataFrame): raise TypeError('first argument must be of type Pandas DataFrame')
        if 'filter' in kwargs: df=df[df[kwargs['filter'][0]].astype(float)<kwargs['filter'][1]]
        
        x_pos=df.columns.get_loc(x_name)
        y_pos=df.columns.get_loc(y_name)
        
        All_x=np.array([float(i) for i in df.iloc[:,x_pos].values])
        All_y=np.array([float(i) for i in df.iloc[:,y_pos].values])
        fig, (ax1,ax2) = plt.subplots(2, 1)
        types=range(6) if not 'types' in kwargs else kwargs['types']
        for i in types:
            to_plot=df.iloc[:,[i,x_pos,y_pos]].dropna()
            x_axis=to_plot.loc[:,x_name].values
            y_axis=to_plot.loc[:,y_name].values
            ax1.plot(x_axis,y_axis, markers[i],label=to_plot.columns[0])  
        
        if not(x_label or y_label): x_label,y_label=x_name,y_name  
        ## Fit plot#####
        pars, cov = curve_fit(f=self.linear_function, xdata=All_x, ydata=All_y)
        stdevs = np.sqrt(np.diag(cov))
        #Plots
        ax1.plot(All_x,self.linear_function(All_x,*pars),label='$'+y_label+'=a'+x_label+'$',linestyle=':',color='red')
        
        if 'xtext' and 'ytext' in kwargs:
            font = {'family': 'serif',
                'color':  'black',
                'weight': 'normal',
                'size': 10
                }
            ax1.text(kwargs['xtext'],kwargs['ytext'],'$a=$'+str(round(pars[0],8))+r'$\pm$'+str(round(stdevs[0],8)),fontdict=font)
            
        ax1.set_xlabel('$'+x_label+'$')
        ax1.set_ylabel('$'+y_label+'$',rotation=0,loc='top')
        ##Residuals plot###############
        res = All_y- self.linear_function(All_x, *pars)
        ax2.plot(All_x,res,'.',color='orange')
        ax2.plot(All_x,np.zeros(len(All_x)),'-',color='grey')
        ax2.set_xlabel('$'+x_label+'$')
        ax2.set_ylabel('$Residuals$')
        
        
        ##Residuals plot###############
        fig.tight_layout()
        loc=kwargs['legend_loc'] if 'legend_loc' in kwargs else 'lower right'
        fig.legend(numpoints=1,loc=loc)
        fig.savefig(self.save_path+self.bar+x_name+'_vs_'+y_name+'_with_fit.png')
