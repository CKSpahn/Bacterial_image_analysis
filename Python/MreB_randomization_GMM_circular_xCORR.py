import numpy as np
import os
import pandas as pd

import matplotlib.pyplot as plt

from sklearn.mixture import GaussianMixture
from scipy.optimize import curve_fit

from scipy import stats 
from scipy.ndimage import gaussian_laplace
from scipy.signal import argrelextrema


#%%

##################################### Methods #####################################

print('Warning: I implemented the normalized cross correlation to deal with amplitude differences')
def ncorr(a,b):
    #normalized correlation computation
    x=np.correlate(a,b,'valid')[0]#I select only the central pic
    div=np.sqrt(np.sum(np.multiply(a,a))) + np.sqrt(np.sum(np.multiply(b,b))) 
    x=np.divide(x,div)
    return x

#circular gaussian model
def gaussian(x,a,mu,sig):
    y = a* np.exp(- (x-mu)**2 / (2* sig**2)) #/ (sig*np.sqrt(2*np.pi))
    y = np.add(y,a* np.exp(- (x-np.max(x)-mu)**2 / (2* sig**2)))#this line allows to deal with circularity right sided (when gaussian overlab)
    y = np.add(y,a* np.exp(- (x+np.max(x)-mu)**2 / (2* sig**2)))#this line allows to deal with circularity left sided (when gaussian overlab)
    return y

def decomposeParameters(param,n_comp):
    a, mu, sig,b = list(param[:n_comp]), list(param[n_comp:2*n_comp]), list(param[2*n_comp:3*n_comp]),param[3*n_comp]
    return a,mu,sig,b

def mergeParameters(a,mu,sig,b):
    param=np.concatenate((a,mu,sig,[b]))
    return param

# GMM with variable number of gaussian
def multi_gaussian(x,n_comp,*args):
    #a, mu, sig,b = list(args[0][:n_comp]), list(args[0][n_comp:2*n_comp]), list(args[0][2*n_comp:3*n_comp]),args[0][3*n_comp]
    a, mu, sig,b = decomposeParameters(args[0],n_comp)
    g=np.ones(len(x))*b #constant background
    for n in range(n_comp):
        g=np.add(g,gaussian(x,a[n],mu[n],sig[n]))
    return g

# same as GMM but returns a list of gaussians instead of the sum ; background b is ignored 
def split_multi_gaussian(x,n_comp,*args):
    a, mu, sig,b = list(args[0][:n_comp]), list(args[0][n_comp:2*n_comp]), list(args[0][2*n_comp:3*n_comp]),args[0][3*n_comp]
    g=[]
    for n in range(n_comp):
        g.append(gaussian(x,a[n],mu[n],sig[n]))
    return g




#%%

##################################### parameters #####################################


sigma=0.15 #initialized size of 1D spots in Mreb signal

number_random_sample=1 #number of simulations per intensity trace

all_plots=False #Activate if all plots for random MreB distribution should be saved; If not activated, the optimized GMM for the first randomization will be saved


folder=r'D:\dSTORM_DATA\GMM\Data'


#%%

p_value_list=dict()




outputFolder=folder+'_randomized_N='+str(number_random_sample)
os.makedirs(outputFolder,exist_ok=True)

plotFolder=outputFolder+'_plots'
os.makedirs(plotFolder,exist_ok=True)

lf=os.listdir(folder)

print(lf)



# Create a file to store the number of components for each file
components_file = os.path.join(plotFolder, 'components.txt')
with open(components_file, 'w') as cf:
    cf.write('Filename\tNumber of Components\tMreB-DNA Correlation\tGMM-DNA Correlation\n')



for i,f in enumerate(lf):
    
    #sample_label='m2'
    #if not f.startswith(sample_label):
    #    continue
    
    print('process ',f)
    path=os.path.join(folder,f)
    
    df=pd.read_csv(path, delimiter="\t")
    col=df.columns.tolist()
    data=dict()
    for c in col:
        data[c]=df[c].tolist()
       
    if i==0:
        plot=True
    else:
        plot=False
        
    if i>=0:
        plot_gmm=True
    else:
        plot_gmm=False
        
        
    #if i>2:
    #    print('WARNING: BREAK')
    #    break
    
    
    x=col[0]
    mreb=col[1]
    dna=col[2]
    

    
    corr_true=ncorr(data[mreb],data[dna])
    

    
    #search local maxima
    gl=-gaussian_laplace(data[mreb],sigma=sigma/data[x][1])#sigma is converted in pixel unit here
    gl[gl<.0]=0
    if all_plots:
        plt.figure()
        plt.plot(data[x],data[mreb],label='Mreb')
        plt.plot(data[x],gl*np.max(data[mreb])/np.max(gl),label=('Detection Laplacian of Gaussian sig='+str(sigma)+' µm'))
        plt.legend()
        plt.xlabel('distance [µm]')  # Add x-axis label
        plt.ylabel('intensity [A.U.]')  # Add y-axis label
        plot_filename = f.replace('.txt', '_LoG_peak_identification_plot.png')
        plt.savefig(os.path.join(plotFolder, plot_filename), bbox_inches='tight', dpi=300)
        plt.close()
    
    local_max=argrelextrema(np.array(gl), np.greater)[0]
    
    n_comp=len(local_max)
    
    
    
    print('i=',i,' nb param = ',n_comp)
    

    mu0=np.array(data[x])[local_max]
    a0=np.array(data[mreb])[local_max] #amplitude init
    sig0=[sigma]*n_comp
    b0=0  #background init
    param=mergeParameters(a0,mu0,sig0,b0)
        
    
    
    if all_plots:
        plt.figure()
        gmm_init=multi_gaussian(data[x],n_comp,param)
        
        plt.plot(data[x],data[mreb],label='MreB')
        plt.plot(data[x],data[dna],label='DNA')
        plt.legend()
        plt.title('correlation(MreB, DNA)='+str(np.round(corr_true,2)))
        plt.xlabel('distance [µm]')  # Add x-axis label
        plt.ylabel('intensity [A.U.]')  # Add y-axis label
        plot_filename = f.replace('.txt', '_MreB_DNA_corr.png')
        plt.savefig(os.path.join(plotFolder, plot_filename), bbox_inches='tight', dpi=300)
        plt.close()
    
    
    
    
    
    if all_plots:
        plt.figure()
        gmm_init=multi_gaussian(data[x],n_comp,param)
        plt.plot(data[x],data[mreb],label='MreB')
        plt.plot(data[x],gmm_init,label='init gmm model')
        plt.legend()
        plt.xlabel('distance [µm]')  # Add x-axis label
        plt.ylabel('intensity [A.U.]')  # Add y-axis label
        plot_filename = f.replace('.txt', '_init_GMM.png')
        plt.savefig(os.path.join(plotFolder, plot_filename), bbox_inches='tight', dpi=300)
        plt.close()
    
    
    
    
    
    
    
    try:
        popt, pcov = curve_fit(lambda x, *params_0: multi_gaussian(x, n_comp, params_0), data[x],data[mreb], p0=param, maxfev=1000,bounds=(0, np.inf))
    except:
        print('WARNING, fit failed, file ignored')
        continue
    
    gmm_opt=multi_gaussian(data[x],n_comp,popt)
    

    
    corr_model=ncorr(gmm_opt,data[dna])
    
    # Save the number of components to the text file
    with open(components_file, 'a') as cf:
        cf.write(f'{f}\t{n_comp}\t{np.round(corr_true, 2)}\t{np.round(corr_model, 2)}\n')

    
    
    if plot_gmm:
        plt.figure()
        
        plt.plot(data[x],data[mreb],label='MreB ; corr='+str(np.round(corr_true,2)))
        plt.plot(data[x],gmm_opt,label='fitted gmm model ; corr='+str(np.round(corr_model,2)))
        g_split=split_multi_gaussian(data[x],n_comp,popt)
        for ii,g in enumerate(g_split): 
            if ii==0: 
                plt.plot(data[x],g,'--',color='gold',alpha=.5,label='gaussians of gmm model')
            else:
                plt.plot(data[x],g,'--',color='gold',alpha=.5)
        plt.legend()
        plt.xlabel('distance [µm]')  # Add x-axis label
        plt.ylabel('intensity [A.U.]')  # Add y-axis label
        plot_filename = f.replace('.txt', '_gmm_plot.png')
        plt.savefig(os.path.join(plotFolder, plot_filename), bbox_inches='tight', dpi=300)
        plt.close()
    
    
    
    
    
    percentageHigher=0
    corr_model_random_list=[] #list of random correlations computed
    a, mu, sig,b = decomposeParameters(popt,n_comp)
    
    
    # Generate and save randomized traces
    for k in range(number_random_sample):
        mu_random = np.random.rand(n_comp) * (np.max(data[x]) - np.min(data[x])) + np.min(data[x])
        p_random = mergeParameters(a, mu_random, sig, b)
        gmm_random = multi_gaussian(data[x], n_comp, p_random)
        
        # Create a new DataFrame for the randomized trace
        randomized_data = pd.DataFrame({
            x: data[x],
            mreb: gmm_random,
            dna: data[dna]
        })
        
        # Save the randomized trace to a separate file
        randomized_filename = f.replace('.txt', f'_random_{k}.txt')
        pathoutput = os.path.join(outputFolder, randomized_filename)
        randomized_data.to_csv(pathoutput, index=False, sep='\t')
        
        if all_plots:
            plt.figure()
            
            plt.plot(randomized_data[x],randomized_data[mreb],label='Random MreB signal')
            plt.plot(randomized_data[x],randomized_data[dna],label='DNA')
            plt.legend()
            plt.xlabel('distance [µm]')  # Add x-axis label
            plt.ylabel('intensity [A.U.]')  # Add y-axis label
            plot_filename = f.replace('.txt', f'_random_{k}.png')
            plt.savefig(os.path.join(plotFolder, plot_filename), bbox_inches='tight', dpi=300)
            plt.close()