#LIBRARIES
import numpy as np
import time
from Fanofunctions_final import alpha2F, descretize_init_phon, descretize_pair_breaking,descretize_phonon_emiss
from simulationfunction import simulation,savedata

#CONSTANTS
Delta = 1 #because we calculate in units of Delta
h = 4.136e-15 #eV s. planck constant
kb = 8.617343e-5 #eV K^-1, Boltzmann constant

#Material
realDelta,wd= 5.75e-4, 200# eV,K Gap energy and bebye temperature of Sn

# realDelta, wd= 15.25e-4, 276 # eV Gap energy and debye temperature of Nb
# realDelta, wd= 1.7e-4, 433 #eV,K gap energy and debye temperature of Al
#realDelta, wd = 0.75e-4, 209 #eV, K, gap energy and debye temperatue of Cd

#general allowed energies of E and Omega
Debye = wd*kb/realDelta


#Simulation parameters
cycles = 1 #number of iterations
spectrum = alpha2F   #phonon spectrum used
control = True #whether or not control variables are printed

stepsize = [1/200]
Q = [10000]
spectra = [alpha2F]

#we use lists to be able to quickly switch between different imput parameters
for j in range(len(spectra)):
    for p in range(len(Q)):
        for i in range(len(stepsize)):

            #initialising and discretisation
            start = stepsize[i]/2
            energy  = np.arange(start,Debye+start,stepsize[i])

            tstart = time.time()
            init_phon = descretize_init_phon(spectra[j],energy)
            pairbreaking = descretize_pair_breaking(energy,stepsize[i],Delta)
            emission = descretize_phonon_emiss(spectra[j],energy,stepsize[i],Delta,)
            print("Finished discretisation in {} s".format(np.round(time.time()-tstart,2)))
        
            #execute simulation
            N, ph, efficiency, E_mean, om_mean = simulation(Q[p],energy,cycles,stepsize[i],pairbreaking,emission,init_phon,control) 
        
            spec = spectra[j].__name__ #give your data a unique name so it is not overwritten
            
            #saving some data, make sure you have a folder named "data" 
            savedata(N,'N',spec,cycles,wd)
            savedata(ph,'ph',spec,cycles,wd)
            savedata(efficiency,'eff',spec,cycles,wd)
            savedata(E_mean,'E',spec,cycles,wd)
            savedata(om_mean,'om',spec,cycles,wd)