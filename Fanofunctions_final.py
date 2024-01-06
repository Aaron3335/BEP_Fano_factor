#DISCLAIMER
#This code was written as part of 
#Author: Aaron Dantuma

#LIBRARIES
import numpy as np

#DISTRIBUTIONS 

#input paramaters are phonon energy (om), quasiparticle energy (E), gap energy (Delta, 1 if you are doing calculations in units of delta)

def alpha2F(om): #quadratic approximation for the phonon spectrum
    return om**2

def pairbreakingdistributions(E,om,Delta): #distribution for the pair breaking interaction
    omminE = np.round(om-E,5)
    dist =  1/(E**2-Delta**2)**0.5   *(E*(omminE)+Delta**2)/((omminE)**2-Delta**2)**0.5
    return dist
   

def phonon_emiss_dist(om,E,spectrum,Delta): #distribution for the phonon emission interaction
    dist = spectrum(om)*(E-om)/((E-om)**2-Delta**2)**0.5   *(1-Delta**2/(E*(E-om)))
    return dist



#DISCRETIZATION
#These fucntions calculate the possible outcomes of all interactions and corresponding chances. 
#The outcomes are indexes in the energy-array, which used as an input to define what energies are possible in the simulation

def descretize_init_phon(spectrum,energy): #the initial phonon spectrum

    om_allowed = energy
    outcomes = np.arange(len(om_allowed))
    chances = spectrum(om_allowed)/sum(spectrum(om_allowed))
    print("Finished initial phonon descretisation")
    return outcomes, chances

def descretize_pair_breaking(energy,stepsize,Delta): #the pair breaking interaction
    om_allowed  =energy
    outcomes = []
    chances = []
    loc_delta = round(1/stepsize)
    indexarray = np.arange(0,len(energy)+1,1)
    for i in range(len(om_allowed)):
        indexes = indexarray[loc_delta:i-loc_delta+1] #define the allowed range of the outcomes
        outcomes.append(indexes)
        E = energy[indexes]
        chance = pairbreakingdistributions(E,om_allowed[i],Delta) #calculate the value of the probability distribution
        if len(chance) > 1:
            chance[-1] = pairbreakingdistributions(E[-1]-stepsize/4,om_allowed[i],Delta)/2 #the last bin needs to be halved in order to prevent an infinite chance, see report
        norm_chance = chance/sum(chance) #normalise the chances
        chances.append(norm_chance)
    print("finished pair breaking descretisation")
    return outcomes, chances

def descretize_phonon_emiss(spectrum,energy,stepsize,Delta):
    E_allowed = energy

    outcomes = []
    chances = []

    indexarray = np.arange(0,len(energy)+1,1)
    loc_delta = round(1/stepsize)

    for i in range(len(E_allowed)):
        indexes = indexarray[0:(i+1)-loc_delta] #define the allowed range of the outcomes
        outcomes.append(indexes)

        om = energy[indexes]
        chance = phonon_emiss_dist(om,E_allowed[i],spectrum,Delta) #calculate the value of the probability distribution
        if len(chance) > 1:
                chance[-1] = phonon_emiss_dist(om[-1]-stepsize/4,E_allowed[i],spectrum,Delta)/2  #the last bin needs to be halved in order to prevent an infinite chance, see report
        norm_chance = chance/sum(chance) #normalise the chances
        chances.append(norm_chance)

    print("finished phonon emission descretisation")
    return outcomes, chances


