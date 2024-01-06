import numpy as np
import time
import json
#this file contains the code that uses the discretisation of the quasiparticle-phonon system from Fanofunctions_final to run the simulation


#PICKING FUNCTIONS
#these functions take the discretisations from Fanofunctions_final.py as input 
#and use them to randomly pick the outcome of the interactions based on the quasiparticle/phonon energy

def pick_init_phon(Q,init_phon,energy,stepsize): #function to initialise the first group of phonons with energy equal to Q, the photon energy
    outcomes,chances = init_phon
    phonons = np.zeros(len(energy))
    
    while sum(phonons*energy) < Q:
        index = np.random.choice(outcomes,p=chances)
        phonons[index] += 1
    phonons[index] += -1
    totalenergy = np.sum(phonons*energy)
    
    lastenergy = int(Q-totalenergy)
    lastindex = int(lastenergy/stepsize)
    phonons[lastindex] += 1
    return phonons


def pick_quasis(index,pairbreaking,energy):
    outcomes, chances = pairbreaking
    try:
        out = outcomes[index]
        probs = chances[index]

        energy1 = np.random.choice(out,p=probs)
        energy2 = index-energy1
    except:
        print(energy[index])
        print(energy[outcomes[index]])
        print(outcomes[index],chances[index])
    quasis  = np.array([energy1,energy2])

    diff = energy[energy1]+energy[energy2]-energy[index]
    if np.abs(diff) > 1/10:
        print("no conservation:", diff,'Delta deviation')
        print(index,energy1,energy2)

    return quasis.astype(int)


def pick_phonon(index,emission,energy): 
    outcomes, chances = emission[0], emission[1]
    try:
        phonon = np.random.choice(outcomes[index],p=chances[index])
    except:
        print(energy[index])
        print(outcomes[index],chances[index])
    return phonon.astype(int)


#main simulation loop
def simulation(Q,energy,cycles,stepsize,pairbreaking,emission,init_phon,control):
    start_time = time.time()
    N = np.zeros(cycles)
    ph = np.zeros(cycles)
    efficiency = np.zeros(cycles)
    E_mean = np.zeros(cycles)
    om_mean = np.zeros(cycles)

    L = len(energy)
    qlimit = round(3/stepsize) #limits where energies of phonons and quasiparticle cannot cause the creation of more quasiparticles
    plimit = round(2/stepsize)
    pairbreaking[1][plimit] = [1] 
    #at a phonon energy 2*Delta, the discretisation does not function properly due to the fact that there is just one outcome (two Delta quasparticles) 
    #possible and this results in infinite probability, which is why we need to set in manually

    for i in range(cycles):
        #some text to make the interface more clear
        print("\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~") #to separate iterations
        print("Cycle number: {}".format(i+1))
        #initialise the first phonons and an empty list for the quasiparticles
        quasiparticles = np.zeros(L)
        phonons = np.zeros(L)
        initphonons = (pick_init_phon(Q,init_phon,energy,stepsize))
        phonons += initphonons

        #cycle that repeats until all interacting particles are exhausted
        while sum(phonons[plimit:])>0 or sum(quasiparticles[qlimit:])>0:
            
            for p in range(plimit,L):
                if phonons[p] > 0:
                    for j in range(int(phonons[p])):
                        phonons[p] += -1
                        indexes = pick_quasis(p,pairbreaking,energy)
                        quasiparticles[indexes] += 1
            
            for q in range(qlimit,L):
                if quasiparticles[q] > 0:
                    for j in range(int(quasiparticles[q])):
                        #print(quasiparticles[q])
                        quasiparticles[q] += -1
                        index = pick_phonon(q,emission,energy)
                        phonons[index] += 1
                        newquasis = q - index
                        quasiparticles[newquasis] += 1

            

        #the sum of the quasiparticle array gives the no. of quasiparticles
        N[i]= sum(quasiparticles)
        ph[i] = sum(phonons)
        efficiency[i]= (sum(quasiparticles)*(1+stepsize/2)/Q)
        E_mean[i] = sum(quasiparticles*energy)/sum(quasiparticles)
        om_mean[i] = sum(phonons*energy)/sum(phonons)

        energycontrol = sum(phonons*energy)+sum(quasiparticles*energy)
        stop_time = time.time()
        Delta_t = stop_time-start_time
        print("Time passed: {} s".format(round(Delta_t,2)))
        if control == True: #print control variables if True
            print("Efficiency: {} ".format(efficiency[i]))
            print("Number of quasiparticles: {}".format(N[i]))
            print("Total energy: {} \u0394".format(round(energycontrol),3))
            print("Phonon energy: {} \u0394".format(round(sum(phonons*energy))))
            print("Quasiparticle energy: {} \u0394".format(round(sum(quasiparticles*energy))))
            print("Mean phonon energy: {} \u0394".format(round(om_mean[i],2)))
            print("Mean quasiparticle energy: {} \u0394".format(round(E_mean[i],2)))
            print("Energy deficiency: {} \u0394".format(sum(initphonons*energy)-energycontrol))

    return N, ph, efficiency, E_mean,om_mean

#SAVING DATA
def savedata(data,name,spec,cycles,wd):
    filename = name+str(spec)+"_"+"wd"+str(int(wd))+"_"+"2024-1-3"+"_"+str(cycles)+"its"
    data = list(data)
    with open("Data\\"+filename, 'w') as file:
        json.dump(data, file)

