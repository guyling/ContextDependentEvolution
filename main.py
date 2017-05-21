'''
Created on Jan 16, 2017

@author: Guyling1
'''
import evoModel
import simulations
from simulations import NUCS
import sys
print "ok"
import time
import pandas as pd
import numpy as np
import numpy.random as nr
from featureExtraction import featureExtraction
global NUCS
NUCS=['A','C','G','T']
import glob
import os

mod=evoModel.evoModel('JC',rateParam={'mu':5*10**-5})
from SimulationGenerator import createSimulationRandomly,analyseSample   
location='/sternadi/home/volume1/guyling/MCMC/dataSimulations'
os.chdir(location)
i=int(sys.argv[1])
files=glob.glob("*.csv")
index=int(i//3)
newNuc=int(i%3)
if newNuc==1:
    newNuc=3
f=str(files[index])
analyseSample(f,newNuc=newNuc)
"""    
location='/sternadi/home/volume1/guyling/MCMC/dataSimulations'
file="2motif_{}.csv".format(index)
analyseSample(file,newNuc=newNuc)
"""
#createSimulationRandomly(index,numOfMotifs=int(i/20.)+1)



"""
if newNuc==1:
    newNuc=3
location='/sternadi/home/volume1/guyling/MCMC/dataSimulations'
file="{}motif_{}.csv".format(motifNumber,simulationNumber)
analyseSample(file,newNuc=newNuc)

createSimulationRandomly("1motif_{}".format(i),numOfMotifs=1)
createSimulationRandomly("2motif_{}".format(i),numOfMotifs=2)
createSimulationRandomly("3motif_{}".format(i),numOfMotifs=3) 

     
sim=simulations.EvolustionSimulations(mod,time=15,sampleNum=10000,chroLen=800)
sim.parseChromosome(r'C:\Users\Guyling1\ContextProject\AggarwalaPaper\sabin2.full_genome.U882C.A2973G.C4905U.C5526U.fasta')
sim.addMotif(([-1,-1,1,2,-1],3),2)
sim.addMotif(([-1,1,1,1,-1],3),-2)   
sim.addMotif(([3,-1,1,-1,-1],3),1)  
#sim.parseChromosome(r'C:\Users\Guyling1\Documents\guyReserch\sabin.txt')
sim.initializeProb()
#sim.setName('testRun')
tic=time.time()
#sim.evolve()
moranMatrix=sim.moranModelByPosition()
sim.toFreqsFile(moranMatrix, r'C:\Users\Guyling1\Documents\guyReserch\moranFreqs.freqs')
fc=featureExtraction(r'C:\Users\Guyling1\Documents\guyReserch\moranFreqs.freqs',[(748,7371)],5)
fc.createRegressionTable()  
fc.regTable.to_csv(r'C:\Users\Guyling1\Documents\guyReserch\moranFreqsTry.csv')
#sim.ouputTofasta(r'C:\Users\Guyling1\Documents\guyReserch\outputSimulator.txt')
tac=time.time()
print tac-tic
"""