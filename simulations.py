'''
Created on Jan 16, 2017

@author: Guyling1
'''
import evoModel
import numpy as np
from numpy import intersect1d
from Bio import SeqIO
from Bio.Seq import Seq
import MoranModel as mm
import scipy.special as sc
from sklearn import preprocessing

global NUCS
NUCS=['A','C','G','T']
class EvolustionSimulations(object):
    '''
    classdocs
    '''


    def __init__(self, evoModel,time,sampleNum,chroLen,k=5,motifs=[],simName='None'):
        '''
        Constructor
        '''
        self.model=evoModel
        self.beta=np.array([])
        self.gamma=np.array([])
        self.time=time
        self.sampleNum=sampleNum
        self.chrolen=chroLen
        self.nullProbMatrix=self.model.getProb(self.time)
        self.n=self.model.getSpaceLen()
        self.outputMatrix=None
        self.motifMatrix=None
        self.chromosomeMatrix=None
        self.motifs=motifs
        self.nullChromosome=None
        self.simulationName=simName
        self.k=k
        
    def setName(self,name):
        self.simulationName=name
        pass
    
    def initializeProb(self):
        chromosome=self.nullChromosome
        self.motifMatrix=self.createMotifMatrix(self.chrolen,chromosome, k=self.k)
        self.chromosomeMatrix=self.createChromosomeTable(self.chrolen,chromosome,self.motifMatrix)
        pass
    
    """
    creating a random chromosome
    """
    def createChromosome(self):
        chrSpace=self.model.getSpaceLen()
        p=self.model.getPiVector()
        chro=np.random.choice(chrSpace,size=self.chrolen,p=p)
        #self.createChromosomeTable(len(chro),chro)
        self.nullChromosome=chro
        print chro
        return chro
    """
    Converting a string of nuc-s in to a np array 
    """
    def nucStringTovector(self,string):
        string=string.upper()
        string=list(string)
        for i in range(len(string)):
            if string[i]=='A':
                string[i]=0
            if string[i]=='C':
                string[i]=1
            if string[i]=='G':
                string[i]=2
            if string[i]=='T':
                string[i]=3
        chro=np.array(string)
        return chro
    """
    #adding the motif to the motif list  and updating the realBeta and gamma values
    #Motif is a tuple (motif,char) where motif is an array of size k over the char space and char is in the char space 
    #example ([-1,-1,2,3,-1],4). -1 is not relevent in the motif
    """
    def addMotif(self,motif,coeff):
        self.motifs.append(motif)
        self.beta=np.append(self.beta,coeff)
        if coeff==0:
            self.gamma=np.append(self.gamma,0)
        else:
            self.gamma=np.append(self.gamma,1)

    """
    taking the list of motifs and creating a binary matrix of (chromlen,motiflen) that indicates which 
    motifs are in effect for each of the positions
    """
    def createMotifMatrix(self,chrolen,chromosome,k=5):
        motifs=self.motifs
        #print chrolen,len(motifs)
        motifMatrix=np.zeros((chrolen,len(motifs)))
        for i in range(len(motifs)):
            motif=np.array(motifs[i][0])
            motifIndices=np.where(motif>=0)[0]
            middleIndex=k//2
            motifLocations=np.where(chromosome==motif[middleIndex])[0]
            for index in motifIndices:
                if index==middleIndex:
                    pass
                else:
                    charLocations=np.where(chromosome==motif[index])[0]
                    shift=index-middleIndex
                    charLocations-=shift
                    motifLocations=np.intersect1d(motifLocations,charLocations)
            motifMatrix[motifLocations,i]=1
        #print np.count_nonzero(motifMatrix)
        #print motifMatrix
        return motifMatrix
        #print np.count_nonzero(self.motifMatrix)
        pass
    

    """
    parsing a fasta file in to a np.vector to work in the evolution evoModel
    """
    def parseChromosome(self,path):
        fasta_sequences = SeqIO.parse(open(path),'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            chro=self.nucStringTovector(sequence)
        self.chrolen=len(chro)
        self.nullChromosome=chro
        pass
            
            
            
    """
    The crhomosome table in a (len,space) table that has the values of change for each of the possible items in 
    the space. Everey row in the matrix corrosponds to a index in the chromosome. This is part of the object so that could be
    created only once per every simulation. after creating the 'null'evoModel in regards to context it's possible to add 
    the context effect
    """
    def createChromosomeTable(self,chrolen,chromosome,motifMatrix):
        probMatrix=self.nullProbMatrix
        #print probMatrix
        n=self.n
        chromMatrix=np.zeros((chrolen,n))
        for i in range(n):#iterating over the space and assinign the evoModel probability
            indices=np.where(chromosome==i)[0] 
            chromMatrix[indices]=probMatrix[i] 
        chromMatrix
        if len(self.motifs)!=0:
            return self.addContextToChromosomeTable(chromosome,chromMatrix,motifMatrix)
        return chromMatrix 
        pass
    """
    finding the index not relevant for the current char mutation
    """
    def motifIndexByChar(self,char):
        indicesOfNotCharMotifs=[]
        motifs=self.motifs
        for i in range(len(motifs)):
            if motifs[i][1]!=char:
                indicesOfNotCharMotifs.append(i)
        return np.array(indicesOfNotCharMotifs)
            
    """
    works in the following way, for each char, cleaving the relevent part of the 
    matrix and realBeta-gamma vectors and then adds to the chromMatrix that has the regular evoModel in 
    it the effects of the motif.te Afr this the table is expit so to make it again a stoch. matrix 
    """
    def addContextToChromosomeTable(self,chromosome,chromMatrix,motifMatrix):
        chromMatrix=sc.logit(chromMatrix)
        #print chromMatrix
        gamma,realBeta=self.gamma,self.beta
        n=self.n
        for char in range(n):
            relevantIndex=self.motifIndexByChar(char)
            charMotifMatrix=np.delete(motifMatrix,np.s_[relevantIndex],1)
            betaBychar=np.delete(realBeta,np.s_[relevantIndex])
            gammaByChar=np.delete(gamma,np.s_[relevantIndex])
            coeffByChar=gammaByChar*betaBychar
            additionOfcontextBychar=np.dot(coeffByChar,charMotifMatrix.T)
            chromMatrix[:,char]+=additionOfcontextBychar
        chromMatrix=sc.expit(chromMatrix)
        chromMatrix=preprocessing.normalize(chromMatrix, norm='l1',axis=1)#Normlizing the matrix back to being a stoch. matrix
        #print np.sum(chromMatrix,axis=1)
        return chromMatrix
        
        
    """
    a no loop vectorized way to mutate the given chromosome. This is much more efficent than the loop version
    but the code is slightly mode complicated 
    """
    def mutateChromosome(self,chromosome):
        chromMatrix=self.chromosomeMatrix
        chromMatrix=np.cumsum(chromMatrix, axis=1)#creating 'buckets' of the probabilities
        randVec=np.random.random((self.chrolen,1))#creating a random vector for all the positions at once
        chromMatrix=np.concatenate((chromMatrix,randVec),axis=1)#concatenating the vector to the table (A,v)
        chromMatrix.sort(axis=1)#sorting the table so we know in which bucket each position is
        ix = np.in1d(chromMatrix.ravel(), randVec).reshape(chromMatrix.shape)#finding the position of the random vectors
        newChrom=np.where(ix)[1]
        return newChrom

    """
    Different kind of simulation, implementing a moran evoModel SSA analysis that is a wright-fisher continous time evoModel
    initFreqs is a (nuc*chromlen) matrix that is the initial probability of each nuc. at each position
    """
    
    def moranModelByPosition(self,initFreqs=None):
        chro=self.nullChromosome
        if initFreqs==None:
            initFreqs=np.zeros((self.chrolen,self.n))
            for i in range(self.n):
                indices=np.where(chro==i)[0]
                #print indices
                if len(indices)>0:
                    initFreqs[indices,i]=10000
            #print initFreqs
        S=mm.createStochMatrixForNucs()
        moranResults=[]
        selection=np.zeros(self.n)
        Nr=S.shape[1]
        t=self.time
        func=mm.changeRate
        k=self.k
        print self.chrolen
        for i in range(self.chrolen)[k//2:-k//2]:
            changeProbs=np.zeros((4,4))#nuc wise  mutation prob matrix for a given position
            for nuc in range(self.n):
                mutChro=np.copy(chro)#a copy of the chromosome, we change the ith nuc to recalculate the motifs for it
                mutChro[i]=nuc
                partOf=mutChro[i-k//2:i+k//2+1]
                partOfMotifMatrix=self.createMotifMatrix(len(partOf),partOf,self.k)
                partOfProbs=self.createChromosomeTable(len(partOf),partOf,partOfMotifMatrix)
                #print partOfProbs
                changeProbs[nuc,:]=partOfProbs[k//2]
            x0=initFreqs[i]
            #print changeProbs
            simResults=mm.ssa(x0, 0, func, t, S, Nr, selection, changeProbs,steps=False)
            moranResults.append(simResults)
            print ("position {},originial_nuc {},results: {}".format(i,chro[i],simResults))
        return np.matrix(moranResults)
        pass
    
    def evolve(self,inMemorySize=10**5):
        outputMatrix=np.zeros((self.sampleNum,self.chrolen))
        originalChro=self.nullChromosome
        print ("original chro:")
        print originalChro
        print("----------------")
        for i in range(self.sampleNum):
            chr=self.mutateChromosome(originalChro)
            outputMatrix[i]=chr
        self.outputMatrix=outputMatrix
    """
    taking the matrix of vectors and creating string to be outputted as fasta files
    This is much faster than the iterative version though slightly more complicated in the code
    """
    def numpyMatrixToNucString(self,matrix):
        outputChr=np.chararray(matrix.shape)#creating a spacesholder of chararray type
        for i in range(self.n):#iterating over the char space
            indices=np.where(matrix==i)# looking for the indices 
            outputChr[indices[0],indices[1]]=NUCS[i]  #putting the right nucs in the indices
        spaces=np.chararray((outputChr.shape[0],1)) #creating a delimiter array 
        spaces[:]="*" #'*' is the delimiter
        outputChr=np.concatenate((outputChr,spaces),axis=1) #concatenating the delimiter in the end of each row
        strings=outputChr.tostring().split('*')[:-1] # returning a list of strings 
        return strings#returning the string 
    
    """
    outputng a moran model simulation to a .Freqs file and format for future inquiry 
    """
    
    def toFreqsFile(self,moranMatrix,ouputName):
        f=open(ouputName,'w')
        f.write('#moran model results\n')
        f.write('Pos    Base    Freq    Ref    Read_count    Rank    Prob\n')
        for i in range(moranMatrix.shape[0]):
            simulationPositionResult=np.asarray(moranMatrix[i])[0]
            #print simulationPositionResult
            for nuc in range(self.n):
                pos=i+1#we only check the ones in this range because they have context 
                ref=NUCS[self.nullChromosome[i+self.k//2]]
                read_count=10000
                base=NUCS[nuc]
                freq=simulationPositionResult[nuc]/read_count
                num=simulationPositionResult[nuc]
                rank=np.where(np.sort(simulationPositionResult)==num)[0][0]
                prob=1.0
                line=[pos,base,freq,ref,read_count,rank,prob]
                line=[str(x) for x in line]
                f.write("    ".join(line))
                f.write("\n")
        f.close()
                
            
    def ouputTofasta(self,path):
        f=open(path,'w')
        matrix=self.outputMatrix
        sampleNum=matrix.shape[0]
        results=self.numpyMatrixToNucString(matrix)
        for i in range(len(results)):
            id=self.simulationName+"_sampleNumber_{}".format(i)
            seq=results[i]
            f.write(id+'>')
            f.write('\n')
            f.write(seq)
            f.write('\n')
        f.close()
            
            
            
        
        
    
        
        
    

            
            
            
            
            
        
        
        