from __future__ import division
import pymc3 as pm
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set() 
from theano import shared
import scipy as sc
import theano.tensor as T
import pandas as pd
from theano import as_op
import featurizeMutationFile as f 
from sets import Set
import time
import gammaPrior as gp
from pymc3.backends.base import merge_traces
from reverseJumpSampler import Metropolis
plt.switch_backend('agg')



"""
size=2000
initNum=np.random.poisson(lam=300,size=size)

pinit=np.array([-2.1972245773362191 for i in range(size)])
X=np.random.randint(0,high=2,size=(size,covars))
#print X.shape
#covarMatrix=np.dot(X.T,X)*size**-1
#print covarMatrix.shape
realBeta=np.random.multivariate_normal(np.zeros(covars),np.identity(covars))
realBeta[0]=5
realBeta[1]=-5
#realBeta[-1]=-5
realGamma=np.zeros(covars)
realGamma[0]=1
realGamma[1]=1
#realGamma[-1]=1
betaGamma=realBeta*realGamma
means=pinit+np.dot(X,betaGamma)
means=sc.special.expit(means)
#print means[:100]
y_=np.random.binomial(initNum,means)
#print y_[0:100]
#Split Data
X_tr, X_te, y_tr, y_te,initNum_tr,initNum_te = train_test_split(X,y_,initNum,test_size=0.25, random_state=0)
"""

def featureListToGraphMatrix(features,sigFeaturesIndices=[]):
    n=len(features)
    if len(sigFeaturesIndices)==0:
        sigFeaturesIndices=[i for i in range(n)]
    graphMatrix=np.zeros((n,n)) #creating an empty graph connectivity matrix
    for i in sigFeaturesIndices:
        feature_i=features[i] #takign only the features that are sig from before
        positions_i,Nucs_i=feature_i.split("__")#spliting to make a list of positions
        positions_i=positions_i.split("_")# list from the string of the motif
        Nucs_i=list(Nucs_i)
        #print positions_i,Nucs_i
        for j in sigFeaturesIndices:
            feature_j=features[j]
            positions_j,Nucs_j=feature_j.split("__")
            positions_j=positions_j.split("_")
            Nucs_j=list(Nucs_j)
            if len(positions_i)==len(positions_j):
                continue
            minSize=min(len(positions_i),len(positions_j))# the size of the smaller motif of the two
            S_i,S_j=Set(positions_i),Set(positions_j) #crating a set from the positions
            if len(S_i&S_j)==minSize:# meaning small set is completly in the big set 
                N_i, N_j=Set(Nucs_i),Set(Nucs_j)#sets of the Nucs 
                if len(N_i&N_j)==minSize and len(positions_i)==minSize:
                    graphMatrix[i,j]=1
    #print np.nonzero(graphMatrix)
    return graphMatrix

@as_op(itypes=[T.dvector], otypes=[T.dvector])
def exp_it(probV):
    #print probV
    return sc.special.expit(probV)

@as_op(itypes=[T.lvector,T.dscalar,T.lscalar,T.lscalar], otypes=[T.dvector])
def GammaCovariance(gamma,sigma,shape,C):
    #print sigma
    indices=np.where(gamma==1)[0]
    sigmaVec=np.ones(shape)*sigma
    sigmaVec[indices]*=C
    #print sigmaVec
    return sigmaVec

def analyze(X,y_,initNum,stepNum=1000,outputFile=None,plot=False,v=0.5,multi=False,chains=4,C_sigma=2,name=None,indices=None,features=None):   
    #Preprocess data for Modeling
    #betaCovars=['{}_beta'.format(i) for i in X.columns]
    #gammaCovars=['{}_gamma'.format(i) for i in X.columns]
    #covarNames=gammaCovars+betaCovars+['alpha']
    #features=X.columns[1:] 
    #X=np.matrix(X)[:,1:]
    #print X
    m=X.shape[1]
    if indices==None:
        indices=range(m) # passing on the indices of the features we are analyzing. this is for the summary report 
    sha_m=shared(m)
    y_=np.array(y_)
    initNum=np.array(initNum)
    shA_X = shared(X)
    #print "data created"
    #Generate Model
    #graphMatrix=featureListToGraphMatrix(features)
    linear_model = pm.Model()
    #print m
    with linear_model:
        # Priors for unknown evoModel parameters
        #print y_.mean(),y_.std()
        obs=np.floor((y_*initNum))
        #defining mean
        alpha = pm.Normal("alpha",mu=sc.special.logit(y_.mean()),sd=sc.special.logit(y_.std()))
        #gamma=gp.GammaaPrior("gamma",v=v,graphMatrix=graphMatrix)
        gamma=pm.Bernoulli("gamma",p=v,shape=m)
        #defining the covaricnce matrix
        sha_C=shared(C_sigma)
        sigma =pm.Uniform('sigma',lower=0.2,upper=1.5)
        nu = pm.Uniform('nu', 0, 5)
        C_triu = pm.LKJCorr('C_triu', nu, 2) 
        C = pm.Deterministic('C', T.fill_diagonal(C_triu[np.zeros((m, m), dtype=np.int64)], 1.))
        #print type(gamma),type(sigma),type(X.shape[1])
        GammaPriorVariance=GammaCovariance(gamma, sigma,sha_m,sha_C)
        sigma_diag = pm.Deterministic('sigma_mat', T.nlinalg.diag(GammaPriorVariance))
        cov = pm.Deterministic('cov', T.nlinalg.matrix_dot(sigma_diag, C, sigma_diag))            
        #sigma = pm.InverseGamma("sigma",alpha=0.05,beta=0.05)#http://www.stat.columbia.edu/~gelman/research/published/taumain.pdf
        betas = pm.MvNormal("betas",0,
                                  cov=cov,shape=m)
        # Expected value of outcome
        mu = exp_it(alpha+np.array([betas[j]*gamma[j]*shA_X[:,j] for j in range(m)]).sum())   
        likelihood =pm.Binomial("likelihood", n=initNum, p=mu, observed=obs)

        stepC=Metropolis([betas,gamma,nu,C_triu,sigma,alpha],sigmaFactor=2,scaling=0.5)#
        #stepD=pm.BinaryGibbsMetropolis([gamma])#,scalingFunc=coolDownFunc
        #step=
        # MCMC
        tic=time.time()
        #trace = merge_traces([pm.sample(stepNum/5,tune=stepNum/25.,step=[stepC],progressbar=True,njobs=1,chain=i) for i in range(5)])
        trace = pm.sample(stepNum,tune=stepNum/5.,step=[stepC],progressbar=True,njobs=1)
        tac=time.time()
        print tac-tic
        summ=pm.stats.df_summary(trace)
    if not(features is None):
        #the point is to make the summary more readable by changing the index of the run to the index of 
        print (summ.index.values)
        variables=summ.index.values #taking the values of the variables
        for i in range(len(variables)):
            if "betas" in variables[i]: #if it's a beta coeff
                num=int(variables[i].split("__")[1]) #take the number
                variables[i]=features[indices[num]]#we are running on the indices. and then taking the relevent features from the list
        summ.index=variables#placing back the index
    #dfSummary=pm.stats.df_summary(trace)
    #Traceplot
    if plot:
        pm.traceplot(trace,varnames=['alpha','betas','sigma'],combined=True)
        plt.title(name)
        plt.show()
    else:
        pm.traceplot(trace,varnames=['alpha','betas','sigma'],combined=True)
        plt.title(name)
        plt.savefig('/sternadi/home/volume1/guyling/MCMC/dataSimulations/plots/{}.png'.format(name))
    #print summ
    return summ 
"""c 
Iterating over the features and outputing the indices of the features that are alpha=0.95 influential to the 
data == 95% of the postirior dist. the beta is>0 or 95%<0. This veresion is taking pairs of features and checking them
and is designed for cases in which the feature number is very big > 15
"""

def featureSelectionPairs(X,y_,initNum,features,featureIndices=[],stepNum=10000):
    sigFeatures=[]
    if len(featureIndices)==0 :
        print "ok"
        featureIndices=range(0,X.shape[1]-1)
    #print features
    X=np.matrix(X)
    X=X[:,1:]#removing 1 col. temp work around
    #print X
    for i in range(0,len(featureIndices),2):
        #print i
        if i==len(featureIndices)-1:#end case in which the len is odd and then we want to take the last one with the prev.
            coupleIndices=featureIndices[i-1:i+1]
            print features[featureIndices[i-1]],features[featureIndices[i]]
        else:    
            print features[featureIndices[i]],features[featureIndices[i+1]]
            coupleIndices=featureIndices[i:i+2]
            featureCouple=X[:,coupleIndices]
            res=analyze(featureCouple, y_, initNum,stepNum)
            print res.loc['betas__0'],res.loc['betas__1']
            if res.loc['betas__0']['hpd_2.5']*res.loc['betas__0']['hpd_97.5']>0.7: #if they are both more than 0 or both less than 0
                sigFeatures.append(featureIndices[i])
                print "{} is a sigFeature at mean {}".format(features[featureIndices[i]],res.loc['betas__0']['mean'])
            else:
                print "{} is a  not a sigFeature at mean {}".format(features[featureIndices[i]],res.loc['betas__0']['mean'])
            if res.loc['betas__1']['hpd_2.5']*res.loc['betas__1']['hpd_97.5']>0.7:
                sigFeatures.append(featureIndices[i+1]) 
                print "{} is a sigFeature at mean {}".format(features[featureIndices[i+1]],res.loc['betas__1']['mean'])
            else:
                print "{} is not a sigFeature at mean {}".format(features[featureIndices[i+1]],res.loc['betas__1']['mean'])
    return sigFeatures

def featureSelection(X,y_,initNum,features,featureIndices,stepNum=100000,name=None):
    sigFeatures=[]
    X=np.matrix(X)
    print name
    stepNum=10**4*len(featureIndices)#the ammount of steps is prop. to the ammount of featuresh
    X=X[:,1:]#removing 1 col. temp work around
    onlyFeatures=X[:,featureIndices]# taking only the features we want to check if are sig.
    res=analyze(onlyFeatures, y_, initNum,stepNum,plot=False,v=1./max(len(featureIndices),25),name=name,indices=featureIndices,features=features)#the v parameter is default to be 1/the number of features. Meaning assuming only 1 feature 
    #print res.head()    
    for i in range(len(featureIndices)):
        if np.abs(res.loc['betas__{}'.format(i)]['mean'])>0.7:
            if res.loc['betas__{}'.format(i)]['sd']<0.4:
                sigFeatures.append((featureIndices[i],res.loc['betas__{}'.format(i)]['mean'])) #inserting the feature name and the coeff
                print "{} is a sigFeature at mean {}".format(features[featureIndices[i]],res.loc['betas__{}'.format(i)]['mean'])
        else:
                print "{} is a sigFeature at mean {}".format(features[featureIndices[i]],res.loc['betas__{}'.format(i)]['mean'])
    return sigFeatures,res    
    
"""
Many times in the data you get false positives for sig. features cause features are correlated with each other.
For example if the real feature is P2_P4_CC the P2_C and P4_C will be sig. creating a grpah where V=features e_ij in E iff
i part of j
"""


            
                    
            
    
def analyzeSigFeatuesByCorrelation(X,y_,initNum,sigFeatures,features):
    X=np.matrix(X)
    X=X[:,1:]#walk around, the 0th row is not a feature.
    realSigFeatures=Set(sigFeatures)
    graphMatrix=featureListToGraphMatrix(sigFeatures, features)   
    for s in sigFeatures:
        correlatedFeatures=list(np.nonzero(graphMatrix[s])[0])+[s]
        #print type(correlatedFeatures)
        print "the feature is {} and correlatedFeatures are {}".format(s,correlatedFeatures)
        if len(correlatedFeatures)>1: # not only idenityty (there is t least one corr. feature)
            print (Set(correlatedFeatures)^Set(featureSelection(X, y_, initNum, features,correlatedFeatures)))
            realSigFeatures=realSigFeatures - (Set(correlatedFeatures)^Set(featureSelection(X, y_, initNum, features,correlatedFeatures)))
        else:
            continue
    
    realSigFeatures=list(realSigFeatures)
    print realSigFeatures
    X=X[:,realSigFeatures]
    return analyze(X, y_, initNum,stepNum=50000,plot=True)
    
    


def main():
    """
    size=256
    initNum=np.random.poisson(lam=1000000,size=size)
    data=pd.read_csv(r'C:\Users\Guyling1\ContextProject\SABIN\regressionCSVOneFileCountsK5P7.csv')
    data=f.getDummyVarPositionForK(data,1)
    data=f.cleanDataMatrix(data,withNull=True)
    covariates=f.createCovarsAndProbVectorForMutationType(data,'CT')[0]
    covariates=np.matrix(covariates)[:,1:]
    pinit=np.ones(size)*-11
    X=covariates
    #print covarMatrix.shape
    realBeta=np.random.multivariate_normal(np.zeros(covars),np.identity(covars))#
    realBeta[0]=3
    realBeta[1]=-3
    #realBeta[-1]=-5
    realGamma=np.zeros(covars)
    realGamma[0]=1
    realGamma[1]=1
    #realGamma[-1]=1
    betaGamma=realBeta*realGamma
    print X
    means=pinit+np.dot(X,betaGamma)
    means=sc.special.expit(means)
    #print means[:100]
    y_=np.random.binomial(initNum,means)
    y_=y_/initNum
    #print y_
    #print y_[0:100]
    #Split Data
    #X_tr, X_te, y_tr, y_te,initNum_tr,initNum_te = train_test_split(X,y_,initNum,test_size=0.25, random_state=0)
    
    """
    #location=r'C:\Users\Guyling1\Documents\guyReserch\\'
    location='/sternadi/home/volume1/guyling/MCMC/'
    data=pd.read_csv(location+'moranFreqsTry.csv')
    data=f.getDummyVarPositionForK(data,1)
    data=f.cleanDataMatrix(data,withNull=True)
    X,y_,initNum,D=f.createCovarsAndProbVectorForMutationType(data,'CT')
    #D.to_csv(r'C:\Users\Guyling1\Documents\guyReserch\testFeatures.csv')
    #outputFile='C:\Users\Guyling1\ContextProject\MCMC\P7K5_{}.csv'.format(mt)
    features=X.columns[1:]
    print len(features)
    #print features
    sigList=[3, 8, 9, 10, 28, 32, 33, 34, 36, 37, 38, 40, 41, 42, 46, 62, 63, 65, 66, 69, 73, 74, 78, 99, 100, 101, 103, 104, 105, 106, 107]
    #sigList=featureSelectionPairs(X, y_, initNum, features) 
    #print sigList
    #featureListToGraphMatrix(sigList, features)
    #featureSelection(X, y_, initNum, features, sigList )
    #analyzeSigFeatuesByCorrelation(X, y_, initNum, sigList, features)
    analyze(X, y_, initNum, stepNum=400000,plot=False,v=0.1)
#main()