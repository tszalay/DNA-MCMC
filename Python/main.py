import numpy as np
import matplotlib.pyplot as plt
import time

import sim_draw
import sim_free
import sim_stretch
import sim_pore

Lp = 50.0 # nm

# 1 kT/nm =  4.114 pN
pNperkT = 4.114
kTperpN = 1/pNperkT

def dnaParams(N,a):
    params = {}
    params['N'] = N
# paper used h = 100kT/L0^2, close to what I have
    #params['ks'] = 100.0/a**2
    # this one is equivalent to 1000pN elastic modulus
    params['ks'] = 1000.0*kTperpN/a
    params['a'] = a
    params['kb'] = Lp/a
    params['L'] = N*a
    return params
    
def plotCorrs(X,params=None):
    # compute directional correlations
    # first, tangent/direction vectors
    ts = X[:,1:,:] - X[:,:-1,:]
    # that we wish to normalize
    ts = ts / np.sqrt(np.sum(ts**2,2).reshape((ts.shape[0],ts.shape[1],1)))
    # now dot the first into the rest, and take expectation (avg)
    tdots = np.mean(np.sum(ts*ts[:,0,:].reshape((ts.shape[0],1,ts.shape[2])),2),0)
    xs = np.arange(len(tdots))
    if params is not None:
        xs *= params['a']
    plt.plot(xs,tdots)
    # plot expected theoretical value
    if params is not None:
        plt.plot(xs,np.exp(-params['a']*np.arange(len(tdots))/Lp),'r')
    plt.xlabel('Distance [nm]')
    plt.ylabel('Correlation')
    plt.title('Tangent Correlation for a = '+str(params['a']))
    plt.show()
    
def plotExt(X,params=None):
    plt.plot(X[:,-1,2])
    # plot expected theoretical value
    #if params is not None:
    #    plt.plot(np.exp(-params['a']*np.arange(len(tdots))/Lp),'r')
    plt.show()

def forceExtension():
    params = dnaParams(101,10)
    # what force values to use, starts in pN
    Fvs = 5*np.logspace(-3,1,10)
    # convert to kT/nm units
    Fs = Fvs*kTperpN
    # and output array
    Es = []
    
    for f in Fs:
        # get samples
        params['force'] = f
        Xs = sim_stretch.sampleConfigs(nsamp=1000,burnin=150,thin=1,params=params)
        Es.append(np.mean(Xs[:,-1,2]/params['L']))
    
    # WLC approx.
    E0s = np.linspace(0.01,0.99,50)
    F0s = (pNperkT/Lp)*( 1/(4*(1-E0s)**2) + E0s - 0.25 )
    plt.semilogy(Es,Fvs)
    plt.semilogy(E0s,F0s,'r')
    plt.xlabel('Extension')
    plt.ylabel('Force [pN]')
    plt.title('Force-Extension Curves of 1um dsDNA')
    plt.show()
    #print Fvs, Es

def getPorePos(X,params):
    pz = params['pore_z']
    pos = np.zeros(X.shape[0])
    
    for i in range(X.shape[0]):
        minind = np.min(np.where(X[i,:,2]>pz))-1
        s = (pz-X[i,minind,2])/(X[i,minind+1,2]-X[i,minind,2])
        pos[i] = (minind+s)*params['a']
        
    return pos

def poreForces():
    params = dnaParams(151,10)
    params['pore_z'] = 400
    params['pore_r'] = 5

    # what force values to use, starts in pN
    Fvs = np.logspace(0,2,10)
    # convert to kT/nm units
    Fs = Fvs*kTperpN
    # and output array
    Sigs = []
    
    for f in Fs:
        # get samples
        params['force'] = f
        Xs = sim_pore.sampleConfigs(nsamp=1000,burnin=100,thin=1,params=params)
        pos=getPorePos(Xs,params)
        Sigs.append(np.std(pos))
    
    plt.semilogx(Fvs,Sigs)
    plt.xlabel('Force [pN]')
    plt.ylabel('Fluctuation Sigma [nm]')
    plt.title('Force-Fluctuation with d = 400nm')
    plt.show()
    
def poreDists():
    params = dnaParams(151,10)
    params['pore_r'] = 5
    params['force'] = 50*kTperpN

    dists = 0.8*np.logspace(1.5,3,10)
    dists = dists[::-1]
    # and output array
    Sigs = []
    
    for d in dists:
        # get samples
        params['pore_z'] = d
        Xs = sim_pore.sampleConfigs(nsamp=500,burnin=100,thin=1,params=params)
        #sim_draw.drawData(Xs,params)
        pos=getPorePos(Xs,params)
        Sigs.append(np.std(pos))
    
    dists = dists[::-1]
    Sigs = Sigs[::-1]
    
    plt.plot(dists,Sigs)
    plt.xlabel('Pore Distance [nm]')
    plt.ylabel('Fluctuation Sigma [nm]')
    plt.title('Force-Fluctuation with F = 50pN')
    plt.show()


#forceExtension()
#poreForces()
#poreDists()


#sim_draw.drawData(Xs,params)

#plotCorrs(Xs,params)
#plotExt(Xs/params['L'])
params = dnaParams(151,10)
params['force'] = 75*kTperpN
params['pore_r'] = 5
params['pore_z'] = 400
Xs = sim_pore.sampleConfigs(nsamp=1000,burnin=100,thin=1,params=params)
#plotCorrs(Xs,params=params)
#np.save('data/run6.npy',Xs)