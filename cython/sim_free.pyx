import numpy as np
cimport numpy as np
    
cdef np.ndarray randRot(double d):
    '''Get randomized rotation matrix, using Arvo 1992'''
    cdef double th = d*np.random.rand()*np.pi/2
    cdef double phi = np.random.rand()*np.pi/2
    cdef double z = d*np.random.rand()*2
    cdef double r = np.sqrt(z)
    cdef np.ndarray V = np.array([np.sin(phi)*r,np.cos(phi)*r,np.sqrt(2-z)])
    cdef np.ndarray R = np.array([[np.cos(th),-np.sin(th),0],[np.sin(th),np.cos(th),0],[0,0,1]])
    return np.dot(np.outer(V,V)-np.identity(3),R)

    
def initConfig(dnaParams):
    '''Initialize DNA configuration using given params'''
    # one row per bead x 3 columns
    N = dnaParams['N']
    X = np.zeros((N,3))
    # they are evenly spaced in z, starting at z=0
    # a is inter-bead spacing, or L/(N-1)
    for i in range(N):
        X[i,2] = i*dnaParams['a']
    return X

cdef double Upartial(np.ndarray[double,ndim=2] X, int i, np.ndarray[double,ndim=1] dnaParams):
    '''Return component of energy associated with ith bead only'''
    # scale down by 0.5 to avoid double-counting pairs
    
    # first calculate stretching and bending
    cdef double U = 0
    cdef double a = dnaParams[1]
    cdef double n1 = 0
    cdef double n2 = 0

    if (i>0):
        d1 = X[i,:]-X[i-1,:]
        n1 = np.sqrt(d1.dot(d1))
        U += 0.5*dnaParams[2]*(n1-a)**2
    if (i<dnaParams[0]-1):
        d2 = X[i+1,:]-X[i,:]
        n2 = np.sqrt(d2.dot(d2))
        U += 0.5*dnaParams[2]*(n2-a)**2
    if (i>0) and (i<dnaParams[0]-1): 
        U += 0.5*dnaParams[3]*np.arccos(np.dot(d1,d2)/(n1*n2))**2
    # this si a start
    return U
    
cdef double Utotal(np.ndarray[double,ndim=2] X, np.ndarray[double,ndim=1] dnaParams):
    '''Return total energy'''
    cdef double U = 0
    for i in range(X.shape[0]):
        U += 0.5*Upartial(X,i,dnaParams)
    return U
    
cdef np.ndarray moveGlobal(np.ndarray X, int i, double dth):
    '''Make a global move, pick pivot and rotate'''
    Xnew = np.copy(X)
    # make a random small-angle rotation matrix
    R = randRot(dth)
    
    if np.random.rand() < 0.5:
        Xnew[:i,:] = Xnew[i,:]+np.dot(Xnew[:i,:]-Xnew[i,:],R)
    else:
        Xnew[(i+1):,:] = Xnew[i,:]+np.dot(Xnew[(i+1):,:]-Xnew[i,:],R)
    return Xnew
    
    
def sampleConfigs(nsamp, burnin, thin, delta, dnaParams):
    '''Sample equilibruim configurations of DNA'''
    # make big 3d array
    cdef int N = dnaParams['N']
    cdef np.ndarray Xs = np.zeros((nsamp, N, 3))
    # and initialize current config
    cdef np.ndarray X = initConfig(dnaParams)
    # now loop through and sample configuration space
    cdef int curind = 0
    # number since last one stored
    cdef int nsince = -burnin*thin*N
    cdef int nacc = 0
    cdef int ntot = 0
    
    cdef double U0 = 0
    cdef double U1 = 0
    
    dp = np.array([dnaParams['N'],dnaParams['a'],dnaParams['ks'],dnaParams['kb']])
    
    while True:
        # propose new configuration
        # pick random bead to perturb
        i = np.random.randint(N)
        Xnew = np.copy(X)
            
        if np.random.rand() < 0.7:
            # and perturb it
            Xnew[i,:] += 2*delta*(np.random.rand(3)-0.5)
            # now calculate and compare energies
            U0 = Upartial(X,i,dp)
            U1 = Upartial(Xnew,i,dp)
        else:
            # make global move instead
            Xnew = moveGlobal(X,i,np.sqrt(2/dp[3]))
            U0 = Utotal(X,dp)
            U1 = Utotal(Xnew,dp)
        
        if np.random.rand() < np.exp(U0-U1):
            # accept proposal
            X = Xnew
            nacc += 1    

        # do we savae current conformation or not?
        if nsince >= thin*N:
            Xs[curind,:,:] = X
            curind += 1
            nsince = 0
            print curind
        
        ntot += 1
        nsince += 1
        if curind == nsamp:
            break
    
    print 'Accept ratio: ' + str(nacc/float(ntot))
    
    return Xs