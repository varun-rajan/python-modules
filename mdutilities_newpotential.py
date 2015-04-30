import numpy as np
import mymath as Mmath
import scipy.optimize as spo

r0 = 2**(1/6) # equilibrium distance for nearest neighbor
npoints = 3000
rmin, rmax = 0.1, 2.0 # range for tabulated potential
rvecmaster = np.linspace(rmin,rmax,npoints)

r1bartry = 1.05 # use morse to 5% strain
r1bartrywarner = 1.1 # use morse to 10% strain
r2bardict = {2: np.sqrt(2), 3: np.sqrt(3/2)}
r3bardict = {2: np.sqrt(3)-0.02, 3: np.sqrt(2)-0.02}
alphabardict = {2: 3.5*r0, 3: 6.5*r0}
rhodict = {2: 2/(np.sqrt(3)*r0**2), 3: np.sqrt(2)/r0**3}
keys2d = ['brittle','inter5','inter4','inter3','inter2','inter1','ductile']
keys3d = ['brittle3d','inter43d','inter33d','inter23d','inter13d','ductile3d']
epvec2d = [0.102,0.160,0.214,0.264,0.308,0.346,0.378]
epvec3d = [0.009,0.077,0.143,0.207,0.266,0.318]
usfvec2d = [0.6333,0.5265,0.4309,0.3422,0.2641,0.1966,0.1397]
usfvec3d = [0.8818,0.7344,0.6024,0.4816,0.3758,0.2824]
epvecdict = {2: dict(zip(keys2d,epvec2d)), 3: dict(zip(keys3d,epvec3d))}
usfvecdict = {2: dict(zip(keys2d,usfvec2d)), 3: dict(zip(keys3d,usfvec2d))}

def genPotential1(dimensions,potentialname,extend=True): # thin wrapper, using default values
    ep = epvecdict[dimensions][potentialname]
    return genPotential(dimensions,ep,extend=extend)
    
def genPotential2(dimensions,usftarget,extend=True): # thin wrapper, attempting to generate potential with desired usf
    ep = estimateEp(usf,dimensions)
    return genPotential(dimensions,ep,extend=extend)  
 
def genPotential(dimensions,ep,rvec=rvecmaster,style='New',extend=False,**kwargs):
    if style == 'Old':
        paramsfun, energyfun = solveEqnsOld, getEnergyOld
    elif style == 'New':
        paramsfun, energyfun = solveEqnsNew, getEnergyNew
    paramsdict = paramsfun(dimensions,ep,**kwargs) # pass kwargs to change defaults of solveEqns (e.g. offset, alphabar, etc.)
    energyvec = np.array([energyfun(rbar,**paramsdict) for rbar in rvec/r0])
    forcevec = Mmath.getDiffVec(rvec,-energyvec)
    if extend:
        rvec, energyvec, forcevec = extendToZero(rvec,energyvec,forcevec)
    numvec = np.arange(np.shape(forcevec)[0])+1
    return np.column_stack((numvec,rvec,energyvec,forcevec)), paramsdict
    
def extendToZero(rvec,energyvec,forcevec,tol=1.e-8):
    rvecnew = np.insert(rvec,0,tol)
    forcevecnew = np.insert(forcevec,0,forcevec[0])
    energy0 = forcevec[0]*(rvec[0]-tol) + energyvec[0] # linear extrapolation
    energyvecnew = np.insert(energyvec,0,energy0)
    return rvecnew, energyvecnew, forcevecnew
    
def solveEqnsNew(dimensions,ep,r1bar=r1bartry,r2bar=None,r3bar=None,alphabar=None,offset=0,step=0.01,tol=1.e-8,r3barfac=0.98,n1=3,n2=3):
    # largest cutoff radius possible is r3bar
    # smallest cutoff radius possible is r2bar
    # tries to solve for largest cutoff radius possible
    # constructs smooth potential so that value of potential at r2bar is -ep
    # n1, n2 are order of splines
    if r2bar is None:
        r2bar = r2bardict[dimensions]
    if r3bar is None:
        r3bar = r3bardict[dimensions]
    if alphabar is None:
        alphabar = alphabardict[dimensions]
    while r3barfac > 0:
        r3barcurr = r2bar + (r3bar - r2bar)*r3barfac
        func = lambda coeff: solveEqnsSub(coeff,ep,r1bar,r2bar,r3barcurr,alphabar,offset,n1,n2)
        coeffsol = spo.fsolve(func,np.zeros(n1+n2+2))
        spline1, spline2 = coeffsol[0:n1+1], coeffsol[n1+1:n1+n2+2]
        maxval2 = Mmath.maxPolyRoot(spline2)
        if maxval2 > tol:
            r3barfac = r3barfac - step
        else:
            keyslist = ['spline1','spline2','r1bar','r2bar','r3bar','alphabar','offset','ep']
            valslist = [spline1,spline2,r1bar,r2bar,r3barcurr,alphabar,offset,ep]
            return dict(zip(keyslist,valslist))
    raise ValueError('Bad input')
    
def solveEqnsSub(coeff,ep,r1bar,r2bar,r3bar,alphabar,offset,n1,n2): # if n1, n2 are not 3, 3, need to add or eliminate conditions below
    spline1, spline2 = coeff[0:n1+1], coeff[n1+1:n1+n2+2]
    dspline1, dspline2 = np.polyder(spline1), np.polyder(spline2)
    ddspline1, ddspline2 = np.polyder(dspline1), np.polyder(dspline2)
    energy1 = morse(r1bar,alphabar,offset)
    denergy1 = morseD(r1bar,alphabar)
    f = np.zeros(n1+n2+2)
    f[0] = np.polyval(spline1,r1bar) - energy1
    f[1] = np.polyval(dspline1,r1bar) - denergy1
    f[2] = np.polyval(spline2,r3bar)
    f[3] = np.polyval(dspline2,r3bar)
    f[4] = np.polyval(spline1,r2bar) - (-ep)
    f[5] = np.polyval(spline2,r2bar) - (-ep)
    f[6] = np.polyval(dspline1,r2bar) - np.polyval(dspline2,r2bar)
    f[7] = np.polyval(ddspline1,r2bar) - np.polyval(ddspline2,r2bar) # unnecessary, but nice to have
    return np.array(f)
    
def getEnergyNew(rbar,spline1,spline2,r1bar,r2bar,r3bar,alphabar,offset,**kwargs):
    if rbar < r1bar: # use morse
        return morse(rbar,alphabar,offset)
    elif rbar < r2bar: # use spline1
        return np.polyval(spline1,rbar)
    elif rbar < r3bar: # use spline2
        return np.polyval(spline2,rbar)
    else:
        return 0
    
def writeToFile(potentialname,data,paramsdict,filename=None,writeoption='w',writeformat='%d %10.8f %10.8f %10.8f'):  # use 'a' for append, 'w' for (over)write
    if filename is None:
        filename = potentialname
    with open(filename,writeoption) as f:
        ndatapoints = np.shape(data)[0]
        writeHeader(potentialname,ndatapoints,f,paramsdict)
    with open(filename,'ab') as fb:
        np.savetxt(fb,data,fmt=writeformat)
        
def writeHeader(potentialname,ndatapoints,f,paramsdict):
    f.write('\n\n') # pad with blank lines for readability if multiple potentials
    writeParams(f,**paramsdict)
    f.write(potentialname + '\n')
    f.write('N %d \n \n' % ndatapoints)

def writeParams(f,ep=np.nan,r1bar=np.nan,r2bar=np.nan,r3bar=np.nan,alphabar=np.nan,offset=np.nan,**kwargs):
    f.write('# Tabulated potential, ep = %f, r1 = %f, r2 = %f, r3 = %f, alpha = %f, offset = %f \n \n' % (ep,r1bar*r0,r2bar*r0,r3bar*r0,alphabar/r0,offset))   

def morse(rbar,alphabar,offset=0):
    return (np.exp(alphabar*(1 - rbar)) - 1)**2 - 1 - offset
    
def morseD(rbar,alphabar):
    return -2*alphabar*np.exp(alphabar*(1 - rbar))*(np.exp(alphabar*(1 - rbar)) - 1)

# analogous functions for Warner potential

def solveEqnsOld(dimensions,ep=None,r1bar=r1bartrywarner,r3bar=None,alphabar=None,offset=0,**kwargs): # as r3bar (cutoff) increases, potential changes from brittle to ductile; brittle - r3bar = r2bartry; ductile - r3bar = 0.99*r3bartry (values outside bounds lead to bad potentials)
    if r3bar is None:
        r3bar = r3bardict[dimensions]
    if alphabar is None:
        alphabar = alphabardict[dimensions]
    func = lambda coeff: solveEqnsSubOld(coeff,r1bar,r3bar,alphabar,offset)
    coeffsol = spo.fsolve(func,np.zeros(4))
    keyslist = ['spline1','r1bar','r3bar','alphabar','offset']
    valslist = [coeffsol,r1bar,r3bar,alphabar,offset]
    return dict(zip(keyslist,valslist))
    
def solveEqnsSubOld(spline1,r1bar,r3bar,alphabar,offset):
    dspline1 = np.polyder(spline1)
    ddspline1 = np.polyder(dspline1)
    energy1 = morse(r1bar,alphabar,offset)
    denergy1 = morseD(r1bar,alphabar)
    f1 = np.polyval(spline1,r1bar) - energy1
    f2 = np.polyval(dspline1,r1bar) - denergy1
    f3 = np.polyval(spline1,r3bar)
    f4 = np.polyval(dspline1,r3bar)
    return np.array([f1,f2,f3,f4])
    
def getEnergyOld(rbar,spline1,r1bar,r3bar,alphabar,offset): # old method (Warner potential)
    if rbar < r1bar: # use morse
        return morse(rbar,alphabar,offset)
    elif rbar < r3bar: # use spline
        return np.polyval(spline1,rbar)
    else:
        return 0
        
# misc

areadict = {2: r0, 3: np.sqrt(3)/2*r0**2}
offsetdict = {2: -0.078, 3: -0.052} # for estimating usf (Equation 5)

def estimateEp(usf,dimensions):        
    offset = offsetdict[dimensions]
    area = areadict[dimensions]
    return -((usf - offset)*area - 1)/2
    
def estimateUsf(ep,dimensions): # takes ep as being positive
    offset = offsetdict[dimensions]
    area = areadict[dimensions]
    return (1 - 2*np.abs(ep))/area + offset 

# summary of properties of potentials
def getPotentialProps(potential,lattice):
    dimensiondict = {'hex': 2, 'fcc': 3, 'hcp': 3}
    dimensions = dimensiondict[lattice]
    alpha = alphabardict[dimensions]/r0
    usf = usfvecdict[dimensions][potential]
    ep = epvecdict[dimensions][potential]
    rho = rhodict[dimensions]
    r2 = r2bardict[dimensions]*r0
    r1 = r1bartry*r0
    # other properties
    if lattice == 'hex':
        symmetry = 'cubic'
        C11 = 3/2*np.sqrt(3)*alpha**2
        Velastic = {'11': C11, '12': C11/3, '44': C11/3}
        Vsurface = {'112': {'surf': 1/r0}}
    elif lattice == 'fcc':
        symmetry = 'cubic'
        C11 = 2*np.sqrt(2)*alpha**2/r0
        Velastic = {'11': C11, '12': C11/2, '44': C11/2}
        Vsurface = {'111': {'surf': np.sqrt(3)/r0**2}, '001': {'surf': 2/r0**2}, '011': {'surf': np.sqrt(9/2)/r0**2}}
    elif lattice == 'hcp':
        symmetry = 'hexagonal'
        C11 = 5/np.sqrt(2)*alpha**2/r0
        Velastic = {'11': C11, '33': C11*16/15, '13': C11*4/15, '44': C11*4/15, '66': C11/3}
        # Velastic = {'11': 128.6, '33': 141.95, '13': 35.49, '44': 35.49, '66': 39.92}
        Vsurface = {'001': {'surf': 1.458}}
    return {'elastic': Velastic, 'surface': Vsurface, 'unstable': {'full': usf, 'partial': usf}, 'symmetry': symmetry, 'rho': rho, 'r1': r1, 'r2': r2, 'ep': ep}

# obsolete stuff
# 2d
# smoothed potential, cutoff -0.05 (new potential #2)
# epsvecnew22d = [0.113,0.169,0.218,0.262,0.301,0.348,0.39]
# usfvecnew22d = [0.6121,0.5067,0.4116,0.3261,0.2504,0.1844,0.1479]
# original potential
# epsvecold2d = [0.1,0.15,0.22,0.27,0.31,0.35,0.38]
# usfvecold2d = [0.6371,0.5087,0.4211,0.3334,0.2631,0.1927,0.1399]
# 3d
# potential with incorrect r2bar = sqrt(17/12)
# epsvecnew23d = [0.01,0.17,0.247,0.313,0.372,0.422]
# usfvecnew23d = [0.9165,0.7667,0.6299,0.5061,0.3952,0.3013]
    