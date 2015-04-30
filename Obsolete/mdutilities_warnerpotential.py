import numpy as np
import mymath as Mmath
import scipy.optimize as spo
import mydictionaries as Mdict

r0 = 2**(1/6) # equilibrium distance for nearest neighbor
npoints = 3000
rmin, rmax = 0.1, 2.0 # range for tabulated potential
rvecmaster = np.linspace(rmin,rmax,npoints)

r1bartry = 1.05 # use morse to 5% strain
r1bartrywarner = 1.1 # use morse to 10% strain
r2bartry = np.sqrt(2) # second-nearest neighbor distance in USF configuration
r3bartry = np.sqrt(3)-0.02 # second-nearest neighbor distance in equilibrium configuration, avoiding compression issue
alphabartry = 3.5*r0

# summary of properties for brittle, ductile potentials (from simple lammps tests)
# these are not necessary for functioning of this module, but can be used by other modules
keys = ['brittle','inter5','inter4','inter3','inter2','inter1','ductile']
# cutoff -0.02 (new potential #1)
epsvecnew1 = [0.102,0.160,0.214,0.264,0.308,0.346,0.378]
usfvecnew1 = [0.6333,0.5265,0.4309,0.3422,0.2641,0.1966,0.1397]
# smoothed potential, cutoff -0.05 (new potential #2)
epsvecnew2 = [0.113,0.169,0.218,0.262,0.301,0.348,0.39]
usfvecnew2 = [0.6121,0.5067,0.4116,0.3261,0.2504,0.1844,0.1479]
# original potential
epsvecold = [0.1,0.15,0.22,0.27,0.31,0.35,0.38]
usfvecold = [0.6371,0.5087,0.4211,0.3334,0.2631,0.1927,0.1399]
epsvecdict = {'new1': epsvecnew1, 'new2': epsvecnew2, 'old': epsvecold}
usfvecdict = {'new1': usfvecnew1, 'new2': usfvecnew2, 'old': usfvecold}
# surface energy
surf = 1.0/r0
 
def getAll(rvec=rvecmaster,style='New',extend=False,**kwargs): # pass ep for New Potential; pass r3bar for Warner
    if style == 'Warner': # old
        paramsfun, energyfun = solveEqnsWarner, getEnergyWarner
    elif style == 'New':
        paramsfun, energyfun = solveEqnsNew, getEnergyNew
    paramsdict = paramsfun(**kwargs)
    energyvec = np.array(getEnergyLoop(rvec/r0,paramsdict,energyfun))
    forcevec = Mmath.getDiffVec(rvec,-energyvec)
    if extend:
        rvec, energyvec, forcevec = extendToZero(rvec,energyvec,forcevec)
    numvec = np.arange(np.shape(forcevec)[0])+1
    return np.column_stack((numvec,rvec,energyvec,forcevec)), paramsdict
    
def getEnergyLoop(rbarvec,paramsdict,fun):
    return [fun(rbar,**paramsdict) for rbar in rbarvec]

def extendToZero(rvec,energyvec,forcevec,tol=1.e-8):
    rvecnew = np.insert(rvec,0,tol)
    forcevecnew = np.insert(forcevec,0,forcevec[0])
    energy0 = forcevec[0]*(rvec[0]-tol) + energyvec[0] # linear extrapolation
    energyvecnew = np.insert(energyvec,0,energy0)
    return rvecnew, energyvecnew, forcevecnew
    
def solveEqnsNew(ep,r1bar=r1bartry,r2bar=r2bartry,r3bar=r3bartry,alphabar=alphabartry,offset=0,step=0.01,tol=1.e-8,r3barfac=0.98,n1=3,n2=3):
    # largest cutoff radius possible is r3bar
    # smallest cutoff radius possible is r2bar
    # tries to solve for largest cutoff radius possible
    # constructs smooth potential so that value of potential at r2bar is -ep
    # n1, n2 are order of splines
    while r3barfac > 0:
        r3barcurr = r2bar + (r3bar - r2bar)*r3barfac
        func = lambda coeff: solveEqnsSub(coeff,ep,r1bar,r2bar,r3barcurr,alphabar,offset,n1,n2)
        coeffsol = spo.fsolve(func,np.zeros(n1+n2+2))
        spline1, spline2 = coeffsol[0:n1+1], coeffsol[n1+1:n1+n2+2]
        maxval2 = Mmath.maxPolyRoot(spline2)
        if maxval2 > tol:
            r3barfac = r3barfac - step
        else:
            keyslist = ['spline1','spline2','r1bar','r2bar','r3bar','alphabar','offset']
            valslist = [spline1,spline2,r1bar,r2bar,r3barcurr,alphabar,offset]
            return Mdict.instantiate(keyslist,valslist)
    print('Bad input')
    
def solveEqnsSub(coeff,ep,r1bar,r2bar,r3bar,alphabar,offset,n1,n2):
    spline1, spline2 = coeff[0:n1+1], coeff[n1+1:n1+n2+2]
    dspline1, dspline2 = np.polyder(spline1), np.polyder(spline2)
    ddspline1, ddspline2 = np.polyder(dspline1), np.polyder(dspline2)
    energy1 = morse(r1bar,alphabar,offset)
    denergy1 = morseD(r1bar,alphabar)
    f = np.zeros(n1+n2+2)
    f[0] = np.polyval(spline1,r1bar) - energy1
    f[1] = np.polyval(dspline1,r1bar) - denergy1
    f[2] = np.polyval(spline1,r2bar) - (-ep)
    f[3] = np.polyval(spline2,r2bar) - (-ep)
    f[4] = np.polyval(dspline1,r2bar) - np.polyval(dspline2,r2bar)
    f[5] = np.polyval(spline2,r3bar)
    f[6] = np.polyval(dspline2,r3bar)
    f[7] = np.polyval(ddspline1,r2bar) - np.polyval(ddspline2,r2bar) # unnecessary, but nice to have
    return np.array(f)
    
def getEnergyNew(rbar,spline1,spline2,r1bar,r2bar,r3bar,alphabar,offset):
    if rbar < r1bar: # use morse
        return morse(rbar,alphabar,offset)
    elif rbar < r2bar: # use spline1
        return np.polyval(spline1,rbar)
    elif rbar < r3bar: # use spline2
        return np.polyval(spline2,rbar)
    else:
        return 0
    
def writeToFile(potentialname,data,paramsdict,filename=None,style='New',writeoption='w',writeformat='%d %10.8f %10.8f %10.8f'):  # use 'a' for append, 'w' for (over)write
    if filename is None:
        filename = potentialname
    with open(filename,writeoption) as f:
        ndatapoints = np.shape(data)[0]
        writeHeader(potentialname,ndatapoints,f,paramsdict,style)
    with open(filename,'ab') as fb:
        np.savetxt(fb,data,fmt=writeformat)
        
def writeHeader(potentialname,ndatapoints,f,paramsdict,style):
    f.write('\n\n') # pad with blank lines for readability if multiple potentials
    if style == 'New':
        writeParamsNew(f,**paramsdict)
    elif style == 'Warner':
        writeParamsWarner(f,**paramsdict)
    f.write(potentialname + '\n')
    f.write('N %d \n \n' % ndatapoints)

def writeParamsNew(f,r1bar,r2bar,r3bar,alphabar,offset,**kwargs):
    f.write('# Tabulated potential, r1 = %f, r2 = %f, r3 = %f, offset = %f, alpha = %f \n \n' % (r1bar*r0,r2bar*r0,r3bar*r0,offset,alphabar/r0))   

def morse(rbar,alphabar,offset=0):
    return (np.exp(alphabar*(1 - rbar)) - 1)**2 - 1 - offset
    
def morseD(rbar,alphabar):
    return -2*alphabar*np.exp(alphabar*(1 - rbar))*(np.exp(alphabar*(1 - rbar)) - 1)
    
# analogous functions for Warner potential...

def solveEqnsWarner(r3bar,r1bar=r1bartrywarner,alphabar=alphabartry,offset=0): # as r3bar (cutoff) increases, potential changes from brittle to ductile; brittle - r3bar = r2bartry; ductile - r3bar = 0.99*r3bartry (values outside bounds lead to bad potentials)
    func = lambda coeff: solveEqnsSubWarner(coeff,r1bar,r3bar,alphabar,offset)
    coeffsol = spo.fsolve(func,np.zeros(4))
    keyslist = ['spline1','r1bar','r3bar','alphabar','offset']
    valslist = [coeffsol,r1bar,r3bar,alphabar,offset]
    return Mdict.instantiate(keyslist,valslist)
    
def solveEqnsSubWarner(spline1,r1bar,r3bar,alphabar,offset):
    dspline1 = np.polyder(spline1)
    ddspline1 = np.polyder(dspline1)
    energy1 = morse(r1bar,alphabar,offset)
    denergy1 = morseD(r1bar,alphabar)
    f1 = np.polyval(spline1,r1bar) - energy1
    f2 = np.polyval(dspline1,r1bar) - denergy1
    f3 = np.polyval(spline1,r3bar)
    f4 = np.polyval(dspline1,r3bar)
    return np.array([f1,f2,f3,f4])
    
def getEnergyWarner(rbar,spline1,r1bar,r3bar,alphabar,offset): # old method (Warner potential)
    if rbar < r1bar: # use morse
        return morse(rbar,alphabar,offset)
    elif rbar < r3bar: # use spline
        return np.polyval(spline1,rbar)
    else:
        return 0

def writeParamsWarner(f,r1bar,r3bar,alphabar,offset,**kwargs):
    f.write('# Tabulated potential, r1 = %f, r3 = %f, offset = %f, alpha = %f \n \n' % (r1bar*r0,r3bar*r0,offset,alphabar/r0))        
    