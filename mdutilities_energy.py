import numpy as np
import mdutilities_io as mduio
import mymath as Mmath

energyfacdict = {'metal': 1.60217657e-19, 'LJ': 1}
lengthfacdict = {'metal': 1.e-10, 'LJ': 1}
directiondict = {'x': 0, 'y': 1, 'z': 2}

def getBlockForce(filename,option=1,units='metal',dimensions=3,energyindex=2,subdir='Log_Files/'):
    xvec, energyvec, npoints = energySub(filename,energyindex,units,dimensions,subdir)
    if option == 1: # get force curve
        return np.column_stack((xvec,Mmath.getDiffVec(xvec,-energyvec)))

def getSurfaceEnergy(filename,option=1,units='metal',dimensions=3,energyindex=2,subdir='Log_Files/'):
    xvec, energyvec, npoints = energySub(filename,energyindex,units,dimensions,subdir)
    if option == 1:
        return np.column_stack((xvec,energyvec))
    elif option == 2: # get surface energy
        return energyvec[-1]/2 # two surfaces

def getGSF(filename,option=1,units='metal',dimensions=3,energyindex=2,subdir='Log_Files/'):
    xvec, energyvec, npoints = energySub(filename,energyindex,units,dimensions,subdir)
    if option == 1:
        return np.column_stack((xvec,energyvec))
    elif option == 2: # get stacking fault energy
        return energyvec[(npoints+1)//2]
    elif option == 3: # get unstable stacking fault energy
        return np.max(energyvec[0:(npoints+1)//2])

def energySub(filename,energyindex,units,dimensions,subdir):
    data = mduio.getOptimizedConfig(subdir + filename)
    area = getArea(subdir+filename,units,dimensions)
    energyfac = energyfacdict[units]
    energyvec = energyfac*(data[:,energyindex] - data[0,energyindex])/area
    npoints = np.shape(energyvec)[0]
    xvec = np.arange(npoints)
    return xvec, energyvec, npoints

def getStressStrain(filename,stressindex,lengthindex,units='metal',subdir='Log_Files/'):
    data = mduio.getOptimizedConfig(subdir + filename)
    stressfac = energyfacdict[units]/lengthfacdict[units]**3
    strain = (data[:,lengthindex] - data[0,lengthindex])/data[0,lengthindex]
    stress = -data[:,stressindex]*stressfac # -1 because LAMMPS outputs pressures
    return np.column_stack((strain,stress))
        
def getArea(filename,units,dimensions,plane='xz'): # default: xz area
    directions = list(plane) # e.g. convert 'xz' to ['x','z']
    boxdims = mduio.getBoxDimLogFile(filename)
    L = boxdims[:,1] - boxdims[:,0]
    L1 = L[directiondict[directions[0]]]
    lengthfac = lengthfacdict[units]
    if dimensions == 2:
        return L1*lengthfac
    elif dimensions == 3:
        L2 = L[directiondict[directions[1]]]
        return L1*L2*lengthfac**2    