import pandas as pd
import numpy as np
import myio as Mio
import mymath as Mmath

def getAndStoreDump(simname,increment,rootdir,dumpdir='Dump_Files/',pickledir='Pickle_Files/',override=False,bounds=None,**kwargs):
    return Mio.getAndStore(readDumpFile,getDumpFilename,override=override,subdirstore=rootdir+pickledir,simname=simname,increment=increment,subdirread=rootdir+dumpdir,bounds=bounds)

def readDumpFile(simname,increment,subdirread,bounds,**kwargs):
    filename = subdirread + getDumpFilename(simname,increment) + '.dump'
    headerlines = getHeaderLinesDumpFile(filename)
    boxdims = getBoxDimDumpFile(filename)
    res = cutDownDumpFile(filename,headerlines,bounds,boxdims)
    return {'boxdims': boxdims, 'array': unscaleCoords(res,boxdims)}

def cutDownDumpFile(filename,headerlines,bounds,boxdims,indices=range(2,5),chunksize=1000000):
    if bounds is not None:
        boundsscaled = scaleBounds(bounds,boxdims)
    reader = pd.read_csv(filename,sep=' ',skiprows=headerlines,iterator=True,chunksize=chunksize)
    chunkall = None
    for chunk in reader:
        values = chunk.values
        indexall = np.ones((values.shape[0],), dtype=bool)
        if bounds is not None:
            for [posmin, posmax], index in zip(boundsscaled,indices):
                valuescurr = values[:,index]
                indexall = indexall & (posmin <= valuescurr) & (valuescurr <= posmax)
        chunknew = values[indexall,:]
        try:
            chunkall = np.concatenate((chunkall,chunknew),axis=0)
        except ValueError:
            chunkall = chunknew
    return chunkall

def unscaleCoords(dumparray,boxdims,indexstart=2):
    for i in range(3):
        dumparray[:,i+indexstart] = Mmath.rescaleCoords(dumparray[:,i+indexstart],[0,1],boxdims[i,:])
    return dumparray
    
def scaleBounds(bounds,boxdims):
    boundsnew = np.empty(np.shape(bounds))
    for i in range(3): # dimensions
        boundsnew[i,:] = Mmath.rescaleCoords(bounds[i,:],boxdims[i,:],[0,1]) # scale bounds
    return boundsnew

def getDumpFilename(simname,increment,**kwargs):
    return simname + '.' + str(increment)
    
def getBoxDimDumpFile(filename):
    lines, _ = Mio.readFileForKey(filename,'ITEM: BOX',4)
    return np.loadtxt(lines[1:])
                
def getHeaderLinesDumpFile(filename):
    _, linenum = Mio.readFileForKey(filename,'ITEM: ATOMS')
    return linenum

def getBoxDimLogFile(filename): # get box dimensions from lammps log file
    lines, _ = Mio.readFileForKey(filename,'Created orthogonal box')
    line = lines[0]
    indicesstart = Mio.findCharIndices(line,'(')
    indicesend = Mio.findCharIndices(line,')')
    bounds = [line[istart+1:iend] for istart, iend in zip(indicesstart,indicesend)]
    return np.transpose(np.loadtxt(bounds))

def getOptimizedConfig(filename): # just get the last line of the MS run (output: numpy array)
    blockall = []
    with open(filename,'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('Loop time'): # stop reading, store results
                blockall.append(lineprev)
            lineprev = line
    return np.loadtxt(blockall)
        
def readLogFile(filename): # get all lines of MS/MD run (output: list of numpy arrays)
    blockall, blockoflines = [], []
    with open(filename,'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('Step'): # start reading, get key, clear previous results
                blockoflines = []
            elif line.startswith('Loop time'): # stop reading, store results
                array = np.loadtxt(blockoflines) # convert to numerical array
                blockall.append(array)
            else:
                blockoflines.append(line)
    return blockall
    
# obsolete    
def readLogFileObs(filename): # reads in all numeric lines from log file
    return Mio.readBlockFile(filename,[],[],fortranoption=False,keyfun=getKeyLogFileObs)
    
def getKeyLogFileObs(line,keystarts,keyends,keyold):
    if 'Step' not in line: # to remove some spurious chunks found by readBlockFile
        return None
    else:
        return keyold + 1