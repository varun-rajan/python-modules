import numpy as np
import myio as Mio
import mydictionaries as Mdict

def readDDFile(inputfile,filetype,increment=None,override=False): # e.g. disloc4.in or init_1.plt or displ100k.plt
    return Mio.getAndStore(readAndProcessDDFile,getPickleName,override,inputfile=inputfile,filetype=filetype,increment=increment,id=increment)

def readAndProcessDDFile(inputfile,filetype,increment,subdir='DD_Files/',**kwargs):
    filename = getDDFilename(inputfile,filetype,increment)
    return readDDFileSub(subdir + filename,filetype)	

def getDDFilename(inputfile,filetype,increment):
    extension = '.plt'
    if increment is not None:
        prefixdict = {'disloc': 'disl_', 'stress': 'stress_', 'displ': 'displ'}
        return incrementToFilename(increment,prefixdict[filetype],extension)
    else:
        filenamedict = {'input': inputfile, 'sn': 'init_1.plt', 'strstr': 'strstrain.dat'}
        return filenamedict[filetype]
        
def readDDFileSub(filename,filetype):
    if filetype == 'input':
        start, end = ['/'],[':',',']
    elif filetype == 'sn':
        start, end = ['"'],['"']
    elif filetype == 'disloc':
        start, end = ['"'],['"']
    else:
        start, end = [],[]
    return Mio.readBlockFile(filename,' ',start,end)
    
def getPickleName(inputfile,filetype,id,**kwargs):
    return Mio.getFilePrefix(inputfile) + '_' + filetype + '_' + str(id)

def filenameToIncrement(filename,prefix):
    factor = 1000
    factorsuffix = 'k'
    filenamenew = Mio.getFilePrefix(filename)
    index1 = len(prefix) 
    if filenamenew[-1] == factorsuffix:
        return factor*int(filenamenew[index1:-1])
    else:
        return int(filenamenew[index1:])
        
def incrementToFilename(increment,prefix,extension):
    factor = 1000
    factorsuffix = 'k'
    cutoff = 100000
    suffix = ''
    if increment >= cutoff:
        (quot, rem) = divmod(increment, factor)
        if rem == 0:
            increment = quot
            suffix = 'k'
    return prefix + str(increment) + suffix + extension
    
