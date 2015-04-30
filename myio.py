import pickle
import scipy.io as spio
import operator as op
import os
import subprocess

def getAndStore(computefun,picklenamefun,override=False,subdirstore='Pickle_Files/',**kwargs):
    # try to load data from file with name picklenamefun(**kwargs)
    # if that fails, compute data using computefun(**kwargs)
    # can recompute data using override, if desired
    picklename = picklenamefun(**kwargs) + '.pkl'
    try:
        if override: 
            raise FileNotFoundError() # dirty
        result = pickle.load(open(subdirstore + picklename,'rb'))
    except FileNotFoundError:
        result = computefun(**kwargs)
        pickle.dump(result,open(subdirstore + picklename,'wb'))
    return result

def readFileForKey(filename,key,numberoflines=1): # goes through file, returns linenumber and next [numberoflines] lines (useful for getting number of header lines, etc.)
    with open(filename, 'r') as f:
        for linenum, line in enumerate(f):
            if line.startswith(key):
                lines = [line]
                for i in range(numberoflines-1):
                    lines.append(next(f))
                return lines, linenum

def findKey(line,keystart,keyend):
    indexstart = str.find(line,keystart)
    line = line[indexstart+len(keystart):]
    indexend = str.rfind(line,keyend)
    return line[:indexend]
                
def fortranReplace(s):
    # because Fortran...
    return multipleReplace(s,['d','D'],['e','E'])
    
def isStringNumeric(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
    
def getFilePrefix(filepath):
    return getFileAttr(filepath)[0]
    
def getFileExt(filepath):
    return getFileAttr(filepath)[1]
    
def getFileAttr(filepath):
    base = os.path.basename(filepath)
    return os.path.splitext(base)

def findCharIndices(line,char):
    return [i for i, letter in enumerate(line) if letter == char]
    
def getMatlabObject(filename):
    return spio.loadmat(filename)
    
def copyReplaceFile(stroldlist,strnewlist,filenameold,filenamenew=None,subdir=''):
    # copy file, replacing strold with strnew in both file and filename
    # if filenew is None, rename new file to filenew
    if filenamenew is None:
        filenamenew = multipleReplace(filenameold,stroldlist,strnewlist)
    if filenamenew != filenameold: # don't overwrite file!
        with open(subdir + filenameold, 'r') as f:
            with open(subdir + filenamenew, 'w') as f2:
                for line in f:
                    f2.write(multipleReplace(line,stroldlist,strnewlist))
                
def multipleReplace(line,stroldlist,strnewlist):
    for strold, strnew in zip(stroldlist,strnewlist):
        line = line.replace(strold,strnew)
    return line
    
def bashBatchRun(bashfile,jobfile,filepref,searchstring='FILE_PREF='):
    with open(jobfile, 'r') as jobs:
        for line in jobs:
            job = line.rstrip()
            with open(bashfile, 'r') as bash:
                with open('temp', 'w') as temp:
                    for line in bash:
                        if line.startswith(searchstring):
                            linenew = searchstring + filepref + job
                            temp.write(linenew + '\n')
                        else:
                            temp.write(line)
            os.remove(bashfile)
            os.rename('temp',bashfile)
            subprocess.call('bash ' + bashfile,shell=True)
    