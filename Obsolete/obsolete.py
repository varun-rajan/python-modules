def readLastNumericLine(filename,key,delimiter=' '):
    # for looking at Lammps log files
    # reads last numeric line before a specific key
    # e.g. reads converged configuration after several iterations
    numericlines = []
    readflag = False
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if readflag:
                if isLineNumeric(line,delimiter):
                    lineprev = line
                else:
                    readflag = False
                    numericlines.append(lineprev)
            if line.startswith(key):
                readflag = True
    return np.loadtxt(numericlines)

def readBlockFile(filename,keystarts,keyends,delimiter=' ',fortranoption=True,keyfun=findKey,keylast=0):
    mydict, numericlines = {}, []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if isLineNumeric(line,delimiter,fortranoption):
                if fortranoption:
                    line = fortranReplace(line)
                numericlines.append(line)
            else:
                if numericlines:
                    if key is not None:
                        mydict[key] = np.loadtxt(numericlines)
                        keylast = key
                    numericlines = []
                key = keyfun(line,keystarts,keyends,keylast)
    if numericlines: # last block
        if key is not None:
            mydict[key] = np.loadtxt(numericlines)
    return mydict
    
def isLineNumeric(line,delimiter,fortranoption=False):
    # just checks first and last whole string (between delimiters)
    linelist = line.split(delimiter)
    s, t = linelist[0], linelist[-1]
    if fortranoption:
        s, t = fortranReplace(s), fortranReplace(t)
    return (isStringNumeric(s) and isStringNumeric(t))

def findKey(line,keystarts,keyends):
    (index, key) = findKeyLeft(line,keystarts)
    line = line[index+len(key):]
    (index, key) = findKeyRight(line,keyends)
    return line[:index]
    
def findKey2(line,keystarts,keyends):
    (index, key) = findKeyLeft(line,keystarts)
    line = line[index+len(key):]
    (index, key) = findKeyLeft(line,keyends)
    return line[:index]

def findKeyLeft(line,keylims): # finds all keys, reading the line from the left, picks the rightmost one
    return findKeySub(line,keylims,op.ge,str.find,0)

def findKeyRight(line,keylims): # finds all keys, reading the line from the right, picks the leftmost one
    return findKeySub(line,keylims,op.le,str.rfind,len(line))
    
def findKeySub(line,keylims,comparefun,findfun,indexdefault):
    index, key = indexdefault, ''
    for keytry in keylims:
        indexcurr = findfun(line,keycurr)
        if indexcurr != -1 and comparefun(indexcurr,index):
            index, key = indexcurr, keycurr
    return (index, key)