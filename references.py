import myio as Mio
import mydictionaries as Mdict

allreferences = ['batteries','compositemechanics','computationalfracture','constitutivemodels','dynamicfracture','hierarchicalbiocomposites','hybridcomposites','metallicglasses','misc','multiscale','notches','quasistaticfracture']

def readRefFiles(suffixes=allreferences):
    dictvec = [readRefFile('references_' + suffix + '.bib') for suffix in suffixes]
    return Mdict.dictUnion(dictvec)
    
def readRefFile(filename):
    referencesdict, blockoflines = {}, []
    with open(filename,'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('@'): # start reading, get key, clear previous results
                articletype, key = getReferenceKey(line)
                blockoflines = []
            elif line.startswith('}'): # stop reading, store results
                referencesdict[key] = getReferenceInfo(blockoflines)
                referencesdict[key]['articletype'] = articletype
            else:
                blockoflines.append(line)
    return referencesdict
                
def getReferenceKey(line):
    articletype = Mio.findKey(line,'@','{')
    key = Mio.findKey(line,'{',',')
    return articletype, key
    
def getReferenceInfo(blockoflines):
    res = {}
    for line in blockoflines:
        key, value = getLineInfo(line)
        res[key] = value
    return res
    
def getLineInfo(line):
    key = Mio.findKey(line,'','=').rstrip()
    value = Mio.findKey(line,'{','}')
    return key, value

def searchAuthor(dict,authorname):
    search(dict,'author',authorname)
            
def searchTitle(dict,keyword):
    search(dict,'title',keyword)
    
def search(dict,entrykey,entrykeyval):
    counter = 0
    entrykeyvalnew = entrykeyval.lower()
    for key, val in dict.items():
        try: # some entries don't have key
            entry = val[entrykey]
            if entrykeyvalnew in entry.lower(): # remove case sensitivity
                counter = counter + 1
                try:
                    title = val['title']
                except KeyError:
                    title = '[No Title]'
                print('{0}) Key: {1}; Title: {2}'.format(counter,key,title))
        except KeyError:
            pass            
    
def printReferences(dict,filename):
    with open(filename,'w') as f:
        for key, val in dict.items():
            articletype = val['articletype']
            f.write('@{0}{{{1},\n'.format(articletype,key))
            for keysub, valsub in val.items():
                if keysub != 'articletype':
                    f.write('{0} = {{{1}}},\n'.format(keysub,valsub))
            f.write('}\n\n')
    