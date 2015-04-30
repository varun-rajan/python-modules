import ddutilities_structureplot as ddusp
import contourplot as cp
import lineplot as lp
import ddutilities_io as dduio
import strains
import numpy as np
import mydictionaries as Mdict

def stressPlot(inputfile,increment,variable,axislims=[],scaletype=[None,1e6],contourspacing=5,datalims=[]):
    (filetype, datakey) = getFileInfo(variable)
    datadict = dduio.readDDFile(inputfile,filetype,increment)
    data = datadict[datakey]
    xvec, yvec, zvec = getData1(variable,data)
    spacing = cp.getSpacingGriddedVecList([xvec,yvec]) # need to fix; broken
    (xmat, ymat, zmat) = cp.interpolateScattered(xvec, yvec, zvec, spacing)
    cp.contourPlot(xmat,ymat,zmat,contourspacing,axislims,datalims,scaletype)

def stressStrainPlot(inputfile,component1,component2):
    datadict = dduio.readDDFile(inputfile,'strstr')
    for data in datadict.values(): # only 1 item
        xy, scale = getData2(component1,component2,data)
        lp.linePlot(data=[xy],scale=scale)
    
def slipPlot(inputfile,increment,phivec,contourspacing=[],axislims=[],datalims=[],scaletype=[None,1.]):
    (filetype, datakey) = getFileInfo('u')
    datadict = dduio.readDDFile(inputfile,filetype,increment)
    data = datadict[datakey]
    (x, y, u) = getData1('u',data)
    (x, y, v) = getData1('v',data)
    spacing = cp.getSpacingGriddedVecList([x,y])
    (xmat, ymat, umat) = cp.interpolateScattered(x, y, u, spacing)
    (xmat, ymat, vmat) = cp.interpolateScattered(x, y, v, spacing)
    strainmat = strains.slipSub(umat,vmat,spacing)
    slipmat = np.zeros(np.shape(strainmat)[0])
    for i, phi in enumerate(phivec):
        e1 = [np.cos(phi), np.sin(phi)]
        e2 = [-np.sin(phi), np.cos(phi)]
        slipmat += np.abs(strains.resolvedStrain(e1,e2,strainmat))
    cp.ContourPlot(xmat,ymat,np.reshape(slipmat,np.shape(xmat)),contourspacing,axislims,datalims,scaletype)

def drawStructure(inputfile,increment,phivec,bounds=[],option=1,phioption=True):
    sndict = dduio.readDDFile(inputfile,'sn')
    dislocdict = dduio.readDDFile(inputfile,'disloc',increment)
    if option == 1:
        objectdict = Mdict.dictUnion([sndict,dislocdict])
        objectlist = ['pos','neg','spos','opos']
    elif option == 2:
        objectdict = sndict
        objectlist = ['spos','opos']
    elif option == 3:
        objectdict = dislocdict
        objectlist = ['pos','neg']
    ddusp.drawAllObjects(phivec,objectdict,objectlist,bounds,phioption)
    
def getFileInfo(variable):
    keys = ['sxx', 'syy', 'sxy', 'u', 'v']
    filetypedict = Mdict.instantiate(keys,['stress']*3 + ['displ']*2) # this depends on what's used in GetDDFileName
    datakeydict = Mdict.instantiate(keys,['Stress']*3 + ['Displacement']*2) # this depends on the key used in .plt file
    return (filetypedict[variable], datakeydict[variable])

def getData1(variable,data):
    keys = ['sxx', 'syy', 'sxy', 'u', 'v']
    indexdict = Mdict.instantiate(keys,[2,3,4,2,3])
    return data[:,0], data[:,1], data[:,indexdict[variable]]
    
def getData2(variable1,variable2,data):
    fudge = 1.1 # Srinath added 10% to wcell??
    keys = ['increment','time','exx','sxx','epamp','sxxT']
    indexdict = Mdict.instantiate(keys,list(range(len(keys))))
    scaledict = Mdict.instantiate(keys,[1,1.e9,1.e2*fudge,1.e6,1.e2*fudge,1.e6])
    xy = data[:,[indexdict[variable1],indexdict[variable2]]]
    scale = [scaledict[variable1],scaledict[variable2]]
    return xy, scale