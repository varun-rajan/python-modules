import numpy as np
import mymath as Mym
import matplotlib.pyplot as plt
import mydictionaries as Mdict
import contourplot as cp

def drawAllObjects(phivec,objectdict,objectlist,bounds=[],phioption=True):
    if bounds:
        objectdict = ddusp.enforceBoundsDict(objectdict,bounds)
    plotdict, markerdict = createPlotDicts()
    for object in objectlist:
        plotdictcurr = plotdict[object]
        markerdictcurr = markerdict[object]
        if phioption: # plot for each slip plane separately
            for i, phi in enumerate(phivec):
                mymarker = myMarker(phi,**markerdictcurr)
                key = object + str(i+1)
                drawObjects(objectdict,key,Mdict.writeNewDict(plotdictcurr,i),mymarker)
        else:
            mymarker = MyMarker(0,**markerdictcurr)
            drawObjects(objectdict,object,plotdictcurr,mymarker)
    cp.postProcessPlot()
    
def drawObjects(objectdict,key,plotdict,mymarker):
    try: # possible that object is undefined (empty)
        objectcurr = objectdict[key]
        plt.scatter(objectcurr[:,0],objectcurr[:,1],marker=mymarker,**plotdict)
    except KeyError:
        print('Object ' + key + ' not found')    
    
def enforceBoundsDict(mydict,bounds):
    for key, array in mydict.items():
        mydict[key] = enforceBoundsXY(array,bounds)
    return mydict

def enforceBoundsXY(xy,bounds):
    indexx = bounds[0] <= xy[:,0] <= bounds[1]
    indexy = bounds[2] <= xy[:,1] <= bounds[3]
    return xy[indexx & indexy,:]
    
def myMarker(phi,markertype,rotation):
    if markertype == 'disl':
        return edgeDislMarker(phi+rotation)
    else:
        return markertype
            
def edgeDislMarker(phi):
    # note to self: reorients the marker w.r.t. visual space (i.e. what appears on screen), not real space (i.e. the physical coordinates). The two only coincide if the aspect ratio correct (i.e. plt.axis('equal'))
    scalefac = 1.5
    markerinit = [[-1/scalefac,0],[0,0],[0,1],[0,0],[1/scalefac,0]]
    rotmatrix = Mym.getRotMatrix(phi,degoption=True)
    return np.dot(np.array(markerinit),rotmatrix)
    
def createPlotDicts():
    colorlist = ['r','g','b','m','y']
    markersizedisl = 100 # in points
    markersizesn = 25 # in points
    plotkeys = ['edgecolor','facecolor','s']
    markerkeys = ['markertype', 'rotation']
    objectkeys = ['spos','opos','pos','neg']
    
    nucdictp = Mdict.instantiate(plotkeys,['k','k',markersizesn])
    obsdictp = Mdict.instantiate(plotkeys,['k','None',markersizesn])
    disldictp = Mdict.instantiate(plotkeys,[colorlist,'None',markersizedisl])
    nucdictm = Mdict.instantiate(markerkeys,['o',0])
    obsdictm = Mdict.instantiate(markerkeys,['o',0])
    dislposdictm = Mdict.instantiate(markerkeys,['disl',0])
    dislnegdictm = Mdict.instantiate(markerkeys,['disl',180])
    
    plotdict = Mdict.instantiate(objectkeys,[nucdictp,obsdictp,disldictp,disldictp])
    markerdict = Mdict.instantiate(objectkeys,[nucdictm,obsdictm,dislposdictm,dislnegdictm])
    return (plotdict, markerdict)