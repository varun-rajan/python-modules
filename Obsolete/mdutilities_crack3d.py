import numpy as np
import mdutilities_io as mduio
import myio as Mio
import mymath as Mmath
import Ccircumradii14 as C
import networkx as nx
import scipy.spatial as spsp
import matplotlib as mpl
import matplotlib.tri as mtri
import matplotlib.pyplot as plt

def parseCrackData(crackdata,cR,option,r0=2**(1/6),timeincrement=0.002,cL=1.0):
    time = crackdata[:,0]*timeincrement
    xposition = crackdata[:,3] - crackdata[0,3]
    velocity = Mmath.getDiffVec(time,xposition)
    timenorm = time/(r0/cL)
    velocitynorm = velocity/cR
    xpositionnorm = xposition/r0
    if option == 1: # position vs. time
        return [time,xposition]
    elif option == 2: # normalized velocity vs. time
        return [time,velocitynorm]
    elif option == 3: # normalized velocity vs. position
        return [xposition,velocitynorm]
    elif option == 4:
        return [timenorm,velocitynorm]
    elif option == 5:
        return [xpositionnorm,velocitynorm]

def loopCrackNodes(simname,increments,rootdir,bounds,centrocutoff=10.0,override1=False,override2=False,pickledir='Pickle_Files/',voidoption=False):
    crackdata = np.zeros((np.shape(increments)[0],5))
    for i, increment in enumerate(increments):
        res = mduio.getAndStoreDump(simname,increment,rootdir=rootdir,override=override1,bounds=bounds)
        badxy = getAndStoreBadXY(res['array'],centrocutoff,dimensions,simname,increment,rootdir,override=override2)
        indexleft, indexright = np.argmin(badxy[:,0]), np.argmax(badxy[:,0])
        crackdata[i,:] = [increment,badxy[indexleft,0],badxy[indexleft,1],badxy[indexright,0],badxy[indexright,1]]
    return crackdata

def getAndStoreBadXY(dumparray,circumradius,dimensions,bounds,simname,increment,rootdir,dumpdir='Dump_Files/',pickledir='Pickle_Files/',override=False,voidoption=False,dimensions=2):
    return Mio.getAndStore(getCrackNodesSub,getBadXYFilename,override=override,subdirstore=rootdir+pickledir,simname=simname,increment=increment,dumparray=dumparray,dimensions=dimensions)
        
def getCrackNodesSub(dumparray,centrocutoff,dimensions,indexstart=2,**kwargs):
    xyzmatold = dumparray[:,indexstart:indexstart+dimensions]
    centro = dumparray[:,indexstart+dimensions+1]
    indexcrack = centro > centrocutoff
    
    tri = getNearCrackTri(xymatold,bounds)
    xymat, triangles, neighbors = tri.points, tri.simplices, tri.neighbors
    ntri = np.shape(triangles)[0]
    triprops = np.asarray(C.CGetCircumradii(xymat,triangles,ntri))
    badtriangles = triprops[triprops[:,1] > circumradius,0].astype('int32')
    badtriclusters = findTriClusters(set(badtriangles),neighbors)
    clusterindexvec = findCentralClusters(badtriclusters,triangles,xymat)
    if voidoption: # take all clusters, including possible voids
        badnodeslist = [triangles[badtriclusters[index],:] for index in clusterindexvec]
        badnodes = np.vstack(badnodeslist)
    else: # just take largest connected cluster (which is the first one)
        badnodes = triangles[badtriclusters[clusterindexvec[0]],:]
    return xymat[np.unique(np.ravel(badnodes)),:]

def getNearCrackTri(xymat,bounds):
    indexx = (bounds[0] <= xymat[:,0]) & (xymat[:,0] <= bounds[1])
    indexy = (bounds[2] <= xymat[:,1]) & (xymat[:,1] <= bounds[3])
    return spsp.Delaunay(xymat[indexx & indexy,:])

def findTriClusters(triangles,neighbors):
    G = nx.Graph()
    G.add_nodes_from(triangles)
    for tri in triangles:
        for i in range(3):
            neighborcurr = neighbors[tri,i]
            if neighborcurr in triangles:
                G.add_edge(tri,neighborcurr)
    return nx.connected_components(G)

def findCentralClusters(badtriclusters,triangles,xymat,fac=100,fac2=0.3):
    xmin, xmax, ymin, ymax = getBounds(xymat)
    xlen, ylen = xmax - xmin, ymax - ymin
    res = np.zeros((len(badtriclusters),1))
    for i, cluster in enumerate(badtriclusters):
        goodcount = 0
        for tri in cluster:
            trinodes = triangles[tri,:]
            tricoords = xymat[trinodes,:]
            tricen = np.mean(tricoords,axis=0)
            # tally data on clusters far away from edge (others are spurious)
            if np.abs(tricen[0] - xmin) > xlen/fac:
                if np.abs(tricen[0] - xmax) > xlen/fac:
                    if np.abs(tricen[1] - ymin) > ylen/fac:
                        if np.abs(tricen[1] - ymax) > ylen/fac:
                            goodcount = goodcount + 1
        res[i] = goodcount/len(cluster) > fac2
    return list(np.where(res==True)[0])
    
def getBounds(xymat):
    xmin = np.min(xymat[:,0])
    xmax = np.max(xymat[:,0])
    ymin = np.min(xymat[:,1])
    ymax = np.max(xymat[:,1])
    return xmin, xmax, ymin, ymax
    
def getBadXYFilename(simname,increment,voidoption,**kwargs):
    return 'bad_xy' + simname + '.' + str(increment) + '.' + str(voidoption)
    