import numpy as np
import mdutilities_io as mduio
import myio as Mio
import mymath as Mmath
import Ccircumradii14 as C
import networkx as nx
import scipy.spatial as spsp

def parseCrackData(crackdata,cR,option,timeincrement=0.002):
    time = crackdata[:,0]*timeincrement
    xposition = crackdata[:,3] - crackdata[0,3]
    velocity = Mmath.getDiffVec(time,xposition)
    velocitynorm = velocity/cR
    if option == 1: # position vs. time
        return [time,xposition]
    elif option == 2: # normalized velocity vs. time
        return [time,velocitynorm]
    elif option == 3: # normalized velocity vs. position
        return [xposition,velocitynorm]

def loopCrackNodes(simname,increments,rootdir,bounds,crackoption=2,override1=False,override2=False,**kwargs):
    crackdata = np.empty((np.shape(increments)[0],5))
    for i, increment in enumerate(increments):
        badxy = getAndStoreBadXY(rootdir,crackoption,override=override2,dumparray=None,simname=simname,increment=increment,bounds=bounds,**kwargs) # just load from pickle file
        if badxy is None: # dumparray is not initialized, so need to get from pickle or read from dump file
            try:
                res = mduio.getAndStoreDump(simname,increment,rootdir,bounds=bounds,override=override1,**kwargs)
                badxy = getAndStoreBadXY(rootdir,crackoption,override=True,dumparray=res['array'],simname=simname,increment=increment,bounds=bounds,**kwargs)
            except FileNotFoundError: # file not found (increment wrong), just return what we have till now
                return crackdata[:i,:]
        indexleft, indexright = np.argmin(badxy[:,0]), np.argmax(badxy[:,0])
        crackdata[i,:] = [increment,badxy[indexleft,0],badxy[indexleft,1],badxy[indexright,0],badxy[indexright,1]]
    return crackdata

def getAndStoreBadXY(rootdir,crackoption,override,voidoption=False,pickledir='Pickle_Files/',**kwargs):
    if crackoption == 1: # my Delaunay method
        crackfun = getCrackNodesSub
    elif crackoption == 2: # use centrosymmetry
        crackfun = getCrackNodesSub2
    return Mio.getAndStore(crackfun,getBadXYFilename,override,rootdir+pickledir,voidoption=voidoption,**kwargs)
    
def getCrackNodesSub2(dumparray,centrocutoff=2.0,indexstart=2,tol=2.0,**kwargs):
    if dumparray is None:
        return None
    xyzmatold = dumparray[:,indexstart:indexstart+3]
    centro = dumparray[:,indexstart+3]
    xmin, xmax, ymin, ymax = getBounds(xyzmatold,dim=2) # don't need 3rd dim, since it's periodic and therefore not a free surface
    xyzmat = xyzmatold[centro > centrocutoff,:]
    cond1 = np.abs((xyzmat[:,1] - ymax)) > tol
    cond2 = np.abs((xyzmat[:,1] - ymin)) > tol
    cond3 = np.abs((xyzmat[:,0] - xmax)) > tol
    cond4 = np.abs((xyzmat[:,0] - xmin)) > tol
    return xyzmat[cond1 & cond2 & cond3 & cond4,:]
    
def getBounds(xymat,dim=2):
    xmin, xmax = np.min(xymat[:,0]), np.max(xymat[:,0])
    ymin, ymax = np.min(xymat[:,1]), np.max(xymat[:,1])
    if dim == 3:
        zmin, zmax = np.min(xymat[:,2]), np.max(xymat[:,2])
        return xmin, xmax, ymin, ymax, zmin, zmax
    else:
        return xmin, xmax, ymin, ymax
    
def getBadXYFilename(simname,increment,voidoption,**kwargs):
    return 'bad_xy' + simname + '.' + str(increment) + '.' + str(voidoption)

# sort of obsolete method based on Delaunay triangulation
    
def getCrackNodesSub(dumparray,bounds,voidoption=False,circumcutoff=1.0,indexstart=2,**kwargs):
    if dumparray is None:
        return None
    xymatold = dumparray[:,indexstart:indexstart+2]
    tri = getNearCrackTri(xymatold,bounds)
    xymat, triangles, neighbors = tri.points, tri.simplices, tri.neighbors
    ntri = np.shape(triangles)[0]
    triprops = np.asarray(C.CGetCircumradii(xymat,triangles,ntri))
    badtriangles = triprops[triprops[:,1] > circumcutoff,0].astype('int32')
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
    xmin, xmax, ymin, ymax = getBounds(xymat,dim=2)
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
    