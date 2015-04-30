import numpy as np
import scipy.linalg as spla

degfac = 180/np.pi
    
def computeAngle(vec1,vec2):
	res = np.dot(normalizeVec(vec1),normalizeVec(vec2))
	res = np.clip(res,-1.0,1.0) # roundoff issues
	return np.arccos(res)

def normalizeVecAll(a,axis): # 0 for columns, 1 for rows
	return np.apply_along_axis(normalizeVec,axis,a)

def normalizeVec(vec):
	return vec/spla.Norm(vec)
    
def getDist(vec1, vec2):
    return spla.Norm(vec1 - vec2)
    
def getAreaTriangle(a,b,c):
    ab = b - a
    ac = c - a
    onevec = np.ones((3,1))
    mat = np.column_stack((ab,ac,onevec))
    return np.abs(1/2*spla.det(mat))

def roundtoEven(number):
    return roundtoIntOffset(number,2,0)
    
def roundtoOdd(number):
    return roundtoIntOffset(number,2,1)
    
def roundtoIntOffset(number,integer,offset):
    return integer*round((number + offset)/integer) - offset
    
def getRotMatrix(phi,degoption=False):
    if degoption: # convert to radians
        phi = phi/degfac
    return np.array([[np.cos(phi), np.sin(phi)],[-np.sin(phi), np.cos(phi)]])
    
def rescaleCoords(data,boundsorig,boundsnew):
    scale = 1/(boundsorig[1] - boundsorig[0])
    numer = data*(boundsnew[1] - boundsnew[0]) 
    return scale*(numer + boundsorig[1]*boundsnew[0] - boundsorig[0]*boundsnew[1])
    
def getDiffVec(vec1,vec2):
    return np.gradient(vec2)/np.gradient(vec1)
    
def maxPolyRoot(coeff): # find maximum of polynomial given its coefficients. Requires at least one real root
    roots = np.roots(np.polyder(coeff))
    return np.max([np.polyval(coeff,root) if not np.iscomplex(root) else -np.Inf for root in roots]) # returns -inf for complex roots, so these are ignored
    
def setDiff(set1,set2):
    return [x for x in set1 if x not in set2]
    
def quadEig(K,C,M):
    # solves quadratic eigenvalue problem (K + lambda*C + lambda^2*M)*a = 0 by reducing
    # the problem to a generalized eigenvalue problem A*x = lambdanew*B*x
    dim = np.shape(K)[0]
    Arow1 = np.hstack((-C,-K))
    Arow2 = np.hstack((np.eye(dim),np.zeros((dim,dim))))
    A = np.vstack((Arow1,Arow2))
    Brow1 = np.hstack((M,np.zeros((dim,dim))))
    Brow2 = np.hstack((np.zeros((dim,dim)),np.eye(dim)))
    B = np.vstack((Brow1,Brow2))
    eigvals, eigvecs = spla.eig(A, B)
    eigvecsnew = normalizeVecAll(eigvecs[:dim,:],0) # convert gen. eig. solution back to quad. eig. problem
    return eigvals, eigvecsnew
        