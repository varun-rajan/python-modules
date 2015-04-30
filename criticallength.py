import numpy as np
import scipy.optimize as spo
import myio as Mio
import mydictionaries as Mdict
import lineplot as lp
    
def getOkabeData(fibertype,matrixtype,Vf,sigmayoption,matfilename,normoption=False):
    mat = Mio.getAndStore(getMatlabObject,Mio.getFilePrefix,override=True,filename=matfilename)
    key = fibertype + matrixtype + str(round(Vf*100)).rjust(2, '0')
    data = mat[key]
    cp = properties(fibertype,matrixtype,Vf)
    sigmay = getSigmay(cp,0,sigmayoption) # initial sigmay (zero load)
    sigmac, deltac = getGLSScalings(cp,sigmay)
    return normalizeData(data,sigmac,deltac,Vf,normoption)

def getMatlabObject(filename,subdir = 'Mat_Files/'):
    return Mio.getMatlabObject(subdir + filename)

def properties(fibertype,matrixtype,Vf):
    # fiber properties: radius, modulus
    fiberdict =	{'GF': {'rf': 6.5e-3, 'Ef': 76e9, 'L0': 24.0, 'rho': 6.34, 's0': 1550.0e6},
                'CF': {'rf': 3.5e-3, 'Ef': 230e9, 'L0': 25.0, 'rho': 5.55, 's0': 4316.0e6},
                'GF2': {'rf': 10.0e-3, 'Ef': 76e9},
                'Fake1': {'rf': 1., 'Ef': 1., 'L0': 1., 'rho': 4, 's0': 1.},
                'Fake2': {'rf': 1., 'Ef': 1., 'L0': 1., 'rho': 10, 's0': 1.},
                }
    # matrix properties: sigmayi, sigmaylong, Em
    matrixdict =	{'Epoxy': {'syi': 73e6, 'syl': 75e6, 'Em': 3.4e9},
                    'PP': {'syi': 11e6, 'syl': 30e6, 'Em': 1.8e9},
                    'PP2': {'syi': 26e6, 'syl': 26e6, 'Em': 1.8e9},
                    'PP3': {'syi': 25e6, 'syl': 25e6, 'Em': 1.8e9},
                    'Fake': {'syi': 0.0001, 'syl': 0.0001, 'Em': 1.8e2}
                    }
    # gls strength: sigmalong, Vflong
    sigmalongdict =	{'GFEpoxy': {'slong': 971e6, 'Vflong': 0.542, 'eta': 1.0},
                    'GFPP': {'slong': 800e6, 'Vflong': 0.542, 'eta': 1.0},
                    'GF2PP2': {'slong': 0.0e6, 'Vflong': 0.542, 'L': 3.1, 'eta': 1.0},
                    'CFPP3': {'slong': 0.0e6, 'Vflong': 0.20, 'eta': 0.375},
                    'Fake1Fake': {'slong': 0.0e6, 'Vflong': 0.20, 'eta': 1.0},
                    'Fake2Fake': {'slong': 0.0e6, 'Vflong': 0.20, 'eta': 1.0}
                    }
    fp = fiberdict[fibertype]
    mp = matrixdict[matrixtype]
    sp = sigmalongdict[fibertype + matrixtype]
    cpold = {'Vf': Vf}
    return Mdict.dictUnion([cpold,fp,mp,sp])
    
def model(fibertype,matrixtype,Vf,fun,sigmayoption,lengthoption,normoption=False):
    cp = properties(fibertype,matrixtype,Vf)
    sigmay = getSigmay(cp,0,sigmayoption) # initial sigmay (zero load)
    sigmac, deltac = getGLSScalings(cp,sigmay)
    data = modelData(cp,fun,sigmay,sigmac,deltac,lengthoption)
    return normalizeData(data,sigmac,deltac,Vf,normoption)
    
def modelData(cp,fun,sigmay,sigmac,deltac,lengthoption,Lmin=0.01,Lmax=7.0,npoints=200): # in mm
    Lvec = np.linspace(Lmin,Lmax,npoints)*deltac
    sigmaappvec = np.fromiter([fun(cp,L,sigmay,sigmac,deltac,lengthoption) for L in Lvec],'float',count=-1)
    return np.column_stack((Lvec,sigmaappvec))
        
def shortFiberGLS(cp,Lcurr,sigmay,sigmac,deltac,lengthoption):
    Vf, rho, eta = cp['Vf'], cp['rho'], cp['eta']
    if lengthoption == 'Exp':
        sigma = lambda Tbar: -sigmaSubSimp(deltac,sigmac,Lcurr,rho,Tbar)
        guessandbounds = [0.001,1.,5.]
        Tbarmin = spo.minimize_scalar(sigma,guessandbounds)
        sigmacr = -Tbarmin.fun*sigmac
    elif lengthoption == 'Unif':
        sigma = lambda Tbar: -sigmaSub(deltac,sigmac,Lcurr,rho,Tbar)
        eps = 0.2
        if Lcurr < eps*deltac:
            sigmacr = (Lcurr/deltac)*sigmac/2
        else:
            guess = 0.99*(Lcurr/deltac)
            guessandbounds = [0,guess,max(5.,guess*2)]
            rhobarmin = spo.minimize_scalar(sigma,guessandbounds)
            sigmacr = -rhobarmin.fun*sigmac
    elif lengthoption == 'Approx':
        sigmacr = ((Lcurr/deltac) - (rho-2)/2/(rho-1)*(Lcurr/deltac)**2)*sigmac
    return eta*(Vf*sigmacr + (1-Vf)*sigmay)

def sigmaSubSimp(deltac,sigmac,Lcurr,rho,Tbar):
    rhobar = lambda Tbar: Tbar**rho + deltac/Lcurr
    sigma = lambda Tbar: (1 - np.exp(-rhobar(Tbar)*Tbar))/rhobar(Tbar)
    eps = 1.e-5
    if Tbar < eps:
        return Tbar*(1 - 0.5*Tbar*deltac/Lcurr)
    else:
        return sigma(Tbar)
    
def sigmaSub(deltac,sigmac,Lcurr,rho,Tbar):
    rhobar = lambda Tbar: Tbar**rho
    interm = lambda Tbar: np.exp(-Lcurr*rhobar(Tbar)/deltac)
    eps = 1.e-9
    if Tbar < eps:
        return 0
    elif Tbar > 1.01*Lcurr/deltac:
        return 0
    else:
        return -(-1 + np.exp(-rhobar(Tbar)*Tbar) + interm(Tbar)*rhobar(Tbar)*Tbar)/(rhobar(Tbar)*(1-interm(Tbar)))
        
def sigmaSub2(deltac,sigmac,Lcurr,rho,Tbar):
    rhobar = lambda Tbar: Tbar**rho
    interm = lambda Tbar: np.exp(-Lcurr*rhobar(Tbar)/deltac)
    eps = 1.e-5
    if rhobar(Tbar) < eps:
        return Tbar*(1 - 0.5*Tbar*deltac/Lcurr)
    else:
        return -(-1 + np.exp(-rhobar(Tbar)*Tbar) + interm(Tbar)*rhobar(Tbar)*Tbar)/(rhobar(Tbar)*(1-interm(Tbar)))
        
def matrixSlip(cp,Lcurr,sigmay,sigmac,deltac,lengthoption): # get 
    rf, Vf, eta = cp['rf'], cp['Vf'], cp['eta']
    A = getA(lengthoption)
    fac = A*(1-0.5*A)
    return eta*(fac*Vf*sigmac*(Lcurr/deltac) + sigmay*(1-Vf))
    
def matrixSlip2(cp,Lcurr,sigmay,sigmac,deltac,lengthoption):
    rho, rf, Vf, eta = cp['rho'], cp['rf'], cp['Vf'], cp['eta']
    A = getA(lengthoption)
    if lengthoption == 'Exp':
        f = sigmaSubSimp
    elif lengthoption == 'Unif':
        f = sigmaSub2
    Tbar = A*Lcurr/deltac
    sbar = f(deltac,sigmac,Lcurr,rho,Tbar)
    return eta*(Vf*sigmac*sbar + sigmay*(1-Vf))

def getSigmay(cp,sigmaapp,sigmayoption):
    sigmayi, sigmaylong = cp['syi'], cp['syl']
    if sigmayoption == 1:
       return sigmayi;
    elif sigmayoption == 2:
       return sigmayi + sigmaapp*(sigmaylong-sigmayi)/sigmaylong;
    elif sigmayoption == 3:
       return sigmaylong
    elif sigmayoption == 4:
        return 0

def getA(lengthoption):
    if lengthoption == 'Exp':
        return 1/2
    elif lengthoption == 'Unif':
        return 1/4
        
def getGLSScalings(cp,sigmay):
    rf, Ef, L0, sigma0, rho, Vf, eta = cp['rf'], cp['Ef'], cp['L0'], cp['s0'], cp['rho'], cp['Vf'], cp['eta']
    tau = sigmay/np.sqrt(3)
    sigmac = (tau*L0/rf)**(1/(rho+1))*sigma0**(rho/(rho+1))
    deltac = sigmac*rf/tau
    return sigmac, deltac
    
def normalizeData(data,sigmac,deltac,Vf,normoption):
    if normoption:
        return data/[deltac,sigmac*Vf]
    else:
        return data