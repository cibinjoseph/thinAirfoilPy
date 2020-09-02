#!/usr/bin/python3

import numpy as np
from scipy.integrate import simps
import sys

_nPoints = 5001
_eps = np.finfo(np.pi).eps

def getNACA4(mpxx, n):
    """ Returns NACA 4-digit airfoil camberline coordinates and slope """
    m = int(mpxx[0])/100.0
    p = int(mpxx[1])/10.0

    thetavec = np.linspace(0.0, np.pi, n)
    xvec = (1-np.cos(thetavec))/2.0
    yvec = []
    dyvec = []
    for x in xvec:
        if (x < p):
            y = (2.0*p*x - x*x)/(p*p)
            dy = (2.0*p - 2.0*x)/(p*p)
        else:
            y = (1.0 - 2.0*p + 2.0*p*x - x*x)/((1 - p)**2.0)
            dy = (2.0*p - 2.0*x)/((1 - p)**2.0)
        yvec.append(y)
        dyvec.append(dy)

    yvec = np.array(yvec)*m
    dyvec = np.array(dyvec)*m

    return thetavec, xvec, yvec, dyvec

def calcFourierCoeffs(alphaDeg, thetavec, dyvec):
    """ Calculates A0, A1, A2 coeffs present in thin airfoil theory """
    oneBypi = 1.0/np.pi
    A0 = alphaDeg*(np.pi/180.0) - (simps(dyvec, thetavec))*oneBypi
    A1 = (2.0*oneBypi)*simps(dyvec*np.cos(thetavec), thetavec)
    A2 = (2.0*oneBypi)*simps(dyvec*np.cos(2.0*thetavec), thetavec)

    return A0, A1, A2

def getCoeffs(mpxx, alphaDeg):
    """ Calculates CL and CM at alphaDeg for NACA 4-digit airfoil """
    th, x, y, dy = getNACA4(mpxx, _nPoints)
    A0, A1, A2 = calcFourierCoeffs(alphaDeg, th, dy)
    CL = 2.0*np.pi*(A0+A1*0.5)
    if np.abs(CL) < _eps:
        xcp = None
    else:
        xcp = (np.pi*0.25*(A1-A2))/CL

    coeffs = {
        'CL': CL,
        'CM_le': -0.5*(A0+A1-A2*0.5),
        'CM_c4': -np.pi*0.25*(A1-A2),
        'xcp': xcp
    }

    return coeffs

def main():
    args = sys.argv
    if len(args) == 1:
        print('Usage: thinAirfoil.py MPXX [alphaDeg]')
        print('Computes airfoil coefficients using thin airfoil theory')

    elif len(args) == 2:
        mpxx = args[1]
        # Print lift curve characteristics for airfoil
        CL0 = getCoeffs(mpxx, 0.0)['CL']
        CL2 = getCoeffs(mpxx, 2.0)['CL']
        CLa = (CL2-CL0)/(2.0*np.pi/180.0)
        alpha0 = -CL0/CLa*(180.0/np.pi)

        print('CL0         = ' + str(CL0))
        print('CLa (1/rad) = ' + str(CLa))
        print('alf (deg)   = ' + str(alpha0))

    else:
        mpxx = args[1]
        for alfStr in args[2:]:
            alphaDeg = float(alfStr)
            coeffs = getCoeffs(mpxx, alphaDeg)
            print('Alpha = ' + str(alphaDeg))
            for i in coeffs:
                print(i + ' = '  + str(coeffs[i]))
            print()


if __name__ == '__main__':
    main()
