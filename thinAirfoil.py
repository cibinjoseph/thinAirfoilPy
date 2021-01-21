#!/usr/bin/python3

import numpy as np
from scipy.integrate import simps
import sys

_nPoints = 5001
_eps = np.finfo(np.pi).eps

def checkAirfoilType(airfoilCode):
    if len(airfoilCode) == 4:
        return 'NACA4'
    elif len(airfoilCode) == 5:
        return 'NACA5'
    else:
        raise ValueError('Invalid airfoil type')

def getNACA4(mpxx, n):
    """ Returns NACA 4-digit airfoil camberline coordinates and slope """
    m = int(mpxx[0])/100.0
    p = int(mpxx[1])/10.0

    thetavec = np.linspace(0.0, np.pi, n)
    xvec = (1-np.cos(thetavec))/2.0

    # Split xvec using position of max camber
    xFore = xvec[xvec < p]
    yFore = (2.0*p*xFore - xFore**2.0)/(p*p)
    dyFore = (2.0*p - 2.0*xFore)/(p*p)

    xRear = xvec[xvec >= p]
    yRear = (1.0 - 2.0*p + 2.0*p*xRear - xRear**2.0)/((1 - p)**2.0)
    dyRear = (2.0*p - 2.0*xRear)/((1 - p)**2.0)

    yvec = np.concatenate((yFore, yRear))*m
    dyvec = np.concatenate((dyFore, dyRear))*m

    return thetavec, xvec, yvec, dyvec

def getNACA5(LPQXX, n):
    """ Returns NACA 5-digit airfoil camberline coordinates and slope """
    l = int(LPQXX[0])/30
    p = int(LPQXX[1])/20
    q = int(LPQXX[2])

    isstandard = False
    if q == 0:
        isstandard = True

    thetavec = np.linspace(0.0, np.pi, n)
    xvec = (1-np.cos(thetavec))/2.0

    if isstandard:
        r = 3.33333333333212*(p**3) + \
                0.700000000000909*(p**2) + \
                1.19666666666638*p - \
                0.00399999999996247;
        k1 = 1514933.33335235*(p**4) - \
                1087744.00001147*(p**3) + \
                286455.266669048*(p**2) - \
                32968.4700001967*p + \
                1420.18500000524;
        k2_k1 = 0

        # Split xvec using position of max camber
        xFore = xvec[xvec < r]
        xRear = xvec[xvec >= r]

        yFore = k1/6*(xFore**3-3*r*xFore**2+r**2*(3-r)*xFore)
        dyFore = k1/6*(3*xFore**2-6*r*xFore+r**2*(3-r))
        yRear = k1*r**3*(1-xRear)/6
        dyRear = -k1*r**3/6 + 0*xRear

    else:  # reflexed airfoil
        r = 10.6666666666861*(p**3) - \
                2.00000000001601*(p**2) + \
                1.73333333333684*p - \
                0.0340000000002413
        k1 = -27973.3333333385*(p**3) + \
                17972.8000000027*(p**2) - \
                3888.40666666711*p + \
                289.076000000022
        k2_k1 = 85.5279999999984*(p**3) - \
                34.9828000000004*(p**2) + \
                4.80324000000028*p - \
                0.21526000000003

        # Split xvec using position of max camber
        xFore = xvec[xvec < r]
        xRear = xvec[xvec >= r]

        yFore = k1/6*((xFore-r)**3-k2_k1*(1-r)**3*xFore+(1-xFore)*r**3)
        dyFore = k1/6*(3*(xFore-r)**2-k2_k1*(1-r)**3-r**3)
        yRear = k1/6*(k2_k1*(xRear-r)**3-k2_k1*(1-r)**3*xRear+r**3*(1-xRear))
        dyRear = k1/6*(3*k2_k1*(xRear-r)**2 - k2_k1*(1-r)**3 -r**3)

    yvec = np.concatenate((yFore, yRear))
    dyvec = np.concatenate((dyFore, dyRear))

    return thetavec, xvec, yvec, dyvec

def calcFourierCoeffs(alphaDeg, thetavec, dyvec):
    """ Calculates A0, A1, A2 coeffs present in thin airfoil theory """
    oneBypi = 1.0/np.pi
    A0 = alphaDeg*(np.pi/180.0) - (simps(dyvec, thetavec))*oneBypi
    A1 = (2.0*oneBypi)*simps(dyvec*np.cos(thetavec), thetavec)
    A2 = (2.0*oneBypi)*simps(dyvec*np.cos(2.0*thetavec), thetavec)

    return A0, A1, A2

def getCoeffs(airfoilCode, alphaDeg):
    """ Calculates CL and CM at alphaDeg for NACA 4,5-digit airfoil """
    if checkAirfoilType(airfoilCode) == 'NACA4':
        th, x, y, dy = getNACA4(airfoilCode, _nPoints)
    elif checkAirfoilType(airfoilCode) == 'NACA5':
        th, x, y, dy = getNACA5(airfoilCode, _nPoints)
    A0, A1, A2 = calcFourierCoeffs(alphaDeg, th, dy)

    CL = 2.0*np.pi*(A0+A1*0.5)
    if np.abs(CL) < _eps:
        xcp = 0.25
    else:
        xcp = 0.25-(np.pi*0.25*(A1-A2))/CL

    coeffs = {
        'CL': CL,
        'CM_le': -0.5*(A0+A1-A2*0.5),
        'CM_c/4': -np.pi*0.25*(A1-A2),
        'xcp': xcp
    }

    return coeffs

def usage():
    print('Usage: thinAirfoil.py MPXX [alphaDeg]')
    print('Computes airfoil coefficients using thin airfoil theory')

def main():
    args = sys.argv
    if len(args) == 1:
        usage()

    elif len(args) == 2:
        if args[1] == '-h' or args[1] == '--help':
            usage()
        else:
            airfoilCode = args[1]
            # Print lift curve characteristics for airfoil
            CL0 = getCoeffs(airfoilCode, 0.0)['CL']
            CM  = getCoeffs(airfoilCode, 0.0)['CM_c/4']
            CL2 = getCoeffs(airfoilCode, 2.0)['CL']
            CLa = (CL2-CL0)/(2.0*np.pi/180.0)
            alpha0 = -CL0/CLa*(180.0/np.pi)

            print('CLa  (1/rad) = ' + '{: f}'.format(CLa))
            print('CL0          = ' + '{: f}'.format(CL0))
            print('CM c/4       = ' + '{: f}'.format(CM))
            print('alf0 (deg)   = ' + '{: f}'.format(alpha0))

    else:
        airfoilCode = args[1]
        print(' Alpha       CL         CM_le      CM_c/4     xcp/c' )
        for alfStr in args[2:]:
            alphaDeg = float(alfStr)
            coeffs = getCoeffs(airfoilCode, alphaDeg)
            print('{: f}'.format(alphaDeg), end='  ')
            print('{: f}'.format(coeffs['CL']), end='  ')
            print('{: f}'.format(coeffs['CM_le']), end='  ')
            print('{: f}'.format(coeffs['CM_c/4']), end='  ')
            print('{: f}'.format(coeffs['xcp']))


if __name__ == '__main__':
    main()
