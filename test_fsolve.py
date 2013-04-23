#!/usr/bin/env python

from scipy.optimize import fsolve
import uncertainties as u

Fa0 = u.ufloat(5.0, 0.05)
v0 = u.ufloat(10., 0.1)

V = u.ufloat(66000.0, 100)  # reactor volume L^3
k = u.ufloat(3.0, 0.2)      # rate constant L/mol/h

def func(Ca):
    "Mole balance for a CSTR. Solve this equation for func(Ca)=0"
    
    Fa = v0 * Ca     # exit molar flow of A
    ra = -k * Ca**2  # rate of reaction of A L/mol/h
    return Fa0 - Fa + V * ra

# CA guess that that 90 % is reacted away
CA_guess = 0.1 * Fa0 / v0

wrapped_fsolve = u.wrap(lambda func, x0: fsolve(func, x0)[0])
CA_sol = wrapped_fsolve(func, CA_guess)

print 'The exit concentration is {0} mol/L'.format(CA_sol)
