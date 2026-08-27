''' 
System dynamics
Author: Guanyu Chen
Date: 07/24/2024
It incorperates two versions of system dynamics
1. A simple 2-dof horizontal motion (should + elbow)
2. A general model for the 4-dof arm (3 dof on shoulder and 1 on elbow)
We do not consider the dofs on the wrist since the wrsit motion will not affect the transmission of the interaction force to the elbow and shoulder joints.  
'''
import sympy as sp
import numpy as np

# Creating indexed symbols
q1 = sp.Symbol('q1')
q2 = sp.Symbol('q2')
q3 = sp.Symbol('q3')
q4 = sp.Symbol('q4')
q1d = sp.Symbol('q1d')
q2d = sp.Symbol('q2d')
q3d = sp.Symbol('q3d')
q4d = sp.Symbol('q4d')
q1dd = sp.Symbol('q1dd')
q2dd = sp.Symbol('q2dd')
q3dd = sp.Symbol('q3dd')
q4dd = sp.Symbol('q4dd')
tau1 = sp.Symbol('tau1')
tau2 = sp.Symbol('tau2')
tau3 = sp.Symbol('tau3')
tau4 = sp.Symbol('tau4')
u1 = sp.Symbol('u1')
u2 = sp.Symbol('u2')
u3 = sp.Symbol('u3')
u4 = sp.Symbol('u4')


Mh = [0.9, 1.1] # mass, kg
Lh = [0.25, 0.35] # length, m
Sh = [0.11, 0.15] # distance from center of mass to joint, m
Ih = [0.065, 0.1] # rotary inertia, kg*m^2
bh = [0.08, 0.08] # coefficient of viscosity, kg*m^2/s 

Mg = []
Lg = []
Sg = []
Ig = []
bg = []

def horizontalDynamics():

    eqn = [(Ih[0] + Ih[1] + 2*Mh[1]*Lh[0]*Sh[1]*sp.cos(q2) + Mh[1]*Lh[0]**2)*q1dd + (Ih[1] + Mh[1]*Lh[0]*Sh[1]*sp.cos(q2))*q2dd - Mh[1]*Lh[0]*Sh[1]*(2*q1d + q2d)*q2d*sp.sin(q2) + bh[0]*q1d - tau1,
     (Ih[1] + Mh[1]*Lh[0]*Sh[1]*sp.cos(q2))*q1dd + Ih[1]*q2dd + Mh[1]*Lh[0]*Sh[1]*q1d**2*sp.sin(q2)+bh[1]*q2d - tau2]

    return eqn

def generalDynamics():

    pass


    return

def main():
    print(horizontalDynamics())

if __name__ == "__main__":
    main()