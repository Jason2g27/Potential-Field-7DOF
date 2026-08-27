import numpy as np
import sympy as sp
import SystemDynamics as sd
from SystemDynamics import u1, u2, u3, u4, q1, q2, q3, q4, q1d, q2d, q3d, q4d, q1dd, q2dd, q3dd, q4dd, tau1, tau2, tau3, tau4
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import math


# initialize
gq0 = [math.radians(10),math.radians(10)]
gqf = [math.radians(90),math.radians(90)]
dim = len(gq0)
gnum = 100
gt0 = 0
gtf = 0.5
galpha = [0]*3*dim
# galpha = [0,0,0,0,0,0]

generaltau = [tau1, tau2, tau3, tau4]
generalq = [q1, q2, q3, q4]
generalqd = [q1d, q2d, q3d, q4d]
generalqdd = [q1dd, q2dd, q3dd, q4dd]
generalu = [u1, u2, u3, u4]
tau = generaltau[:dim]
q = generalq[:dim]
qd = generalqd[:dim]
qdd = generalqdd[:dim]
u = generalu[:dim]

eqn = sd.horizontalDynamics()
eqnSystem = []
for i in range(0, len(eqn)):
    eqnSystem.append(sp.Eq(eqn[i], 0))

dydt = sp.solve(tuple(eqnSystem), tuple(qdd))
dydt = [dydt.get(dy_var, 0) for dy_var in generalqdd[:dim]]
fxu = qd + dydt + u
dXdtMatrix = sp.Matrix(fxu)
dpdt = -dXdtMatrix.jacobian(sp.Matrix(q + qd + tau)).T

def solveSystemDynamics(t0, tf, num, q0, qf, alpha):
    # define boundary conditions
    specificTau = sp.solve(eqnSystem, tuple(tau))
    specificValues = {q1d: 0, q2d: 0, q3d: 0, q4d: 0, q1dd: 0, q2dd: 0, q3dd: 0, q4dd: 0}
    initialSpecificValues = {q[i]: q0[i] for i in range(dim)}
    finalSpecificValues = {q[i]: qf[i] for i in range(dim)}
    tau0 = {var: val.subs(initialSpecificValues|specificValues) for var, val in specificTau.items()}
    tauf= {var: val.subs(finalSpecificValues|specificValues) for var, val in specificTau.items()}
    X0 = [0] * 3 * dim
    Xf = [0] * 3 * dim
    X0[:dim] = q0
    Xf[:dim] = qf
    X0[2*dim:] = [tau0.get(tau_var, 0) for tau_var in generaltau[:dim]]
    Xf[2*dim:] = [tauf.get(tau_var, 0) for tau_var in generaltau[:dim]]
    u0 = alpha[2*dim:]
    z0 = np.concatenate((X0, alpha, u0))
    sol = solve_ivp(system, [t0, tf], z0, t_eval=np.linspace(t0, tf, num))
    return sol

def f(X,u):
    replacement_dict = {}
    symbols = [q,qd,tau]
    varSub = [item for sublist in symbols for item in sublist]
    for symbol, value in zip(varSub, X):
        replacement_dict[symbol] = value
    specificdydt = [expr.subs(replacement_dict) for expr in dydt]
    dXdt = np.concatenate((X[dim:2*dim], specificdydt, u))
    return dXdt

def dphidt(t, phi, X):
    replacement_dict = {}
    symbols = [q,qd,tau]
    varSub = [item for sublist in symbols for item in sublist]
    for symbol, value in zip(varSub, X):
        replacement_dict[symbol] = value
    specificdpdt = dpdt.subs(replacement_dict)
    return specificdpdt @ phi
    

def system(t, z):
    X = z[:3*dim]
    phi = z[3*dim:6*dim]
    u = z[6*dim:]
    dXdt = f(X,u)
    dpdt = dphidt(t, phi, X)
    dudt = dpdt[2*dim:]
    return np.concatenate((dXdt, dpdt, dudt))
    

def main():
    sol = solveSystemDynamics(gt0,gtf,gnum,gq0,gqf,galpha)
    # sol contains X and phi
    # what we want is pos, vel, and torque
    XPhi = sol.y
    pos = XPhi[:dim,:]
    vel = XPhi[dim:2*dim,:]
    torque = XPhi[2*dim:,:]
    time = sol.t
    plt.figure(figsize=(10, 6))
    plt.plot(time, 180/math.pi*pos[0,:], label=r'$\theta_1$')
    plt.plot(time, 180/math.pi*pos[1,:], label=r'$\theta_2$')
    plt.xlabel('Time')
    plt.ylabel('Joint positions [deg]')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()