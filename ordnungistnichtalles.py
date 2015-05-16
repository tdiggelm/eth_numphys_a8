#!/usr/bin/env python
#############################################################################
# course:   Numerische Methoden D-PHYS
# exercise: assignment 8
# author:   Thomas Diggelmann <thomas.diggelmann@student.ethz.ch>
# date:     16.04.2015
#############################################################################
from numpy import *
from matplotlib.pylab import *
from rk import *

def rhs(t, y):
    """Berechne rechte Seite der Diff.-Gl.

    Input: y ... y(t)
    Output: dydt ...  rechte Seite der Diff.-Gl.
    """
    dydt = abs(1.1 - y) + 1.0
    return dydt


# Exakte Loesungen
yExactT009 = 1.094675696201649
yExactT010 = 1.104700834614225

# Anfangswerte
t0 = 0.0
y0 = 1.0


################################################################
#                                                              #
# TODO: Implementieren Sie hier ein explizites Eulerverfahren. #
#                                                              #
################################################################

def odeint_ee(rhs, y0, t0, T, N, flag=False):
    r"""Integrate ODE with explicit Euler method

    Input: y0     ... initial condition
           t0     ... start x
           T      ... end   x
           N      ... number of steps (h = (xEnd - xStart)/N)
           flag   ... flag == False return complete solution: (phi, phi', t)
                      flag == True  return solution at endtime only: phi(tEnd)

    Output: t ... variable
            y ... solution
    """
    t, h = linspace(t0, T, N, retstep=True)
    y0 = atleast_1d(y0)
    y = zeros((size(y0), N))
    y[:,0] = y0
    for k in xrange(N-1):
        y[:,k+1] = y[:,k] + h * rhs(t[k], y[:,k])
    if flag:
        return t[-1], y[:,-1]
    else:
        return t, y

"""
# implementation of EE with RungeKutta class for testing
def odeint_ee2(rhs, y0, t0, T, N, flag=False):    
    # Butcher Schema
    A = array([[0]])
    b = array([1])
    c = array([0])
    rk = RungeKutta(A, b, c)
    rk.set_rhs(rhs)
    rk.set_iv(t0, atleast_1d(y0))
    t, y = rk.integrate(T, N)
    if flag:
        return t[-1], y[:,-1]
    else:
        return t, y
"""

#####################################################
#                                                   #
# TODO: Implementieren Sie hier das Heun Verfahren. #
#                                                   #
#####################################################


def odeint_he(rhs, y0, t0, T, N, flag=False):    
    # Butcher Schema
    A = array([[0, 0],
               [1, 0]])
    b = array([0.5, 0.5])
    c = array([0, 1])
    rk = RungeKutta(A, b, c)
    rk.set_rhs(rhs)
    rk.set_iv(t0, atleast_1d(y0))
    t, y = rk.integrate(T, N)
    if flag:
        return t[-1], y[:,-1]
    else:
        return t, y

#####################################################################
#                                                                   #
# TODO: Implementieren Sie hier das gegebene Runge-Kutta Verfahren. #
#                                                                   #
#####################################################################

def odeint_rk(rhs, y0, t0, T, N, flag=False):    
    # Butcher Schema
    A = array([[ 0.0, 0.0, 0.0],
               [ 0.5, 0.0, 0.0],
               [-1.0, 2.0, 0.0]])
    b = array([1.0, 4.0, 1.0]) / 6.0
    c = array([0.0, 0.5, 1.0])
    rk = RungeKutta(A, b, c)
    rk.set_rhs(rhs)
    rk.set_iv(t0, atleast_1d(y0))
    t, y = rk.integrate(T, N)
    if flag:
        return t[-1], y[:,-1]
    else:
        return t, y


"""
# code for testing the integrators
def testing():
    figure()
    y0 = -6
    for name, integrator in [("ee2", odeint_ee2), ("ee", odeint_ee), ("he", odeint_he), ("rk", odeint_rk)]:
        t, y = integrator(rhs, y0, -2, 2, 10)
        plot(t, y[0], label=name)
    grid()
    legend(loc='best')
    savefig('test.pdf')
testing()    
from sys import exit
exit()
"""


# Anzahl Schritte fuer Konvergenzstudien
NN = 2**arange(1,10)


# Konvergenzstudie fuer T = 0.09
tEnd   = 0.09

error_ee_a = []
error_he_a = []
error_rk_a = []
for i, N in enumerate(NN):
    ############################################################
    #                                                          #
    # TODO: Berechnen Sie hier die Loesung mit jedem Verfahren #
    #       und anschliessend den Fehler zum Endzeitpunkt fuer #
    #       die verschiedenen Anzahlen an Schritten.           #
    #                                                          #
    ############################################################
    #for name, odeint in integrators:
    t, y = odeint_ee(rhs, y0, t0, tEnd, N, flag=True)
    error_ee_a.append(abs(yExactT009-y[0]))
    
    t, y = odeint_rk(rhs, y0, t0, tEnd, N, flag=True)
    error_rk_a.append(abs(yExactT009-y[0]))
    
    t, y = odeint_he(rhs, y0, t0, tEnd, N, flag=True)
    error_he_a.append(abs(yExactT009-y[0]))

error_ee_a = array(error_ee_a)
error_he_a = array(error_he_a)
error_rk_a = array(error_rk_a)


# Konvergenzstudie fuer T = 0.1
tEnd   = 0.1

error_ee_b = []
error_he_b = []
error_rk_b = []
for i, N in enumerate(NN):
    ############################################################
    #                                                          #
    # TODO: Berechnen Sie hier die Loesung mit jedem Verfahren #
    #       und anschliessend den Fehler zum Endzeitpunkt fuer #
    #       die verschiedenen Anzahlen an Schritten.           #
    #                                                          #
    ############################################################
    t, y = odeint_ee(rhs, y0, t0, tEnd, N, flag=True)
    error_ee_b.append(abs(yExactT010-y[0]))
    
    t, y = odeint_rk(rhs, y0, t0, tEnd, N, flag=True)
    error_rk_b.append(abs(yExactT010-y[0]))
    
    t, y = odeint_he(rhs, y0, t0, tEnd, N, flag=True)
    error_he_b.append(abs(yExactT010-y[0]))

error_rk_b = array(error_rk_b)
error_he_b = array(error_he_b)
error_ee_b = array(error_ee_b)


# Plot Konvergenzstudie
figure()

##########################################################
#                                                        #
# TODO: Plotten Sie den Fehler gegen die Anzahl Schritte #
#       fuer alle Verfahren und beide Endzeiten doppelt  #
#       logarithmisch.                                   #
#                                                        #
##########################################################

# Konvergenzraten

rate_ee_a = -polyfit(log(NN), log(error_ee_a), 1)[0]
rate_he_a = -polyfit(log(NN), log(error_he_a), 1)[0]
rate_rk_a = -polyfit(log(NN), log(error_rk_a), 1)[0]

print("Konvergenzrate Euler T=0.09 : %.4f" % rate_ee_a)
print("Konvergenzrate Heun  T=0.09 : %.4f" % rate_he_a)
print("Konvergenzrate RK    T=0.09 : %.4f" % rate_rk_a)

rate_ee_b = -polyfit(log(NN), log(error_ee_b), 1)[0]
rate_he_b = -polyfit(log(NN), log(error_he_b), 1)[0]
rate_rk_b = -polyfit(log(NN), log(error_rk_b), 1)[0]

print("Konvergenzrate Euler T=0.1 : %.4f" % rate_ee_b)
print("Konvergenzrate Heun  T=0.1 : %.4f" % rate_he_b)
print("Konvergenzrate RK    T=0.1 : %.4f" % rate_rk_b)

loglog(NN, error_rk_a, "b+-", label="RK $T=0.09$ $p=%.4f$" % rate_rk_a)
loglog(NN, error_rk_b, "b+--", label="RK $T=0.10$ $p=%.4f$" % rate_rk_b)
loglog(NN, error_he_a, "r+-", label="HE $T=0.09$ $p=%.4f$" % rate_he_a)
loglog(NN, error_he_b, "r+--", label="HE $T=0.10$ $p=%.4f$" % rate_he_b)
loglog(NN, error_ee_a, "g+-", label="EE $T=0.09$ $p=%.4f$" % rate_ee_a)
loglog(NN, error_ee_b, "g+--", label="EE $T=0.10$ $p=%.4f$" % rate_ee_b)
title("Konvergenzordnungen")
grid(True)
legend(loc='best', prop={'size': 10})
xlabel(r'Anzahl Schritte $N$')
ylabel('Fehler')
savefig('ordnungistnichtalles.pdf')


