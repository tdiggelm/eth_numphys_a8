from numpy import *
from matplotlib.pylab import *

########################################################
#                                                      #
# TODO (optional): Sie koennen den Runge-Kutta Code    #
#                  aus Serie 5 auf explizite Verfahren #
#                  anpassen und hier benutzen.         #
#                                                      #
# from rk import *                                     #
#                                                      #
########################################################

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


# Butcher Schema
A = array([[0, 0],
           [1, 0]])
b = array([0.5, 0.5])
c = array([0, 1])

#####################################################
#                                                   #
# TODO: Implementieren Sie hier das Heun Verfahren. #
#                                                   #
#####################################################


# Butcher Schema
A = array([[ 0.0, 0.0, 0.0],
           [ 0.5, 0.0, 0.0],
           [-1.0, 2.0, 0.0]])
b = array([1.0, 4.0, 1.0]) / 6.0
c = array([0.0, 0.5, 1.0])

#####################################################################
#                                                                   #
# TODO: Implementieren Sie hier das gegebene Runge-Kutta Verfahren. #
#                                                                   #
#####################################################################


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
    pass

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
    pass

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

grid(True)
legend(loc='best')
xlabel(r'Anzahl Schritte $N$')
ylabel('Fehler')
savefig('ordnungistnichtalles.png')


# Konvergenzraten
rate_ee_a = 0.0
rate_he_a = 0.0
rate_rk_a = 0.0

##########################################################
#                                                        #
# TODO: Berechnen Sie hier die Konvergenzraten zu T=0.09 #
#                                                        #
##########################################################

print("Konvergenzrate Euler T=0.09 : %.4f" % rate_ee_a)
print("Konvergenzrate Heun  T=0.09 : %.4f" % rate_he_a)
print("Konvergenzrate RK    T=0.09 : %.4f" % rate_rk_a)

rate_ee_b = 0.0
rate_he_b = 0.0
rate_rk_b = 0.0

#########################################################
#                                                       #
# TODO: Berechnen Sie hier die Konvergenzraten zu T=0.1 #
#                                                       #
#########################################################

print("Konvergenzrate Euler T=0.1 : %.4f" % rate_ee_b)
print("Konvergenzrate Heun  T=0.1 : %.4f" % rate_he_b)
print("Konvergenzrate RK    T=0.1 : %.4f" % rate_rk_b)
