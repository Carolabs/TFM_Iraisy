from LO_subsystems import mass_yF_uS, mass_yS_uF
from CosimManager import JacobiManagerSR
import numpy as np
from numpy import linalg as LA
from LO_properties import sys
from LO_analytical import analyticalSolution

# Parámetros
# Masa
m1 = 1.0
m2 = 1.0

# Rigidez
k1 = 10.0
k2 = 1000.0
kc = 100.0

# Amortiguamiento
c1 = 0.0
c2 = 0.0
cc = 0.0

# Condiciones iniciales
s10 = 0
s20 = 0
s1d0 = 100
s2d0 = -100

# Crea un sys con todos los parámetros
properties = sys(m1, m2, k1, k2, kc, c1, c2, cc, s10, s1d0, s20, s2d0)

# Pasos de tiempo
h1 = 1e-3       # Macropasos
h2 = 1e-3
H = 1e-3        # Paso del manager

# Tiempo de simulación
tEnd = 10
steps = 1 + np.round(tEnd/H)

# Solución analítica
t, s1, s1d, s1dd, s2, s2d, s2dd, fc, _, _, energy, _, _, _, _ = analyticalSolution(properties, H, tEnd)


# Esta función se llama para hacer la cosimulación.
# Instancia los subsistemas, los inicializa y cosimula
# Por último, devuelve los resultados

def cosimulate(individual):
    mu, k, limit = individual

    # Instancias de ambos sistemas
    ss1 = mass_yF_uS('SS1', properties, 1)
    ss2 = mass_yS_uF('SS2', properties, 2)

    # Manager single-rate
    manager = JacobiManagerSR(properties, H, h1, h2, tEnd) # Crea el manager
    manager.assignParams(mu, k, limit) # Aquí asignamos los parámetros de la corrección

    # Assigning subsystems
    manager.assignSS1(ss1)
    manager.assignSS2(ss2)

    # Inicialización
    _, _ = manager.initialize()

    # Loop
    manager.run()

    # Resultados
    cosim_t, cosim_fc, cosim_s1, cosim_s1d, cosim_s1dd, cosim_s2, cosim_s2d, cosim_s2dd, cosim_T, cosim_V, cosim_E, cosim_Balance = manager.getSTORE()
    cosim_dP, cosim_dE, cosim_dEAcc = manager.evalResidual()
    indicador = manager.getIndicator()

    # Terminación
    manager.terminate()

    # Usar el residuo de potencia modificado
    return indicador[-2]

# Comprueba la energía para las soluciones finales
def check(individual):
    mu, k, limit = individual

    # Instancias de ambos sistemas
    ss1 = mass_yF_uS('SS1', properties, 1)
    ss2 = mass_yS_uF('SS2', properties, 2)

    # Manager single-rate
    manager = JacobiManagerSR(properties, H, h1, h2, tEnd) # Crea el manager
    manager.assignParams(mu, 0*k, limit) # Aquí asignamos los parámetros de la corrección

    # Assigning subsystems
    manager.assignSS1(ss1)
    manager.assignSS2(ss2)

    # Inicialización
    _, _ = manager.initialize()

    # Loop
    manager.run()

    # Soluciones
    cosim_t, cosim_fc, cosim_s1, cosim_s1d, cosim_s1dd, cosim_s2, cosim_s2d, cosim_s2dd, cosim_T, cosim_V, cosim_E, cosim_Balance = manager.getSTORE()
    cosim_dP, cosim_dE, cosim_dEAcc = manager.evalResidual()
    indicador = manager.getIndicator()

    # Terminación
    manager.terminate()

    # Devuelve residuo de energía, enerǵia del sistema cosimulado
    # energía de la solución analítica y balance energético
    return cosim_t, cosim_fc, cosim_s1, cosim_s1d, cosim_s1dd, cosim_s2, cosim_s2d, cosim_s2dd, cosim_T, cosim_V, cosim_E, cosim_dP, cosim_dE, cosim_dEAcc, indicador, energy, cosim_Balance
