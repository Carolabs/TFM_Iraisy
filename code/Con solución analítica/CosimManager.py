import numpy as np


# Residual power
# tEnd          Final time
# H:        	Macro step-size
class ResidualPower:

    # _______________________________________________________ CONSTRUCTOR
    def __init__(self, tEnd, H):
        # Properties
        self.__H                    = H
        self.__tEnd                 = tEnd
        self.__power                = np.zeros(round(self.__tEnd / self.__H) + 1)   # Vector de residuos de potencia
        self.__inputCounter         = {}                                            # Diccionario que contiene el contador que apunta a la posición del vector de inputs donde se graba el siguiente input
        self.__outputCounter        = {}                                            # Diccionario que contiene el contador que apunta a la posición del vector de outputs donde se graba el siguiente output
        self.__inputs               = {}                                            # Diccionario que contiene el vector con los inputs
        self.__outputs              = {}                                            # Diccionario que contiene el vector con los outputs
        self.__names                = []                                            # Vector con los nombres

    def addSubsystem(self, name):

        # Al añadir un subsistema primero crea la clave 'name' e inicializa a 0 los contadores
        self.__inputCounter[name] = 0
        self.__outputCounter[name] = 0

        # Crea los vectores de inputs y outputs y los rellena de 0
        self.__inputs[name] = np.zeros(round(self.__tEnd / self.__H) + 1)
        self.__outputs[name] = np.zeros(round(self.__tEnd / self.__H) + 1)
        
        # Añade el nombre al vector de nombres
        self.__names.append(name)

    def updateInputCounter(self, name):
        # Actualiza la posición donde se graba el siguiente input
        self.__inputCounter[name] += 1

    def updateOutputCounter(self, name):
        # Actualiza la posición donde se graba el siguiente output
        self.__outputCounter[name] += 1

    def assignInput(self, name, input):
        # Añade el input al vector de inputs para este subsistema
        self.__inputs[name][self.__inputCounter[name]] = input

    def assignOutput(self, name, output):
        # Añade el output al vector de outputs para este subsistema
        self.__outputs[name][self.__outputCounter[name]] = output
    
    def evaluatePower(self):
        # Comprueba que los contadores apuntan a la misma posición
        if self.__inputCounter['SS1'] == self.__outputCounter['SS1'] and self.__inputCounter['SS2'] == self.__outputCounter['SS2']:
            # Evalúa el residuo de potencia para la última posición dP = -u1y1 - u2y2
            self.__power[self.__inputCounter['SS1'] - 1] = -self.__inputs['SS1'][self.__inputCounter['SS1'] - 1] * self.__outputs['SS1'][self.__outputCounter['SS1'] - 1] - self.__inputs['SS2'][self.__inputCounter['SS2'] - 1] * self.__outputs['SS2'][self.__outputCounter['SS2'] - 1] 
        else:
            raise Exception('Error')

        return self.__power[self.__inputCounter['SS1'] - 1]

    def getPower(self):
        # Devuelve todo el array con los residuos de potencia
        return self.__power


# Create the manager
# SYS:          System properties class
# H:        	Macro step-size
# h1:        	Macro step-size SS1
# h2:        	Macro step-size SS2
# tEnd          Final time
class JacobiManagerSR:

    # _______________________________________________________ CONSTRUCTOR
    def __init__(self, SYS, H, h1, h2, tEnd):
        # Properties
        self.__H          = H           # Paso del maanger
        self.__h1         = h1          # Paso SS1
        self.__h2         = h2          # Paso SS2
        self.__tEnd       = tEnd        # Tiempo final
        self.__storeIdx   = 0           # Índice donde grabar
        self.__time       = 0.0         # Tiempo manager
        self.__time1      = 0.0         # Tiempo SS1
        self.__time2      = 0.0         # Tiempo SS2
        self.__dP         = ResidualPower(tEnd, H)
        self.__dPcorr     = ResidualPower(tEnd, H)
        
        
        # Create structures for storage
        self.__deltaP       = 0.0   # Residual power
        self.__deltaE       = 0.0   # Residual energy
        self.__deltaPcorr   = 0.0   # Residual power
        self.__deltaEcorr   = 0.0   # Residual energy
        self.__fCorr        = 0.0   # Corrective force
        self.__WfCorr       = 0.0   # Work of the corrective force during a step
        self.__Balance      = 0.0   # Energy balance: (deltaE/2.0 - dWfCorr)


        # Mass, damping, and stiffness matrices
        self.__M    = SYS.M
        self.__K    = SYS.K
        self.__C    = SYS.C

        # _________________________________________________________________________
        #                                                                     Store
        self.__t        = np.zeros(round(tEnd / H) + 1)         # Tiempo
        self.__fc       = np.zeros(round(tEnd / H) + 1)         # Fuerza de acoplamiento
        self.__x1       = np.zeros(round(tEnd / H) + 1)         # Posición 1
        self.__x1d      = np.zeros(round(tEnd / H) + 1)         # Velocidad 1
        self.__x1dd     = np.zeros(round(tEnd / H) + 1)         # Aceleración 1
        self.__x2       = np.zeros(round(tEnd / H) + 1)
        self.__x2d      = np.zeros(round(tEnd / H) + 1)
        self.__x2dd     = np.zeros(round(tEnd / H) + 1)
        self.__T        = np.zeros(round(tEnd / H) + 1)         # Energía cinética
        self.__V        = np.zeros(round(tEnd / H) + 1)         # Energía potencial
        self.__E        = np.zeros(round(tEnd / H) + 1)         # Energía total
        self.__dWfCorr  = np.zeros(round(tEnd / H) + 1)         # Trabajo de la fuerza de corrección
        self.__fcorr    = np.zeros(round(tEnd / H) + 1)         # Fuerza de corrección
        self.__dE       = np.zeros(round(tEnd / H) + 1)         # Residuo de energía
        self.__dEcorr   = np.zeros(round(tEnd / H) + 1)
        self.__dEAcc    = np.zeros(round(tEnd / H) + 1)         # Acumulado del residuo de energía
        self.__dEcorrAcc= np.zeros(round(tEnd / H) + 1)
        self.__Ebalance = np.zeros(round(tEnd / H) + 1)         # Balance de energía
        self.__u        = np.zeros((round(tEnd / H) + 1, 3))    # Inputs


       

        # _________________________________________________________________________
        #                                                       Connectivity matrix
        self.__L    = np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]])
        self.__yEx  = np.zeros(3)       # Outputs
        self.__uEx  = np.zeros(3)       # Inputs
        self.__y1   = np.zeros(4)       # Outputs SS1
        self.__y2   = np.zeros(3)       # Outputs SS2


    # _______________________________________________________ ASSIGN SUBSYSTEMS
    def assignSS1(self, SS1):
        # Asigna el subsistema como SS1 con su nombre y paso de tiempo
        self.__SS1 = SS1
        self.__SS1.setStepSize(self.__h1)
        self.__dP.addSubsystem('SS1')
        self.__dPcorr.addSubsystem('SS1')
        #print('SS1 assigned correctly')
    
    def assignSS2(self, SS2):
        self.__SS2 = SS2
        self.__SS2.setStepSize(self.__h2)
        self.__dP.addSubsystem('SS2')
        self.__dPcorr.addSubsystem('SS2')
        #print('SS2 assigned correctly')

    # _______________________________________________________ POWER RESIDUAL

    def assignParams(self, mu, k, limit):
        self.__mu       = mu            # Proporcional
        #self.__mu       = 0.2 + 0.6 * mu
        self.__k        = k             # Integral
        #self.__k        = 0.05*(k -0.5)
        self.__limit    = limit         # Límite de corrección
        
    # _______________________________________________________ EVALUATION
    def evalManager(self):
        # y1: force, x1, x1, x1dd
        # y2: x2, x2d, x2dd

        # _________________________________________________________________________
        #                                               Update inputs to subsystems

        # Outputs actually exchanged between subsystems
        # These do not include any correction
        self.__yEx = np.array((self.__y1[0], self.__y2[0], self.__y2[1]))
        self.__uEx = np.matmul(self.__L, self.__yEx)

    # _______________________________________________________ STORAGE
    def store(self):

        # Guarda las variables más relevantes de la cosimulación en cada paso de tiempo

        # Store system outputs
        self.__t[self.__storeIdx]        = self.__time
        self.__fc[self.__storeIdx]       = -self.__y1[0]
        self.__x1[self.__storeIdx]       = self.__y1[1]
        self.__x1d[self.__storeIdx]      = self.__y1[2]
        self.__x1dd[self.__storeIdx]     = self.__y1[3]
        self.__x2[self.__storeIdx]       = self.__y2[0]
        self.__x2d[self.__storeIdx]      = self.__y2[1]
        self.__x2dd[self.__storeIdx]     = self.__y2[2]
        self.__u[self.__storeIdx, :]     = self.__uEx

        self.__dWfCorr[self.__storeIdx]  = self.__WfCorr
        #self.__dE[self.__storeIdx]       = self.__deltaE
        #self.__dEcorr[self.__storeIdx]   = self.__deltaEcorr
        #self.__Ebalance[self.__storeIdx] = self.__Balance
        self.__fcorr[self.__storeIdx]    = self.__fCorr


        # Evaluate magnitudes
        # Assemble arrays
        x  = np.array((self.__y1[1], self.__y2[0]))
        xd = np.array((self.__y1[2], self.__y2[1]))

        # Evaluation of kinetic and potential energy
        self.__T[self.__storeIdx] = 0.5 * np.matmul(xd.T, np.matmul(self.__M, xd))   # T = 0.5*vt*M*v  energía cinética
        self.__V[self.__storeIdx] = 0.5 * np.matmul(x.T, np.matmul(self.__K, x))     # V = 0.5*xt*k*x  energía potencial
        self.__E[self.__storeIdx] = self.__T[self.__storeIdx] + self.__V[self.__storeIdx]
        
        # Update index
        self.__storeIdx += 1
    
    # _______________________________________________________ INITIALIZATION
    def initialize(self):
        # Iterative initialization
        y1It    = np.zeros(4)
        error   = 10
        tol     = 1e-10

        # Inicializa los subsistemas en un proceso iterativo
        # El proceso finaliza con la convergencia de los valores de los outputs
        # Para el LO no sería necesario
        while error > tol:

            # Initialize units
            self.__SS1.initialize(0.0)
            self.__SS2.initialize(0.0)

            # Retrieve outputs of subsystems
            self.__y1 = self.__SS1.readOutputs()
            self.__y2 = self.__SS2.readOutputs()
        
            # Evaluation
            self.evalManager() # Aplica la matriz de conectividad uEx = L * yEx, donde yEx se obtiene de y1 e y2

            # Send inputs
            self.__SS1.setInputs(self.__uEx[0:2]) # Aquí envía posición y velocidad del SS2 como input del SS1
            self.__SS2.setInputs(self.__uEx[2:]) # Aquí envía fuerza del SS1 como input del SS2

            # Evaluate error
            self.__deltaP = 0
            self.__deltaPcorr = 0
            self.__deltaE = 0
            self.__deltaEcorr = 0
            error = np.linalg.norm(y1It - self.__y1, ord = 2) # Comprueba que la fuerza de acoplamiento entre iteraciones no cambie
            y1It = self.__y1

        # Store initial equilibrium
        self.store()

        return self.__SS1.readOutputs(), self.__SS2.readOutputs()

    # _______________________________________________________ MAIN LOOP
    def run(self):
        # The simulation loop follows a Jacobi non-iterative scheme
        while self.__tEnd - self.__time > self.__H / 2.0:

            # Report when required
            #if self.__storeIdx % 1000 == 0:
                #print('Macro time: ', round(self.__time, 4), 's')
           
            # Assign inputs
            self.__dP.assignInput('SS1', self.__uEx[1]) # Velocidad
            self.__dP.updateInputCounter('SS1')

            self.__dP.assignInput('SS2', -self.__y1[0]) # Fuerza
            self.__dP.updateInputCounter('SS2')

            self.__dPcorr.assignInput('SS1', self.__uEx[1]) # Velocidad
            self.__dPcorr.updateInputCounter('SS1')

            self.__dPcorr.assignInput('SS2', -self.__y1[0] + self.__fCorr) # Fuerza + fcorr
            self.__dPcorr.updateInputCounter('SS2')

            # Exchange coupling variables and call manager
            self.__y1 = self.__SS1.readOutputs()
            self.__y2 = self.__SS2.readOutputs()

            # Assign outputs
            self.__dP.assignOutput('SS1', self.__y1[0])
            self.__dP.updateOutputCounter('SS1')
            self.__dP.assignOutput('SS2', self.__y2[1])
            self.__dP.updateOutputCounter('SS2')

            self.__dPcorr.assignOutput('SS1', self.__y1[0])
            self.__dPcorr.updateOutputCounter('SS1')
            self.__dPcorr.assignOutput('SS2', self.__y2[1])
            self.__dPcorr.updateOutputCounter('SS2')

            # Eval correction terms (copiado del ejemplo de Matlab)
            self.__deltaP = self.__dP.evaluatePower()
            self.__deltaE = self.__deltaP * self.__H

            self.__deltaPcorr = self.__dPcorr.evaluatePower()
            self.__deltaEcorr = self.__deltaPcorr * self.__H


            # Evalúa dE acumulado
            self.__dE[self.__storeIdx] = self.__deltaE
            self.__dEAcc[self.__storeIdx] = self.__dEAcc[self.__storeIdx-1] + self.__deltaE

            self.__dEcorr[self.__storeIdx] = self.__deltaEcorr
            self.__dEcorrAcc[self.__storeIdx] = self.__dEcorrAcc[self.__storeIdx-1] + abs(self.__deltaEcorr)

            # Evaluation
            self.evalManager() # Aplica la matriz de conectividad uEx = L * yEx, donde yEx se obtiene de y1 e y2

            # Real work of the corrective force
            self.__WfCorr = -self.__fCorr * self.__uEx[1] * self.__H # fCorr is the correction force from the previous step

            # Evaluar fuerza de corrección
            #self.__fCorr = (self.__mu * self.__deltaE + self.__k * self.__Balance) / (self.__uEx[1] * self.__H)
            self.__fCorr = (self.__mu * self.__deltaP) / (self.__uEx[1])


            # Limita el valor de la corrección
            if self.__fCorr > abs(self.__limit) * abs(self.__uEx[2:]):
                self.__fCorr = abs(self.__limit) * abs(self.__uEx[2:])

            if self.__fCorr < -abs(self.__limit) * abs(self.__uEx[2:]):
                self.__fCorr = -abs(self.__limit) * abs(self.__uEx[2:])

            # Store system states
            if self.__time > self.__H / 2.0:
                self.store()
    

            # Send inputs
            self.__SS1.setInputs(self.__uEx[0:2]) # Aquí envía posición y velocidad del SS2 como input del SS1
            self.__SS2.setInputs(self.__uEx[2:] + self.__fCorr) # Envía la fuerza obtenida en SS1 al SS2, añadiendo la corrección
            
            # Increase goal time
            self.__time += self.__H
            
            # Advance integrations
            while abs(self.__time1 - self.__time) > self.__h1 / 2.0:
                self.__time1 += self.__h1
                self.__SS1.doStep(self.__time1)

            while abs(self.__time2 - self.__time) > self.__h2 / 2.0:
                self.__time2 += self.__h2
                self.__SS2.doStep(self.__time2)


        
        # Save last step
        # Report when required
        #if self.__storeIdx % 1000 == 0:
            #print('Macro time: ', round(self.__time, 4), 's')
        
        # Exchange coupling variables and call manager
        self.__y1 = self.__SS1.readOutputs()
        self.__y2 = self.__SS2.readOutputs()
        
        # Evaluation
        self.evalManager()

        # Store system states
        if self.__time > self.__H / 2.0 and round(self.__time * 1000) % 10 == 0:
            self.store()

    # _______________________________________________________ TERMINATE
    def terminate(self):
        self.__SS1.terminate()
        self.__SS2.terminate()

    # _______________________________________________________ ACCESSORS
    def getSTORE(self):
        return  self.__t, self.__fc, self.__x1, self.__x1d, self.__x1dd, self.__x2, self.__x2d, self.__x2dd, self.__T, self.__V, self.__E, self.__Ebalance

    def evalResidual(self):
        return self.__dP, self.__dE, self.__dEAcc

    def getError(self):
        return self.__dEcorrAcc[-2]
