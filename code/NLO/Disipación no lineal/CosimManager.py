import numpy as np

# Create an extrapolator for the input extrapolation within the cosim manager
def extrapolate(idx, order, t, extrBuffx, extrBuffy):

    # Determine extrapolation method - we need enough stored points to approximate
    __order = min(idx, order)

    return extrapolationDetail(extrBuffx, extrBuffy, t, __order)

def extrapolationDetail(extrBuffx, extrBuffy, t, order):

    yex = np.zeros(extrBuffy.shape[0])

    if order == 0:
        yex = extrBuffy[:, -1]

    elif order == 1:

        for var in range(0, extrBuffy.shape[0]):
            poly = np.polyfit(extrBuffx[[-2, -1]], extrBuffy[var,[-2, -1]], deg = 1)
            yex[var]  = np.polyval(poly, t)

    elif order == 2:

        for var in range(0, extrBuffy.shape[0]):
            poly = np.polyfit(extrBuffx[[-3, -2, -1]], extrBuffy[var,[-3, -2, -1]], deg = 2)
            yex[var]  = np.polyval(poly, t)

    else:
        raise ValueError('Incorrect extrapolation order')

    return yex

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
        self.__H          = H
        self.__h1         = h1
        self.__h2         = h2
        self.__tEnd       = tEnd
        self.__storeIdx   = 0
        self.__time       = 0.0  # Manager
        self.__time1      = 0.0  # SS1
        self.__time2      = 0.0  # SS2
        self.__corr       = False
        self.__SYS        = SYS



        # _________________________________________________________________________
        #                                                                     Store
        self.__t        = np.zeros(round(tEnd / H) + 1)
        self.__fc       = np.zeros(round(tEnd / H) + 1)
        self.__x1       = np.zeros(round(tEnd / H) + 1)
        self.__x1d      = np.zeros(round(tEnd / H) + 1)
        self.__x1dd     = np.zeros(round(tEnd / H) + 1)
        self.__x2       = np.zeros(round(tEnd / H) + 1)
        self.__x2d      = np.zeros(round(tEnd / H) + 1)
        self.__x2dd     = np.zeros(round(tEnd / H) + 1)
        self.__T        = np.zeros(round(tEnd / H) + 1)
        self.__V        = np.zeros(round(tEnd / H) + 1)
        self.__E        = np.zeros(round(tEnd / H) + 1)
        self.__u        = np.zeros((round(tEnd / H) + 1, 3))
        

        # _________________________________________________________________________
        #                                                       Connectivity matrix
        self.__L    = np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]])
        self.__yEx  = np.zeros(3)
        self.__uEx  = np.zeros(3)
        self.__y1   = np.zeros(4)
        self.__y2   = np.zeros(3)


    # _______________________________________________________ ASSIGN SUBSYSTEMS
    def assignSS1(self, SS1):
        self.__SS1 = SS1
        self.__SS1.setStepSize(self.__h1)
        print('SS1 assigned correctly')
    
    def assignSS2(self, SS2):
        self.__SS2 = SS2
        self.__SS2.setStepSize(self.__h2)
        print('SS2 assigned correctly')

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

    # _______________________________________________________ CORRECTIONS
    def correctionCoeffs(self, coeff, limit):
        self.__corrCoeffs               = coeff
        self.__corr                     = True
        self.__limit                    = abs(limit)
        self.__vec                      = np.ones(coeff.shape)

    # _______________________________________________________ STORAGE
    def store(self):
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

        # Evaluate magnitudes
        # Assemble arrays
        x  = np.array((self.__y1[1], self.__y2[0]))
        xd = np.array((self.__y1[2], self.__y2[1]))

        # Evaluation of kinetic and potential energy
        self.__T[self.__storeIdx] = 0.5 * self.__SYS.m[0] * (xd[0]**2) + 0.5 * self.__SYS.m[1] * (xd[1]**2)
        self.__V[self.__storeIdx] = 0.5 * self.__SYS.k[0] * (x[0]**2) + 0.5 * self.__SYS.k[1] * ((x[1])**2) + 0.5 * self.__SYS.kc * ((x[0] - x[1])**2)      
        self.__E[self.__storeIdx] = self.__T[self.__storeIdx] + self.__V[self.__storeIdx]
        
        # Update index
        self.__storeIdx += 1
    
    # _______________________________________________________ INITIALIZATION
    def initialize(self):
        # Iterative initialization
        y1It    = np.zeros(4)
        error   = 10
        tol     = 1e-10

        while error > tol:

            # Initialize units
            self.__SS1.initialize(0.0)
            self.__SS2.initialize(0.0)

            # Retrieve outputs of subsystems
            self.__y1 = self.__SS1.readOutputs()
            self.__y2 = self.__SS2.readOutputs()
        
            # Evaluation
            self.evalManager()

            # Send inputs
            self.__SS1.setInputs(self.__uEx[0:2])
            self.__SS2.setInputs(self.__uEx[2:])

            # Evaluate error
            error = np.linalg.norm(y1It - self.__y1, ord = 2)
            y1It = self.__y1

        # Store initial equilibrium
        self.store()

        # Update correction
        if self.__corr:
            self.__vec[2] = self.__uEx[2:]
            self.__vec[1] = self.__vec[2]

        return self.__SS1.readOutputs(), self.__SS2.readOutputs()

    # _______________________________________________________ MAIN LOOP
    def run(self):
        # The simulation loop follows a Jacobi non-iterative scheme
        while self.__tEnd - self.__time > self.__H / 2.0:

            # Report when required
            if self.__storeIdx % 1000 == 0:
                print('Macro time: ', round(self.__time, 4), 's')
           
           # Exchange coupling variables and call manager
            self.__y1 = self.__SS1.readOutputs()
            self.__y2 = self.__SS2.readOutputs()
            
            # Evaluation
            self.evalManager()

            # Update correction
            if self.__corr:
                self.__vec[2] = self.__vec[1]
                self.__vec[1] = self.__uEx[2:]
                corr = np.matmul(self.__corrCoeffs, self.__vec)
                
                # Correction limits
                #if corr > self.__limit:
                    #self.__uEx[2:] = self.__limit

                #elif corr < -self.__limit:
                    #self.__uEx[2:] = -self.__limit
                if corr - self.__uEx[2:] > self.__limit:
                    self.__uEx[2:] += self.__limit

                elif corr - self.__uEx[2:] < -self.__limit:
                    self.__uEx[2:] -= self.__limit

                else:
                    self.__uEx[2:] = corr

            # Store system states
            if self.__time > self.__H / 2.0:
                self.store()
                
            # Send inputs
            self.__SS1.setInputs(self.__uEx[0:2])
            self.__SS2.setInputs(self.__uEx[2:])
            
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
        if self.__storeIdx % 1000 == 0:
            print('Macro time: ', round(self.__time, 4), 's')
        
        # Exchange coupling variables and call manager
        self.__y1 = self.__SS1.readOutputs()
        self.__y2 = self.__SS2.readOutputs()
        
        # Evaluation
        self.evalManager()

        # Store system states
        if self.__time > self.__H / 2.0:
            self.store()

    # _______________________________________________________ TERMINATE
    def terminate(self):
        self.__SS1.terminate()
        self.__SS2.terminate()

    # _______________________________________________________ ACCESSORS
    def getSTORE(self):
        return  self.__t, self.__fc, self.__x1, self.__x1d, self.__x1dd, self.__x2, self.__x2d, self.__x2dd, self.__T, self.__V, self.__E, self.__u


