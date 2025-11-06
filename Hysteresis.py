import numpy as np
import matplotlib.pyplot as plt

class AFE_FET:
    current = 0
    voltage = 0
    stateOn = False
    def __init__(self, Vstart, Vsat, Isat, Roff, Ron, Vhl, Vlh):
        self.Vstart = Vstart
        self.Vsat = Vsat
        self.Isat = Isat
        self.Roff = Roff
        self.Ron = Ron
        self.Vhl = Vhl
        self.Vlh = Vlh
        if self.Vstart > self.Vhl:
            raise ValueError("Vstart should be smaller than Vhl!")
        if self.Vhl > self.Vsat:
            raise ValueError("Vhl should be smaller than Vsat!")
        if self.Vlh > self.Vsat:
            raise ValueError("Vlh should be smaller than Vsat!")

    def getSlopeUp(self):
        return (self.Isat - self.Vhl/self.Roff)/(self.Vsat-self.Vhl)
    
    def getSlopeDown(self):
        return (self.Line(self.Vlh-self.Vsat, 1/self.Ron, self.Isat)-self.Vstart/self.Roff)/(self.Vlh-self.Vstart)
    
    
    
    # ------------------------------------------------------------------------------
    #   Functions which define behaviour in different parts of the hysteresis loop
    # ------------------------------------------------------------------------------

    # Gives a straight line with variable slope and intercept
    def Line(self, V, slope, intersect):
        return slope*V + intersect
    

    def Update(self, Vnew):
        self.Vold = self.voltage
        self.voltage = Vnew
        # Corresponding to the baseline through the origin
        if Vnew <= self.Vstart:
            self.current = self.Line(Vnew, 1/self.Roff, 0)
            if not self.stateOn:
                self.stateOn = False
                # One loop is done! Update V and I!
        # Correponding to the saturated current after Vsat
        elif Vnew >= self.Vsat:
            self.current = self.Line(Vnew-self.Vsat, 1/self.Ron, self.Isat)
            self.stateOn = True
        # Corresponding to the loop
        else:
            if not self.stateOn:
                if Vnew <= self.Vhl:
                    self.current = self.Line(Vnew, 1/self.Roff, 0)
                else:
                    if Vnew >= self.Vold:
                        self.current = self.Line(Vnew-self.Vhl, (self.Isat - self.Vhl/self.Roff)/(self.Vsat-self.Vhl), self.Vhl/self.Roff)
                    else:
                        self.current = self.Line(Vnew-self.Vsat, 1/self.Ron, self.Isat)
                        self.stateOn = True
            else:
                if Vnew >= self.Vlh:
                    self.current = self.Line(Vnew-self.Vsat, 1/self.Ron, self.Isat)
                else:
                    if Vnew <= self.Vold:
                        self.current = self.Line(Vnew-self.Vstart, (self.Line(self.Vlh-self.Vsat, 1/self.Ron, self.Isat)-self.Vstart/self.Roff)/(self.Vlh-self.Vstart), self.Vstart/self.Roff)
                        self.stateOn = True

        



x = AFE_FET(1, 7, 100, 1e3, 1e3, 5, 3)
TP = 8
Alist = np.linspace(0, TP, 100)
Blist = np.linspace(TP, 0, 100)
Vlist = np.concatenate((Alist, Blist, Alist, Blist))
#Vlist = [0,1,2,3,4,5,6,7,8,7,6,5,4,3,2,1,0]
#Vlist = [1,3,2,8,4,6,3,7,3,1]
Ilist = []
for V in Vlist:
    x.Update(V)
    Ilist.append(x.current)

#print(Ilist)
plt.close()
plt.scatter(Vlist, Ilist)
plt.show()