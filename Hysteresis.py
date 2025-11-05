import numpy as np
import matplotlib.pyplot as plt

class AFE_FET:
    current = 0
    voltage = 0
    def __init__(self, Vstart, Roff, Ron, Vhl, Vlh, SlopeUp, SlopeDown):
        self.Vstart = Vstart
        self.Roff = Roff
        self.Ron = Ron
        self.Vhl = Vhl
        self.Vlh = Vlh
        self.SlopeUp = SlopeUp
        self.SlopeDown = SlopeDown
    
    # ------------------------------------------------------------------------------
    #   Functions which define behaviour in different parts of the hysteresis loop
    # ------------------------------------------------------------------------------

    # A linear stepper with variable slope
    def LinStep(self, Vold, Vnew, Iold, slope):
        return slope*(Vnew-Vold)  + Iold

    def Update(self, Vnew):
        self.Vold = self.voltage
        self.Iold = self.current
        self.voltage = Vnew
        if Vnew <= self.Vstart and abs(self.Iold) < 1e-1:
            self.current = 0
        else:
            if Vnew > self.Vold:
                if Vnew <= self.Vhl:
                    self.current = self.LinStep(self.Vold, Vnew, self.Iold, 1/self.Roff)
                else:
                    self.current = self.LinStep(self.Vold, Vnew, self.Iold, self.SlopeUp)
            else:
                if Vnew >= self.Vlh:
                    self.current = self.LinStep(self.Vold, Vnew, self.Iold, 1/self.Ron)
                else:
                    self.current = self.LinStep(self.Vold, Vnew, self.Iold, self.SlopeDown)
        



x = AFE_FET(1, 1, 1, 5, 2, 10, 10)
#Vlist = [0,1,2,3,4,5,6,7,8,7,6,5,4,5,6,7,8,7,6,5,4,3,2,1,2,1,0]
Vlist = [0,1,2,3,4,5,6,7,8,7,6,5,4,3,2,1,0.7,0.5,0.3,0]
Ilist = []
for V in Vlist:
    x.Update(V)
    Ilist.append(x.current)

print(Ilist)
plt.close()
plt.plot(Vlist, Ilist)
plt.show()