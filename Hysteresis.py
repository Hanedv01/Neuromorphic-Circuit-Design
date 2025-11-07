import numpy as np
import matplotlib.pyplot as plt

class AFE_FET:
    current = 0
    voltage = 0
    stateOn = False
    SpikeTimes = []
    def __init__(self, Vstart, Vsat, Isat, Roff, Ron, Vhl, Vlh, tau, deltaV, deltaI):
        self.Vstart0 = Vstart
        self.Vsat0 = Vsat
        self.Isat0 = Isat
        self.Roff = Roff
        self.Ron = Ron
        self.Vhl0 = Vhl
        self.Vlh0 = Vlh
        self.tau = tau
        self.deltaV = deltaV
        self.deltaI = deltaI

        if self.Vstart0 > self.Vhl0:
            raise ValueError("Vstart should be smaller than Vhl!")
        if self.Vhl0 > self.Vsat0:
            raise ValueError("Vhl should be smaller than Vsat!")
        if self.Vlh0 > self.Vsat0:
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
    
    def AddSpike(self, tSpike):
        self.SpikeTimes = np.append(self.SpikeTimes, tSpike)

    def Shift(self, t, delta, tau):
        output = 0
        for tSpike in self.SpikeTimes:
            if t >= tSpike:
                output += delta * np.exp(-(t-tSpike) / tau)
        return output

    # Gives the new current, given past state and new voltage
    def Update(self, Vnew, t):
        self.Vold = self.voltage
        self.voltage = Vnew
        # Update voltages and currents according to shift
        Vstart = self.Vstart0 + self.Shift(t, self.deltaV, self.tau)
        Vsat   = self.Vsat0   + self.Shift(t, self.deltaV, self.tau)
        Isat   = self.Isat0   + self.Shift(t, self.deltaI, self.tau)
        Vhl    = self.Vhl0    + self.Shift(t, self.deltaV, self.tau)
        Vlh    = self.Vlh0    + self.Shift(t, self.deltaV, self.tau)
        print(Vhl)
        # Corresponding to the baseline through the origin
        if Vnew <= Vstart:
            self.current = self.Line(Vnew, 1/self.Roff, 0)
            if self.stateOn:
                self.stateOn = False
                self.AddSpike(t)
        # Correponding to the saturated current after Vsat
        elif Vnew >= Vsat:
            self.current = self.Line(Vnew-Vsat, 1/self.Ron, Isat)
            self.stateOn = True
        # Corresponding to the loop
        else:
            if not self.stateOn:
                if Vnew <= Vhl:
                    self.current = self.Line(Vnew, 1/self.Roff, 0) + self.Shift(t, self.deltaI, self.tau)
                else:
                    if Vnew >= self.Vold:
                        self.current = self.Line(Vnew-Vhl, (Isat - Vhl/self.Roff)/(Vsat-Vhl), Vhl/self.Roff)
                    else:
                        self.current = self.Line(Vnew-Vsat, 1/self.Ron, Isat)
                        self.stateOn = True
            else:
                if Vnew >= Vlh:
                    self.current = self.Line(Vnew-Vsat, 1/self.Ron, Isat)
                else:
                    if Vnew <= self.Vold:
                        self.current = self.Line(Vnew-Vstart, (self.Line(Vlh-Vsat, 1/self.Ron, Isat)-Vstart/self.Roff)/(Vlh-Vstart), Vstart/self.Roff)
                        self.stateOn = True

        
"""
class TestCicuit:
    current = 0
    time = 0
    shift = 0
    SpikeTimes = np.array([])
    def __init__(self, tau):
        self.tau = tau

    def AddSpike(self, tSpike):
        self.SpikeTimes = np.append(self.SpikeTimes, tSpike)

    def Shift(self, t, delta):
        output = 0
        for tSpike in self.SpikeTimes:
            if t >= tSpike:
                output += delta * np.exp(-(t-tSpike) / self.tau)
        return output
         

tlist = np.linspace(0,1,1000)

circuit = TestCicuit(0.1)
Ilist = []
for ind in range(len(tlist)):
    if ind%200 == 0 and ind<=500:
        circuit.AddSpike(tlist[ind])
    Ilist.append(circuit.Shift(tlist[ind], 2))
    print(f"{tlist[ind]},   {Ilist[-1]},   {circuit.SpikeTimes}")

plt.plot(tlist, Ilist)
plt.show()

"""

x = AFE_FET(1, 7, 100, 1e3, 1e3, 5, 3, 100, 0.1, 5)

TP = 8
Alist = np.linspace(0, TP, 100)
Blist = np.linspace(TP, 0, 100)
Vlist = np.concatenate((Alist, Blist, Alist, Blist))
#Vlist = [0,1,2,3,4,5,6,7,8,7,6,5,4,3,2,1,0]
#Vlist = [1,3,2,8,4,6,3,7,3,1]
Ilist = []
tlist = [0]
for V in Vlist:
    t = tlist[-1] + 0.1
    x.Update(V, t)
    Ilist.append(x.current)

print(x.SpikeTimes)

#print(Ilist)
plt.close()
plt.plot(Vlist, Ilist)
plt.show()
