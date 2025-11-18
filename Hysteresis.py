import numpy as np


class AFE_FET:
    current = 0
    voltage = 0
    stateOn = False
    SpikeTimes = []
    def __init__(self, Vstart, Vsat, Isat, Roff, Ron, Vhl, Vlh, tau, deltaVon, deltaVoff, deltaI):
        """
        Represents a device with the I-V characteristics of a hysteresis loop. 
        
        """
        self.Vstart0 = Vstart
        self.Vsat0 = Vsat
        self.Isat0 = Isat
        self.Roff = Roff
        self.Ron = Ron
        self.Vhl0 = Vhl
        self.Vlh0 = Vlh
        self.tau = tau
        self.deltaVon = deltaVon
        self.deltaVoff = deltaVoff
        self.deltaI = deltaI

        if self.Vstart0 > self.Vhl0:
            raise ValueError("Vstart should be smaller than Vhl!")
        if self.Vhl0 > self.Vsat0:
            raise ValueError("Vhl should be smaller than Vsat!")
        if self.Vlh0 > self.Vsat0:
            raise ValueError("Vlh should be smaller than Vsat!")

    
    # ------------------------------------------------------------------------------
    #   Functions which define behaviour in different parts of the hysteresis loop
    # ------------------------------------------------------------------------------

    # Gives a straight line with variable slope and intercept
    def Line(self, V, slope, intersect):
        return slope*V + intersect
    
    # Gives smooth growth between two given points and given slope at V0
    def Smooth(self, V, V0, V1, I0, I1, slope):
        if not V0 == V1:
            a = (-slope*(V0-V1) + I0 - I1)/(V0*V1 - V1**2)
            b = slope - a*V0
            c = I0 - a*V0**2 - b*V0
            return a*V**2 + b*V + c
        else:
            raise ValueError("Make sure that V0 and V1 differ!")
    def AddSpike(self, tSpike):
        self.SpikeTimes = np.append(self.SpikeTimes, tSpike)

    def Shift(self, t, delta, tau):
        output = 0
        for tSpike in self.SpikeTimes:
            if t >= tSpike:
                output += delta * np.exp(-(t-tSpike) / tau)
        return output
    
    def ShiftNotLast(self, t, delta, tau):
        output = 0
        for tSpike in self.SpikeTimes[:-1]:
            if t >= tSpike:
                output += delta * np.exp(-(t-tSpike) / tau)
        return output


    def Update(self, Vnew, t):
        """
        Updates the current, given the past state and new voltage
        """
        Vold = self.voltage
        self.voltage = Vnew
        # Update voltages and currents according to shift
        Vstart = self.Vstart0 + self.Shift(t, self.deltaVoff, self.tau)
        Vsat   = self.Vsat0   + self.Shift(t, self.deltaVon, self.tau)
        Isat   = self.Isat0   + self.Shift(t, self.deltaI, self.tau)
        Vhl    = self.Vhl0    + self.Shift(t, self.deltaVon, self.tau)
        Vlh    = self.Vlh0    + self.Shift(t, self.deltaVoff, self.tau)
        # Create an upper current ceiling
        UpperLimit = self.Line(Vnew-Vsat, 1/self.Ron, Isat)
        # Corresponding to the baseline through the origin
        if Vnew <= Vstart:
            if self.stateOn:
                    self.stateOn = False
                    self.AddSpike(t)
            if Vnew <= self.Vstart0:
                self.current = self.Line(Vnew, 1/self.Roff, 0)
            # If shifted, connect smoothly between Vstart0 and Vstart
            else:
                if Vnew >= Vold:
                    self.current = self.Smooth(Vnew, self.Vstart0, Vstart, self.Line(self.Vstart0, 1/self.Roff, 0), self.Line(Vstart, 1/self.Roff, 0) + self.Shift(t, self.deltaI, self.tau), 1/self.Roff)
                else:
                    self.current = self.Smooth(Vnew, self.Vstart0, self.Vstart0+self.ShiftNotLast(t, self.deltaVoff, self.tau), self.Line(self.Vstart0, 1/self.Roff, 0), self.Line(self.Vstart0+self.ShiftNotLast(t, self.deltaVoff, self.tau), 1/self.Roff, 0) + self.ShiftNotLast(t, self.deltaI, self.tau), 1/self.Roff)
        # Correponding to the saturated current after Vsat
        elif Vnew >= Vsat:
            self.current = UpperLimit
            self.stateOn = True
        # Corresponding to the loop
        else:
            if not self.stateOn:
                if Vnew <= Vhl:     # Segment 1
                    self.current = self.Line(Vnew, 1/self.Roff, 0) + self.Shift(t, self.deltaI, self.tau)
                else:               # Segment 2
                    self.current = self.Line(Vnew-Vhl, (Isat-Vhl/self.Roff)/(Vsat-Vhl), Vhl/self.Roff + self.Shift(t, self.deltaI, self.tau))
                    if self.current >= UpperLimit:      #Sanity check
                        self.current = UpperLimit
            else:
                if Vnew >= Vlh:     # Segment 3
                    self.current = self.Line(Vnew-Vsat, 1/self.Ron, Isat)
                else:               # Segment 4
                    self.stateOn = True
                    self.current = self.Line(Vnew-Vstart, (self.Line(Vlh-Vsat, 1/self.Ron, Isat)-Vstart/self.Roff)/(Vlh-Vstart), Vstart/self.Roff+ self.Shift(t, self.deltaI, self.tau))
                    if self.current >= UpperLimit:      #Sanity check
                        self.current = UpperLimit
                        


# Flytta turn on till ett senare segment
# Alternativt flytta den till Vstart0 och brute forcea om vi har en switch mellan Vstart0 och Vstart


def main():
    import matplotlib.pyplot as plt
    x = AFE_FET(1, 7, 10, 1e1, 1e1, 5, 3, 1000, 1, 0, 0)
    """
    for tSpike in [0.3, 0.4, 0.8]:
        x.AddSpike(tSpike)

    tlist = np.linspace(0,1,100)
    Ilist = []
    for t in tlist:
        Ilist.append(x.Shift(t, 1, 0.1))
    
    plt.plot(tlist, Ilist)
    plt.show()
    print(x.SpikeTimes)
    """


    """
    Vlist = [0,1,2,3,4,5,6,5,6,7,8,7,6,5,4,3,2,1,0,-1,-2,-1,0,1,2,3,4,5,6,7,8,9,8,7,6,5,4,3,2,1,0]
    Ilist = []
    tlist = [0]
    for V in Vlist:
        t = tlist[-1] + 0.1
        tlist.append(t)
        x.Update(V, t)
        Ilist.append(x.current)

    plt.plot(Vlist, Ilist, color="blue")
    plt.show()
    """
    """
    def smallTest():
        A = AFE_FET(0,1,1,1e3,1e3,1,0,1000, 0,0)
        TP = 1.3
        Alist = np.linspace(0, TP, 100)
        Blist = np.linspace(TP, 0, 100)
        Vlist = np.concatenate((Alist, Blist))
        Ilist = []
        tlist = [0]
        for V in Vlist:
            t = tlist[-1] + 0.1
            tlist.append(t)
            A.Update(V, t)
            Ilist.append(A.current)
        plt.plot(Vlist, Ilist)
        plt.show()
    smallTest()

    """
    TP = 10
    Alist = np.linspace(0, TP, 100)
    Blist = np.linspace(TP, 0, 100)
    Vlist = np.concatenate((Alist, Blist))
    Vlist2 = np.concatenate((Alist, Blist, Alist, Blist))
    #Vlistbad = [0,1,1.2,1.4,1.6,1.8,1.6,1.4,1.2,1,0.9,1,1.2,1.4,1.6,1.8,2,3,4,5,6,7,8,7,6,5,4,3,2,1.8,1.6,1.4,1.2,1,0.9,1,1.2,1.4,1.6,1.8,2,2.1,2,1.8,1.6,1.4,1.2,1,0]
    #Vlist = [1,3,2,8,4,6,3,7,3,1]
    Ilist1 = []
    Ilist2 = []
    Ilist3 = []
    tlist = [0]
    for V in Vlist:
        t = tlist[-1] + 0.1
        tlist.append(t)
        x.Update(V, t)
        Ilist1.append(x.current)
    print(x.SpikeTimes)

    for V in Vlist:
        t = tlist[-1] + 0.1
        tlist.append(t)
        x.Update(V, t)
        Ilist2.append(x.current)
    print(x.SpikeTimes)

    for V in Vlist:
        t = tlist[-1] + 0.1
        tlist.append(t)
        x.Update(V, t)
        Ilist3.append(x.current)
    print(x.SpikeTimes)
    print(t)

    plt.plot(Vlist, Ilist1, color="red")
    plt.plot(Vlist, Ilist2, color="blue")
    plt.plot(Vlist, Ilist3, color="green")
    plt.show()

    

if __name__ == "__main__":
    main()