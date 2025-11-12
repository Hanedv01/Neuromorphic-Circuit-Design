import numpy as np
import matplotlib.pyplot as plt
from Hysteresis import AFE_FET

class Resistor:
    def __init__(self, R):
        """Represents an ideal resistor."""
        self.R = R

    def current(self, V):
        return V/self.R

class Capacitor:
    voltage = 0
    def __init__(self, C):
        """Represents an ideal capacitor."""
        self.C = C

    def Charge(self):
        """Gives the charge on the capacitor."""
        return self.voltage * self.C
    
    def StepUpdate(self, dV):
        """
        Updates the voltage over the capacitor
        """
        self.voltage += dV

class FET:
    def __init__(self, Vth, lam, K):
        self.Vth = Vth
        self.lam = lam
        self.K = K

    #Based on Shichman-Hodges
    def GetIds(self, Vgs, Vds):
        """Returns the current Ids"""
        if Vds < 0:
            raise ValueError("Vds should be non-negative!")
        if Vds <= self.Vth:     # Cut-off region
            return 0
        elif 0 < Vds and Vds < Vgs - self.Vth:
            return self.K * ((Vgs-self.Vth)*Vds - (Vds**2)/2)*(1 + self.lam*Vds)
        elif 0 < Vgs - self.Vth and Vgs - self.Vth < Vds:
            return (self.K/2)*((Vgs-self.Vth)**2)*(1 + self.lam*Vds)
        else:
            raise ValueError("something is wrong with the arguments!")


def main():
    def Stepfunction(t):
        if t >= 1:
            return 1
        else:
            return 0

    def CapacitorCheck(C, N, T, VinFunc):
        C1 = Capacitor(C)
        dt = T/(N-1)
        tlist = np.linspace(0, T, N)
        Ilist = []
        for t in tlist:
            Vin = VinFunc(t)
            C1.Update(Vin,dt)
            Ilist.append(C1.current)

        plt.plot(tlist, Ilist)
        plt.show()
    #CapacitorCheck(1e-12, 201, 10, np.sin)

    def RC_CircuitCheck(C, R, N, T, VinFunc):
        Cap = Capacitor(C)
        Res = Resistor(R)
        dt = T/(N-1)
        tlist = np.linspace(0, T, N)
        Ilist = []
        Vclist = []
        for t in tlist:
            Vin = VinFunc(t)
            dVc = (Vin-Cap.voltage)*dt / (R*C)
            Cap.StepUpdate(dVc)
            Ilist.append(C* dVc/dt)
            Vclist.append(Cap.voltage)
        plt.plot(tlist, Vclist)
        plt.show()
    #RC_CircuitCheck(1e-6, 1e6, 1001, 10, Stepfunction)

    def FETCheck(Vth, l, K):
        TestFET = FET(1, 0, 1e-4)




if __name__ == "__main__":
    main()
