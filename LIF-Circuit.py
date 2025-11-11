import numpy as np
import matplotlib.pyplot as plt
from Hysteresis import AFE_FET

class Resistor:
    def __init__(self, R):
        """Represents an ideal resistor."""
        self.R = R

    def current(self, V):
        return V/self.R

# Behöver sannolikt ha spänning som in- och utdata
class Capacitor:
    voltage = 0
    current = 0
    def __init__(self, C):
        """Represents an ideal capacitor."""
        self.C = C

    def Charge(self):
        """Gives the charge on the capacitor."""
        return self.voltage * self.C
    
    def Update(self, Vnew, dt):
        """
        Updates the current through the capacitor, given a voltage.
        """
        Vold = self.voltage
        dV = Vnew-Vold
        self.current = self.C * dV/dt
        self.voltage = Vnew


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
        for t in tlist:
            Vin = VinFunc(t)
            dVc = (Vin-Cap.voltage)*dt / (R*C)
            Cap.Update(Vin, dt)
            Ilist.append(C*dVc/dt)
        plt.plot(tlist, Ilist)
        plt.show()
    RC_CircuitCheck(1e-6, 1e6, 1001, 2, Stepfunction)

    """

    def StepFunction(t):
        if t > 2:
            return 1
        else:
            return 0
        
    Circuit = RC_Circuit(1e3, 1e-9)
    T = 2.1
    dt = 0.1
    N = int(T/dt)
    tlist = np.linspace(0,T,N)
    Ilist = []
    Vlist = []
    for t in tlist:
        V = StepFunction(t)
        Vlist.append(V)
        Circuit.Update(V, dt)
        Ilist.append(Circuit.C.current)
    
    plt.plot(tlist, Ilist)
    plt.show()
    """



if __name__ == "__main__":
    main()
