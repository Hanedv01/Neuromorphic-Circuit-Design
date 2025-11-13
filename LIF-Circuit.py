import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
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
        if Vgs <= self.Vth:     # Cut-off region
            return 0
        elif 0 <= Vds and Vds <= Vgs - self.Vth:    # Linear region
            return self.K * ((Vgs-self.Vth)*Vds - (Vds**2)/2)*(1 + self.lam*Vds)
        elif 0 <= Vgs - self.Vth and Vgs - self.Vth <= Vds:     # Saturated region
            return (self.K/2)*((Vgs-self.Vth)**2)*(1 + self.lam*Vds)
        else:
            raise ValueError(f"Something is wrong with the arguments! Vgs={Vgs}V, Vds={Vds}V")
        
    
class LIF_Circuit:
    def __init__(self, Rin, Rd, C, HystComp, MOSFET):
        """
        Components:
        -------------------
        Rin and Rd: float [Ohm]

        C: float [Farad]

        HystComp: AFE_FET object from Hysteresis.py

        FET: FET object
        """
        self.Rin = Rin
        self.C = C
        self.Rd = Rd
        self.HystComp = HystComp
        self.MOSFET = MOSFET
    
    Vmem = 0
    Vout = 0

    def SolveRecursiveIds(self, Vout, Vmem):
        def residual(Ids):
            Vgs = Vout - Vmem
            Vds = Vmem - Ids*self.Rd
            print(Vgs, Vds)
            print(self.MOSFET.GetIds(Vgs, Vds))
            return Ids - self.MOSFET.GetIds(Vgs, Vds)
        
        Imin, Imax = 0, 10
        sol = root_scalar(residual, bracket=[Imin, Imax], method='brentq')

        if not sol.converged:
            raise RuntimeError("Self-consistent Ids solution did not converge.")
        return sol.root
    
    def Step(self, t, dt, Vin):
        self.HystComp.Update(self.Vmem, t)
        self.Vout = self.HystComp.current
        print(self.Vout)
        print(self.Vmem)
        Ids = self.SolveRecursiveIds(self.Vout, self.Vmem)
        dVmem = ((Vin-self.Vmem)/(self.Rin*self.C) + Ids/self.C)*dt
        self.Vmem += dVmem



def main():
    def Stepfunction(t):
        if t >= 1:
            return 1
        else:
            return 0
        
    def CircuitTest(T, N):
        MOSFET = FET(1,0,1e-4)
        HystDevice = AFE_FET(1,8,10,1e3,1e3,2,7,100,1,1)
        Circuit = LIF_Circuit(1e3, 1e3, 1e-12, HystDevice, MOSFET)
        print(Circuit.SolveRecursiveIds(0,0))
        dt = T/(N-1)
        tlist = np.linspace(0,T,N)
        Ilist = []
        for t in tlist:
            Circuit.Step(t, dt, Stepfunction(t))
            Ilist.append(Circuit.Vout)
        plt.plot(tlist, Ilist)
        plt.show()
        
    CircuitTest(10, 501)

    def CapacitorCheck(C, N, T, VinFunc):
        C1 = Capacitor(C)
        dt = T/(N-1)
        tlist = np.linspace(0, T, N)
        Clist = []
        for t in tlist:
            dVin = VinFunc(t) - VinFunc(t-dt)
            C1.StepUpdate(dVin)
            Clist.append(C1.Charge())

        plt.plot(tlist, Clist)
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
        TestFET = FET(Vth, l, K)
        Vlist = np.linspace(0,10,201)
        Vrange = [0,1,2,4,6]
        Ilist = []
        for Vgs in Vlist:
            Ilist.append(TestFET.GetIds(Vgs, 5))
        #plt.plot(Vlist, Ilist)
        #plt.xlabel("Vgs")
        #plt.ylabel("Ids")
        #plt.show()

        for Vgs in Vrange:
            Ilist = []
            for Vds in Vlist:
                Ilist.append(TestFET.GetIds(Vgs, Vds))
            plt.plot(Vlist, Ilist)
            plt.title(f"Vgs = {str(Vgs)}")
            plt.xlabel("Vds")
            plt.ylabel("Ids")
            plt.show()
    #FETCheck(1.7,0,1e-4)



if __name__ == "__main__":
    main()
