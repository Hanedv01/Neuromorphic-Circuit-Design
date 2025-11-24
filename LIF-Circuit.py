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
            #return 0
            raise ValueError(f"Vds should be non-negative! It is {Vds}")
        if Vgs <= self.Vth:                                             # Cut-off region
            return 0
        elif 0 <= Vds and Vds <= Vgs - self.Vth:                        # Linear region
            return self.K * ((Vgs-self.Vth)*Vds - (Vds**2)/2)*(1 + self.lam*Vds)
        elif 0 <= Vgs - self.Vth and Vgs - self.Vth <= Vds:             # Saturated region
            return (self.K/2)*((Vgs-self.Vth)**2)*(1 + self.lam*Vds)
        else:
            raise ValueError(f"Something is wrong with the arguments! Vgs={Vgs}V, Vds={Vds}V")
        
    
class LIF_Circuit:
    def __init__(self, Rin, Rl, Rd, C, HystComp, MOSFET):
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
        self.Rl = Rl
        self.Rd = Rd
        self.HystComp = HystComp
        self.MOSFET = MOSFET
    
    Vmem = 0
    Vout = 0
    
    spikeWidths = []
    timeBetweenSpikes = []
    isInSpike = False
    t_temp = 0

    def SolveRecursiveIds(self, Vout, Vmem):
        def FET_SignFixer(Ids):
            """
            Swaps source and drain if Vds < 0
            """
            Vd = Vmem
            Vs = Ids*self.Rd
            Vg = Vout
            Vds = Vd - Vs
            print(Vds)
            if Vds >= 0:    # If everything is as it should be
                Vgs = Vg - Vs
                return self.MOSFET.GetIds(Vgs, Vds)
            else:           # If the formula for Ids would get a negative argument
                Vgd = Vg - Vd
                Vsd = Vs - Vd
                Isd = self.MOSFET.GetIds(Vgd, Vsd)
                return -Isd
            
        def Ids_optionA(Ids):
            Vgs = self.Vout - Ids*self.Rd
            Vds = self.Vmem - Ids*self.Rd
            if abs(Vds) <= 1e-15:       # Necessary for stability
                Vds = 0
            #print(f"Vds = {Vds}")
            return self.MOSFET.GetIds(Vgs,Vds)

        def Ids_optionB(Ids):
            Vgs = self.Vout - self.Vmem
            Vds = Ids*self.Rd - self.Vmem
            #print(f"Vds = {Vds}")
            return self.MOSFET.GetIds(Vgs,Vds)
            
        def residual(Ids):
            return Ids - Ids_optionA(Ids)
            #return Ids - FET_SignFixer(Ids)
        
        Imax = Vmem/self.Rd
        #print(f"Imax = {Imax}")
        Ilow, Ihigh = -Imax, Imax
        sol = root_scalar(residual, bracket=[Ilow,Ihigh], method='brentq')

        if not sol.converged:
            raise RuntimeError("Ids solution did not converge.")
        return sol.root
    
    def Step(self, t, dt, Vin):
        def WidthUpdate(t):
            if self.isInSpike == False and self.Vout >= self.HystComp.Isat0:
                self.isInSpike = True
                self.timeBetweenSpikes.append(t - self.t_temp)
                self.t_temp = t
            if self.isInSpike == True and self.Vout < 0.5*self.HystComp.Isat0:
                self.isInSpike = False
                self.spikeWidths.append(t - self.t_temp)
                self.t_temp = t

        Ids = self.SolveRecursiveIds(self.Vout, self.Vmem)
        #print(f"Ids = {Ids}")
        if Vin < self.Vmem:
            dVmem = (-Ids/self.C - self.Vmem/self.Rl)*dt
        else:
            dVmem = ((Vin-self.Vmem)/(self.Rin*self.C) - self.Vmem/self.Rl - Ids/self.C)*dt
        self.Vmem += dVmem
        self.HystComp.Update(self.Vmem, t)
        self.Vout = self.HystComp.current
        WidthUpdate(t)
        #print(f"Vout = {self.Vout}")
        #print(f"Vmem = {self.Vmem}")
        #print(f"Vs = {Ids*self.Rd}")
        #print(f"Vds = {self.Vmem - Ids*self.Rd}")
        #print(f"Vgs = {self.Vout - Ids*self.Rd}")
        #print(f"Ids = {Ids}")
        return Ids



def main():
    def Stepfunction(t):
        if t >= 0:
            return 1
        else:
            return 0
        
    def SpikeTrain(t, tSpike, tPause):
        tTot = tSpike + tPause
        if t%tTot < tPause:
            return 0
        else:
            return 1
        
    def CircuitTest(T, N, plot=False):
        MOSFET = FET(0.6,0,1e-4)
        Vstart = 0.3
        Vsat = 0.6
        Isat = 1
        Vhl = 0.6
        Vlh = 0.3
        HystDevice = AFE_FET(Vstart,Vsat,Isat,1e3,1e3,Vhl,Vlh,1e-1,0,0,0)
        Circuit = LIF_Circuit(5e5, 1e12, 1e1, 1e-12, HystDevice, MOSFET)
        #print(Circuit.SolveRecursiveIds(0,0))
        dt = T/(N-1)
        tlist = np.linspace(0,T,N)
        Vinlist  = []
        Voutlist = []
        Vmemlist = []
        Vdslist  = []
        for t in tlist:
            #Vdslist.append(Circuit.Vmem - Circuit.Rd*Circuit.Step(t, dt, Stepfunction(t)))
            #print(t)
            #Vinlist.append(Stepfunction(t))
            Circuit.Step(t, dt, SpikeTrain(t, 1.57e-7, 5.6e-8))
            Vinlist.append(SpikeTrain(t, 1.57e-7, 5.6e-8))
            Voutlist.append(Circuit.Vout)
            Vmemlist.append(Circuit.Vmem)

        
        if plot:
            plt.plot(tlist, Vinlist,  label="Vin")
            plt.plot(tlist, Voutlist, label="Vout")
            plt.plot(tlist, Vmemlist, label="Vmem")
            #plt.plot(tlist, Vdslist,  label="Ids")
            plt.xlabel("t [s]")
            plt.ylabel("V [V]")
            plt.legend()
            plt.show()

        print(f"Width of spikes: {Circuit.spikeWidths}")
        print("")
        print(f"Time between spikes: {Circuit.timeBetweenSpikes}")
        print("")
        print(f"Frequency: {1/(Circuit.spikeWidths[-1] + Circuit.timeBetweenSpikes[-1]):.0f} Hz")
        
    CircuitTest(0.000002, 40001, plot = True)


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
