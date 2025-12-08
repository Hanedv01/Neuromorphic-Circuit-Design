import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from Hysteresis import AFE_FET

# För logaritmisk z-skala
from matplotlib.colors import LogNorm

class FET:
    def __init__(self, Vth, lam, K):
        self.Vth = Vth
        self.lam = lam
        self.K = K

    #Based on Shichman-Hodges
    def GetIds(self, Vgs, Vds):
        """Returns the current Ids"""
        if Vds < 0:
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
    def __init__(self, Rin, Rl, Rs, C, HystComp, MOSFET):
        """
        Components:
        -------------------
        Rin and Rs: float [Ohm]

        C: float [Farad]

        HystComp: AFE_FET object from Hysteresis.py

        FET: FET object
        """
        self.Rin = Rin
        self.C = C
        self.Rl = Rl
        self.Rs = Rs
        self.HystComp = HystComp
        self.MOSFET = MOSFET
    
    Vmem = 0
    Vout = 0
    
    spikeWidths = []
    timeBetweenSpikes = []
    isInSpike = False
    t_temp = 0

    def SolveRecursiveIds(self, Vmem, Vout):          
        def GetCircuitIds(Ids):
            Vgs = Vout - Ids*self.Rs
            Vds = Vmem - Ids*self.Rs
            if abs(Vds) <= 1e-15:       # Necessary for stability
                Vds = 0
            #print(f"Vds = {Vds}")
            return self.MOSFET.GetIds(Vgs,Vds)
            
        def residual(Ids):
            return Ids - GetCircuitIds(Ids)
            #return Ids - FET_SignFixer(Ids)
        
        Imax = Vmem/self.Rs
        #print(f"Imax = {Imax}")
        Ilow, Ihigh = -Imax, Imax
        sol = root_scalar(residual, bracket=[Ilow,Ihigh], method='brentq')

        if not sol.converged:
            raise RuntimeError("Ids solution did not converge.")
        return sol.root
    
    def GetAsymptote(self, Vin):
        """
        Requires relatively large Ron!
        """
        storedTau = self.HystComp.tau
        self.HystComp.tau = 1e20

        def FuncVmem(Vmem):
            return (Vin/self.Rin - self.SolveRecursiveIds(Vmem, self.HystComp.Isat0))/(1/self.Rin + 1/self.Rl)

        def residual(Vmem):
            return Vmem - FuncVmem(Vmem)
        
        sol = root_scalar(residual, bracket=[0,Vin], method='brentq')
        if not sol.converged:
            raise RuntimeError("Asymptotic Vmem solution did not converge.")
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

        Ids = self.SolveRecursiveIds(self.Vmem, self.Vout)
        #print(f"Ids = {Ids}")
        if Vin < self.Vmem:
            dVmem = (-Ids/self.C - self.Vmem/(self.Rl*self.C))*dt
        else:
            dVmem = ((Vin-self.Vmem)/(self.Rin*self.C) - self.Vmem/(self.Rl*self.C) - Ids/self.C)*dt
        self.Vmem += dVmem
        self.HystComp.Update(self.Vmem, t)
        self.Vout = self.HystComp.current
        WidthUpdate(t)
        #print(f"Vout = {self.Vout}")
        #print(f"Vmem = {self.Vmem}")
        #print(f"Vs = {Ids*self.Rs}")
        #print(f"Vds = {self.Vmem - Ids*self.Rs}")
        #print(f"Vgs = {self.Vout - Ids*self.Rs}")
        #print(f"Ids = {Ids}")
        return Ids



def main():
    def Stepfunction(t):
        if t >= 0:
            return 1
        else:
            return 0
        
    def SpikeTrain(t):
        tSpike = 3e-7
        tPause = 27e-7
        tTot = tSpike + tPause
        if t%tTot < tPause:
            return 0
        else:
            return 1
        
    def CircuitTest(Circuit, T, N, VinFunc, plot=False, sideBySide=False):
        dt = T/(N-1)
        tlist = np.linspace(0,T,N)
        Vinlist  = []
        Voutlist = []
        Vmemlist = []
        Vdslist  = []

        for t in tlist:
            Vin = VinFunc(t)
            #if t > 3e-6 and t < 10e-6:
            #    Vin = 0
            Vdslist.append(Circuit.Vmem - Circuit.Rs*Circuit.Step(t, dt, Vin))
            Vinlist.append(Vin)
            Voutlist.append(Circuit.Vout)
            Vmemlist.append(Circuit.Vmem)

        if plot:
            tlist = tlist * 1e6
            if not sideBySide:
                plt.plot(tlist, Vinlist,  label="Vin")
                plt.plot(tlist, Voutlist, label="Vout")
                plt.plot(tlist, Vmemlist, label="Vmem")
                #plt.plot(tlist, Vdslist,  label="Ids")
                plt.xlabel("t [μs]", fontsize=16)
                plt.ylabel("V [V]", fontsize=16)
                plt.xticks(fontsize=14)
                plt.yticks(fontsize=14)
                plt.subplots_adjust(bottom=0.125)
                plt.title("Neuron dynamics: hysteresis device never turns off", fontsize=14)
                #plt.legend(loc=(0.04,0.68), fontsize=14)
                plt.legend(fontsize=14)
            else:
                fig, axs = plt.subplots(3, sharey=True)
                fig.suptitle('Neuron dynamics: simulated adaptivity')
                axs[0].plot(tlist, Vinlist)
                axs[1].plot(tlist, Vmemlist, 'tab:orange')
                axs[2].plot(tlist, Voutlist, 'tab:green')
                for ax in axs.flat:
                    ax.set(xlabel='t [μs]', ylabel='V [V]')

                # Hide x labels and tick labels for top plots and y ticks for right plots.
                for ax in axs.flat:
                    ax.label_outer()
            plt.show()

        if not len(Circuit.spikeWidths) == 0 and not len(Circuit.timeBetweenSpikes) == 0:
            #print(f"Width of spikes: {Circuit.spikeWidths}")
            #print("")
            #print(f"Time between spikes: {Circuit.timeBetweenSpikes}")
            #print("")
            frequency = 1/(Circuit.spikeWidths[-1] + Circuit.timeBetweenSpikes[-1])
            #print(f"Frequency: {frequency:.0f} Hz")
            return [Circuit.spikeWidths[-1], Circuit.timeBetweenSpikes[-1], frequency]
        else:
            return [0,0,0]
        
    
    MOSFET = FET(0.6,0,3e-5)
    Vstart = 0.1
    Vsat = 0.6
    Isat = 1
    Vhl = 0.6
    Vlh = 0.1
    HystDevice = AFE_FET(Vstart,Vsat,Isat,1e3,1e3,Vhl,Vlh,1e-2,0.02,0.01,0)
    Circuit = LIF_Circuit(4e5, 4e6, 1e1, 1e-12, HystDevice, MOSFET)
    #print(CircuitTest(Circuit, 0.000032, 40001, plot = True, sideBySide=True))
    #print(Circuit.GetAsymptote(1))
    

    def AsymptoteSearch(parameter):
        Vstart = 0.1
        Vsat = 0.7
        Isat = 1
        Vhl = 0.7
        Vlh = 0.1
        HystDevice = AFE_FET(Vstart,Vsat,Isat,1e3,1e3,Vhl,Vlh,1e-1,0,0,0)
        resultlist = []
        MOSFET = FET(0.7,0,2e-5)
        if parameter == "Rin":
            Rlist = np.logspace(2, 13, 100)
            for R in Rlist:
                Circuit = LIF_Circuit(R, 1e32, 1e2, 1e-12, HystDevice, MOSFET)
                resultlist.append(Circuit.GetAsymptote(1))
            plt.plot(Rlist, resultlist)
            plt.title("On-state asymptote: dependence on " + r"$\mathregular{R_{in}}$", fontsize=14)
            plt.xlabel(r"$\mathregular{R_{in}}$" +" [Ω]", fontsize=16)
            plt.subplots_adjust(bottom=0.125)
            plt.xscale("log")
        elif parameter == "Rl":
            Rlist = np.logspace(1, 18, 100)
            for R in Rlist:
                Circuit = LIF_Circuit(3e6, R, 1e3, 1e-12, HystDevice, MOSFET)
                resultlist.append(Circuit.GetAsymptote(1))
            plt.plot(Rlist, resultlist)
            plt.xlabel("Rl [Ohm]")
            plt.xscale("log")
        elif parameter == "C":      #No effect
            Clist = np.linspace(1e-11, 1e-13, 100)
            for C in Clist:
                Circuit = LIF_Circuit(3e6, 1e12, 1e3, C, HystDevice, MOSFET)
                resultlist.append(Circuit.GetAsymptote(1))
            plt.plot(Clist, resultlist)
            plt.xlabel("C [F]")
        elif parameter == "K":
            Klist = np.logspace(-9, -1, 100)
            for K in Klist:
                MOSFET = FET(0.6,0,K)
                Circuit = LIF_Circuit(3e6, 1e12, 1e3, 1e-12, HystDevice, MOSFET)
                resultlist.append(Circuit.GetAsymptote(1))
            plt.plot(Klist, resultlist)
            plt.title("On-state asymptote: dependence on K", fontsize=14)
            plt.xlabel("K [A/V$^2$]", fontsize=16)
            plt.subplots_adjust(bottom=0.125)
            plt.xscale("log")
        elif parameter == "Vth":
            Vlist = np.linspace(0.1, 0.9, 100)
            for Vth in Vlist:
                MOSFET = FET(Vth,0,2e-5)
                Circuit = LIF_Circuit(3e6, 1e12, 1e3, 1e-12, HystDevice, MOSFET)
                resultlist.append(Circuit.GetAsymptote(1))
            plt.plot(Vlist, resultlist)
            plt.xlabel("Vth [V]")
        else:
            print(f"Please enter a valid parameter! You entered {parameter}")
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.ylabel("Asymptote [V]", fontsize=16)
        plt.show()
    #AsymptoteSearch("Rin")



    def GridSearch():
        def Ratio(result):
            if abs(result[0]) > 1e-12:
                #print(result[0])
                ratio = result[1]/result[0]
            else:
                ratio = 0
            return ratio
        VlhList = np.linspace(0, 1, 21)
        VhlList = np.linspace(0, 1, 21)
        MOSFET = FET(0.7,0,2e-5)
        reslist = []
        for Vlh in VlhList:
            for Vhl in VhlList:
                if Vhl > Vlh:
                    print(f"Vlh = {Vlh}, Vhl = {Vhl}")
                    HystDevice = AFE_FET(Vlh,Vhl,1,1e3,1e3,Vhl,Vlh,1e-1,0,0,0)
                    Circuit = LIF_Circuit(7e6, 1e32, 1e2, 1e-11, HystDevice, MOSFET)
                    Circuit.spikeWidths = []
                    Circuit.timeBetweenSpikes = []
                    result = CircuitTest(Circuit, 0.001, 40001, plot = False)
                    print(result)
                    ratio = Ratio(result)
                    if not ratio == 0:
                        if ratio < 8:
                            Rin = Circuit.Rin*9/ratio
                            print(f"Ratio was < 8 ({ratio}). Attempting Rin = {Rin}.")
                            Circuit = LIF_Circuit(Rin, 1e32, 1e2, 1e-11, HystDevice, MOSFET)
                            Circuit.spikeWidths = []
                            Circuit.timeBetweenSpikes = []
                            result = CircuitTest(Circuit, 0.001, 40001, plot = False)
                            ratio = Ratio(result)
                        while ratio > 10:
                            Rin = Circuit.Rin*9/ratio
                            print(f"Ratio was > 10 ({ratio}). Attempting Rin = {Rin}.")
                            Circuit = LIF_Circuit(Rin, 1e32, 1e2, 1e-11, HystDevice, MOSFET)
                            Circuit.spikeWidths = []
                            Circuit.timeBetweenSpikes = []
                            result = CircuitTest(Circuit, 0.001, 40001, plot = False)
                            ratio = Ratio(result)
                    reslist.append([Vlh, Vhl, result[2]])
                    print(f"The final ratio was {ratio}.")
                    
        reslist = np.array(reslist)
        #print(reslist)
        x = reslist[:,0]
        y = reslist[:,1]
        z = reslist[:,2]   

        x_unique = np.unique(x)
        y_unique = np.unique(y)
        # create a grid initialized with NaN
        heat = np.full((len(y_unique), len(x_unique)), np.nan)

        # fill the grid
        for xi, yi, zi in reslist:
            ix = np.where(x_unique == xi)[0][0]
            iy = np.where(y_unique == yi)[0][0]
            heat[iy, ix] = zi

        # plot heatmap
        plt.imshow(heat, origin='lower', 
                extent=[x_unique.min(), x_unique.max(),
                        y_unique.min(), y_unique.max()],
                aspect='auto')
                #,norm=LogNorm())

        plt.colorbar(label="Frequency [Hz]")
        plt.xlabel("Vlh [V]")
        plt.ylabel("Vhl [V]")
        plt.title("Frequency heatmap")
        plt.show()
    #GridSearch()



    def FET_VgsSweep():
        Vgslist = [0, 2, 4, 6]
        Vdslist = np.linspace(0, 10, 100)
        TestFET = FET(0.6, 0, 2e-5)
        for Vgs in Vgslist:
            Ilist = []
            for Vds in Vdslist:
                Ilist.append(TestFET.GetIds(Vgs, Vds)*1000)
            plt.plot(Vdslist, Ilist, label=f"Vgs = {Vgs}")
            textposI = Ilist[-1] + 0.002
            plt.text(Vdslist[-1], textposI, '$V_{GS}$'+f' = {Vgs}', fontsize=12, horizontalalignment='right')
        plt.xlabel(r"$\mathregular{V_{DS}}$ [V]", fontsize=14)
        plt.ylabel(r"$\mathregular{I_{DS}}$ [μA]", fontsize=14)
        plt.title("$I_{DS}$ dependence on $V_{DS}$ for different $V_{GS}$")
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.show()
    #FET_VgsSweep()
            


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
    #FETCheck(0.6,0,2e-5)

    def VinPlot():
        MOSFET = FET(0.7,0,2e-5)
        Vstart = 0.1
        Vsat = 0.7
        Isat = 1
        Vhl = 0.7
        Vlh = 0.1
        HystDevice = AFE_FET(Vstart,Vsat,Isat,1e3,1e3,Vhl,Vlh,1e-2,0,0,0)
        Circuit = LIF_Circuit(7e6, 1e32, 1e2, 1e-11, HystDevice, MOSFET)

        VinList = np.linspace(0.6, 0.8, 11)
        freqList = []
        for Vin in VinList:
            Circuit = LIF_Circuit(7e6, 1e32, 1e2, 1e-11, HystDevice, MOSFET)
            Circuit.spikeWidths = []
            Circuit.timeBetweenSpikes = []
            def StepScaled(t):
                return Vin*Stepfunction(t)
            result = CircuitTest(Circuit, 0.1, 40001, StepScaled, plot = False)
            if result[2] == 0:
                freq = None
            else:
                freq = result[2]/1e3
            freqList.append(freq)
            print(f"Finished Vin = {Vin}")
        print(f"Vinlist = {VinList}")
        print(f"Frequencies = {freqList}")
        plt.plot(VinList, freqList, '.')
        plt.ylim(bottom=0)
        plt.xlim(left=0)
        plt.ylabel("Frequency [kHz]")
        plt.xlabel(r"$\mathregular{V_{in}}$"+" [V]")
        plt.title("Frequency as a function of the in-signal amplitude")
        plt.show()
    #VinPlot()




if __name__ == "__main__":
    main()
