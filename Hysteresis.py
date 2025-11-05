import numpy as np
import matplotlib.pyplot as plt


class Hysteresis:
    def __init__(self, V_start, V_HL, V_LH, R_down, R_up, R_off, R_on):
        # Store parameters
        self.V_start = V_start
        self.V_HL = V_HL
        self.R_up = R_up
        self.R_off = R_off
        self.R_on = R_on
        self.R_down = R_down
        self.V_LH = V_LH
        self.V_old = 0
        self.I_old = 0
        self.V = 0
        self.I = 0
        self.state = "off"

    def Update(self, V):
        """Simulate the hysteresis loop and store V and I."""
        self.V_old = self.V
        self.V = V
        self.I_old = self.I


        if self.state == "off" and self.V <= self.V_HL:
            self.I = (self.V-self.V_start)/self.R_off

        elif self.state == "off" and self.V > self.V_HL:
            self.I = (self.V_HL-self.V_start)/self.R_off + (self.V - self.V_HL)/self.R_up
            self.state = "up"

        elif self.state == "up" and self.V > self.V_old:
            self.I = self.I_old + (self.V-self.V_old)/self.R_up

        elif self.state == "up" and self.V < self.V_old:
            self.I = (self.V_HL-self.V_start)/self.R_off + (self.V_old - self.V_HL)/self.R_up + (self.V-self.V_old)/self.R_on
            self.state = "on"



        elif self.state == "on" and self.V >= self.V_HL and self.V >= self.V_LH and self.V > self.V_old:
            if self.V_old >= self.V_HL:
                self.I = self.I_old + (self.V-self.V_old)/self.R_up
                self.state = "up"
            else:
                self.I = self.I_old + (self.V_HL-self.V_old)/self.R_on + (self.V-self.V_HL)/self.R_up
                self.state = "up"

        elif self.state == "on" and self.V >= self.V_LH and self.V <= self.V_old:
            self.I = self.I_old + (self.V-self.V_old)/self.R_on
        
        elif self.state == "on" and self.V <= self.V_LH:
            if self.V_old >= self.V_LH:
                self.I = self.I_old + (self.V_LH - self.V_old)/self.R_on + (self.V - self.V_LH)/self.R_down
                self.state = "down"

        elif self.state == "down" and self.V <= self.V_old:
            if self.I_old + (self.V-self.V_old)/self.R_down > (self.V-self.V_start)/self.R_off:
                self.I = self.I_old + (self.V-self.V_old)/self.R_down
            else:
                self.state = "off"
                self.I = (self.V-self.V_start)/self.R_off

        elif self.state == "down" and self.V >= self.V_old:
            self.I = self.I_old + (self.V-self.V_old)/self.R_off
            self.state = "off" 
        else:
            print("Hoppsan!")


x = Hysteresis(0, 5, 2, 0.1, 0.1, 1, 1)
#Vlist = [-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,7,6,5,4,3,2,1,0,-1,-2,-3,-4]
Vlist = [0, 5, 10, 2, -3, 1, 3]
Ilist = []
for V in Vlist:
    x.Update(V)
    Ilist.append(x.I)
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
        if Vnew <= self.Vstart:
            self.current = 0
        else:
            if Vnew > self.Vold:
                if Vnew <= self.Vhl:
                    self.current = self.LinStep(self.Vold, Vnew, self.Iold, 1/self.Roff)
                else:
                    self.current = self.LinStep(self.Vold, Vnew, self.Iold, self.SlopeUp)
            else:
                if Vnew > self.Vlh:
                    self.current = self.LinStep(self.Vold, Vnew, self.Iold, 1/self.Ron)
                else:
                    self.current = self.LinStep(self.Vold, Vnew, self.Iold, self.SlopeDown)
        



x = AFE_FET(0, 1, 1, 5, 2, 10, 10)
Vlist = [0,1,2,3,4,5,6,7,8,7,6,5,4,3,2,1,0]
Ilist = []
for V in Vlist:
    x.Update(V)
    Ilist.append(x.current)

plt.close()
plt.plot(Vlist, Ilist)
plt.show()


