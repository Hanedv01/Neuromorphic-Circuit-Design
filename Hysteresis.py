import numpy as np
import matplotlib.pyplot as plt



class Hysteresis:
    def __init__(self, V_start, V_HL, R_down, R_up, R_off, R_on, I_sat, delta_down, delta_up, delta_right, delta_left, tau, k = 3):
        # Store parameters
        self.V_start = V_start #where the left side meets the off state
        self.V_HL = V_HL #Default threshold voltage
        self.R_up = R_up #resistance on right side
        self.R_off = R_off #resistance on bottom side
        self.R_on = R_on #resistance on top side
        self.R_down = R_down #resistance on left side
        self.V = 0 #Voltage at this time step
        self.I = 0 #current at this time step
        self.state = "off" #which of the four states the Neuron is in
        self.I_sat = I_sat #default saturation current
        self.delta_down = delta_down #pushing the bottom upwards (changes off current)
        self.delta_up = delta_up #pusing the top upwards (changes on current)
        self.delta_right = delta_right #pushing right side to the right (changes V_HL)
        self.delta_left = delta_left  #pushing left side to the right (changes V_start)
        self.decay = 0  #amounts of decays
        self.tau = tau
        self.k = k

    def Update(self, V, dt):
        """Simulate the hysteresis loop and store V and I."""
        self.V_old = self.V
        self.V = V
        self.I_old = self.I
        self.decay = self.decay*np.exp(-1*dt/self.tau)
        V_start = self.V_start + self.decay*self.delta_left + self.decay*self.delta_down*self.R_down
        I_start = (self.V_start + self.decay*self.delta_left)/self.R_off + self.decay*self.delta_down
        V_HL = self.V_HL + self.decay*self.delta_right + self.decay*self.delta_down*self.R_up
        I_HL = V_HL/self.R_off + self.decay*self.delta_down
        I_sat = self.I_sat + self.decay*self.delta_up + self.decay*self.delta_right/self.R_on
        V_sat = V_HL + (I_sat - I_HL) * self.R_up #V_HL + (I_sat - I_HL) * self.R_up
        #I = V/R_on + I_sat - V_sat/R_on = I_start - V_start/R.down + V/R_down
        # V(1/R_on - 1/R_down) = V_sat/R_on - I_sat + I_start - V_start/R_down
        V_LH = (V_sat/self.R_on - I_sat + I_start - V_start/self.R_down)/(1/self.R_on - 1/self.R_down)
        if V <= V_start:
            if (self.state == "down" or self.state == "on") and V <= self.V_start:
                self.decay += 1
                print(V)
                self.state = "off"
            if V >= self.V_start:
                self.I = V/self.R_off + (I_start-V/self.R_off)*(np.exp(self.k*(V-self.V_start))-1)/(np.exp(self.k*(V_start-self.V_start))-1)
            elif V < self.V_start:
                self.I = V/self.R_off
        
        elif self.state == "off" or self.state == "up":
            if V < V_HL:        #check if smaller than V_HL
                self.I = I_start + (V-V_start)/self.R_off #decay after V_start
                self.state = "off"
            else:                       #V is bigger or equal to V_HL
                if V > V_sat: #check if V saturates current
                    self.I = I_sat + (V-V_sat)/self.R_on #propogate the on state
                    self.state = "on" #state is "on"
                else:
                    self.I = V_HL/self.R_off + self.decay*self.delta_down + (V-V_HL)/self.R_up #current up to V_HL and then state change up to V
                    self.state = "up" #state is up but not fully saturated

        elif self.state == "on" or self.state == "down":
            if V > V_LH:
                self.I = I_sat + (V-V_sat)/self.R_on
                self.state = "on"
            else:
                self.I = I_start + (V-V_start)/self.R_down
                self.state = "down"


"""
x = Hysteresis(2, 5, 0.1, 0.1, 1, 1, 15, 1, 1 , 1, 1, 1e32)
Vlist = np.append(np.linspace(0,8,1000), np.linspace(8,0,1000))
Vlist2 = np.append(np.linspace(0,8,1000), np.linspace(8,0,1000))
Vlist3 = np.append(np.linspace(0,10,1000), np.linspace(10,0,1000))
#Vlist = [-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,7,6,5,4,3,2,1,0,-1,-2,-3,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,7,6,5,4,3,2,1,0,-1,-2,-3,-4]
#Vlist = [0, 5, 10, 2, -3, 1, 3]
#Vlist = [0,1,2,3,4,5,6,7,8,7,6,5,4,3,4,5,6,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,-1,-2,-3,-4]
Ilist = []
Ilist2 = []
Ilist3 = []


for V in Vlist:
    x.Update(V, 0.001)
    Ilist.append(x.I)

for V in Vlist2:
    x.Update(V, 0.001)
    Ilist2.append(x.I)

for V in Vlist3:
    x.Update(V, 0.001)
    Ilist3.append(x.I)

plt.close()
plt.plot(Vlist, Ilist)
plt.plot(Vlist2, Ilist2)
plt.plot(Vlist3, Ilist3)
plt.show()
"""




"""
time = np.linspace(0,1.4,1400)
for V in Vlist:
    x.Update(V, 0.001)
    Ilist.append(x.decay)
for V in Vlist:
    x.Update(V, 0.001)
    Ilist.append(x.decay)

for V in Vlist:
    x.Update(V, 0.001)
    Ilist.append(x.decay)

for V in Vlist:
    x.Update(V, 0.001)
    Ilist.append(x.decay)

for V in Vlist:
    x.Update(np.zeros_like(V), 0.001)
    Ilist.append(x.decay)
for V in Vlist:
    x.Update(np.zeros_like(V), 0.001)
    Ilist.append(x.decay)
for V in Vlist:
    x.Update(np.zeros_like(V), 0.001)
    Ilist.append(x.decay)


plt.close()
plt.plot(time, Ilist)
plt.show()
"""

def V_GS(hyst_obj_copy, V_mem, I_DS, R_load)_
    hyst_obj_copy.update(V_mem)
    V_hyst = hyst_obj_copy.V
    V_load = I_DS*R_load
    return V_hyst - V_load


def I_DS(K, V_GS, V_DS, V_th):
    if V_GS < V_th:
        return 0
    elif V_DS > 0 and V_DS < V_GS - V_th:
        return K*((V_GS - V_th)*V_DS - V_DS**2/2)
    elif V_DS > V_GS - V_th and V_GS - V_th > 0:
        return (K/2)*(V_GS - V_th)**2
    
def V_DS(V_mem, I_DS, R_load):
    return V_mem - I_DS*R_load