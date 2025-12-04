import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import copy


class Hysteresis:
    def __init__(self, V_start, V_HL, R_down, R_up, R_off, R_on, I_sat, delta_down, delta_up, delta_right, delta_left, tau, k = 10):
        # Store parameters
        self.V_start = V_start #where the left side meets the off state
        self.V_HL = V_HL #Default threshold voltage
        self.R_up = R_up #resistance on right side
        self.R_off = R_off #resistance on bottom side
        self.R_on = R_on #resistance on top side
        self.R_down = R_down #resistance on left side
        self.V = 0.0 #Voltage at this time step
        self.I = 0.0 #current at this time step
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
        self.current_V_start = self.V_start + self.decay*self.delta_left + self.decay*self.delta_down*self.R_down
        I_start = (self.V_start + self.decay*self.delta_left)/self.R_off + self.decay*self.delta_down
        V_HL = self.V_HL + self.decay*self.delta_right + self.decay*self.delta_down*self.R_up
        I_HL = V_HL/self.R_off + self.decay*self.delta_down
        I_sat = self.I_sat + self.decay*self.delta_up + self.decay*self.delta_right/self.R_on
        V_sat = V_HL + (I_sat - I_HL) * self.R_up #V_HL + (I_sat - I_HL) * self.R_up
        #I = V/R_on + I_sat - V_sat/R_on = I_start - V_start/R.down + V/R_down
        # V(1/R_on - 1/R_down) = V_sat/R_on - I_sat + I_start - V_start/R_down
        V_LH = (V_sat/self.R_on - I_sat + I_start - self.current_V_start/self.R_down)/(1/self.R_on - 1/self.R_down)
        if V <= self.current_V_start:
            if (self.state == "down" or self.state == "on"):
                self.decay += 1
                self.state = "off"
            if V >= self.V_start:
            #    self.I = V/self.R_off + (I_start-V/self.R_off)*(np.exp(self.k*(V-self.V_start))-1)/(np.exp(self.k*(self.current_V_start-self.V_start))-1)
                self.I = self.I = V/self.R_off
            elif V < self.V_start:
                self.I = V/self.R_off
        
        elif self.state == "off" or self.state == "up":
            if V < V_HL:        #check if smaller than V_HL
                self.I = I_start + (V-self.current_V_start)/self.R_off #decay after V_start
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
                self.I = I_start + (V-self.current_V_start)/self.R_down
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

"""
Vlist = np.append(np.linspace(0,10,1600), np.linspace(10,0,1600))
tlist = np.linspace(0,10,6400)
Device = Hysteresis(0.1, 7, 1e-6, 1e-6, 100, 100, 1, 1, 1 , 0, 1, 1e32)
dt = tlist[1] - tlist[0]

Ilist = []
for n in range(3):
    Ilist = []
    for i in range(len(Vlist)):
        V = Vlist[i]
        Device.Update(V, dt)
        Ilist.append(Device.I)
    plt.plot(Vlist,Ilist)

plt.xlabel("Voltage [V]")
plt.ylabel("Current [A]")
plt.show()
"""


def f_V_GS(hyst_obj_copy, V_mem, V_DS, dt):
    hyst_obj_copy.Update(V_mem, dt)
    V_hyst = hyst_obj_copy.I
    if np.isnan(V_hyst):
        V_hyst = 0.0
    V_hyst = max(V_hyst, 0)

    D_V = (V_mem - V_DS)
    if np.isnan(D_V):
        D_V = 0.0
    #if D_V < -1e-3:
    #    print(V_hyst, D_V)
    D_V = max(D_V, 0)
    return V_hyst - D_V


def f_V_DS(K, V_GS, V_mem, V_th, R_S, V_DS):
    if V_GS < V_th:
        return V_DS - V_mem
    elif V_DS > 0 and V_DS < V_GS - V_th:
        expr = V_mem - R_S*K*((V_GS - V_th)*V_DS - V_DS**2/2)
        return V_DS - expr
    elif V_DS > V_GS - V_th and V_GS - V_th > 0:
        expr = V_mem - R_S*(K/2)*(V_GS - V_th)**2
        return V_DS - expr
    return V_DS - V_mem

def Delta_V_mem(dt, V_in, V_mem, C_mem, R_in, R_S, V_DS, R_L):
        RC_Voltage = max((V_in - V_mem)/(C_mem*R_in), 0)
        return dt*(RC_Voltage - (V_mem-V_DS)/(C_mem*R_S) - V_mem/(C_mem*R_L))

def equations(x, hyst_obj, V_in, V_mem_old, V_th, R_S, R_in, C_mem, K, dt, R_L):
    V_GS, V_DS, V_mem_new = x

    # 1. Hysteresis relation
    hyst_obj_copy = copy.deepcopy(hyst_obj)   # Make a temporary copy
    f1 = f_V_GS(hyst_obj_copy, V_mem_old, V_DS, dt) - V_GS

    # 2. MOSFET DS equation
    f2 = f_V_DS(K, V_GS, V_mem_old, V_th, R_S, V_DS)

    # 3. Membrane capacitor update
    f3 = V_mem_new - (V_mem_old + Delta_V_mem(dt, V_in, V_mem_old, C_mem, R_in, R_S, V_DS, R_L))

    return [f1, f2, f3]

#E_spike_list = [] #left overs from running for loop
#P_list = []


Vlist = np.ones(81920)*1
Vlist = np.append(Vlist, Vlist)
tlist = np.linspace(0,0.0008,163840)
R_in = 7e6 #3.3e6
R_S = 1e2 #1e2
R_L = 1e32
R_on = 1e3
R_off = 1e3
V_LH = 0.1
V_HL = 0.7
C_mem = 1e-11 #works well for 1.5e-14
Device = Hysteresis(V_LH, V_HL, 1e-32, 1e-32, R_off, R_on, 1, 0 ,0 , 0, 0, 0.2)
V_mem = 0
V_DS = 0
V_GS = 0
V_th = 0.7
K = 2e-5
dt = tlist[1] - tlist[0]

V_mem_list = []
V_DS_list = []
V_GS_list = []
V_out_list = []

for i in range(len(tlist)):
    V_in = Vlist[i]

    # initial guess
    x0 = [V_GS, V_DS, V_mem]

    sol = fsolve(
        equations, 
        x0,
        args=(Device, V_in, V_mem, V_th, R_S, R_in, C_mem, K, dt, R_L)
    )

    V_GS, V_DS, V_mem = sol
    V_mem = max(V_mem, 0)
    #V_DS = max(V_DS, 0)
    #V_GS = max(V_GS, 0)
    Device.Update(V_mem, dt)
    V_out = Device.I
    V_out = max(V_out, 0)
    
    V_GS_list.append(V_GS)
    V_DS_list.append(V_DS)
    V_mem_list.append(V_mem)
    V_out_list.append(V_out)
    #print("V_GS - V_th: ", V_GS - V_th)
    #print("V_DS: ", V_DS)
    t = tlist[i]
    #if t > 0.1339 and t < 0.1374:
    #    print("V_GS, V_G, V_mem - V_DS", V_GS, V_out, V_mem - V_DS)

on = "False"
t_start = 0 #start of spike
t_end = 0 #end of spike
spike_interval = [] #list of time between spikes
spike_width = [] #list of time of spikes
a_list = [] #list of voltage slopes of spikes
V_end = 0 #Voltage at the end of spike
V_start = 0 #Voltage at the beginning of spike

for i in range(len(V_out_list)):
    if i == 0 or i == 1: #No spikes will occur at the first two time steps so this is okay
        pass
    elif V_out_list[i-1] > 0.5: #When the output voltage increases over 0 (0.5 for saftey) we have a spike
        if V_out_list[i] < 0.5: # If the next timestep have a lower voltage we are at the end of the spike
            t_end = tlist[i-1]
            V_end = V_out_list[i-1]
            spike_width.append(t_end - t_start)
            a_list.append((V_end-V_start)/(t_end - t_start)) #the formula for the slope of V(t) in the spike
        elif V_out_list[i-2] <= 0.5: #if i-1 is above 0.5V and the step before that is below --> i-1 is at the beginning of the spike
            t_start = tlist[i-1]
            V_start = V_out_list[i-1]
            spike_interval.append(t_start - t_end) #spike interval is the difference between the start of new spike and the end of previous spike



plt.figure()
tlist = tlist*1e9
plt.plot(tlist,Vlist, label = "V_in", color = "k")
plt.plot(tlist, V_mem_list, label = "V_mem", color = "r")
#plt.plot(tlist, V_DS_list, label = "V_DS")
#plt.plot(tlist, V_GS_list, label = "V_GS")
plt.plot(tlist, V_out_list, label = "V_out", color = "b")
plt.xlabel(u'Time [\xb5s]')
plt.ylabel("Voltage [V]")
plt.legend()
plt.show()


a = a_list[-1]
V0 = V_end
T = spike_width[-1]

Frequency = 1/(spike_width[-1]+spike_interval[-1])
E_per_spike = 1/R_on*(V0**2*T + a*V0*T**2 + a**2*T**3/3)
Power_consumption = E_per_spike*Frequency

print("V_th = ", V_th, " completed")
print(u"spike width [\xb5s]: ", spike_width[-1]*1e6)
print(u"spike distance [\xb5s]: ", spike_interval[-1]*1e6)
print("Frequency kHz", Frequency/1e3)
print("Energy per spike [nJ]: ", E_per_spike*1e9)
print(u"Power consumption [\xb5W]: ", Power_consumption*1e6)
print("time step [ns]: ", dt*1e9)


