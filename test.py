"""
Stable simulation of V_mem, V_DS, V_GS using the user's Hysteresis model.

Changes from the previous draft:
- Uses the exact Hysteresis class you provided (unchanged logic).
- Interprets your statement "The current sets the voltage" as a linear mapping
  V_hyst = R_gate * I + V_offset; you can change R_gate or the mapping easily.
- Calls Device.Update(V_mem, dt) exactly once per time step (so the hysteresis
  internal state is advanced only once per dt).
- Uses backward-Euler closed-form for the membrane node for stability.
- Solves V_DS self-consistently while treating the hysteresis current as fixed
  during the inner solve (that matches the physical picture: the hysteresis
  responds to V_mem, not to instantaneous V_DS).

To change the V_GS coupling from current->voltage, edit the R_gate and V_offset
parameters near the top of the simulation.
"""

import numpy as np
import matplotlib.pyplot as plt

# -------------------------
# User-provided Hysteresis class (exactly as given)
# -------------------------
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
        V_sat = V_HL + (I_sat - I_HL) * self.R_up
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

# -------------------------
# Device / circuit functions
# -------------------------

def f_V_DS_solve_quadratic(K, V_GS, V_mem, V_th, R_load):
    """Solve for V_DS in the linear region by solving the quadratic.

    Equation derived from:
        V_mem - V_DS = R_load * I_D
    with linear-region current
        I_D = K*((V_GS - V_th)*V_DS - 0.5 * V_DS**2)
    """
    if V_GS <= V_th:
        return 0.0

    a = -0.5 * R_load * K
    b = R_load * K * (V_GS - V_th) + 1.0
    c = -V_mem

    disc = b * b - 4 * a * c
    if disc < 0:
        return max(0.0, min(V_mem, V_GS - V_th))

    sqrt_disc = np.sqrt(disc)
    v1 = (-b + sqrt_disc) / (2 * a) if a != 0 else -c / b
    v2 = (-b - sqrt_disc) / (2 * a) if a != 0 else -c / b

    candidates = [v for v in (v1, v2) if np.isfinite(v)]
    lower = 0.0
    upper = V_GS - V_th
    valid = [v for v in candidates if lower - 1e-12 <= v <= upper + 1e-12]

    if valid:
        return float(sorted(valid)[0])

    guess = max(lower, min(V_mem, upper))
    return guess


def f_V_DS(K, V_GS, V_mem, V_th, R_load):
    if V_GS <= V_th:
        return V_mem

    I_sat = 0.5 * K * (V_GS - V_th) ** 2
    V_DS_sat_candidate = V_mem - R_load * I_sat

    if V_DS_sat_candidate >= V_GS - V_th:
        return max(0.0, V_DS_sat_candidate)

    V_DS_linear = f_V_DS_solve_quadratic(K, V_GS, V_mem, V_th, R_load)
    V_DS_linear = float(max(0.0, min(V_DS_linear, V_GS - V_th)))
    return V_DS_linear


def Delta_V_mem_backward_euler(dt, V_in, V_mem, C_mem, R_in, R_load, V_DS):
    A = V_in / (C_mem * R_in) + V_DS / (C_mem * R_load)
    B = 1.0 / (C_mem * R_in) + 1.0 / (C_mem * R_load)
    V_new = (V_mem + dt * A) / (1.0 + dt * B)
    return V_new

# -------------------------
# Simulation: parameters, loop
# -------------------------
if __name__ == "__main__":
    # Input waveform
    Vlist = np.append(np.linspace(0, 100, 100), np.linspace(100, 0, 100))
    Vlist = np.append(Vlist, np.zeros_like(Vlist))

    # time step and time vector
    dt = 1.0
    tlist = np.arange(len(Vlist)) * dt

    # circuit/device params (change these to experiment)
    R_in = 1e14
    R_load = 100.0
    C_mem = 1e-12
    Device = Hysteresis(2, 5, 0.1, 0.1, 1, 1, 15, 0, 0 , 0, 0, 1e32)

    V_mem = 0.0
    V_DS = 0.0
    V_GS = 0.0
    V_th = 8.0
    K = 1.0

    # Mapping from hysteresis current to gate voltage: V_hyst = R_gate * I + V_offset
    R_gate = 1.0    # ohms (tune as needed)
    V_offset = 0.0  # volts

    CLAMP_NONNEGATIVE = True

    V_mem_list = [V_mem]
    V_DS_list = [V_DS]
    V_GS_list = [V_GS]

    for idx, V_in in enumerate(Vlist):
        # Update hysteresis once using the membrane voltage (advance its internal state)
        Device.Update(V_mem, dt)

        # initial guesses
        V_DS_guess = V_DS
        V_GS_guess = V_GS

        # Precompute hysteresis-driven gate voltage contribution (treated as fixed during inner solve)
        V_hyst_from_I = R_gate * Device.I + V_offset

        # fixed-point iterate to make V_DS and V_GS consistent
        for it in range(8):
            # V_GS depends on hysteresis (current->voltage) and on V_DS via subtraction term
            V_GS_guess = V_hyst_from_I - (V_mem - V_DS_guess)

            # compute new V_DS based on updated V_GS and current V_mem
            V_DS_new = f_V_DS(K, V_GS_guess, V_mem, V_th, R_load)

            # damping to help convergence
            V_DS_guess = 0.6 * V_DS_guess + 0.4 * V_DS_new

        V_GS = V_hyst_from_I - (V_mem - V_DS_guess)
        V_DS = V_DS_guess

        # update membrane using backward-Euler closed form (stable)
        V_mem = Delta_V_mem_backward_euler(dt, V_in, V_mem, C_mem, R_in, R_load, V_DS)

        if CLAMP_NONNEGATIVE:
            V_mem = max(0.0, V_mem)
            V_DS = max(0.0, V_DS)

        V_mem_list.append(V_mem)
        V_DS_list.append(V_DS)
        V_GS_list.append(V_GS)

    # -------------------------
    # Plot results
    # -------------------------
    plt.figure(figsize=(8, 4))
    plt.plot(tlist, V_mem_list[:-1], label="V_mem")
    plt.plot(tlist, V_DS_list[:-1], label="V_DS")
    plt.plot(tlist, V_GS_list[:-1], label="V_GS")
    plt.xlabel('time (s)')
    plt.ylabel('Voltage (V)')
    plt.legend()
    plt.tight_layout()
    plt.show()
