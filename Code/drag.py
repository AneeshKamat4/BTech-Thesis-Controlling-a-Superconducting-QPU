"""
An experiment to calibrate the DRAG coefficient: drag_coef
This protocol is described in Reed's thesis (Fig. 5.8) https://rsl.yale.edu/sites/default/files/files/RSL_Theses/reed.pdf
This protocol was also cited in: https://doi.org/10.1103/PRXQuantum.2.040202
"""
from qm.qua import *
from qm.QuantumMachinesManager import QuantumMachinesManager
from configuration_4qubitsv4 import *
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import numpy as np
from qm import SimulationConfig
from qualang_tools.loops import from_array
from scipy.optimize import curve_fit

from qualang_tools.units import unit
u = unit()

q_no = 1
qe = f"q{q_no}"
rr = f"rr{q_no}"
con = f"con{dac_mapping[qe][0]}"
ro_len = ro_len_clk[str(q_no)]
out = adc_mapping[rr]
drag_coef = drag_dict[f"{q_no}"]["alpha"]

# ----------------------------
# Time-Rabi to find the levels
# ----------------------------

t_min_ns = 0
t_max_ns = 1200
dt_ns = 8 # minimum 4ns
n_avg_rabi = 300
rabi_pi_amp = 0.5

t_min = int(t_min_ns/4) # in clocks
t_max = int(t_max_ns/4)
dt = int(dt_ns/4)
t_list = np.arange(t_min, t_max, dt)

rep_rate_clk = 250000
wait_rr = 16

with program() as rabi:
    n = declare(int)
    I = declare(fixed)
    I_st = declare_stream()
    Q = declare(fixed)
    Q_st = declare_stream()
    t = declare(int)

    with for_(n, 0, n < n_avg_rabi, n + 1):
        with for_(t, t_min, t < t_max + 0.1, t + dt):
            wait(rep_rate_clk - t - wait_rr - ro_len)
            play("grft"*amp(rabi_pi_amp), qe, t)
            # play("X180"*amp(0.8), qe, t)
            align(qe, rr)
            wait(wait_rr, rr)
            measure("readout", rr, None,
                    demod.full("integW_cos", I, out),
                    demod.full("integW_minus_sin", Q, out))
            save(I, I_st)
            save(Q, Q_st)

    with stream_processing():
        I_st.buffer(len(t_list)).average().save('I')
        Q_st.buffer(len(t_list)).average().save('Q')


qmm = QuantumMachinesManager(qm_ip)

# ------------------------
#       Execute on the OPX
# ------------------------
qm1 = qmm.open_qm(config)
job = qm1.execute(rabi)
res_handles = job.result_handles
I_handle = job.result_handles.get("I")
Q_handle = job.result_handles.get("Q")

t_list = 4 * t_list
plt.figure()
plt.title("Rabi")
I_handle.wait_for_values(1)
Q_handle.wait_for_values(1)

while res_handles.is_processing():
    I = I_handle.fetch_all()
    Q = Q_handle.fetch_all()
    sig = I + 1j * Q

    plt.clf()
    plt.plot(t_list, Q, marker='.', label="Q")
    plt.plot(t_list, I, marker='.', label="I")
    plt.xlabel("Time (ns)")
    plt.ylabel("Rabi Amplitude")
    plt.title(f"Time Rabi on Qubit {q_no}")
    plt.legend()
    plt.grid()
    # plt.ylim((-0.00015, 0))
    plt.pause(0.25)

I = job.result_handles.get("I").fetch_all()
Q = job.result_handles.get("Q").fetch_all()

qm1.close()

# ############
# # analysis #
# ############
sig = I + 1j*Q

def rabi_fit(t, A, f, d, p, c):

    return A * np.exp(-t/d) * np.sin(2*np.pi*f*t + p) + c

pars, cov = curve_fit(f=rabi_fit, xdata=t_list, ydata=I, p0=[3e-3,0.01,100,0, 1e-5], bounds=(-np.inf, np.inf), maxfev=2000)
# init_pars = p0 = [A, f, d, p, c]

print('######################### \n### Fitted Parameters ### \n######################### ')
print(f"Rabi frequency = {np.round(1e3*pars[1],2)} MHz")
print(f"Pi pulse = {np.round(0.5/pars[1],3)} ns")
print(f"Rabi amplitude = {pars[0]}")
print(f"Rabi decay constant = {pars[2]*1e-3} us")

I_eq = pars[4]
rabi_amp = abs(pars[0])
I_g = I_eq - rabi_amp
I_e = I_eq + rabi_amp

print(f"Rabi DC offset = {I_eq}")
print(f"In arbitrary units \n Ground state [-1] = {I_g} \n Equator [0] = {I_eq} \n Excited state [+1] = {I_e}")

###################
# The QUA program #
###################

# set the drag_coef in the configuration
# drag_coef = 1

n_avg = 1e4

best_qubit_T1 = 40 * u.us
cooldown_time = int(5 * best_qubit_T1 // 4) # in clock cycles

a_min = -1.0
a_max = 1.0
da = 0.01
amps = np.arange(a_min, a_max + da / 2, da)  # + da/2 to add a_max to amplitudes

with program() as drag:
    n = declare(int)
    n_st = declare_stream()
    a = declare(fixed)
    I = declare(fixed)
    Q = declare(fixed)
    I1_st = declare_stream()
    Q1_st = declare_stream()
    I2_st = declare_stream()
    Q2_st = declare_stream()

    with for_(n, 0, n < n_avg, n + 1):
        with for_(*from_array(a, amps)):
            play("d_X180" * amp(1, 0, 0, a), qe)
            play("d_Y90" * amp(a, 0, 0, 1), qe)

            align(qe, rr)
            measure("readout", rr, None,
                demod.full("integW_cos", I, out),
                demod.full("integW_minus_sin", Q, out))
            
            save(I, I1_st)
            save(Q, Q1_st)
            wait(cooldown_time, rr)
            align()

            play("d_Y180" * amp(a, 0, 0, 1), qe)
            play("d_X90" * amp(1, 0, 0, a), qe)

            align(qe, rr)
            measure("readout", rr, None,
                demod.full("integW_cos", I, out),
                demod.full("integW_minus_sin", Q, out))
            
            save(I, I2_st)
            save(Q, Q2_st)
            wait(cooldown_time, rr)
        save(n, n_st)
    
    with stream_processing():
        I1_st.buffer(len(amps)).average().save("I1")
        Q1_st.buffer(len(amps)).average().save("Q1")
        I2_st.buffer(len(amps)).average().save("I2")
        Q2_st.buffer(len(amps)).average().save("Q2")
        n_st.save("Iteration")

#####################################
#  Open Communication with the QOP  #
#####################################
qmm = QuantumMachinesManager(qm_ip)

simulate = False

if simulate:
    simulation_config = SimulationConfig(duration=1000)  # in clock cycles
    job = qmm.simulate(config, drag, simulation_config)
    samples = job.get_simulated_samples()
    # plot all ports:
    sim_output = getattr(samples, con)
    sim_output.plot()
    plt.legend("")
    plt.show()
    raise(Halted)

qm2 = qmm.open_qm(config)

job = qm2.execute(drag)
# Get results from QUA program
results = fetching_tool(job, data_list=["I1", "I2", "Iteration"], mode="live")

# Live plotting
fig = plt.figure()
interrupt_on_close(fig, job)  # Interrupts the job when closing the figure

while results.is_processing():
    # Fetch results
    I1, I2, iter = results.fetch_all()

    # Rescale results
    I1 = np.array(I1)
    I2 = np.array(I2)
    I_pair = np.array([I1, I2])
    z_pair = (I_pair - I_eq) / rabi_amp
    state_pair = (z_pair + 1) / 2

    # Progress bar
    progress_counter(iter, n_avg, start_time=results.get_start_time())
    # Plot results
    plt.cla()
    plt.plot(amps * drag_coef, state_pair[0], label="x180y90")
    plt.plot(amps * drag_coef, state_pair[1], label="y180x90")
    plt.xlabel("DRAG coefficient")
    plt.ylabel("State probability")
    plt.legend()
    plt.tight_layout()
    plt.pause(0.5)
    plt.show()

# Fetch results
I1, I2, iter = results.fetch_all()

# Rescale results
I1 = np.array(I1)
I2 = np.array(I2)
I_pair = np.array([I1, I2])
z_pair = (I_pair - I_eq) / rabi_amp
state_pair = (z_pair + 1) / 2

qm2.close()

# Analysis

def linear_fit(x, m, c):
    return m * x + c

pars1, cov1 = curve_fit(f=linear_fit, xdata=amps, ydata=state_pair[0], p0=[-0.4, 0.5], bounds=(-np.inf, np.inf), maxfev=2000)
pars2, cov2 = curve_fit(f=linear_fit, xdata=amps, ydata=state_pair[1], p0=[0.4, 0.5], bounds=(-np.inf, np.inf), maxfev=2000)

opt_drag_coef = (pars1[1] - pars2[1])/(pars2[0] - pars1[0])

print(f"The optimal DRAG coefficient is {opt_drag_coef}")

plt.figure()
plt.plot(amps * drag_coef, state_pair[0], label="x180y90")
plt.plot(amps * drag_coef, state_pair[1], label="y180x90")
plt.plot(amps * drag_coef, linear_fit(amps, *pars1), label="x180y90")
plt.plot(amps * drag_coef, linear_fit(amps, *pars2), label="y180x90")
plt.axvline(x=opt_drag_coef)
plt.xlabel("DRAG Coefficient")
plt.ylabel("State probability")
plt.title("Yale method to find the DRAG Coefficient")
plt.legend()
plt.tight_layout()
plt.show()
