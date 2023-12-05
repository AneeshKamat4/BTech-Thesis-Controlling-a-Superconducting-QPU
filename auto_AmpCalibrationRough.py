from allxy_macro import *
from qm.QuantumMachinesManager import QuantumMachinesManager
from scipy.optimize import minimize
from qualang_tools.units import unit
u = unit()

qop_ip = '192.168.0.106'
qmm = QuantumMachinesManager(qm_ip)

show_Rabi = True
print_Rabi = False
simulate = False
init_Rabi = True

q_no = 1
qe = f"q{q_no}"
rr = f"rr{q_no}"
# con = f"con{dac_mapping[qe][0]}"
# ro_len = ro_len_clk[str(q_no)]
out = adc_mapping[rr]

n_avg = 1e2
wait_init = 250000 # in clocks ~ 1 ms

t_min_ns = 16
t_max_ns = 1600
dt_ns = 8 # minimum 4ns
rabi_pi_amp = 0.5

t_min = int(t_min_ns/4) # in clocks
t_max = int(t_max_ns/4)
dt = int(dt_ns/4)
t_list_clk = np.arange(t_min, t_max, dt)

ensemble_params = (n_avg, t_min, t_max, dt, t_list_clk, wait_init, qe, rr, out, rabi_pi_amp)

initial_amps = initialize_Rabi(q_no, qmm, *ensemble_params, show_Rabi)

if init_Rabi is False:
    initial_amps = np.array([0.5]*2 + [0.25]*2)

print(f"The roughly calibrated amplitudes are: {initial_amps}.")

allxy_fid = opt_allxy(initial_amps, q_no, qmm, *ensemble_params)

print(f"The ALLXY fidelity is: {allxy_fid}.")


amp_bounds = [(max(0, 0.5*initial_amps[i]), min(0.9, 2*initial_amps[i])) for i in range(4)]

# Optimize over the amplitude for the single-qubit gate pulses
res = minimize(opt_allxy, x0=initial_amps, args=(q_no, qmm, *ensemble_params), method='L-BFGS-B',
               bounds=amp_bounds, options={"maxfun": 50, "disp": True, "ftol": 5e-2})
opt_amplitude = res.x
print("The optimal amplitudes for the single-qubit gateset are \n"
       f"X180 - {opt_amplitude[0]},\n Y180 - {opt_amplitude[1]},\n"
       f"X90 - {opt_amplitude[2]},\n Y90 - {opt_amplitude[3]}].")