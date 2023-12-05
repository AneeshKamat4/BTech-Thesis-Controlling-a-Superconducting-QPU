from qm import SimulationConfig
from qm.QuantumMachinesManager import QuantumMachinesManager

from configuration_4qubitsv4 import *
from macros import *
from allxy_macro import *
from helper_functions_RB import *

import numpy as np
from qm.qua import *
from scipy.optimize import curve_fit, minimize
import matplotlib

matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt

qmm = QuantumMachinesManager(qm_ip)

simulate = False
lsb = False
use_drag = False
use_pi12 = False
save_data = False

q_no = 1
qe = f"q{q_no}"
rr = f"rr{q_no}"
# con = f"con{dac_mapping[qe][0]}"
# ro_len = ro_len_clk[str(q_no)]
out = adc_mapping[rr]

show_Rabi = False
init_Rabi = True

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

amp_bounds = [(max(0, 0.5*initial_amps[i]), min(0.9, 2*initial_amps[i])) for i in range(4)]

seed = 345323  # 345324
''' 
Choose appropriate values for ORBIT and Simultaneous RB parameters!
'''

# For ORBIT measurement
RB_length = 3 #00
RB_number = 10
shots = 2 #00

# For ORBIT optimization
# fatol = 1e-4
maxiter = 15
# maxfev = 10

# For Simultaneous RB on all Qubits under consideration
max_circuit_depth = 3 #00
delta_depth = 1  # must be 1!!
num_of_sequences = 10
avgs = 2 #00

if simulate:
	wait_init = 100
	shots = 3
	avgs = 5

with program() as orbit:
	depth = declare(int)
	m = declare(int)
	n = declare(int)

	Ig = declare(fixed)
	Qg = declare(fixed)
	Ie = declare(fixed)
	Qe = declare(fixed)
	I = declare(fixed)
	Q = declare(fixed)

	Ig_st = declare_stream()
	Ie_st = declare_stream()
	I_st = declare_stream()

	with for_(m, 0, m < RB_number, m + 1):
		sequence_list, inv_gate_list = generate_sequence(RB_length, seed)
		assign(sequence_list[RB_length], inv_gate_list[RB_length - 1])

		with for_(n, 0, n < shots, n + 1):
			# Start of Rabi level test
			reset_phase(rr)
			wait(wait_init, qe, rr)

			measure_macro(qe, rr, out, Ig, Qg, pi_12=use_pi12)
			save(Ig, Ig_st)

			reset_phase(rr)
			wait(wait_init, qe, rr)

			play("X180", qe)
			measure_macro(qe, rr, out, Ie, Qe, pi_12=use_pi12)
			save(Ie, Ie_st)
			# End of Rabi level test

			reset_phase(rr)
			wait(wait_init, qe, rr)

			play_sequence(qe, sequence_list, RB_length, drag=use_drag)
			measure_macro(qe, rr, out, I, Q, pi_12=use_pi12)
			save(I, I_st)

	with stream_processing():
		for q_i in q_no_list:
			Ig_st.average().save("Ig")
			Ie_st.average().save("Ie")
			I_st.buffer(shots).map(FUNCTIONS.average()).buffer(RB_number).save("I_avg")


def config_updater(params_array, config, q_no, drag: bool, det: bool):
	'''
	Updates the config with the new parameters

	params_array = [amp_X_pi, amp_Y_pi, amp_X_piby2, amp_Y_piby2]
	'''
	new_amp_X_pi = params_array[0]
	new_amp_Y_pi = params_array[1]
	new_amp_X_piby2 = params_array[2]
	new_amp_Y_piby2 = params_array[3]

	amp_scale[f"{q_no}"]["X180"] = new_amp_X_pi
	amp_scale[f"{q_no}"]["Y180"] = new_amp_Y_pi
	amp_scale[f"{q_no}"]["X90"] = new_amp_X_piby2
	amp_scale[f"{q_no}"]["Y90"] = new_amp_Y_piby2

	rr_no = q_no

	new_config = config_add_elements_q_rr(config, q_no, rr_no, dac_mapping, q_LO, q_IF, q_anh, rr_LO, rr_IF,
	                                      pi_len_ns, piby2_len_ns, pi_rise_grft_ns, amp_scale, drag_dict,
	                                      mixers, mixer_corrections, ro_amp, ro_len_clk, tof, integ_len_clk,
	                                      optimal_readout_phase, smearing=0)

	return new_config


with program() as rb:
	depth = declare(int)
	saved_gate = declare(int)
	m = declare(int)
	n = declare(int)

	I = declare(fixed)
	Q = declare(fixed)

	I_st = declare_stream()
	Q_st = declare_stream()

	with for_(m, 0, m < num_of_sequences, m + 1):
		sequence_list, inv_gate_list = generate_sequence(max_circuit_depth, seed)

		with for_(depth, 1, depth <= max_circuit_depth, depth + delta_depth):
			with for_(n, 0, n < avgs, n + 1):
				assign(saved_gate, sequence_list[depth])
				assign(sequence_list[depth], inv_gate_list[depth - 1])

				reset_phase(rr)
				wait(wait_init, qe, rr)

				play_sequence(qe, sequence_list, depth, drag=use_drag)
				measure_macro(qe, rr, out, I, Q, pi_12=use_pi12)

				save(I, I_st)
				save(Q, Q_st)

				assign(sequence_list[depth], saved_gate)

	with stream_processing():
		I_st.buffer(avgs).map(FUNCTIONS.average()).buffer(num_of_sequences, max_circuit_depth).save("I_avg")
		Q_st.buffer(avgs).map(FUNCTIONS.average()).buffer(num_of_sequences, max_circuit_depth).save("Q_avg")


def Standard_RB(config, optimized=False):
	###########
	# execute #
	###########
	qmm = QuantumMachinesManager(qm_ip)
	qm = qmm.open_qm(config)

	if simulate:
		job = qmm.simulate(config, rb, SimulationConfig(int(60000)))
		# get DAC and digital samples
		samples = job.get_simulated_samples()
		# plot all ports:
		samples.con1.plot()
		samples.con2.plot()
		plt.legend("")
		raise Halted()

	job = qm.execute(rb, duration_limit=0, data_limit=0)

	############
	# analysis #
	############

	res_handles = job.result_handles
	res_handles.wait_for_all_values()

	def power_law(m, a, b, p):
		return a * (p ** m) + b

	avg_trace_values = []
	x = np.linspace(1, max_circuit_depth, max_circuit_depth)

	Ivalue = getattr(res_handles, "I_avg").fetch_all()
	Qvalue = getattr(res_handles, "Q_avg").fetch_all()

	avg_trace = np.average(Ivalue, axis=0)
	avg_trace_values.append(avg_trace)
	init_vals = [-6e-5, 6e-5, 0.99]

	pars, cov = curve_fit(f=power_law, xdata=x, ydata=avg_trace, p0=init_vals,
	                      bounds=(-np.inf, np.inf), maxfev=2000)
	stdevs = np.sqrt(np.diag(cov))
	one_minus_p = 1 - pars[2]
	r_c = one_minus_p * (1 - 1 / 2 ** 1)
	r_g = r_c / 1.875
	r_c_std = stdevs[2] * (1 - 1 / 2 ** 1)
	r_g_std = r_c_std / 1.875

	fid = 1 - r_c

	plt.figure()
	plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

	if optimized:
		plt.title(f"Optimized RB \n Qubit {q_no} Fidelity = {fid:.3%} ", fontsize=14)
	else:
		plt.title(f"Initial RB \n Qubit {q_no} Fidelity = {fid:.3%} ", fontsize=14)

	plt.ylabel("Voltage (a.u.)", fontsize=16)
	plt.xlabel("No. of Cliffords", fontsize=16)
	plt.plot(avg_trace, ".r", markersize=6, alpha=0.7)  # plot averaged trace

	for j in range(Ivalue.shape[0]):
		plt.plot(Ivalue[j], '.', alpha=0.4, markersize=3)  # plot individual traces in 4k colour

	plt.plot(x, power_law(x, *pars), '-r')
	plt.grid()
	plt.show()

	print(f'~~~~~~~~~~~~~~ FOR QUBIT {q_i} ~~~~~~~~~~~~~~~~~')
	print('#########################')
	print('### Fitted Parameters ###')
	print('#########################')
	print(
		f'A = {pars[0]:.3} ({stdevs[0]:.1}), B = {pars[1]:.3} ({stdevs[1]:.1}), p = {pars[2]:.3} ({stdevs[2]:.1})')
	print('Covariance Matrix')
	print(cov)

	print('#########################')
	print('### Useful Parameters ###')
	print('#########################')
	print(f'1-p = {np.format_float_scientific(one_minus_p, precision=2)} ({stdevs[2]:.1}), '
	      f'r_c = {np.format_float_scientific(r_c, precision=2)} ({r_c_std:.1}), '
	      f'r_g = {np.format_float_scientific(r_g, precision=2)}  ({r_g_std:.1})'
	      f'fid = {fid:.3%}')

	return f"{fid:.3%}"


def ORBIT(params_list, config, orbit, qmm, shape, drag: bool, det: bool):
	params_array = params_list.reshape(shape)
	mod_config = config_updater(params_array, config, q_no_list, drag, det)

	###########
	# execute #
	###########

	qm = qmm.open_qm(mod_config)
	job = qm.execute(orbit, duration_limit=0, data_limit=0)

	############
	# analysis #
	############

	res_handles = job.result_handles
	res_handles.wait_for_all_values()

	goal_list = []

	Ig = getattr(res_handles, "Ig").fetch_all()
	Ie = getattr(res_handles, "Ie").fetch_all()
	I_avg = getattr(res_handles, "I_avg").fetch_all()

	pop = np.divide(I_avg - Ig, Ie - Ig)
	goal = np.mean(pop)

	goal_list.append(goal)
	goal_avg = np.mean(goal_list)

	return goal_avg


# ORBIT Optimizer
''' 
   	params_array = [amp_X_pi, amp_Y_pi, amp_X_piby2, amp_Y_piby2]

'''
drag = use_drag
det = False

initial_amps, amp_bounds


print(f'The initial parameters are {initial_amps}')

options = {'disp': True, 'return_all': True, 'maxiter': maxiter} #'fatol': fatol,

# Formatting the initial array for config updates
init_config = config_updater(initial_amps, config, q_no, drag, det)
init_fid = Standard_RB(init_config)
print(f"The initial fidelities for the qubits under consideration with RB are: \n{init_fid_qpu}")

result = minimize(fun=ORBIT, x0=init_params, args=(config, orbit, qmm, shape, drag, det), bounds=bounds_params,
                  method='Nelder-Mead', tol=1e-1, options=options)

all_params_list = result.allvecs
iter_array = []

for i in range(len(all_params_list)):
	iter_array.append(i)

	one_param_list = np.array(all_params_list[i])
	one_param_array = one_param_list.reshape(shape)

	par_list = ["amp_X_pi", "amp_Y_pi", "amp_X_piby2", "amp_Y_piby2"]

# # Comparing one parameter for all qubits
for par in par_list:
	plt.figure()
	for q_i in q_no_list:
		plt.plot(iter_array, eval(f"p_{par}_{q_i}"), '.', markersize=5, linestyle='-', label=f'Q{q_i} {par} amp')

	plt.xlabel('Iterations')
	plt.ylabel('Parameter Value')
	plt.legend()
	plt.tight_layout()
	plt.show()

# Comparing all parameters for one qubit
for q_i in q_no_list:
	plt.figure()
	for par in par_list:
		plt.plot(iter_array, eval(f"p_{par}_{q_i}"), '.', markersize=5, linestyle='-', label=f'Q{q_i} {par} amp')

	plt.xlabel('Iterations')
	plt.ylabel('Parameter Value')
	plt.legend()
	plt.tight_layout()
	plt.show()

opt_params_list = np.array(result.x)

# Formatting the optimized array for config updates
opt_params_array = opt_params_list.reshape(shape)
opt_config = config_updater(opt_params_array, config, q_no_list, drag, det)
opt_fid_list = QPU_Simultaneous_RB(opt_config, True)
opt_fid_qpu = {}
for q_i in q_no_list:
	opt_fid_qpu[f"Qubit {q_i}"] = opt_fid_list[q_index[f"{q_i}"]]
print(f"The optimal fidelities for the qubits under consideration with RB are: \n{opt_fid_qpu}")
