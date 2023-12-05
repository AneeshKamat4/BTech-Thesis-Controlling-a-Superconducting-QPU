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

n_qubits = 4
wait_init = 250000

q_no_list = [i + 1 for i in range(n_qubits)]
q_no_list = [1]

# [q_i-1] -> q_index[f"{q_i}"]
q_index = {}
qe_list, rr_list, out_list = [], [], []
# qe, rr, out,
for q_i in q_no_list:
	q_index[f"{q_i}"] = q_no_list.index(q_i)
	qe_list.append(f"q{q_i}")
	rr = f"rr{q_i}"
	rr_list.append(rr)
	out_list.append(adc_mapping[rr])

seed = 345323  # 345324


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

	Ig = declare(fixed, size=len(q_no_list))
	Qg = declare(fixed, size=len(q_no_list))
	Ie = declare(fixed, size=len(q_no_list))
	Qe = declare(fixed, size=len(q_no_list))
	I = declare(fixed, size=len(q_no_list))
	Q = declare(fixed, size=len(q_no_list))

	for q_i in q_no_list:
		vars()[f"Ig_{q_i}_st"] = declare_stream()
		vars()[f"Ie_{q_i}_st"] = declare_stream()
		vars()[f"I_{q_i}_st"] = declare_stream()

	with for_(m, 0, m < RB_number, m + 1):
		sequence_list, inv_gate_list = generate_sequence(RB_length, seed)
		assign(sequence_list[RB_length], inv_gate_list[RB_length - 1])

		with for_(n, 0, n < shots, n + 1):
			for q_i in q_no_list:  # Python loop
				# Start of Rabi level test
				reset_phase(rr_list[q_index[f"{q_i}"]])
				wait(wait_init, qe_list[q_index[f"{q_i}"]], rr_list[q_index[f"{q_i}"]])

				measure_macro(qe_list[q_index[f"{q_i}"]], rr_list[q_index[f"{q_i}"]], out_list[q_index[f"{q_i}"]],
				              Ig[q_index[f"{q_i}"]], Qg[q_index[f"{q_i}"]], pi_12=use_pi12)
				save(Ig[q_index[f"{q_i}"]], eval(f'Ig_{q_i}_st'))

				reset_phase(rr_list[q_index[f"{q_i}"]])
				wait(wait_init, qe_list[q_index[f"{q_i}"]], rr_list[q_index[f"{q_i}"]])

				play("X180", qe_list[q_index[f"{q_i}"]])
				measure_macro(qe_list[q_index[f"{q_i}"]], rr_list[q_index[f"{q_i}"]], out_list[q_index[f"{q_i}"]],
				              Ie[q_index[f"{q_i}"]], Qe[q_index[f"{q_i}"]], pi_12=use_pi12)
				save(Ie[q_index[f"{q_i}"]], eval(f'Ie_{q_i}_st'))
				# End of Rabi level test

				reset_phase(rr_list[q_index[f"{q_i}"]])
				wait(wait_init, qe_list[q_index[f"{q_i}"]], rr_list[q_index[f"{q_i}"]])

				play_sequence(qe_list[q_index[f"{q_i}"]], sequence_list, RB_length, drag=use_drag)
				measure_macro(qe_list[q_index[f"{q_i}"]], rr_list[q_index[f"{q_i}"]], out_list[q_index[f"{q_i}"]],
				              I[q_index[f"{q_i}"]], Q[q_index[f"{q_i}"]], pi_12=use_pi12)
				save(I[q_index[f"{q_i}"]], eval(f"I_{q_i}_st"))

	with stream_processing():
		for q_i in q_no_list:
			eval(f"Ig_{q_i}_st").average().save(f"Ig_{q_i}")
			eval(f"Ie_{q_i}_st").average().save(f"Ie_{q_i}")
			eval(f"I_{q_i}_st").buffer(shots).map(FUNCTIONS.average()).buffer(RB_number).save(f"I_{q_i}_avg")


def config_updater(params_array, config, q_no_list, drag: bool, det: bool):
	'''
	Updates the config with the new parameters

	params_{q_i} = [amp_pi_{q_i}, amp_piby2_{q_i}, alpha_{q_i}, det_{q_i}]
	params_array = [params1, params2, ...]

	params_1 = [amp_pi_1, amp_piby2_1, 0, 0]
	'''
	for q_i in q_no_list:  # Python loop
		vars()[f"new_amp_pi_{q_i}"] = params_array[q_index[f"{q_i}"]][0]
		vars()[f"new_amp_piby2_{q_i}"] = params_array[q_index[f"{q_i}"]][1]

		amp_scale[f"{q_i}"]["X180"] = eval(f"new_amp_pi_{q_i}")
		amp_scale[f"{q_i}"]["Y180"] = eval(f"new_amp_pi_{q_i}")
		amp_scale[f"{q_i}"]["X90"] = eval(f"new_amp_piby2_{q_i}")
		amp_scale[f"{q_i}"]["Y90"] = eval(f"new_amp_piby2_{q_i}")

		if drag:
			vars()[f"new_alpha_{q_i}"] = params_array[q_index[f"{q_i}"]][2]
			drag_dict[f"{q_i}"]["alpha"] = eval(f"new_alpha_{q_i}")
			if det:
				vars()[f"new_det_{q_i}"] = params_array[q_index[f"{q_i}"]][3]
				drag_dict[f"{q_i}"]["det"] = eval(f"new_det_{q_i}") * u.MHz

		rr_i = q_i

		new_config = config_add_elements_q_rr(config, q_i, rr_i, dac_mapping, q_LO, q_IF, q_anh, rr_LO, rr_IF,
		                                      pi_len_ns, piby2_len_ns, pi_rise_grft_ns, amp_scale, drag_dict,
		                                      mixers, mixer_corrections, ro_amp, ro_len_clk, tof, integ_len_clk,
		                                      optimal_readout_phase, smearing=0)

	return new_config


with program() as rb:
	depth = declare(int)
	saved_gate = declare(int)
	m = declare(int)
	n = declare(int)

	I = declare(fixed, size=len(q_no_list))
	Q = declare(fixed, size=len(q_no_list))

	for q_i in q_no_list:
		vars()[f"I{q_i}_st"] = declare_stream()
		vars()[f"Q{q_i}_st"] = declare_stream()

	with for_(m, 0, m < num_of_sequences, m + 1):
		sequence_list, inv_gate_list = generate_sequence(max_circuit_depth, seed)

		with for_(depth, 1, depth <= max_circuit_depth, depth + delta_depth):
			with for_(n, 0, n < avgs, n + 1):
				assign(saved_gate, sequence_list[depth])
				assign(sequence_list[depth], inv_gate_list[depth - 1])
				for q_i in q_no_list:
					reset_phase(rr_list[q_index[f"{q_i}"]])
					wait(wait_init, qe_list[q_index[f"{q_i}"]], rr_list[q_index[f"{q_i}"]])

					play_sequence(qe_list[q_index[f"{q_i}"]], sequence_list, depth, drag=use_drag)
					measure_macro(qe_list[q_index[f"{q_i}"]], rr_list[q_index[f"{q_i}"]], out_list[q_index[f"{q_i}"]],
					              I[q_index[f"{q_i}"]], Q[q_index[f"{q_i}"]], pi_12=use_pi12)

					save(I[q_index[f"{q_i}"]], eval(f"I{q_i}_st"))
					save(Q[q_index[f"{q_i}"]], eval(f"Q{q_i}_st"))

				assign(sequence_list[depth], saved_gate)

	with stream_processing():

		for q_i in q_no_list:
			eval(f"I{q_i}_st").buffer(avgs).map(FUNCTIONS.average()).buffer(num_of_sequences, max_circuit_depth).save(
				f"I{q_i}_avg")
			eval(f"Q{q_i}_st").buffer(avgs).map(FUNCTIONS.average()).buffer(num_of_sequences, max_circuit_depth).save(
				f"Q{q_i}_avg")


def QPU_Simultaneous_RB(config, optimized=False):
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
	fid_list = []

	for q_i in q_no_list:

		vars()[f"I{q_i}value"] = getattr(res_handles, f"I{q_i}_avg").fetch_all()
		vars()[f"Q{q_i}value"] = getattr(res_handles, f"Q{q_i}_avg").fetch_all()

		avg_trace = np.average(eval(f"I{q_i}value"), axis=0)
		avg_trace_values.append(avg_trace)
		init_vals = [-6e-5, 6e-5, 0.99]

		pars, cov = curve_fit(f=power_law, xdata=x, ydata=avg_trace, p0=init_vals, bounds=(-np.inf, np.inf),
		                      maxfev=2000)
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
			plt.title(f"Optimized RB \n Qubit {q_i} Fidelity = {fid:.3%} ", fontsize=14)
		else:
			plt.title(f"Initial RB \n Qubit {q_i} Fidelity = {fid:.3%} ", fontsize=14)

		plt.ylabel("Voltage (a.u.)", fontsize=16)
		plt.xlabel("No. of Cliffords", fontsize=16)
		plt.plot(avg_trace, ".r", markersize=6, alpha=0.7)  # plot averaged trace

		for j in range(eval(f"I{q_i}value").shape[0]):
			plt.plot(eval(f"I{q_i}value")[j], '.', alpha=0.4, markersize=3)  # plot individual traces in 4k colour

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
		fid_list.append(f"{fid:.3%}")

	return fid_list


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
	for q_i in q_no_list:
		vars()[f"Ig_{q_i}"] = getattr(res_handles, f"Ig_{q_i}").fetch_all()
		vars()[f"Ie_{q_i}"] = getattr(res_handles, f"Ie_{q_i}").fetch_all()
		vars()[f"I_{q_i}_avg"] = getattr(res_handles, f"I_{q_i}_avg").fetch_all()

		vars()[f"pop_{q_i}"] = np.divide(eval(f"I_{q_i}_avg") - eval(f"Ig_{q_i}"),
		                                 eval(f"Ie_{q_i}") - eval(f"Ig_{q_i}"))
		vars()[f"goal_{q_i}"] = np.mean(eval(f"pop_{q_i}"))

		goal_list.append(eval(f"goal_{q_i}"))
	goal_avg = np.mean(goal_list)

	return goal_avg


# ORBIT Optimizer
''' 
    params_{q_i} = [amp_pi_{q_i}, amp_piby2_{q_i}, alpha_{q_i}, det_{q_i}]
    params_list = [params1, params2, ...]
'''
drag = use_drag
det = False

params_list = []

ensemble_params = (n_avg, t_min, t_max, dt, t_list_clk, wait_init, qe_list[0], rr_list[0], out_list[0], rabi_pi_amp)
show_Rabi = True
initial_amps = initialize_Rabi(q_no, qmm, *ensemble_params, show_Rabi)

for q_i in q_no_list:
	# vars()[f"init_amp_pi_{q_i}"] = 0.5
	# vars()[f"init_amp_piby2_{q_i}"] = 0.25

	vars()[f"init_alpha_{q_i}"] = 0
	vars()[f"init_det_{q_i}"] = 0

	# vars()[f"init_params_{q_i}"] = [eval(f"amp_scale['{q_i}']['X180']"), eval(f"amp_scale['{q_i}']['X90']")]
	vars()[f"init_params_{q_i}"] = [initial_amps[0], initial_amps[2]]

	if drag:
		eval(f"init_params_{q_i}").append(eval(f"init_alpha_{q_i}"))
		if det:
			eval(f"init_params_{q_i}").append(eval(f"init_det_{q_i}"))

	params_list.append(eval(f"init_params_{q_i}"))

params_array = np.array(params_list)
init_params = params_array.reshape(params_array.size)
shape = params_array.shape
n_params = params_array.size / len(q_no_list)

print(f'The initial parameters are {init_params}')

bounds_list = params_array.tolist()
for q_i in q_no_list:
	bounds_list[q_index[f"{q_i}"]][0] = (0.0, min(2 * eval(f"amp_scale['{q_i}']['X180']"), 0.9))
	bounds_list[q_index[f"{q_i}"]][1] = (0.0, min(2 * eval(f"amp_scale['{q_i}']['X90']"), 0.45))
	if drag:
		bounds_list[q_index[f"{q_i}"]][2] = (-1.0, 1.0)
		if det:
			bounds_list[q_index[f"{q_i}"]][3] = (-3.0, 3.0)

bounds_array = np.array(bounds_list, dtype=object)
bounds_params = bounds_array.reshape(params_array.size, 2)

options = {'disp': True, 'return_all': True, 'maxiter': maxiter} #'fatol': fatol,

# Formatting the initial array for config updates
init_params_array = params_array
init_config = config_updater(init_params_array, config, q_no_list, drag, det)
init_fid_list = QPU_Simultaneous_RB(init_config)
init_fid_qpu = {}
for q_i in q_no_list:
	init_fid_qpu[f"Qubit {q_i}"] = init_fid_list[q_index[f"{q_i}"]]
print(f"The initial fidelities for the qubits under consideration with RB are: \n{init_fid_qpu}")

result = minimize(fun=ORBIT, x0=init_params, args=(config, orbit, qmm, shape, drag, det), bounds=bounds_params,
                  method='Nelder-Mead', tol=1e-1, options=options)

all_params_list = result.allvecs
iter_array = []

for i in range(len(all_params_list)):
	iter_array.append(i)

	one_param_list = np.array(all_params_list[i])
	one_param_array = one_param_list.reshape(shape)

	for q_i in q_no_list:
		par_list = ["pi", "piby2"]
		vars()[f"p_pi_{q_i}"], vars()[f"p_piby2_{q_i}"] = [], []
		eval(f"p_pi_{q_i}").append(one_param_array[q_index[f"{q_i}"]][0])
		eval(f"p_piby2_{q_i}").append(one_param_array[q_index[f"{q_i}"]][1])
		if drag:
			par_list.append("alpha")
			vars()[f"p_alpha_{q_i}"] = []
			eval(f"p_alpha_{q_i}").append(one_param_array[q_index[f"{q_i}"]][2])
			if det:
				par_list.append("det")
				vars()[f"p_det_{q_i}"] = []
				eval(f"p_det_{q_i}").append(one_param_array[q_index[f"{q_i}"]][3])

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