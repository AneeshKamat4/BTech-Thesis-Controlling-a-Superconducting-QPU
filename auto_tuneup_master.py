from configuration_4qubitsv4 import *
from auto_helper_functions_spectro import *
from auto_helper_functions_init import *

"""
To be done before using this script:

1. Readout spectroscopy
2. Mixer calibration (offset and sideband) for qe, rr as needed
3. Time of flight calibration
4. Use 20dB attenuator, change RO amp to 0.5
5. Run readout spectroscopy and update rr_LO ensuring rr_IF = 20 MHz
"""

q_no = 1

n_avgs = 5e2 # SNR dependent, upper bound
n_samples = 500
anharm = 250 # MHz
qmm = QuantumMachinesManager(qm_ip)

q_props, zero2by2_props, anharm = auto_spectro(q_no, qmm, config, adc_mapping, q_LO, 
                                               n_samples, n_avgs, anharm, pi_12=False)

[q_freq, q_fwhm] = q_props
[zero2by2_freq, zero2by2_fwhm] = zero2by2_props


def config_update_freq(config, q_no, q_freq):
	'''
	Updates the config with the new parameters.
	
	'''
	rr_no = q_no
	NEW_q_IF = q_freq * u.GHz - q_LO
	new_config = config_add_elements_q_rr(config, q_no, rr_no, dac_mapping, q_LO, NEW_q_IF, q_anh, rr_LO, rr_IF,
										  pi_len_ns, piby2_len_ns, pi_rise_grft_ns, amp_scale, drag_dict,
										  mixers, mixer_corrections, ro_amp, ro_len_clk, tof, integ_len_clk,
										  optimal_readout_phase, smearing=0)

	return new_config

mod1_config = config_update_freq(config, q_no, q_freq)

# Time Rabi initialization followed by ALLXY diagnosis

qe = f"q{q_no}"
rr = f"rr{q_no}"
out = adc_mapping[rr]

wait_init = 250000 # in clocks ~ 1 ms

t_min_ns = 16
t_max_ns = 1600
dt_ns = 8 # minimum 4ns
rabi_pi_amp = 0.5

show_Rabi = True

t_min = int(t_min_ns/4) # in clocks
t_max = int(t_max_ns/4)
dt = int(dt_ns/4)
t_list_clk = np.arange(t_min, t_max, dt)

ensemble_params = (n_avgs, t_min, t_max, dt, t_list_clk, wait_init, qe, rr, out, rabi_pi_amp)

initial_amps = initialize_Rabi(mod1_config, q_no, qmm, *ensemble_params, show_Rabi=show_Rabi)

def config_update_amps(config, q_no, initial_amps):
	'''
	Updates the config with the new parameters.
	'''
	rr_no = q_no
	
	NEW_amp_scale = amp_scale
	NEW_amp_scale[f"{q_no}"]["X180"] = initial_amps[0]
	NEW_amp_scale[f"{q_no}"]["Y180"] = initial_amps[1]
	NEW_amp_scale[f"{q_no}"]["X90"] = initial_amps[2]
	NEW_amp_scale[f"{q_no}"]["Y90"] = initial_amps[3]

	new_config = config_add_elements_q_rr(config, q_no, rr_no, dac_mapping, q_LO, q_IF, q_anh, rr_LO, rr_IF,
										  pi_len_ns, piby2_len_ns, pi_rise_grft_ns, NEW_amp_scale, drag_dict,
										  mixers, mixer_corrections, ro_amp, ro_len_clk, tof, integ_len_clk,
										  optimal_readout_phase, smearing=0)

	return new_config

mod2_config = config_update_amps(mod1_config, q_no, initial_amps)


allxy_fid = opt_allxy(initial_amps, mod2_config, q_no, qmm, n_avgs, t_min, t_max, dt, t_list_clk, 
					  wait_init, qe, rr, out, rabi_pi_amp, show_Rabi=False, print_Rabi=False)

