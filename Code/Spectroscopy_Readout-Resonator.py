from qm import SimulationConfig
from qm.qua import *
from qm import LoopbackInterface
from qm.QuantumMachinesManager import QuantumMachinesManager
from configuration_4qubitsv4 import *
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
save_data = False
simulate = False

###################
# The QUA program #
###################

check_e_delay = False

rr_no = 1
rr = f"rr{rr_no}"
out = adc_mapping[rr]
ro_len = ro_len_clk[str(rr_no)]
rep_rate_clk = 2500
rr_LO = config["elements"][rr]["mixInputs"]["lo_frequency"]

f_min = -10e6
f_max = 50e6
df = 0.1e6

if check_e_delay:
    f_min = -250e6
    f_max = 250e6
    df = 1e6

freq_list = np.arange(f_min, f_max, df)
zeros = np.where(freq_list == 0)
if not len(zeros[0]) == 0:
    zero_i = zeros[0][0]

with program() as rr_spec:
    n = declare(int)
    I = declare(fixed)
    I_st = declare_stream()
    Q = declare(fixed)
    Q_st = declare_stream()
    f = declare(int)

    with for_(n, 0, n < 1500, n + 1):
        with for_(f, f_min, f < f_max, f + df):

            update_frequency(rr, f)
            wait(rep_rate_clk, rr)
            measure("readout"*amp(1.0), rr, None,
                    demod.full("integW_cos", I, out),
                    demod.full("integW_minus_sin", Q, out))

            save(I, I_st)
            save(Q, Q_st)

    with stream_processing():
        I_st.buffer(len(freq_list)).average().save('I')
        Q_st.buffer(len(freq_list)).average().save('Q')

######################################
# Open Communication with the Server #
######################################
qmm = QuantumMachinesManager(host=qm_ip)

####################
# Simulate Program #
####################
if simulate:
    simulation_config = SimulationConfig(
        duration=200000,
        simulation_interface=LoopbackInterface([("con1", 9, "con1", 1), ("con1", 10, "con1", 2)]))
    job = qmm.simulate(config, rr_spec, simulation_config)
    # get DAC and digital samples
    samples = job.get_simulated_samples()
    # plot all ports:
    samples.con1.plot()
    raise Halted()

#############
# execution #
#############
qm = qmm.open_qm(config)
job = qm.execute(rr_spec)
job.result_handles.wait_for_all_values()
I = job.result_handles.get("I").fetch_all()
Q = job.result_handles.get("Q").fetch_all()

############
# analysis #
############
freq_list = 1e-9*(rr_LO + freq_list)
sig = I + 1j*Q
if not len(zeros[0]) == 0:
    freq_list = np.delete(freq_list, zero_i)
    sig = np.delete(sig, zero_i)

e_delay = elec_delay_ns[str(rr_no)]
p_offset =phase_offset_rad[str(rr_no)]
# e_delay = 254.9
# p_offset = 4
sig_corrected = sig*np.exp(1j*2*np.pi*freq_list*e_delay + 1j*p_offset)
phase = np.angle(sig_corrected)
real = np.real(sig_corrected)
imag = np.imag(sig_corrected)
f_res_i = np.argmin(abs(sig))
f_res = freq_list[f_res_i]
plt.plot(freq_list, phase, label = "phase")
#plt.plot(freq_list, imag, label = "Imag")
plt.axvline(x=f_res, linestyle="--")
plt.legend()
plt.figure()
plt.plot(freq_list, phase)
plt.legend()
plt.xlabel("Frequency (GHz)")
plt.title(f'Cavity Spectroscopy (Phase) : Cavity Freqeuncy = {f_res} GHz')
plt.axvline(x=f_res, linestyle = "--")
plt.grid()
plt.show()

if not check_e_delay:
    plt.figure()
    plt.plot(freq_list, real)
    plt.title(f'Cavity Spectroscopy (Real) : Cavity Freqeuncy = {f_res} GHz')
    plt.xlabel("Frequency (GHz)")
    plt.axvline(x=f_res, linestyle = "--")
    plt.grid()
    plt.show()

    plt.figure()
    plt.plot(freq_list, np.abs(sig))
    plt.title(f'Cavity Spectroscopy (Magnitude) : Cavity Freqeuncy = {f_res} GHz')
    plt.xlabel("Frequency (GHz)")
    plt.axvline(x=f_res, linestyle="--")
    plt.grid()
    plt.show()

def extract_params(z1, z2, f1, f2):

    e_delay = -(np.angle(z1) - np.angle(z2))/(2*np.pi*(f1 - f2))
    return e_delay

data = np.transpose([freq_list, np.abs(sig), phase, real, imag])
if save_data :
    file_saver_(data, file_name=__file__, suffix=f"{rr}_amp_1.0", master_folder=ExpName, header_string="Frequency (GHz), Magnitude, Phase, Real, Imag", time_stamp=False)
