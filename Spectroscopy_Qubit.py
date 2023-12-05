from qm import LoopbackInterface, SimulationConfig
from qm.qua import *
from qm.QuantumMachinesManager import QuantumMachinesManager
from configuration_4qubitsv4 import *
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
save_data = False

###################
# The QUA program #
###################
q_no = 1
qe = f"q{q_no}"
rr = f"rr{q_no}"
#rr ="rr2"
out = adc_mapping[rr]

f_LO = q_LO[qe[-1]]
# temp_f = 46.5
# f_min = (temp_f-0.2)*u.MHz # 188.3* u.MHz
# f_max = (temp_f+0.2)*u.MHz #188.4* u.MHz
f_min = -12 * u.MHz # 188.3* u.MHz
f_max = 12 * u.MHz #188.4* u.MHz
df = 0.03 * u.MHz #0.001 * u.MHz
q_amp = 0.7
freqs = np.arange(f_min, f_max, df)

with program() as qubit_spec:
    n = declare(int)
    I = declare(fixed)
    I_st = declare_stream()
    Q = declare(fixed)
    Q_st = declare_stream()
    f = declare(int)

    with for_(n, 0, n < 1500, n + 1):
        with for_(f, f_min, f < f_max, f + df):
            wait(250000, qe)
            update_frequency(qe, f)
            play("const"*amp(q_amp), qe, duration=50000)  #duration=20000
            align(rr, qe)
            measure("readout", rr, None,
                    demod.full("integW_cos", I, out),
                    demod.full("integW_minus_sin", Q, out))
            save(I, I_st)
            save(Q, Q_st)

    with stream_processing():
        I_st.buffer(len(freqs)).average().save('I')
        Q_st.buffer(len(freqs)).average().save('Q')

######################################
# Open Communication with the Server #
######################################
qmm = QuantumMachinesManager(qm_ip)

####################
# Simulate Program #
####################
simulate = False
if simulate:
    simulation_config = SimulationConfig(
        duration=200000,
        simulation_interface=LoopbackInterface([("con1", 9, "con1", 1), ("con1", 10, "con1", 2)]))
    job = qmm.simulate(config, qubit_spec, simulation_config)
    # get DAC and digital samples
    samples = job.get_simulated_samples()
    # plot all ports:
    samples.con1.plot()
    raise Halted()

#############
# execution #
#############
qm = qmm.open_qm(config)
job = qm.execute(qubit_spec)
res_handles = job.result_handles
I_handle = job.result_handles.get("I")
Q_handle = job.result_handles.get("Q")
#job.result_handles.wait_for_all_values()

plt.figure()
plt.title("Qubit Spectroscopy")
I_handle.wait_for_values(1)
Q_handle.wait_for_values(1)
while res_handles.is_processing():

    I = I_handle.fetch_all()
    Q = Q_handle.fetch_all()
    sig = I + 1j * Q
    plt.clf()
    # plt.plot(1e-6*(freqs), np.abs(sig))
    plt.plot(1e-6*(freqs), I, label="I")
    # plt.plot(1e-6*(freqs), Q, label="Q")
    plt.grid()
    plt.legend()
    plt.pause(1)

I = job.result_handles.get("I").fetch_all()
Q = job.result_handles.get("Q").fetch_all()

############
# analysis #
############

freqs = 1e-9*(f_LO + freqs)
plt.figure()
sig = I + 1j * Q

plt.title('qubit spectroscopy')
# plt.plot(freqs, np.abs(sig))
plt.plot(freqs, I)
#plt.plot(freqs, Q)
plt.xlabel("Frequency (GHz)")
plt.grid()
plt.show()

data = np.transpose([freqs, I, Q])


vol = -30.34
if save_data :
    file_saver_(data,file_name=__file__, suffix= f"{q_amp}_Vol-{vol}", master_folder= ExpName,header_string="Frequency (GHz), I, Q",time_stamp=False)