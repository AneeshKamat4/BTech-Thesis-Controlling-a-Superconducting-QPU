from qm import SimulationConfig
from qm.qua import *
from qm import LoopbackInterface
from qm.QuantumMachinesManager import QuantumMachinesManager
from configuration_4qubitsv4 import *
import numpy as np
import matplotlib
#matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
from qualang_tools.analysis.discriminator import two_state_discriminator

save_data = False

qmm = QuantumMachinesManager(qm_ip)
q_no = 1
qe = f"q{q_no}"
rr = f"rr{q_no}"
#rr="rr2"
out = adc_mapping[rr]

rep_rate_clk = 250000
wait_rr = 16
pi_len = pi_len_ns[str(q_no)]
ro_len = ro_len_clk[str(q_no)]
###################
# The QUA program #
###################

dem = -1.742e-3
a_rr = 1.0
with program() as IQ_blobs:

    n = declare(int)
    I0 = declare(fixed)
    I0_st = declare_stream()
    Q0 = declare(fixed)
    Q0_st = declare_stream()
    I1 = declare(fixed)
    I1_st = declare_stream()
    Q1 = declare(fixed)
    Q1_st = declare_stream()

    with for_(n, 0, n < 5000, n + 1):

        # wait(rep_rate_clk - pi_len - wait_rr - ro_len, qe)
        wait(rep_rate_clk, qe)
        play("const" * amp(1), qe_c, cool_pulse_time)
        align(qe_c, qe)
        wait(cool_wait_time_ns, qe)
        play("I", qe)
        align(qe, rr)
        wait(wait_rr, rr)
        measure("readout" * amp(a_rr), rr, None,
                demod.full("integW_cos", I0, out),
                demod.full("integW_minus_sin", Q0, out))
        save(I0, I0_st)
        save(Q0, Q0_st)

        align(rr, qe)
        # wait(rep_rate_clk - pi_len - wait_rr - ro_len, qe)
        wait(rep_rate_clk, qe)
        play("X180", qe)
        align(qe, rr)
        wait(wait_rr, rr)
        measure("readout" * amp(a_rr), rr, None,
                demod.full("integW_cos", I1, out),
                demod.full("integW_minus_sin", Q1, out))
        save(I1, I1_st)
        save(Q1, Q1_st)

    with stream_processing():
        I0_st.save_all('I0')
        Q0_st.save_all('Q0')
        I1_st.save_all('I1')
        Q1_st.save_all('Q1')

qm = qmm.open_qm(config)
job = qm.execute(IQ_blobs)

I0_handle = job.result_handles.get("I0")
Q0_handle = job.result_handles.get("Q0")
I1_handle = job.result_handles.get("I1")
Q1_handle = job.result_handles.get("Q1")
job.result_handles.wait_for_all_values()

I0 = I0_handle.fetch_all()["value"]
Q0 = Q0_handle.fetch_all()["value"]
I1 = I1_handle.fetch_all()["value"]
Q1 = Q1_handle.fetch_all()["value"]


angle, threshold, fidelity, gg, ge, eg, ee = two_state_discriminator(I0, Q0, I1, Q1, b_print=True, b_plot=True)

a = ro_amp[str(q_no)]
if save_data is True:
    file_saver_(I0, file_name=__file__, suffix=f"flux_-38.020_wo_20dB_I0_Amp_{a}", master_folder=ExpName, time_stamp=False)  # _wo_20dB_
    file_saver_(Q0, file_name=__file__, suffix=f"flux_-38.020_wo_20dB_Q0_Amp_{a}", master_folder=ExpName, time_stamp=False)
    file_saver_(I1, file_name=__file__, suffix=f"flux_-38.020_wo_20dB_I1_Amp_{a}", master_folder=ExpName, time_stamp=False)
    file_saver_(Q1, file_name=__file__, suffix=f"flux_-38.020_wo_20dB_Q1_Amp_{a}", master_folder=ExpName, time_stamp=False)
