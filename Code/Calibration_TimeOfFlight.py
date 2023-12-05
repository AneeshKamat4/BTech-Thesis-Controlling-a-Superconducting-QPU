from qm import SimulationConfig
from qm.qua import *
from qm import LoopbackInterface
from qm.QuantumMachinesManager import QuantumMachinesManager
from configuration_4qubitsv4 import *
import numpy as np
from matplotlib import pyplot as plt

###################
# The QUA program #
###################
simulate = False

update_adc_offsets = False

rr_no = 1
rr = f"rr{rr_no}"
ro_len = ro_len_clk[str(rr_no)]
rep_rate_clk = 250000 # 1 ms
# rep_rate_clk = 10000 # 1 ms


with program() as tof_cal:

    n = declare(int)
    adc_st = declare_stream(adc_trace=True)

    with for_(n, 0, n < 500, n+1):

        reset_phase(rr)
        wait(rep_rate_clk - ro_len, rr)

        measure("readout" * amp(1.0), rr, adc_st)


    with stream_processing():
        adc_st.input1().average().save("adc1")
        adc_st.input2().average().save("adc2")

######################################
# Open Communication with the Server #
######################################
qmm = QuantumMachinesManager(host=qm_ip)

####################
# Simulate Program #
####################

if simulate:
    simulation_config = SimulationConfig(
        duration=2000,
        simulation_interface=LoopbackInterface([("con1", 9, "con1", 1), ("con1", 10, "con1", 2)]))
    job = qmm.simulate(config, tof_cal, simulation_config)
    job.result_handles.wait_for_all_values()
    adc1 = job.result_handles.get("adc1").fetch_all()
    adc2 = job.result_handles.get("adc2").fetch_all()

    raise Halted()

#############
# execution #
#############
qm =qmm.open_qm(config)
job = qm.execute(tof_cal)
res_handles = job.result_handles
adc1_handle = res_handles.get("adc1")
adc2_handle = res_handles.get("adc2")
res_handles.wait_for_all_values()
adc1 = job.result_handles.get("adc1").fetch_all()
adc2 = job.result_handles.get("adc2").fetch_all()

############
# analysis #
############
adc1_offset = -np.mean(adc1)*2**-12
adc2_offset = -np.mean(adc2)*2**-12

con = f"con{dac_mapping[rr][0]}"

with open('Calibrations/adc_offsets.json') as f:
    adc_offsets = json.load(f)

adc_offsets[con]["1"] = adc_offsets[con]["1"] + adc1_offset
adc_offsets[con]["2"] = adc_offsets[con]["2"] + adc2_offset

if update_adc_offsets:
    with open('Calibrations/adc_offsets.json', "w") as outfile:
        json.dump(adc_offsets, outfile, indent=4)

plt.figure()
plt.title('time-of-flight calibration analysis')
plt.plot(adc1)
plt.plot(adc2)
plt.legend(
    [f"adc1: offset = {adc1_offset:.4f}",
     f"adc2: offset = {adc2_offset:.4f}"])
plt.show()




