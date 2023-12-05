from qm.QuantumMachinesManager import QuantumMachinesManager
import time
from configuration_4qubitsv4 import *
from auto_mixer_tools_visa import RohdeSchwarzFPC1000
from qm.qua import *
from helper_functionsv4 import *
import json


qmm = QuantumMachinesManager(qm_ip)
qm = qmm.open_qm(config)

address = "TCPIP0::192.168.0.101::inst0::INSTR"
calib = RohdeSchwarzFPC1000(address, qm)

qe_freq = {}
for i in range(1,n_qubits): #use range(1, 5) if 4 quantum elements
    qe_freq[f"q{i}"] = [q_LO[f"{i}"], dac_mapping[f"q{i}"]]
    qe_freq[f"rr{i}"] = [rr_LO[f"{i}"], dac_mapping[f"rr{i}"]]


qe = "q1"  #RF switch A pos 1 => q2, rr2 else q1, rr1

with program() as mixer_calibration_pulse:
    with infinite_loop_():
        play("const" * amp(0.0), qe)
        # for qe_temp in qe_freq.keys():
        #     play("const" * amp(0.0), qe_temp)
        
job = qm.execute(mixer_calibration_pulse)
center_freq, qe_ch = qe_freq[qe] 
span = 50
calib.set_bandwidth(5)
calib.set_sweep_points(501)

calib.set_center_freq(center_freq)
calib.set_span(span)
calib.active_marker(1)
calib.set_marker_freq(1, center_freq)

def set_dc_offset_get_leakage(offset_arr, qe, verbose=True):
    
    I_offset, Q_offset = offset_arr[0], offset_arr[1]
    
    qm.set_output_dc_offset_by_element(qe, "I", I_offset)
    qm.set_output_dc_offset_by_element(qe, "Q", Q_offset)
    
    time.sleep(2)
    
    leakage = calib.query_marker(1)
    if verbose:
        print("Current leakage is {0} dBm".format(leakage))
    
    if leakage < -95:
         return -300
        
    return leakage



from scipy.optimize import minimize

init_vals = [0.0, 0.0]
bnds = ((-0.3, 0.3), (-0.3, 0.3))

xatol = 1e-4  # 1e-4 change in DC offset or gain/phase
#fatol = 3  # dB change tolerance
maxiter = 50  # 50 iterations should be more then enough, but can be changed.
initial_simplex = np.zeros([3, 2])
# initial_simplex[0, :] = [0, 0]
# initial_simplex[1, :] = [0, -0.2]
# initial_simplex[2, :] = [-0.1, 0]

initial_simplex[0, :] = [0, -0.1] #[0, -0.2]
initial_simplex[1, :] = [0.1, 0.1] #[0.2, 0.2]
initial_simplex[2, :] = [-0.1, 0.1] #[-0.2, 0.2]

res = minimize(set_dc_offset_get_leakage, 
               np.array(init_vals),
               args=qe, 
               method="nelder-mead", 
               bounds=bnds, 
               options={
                "xatol": xatol,
#                "fatol": fatol,
                "initial_simplex": initial_simplex,
                "maxiter": maxiter,
               }
              )
print("Final leakage is {0} dBm".format(calib.query_marker(1)))
print(f"Calibrated Mixer offsets for {qe} are {res.x}")


with open('Calibrations/dc_offsets.json') as f:
    dc_offsets = json.load(f)

dc_offsets[f"con{qe_ch[0]}"][f"{qe_ch[1][0]}"] = res.x[0]
dc_offsets[f"con{qe_ch[0]}"][f"{qe_ch[1][1]}"] = res.x[1]

with open('Calibrations/dc_offsets.json', "w") as outfile:
    json.dump(dc_offsets, outfile, indent=4)
