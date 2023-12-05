from qm.QuantumMachinesManager import QuantumMachinesManager
import time
from configuration_4qubitsv4 import *
from auto_mixer_tools_visa import RohdeSchwarzFPC1000
from qm.qua import *
from scipy.optimize import minimize
import json


qmm = QuantumMachinesManager(qm_ip)
qm = qmm.open_qm(config)

address = "TCPIP0::192.168.0.101::inst0::INSTR"
calib = RohdeSchwarzFPC1000(address, qm)

qe_freq = {}
for i in range(1,n_qubits): #use range(1, 5) if 4 quantum elements
    qe_freq[f"q{i}"] = [q_LO[f"{i}"], q_IF[f"{i}"], f"mixer_q{i}"]
    qe_freq[f"rr{i}"] = [rr_LO[f"{i}"], rr_IF[f"{i}"], f"mixer_rr{i}"]
    # qe_freq[f"q12_{i}"] = [q_LO[f"{i}"], q12_IF[f"{i}"], f"mixer_q12_{i}"]

# qe_freq["cr_c2t1"] = [q_LO["2"], q_IF["1"], "mixer_cr_c2t1"]
# qe_freq["cr_c4t1"] = [q_LO["4"], q_IF["1"], "mixer_cr_c4t1"]
# qe_freq["cr_c2t3"] = [q_LO["2"], q_IF["3"], "mixer_cr_c2t3"]
# qe_freq["cr_c4t3"] = [q_LO["4"], q_IF["3"], "mixer_cr_c4t3"]

qe = "q1"

with open('Calibrations/iq_imbalance.json') as f:
    iq_imbalance = json.load(f)


def set_IQ_imbalance_get_leakage(imbalance_arr, qe, qe_freq,verbose=True):
    
    a_imb, p_imb = imbalance_arr
    qe_LO, qe_IF, mixer = qe_freq[qe]
    qm.set_mixer_correction(mixer, int(qe_IF), int(qe_LO), IQ_imbalance(a_imb, p_imb))
    time.sleep(2) ##
    leakage = calib.query_marker(1)
    if verbose:
        print("Current leakage is {0} dBm".format(leakage))
    if leakage < -93:
         return -300    # return low constant values so that the optimization loop quits
    return leakage


with program() as mixer_calibration_pulse:
    with infinite_loop_():
        
        play("const"*amp(1.0), qe)

job = qm.execute(mixer_calibration_pulse)

qe_LO, qe_IF, mixer = qe_freq[qe]
center_freq = qe_LO - qe_IF  # rr_LO
span = 50
calib.set_bandwidth(5)
calib.set_sweep_points(501)

calib.set_center_freq(center_freq)
calib.set_span(span)
calib.active_marker(1)
calib.set_marker_freq(1, center_freq)


init_vals = [iq_imbalance[qe]["a"], iq_imbalance[qe]["p"]]
#init_vals = 0, 0
print(init_vals)
# lower_bound = -60
bnds = ((-0.3, 0.3), (-0.3, 0.3))

xatol = 1e-4  # 1e-4 change in DC offset or gain/phase
# fatol = 3  # dB change tolerance
maxiter = 50  # 50 iterations should be more than enough, but can be changed.
initial_simplex = np.zeros([3, 2])

initial_simplex[0, :] = [0.1, -0.1]  #Srijita's observation: decreasing simplex size helped with convergence
initial_simplex[1, :] = [0.1, 0.1]
initial_simplex[2, :] = [-0.1, 0.1]   #[0.2, -0.2] [-0.1, 0.1]

args = [qe, qe_freq]
print(np.array(init_vals))
res = minimize(set_IQ_imbalance_get_leakage,
               np.array(init_vals),
               args= (qe, qe_freq),
               method="nelder-mead",
               bounds=bnds,
               options={
                "xatol": xatol,
#                "fatol": fatol,
                "initial_simplex": initial_simplex,
                "maxiter": maxiter,
               }
              )

print("Final reject band leakage is {0} dBm".format(calib.query_marker(1)))

center_freq = qe_LO + qe_IF
calib.set_center_freq(center_freq)
calib.set_marker_freq(1, center_freq)
time.sleep(2)
print("Upconverted sideband power is {0} dBm".format(calib.query_marker(1)))
print(f"Calibrated IQ imbalance tuple is {res.x}")


a, p = res.x
iq_imbalance[qe]["a"] = a
iq_imbalance[qe]["p"] = p

with open('Calibrations/iq_imbalance.json', "w") as outfile:
    json.dump(iq_imbalance, outfile, indent=4)

#qm.close()



