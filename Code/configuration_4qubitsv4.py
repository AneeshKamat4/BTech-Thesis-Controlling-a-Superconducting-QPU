from config_builderv4 import *
# from helper_functionsv4 import *
from qualang_tools.units import unit
from qualang_tools.plot import interrupt_on_close
from qualang_tools.results import progress_counter, fetching_tool
from instrumentlib import RhodeandSchwarz_SGS100A
import matplotlib
import json

# ========== Added by Aneesh ======
import numpy as np
qop_ip = "192.168.0.106"
qm_ip = "192.168.0.47"

# ===== End of Additions =========
u = unit()
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

TOF_testing = False
initialize_RF_sources = False

path = r"D:\Experimental Data\2023-11-10 FLX_dev_Ta_03, 3DT_Ta_05AC"
ExpName = r"2023-11-10 FLX_dev_Ta_03, 3DT_Ta_05AC\Aneesh_autotuneup"
# path = r"D:\Experiments\2023-05-20 Ringv1_4Q_INDIQ05_Run2"
# ExpName = "2023-05-20 Ringv1_4Q_INDIQ05_Run2"

n_qubits = 4

# ====================== DAC Mapping ========================================
dac_mapping = {"q1": [1, [1, 2]], "q2": [1, [3, 4]], "q3": [2, [1, 2]], "q4": [2, [3, 4]],
               "rr1": [1, [7, 8]], "rr2": [1, [9, 10]], "rr3": [2, [7, 8]], "rr4": [2, [9, 10]]}
               # "cr_c2t1": [1, [3, 4]], "cr_c4t1": [2, [3, 4]], "cr_c2t3": [1, [3, 4]], "cr_c4t3": [2, [3, 4]]}

adc_mapping = {"rr1": "out1", "rr2": "out2", "rr3": "out1", "rr4": "out2"}

# =======================Frequencies==========================================

qe_list = ["q1", "q2", "q3", "q4",
           "rr1", "rr2", "rr3", "rr4"]
           # "cr_c2t1", "cr_c4t1", "cr_c2t3", "cr_c4t3",
           # "q12_1", "q12_2", "q12_3", "q12_4"]

CrossKerr = {"1": 0*0.1174 * u.MHz,              # Calculate and update later using simulation/ conditional Ramsey
             "2": 0*0.1000 * u.MHz,
             "3": 0*0.1118 * u.MHz,
             "4": 0*0.1292 * u.MHz,
             }

atten = {"1": True,
        "2": False,
        "3": True,
        "4": False,
        }

q_LO = {"1": 5.3 * u.GHz,
        "2": 5 * u.GHz,
        "3": 4.9 * u.GHz,
        "4": 4.9 * u.GHz,
        }
# 5.3413436472942495
q_IF = {"1": 41.343 * u.MHz + CrossKerr["1"],    #Ring 1
        "2": 88.685 * u.MHz + CrossKerr["2"],   #Ring 11
        "3": -46.893 * u.MHz + CrossKerr["3"],   #Ring 5
        "4": 104.179 * u.MHz + CrossKerr["4"],  #Ring 7 104.23544
        }

q_anh = {"1": -260 * u.MHz,
         "2": -309 * u.MHz,
         "3": -312.6 * u.MHz,
         "4": -308 * u.MHz,
        }

rr_LO = {"1": 7.1442 * u.GHz,
         "2": 7.1442 * u.GHz, #7.43615 + 20 IF
         "3": 7.36243 * u.GHz,
         "4": 7.36243 * u.GHz,
         }

rr_IF = {"1": 20 * u.MHz,
         "2": 20 * u.MHz,
         "3": 40 * u.MHz, #40
         "4": 36.755 * u.MHz, # 36.985  * u.MHz
         }

pi_12_readout = {"1": False,
                 "2": False,
                 "3": False,
                 "4": False,
                }

# ======================Readout Parameters===========================================

tof = {"1": 300,
       "2": 300,
       "3": 300,
       "4": 300,
       }

gain = {"1" : {"1": 20, "2": 20},  #Set ADC gain here
        "2" : {"1": 10, "2": 10}}

ro_amp = {"1": 0.5, # 0.05 wo20dB or 0.5 with 20dB for first coherence set
          "2": 2*0.08,
          "3": 0.192, # 0.4,
          "4": 0.07,
          }

# For TOF Testing
if TOF_testing:

    ro_amp = {"1": 1.0,
              "2": 1.0,
              "3": 1.0,
              "4": 1.0,
              }

    rr_IF = {"1": 12 * u.MHz,
             "2": 12 * u.MHz,
             "3": 12 * u.MHz,
             "4": 12 * u.MHz,
             }

ro_len_clk = {"1": 500,
              "2": 500,
              "3": 500,
              "4": 500,
              }

integ_len_clk = {"1": 400,
                 "2": 400,
                 "3": 500,
                 "4": 500,
                 }

optimal_readout_phase = {"rr1": -72.9*(np.pi/180), #
                         "rr2": -36.2*(-np.pi/180),
                         "rr3": -28.3*(np.pi/180),
                         "rr4": 57.1*(np.pi/180),
                         }

demarcations = {"1": 2.515e-05,
                "2": -1.997e-04,
                "3": 5.109e-05,
                "4": -1.190e-05,
                }

elec_delay_ns = {"1": 254,
                "2": 273,
                "3": 277.5,
                "4": 277.5,
                }

phase_offset_rad = {"1": 2,
                "2": 2.199,
                "3": 3.44,
                "4": -2.074-1.31,
                }
# ==================Control Parameters===============================================

pi_rise_grft_ns = 10
pi_len_ns = {"1": 312,
             "2": 52,
             "3": 72,
             "4": 52, #124
             }

piby2_rise_grft_ns = 10
piby2_len_ns = {"1": 312,
                "2": 52,
                "3": 72,
                "4": 52, #124
                }
#
# cr_tail_ns = 16
# cr_len_ns = {"cr_c2t1": 300,                 # CR gate length except 16ns rise and 16ns fall
#              "cr_c4t1": 204,
#              "cr_c2t3": 284,
#              "c4_c4t3": 244,
#              }
#
# cr_amp = {"cr_c2t1": 0.3, "cr_ac_c2t1": 0.03,
#           "cr_c4t1": 0.05, "cr_ac_c4t1": 0.03, # 0.25/2  #0.25/3 tested, good RB and echo
#           "cr_c2t3": 0.6, "cr_ac_c2t3": 0.03}
#
# cr_phase = {"cr_c2t1":  0.43843843843843844, "cr_ac_c2t1": 0.16716716716716717,
#             "cr_c4t1":  0.3743743743743743, "cr_ac_c4t1": 0.1701701701701701,
#             "cr_c2t3":   0.4434434434434434, "cr_ac_c2t3": 0.07107107107107108}



calib_vals = {"1": {"amin": 0.4917043021510755, "amax": 0.5222, "da": 0.0005, "n_pulses": 11},
              "2": {"amin": 0.3, "amax": 0.315, "da": 0.0001, "n_pulses": 17},
              "3": {"amin": 0.59, "amax": 0.61, "da": 0.0001, "n_pulses": 17},
              "4": {"amin": 0.12, "amax": 0.1275, "da": 0.00005, "n_pulses": 17}
              }

amp_scale = {}
drag_dict = {}
q_no_list = [i+1 for i in range(n_qubits)]

for q_i in q_no_list:
    vars()[f"amp_pi_{q_i}"] = 0.5
    vars()[f"amp_piby2_{q_i}"] = 0.25
    vars()[f"alpha_{q_i}"] = 1
    vars()[f"det_{q_i}"] = 0

# Custom initialization from one Power Rabi
amp_pi_1 = 1 #0.4
amp_piby2_1 = 1
alpha_1 = -0.553
#[0.43999519 0.43999519 0.2199976  0.2199976 ]
'''
"1": {"X180": 0.77968984, "Y180": 0.77968984, "X90": 0.40955477, "Y90": 0.40955477},#
             "2": {"X180": 0.70502751, "Y180": 0.70502751, "X90": 0.3638069, "Y90": 0.3638069},  #
'''

for q_i in q_no_list:
    if not f"{q_i}" in amp_scale: 
        amp_scale[f"{q_i}"] = {}
    amp_scale[f"{q_i}"]["X180"] = eval(f"amp_pi_{q_i}")
    amp_scale[f"{q_i}"]["Y180"] = eval(f"amp_pi_{q_i}")
    amp_scale[f"{q_i}"]["X90"] = eval(f"amp_piby2_{q_i}")
    amp_scale[f"{q_i}"]["Y90"] = eval(f"amp_piby2_{q_i}")

    if not f"{q_i}" in drag_dict:
        drag_dict[f"{q_i}"] = {}
    drag_dict[f"{q_i}"]["alpha"] = eval(f"alpha_{q_i}")
    drag_dict[f"{q_i}"]["det"] = eval(f"det_{q_i}") * u.MHz



'''
amp_scale = {"1": {"X180": 1, "Y180":  1, "X90": 1, "Y90": 1},
             "2": {"X180": 0.30521710, "Y180": 0.30507253626, "X90":  0.15164707353, "Y90":  0.15164707353},
             "3": {"X180": 1 , "Y180": 1 , "X90": 1, "Y90": 1},
             "4": {"X180": 1, "Y180": 0.123298649, "X90": 0.06149624812, "Y90": 0.06149624812}, #0.123298649 - X180
            }'''
'''
drag_dict = {"1": {"alpha": 0, "det": 0 * u.MHz}, #0.6718029715665588, 0.4237131091071187, -0.22100726252633543
             "2": {"alpha": 0, "det": 0 * u.MHz},
             "3": {"alpha": -0.5616, "det": 0 * u.MHz}, #-0.5614887526086284
             "4": {"alpha": -1.74, "det": 0 * u.MHz} #0.9576946217201819
             }'''

# 0.6660086451542028
# 0.668641164169927

# drag_coef = 0
# anharmonicity = -308 * u.MHz
# AC_stark_detuning = 0 * u.MHz

# end of additions by Aneesh


#
# # =================== Dictionaries for 1-2 transition =================================
#
# q12_IF = {"1": (60.98 - 303.8) * u.MHz,    #Ring 1
#           "2": (88.67 - 309.28) * u.MHz,   #Ring 11
#           "3": -359.28 * u.MHz,   #Ring 5
#           "4": -204.36 * u.MHz,  #Ring 7
#          }
#
# piby2_12_len_ns = {"1": 52,
#                 "2": 52,
#                 "3": 124,
#                 "4": 124,
#                 }
#
# pi_12_len_ns = {"1": 52,
#                 "2": 52,
#                 "3": 124,
#                 "4": 124,
#                 }
#
# amp_12_scale = {"1": {"X180": 0.4695898, "Y180":  0.4896248124, "X90": 0.243341670, "Y90": 0.2433266633},
#                 "2": {"X180": 0.25472, "Y180": 0.30507253626, "X90":  0.15164707353, "Y90":  0.15164707353},
#                 "3": {"X180": 0.503266, "Y180": 0.211735, "X90": 0.3202251, "Y90": 0.3202251},  #0.5981540770
#                 "4": {"X180": 0.0409754, "Y180": 0.1104652326, "X90": 0.05494747, "Y90": 0.05494747}, #0.1233579289
#                 }
#
#

# ==================Mixer Parameters ===============================================

# =================== DC Offsets ===============================================
with open('./Calibrations/dc_offsets.json') as f:
    dc_offsets = json.load(f)

# ================ List of Mixers ===============================================
mixers = {}
for qe in qe_list:
    mixers[qe] = "mixer_" + qe

# ================ IQ Imbalance Correction Matrices ===============================================
with open('./Calibrations/iq_imbalance.json') as f:
    iq_imbalance = json.load(f)

mixer_corrections = {}
for qe in qe_list:
    a = iq_imbalance[qe]["a"]
    p = iq_imbalance[qe]["p"]
    mixer_corrections[qe] = IQ_imbalance(a, p)

# ============== ADC Offets ===============================================
with open('./Calibrations/adc_offsets.json') as f:
    adc_offsets = json.load(f)


# # ======================= Setting up RF sources with above parameters =========================

def setup_RF_Source(rf_ip, rf_f, rf_p):

    rf = RhodeandSchwarz_SGS100A(rf_ip)
    rf.set_freq_Hz(rf_f)
    rf.set_power_dB(rf_p)
    rf.output_on()
    rf.close_connection()


rf_rr1_ip = "TCPIP0::192.168.0.121::inst0::INSTR"
rf_rr2_ip = "TCPIP0::192.168.0.108::inst0::INSTR"
rf_rr34_ip = "TCPIP0::192.168.0.200::inst0::INSTR"

rf_q12_ip = "TCPIP0::192.168.0.104::inst0::INSTR"
rf_q34_ip = "TCPIP0::192.168.0.103::inst0::INSTR"

rr_LO_dBm = 19
q_LO_dBm = 17

if initialize_RF_sources:
    setup_RF_Source(rf_rr1_ip, rr_LO["1"], rr_LO_dBm)
    setup_RF_Source(rf_rr2_ip, rr_LO["2"], rr_LO_dBm)
    setup_RF_Source(rf_rr34_ip, rr_LO["3"], rr_LO_dBm+3)
    setup_RF_Source(rf_q12_ip, q_LO["1"], q_LO_dBm)





# ======================== QM Config Starts Here =========================================

config = {
    "version": 1,
}

config = config_add_controller(config, 1, dc_offsets, adc_offsets)
config = config_add_controller(config, 2, dc_offsets, adc_offsets)
# config = config_add_controller(config, 3, dc_offsets, adc_offsets)

config = config_add_common_elements(config)

for i in range(1, n_qubits+1):
    q_no = i
    rr_no = i

    config = config_add_elements_q_rr(config, q_no, rr_no, dac_mapping, q_LO, q_IF, q_anh, rr_LO, rr_IF, 
                                      pi_len_ns, piby2_len_ns, pi_rise_grft_ns, amp_scale, drag_dict,
                                      mixers, mixer_corrections, ro_amp, ro_len_clk, tof, integ_len_clk, 
                                      optimal_readout_phase, smearing=0)


# config = config_add_rise_fall(config, cr_tail_ns)
#
# # Check CR config elements
# config = config_add_crgate(config, 2, 1, dac_mapping, q_LO, q_IF, mixers, mixer_corrections)
# config = config_add_crgate(config, 4, 1, dac_mapping, q_LO, q_IF, mixers, mixer_corrections)
# config = config_add_crgate(config, 2, 3, dac_mapping, q_LO, q_IF, mixers, mixer_corrections)
# config = config_add_crgate(config, 4, 3, dac_mapping, q_LO, q_IF, mixers, mixer_corrections)

#
# #Add quantum elements for 1-2 transition
# for i in range(1, n_qubits+1):
#
#     q_no = i
#     config = config_add_q12(config, q_no, dac_mapping, q_LO, q12_IF, pi_12_len_ns, piby2_12_len_ns,
#                             pi_rise_grft_ns, amp_12_scale, mixers, mixer_corrections)

# import pprint
# pprint.pprint(config, width=2)


