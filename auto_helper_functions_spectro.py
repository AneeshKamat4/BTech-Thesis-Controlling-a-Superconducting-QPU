# Master helper functions for automated spectroscopy
from qm.qua import *
from qm.QuantumMachinesManager import QuantumMachinesManager
from qualang_tools.units import unit
u = unit()
from macros import *

import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from scipy import signal as sgn


def q_spectro(q_no, adc_mapping, q_LO, 
              n_samples, n_avgs, q_amp, 
              f_min, f_max, pi_12=False):
    
    
    qe = f"q{q_no}"
    rr = f"rr{q_no}"
    out = adc_mapping[rr]
    f_LO = q_LO[f"{q_no}"]
    freqs = np.linspace(f_min, f_max, n_samples)
    df = freqs[2] - freqs[1]
    # freqs = np.arange(f_min, f_max, df)

    with program() as qubit_spec:
        n = declare(int)
        I = declare(fixed)
        I_st = declare_stream()
        Q = declare(fixed)
        Q_st = declare_stream()
        f = declare(int)

        with for_(n, 0, n < n_avgs, n + 1):
            with for_(f, f_min, f < f_max+1, f + df):
                wait(20000, qe)
                update_frequency(qe, f)
                play("const"*amp(q_amp), qe, duration=20000)
                align(rr,qe)
                # measure("readout", rr, None,
                #         demod.full("integW_cos", I, out),
                #         demod.full("integW_minus_sin", Q, out))
                
                measure_macro(qe, rr, out, I, Q, pi_12)
                save(I, I_st)
                save(Q, Q_st)

        with stream_processing():
            I_st.buffer(len(freqs)).average().save('I')
            Q_st.buffer(len(freqs)).average().save('Q')
    
    return qubit_spec, f_LO, freqs


def smooth_filter(data, window): # moving_average_numpy
    weights = np.repeat(1.0, window) / window
    return np.convolve(data, weights, mode='valid')


# Find the Signal to Noise ratio
def S2N(norm_sig): # normalized signal input
    window_size = int(max(len(norm_sig)//100, 10)) # 10 works for 500 samples
    fil_sig = smooth_filter(norm_sig, window_size)
    prefix = np.zeros(window_size//2)
    suffix = np.zeros(window_size//2 - 1)
    
    signal = np.append(np.append(prefix, fil_sig), suffix)
    noise = norm_sig - signal
    noise_std = np.std(noise)
    fil_noise = [min(2*np.abs(noise_std), np.abs(noise[i])) for i in range(len(noise))]

    P_signal = np.sum(np.square(signal))
    P_noise = np.sum(np.square(fil_noise))

    S2N = P_signal/P_noise
    return S2N, [signal, fil_noise]


def run_spectro(q_no, qmm, config, adc_mapping, q_LO, 
                n_samples, n_avgs, q_amp, f_min, f_max, 
                pi_12=False):
        
    qua_prog, f_LO, freqs = q_spectro(q_no, adc_mapping, q_LO, 
                                      n_samples, n_avgs, q_amp, 
                                      f_min, f_max, pi_12=pi_12)
    
    qm = qmm.open_qm(config)
    job = qm.execute(qua_prog)
    res_handles = job.result_handles
    I_handle = job.result_handles.get("I")
    Q_handle = job.result_handles.get("Q")
    #job.result_handles.wait_for_all_values()

    plt.figure()
    I_handle.wait_for_values(1)
    Q_handle.wait_for_values(1)

    snr = 0
    while (snr < 1 and res_handles.is_processing()):

        I = I_handle.fetch_all()
        Q = Q_handle.fetch_all()
        signal = I + 1j * Q
        phase = np.angle(signal)
        amplitude = np.abs(signal)
        sig = amplitude #phase

        # Ensure signal has peaks not troughs
        if (np.max(sig)-np.median(sig)) < (np.median(sig)-np.min(sig)):
            sig *= -1

        # Normalized signal, baseline oscillates around 0
        norm_sig = (sig - np.median(sig))/(np.max(sig) - np.median(sig))

        snr, [pure_sig, nois_sig] = S2N(norm_sig)
        # snr = np.mean(sig)/np.std(sig)

        plt.clf()
        plt.plot(1e-6*(freqs), norm_sig,label="Data") # c='C1',
        # plt.plot(1e-6*(freqs), pure_sig, linestyle='--', c='C2', alpha=0.8, label="Signal")
        # plt.plot(1e-6*(freqs), nois_sig, linestyle='--', c='C3', alpha=0.8, label="Noise")
        plt.xlabel("Frequency offset from LO (MHz)")
        plt.title(f"Spectroscopy Qubit {q_no}, SNR = {snr}")
        plt.grid()
        plt.legend()
        plt.show()
        plt.pause(1)

    I = I_handle.fetch_all()
    Q = Q_handle.fetch_all()
    qm.close()

    signal = I + 1j * Q
    phase = np.angle(signal)
    amplitude = np.abs(signal)
    sig = amplitude #phase


    # Ensure signal has peaks not troughs
    if (np.max(sig)-np.median(sig)) < (np.median(sig)-np.min(sig)):
        sig *= -1

    # Normalized signal, baseline oscillates around 0
    norm_sig = (sig - np.median(sig))/(np.max(sig) - np.median(sig))

    snr, [pure_sig, nois_sig] = S2N(norm_sig)
    # snr = np.mean(sig)/np.std(sig)

    plt.clf()
    plt.plot(1e-6*(freqs), norm_sig, label="Data") #, c='C1'
    # plt.plot(1e-6*(freqs), pure_sig, linestyle='--', c='C2', alpha=0.8, label="Signal")
    # plt.plot(1e-6*(freqs), nois_sig, linestyle='--', c='C3', alpha=0.8, label="Noise")
    plt.xlabel("Frequency offset from LO (MHz)")
    plt.title(f" Final Spectroscopy Qubit {q_no}, SNR = {snr}")
    plt.grid()
    plt.legend()
    plt.show()
    
    freq_array = 1e-9*(f_LO + freqs)
    return [norm_sig, pure_sig, nois_sig], freq_array

def broad_sweep(q_no, qmm, config, adc_mapping, q_LO, 
                n_samples, n_avgs, anharm, pi_12=False,
                q_amp=0.5, f_min_MHz=-350, f_max_MHz=350):
    
    # n_avgs = 5e2 but for good SNR - 100, bad SNR - 1000
    f_min = f_min_MHz * u.MHz
    f_max = f_max_MHz * u.MHz
    df = (f_max - f_min)/(n_samples-1)
    # df = df_MHz * u.MHz
    min_anharm_MHz = anharm-10
    min_zero2by2_det = (min_anharm_MHz/2) * u.MHz # zero2by2 detuning from qubit
    spacing = min_zero2by2_det/df

    spectro_params = (q_no, qmm, config, adc_mapping, q_LO, 
                      n_samples, n_avgs, q_amp, f_min, f_max, pi_12)

    [norm_sig, pure_sig, nois_sig], freqs = run_spectro(*spectro_params) # sig is returned normalized
    sig = norm_sig
    sig_peaks, peak_props = sgn.find_peaks(sig, prominence=0.6, width=1, rel_height=0.5, distance=spacing)

    # print(f'The signal peaks are found at the indices {sig_peaks}.')
    # print(f'The peak properties are {peak_props}.')

    widths = peak_props['widths']
    q_index = widths.tolist().index(max(widths))
    q_peak = sig_peaks[q_index]
    zero2by2_index = widths.tolist().index(min(widths))
    zero2by2_peak = sig_peaks[zero2by2_index]

    i = 0
    hlines_array = []
    for key, val in peak_props.items():
        # print(f'i is {i}')
        # print(f'key is {key}')
        # print(f'val is {val}')
        
        if i == 3:
            widths = val
        if i==4:
            hlines_array.append(val)
        if i>4:
            hlines_array.append(freqs[0] + val*df*1e-9) # df is in MHz
            
        i += 1
    hlines_list = (hlines_array[0], hlines_array[1], hlines_array[2])

    peaks = np.array([zero2by2_peak, q_peak])
    plt.figure()
    plt.title('Qubit Spectroscopy with 0-2/2 line')
    plt.plot(freqs, sig)
    plt.plot(freqs[sig_peaks], sig[sig_peaks], "x", color='C1')
    plt.plot(freqs[peaks], sig[peaks], "o", color='C2')
    # plt.vlines(x=freqs[sig_peaks], ymin=baseline, ymax=sig[sig_peaks], colors='C2')
    plt.hlines(y=hlines_list[0], xmin=hlines_list[1], xmax=hlines_list[2], colors='C3')
    plt.xlabel('Frequency sweep (GHz)')
    plt.ylabel('Normalized intensity (a.u.)')
    plt.grid()
    plt.show()

    q_freq = freqs[q_peak]
    zero2by2_freq = freqs[zero2by2_peak]
    q_fwhm = hlines_list[2][q_index] - hlines_list[1][q_index]  
    zero2by2_fwhm = hlines_list[2][zero2by2_index] - hlines_list[1][zero2by2_index]  

    print('##### Broad Qubit Sweep Outputs #####')
    print(f'The qubit frequency is {q_freq} GHz.')
    print(f'The FWHM of the qubit peak is {1e3 * q_fwhm} MHz.')
    print(f'The zero2by2 frequency is {zero2by2_freq} GHz.')
    # print(f'The FWHM of the zero2by2 peak is {1e3 * zero2by2_fwhm} MHz.')
    
    zero2by2_det = zero2by2_freq - q_freq
    anharm = 2*(zero2by2_det)
    print(f'The 0-2/2 line is detuned by {1e3 * zero2by2_det} MHz.')
    print(f'The anharmonicity of the qubit is {1e3 * anharm} MHz.')

    rerun_sweep = -zero2by2_det > (50*1e-3 + min_zero2by2_det*1e-9) # Not within 50 MHz of expected 0-2/2 detuning
    print(rerun_sweep)
    return [q_amp], [q_freq, q_fwhm], [zero2by2_freq, zero2by2_fwhm], [rerun_sweep]


def fine_qubit_sweep(q_no, qmm, config, adc_mapping, q_LO, 
                      n_samples, n_avgs, old_q_amp, 
                      old_q_freq, old_q_fwhm, 
                      pi_12=False, into_amp=1/3):
    
    f_LO = q_LO[f"{q_no}"]
    # n_avgs = 5e2 but for good SNR - 100, bad SNR - 1000
    q_offLO = old_q_freq*u.GHz - f_LO
    f_min = q_offLO - 3*old_q_fwhm*u.GHz
    f_max = q_offLO + 3*old_q_fwhm*u.GHz
    new_df = (f_max - f_min)/(n_samples-1)

    # new_df = old_df * into_df
    q_amp = old_q_amp * into_amp

    print('###### Fine Qubit Sweep Outputs')
    print(f'The qubit is offset from LO by {q_offLO}.\n'
          f'f_min, f_max, df are now {f_min, f_max, new_df}\n'
          f'The control amplitude is now {q_amp}.')

    spectro_params = (q_no, qmm, config, adc_mapping, q_LO, 
                      n_samples, n_avgs, q_amp, f_min, f_max,
                      pi_12)

    [norm_sig, pure_sig, nois_sig], freqs = run_spectro(*spectro_params) # sig is returned normalized
    sig = norm_sig

    
    sig -= min(sig)
    sig /= max(sig)

    # plt.plot(sig)

    sig_peaks, peak_props = sgn.find_peaks(sig, prominence=0.8, width=10, rel_height=0.5)

    widths = peak_props['widths']
    q_index = widths.tolist().index(max(widths))
    q_peak = sig_peaks[q_index]
    
    # print(f'The peak is found at the index: {q_peak}, at a frequency offset of {q_peak*new_df*1e-6} MHz.')

    for key, val in peak_props.items():
        vars()[f'q_{key[:-1]}'] = val[q_index]


    baseline = 0*sig[q_peak]

    q_freq = freqs[q_peak]
    q_fwhm = max(widths)*new_df*1e-9 # in GHz

    print('###########')
    print(f'The qubit frequency is {q_freq} GHz.')
    print(f'The FWHM of the qubit peak is {1e3 * q_fwhm} MHz.')
    
    plt.clf()
    plt.title('Qubit Fine Spectroscopy')
    plt.plot(freqs, sig)
    plt.plot(freqs[q_peak], sig[q_peak], "o", color="C1")
    plt.vlines(x=freqs[q_peak], ymin=baseline, ymax=sig[q_peak], colors='C2')
    plt.hlines(eval(f'q_width_height'), freqs[0] + eval(f'q_left_ip')*new_df*1e-9, freqs[0] + eval(f'q_right_ip')*new_df*1e-9, colors="C3")
    plt.show()    

    return q_freq, q_fwhm, q_amp


def auto_spectro(q_no, qmm, config, adc_mapping, q_LO, 
                 n_samples, n_avgs, anharm, pi_12=False):

    prog_props, q_props, zero2by2_props, rerun_sweep = broad_sweep(q_no, qmm, config, adc_mapping, q_LO, 
                                                                   n_samples, n_avgs, anharm, pi_12=pi_12)
    print("------------------\n"
         f"Broad qubit sweep rerun is needed? {rerun_sweep[0]}\n"
          "------------------")
    while rerun_sweep[0]:
        print("--------------------------------------------------------")
        print("Rerunning broad frequency sweep with more power and finer binning.")
        print("--------------------------------------------------------")

        old_amp = prog_props[0]
        # old_df_MHz = prog_props[1]*1e-6
        temp_amp = old_amp*5/4
        # temp_df_MHz = old_df_MHz/2

        prog_props, q_props, zero2by2_props, rerun_sweep  = broad_sweep(q_no, qmm, config, adc_mapping, q_LO, 
                                                                        n_samples, n_avgs, anharm, pi_12=pi_12, q_amp=temp_amp)

    old_q_amp = prog_props[0]
    old_q_freq = q_props[0]
    old_q_fwhm = q_props[1]
    old_zero2by2_freq = zero2by2_props[0]
    old_zero2by2_fwhm = zero2by2_props[1]

        
    q_fwhm = old_q_fwhm * u.GHz
    print(f'Do we need a rerun_sweep of fine sweep? A = {q_fwhm > 50 * u.kHz}')
    plt.figure()
    while q_fwhm > (50 * u.kHz):
    # for i in range(1,4):
        q_freq, q_fwhm, q_amp = fine_qubit_sweep(q_no, qmm, config, adc_mapping, q_LO, 
                                                 n_samples, n_avgs, old_q_amp, 
                                                 old_q_freq, old_q_fwhm, pi_12=pi_12)
        old_q_freq = q_freq
        old_q_fwhm = q_fwhm
        old_q_amp = q_amp
        q_fwhm = q_fwhm * u.GHz


    print('######## Full Spectro Output ########')
    print(f'The final qubit frequency is {q_freq} GHz.')
    print(f'The final FWHM of the qubit peak is {1e-6 * q_fwhm} MHz.')

    zero2by2_fwhm = old_zero2by2_fwhm * u.GHz
    zero2by2_freq = old_zero2by2_freq

    print(f'The final zero2by2 frequency is {zero2by2_freq} GHz.')
    # print(f'The final FWHM of the zero2by2 peak is {1e3 * zero2by2_fwhm} MHz.')
    print('##################################')

    zero2by2_det = zero2by2_freq - q_freq
    anharm = 2*zero2by2_det

    print(f'The 0-2/2 line is detuned by {1e3 * zero2by2_det} MHz.')
    print(f'The final anharmonicity of the qubit is {1e3 * anharm} MHz.')

    return [q_freq, q_fwhm], [zero2by2_freq, zero2by2_fwhm], anharm

