import numpy as np
from numpy import loadtxt
import sys
from scipy.signal.windows import gaussian

#Code to allow halting the program if simulating
class Halted(Exception) :
    def __init__(self): sys.tracebacklimit = 0
if hasattr(sys,'tracebacklimit'): del sys.tracebacklimit
# use as raise Halted()


def IQ_imbalance(g, phi):
    c = np.cos(phi)
    s = np.sin(phi)
    N = 1 / ((1-g**2)*(2*c**2-1))
    return [float(N * x) for x in [(1-g)*c, (1+g)*s, (1-g)*s, (1+g)*c]]

# DO NOT MODIFY ; CALL/INFORM Jay IF NECESSARY
def file_saver_(data, file_name = None, header_string="",suffix="",delimiter=",",master_folder="ExpName",init_path="E:/Experiments"):
    '''
    saves data as a file -> init_path/masterfolder/exp_type/date/file-time.csv

    Arguments : 
    'init_path' -> prefix path like "E:/Experiments"
    'exp_type' -> example "Spectroscopy"
    'master_folder" -> example "2qubit_RB_experiment"
    'header_string' -> string to passed on as header to store useful parameters about exp.
    '''
    from datetime import datetime
    import os

    if file_name is None :
        import inspect
        file_name = inspect.stack()[1].filename
    now = datetime.now()
    current_date = now.strftime("%y-%m-%d")
    current_time = now.strftime("%H-%M-%S")
    #Create relevant Folders
    master_path = init_path+"/"+master_folder
    exp_type = file_name.split("\\")[-1][:-3]
    os.makedirs(master_path+"/"+exp_type+"/"+current_date,exist_ok=True)
    final_path = master_path+"/"+exp_type+"/"+current_date

    #Create file and save DATA
    filename = file_name.split("\\")[-1][:-3]+f"-{suffix}-"+current_time+".csv"
    np.savetxt(final_path+"/"+filename,data,delimiter=delimiter,header=header_string)

    return print("File saved as: ",final_path+"/"+filename)

def integration_weight_update(b,config,integ_len_clk):
    w_plus_cos = [(np.cos(b),  integ_len_clk*4)]
    w_minus_cos = [(-np.cos(b),  integ_len_clk*4)]
    w_plus_sin = [(np.sin(b),  integ_len_clk*4)]
    w_minus_sin = [(-np.sin(b),  integ_len_clk*4)]
    config['integration_weights']['integW_cos']['cosine'] = w_plus_cos
    config['integration_weights']['integW_cos']['sine'] = w_minus_sin
    config['integration_weights']['integW_sin']['cosine'] = w_plus_sin
    config['integration_weights']['integW_sin']['sine'] = w_plus_cos
    config['integration_weights']['integW_minus_sin']['cosine'] = w_minus_sin
    config['integration_weights']['integW_minus_sin']['sine'] = w_minus_cos
    return 1

def _len_guards(M):
    """Handle small or incorrect window lengths"""
    if int(M) != M or M < 0:
        raise ValueError('Window length M must be a non-negative integer')
    return M <= 1

def _extend(M, sym):
    """Extend window by 1 sample if needed for DFT-even symmetry"""
    if not sym:
        return M + 1, True
    else:
        return M, False

def _truncate(w, needed):
    """Truncate window by 1 sample if needed for DFT-even symmetry"""
    if needed:
        return w[:-1]
    else:
        return w


def gauss_window(M, std, sym=True):
    r"""Return a Gaussian window.

    Parameters
    ----------
    M : int
        Number of points in the output window. If zero, an empty array
        is returned. An exception is thrown when it is negative.
    std : float
        The standard deviation, sigma.
    sym : bool, optional
        When True (default), generates a symmetric window, for use in filter
        design.
        When False, generates a periodic window, for use in spectral analysis.

    Returns
    -------
    w : ndarray
        The window, with the maximum value normalized to 1 (though the value 1
        does not appear if `M` is even and `sym` is True).

    Notes
    -----
    The Gaussian window is defined as

    .. math::  w(n) = e^{ -\frac{1}{2}\left(\frac{n}{\sigma}\right)^2 }

    """
    if _len_guards(M):
        return np.ones(M)
    M, needs_trunc = _extend(M, sym)

    n = np.arange(0, M) - (M - 1.0) / 2.0
    sig2 = 2 * std * std
    w = np.exp(-n ** 2 / sig2)

    return _truncate(w, needs_trunc)

def gauss_der_window(M, std, sym=True):
    r"""Return a Gaussian window.

    Parameters
    ----------
    M : int
        Number of points in the output window. If zero, an empty array
        is returned. An exception is thrown when it is negative.
    std : float
        The standard deviation, sigma.
    sym : bool, optional
        When True (default), generates a symmetric window, for use in filter
        design.
        When False, generates a periodic window, for use in spectral analysis.

    Returns
    -------
    w : ndarray
        The window, with the maximum value normalized to 1 (though the value 1
        does not appear if `M` is even and `sym` is True).

    Notes
    -----
    The Gaussian window is defined as

    .. math::  w(n) = e^{ -\frac{1}{2}\left(\frac{n}{\sigma}\right)^2 }

    """
    if _len_guards(M):
        return np.ones(M)
    M, needs_trunc = _extend(M, sym)

    n = np.arange(0, M) - (M - 1.0) / 2.0
    sig2 = 2 * std * std
    w = -2 * (1e9 * n / sig2) * np.exp(-n ** 2 / sig2) #1e9 because pulse length is M in ns

    return _truncate(w, needs_trunc)


def grft_pulse(pi_len_ns, pi_rise_ns):
    risefall = [float(arg) for arg in gauss_window(2 * pi_rise_ns, 2 * pi_rise_ns // 6)]
    pulse = []
    for i in range(pi_len_ns):
        if i < pi_rise_ns:
            pulse.append(risefall[i])

        elif i >= pi_len_ns - pi_rise_ns:
            pulse.append(risefall[i - pi_len_ns + 2 * pi_rise_ns])
            
        else:
            pulse.append(1)

    return np.array(pulse)


def grft_arr_gen(arguments,scales=[1.0]):
    scaling = np.prod(scales) * 0.4
    slist = np.array([float(arg) for arg in scaling * grft_pulse(*arguments)], dtype=complex)
    return slist

def grft_der_pulse(pi_len_ns, pi_rise_ns):
    risefall = [float(arg) for arg in gauss_der_window(2 * pi_rise_ns, 2 * pi_rise_ns // 6)]
    pulse = []
    for i in range(pi_len_ns):
        if i < pi_rise_ns:
            pulse.append(risefall[i])

        elif i >= pi_len_ns - pi_rise_ns:
            pulse.append(risefall[i - pi_len_ns + 2 * pi_rise_ns])

        else:
            pulse.append(0)

    return np.array(pulse)


def grft_der_arr_gen(arguments, scales=[1.0]):
    scaling = np.prod(scales) * 0.4
    slist = np.array([float(arg) for arg in scaling * grft_der_pulse(*arguments)], dtype=complex)
    return slist


def gauss_arr_gen(arguments,scales=[1.0]):
    scaling = np.prod(scales) * 0.4
    slist = [float(arg) for arg in scaling * gaussian(*arguments)]
    return slist

def rise_arr(rise_ns):

    risefall = [float(0.4*arg) for arg in gaussian(2 * rise_ns, 2 * rise_ns // 6)]
    l = int(len(risefall)/2)
    rise = risefall[:l]

    return rise


def fall_arr(rise_ns):
    risefall = [float(0.4*arg) for arg in gaussian(2 * rise_ns, 2 * rise_ns // 6)]
    l = int(len(risefall) / 2)
    fall = risefall[l:]

    return fall

def drag_grft_pulse_waveforms(
    amplitude,
    length,
    rise,
    anharmonicity,
    alpha,
    detuning=0.0,
    **kwargs
):
    """
    Creates Gaussian Rise Flat Top (GRFT) based DRAG waveforms that compensate for the leakage and for the AC stark shift.

    These DRAG waveforms has been implemented following the next Refs.:
    Chen et al. PRL, 116, 020501 (2016)
    https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.116.020501
    and Chen's thesis
    https://web.physics.ucsb.edu/~martinisgroup/theses/Chen2018.pdf

    :param float amplitude: The amplitude in volts. (Will be scaled down by 0.4)
    :param int length: The pulse length in ns.
    :param int rise: The rise-fall pulse length in ns
    :param float alpha: The DRAG coefficient.
    :param float anharmonicity: f_21 - f_10 - The differences in energy between the 2-1 and the 1-0 energy levels, in Hz.
    :param float detuning: The frequency shift to correct for AC stark shift, in Hz.
    :return: Returns a tuple of two lists. The first list is the I waveform (real part) and the second is the
        Q waveform (imaginary part)
    """
    if alpha != 0 and anharmonicity == 0:
        raise Exception("Cannot create a DRAG pulse with `anharmonicity=0`")

    t = np.arange(length, dtype=int)  # An array of size pulse length in ns

    grft_wave = grft_arr_gen((length, rise), [amplitude])
    grft_der_wave = grft_der_arr_gen((length, rise), [amplitude])

    z = grft_wave + 1j * 0
    if alpha != 0:
        # The complex DRAG envelope:
        z += (1j * grft_der_wave * (alpha / (2 * np.pi * anharmonicity - 2 * np.pi * (detuning * 1e6))))
        # The complex detuned DRAG envelope:
        z *= np.exp(1j * 2 * np.pi * (detuning * 1e6) * t * 1e-9)
    I_wf = z.real.tolist()  # The `I` component is the real part of the waveform
    Q_wf = (z.imag.tolist())  # The `Q` component is the imaginary part of the waveform
    return [I_wf, Q_wf]

def drag_gaussian_pulse_waveforms(
    amplitude,
    length,
    sigma,
    alpha,
    anharmonicity,
    detuning=0.0,
    subtracted=True,
    **kwargs
):
    """
    Creates Gaussian based DRAG waveforms that compensate for the leakage and for the AC stark shift.

    These DRAG waveforms has been implemented following the next Refs.:
    Chen et al. PRL, 116, 020501 (2016)
    https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.116.020501
    and Chen's thesis
    https://web.physics.ucsb.edu/~martinisgroup/theses/Chen2018.pdf

    :param float amplitude: The amplitude in volts.
    :param int length: The pulse length in ns.
    :param float sigma: The gaussian standard deviation.
    :param float alpha: The DRAG coefficient.
    :param float anharmonicity: f_21 - f_10 - The differences in energy between the 2-1 and the 1-0 energy levels, in Hz.
    :param float detuning: The frequency shift to correct for AC stark shift, in Hz.
    :param bool subtracted: If true, returns a subtracted Gaussian, such that the first and last points will be at 0
        volts. This reduces high-frequency components due to the initial and final points offset. Default is true.
    :return: Returns a tuple of two lists. The first list is the I waveform (real part) and the second is the
        Q waveform (imaginary part)
    """

    if alpha != 0 and anharmonicity == 0:
        raise Exception("Cannot create a DRAG pulse with `anharmonicity=0`")
    t = np.arange(length, dtype=int)  # An array of size pulse length in ns
    center = (length - 1) / 2
    gauss_wave = amplitude * np.exp(
        -((t - center) ** 2) / (2 * sigma**2)
    )  # The gaussian function
    gauss_der_wave = (
        amplitude
        * (-2 * 1e9 * (t - center) / (2 * sigma**2))
        * np.exp(-((t - center) ** 2) / (2 * sigma**2))
    )  # The derivative of gaussian
    if subtracted:
        gauss_wave = gauss_wave - gauss_wave[-1]  # subtracted gaussian
    z = gauss_wave + 1j * 0
    if alpha != 0:
        # The complex DRAG envelope:
        z += (
            1j
            * gauss_der_wave
            * (alpha / (2 * np.pi * anharmonicity - 2 * np.pi * detuning))
        )
        # The complex detuned DRAG envelope:
        z *= np.exp(1j * 2 * np.pi * detuning * t * 1e-9)
    I_wf = z.real.tolist()  # The `I` component is the real part of the waveform
    Q_wf = (
        z.imag.tolist()
    )  # The `Q` component is the imaginary part of the waveform
    return I_wf, Q_wf

def c3_waveforms(amplitude, I_file_name, Q_file_name):
    '''
    Reads the I and Q waveform samples from their respective .csv files
    '''
    Ihold_wf = loadtxt(I_file_name, delimiter=',')
    Qhold_wf = loadtxt(Q_file_name, delimiter=',')

    I_wf = amplitude*Ihold_wf
    Q_wf = amplitude*Qhold_wf

    return [I_wf, Q_wf]            


