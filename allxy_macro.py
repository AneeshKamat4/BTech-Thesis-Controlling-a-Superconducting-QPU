from qm.qua import *
from macros import *
from configuration_4qubitsv4 import *
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt


# Rabi Analysis
def rabi_fit(t, A, f, d, p, c):
    return A * np.exp(-t/d) * np.sin(2*np.pi*f*t + p) + c


# Time Rabi to guess initial amplitudes
def initialize_Rabi(q_no, qmm, n_avg, t_min, t_max, dt, t_list_clk, wait_init, qe, rr, out, rabi_pi_amp, show_Rabi=False):
    with program() as rabi:
        n = declare(int)
        t = declare(int)

        Irabi = declare(fixed)
        Qrabi = declare(fixed)
        Irabi_st = declare_stream()
        Qrabi_st = declare_stream()

        with for_(n, 0, n < n_avg, n + 1):
            with for_(t, t_min, t < t_max + 0.1, t + dt):
                wait(wait_init)
                play("grft"*amp(rabi_pi_amp), qe, t)
                measure_macro(qe, rr, out, 
                            Irabi, Qrabi, pi_12=False)
                save(Irabi, Irabi_st)
                save(Qrabi, Qrabi_st)

        with stream_processing():
            Irabi_st.buffer(len(t_list_clk)).average().save('Irabi')
            Qrabi_st.buffer(len(t_list_clk)).average().save('Qrabi')

    # ------------------------
    #   Execute on the OPX
    # ------------------------

    qm1 = qmm.open_qm(config)
    job = qm1.execute(rabi)
    res_handles = job.result_handles
    # ------
    I_handle = job.result_handles.get("Irabi")
    Q_handle = job.result_handles.get("Qrabi")

    t_list = 4 * t_list_clk
    
    if show_Rabi is True:
        plt.figure(1)
        plt.title("Initializing Time Rabi")
        I_handle.wait_for_values(1)
        Q_handle.wait_for_values(1)

        while res_handles.is_processing():
            Irabi = I_handle.fetch_all()
            Qrabi = Q_handle.fetch_all()
            sig = Irabi + 1j * Qrabi

            plt.clf()
            plt.plot(t_list, Qrabi, marker='.', label="Q")
            plt.plot(t_list, Irabi, marker='.', label="I")
            plt.xlabel("Time (ns)")
            plt.ylabel("Rabi Amplitude")
            plt.title(f"Time Rabi on Qubit {q_no}")
            plt.legend()
            plt.grid()
            plt.pause(0.25)
        plt.close(1)

    Irabi = job.result_handles.get("Irabi").fetch_all()
    Qrabi = job.result_handles.get("Qrabi").fetch_all()

    qm1.close()

    # ############
    # # analysis #
    # ############
    sig = Irabi + 1j*Qrabi

    pars, cov = curve_fit(f=rabi_fit, xdata=t_list, ydata=Irabi, p0=[3e-3,0.01,100,0, 1e-5], bounds=(-np.inf, np.inf), maxfev=2000)
    # init_pars = p0 = [A, f, d, p, c]

    rabi_pi_len = np.round(0.5/pars[1],3)
    # rabi_pi_amp = pars[0]

    print('######################### \n### Fitted Parameters ### \n######################### ')
    print(f"Rabi frequency = {np.round(1e3*pars[1],2)} MHz")
    print(f"Pi pulse = {rabi_pi_len} ns")
    print(f"Rabi amplitude = {np.abs(pars[0])}")
    print(f"Rabi decay constant = {pars[2]*1e-3} us")

    init_pi_amp = rabi_pi_amp*rabi_pi_len/pi_len_ns[f"{q_no}"]
    init_pi_amp = min(0.9, max(0, init_pi_amp))
    initial_amps = np.array([init_pi_amp]*2 + [init_pi_amp/2]*2)
    # initial_amps = np.fmin(np.array([0.9]*4), np.fmax(np.array([0]*4), initial_amps))

    return initial_amps

# ------------------------------
# All XY sequences. The sequence names must match corresponding operation in the config
sequence = [  # based on https://rsl.yale.edu/sites/default/files/physreva.82.pdf-optimized_driving_0.pdf
    ("I", "I"),
    ("X180", "X180"),
    ("Y180", "Y180"),
    ("X180", "Y180"),
    ("Y180", "X180"),
    ("X90", "I"),
    ("Y90", "I"),
    ("X90", "Y90"),
    ("Y90", "X90"),
    ("X90", "Y180"),
    ("Y90", "X180"),
    ("X180", "Y90"),
    ("Y180", "X90"),
    ("X90", "X180"),
    ("X180", "X90"),
    ("Y90", "Y180"),
    ("Y180", "Y90"),
    ("X180", "I"),
    ("Y180", "I"),
    ("X90", "X90"),
    ("Y90", "Y90"),
]

# Pulses are all pair-wise combinations of I, X180, Y180, X90, Y90
label_seq = [  # for aesthetic purposes
    "IdId",
    "XpXp",
    "YpYp",
    "XpYp",
    "YpXp",
    "X9Id",
    "Y9Id",
    "X9Y9",
    "Y9X9",
    "X9Yp",
    "Y9Xp",
    "XpY9",
    "YpX9",
    "X9Xp",
    "XpX9",
    "Y9Yp",
    "YpY9",
    "XpId",
    "YpId",
    "X9X9",
    "Y9Y9"
]

# All XY macro generating the pulse sequences from a python list.
def run_allXY(pulses, amp_scaling, wait_init, qe, rr, out):
    """
    Generate a QUA sequence based on the two operations written in pulses. Used to generate the all XY program.
    **Example:** I, Q = allXY(['I', 'Y90'])

    :param pulses: tuple containing a particular set of operations to play. The pulse names must match corresponding
        operations in the config except for the identity operation that must be called 'I'.
    :return: two QUA variables for the 'I' and 'Q' quadratures measured after the sequence.
    """
    I_xy = declare(fixed)
    Q_xy = declare(fixed)
    
    wait(wait_init, qe)
    play(pulses[0] * amp(amp_scaling[0]), qe)
    play(pulses[1] * amp(amp_scaling[1]), qe)
    measure_macro(qe, rr, out, I_xy, Q_xy, pi_12=False)

    return I_xy, Q_xy


###################
# The QUA program #
###################
def qua_allxy(opt_amp_dict, n_avg, t_min, t_max, dt, t_list_clk, wait_init, qe, rr, out, rabi_pi_amp):
    
    with program() as ALLXY:
        n = declare(int)
        i = declare(int)
        t = declare(int)

        Irabi = declare(fixed)
        Qrabi = declare(fixed)
        Irabi_st = declare_stream()
        Qrabi_st = declare_stream()

        Ixy_st = [declare_stream() for _ in range(21)]
        Qxy_st = [declare_stream() for _ in range(21)]
        
        # sigmaz = declare(float)
        # z_st = [declare_stream() for _ in range(21)]
        # z_exp = declare_stream()

        with for_(n, 0, n < n_avg, n + 1):
            with for_(t, t_min, t < t_max + 0.1, t + dt):
                # Find Rabi levels
                wait(wait_init)
                play("grft"*amp(rabi_pi_amp), qe, t)
                measure_macro(qe, rr, out, 
                              Irabi, Qrabi, pi_12=False)
                
                save(Irabi, Irabi_st)
                save(Qrabi, Qrabi_st)

            for i in range(len(sequence)):
                # Scale the pulse amplitudes in the sequences 
                pulses = sequence[i]
                amp_scaling = [opt_amp_dict[pulses[0]], opt_amp_dict[pulses[1]]]
                Ixy, Qxy = run_allXY(pulses, amp_scaling, wait_init, qe, rr, out)
                save(Ixy, Ixy_st[i])
                save(Qxy, Qxy_st[i])

        with stream_processing():
            Irabi_st.buffer(len(t_list_clk)).average().save('Irabi')
            Qrabi_st.buffer(len(t_list_clk)).average().save('Qrabi')

            for j in range(21):
                Ixy_st[j].average().save(f"Ixy_{j}")
                Qxy_st[j].average().save(f"Qxy_{j}")
    return ALLXY


def opt_allxy(opt_amps, q_no, qmm, n_avg, t_min, t_max, dt, t_list_clk, wait_init, qe, rr, out, rabi_pi_amp, show_Rabi=False, print_Rabi=False):
    
    opt_amp_dict = {"I": 1,
                    "X180": opt_amps[0],
                    "Y180": opt_amps[1],
                    "X90": opt_amps[2],
                    "Y90": opt_amps[3]}
    
    qm = qmm.open_qm(config)
    job = qm.execute(qua_allxy(opt_amp_dict, n_avg, t_min, t_max, dt, t_list_clk, wait_init, qe, rr, out, rabi_pi_amp))
    res_handles = job.result_handles
    # -------
    I_handle = job.result_handles.get("Irabi")
    Q_handle = job.result_handles.get("Qrabi")

    t_list = 4 * t_list_clk
    # print("Hiiiiiii")
    I_handle.wait_for_values(1)
    Q_handle.wait_for_values(1)
    if show_Rabi is True:
        plt.figure(2)
        plt.title("OptALLXY Rabi")

        while res_handles.is_processing():
            Irabi = I_handle.fetch_all()
            Qrabi = Q_handle.fetch_all()
            sig = Irabi + 1j * Qrabi

            plt.clf()
            plt.plot(t_list, Qrabi, marker='.', label="Rabi - Q")
            plt.plot(t_list, Irabi, marker='.', label="Rabi - I")
            plt.xlabel("Time (ns)")
            plt.ylabel("Rabi Amplitude")
            plt.title(f"OptALLXY Time Rabi on Qubit {q_no}")
            plt.legend()
            plt.grid()
            plt.pause(0.25)
        plt.close(2)

    Irabi = job.result_handles.get("Irabi").fetch_all()
    Qrabi = job.result_handles.get("Qrabi").fetch_all()
    # -------
    # job.result_handles.wait_for_all_values()

    Ixy = []
    Qxy = []

    for j in range(21):
        Ixy.append(job.result_handles.get(f"Ixy_{j}").fetch_all())
        Qxy.append(job.result_handles.get(f"Qxy_{j}").fetch_all())

    I_xy = np.array(Ixy)
    Q_xy = np.array(Qxy)

    qm.close()

    pars, cov = curve_fit(f=rabi_fit, xdata=t_list, ydata=Irabi, p0=[3e-3,0.01,100,0, 1e-5], bounds=(-np.inf, np.inf), maxfev=2000)
    # init_pars = p0 = [A, f, d, p, c]
    
    I_eq = pars[4]
    rabi_amp = abs(pars[0])
    I_g = I_eq - rabi_amp
    I_e = I_eq + rabi_amp

    if print_Rabi is True:
        print('######################### \n### Fitted Parameters ### \n######################### ')
        print(f"Rabi frequency = {np.round(1e3*pars[1],2)} MHz")
        print(f"Pi pulse = {np.round(0.5/pars[1],3)} ns")
        print(f"Rabi amplitude = {pars[0]}")
        print(f"Rabi decay constant = {pars[2]*1e-3} us")
        print(f"Rabi DC offset = {I_eq}")
        print(f"In arbitrary units \n Ground state [-1] = {I_g} \n Equator [0] = {I_eq} \n Excited state [+1] = {I_e}")

    z_xy = (I_xy - I_eq)/rabi_amp # Projecting ALLXY to -1:1 Z-space


    '''fig, ax = plt.subplots(nrows=2, ncols=1)

    ax[0].plot(I_xy, marker='o', linestyle=':', markersize=5, c='green')
    ax[0].axhline(y=I_g, color='r', linewidth=0.5)
    ax[0].axhline(y=I_eq, color='r', linewidth=0.5)
    ax[0].axhline(y=I_e, color='r', linewidth=0.5)
    ax[0].set_xticks(ticks=range(21), labels=[str(el) for el in label_seq], rotation=90)

    ax[1].plot(z_xy, marker='o', linestyle=':', markersize=5, c='blue')
    ax[1].axhline(y=-1, color='r', linewidth=0.5)
    ax[1].axhline(y=0, color='r', linewidth=0.5)
    ax[1].axhline(y=1, color='r', linewidth=0.5)
    ax[1].set_xticks(ticks=range(21), labels=[str(el) for el in label_seq], rotation=90)

    ax[0].set_title('I quadrature [a.u.]')
    ax[1].set_title('Z projection')

    fig.suptitle(f"All XY Calibration on Qubit {q_no}")'''

    plt.figure(11)
    plt.plot(z_xy, marker='o', linestyle=':', markersize=5, c='blue')
    plt.axhline(y=-1, color='r', linewidth=0.5)
    plt.axhline(y=0, color='r', linewidth=0.5)
    plt.axhline(y=1, color='r', linewidth=0.5)
    plt.xticks(ticks=range(21), labels=[str(el) for el in label_seq], rotation=90)

    plt.ylabel('Z projection')
    plt.title(f"All XY Calibration on Qubit {q_no}")
    plt.show()

    '''plt.figure(12)
    plt.plot(I_xy, marker='o', linestyle=':', markersize=5, c='blue')
    plt.axhline(y=I_g, color='r', linewidth=0.5)
    plt.axhline(y=I_eq, color='r', linewidth=0.5)
    plt.axhline(y=I_e, color='r', linewidth=0.5)
    plt.xticks(ticks=range(21), labels=[str(el) for el in label_seq], rotation=90)

    plt.ylabel('I quadrature [a.u.]')
    plt.title(f"All XY Calibration on Qubit {q_no}")
    
    plt.tight_layout()
    plt.show()'''

    target = np.array([-1] * 5 + [0] * 12 + [1] * 4) # Expected sigma z averaged over all pulse sequences
    # target = np.array([-1] * 5 + [0] * 12 + [1] * 4) # Expected sigma z averaged over all pulse sequences
    # mse = np.sum((target[i]-z_xy[i])**2 for i in len(target))/len(target)
    mse = np.mean((target-z_xy)**2)
    # fit = np.linalg.norm(target - z_xy)
    # print(f'The fit is {fit}')
    print(f'The Mean Squared Error (MSE) is {mse}')
    return mse

