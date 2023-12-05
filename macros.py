from qm.qua import *
from configuration_4qubitsv4 import cr_amp, cr_len_ns


def play_flat_top(qe, a, t = 5):

    play("rise"*amp(a), qe)
    play("const"*amp(a), qe, duration=t)
    play("fall"*amp(a), qe)


def Hadamard(qe):

    # frame_rotation_2pi(-0.5, qe)
    play("Y90", qe)
    play("X180", qe)
    play("Y180", qe)

def cooldown(time=250000, active_reset=False, qe=None, qe_12=None, rr=None, out=None, I=None, Q=None, pi_12=False, dem=None):

    if active_reset:

        measure_macro(qe, rr, out, I, Q, pi_12=pi_12)
        if pi_12:
            align(rr, qe_12)
            play("X180", qe_12, condition=I > dem)
            align(qe_12, qe)
            play("X180", qe, condition=I > dem)

        else:
            align(rr, qe)
            play("X180", qe, condition=I > dem)

    else:
        wait(time)

def measure_macro(qe, rr, out, I, Q, pi_12=False):

    wait_rr = 16
    qe12 = f"q12_{rr[-1]}"

    if pi_12:
        align(qe, qe12)
        wait(4, qe12)
        play("X180", qe12)
        align(qe12, rr)
        wait(wait_rr, rr)

    else:
        align(qe, rr)
        wait(wait_rr, rr)

    measure("readout", rr, None,
            demod.full("integW_cos", I, out),
            demod.full("integW_minus_sin", Q, out))

def ZXby4(qe_cr, qe_c, a=1.0, t=28):

    tby2 = t // 2
    play_flat_top(qe_cr, a, tby2)
    align(qe_c, qe_cr)


def ZXby2_echo_noAC(qe_cr, qe_c, qe_t):

    a = cr_amp[qe_cr]
    t = cr_len_ns[qe_cr]//4

    tby2 = t // 2

    align(qe_cr, qe_c, qe_t)
    play_flat_top(qe_cr, a, tby2)
    align(qe_c, qe_cr)
    wait(4, qe_c)
    play("X180", qe_c)
    align(qe_c, qe_cr)
    play_flat_top(qe_cr, -a, tby2)
    align(qe_c, qe_cr)
    wait(4, qe_c)
    play("X180", qe_c)


def CNOT_macro(qe_c, qe_t):

    c, t = qe_c[-1], qe_t[-1]

    if qe_c == "q2" or qe_c == "q4":

        qe_cr = f"cr_c{c}t{t}"
        ZXby2_echo_noAC(qe_cr, qe_c, qe_t)
        align(qe_cr, qe_t)
        wait(4, qe_t)
        play("X90", qe_t)  #fixed for originally mX90
        wait(4, qe_t)
        align(qe_c, qe_t)
        play("X90", qe_c)
        play("Y90", qe_c)
        play("mX90", qe_c)

    if qe_c == "q1" or qe_c == "q3":

        qe_cr = f"cr_c{t}t{c}"
        align(qe_cr, qe_c, qe_t)
        Hadamard(qe_c)
        Hadamard(qe_t)
        wait(4, qe_c)
        align(qe_cr, qe_c, qe_t)
        CNOT_macro(qe_t, qe_c)
        wait(4, qe_c)
        align(qe_cr, qe_c, qe_t)
        Hadamard(qe_c)
        Hadamard(qe_t)
