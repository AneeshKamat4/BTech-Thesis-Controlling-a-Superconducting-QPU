from qm.qua import *

cayley_table = np.int_(np.genfromtxt('Resources/c1_cayley_table.csv', delimiter=','))[1:, 1:]
inv_gates = [int(np.where(cayley_table[i, :] == 0)[0][0]) for i in range(24)]

def generate_sequence(max_circuit_depth, seed):

    cayley = declare(int, value=cayley_table.flatten().tolist())
    inv_list = declare(int, value=inv_gates)
    current_state = declare(int)
    step = declare(int)
    sequence = declare(int, size=max_circuit_depth+1)
    inv_gate = declare(int, size=max_circuit_depth+1)
    i = declare(int)
    rand = Random(seed=seed)

    assign(current_state, 0)
    with for_(i, 0, i < max_circuit_depth, i+1):
        assign(step, rand.rand_int(24))
        assign(current_state, cayley[current_state*24+step])
        assign(sequence[i], step)
        assign(inv_gate[i], inv_list[current_state])

    return sequence, inv_gate

def play_sequence(qe, sequence_list, depth, drag=False):
    i = declare(int)
    if drag:
        pre = 'd_'
    else:
        pre = ''
    with for_(i, 0, i <= depth, i+1):

        with switch_(sequence_list[i], unsafe=True):

            with case_(0):
                play('I', qe)
            with case_(1):
                play(f'{pre}X180', qe)
            with case_(2):
                play(f'{pre}Y180', qe)
            with case_(3):
                play(f'{pre}Y180', qe)
                play(f'{pre}X180', qe)
            with case_(4):
                play(f'{pre}X90', qe)
                play(f'{pre}Y90', qe)
            with case_(5):
                play(f'{pre}X90', qe)
                play(f'{pre}mY90', qe)
            with case_(6):
                play(f'{pre}mX90', qe)
                play(f'{pre}Y90', qe)
            with case_(7):
                play(f'{pre}mX90', qe)
                play(f'{pre}mY90', qe)
            with case_(8):
                play(f'{pre}Y90', qe)
                play(f'{pre}X90', qe)
            with case_(9):
                play(f'{pre}Y90', qe)
                play(f'{pre}mX90', qe)
            with case_(10):
                play(f'{pre}mY90', qe)
                play(f'{pre}X90', qe)
            with case_(11):
                play(f'{pre}mY90', qe)
                play(f'{pre}mX90', qe)
            with case_(12):
                play(f'{pre}X90', qe)
            with case_(13):
                play(f'{pre}mX90', qe)
            with case_(14):
                play(f'{pre}Y90', qe)
            with case_(15):
                play(f'{pre}mY90', qe)
            with case_(16):
                play(f'{pre}mX90', qe)
                play(f'{pre}Y90', qe)
                play(f'{pre}X90', qe)
            with case_(17):
                play(f'{pre}mX90', qe)
                play(f'{pre}mY90', qe)
                play(f'{pre}X90', qe)
            with case_(18):
                play(f'{pre}X180', qe)
                play(f'{pre}Y90', qe)
            with case_(19):
                play(f'{pre}X180', qe)
                play(f'{pre}mY90', qe)
            with case_(20):
                play(f'{pre}Y180', qe)
                play(f'{pre}X90', qe)
            with case_(21):
                play(f'{pre}Y180', qe)
                play(f'{pre}mX90', qe)
            with case_(22):
                play(f'{pre}X90', qe)
                play(f'{pre}Y90', qe)
                play(f'{pre}X90', qe)
            with case_(23):
                play(f'{pre}mX90', qe)
                play(f'{pre}Y90', qe)
                play(f'{pre}mX90', qe)

        # wait(4, qe)
