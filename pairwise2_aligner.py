from tkinter import *
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62

############### start to set env ################
len_str = 50
window = None
seq_a_label = None
seq_b_label = None
input_seq_a = None
input_seq_b = None

res_text = None

############### end setting env #################


def get_pairwise2_globalds_result(asequence, bsequence, gap_open_penalty=10, extension_penalty=0.5,
                                  matrx=blosum62):
    alignments = pairwise2.align.globalds(asequence.upper().replace(" ", ""), bsequence.upper().replace(" ", ""),
                                          matrx, -gap_open_penalty, -extension_penalty)
    alignments_result = pairwise2.format_alignment(*alignments[0])
    align_arr = alignments_result.split("\n")
    return align_arr[0], align_arr[1], align_arr[2], alignments_result


def clear_res_text(res_txt):
    if res_txt.get('1.0', END) != '':
        res_txt.delete('1.0', END)


def valid_input(seq_a, seq_b, res_txt):
    if seq_a.replace(' ', '') == '' and seq_b.replace(' ', '') == '':
        res_txt.insert(CURRENT, 'check [seq A] and [seq B]')
        return False
    elif seq_a.replace(' ', '') == '':
        res_txt.insert(CURRENT, 'check [seq A]')
        return False
    elif seq_b.replace(' ', '') == '':
        res_txt.insert(CURRENT, 'check [seq B]')
        return False
    else:
        return True


def setupGUI():
    global window
    global seq_a_label
    global seq_b_label
    global input_seq_a
    global input_seq_b
    global res_text

    window = Tk()
    window.title('pairwise2 aligner')

    seq_a_label = Label(window, text='seq A ', font='Courier 10 bold', relief=RAISED)
    seq_a_label.grid(row=0, column=0, padx=5, pady=5)
    input_seq_a = Entry(window, font='Terminal 10', width=50)
    input_seq_a.grid(row=0, column=1, padx=3, pady=3)
    # input_seq_a.bind('<Return>', aligner)

    seq_b_label = Label(window, text='seq B ', font='Courier 10 bold', relief=RAISED)
    seq_b_label.grid(row=1, column=0, padx=5, pady=5)
    input_seq_b = Entry(window, font='Terminal 10', width=50)
    input_seq_b.grid(row=1, column=1, padx=3, pady=3)
    input_seq_b.bind('<Return>', aligner)

    aligner_btn = Button(window, text='align', font='Courier 20 bold', command=aligner, height=10)
    aligner_btn.grid(rowspan=1, column=2, padx=3, pady=3)

    res_text = Text(window, font='Terminal 10', relief=RAISED, width=60, height=30)
    res_text.grid(row=2, columnspan=2, padx=5, pady=5)


def aligner():
    global input_seq_a
    global input_seq_b
    global res_text

    seq_a = input_seq_a.get()
    seq_b = input_seq_b.get()

    clear_res_text(res_text)

    if valid_input(seq_a, seq_b, res_text):
        seq_a_n, seq_n, seq_b_n, tot_re = get_pairwise2_globalds_result(seq_a.upper(), seq_b.upper())
        for i in range((len(seq_a_n) // len_str) + 1):
            res_text.insert(CURRENT, seq_a_n[i*len_str: (i + 1)*len_str] + '\n')
            res_text.insert(CURRENT, seq_n[i*len_str: (i + 1)*len_str] + '\n')
            res_text.insert(CURRENT, seq_b_n[i*len_str: (i + 1)*len_str] + '\n\n')


if __name__ == '__main__':
    setupGUI()
    window.mainloop()
