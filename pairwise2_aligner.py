from tkinter import *
from tkinter import ttk
from tkinter import messagebox
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62

############### start to set env ################
len_str = 50
window = None
input_seq_a = None
input_seq_b = None

align_radio = None
align_opt = None
gap_opn_pnlty = None
extensn_pnlty = None

res_text = None

############### end setting env #################


def get_pairwise2_globalds_result(asequence, bsequence, gap_open_penalty=10.0, extension_penalty=0.5,
                                  matrx=blosum62):
    alignments = pairwise2.align.globalds(asequence.upper().replace(" ", ""), bsequence.upper().replace(" ", ""),
                                          matrx, -gap_open_penalty, -extension_penalty)
    alignments_result = pairwise2.format_alignment(*alignments[0])
    align_arr = alignments_result.split("\n")
    return align_arr[0], align_arr[1], align_arr[2], alignments_result


def get_pairwise2_localds_result(asequence, bsequence, gap_open_penalty=10.0, extension_penalty=1.0, matrx=blosum62):
    alignments = pairwise2.align.localds(asequence.upper().replace(" ", ""), bsequence.upper().replace(" ", ""),
                                          matrx, -gap_open_penalty, -extension_penalty)
    alignments_result = pairwise2.format_alignment(*alignments[0])
    align_arr = alignments_result.split("\n")
    return ''.join([i for i in align_arr[0] if not i.isdigit()]), ''.join(
        [i for i in align_arr[1] if not i.isdigit()]), ''.join(
        [i for i in align_arr[2] if not i.isdigit()]), alignments_result


def set_align_opt(val):
    global align_opt
    global gap_opn_pnlty
    global extensn_pnlty

    gap_opn_pnlty.delete(0, 'end')
    gap_opn_pnlty.insert(0, 10.0)

    extensn_pnlty.delete(0, 'end')
    if val == 'local':
        extensn_pnlty.insert(0, 1.0)
    else:
        extensn_pnlty.insert(0, 0.5)

    align_opt = val


def get_idx_seq(input_seq, st_idx, ignr_char='-'):
    result_seq = ''
    for tmp_ch in input_seq:
        if tmp_ch == ignr_char:
            result_seq += ' '
        else:
            result_seq += str(st_idx)
            st_idx += 1
            if st_idx >= 10:
                st_idx = 0
    return result_seq, st_idx


def clear_res_text(res_txt):
    if res_txt.get('1.0', END) != '':
        res_txt.delete('1.0', END)


def valid_input(seq_a, seq_b, res_txt):
    if seq_a.replace(' ', '') == '' and seq_b.replace(' ', '') == '':
        messagebox.showerror('error', 'check [seq A] and [seq B]')
        return False
    elif seq_a.replace(' ', '') == '':
        messagebox.showerror('error', 'check [seq A]')
        return False
    elif seq_b.replace(' ', '') == '':
        messagebox.showerror('error', 'check [seq B]')
        return False
    else:
        return True


def reset():
    global input_seq_a
    global input_seq_b
    global gap_opn_pnlty
    global extensn_pnlty
    global align_radio
    global res_text

    input_seq_a.delete(0, 'end')
    input_seq_b.delete(0, 'end')
    gap_opn_pnlty.delete(0, 'end')
    extensn_pnlty.delete(0, 'end')
    clear_res_text(res_text)

    gap_opn_pnlty.insert(0, 10.0)
    extensn_pnlty.insert(0, 0.5)
    align_radio.set('global')


def setupGUI():
    global window
    global input_seq_a
    global input_seq_b
    global align_radio
    global gap_opn_pnlty
    global extensn_pnlty
    global res_text

    window = Tk()
    window.title('pairwise2 aligner')

    seq_a_label = Label(window, text='seq A ', font='Courier 10 bold', relief=RAISED)
    seq_a_label.grid(row=1, column=0, padx=5, pady=5)
    input_seq_a = Entry(window, font='Terminal 10', width=50)
    input_seq_a.grid(row=1, column=1, padx=3, pady=3)
    # input_seq_a.bind('<Return>', aligner)

    seq_b_label = Label(window, text='seq B ', font='Courier 10 bold', relief=RAISED)
    seq_b_label.grid(row=2, column=0, padx=5, pady=5)
    input_seq_b = Entry(window, font='Terminal 10', width=50)
    input_seq_b.grid(row=2, column=1, padx=3, pady=3)
    input_seq_b.bind('<Return>', aligner)

    align_opt_label1 = Label(window, relief=RAISED)
    align_opt_label1.grid(row=0, column=2, padx=5, pady=5)
    align_radio = StringVar()
    global_opt = ttk.Radiobutton(align_opt_label1, text='global', value='global', variable=align_radio, command=lambda: set_align_opt('global'))
    global_opt.grid(row=0, column=0)
    local_opt = ttk.Radiobutton(align_opt_label1, text=' local ', value='local', variable=align_radio, command=lambda: set_align_opt('local'))
    local_opt.grid(row=1, column=0)
    align_radio.set('global')

    align_opt_label2 = Label(window, relief=RAISED)
    align_opt_label2.grid(row=0, columnspan=2, padx=5, pady=5)
    gap_opn_label = Label(align_opt_label2, text='gap open penalty :', font='Courier 10 bold', relief=FLAT)
    gap_opn_label.grid(row=0, column=0, padx=5, pady=5)
    gap_opn_pnlty = Entry(align_opt_label2, font='Terminal 10', width=5)
    gap_opn_pnlty.grid(row=0, column=1, padx=3, pady=3)
    gap_opn_pnlty.insert(0, 10.0)
    extensn_label = Label(align_opt_label2, text=' extension penalty :', font='Courier 10 bold', relief=FLAT)
    extensn_label.grid(row=0, column=2, padx=5, pady=5)
    extensn_pnlty = Entry(align_opt_label2, font='Terminal 10', width=5)
    extensn_pnlty.grid(row=0, column=3, padx=3, pady=3)
    extensn_pnlty.insert(0, 0.5)

    aligner_btn = Button(window, text='align', font='Courier 20 bold', command=aligner, height=10)
    aligner_btn.grid(rowspan=1, column=2, padx=3, pady=3)

    reset_btn = Button(window, text='reset', font='Courier 10 bold', fg='red', command=reset)
    reset_btn.grid(row=2, column=2, padx=3, pady=3)

    res_text = Text(window, font='Terminal 10', relief=RAISED, width=65, height=30)
    res_text.grid(row=3, columnspan=2, padx=5, pady=5)


def aligner():
    global input_seq_a
    global input_seq_b
    global res_text
    global align_opt
    global gap_opn_pnlty
    global extensn_pnlty

    seq_a = input_seq_a.get()
    seq_b = input_seq_b.get()
    gap_opn = float(gap_opn_pnlty.get())
    extensn = float(extensn_pnlty.get())

    clear_res_text(res_text)

    if valid_input(seq_a, seq_b, res_text):
        try:
            if align_opt == 'local':
                # print(gap_opn, 'gap_opn', extensn, 'extensn', 'local')
                seq_a_n, seq_n, seq_b_n, tot_re = get_pairwise2_localds_result(seq_a, seq_b, gap_opn, extensn)
                res_text.insert(CURRENT, tot_re)

            else:
                # print(gap_opn, 'gap_opn', extensn, 'extensn', 'global')
                st_seq_a = 1
                st_seq_b = 1
                seq_a_n, seq_n, seq_b_n, tot_re = get_pairwise2_globalds_result(seq_a, seq_b, gap_opn, extensn)
                for i in range((len(seq_a_n) // len_str) + 1):
                    tmp_seq_a = seq_a_n[i*len_str: (i + 1)*len_str]
                    tmp_seq_b = seq_b_n[i*len_str: (i + 1)*len_str]
                    seq_a_idx, st_seq_a = get_idx_seq(tmp_seq_a, st_seq_a)
                    seq_b_idx, st_seq_b = get_idx_seq(tmp_seq_b, st_seq_b)
                    res_text.insert(CURRENT, seq_a_idx + '\n')
                    res_text.insert(CURRENT, tmp_seq_a + '\n')
                    res_text.insert(CURRENT, seq_n[i*len_str: (i + 1)*len_str] + '\n')
                    res_text.insert(CURRENT, tmp_seq_b + '\n')
                    res_text.insert(CURRENT, seq_b_idx + '\n')
                    res_text.insert(CURRENT, '\n')
        except Exception as err:
            # messagebox.showerror('error', str(err))
            res_text.insert(CURRENT, 'error : \n' + str(err) + '\n')


if __name__ == '__main__':
    setupGUI()
    window.mainloop()
