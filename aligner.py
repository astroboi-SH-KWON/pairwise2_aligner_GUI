from tkinter import *
from tkinter import ttk
from tkinter import messagebox
from tkinter import filedialog
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62
import pandas as pd
import os

############### start to set env ################
len_str = 50
window = None
input_seq_a = None
input_seq_b = None
input_file = None
image_label = None
load_image = None

align_radio = None
align_radio2 = None
align_opt = None
align_opt2 = None
gap_opn_pnlty = None
extensn_pnlty = None
gap_opn_pnlty2 = None
extensn_pnlty2 = None

res_text = None
log_text = None

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
    return align_arr[0], align_arr[1], align_arr[2], alignments_result
    # return ''.join([i for i in align_arr[0] if not i.isdigit()]), ''.join(
    #     [i for i in align_arr[1] if not i.isdigit()]), ''.join(
    #     [i for i in align_arr[2] if not i.isdigit()]), alignments_result


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


def set_align_opt2(val):
    global align_opt2
    global gap_opn_pnlty2
    global extensn_pnlty2

    gap_opn_pnlty2.delete(0, 'end')
    gap_opn_pnlty2.insert(0, 10.0)

    extensn_pnlty2.delete(0, 'end')
    if val == 'local':
        extensn_pnlty2.insert(0, 1.0)
    else:
        extensn_pnlty2.insert(0, 0.5)

    align_opt2 = val


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


def valid_input(seq_a, seq_b):
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


def reset2():
    global gap_opn_pnlty2
    global extensn_pnlty2
    global align_radio2
    global log_text
    global input_file

    gap_opn_pnlty2.delete(0, 'end')
    extensn_pnlty2.delete(0, 'end')
    input_file.delete(0, 'end')
    clear_res_text(log_text)

    gap_opn_pnlty2.insert(0, 10.0)
    extensn_pnlty2.insert(0, 0.5)
    align_radio2.set('global')


def upload_file():
    global input_file

    input_file.delete(0, 'end')
    input_file.insert(0, filedialog.askopenfilename())


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

    if valid_input(seq_a, seq_b):
        try:
            if align_opt == 'local':
                # print(gap_opn, 'gap_opn', extensn, 'extensn', 'local')
                seq_a_n, seq_n, seq_b_n, tot_re = get_pairwise2_localds_result(seq_a, seq_b, gap_opn, extensn)
                res_text.insert(CURRENT, seq_a_n + '\n')
                res_text.insert(CURRENT, seq_n + '\n')
                res_text.insert(CURRENT, seq_b_n + '\n')

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


def valid_file(df):
    len_seq_a = len(df[df.columns[0]])
    len_seq_b = len(df[df.columns[1]])

    if len_seq_a == len_seq_b:
        return True
    else:
        messagebox.showerror('error', 'check input file')
        return False


def get_align_file():
    global input_file
    global log_text
    global align_opt2
    global gap_opn_pnlty2
    global extensn_pnlty2

    gap_opn = float(gap_opn_pnlty2.get())
    extensn = float(extensn_pnlty2.get())

    input_f = input_file.get()
    df = pd.read_excel(input_f)
    clear_res_text(log_text)
    log_text.insert(CURRENT, 'loading input file :\n ' + input_f + '\n\n')

    if valid_file(df):
        result_dict = {'seq_A': [], 'alignment': [], 'seq_B': []}
        len_df_c = len(df[df.columns[0]])
        dv_num = 10
        len_unit = len_df_c // dv_num
        pcnt = 0
        try:
            for i in range(len_df_c):
                seq_a = df.loc[i][0].replace('\n', '').replace('\t', '')
                seq_b = df.loc[i][1].replace('\n', '').replace('\t', '')
                if align_opt2 == 'local':
                    seq_a_n, seq_n, seq_b_n, tot_re = get_pairwise2_localds_result(seq_a, seq_b, gap_opn, extensn)
                else:
                    seq_a_n, seq_n, seq_b_n, tot_re = get_pairwise2_globalds_result(seq_a, seq_b, gap_opn, extensn)

                result_dict['seq_A'].append(seq_a_n)
                result_dict['alignment'].append(seq_n)
                result_dict['seq_B'].append(seq_b_n)

                if len_df_c > dv_num and (i % len_unit == 0):
                    log_text.insert(CURRENT, '%3d %%\n' % pcnt)
                    pcnt += (100 // dv_num)

            result_df = pd.DataFrame.from_dict(result_dict)
            fn, f_ex = os.path.splitext(input_f)
            output_f = fn + '_align_result' + f_ex
            result_df.to_excel(output_f)
            log_text.insert(CURRENT, '\nDONE!\n\ncheck output file :\n ' + output_f + '\n')
        except Exception as err:
            # messagebox.showerror('error', str(err))
            log_text.insert(CURRENT, '\nerror : \n' + str(err) + '\n')


def setupGUI():
    global window
    global input_seq_a
    global input_seq_b
    global align_radio
    global align_radio2
    global gap_opn_pnlty
    global gap_opn_pnlty2
    global extensn_pnlty
    global extensn_pnlty2
    global input_file
    global image_label
    global load_image
    global res_text
    global log_text

    window = Tk()
    window.title('pairwise2 aligner')

    notebk = ttk.Notebook(window)
    notebk.pack()

    # first tab
    frame1 = Frame(window)
    notebk.add(frame1, text='1:1 aligner')

    # st align options
    align_opt_label1 = Label(frame1, relief=RAISED)
    align_opt_label1.grid(row=0, column=2, padx=5, pady=5)
    align_radio = StringVar()
    global_opt = ttk.Radiobutton(align_opt_label1, text='global', value='global', variable=align_radio, command=lambda: set_align_opt('global'))
    global_opt.grid(row=0, column=0)
    local_opt = ttk.Radiobutton(align_opt_label1, text=' local ', value='local', variable=align_radio, command=lambda: set_align_opt('local'))
    local_opt.grid(row=1, column=0)
    align_radio.set('global')

    align_opt_label2 = Label(frame1, relief=RAISED)
    align_opt_label2.grid(row=0, columnspan=2, padx=5, pady=5)
    gap_opn_label = Label(align_opt_label2, text='gap open penalty :', font='Courier 10 bold', relief=FLAT)
    gap_opn_label.grid(row=0, column=0, padx=5, pady=5)
    # input of gap open penalty
    gap_opn_pnlty = Entry(align_opt_label2, font='Terminal 10', width=5)
    gap_opn_pnlty.grid(row=0, column=1, padx=3, pady=3)
    gap_opn_pnlty.insert(0, 10.0)
    extensn_label = Label(align_opt_label2, text=' extension penalty :', font='Courier 10 bold', relief=FLAT)
    extensn_label.grid(row=0, column=2, padx=5, pady=5)
    # input of extension penalty
    extensn_pnlty = Entry(align_opt_label2, font='Terminal 10', width=5)
    extensn_pnlty.grid(row=0, column=3, padx=3, pady=3)
    extensn_pnlty.insert(0, 0.5)
    # en align options

    # st seq A label & its input
    seq_a_label = Label(frame1, text='seq A ', font='Courier 10 bold', relief=RAISED)
    seq_a_label.grid(row=1, column=0, padx=5, pady=5)
    input_seq_a = Entry(frame1, font='Terminal 10', width=50)
    input_seq_a.grid(row=1, column=1, padx=3, pady=3)
    # input_seq_a.bind('<Return>', aligner)
    # en seq A label & its input

    # st seq B label & its input
    seq_b_label = Label(frame1, text='seq B ', font='Courier 10 bold', relief=RAISED)
    seq_b_label.grid(row=2, column=0, padx=5, pady=5)
    input_seq_b = Entry(frame1, font='Terminal 10', width=50)
    input_seq_b.grid(row=2, column=1, padx=3, pady=3)
    input_seq_b.bind('<Return>', aligner)
    # en seq B label & its input

    # buttons
    aligner_btn = Button(frame1, text='align', font='Courier 20 bold', command=aligner, height=10)
    aligner_btn.grid(rowspan=1, column=2, padx=3, pady=3)

    reset_btn = Button(frame1, text='reset', font='Courier 10 bold', fg='red', command=reset)
    reset_btn.grid(row=1, column=2, padx=3, pady=3)

    # Text for result
    res_text = Text(frame1, font='Terminal 10', relief=RAISED, width=65, height=30)
    res_text.grid(row=3, columnspan=2, padx=5, pady=5)

    # second tab
    frame2 = Frame(window)
    notebk.add(frame2, text='n:n aligner')

    # st align options
    align_opt_label3 = Label(frame2, relief=RAISED)
    align_opt_label3.grid(row=0, column=2, padx=5, pady=5)
    align_radio2 = StringVar()
    global_opt2 = ttk.Radiobutton(align_opt_label3, text='global', value='global', variable=align_radio2,
                                 command=lambda: set_align_opt2('global'))
    global_opt2.grid(row=0, column=0)
    local_opt2 = ttk.Radiobutton(align_opt_label3, text=' local ', value='local', variable=align_radio2,
                                command=lambda: set_align_opt2('local'))
    local_opt2.grid(row=1, column=0)
    align_radio2.set('global')

    align_opt_label4 = Label(frame2, relief=RAISED)
    align_opt_label4.grid(row=0, columnspan=2, padx=5, pady=5)
    gap_opn_label = Label(align_opt_label4, text='gap open penalty :', font='Courier 10 bold', relief=FLAT)
    gap_opn_label.grid(row=0, column=0, padx=5, pady=5)
    # input of gap open penalty
    gap_opn_pnlty2 = Entry(align_opt_label4, font='Terminal 10', width=5)
    gap_opn_pnlty2.grid(row=0, column=1, padx=3, pady=3)
    gap_opn_pnlty2.insert(0, 10.0)
    extensn_label = Label(align_opt_label4, text=' extension penalty :', font='Courier 10 bold', relief=FLAT)
    extensn_label.grid(row=0, column=2, padx=5, pady=5)
    # input of extension penalty
    extensn_pnlty2 = Entry(align_opt_label4, font='Terminal 10', width=5)
    extensn_pnlty2.grid(row=0, column=3, padx=3, pady=3)
    extensn_pnlty2.insert(0, 0.5)
    # en align options

    # st file dialog
    input_file = Entry(frame2, font='Terminal 10', width=50)
    input_file.grid(row=1, column=0, padx=3, pady=3)
    file_button = Button(frame2, text='open',  font='Courier 10 bold', command=upload_file, width=10)
    file_button.grid(row=1, column=1, padx=3, pady=3)
    # en file dialog

    # reset2 btn
    reset_btn2 = Button(frame2, text='reset', font='Courier 10 bold', fg='red', command=reset2)
    reset_btn2.grid(row=1, column=2, padx=3, pady=3)

    # example image
    # os.chdir(os.getcwd())
    # load_image = PhotoImage(file=os.getcwd() + '/example.png')
    # load_image = PhotoImage(file='example.png')
    # image_label = Label(frame2, image=load_image, width=400, height=95)
    image_label = Label(frame2, width=400, height=95)
    image_label.grid(row=2, columnspan=2, padx=3, pady=3)
    example = [
            [' ', 'A', 'B']
            , [1, 'seq_A', 'seq_B']
            , [2, 'ATCGATGCTGATGCGCGACTAGTC', 'agtcgatcgtagcggtacgta']
            , [3, 'ATCGATGCTGATGCGCGACTAGTC', 'AGATCGTAGCTAGTCGATCTAG']
            , [4, 'atcgatgctgatgcgcgactagtc', 'agtcgatcgtagcggtacgta']
            , [5, 'atcgatgctgatgcgcgactagtc', 'AGTCGATCGTAGCGGTACGTA']
            , [6, '', '']
            ]
    for i in range(len(example)):
        for j in range(len(example[i])):
            bg = 'white'
            wid = 30
            font = 'Arial'
            if j == 0 or i == 0:
                bg = 'light gray'
                font = 'Arial bold'
            if j == 0:
                wid = 3

            e = Entry(image_label, width=wid, bg=bg, font=(font, 8))
            e.grid(row=i, column=j)
            e.insert(END, example[i][j])

    example_label = Label(frame2, text='sample\nexcel\ninput file', font='Courier 13 bold', fg='gray', relief=FLAT)
    example_label.grid(row=2, column=2, padx=1, pady=1)

    # Text for log
    log_text = Text(frame2, font='Terminal 10', relief=RAISED, width=65, height=23)
    log_text.grid(row=3, columnspan=2, padx=3, pady=3)

    aligner_btn2 = Button(frame2, text='align', font='Courier 20 bold', command=get_align_file, height=7)
    aligner_btn2.grid(row=3, column=2, padx=3, pady=3)


    # third tab
    frame3 = Frame(window)
    notebk.add(frame3, text='about')
    about = """
THIS APP WAS MADE WITH pairwise2 FROM biopython.
IT IS A FREE APP LIKE biopython.
YOU CAN DOWNLOAD AND EDIT THE ORIGINAL 
BY FOLLOWING THE LINK BELOW.
WE ARE NOT RESPONSIBLE FOR ANY DAMAGE AFTER DOWNLOAD AND USE.
_________________________________________________________
https://github.com/astroboi-SH-KWON/pairwise2_aligner_GUI
_________________________________________________________

이 app은 biopython 기반의 pairwise2를 사용하였습니다.
biopython처럼 이 app은 완전 무료 제품입니다.
github에 원래 소스를 받아서 수정하셔도 됩니다.
우리는 이 app에 의해 생기는 
사용자의 어떠한 손해도 책임지지 않습니다.
    """
    about_label = Label(frame3, font='Terminal 15', text=about, relief=RAISED, width=65, height=30)
    about_label.grid(row=1, columnspan=2, padx=5, pady=5)


if __name__ == '__main__':
    setupGUI()
    window.mainloop()
