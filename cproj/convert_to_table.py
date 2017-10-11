import numpy as np
from oct2py import octave as oc
# The form we will have to send to this function:
# a = x_reform,y_reform
# then input = [a1,a2,a3,a4,a5,a6]
Methods=['MC','MC-S','QMC-S','MLMC','MLMC-A',
    'FFT','FGL','COS',
    'FD','FD-NU','FD-AD',
    'RBF','RBF-FD','RBF-PUM','RBF-LSML','RBF-AD','RBF-MLT']
def makeTable(input):
    temp_t = [None]*6
    temp_err = [None]*6
    j = 0
    for i in input:
        temp_t[j] = oc.transpose(i[0])
        temp_err[j] = oc.transpose(i[1])
        j = j + 1
    if not oc.pwd() == "/proj/bench":
        oc.chdir("bench/")
    table = oc.feval("setupTable",temp_t,temp_err, nout = 1)
    latex = oc.feval("latexTable",table,nout = 1)
    return latex
