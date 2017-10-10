import numpy as np
from oct2py import octave as oc

Methods=['MC','MC-S','QMC-S','MLMC','MLMC-A',
    'FFT','FGL','COS',
    'FD','FD-NU','FD-AD',
    'RBF','RBF-FD','RBF-PUM','RBF-LSML','RBF-AD','RBF-MLT']
def makeTable(input):
    temp_t = [None]*6
    temp_err = [None]*6
    for i in input:
        temp_t = np.array(i[0])
        temp_err = np.array(i[1])
    if not oc.pwd() == "/proj/bench":
        oc.chdir("bench/")
    table = oc.feval("setupTable",temp_t,temp_err, nout = 1)
    latex = oc.feval("latexTable",table,nout = 1)
    return latex
