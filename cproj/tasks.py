from __future__ import absolute_import, unicode_literals
from .celery import app
from oct2py import Oct2Py
import numpy as np

@app.task
def solveproblem(pn):
    oc = Oct2Py()
    if not oc.pwd() == "/proj/bench":
        oc.chdir("bench/")
    rel_err,time = oc.feval(pn,nout=2)
    rel_err_reform = rel_err.tolist()
    time_reform = time.tolist()
    return rel_err_reform,time_reform


@app.task
def solveproblem_par(pn,par_dic):
    oc = Oct2Py()
    if not oc.pwd() == "/proj/bench":
        oc.chdir("bench/")
    rel_err,time = oc.feval(pn +"par",np.array(par_dic['S']),
                float(par_dic['K']),float(par_dic['T']),par_dic['r'],
                par_dic['sig'],np.array(par_dic['U']),nout=2)
    rel_err_reform = rel_err.tolist()
    time_reform = time.tolist()
    return rel_err_reform,time_reform
