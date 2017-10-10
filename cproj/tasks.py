from __future__ import absolute_import, unicode_literals
from .celery import app
from oct2py import Oct2Py

@app.task
def solveproblem(pn):
    oc = Oct2Py()
    if not oc.pwd() == "/proj/bench":
        oc.chdir("bench/")
    x,y = oc.feval(pn,nout=2)
    x_reform = [e for e in x[0]]
    y_reform = [e for e in y[0]]
    return x_reform,y_reform
