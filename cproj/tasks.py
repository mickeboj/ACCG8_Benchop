from __future__ import absolute_import, unicode_literals
from .celery import app
from oct2py import octave as oc

@app.task
def solveproblem(pn):
    
