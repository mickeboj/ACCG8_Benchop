from oct2py import octave as oc
import os
if __name__ =="__main__":
    os.chdir('bench/')
    [time,relerr] = oc.feval("prob1I")
