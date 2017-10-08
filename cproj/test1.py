from oct2py import octave as oc
import os






def test1():
    print "Testing OS.dir and feval"

    x,y = oc.feval("prob1aI",nout=2)
    print "1st", x
    print "2st", y

def test2():
    print "Testing OS.dir and eval"
    
    x,y = oc.eval("prob1aI()",nout=2)
    print "1st", x
    print "2st", y



if __name__ == "__main__":
    oc.chdir("bench/")
    test1()
    test2()
