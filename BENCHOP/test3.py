from oct2py import octave as oc
def test1():
    print "Testing OS.dir and feval"
    x,y = oc.feval("prob1I",nout=2)
    print "1st", x
    print "2st", y

def test2():
    print "Testing OS.dir and eval"
    x,y = oc.eval("prob1I()",nout=2)
    print "1st", x
    print "2st", y



if __name__ == "__main__":
    test1()
    test2()
