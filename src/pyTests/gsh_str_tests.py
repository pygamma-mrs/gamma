from __future__ import division
from __future__ import print_function

import pygamma as pg 

def printsep(val, f):
    sep = '-----------------------------------------------------------'
    print(val, file=f)
    print(sep+'\n', file=f)    

infile = 'gsh_test.sys' 
outfile = "gsh_str_tests.txt" 
h1 = "Testing object print functionality"
h2 = "using input sys file: " + infile 


with open(outfile, 'w') as f:
    print(h1+'\n', file=f)
    printsep(h2+'\n', f)

    # Basics.Isotope
    iso = pg.Isotope()
    print(iso, file=f)
    iso = pg.Isotope('31P')
    printsep(iso, f)

    # Basics.IsotopeData
    iso = pg.IsotopeData()
    print(iso, file=f)
    iso = pg.IsotopeData('1H')
    printsep(iso, f)

    sys = pg.spin_system()
    print(sys, file=f)          # empty
    sys.read(infile)
    printsep(sys, f)  

    # HSLib.GenOp
    # nb. empty gen_op() throws exception
    gop = pg.Hcs(sys)    
    printsep(gop, f)
    gop_list = gop.toList()
    printsep(gop_list, f)
    gop_np = gop.toNParray()
    printsep(gop_np, f)

    # HSLib.SpinOp
    ax = pg.spin_sys(2)
    sop = pg.spin_op()    
    printsep(sop, f)
    sop = pg.Fx(ax)
    printsep(sop, f)
    sop = pg.Fx(ax) + pg.Fy(ax)
    printsep(sop, f)
    mx = pg.matrix()
    printsep(mx, f)
    mx = sop
    printsep(mx, f)

    # Matrix.col_vector
    col = pg.col_vector()
    printsep(col, f)
    col = pg.col_vector(4)
    printsep(col, f)
    col[0] = 5
    printsep(col, f)
    col[1] = 9.2
    col[3] = 11.6
    printsep(col, f)
    printsep(col+col, f)
    printsep(col*3, f)
     
    # Matrix.row_vector
    row = pg.row_vector()
    printsep(row, f)
    row = pg.row_vector(4)
    printsep(row, f)
    row[0] = 5
    printsep(row, f)
    row[1] = 9.2
    row[3] = 11.6
    printsep(row, f)
    printsep(row+row, f)
    printsep(row*3, f)
     
    # Matrix.complex 
    val = pg.complex(0)         # nb. empty complex() has randon numbers in it
    printsep(val, f)
    val = pg.complex(4, 3.2)
    printsep(val, f)
    val2 = pg.complex(val)
    printsep(val2, f)
    val3 = pg.complex(val+1)
    printsep(val3, f)
    

    
    
    

