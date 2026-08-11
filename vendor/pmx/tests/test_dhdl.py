#!/usr/bin/env python
import pytest
from pmx.scripts.analyze_dhdl import entry_point
import os
import numpy as np
from numpy.random import normal,seed
from numpy.testing import assert_almost_equal
import sys


def test_dhdl(tmpdir):
    #make fake dHdl files in a temp directory
    trjlen=51
    dif=3.0
    W_scale=100.0
    dG=100.0
    noise_scale=10.0
    nruns=80
    l=np.linspace(0., 1., num=trjlen, endpoint=True)
    signal=(1/(1+np.exp((-l+0.5)*10))-0.5)*W_scale + dG
    
    np.random.seed(123456) #force reproducible random numbers
    
    fns_f=[]
    fns_b=[]
    for i in range(nruns):
        f=signal + 0.5*dif + np.random.normal(loc=0.0, scale=noise_scale, size=signal.shape)
        b=signal - 0.5*dif + np.random.normal(loc=0.0, scale=noise_scale, size=signal.shape)
        
        dataf=np.transpose(np.vstack((l,f)))
        datab=np.transpose(np.vstack((l,np.flip(b)))) # don't flip l. it's used as time here, not lambda value
        fh_f=tmpdir.join("test_temp_f%d.xvg"%i)
        fh_b=tmpdir.join("test_temp_b%d.xvg"%i)
        np.savetxt(fh_f, dataf, fmt='%16.8f')
        np.savetxt(fh_b, datab, fmt='%16.8f')
        
        fns_f.append(str(fh_f))
        fns_b.append(str(fh_b))
        
    #run analysis on our fake files
    orig_argv = sys.argv
    orig_dir = os.getcwd()
    os.chdir(tmpdir)
    sys.argv = [['analyze_dhdl.py'],\
                ['-fA'], [x for x in fns_f],\
                ['-fB'], [x for x in fns_b]]
    sys.argv = [item for sublist in sys.argv for item in sublist]
    entry_point()
    sys.argv = orig_argv #reset argv
    os.chdir(orig_dir)   #reset cwd

    #read results
    bar=-1
    barerr=-1
    print("#"*40)
    with open(str(tmpdir)+"/results.txt","rb") as fp:
        lines = fp.readlines()
        for line in lines:
            line=line.decode()
            print(line)
            if("BAR: dG =" in line):
                s=line.split()
                print(s)
                bar=float(s[-2])
            elif("BAR: Std Err (analytical) =" in line):
                s=line.split()
                print(s)
                barerr=float(s[-2])
    
    assert_almost_equal(bar, dG*0.98, decimal=2)
    assert_almost_equal(barerr, 0.17, decimal=2)
    