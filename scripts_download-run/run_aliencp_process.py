#!/usr/bin/env python3.6

import os
import sys
import subprocess
from time import time
from aliencp import aliencp
from process_data import process_data


#main_source_dir = '/alice/sim/2017/LHC17f8g/5/255539'

def run_aliencp_process(main_source_dir):
    res = subprocess.check_output(['alien_ls', main_source_dir]).decode('ascii')
    for subdir in res.split('\n'):
        if main_source_dir.endswith(subdir): continue
    
        a = subprocess.check_output(['df', '.','-h']).decode('utf-8')
        disc_left = int(a[a.rfind('G ')-3:a.rfind('G ')])
        if  disc_left < 25: sys.exit('program aborted: too little memory ({})'.format(disc_left)) 
     
        os.system('date')
        start = time()
        acp_arg = os.path.join(main_source_dir, subdir)
        print('\nEXEC: rap: aliencp({})'.format(acp_arg))
        aliencp( acp_arg )
    
        pd_arg = acp_arg.replace('/alice', 'DATA')
        print('\nEXEC: rap: process_data({})'.format(pd_arg))
        process_data( pd_arg )    
    
        print('\t--- exec. time: {:.1f} min ---'.format((time()-start)/60.))


if __name__ == '__main__':
    for pt_bin in ['2','4']:
        run_aliencp_process('/alice/sim/2017/LHC17f8g/'+pt_bin+'/255539')
