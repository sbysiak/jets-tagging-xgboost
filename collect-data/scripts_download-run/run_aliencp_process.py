#!/usr/bin/env python2.7
from __future__ import division
import os
import sys
import subprocess
from time import time
from aliencp import aliencp
from process_data import process_data_kinematics, process_data_aod
from math import ceil


#main_source_dir = '/alice/sim/2017/LHC17f8g/5/255539'
def check_disc_space(limit=25):
    """checks space left on disc and exits if is smaller than limit"""
    a = subprocess.check_output(['df', '.','-h']).decode('utf-8')
    disc_left = int(a[a.rfind('G ')-3:a.rfind('G ')])
    if  disc_left < limit: sys.exit('program aborted: too little memory ({})'.format(disc_left))



def run_aliencp_process(main_source_dir, input_type, limit_low=0, limit_high=99999, part_size=50):
    print '\n\nstarting run_aliencp_process ...'
    print 'parameters: {}'.format([main_source_dir, input_type, limit_low, limit_high, part_size])
    # adds '/' to the end if needed
    main_source_dir = os.path.join(main_source_dir, '')

    res = subprocess.check_output(['alien_ls', main_source_dir]).decode('ascii')
    res = res.split('\n')

    # define files names according to input_type
    if input_type == 'kinematics':
        fnames = {'download':'root_archive.zip',
                  'input':['Kinematics.root','galice.root'],
                  'output':['jetsTree.root', 'histos.root']}
    elif input_type == 'aod':
        f = main_source_dir
        suffix = f[f.find('LHC'):f.rfind(os.sep)].replace(os.sep,'_')
        if limit_high != 99999:
            suffix += '_part'+str(int(limit_low/part_size+1))
        idx = f.find('LHC')
    	lhc_dataset = f[idx : f.find(os.sep, idx)]

        fnames = {'download':'aod_archive.zip',
                  'input':['AliAOD.root', 'pyxsec_hists.root'], # must be in list
                  #'input':['AliAOD.root'],
                  'output':[suffix, ], # suffix
                  'dataset':lhc_dataset
                  }

    fnames['ignore_rm'] = fnames['output']
    save_inputs = False
    if save_inputs:
        for fn in fnames['input']: fnames['ignore_rm'].append(fn)


    if input_type == 'kinematics':
        for subdir in res:
            if main_source_dir.endswith(subdir): continue
            if  check_disc_space(): sys.exit('program aborted: too little memory ({})'.format(disc_left))

            os.system('date')
            start = time()

            acp_arg = os.path.join(main_source_dir, subdir)
            print('\nEXEC: rap: aliencp({}, {})'.format(acp_arg, fnames))
            aliencp( acp_arg, fnames )

            pd_arg = acp_arg.replace('/alice', 'DATA')
            print('\nEXEC: rap: process_data_kinematics({}, {})'.format(pd_arg, fnames))
            process_data_kinematics( pd_arg, fnames )

            print('\t--- exec. time of aliencp + process_data: {:.1f} min ---'.format((time()-start)/60.))

    elif input_type == 'aod':
        print 'entering AOD part...'
        if len(res[limit_low:limit_high]) > 2*part_size+12:  # 12 = usual no. additional, useless files
            # if there are too many files to process (in reasonable time), then split it in parts
            n_parts = int(ceil(len(res[limit_low:limit_high])/part_size))
            for i_part in range(n_parts):
                print 'about to run \"run_aliencp_process\" again...'
                run_aliencp_process(main_source_dir=main_source_dir,
                                    input_type=input_type,
                                    limit_low=i_part*part_size+limit_low,
                                    limit_high=(i_part+1)*part_size+limit_low,
                                    part_size=part_size)
        else:
            for subdir in res:
                try:
                    chunk_num = int(subdir.split(os.sep)[-1])
                except ValueError:
                    continue
                if chunk_num < limit_low or chunk_num >= limit_high: continue
                #print 'RAP: subdir=',subdir, ' main_source_dir=',main_source_dir
                if main_source_dir.endswith(subdir): continue
                acp_arg = os.path.join(main_source_dir, subdir)
                if 'AOD' in acp_arg: continue
                os.system('date')
                start = time()
                print('\nEXEC: rap: aliencp({}, {})'.format(acp_arg, fnames))
                aliencp( acp_arg, fnames )
                print('\t--- exec. time of aliencp: {:.1f} min ---'.format((time()-start)/60.))

            os.system('date')
            start = time()

            pd_arg = main_source_dir.replace('/alice', 'DATA')
            print('\nEXEC: rap: process_data({}, {})'.format(pd_arg, fnames))
            process_data_aod( pd_arg, fnames, (limit_low,limit_high) )
            print('\t--- exec. time of process_data: {:.1f} min ---'.format((time()-start)/60.))






if __name__ == '__main__':
    # run_aliencp_process('/alice/sim/2017/LHC17f8g/5/255539/', 'aod', limit_low=300, limit_high=350, part_size=50)
    # run_aliencp_process('/alice/sim/2017/LHC17f8g/15/255539/', 'aod', limit_low=150)
    # run_aliencp_process('/alice/sim/2017/LHC17f8g/7/255539/', 'aod')
    # run_aliencp_process('/alice/sim/2017/LHC17f8g/8/255539/', 'aod', limit_low=300, limit_high=350)
    run_aliencp_process('/alice/sim/2017/LHC17f8g/11/255539/', 'aod', limit_low=200, limit_high=500)
    # run_aliencp_process('/alice/sim/2017/LHC17f8g/11/255539/', 'kinematics')
    run_aliencp_process('/alice/sim/2017/LHC17f8g/11/255540/', 'aod')
    run_aliencp_process('/alice/sim/2017/LHC17f8g/11/255541/', 'aod')
    run_aliencp_process('/alice/sim/2017/LHC17f8g/11/255542/', 'aod')
    run_aliencp_process('/alice/sim/2017/LHC17f8g/11/255543/', 'aod')

    run_aliencp_process('/alice/sim/2017/LHC17f8g/11/255577/', 'aod')
    run_aliencp_process('/alice/sim/2017/LHC17f8g/11/255582/', 'aod')
    run_aliencp_process('/alice/sim/2017/LHC17f8g/11/255583/', 'aod')
