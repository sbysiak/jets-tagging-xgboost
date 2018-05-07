#!/usr/bin/env python3.6

import os
import sys
import subprocess




def aliencp(main_source_dir='/alice/sim/2017/LHC17f8g/1/255618/001',
                fnames={},
                main_target_dir='DATA',
                max_depth=999, depth=0,
                counter=0):

    # prevent infinite looping
    if depth > max_depth: return counter

    if not all(k in fnames.keys() for k in ['download', 'input', 'output']):
        sys.exit('program aborted: incorrect fnames argument ({}), should contain \'download\', \'input\' and \'output\' keys'.format(fnames))

    fname_download = fnames['download']
    fname_input = fnames['input']
    fname_output = fnames['output']

    # fname_download = 'root_archive.zip'
    # fname_input = ['Kinematics.root', 'galice.root']
    # fname_output = ['jetsTree.root', 'histos.root']
    #
    #


    cmd_ls = 'alien_ls'
    cmd_cp = 'alien_cp'

    # call alien_ls
    res = subprocess.check_output([cmd_ls, main_source_dir]).decode('ascii')

    if fname_download in res:  # res = result from alien_ls

        target = main_target_dir + main_source_dir.replace('/alice', '')
        if not os.path.isdir(target):
            subprocess.call(['mkdir', '--parents', target])

        source_target = 'alien:///'+main_source_dir+'/'+fname_download+' '+target

        if (os.path.isfile(target+'/'+fname_download)):
            # print('ls '+target+' : ')
            # os.system('ls '+target)
            print('\tacp: file to download ({}) already in dir'.format(fname_download))
            pass
        elif all( [os.path.isfile(target+'/'+f) for f in fname_output] ) :
            print('\tacp: output files ({}) already in dir'.format(fname_output))
            pass
        elif all( [os.path.isfile(target+'/'+f) for f in fname_input] ) :
            print('\tacp: input files ({}) already in dir'.format(fname_input))
            pass
        else:
            print('EXEC: acp: $ '+cmd_cp+' '+source_target)
            # calling alien_cp !!!
            os.system(cmd_cp+' '+source_target)
            counter += 1
            #subprocess.call([cmd_cp, source_target])
    else:
        for subdir in res.split('\n'):
            #print 'ACP: subdir=',subdir, ' main_source_dir=',main_source_dir
            if 'no such file or directory' in subdir:
                continue
            if main_source_dir.endswith(subdir):
                continue
            if '/AOD' in subdir:
                continue
            if '.' in subdir:
                if subdir.endswith('.zip'):
                    print("\n\tacp: inner:::alien_cp alien:///"+main_source_dir+'/'+subdir+"*.zip PATH'")
                continue

            counter = aliencp(main_source_dir=main_source_dir+'/'+subdir,
                              fnames=fnames,
                              main_target_dir=main_target_dir,
                              max_depth=max_depth,
                              depth=depth+1,
                              counter=counter
                              )

    # logging
    return counter


if __name__ == '__main__':

    #alien_cp('/alice/sim/2017/LHC17f8g/1/255618/')
    counter = aliencp('/alice/sim/2017/LHC17f8g/1/255618/')
    print('there were {} *.zip files of each type downloaded'.format(counter))

#alien_file = ''
#destination = ''
#if len(sys.argv) > 1:
#    alien_file = sys.argv[1]
#if len(sys.argv) > 2:
#    destination = sys.argv[2]




#for dir in dirlst:
#    if os.path.exists(dir):
#        print(dir+' exists')
#    else:
#        print(dir+' doesnt exist')
#        print('mkdir --parents '+dir)
#        #os.system('mkdir --parents '+dir)


#print 'alien:',alien_file, 'dest: ', destination

#command = 'alien_cp alien:///alice/sim/2017/LHC17f8g/20/255618/001/*.zip ~/Pulpit/MASTER_THESIS/DATA'
#os.system('alien_cp alien:///alice/sim/2017/LHC17f8g/20/255618/001/*.zip ~/Pulpit/MASTER_THESIS/DATA')
