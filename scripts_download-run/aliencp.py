#!/usr/bin/env python3.6

import os
import sys
import subprocess




def aliencp(main_source_dir='/alice/sim/2017/LHC17f8g/1/255618/001', 
                main_target_dir='DATA', 
                max_generation=999, generation=0, counter=0):

    if generation > max_generation: return counter

    cmd_ls = 'alien_ls'
    cmd_cp = 'alien_cp'

    #print('\t*'*generation+'calling: $'+cmd_ls+' '+main_source_dir+' cnt=',counter)
    # call alien_ls
    res = subprocess.check_output([cmd_ls, main_source_dir]).decode('ascii')

    if 'root_archive.zip' in res:  # res = result from alien_ls

        target = main_target_dir + main_source_dir.replace('/alice', '')
        if not os.path.isdir(target): 
            subprocess.call(['mkdir', '--parents', target])

        # temporary !!! only for stopped jobs. Noramally indented below
        source_target = 'alien:///'+main_source_dir+'/*.zip '+target
        if (os.path.isfile(target+'/root_archive.zip') 
        and os.path.isfile(target+'/aod_archive.zip')
        and os.path.isfile(target+'/QA_archive.zip')):
            print('\tacp: all zip files already in dir')
            pass
        elif (os.path.isfile(target+'/jetsTree.root') 
        and os.path.isfile(target+'/histos.root')):
            print('\tacp: jetsTree and histos already in dir')
            pass
        else:
            print('EXEC: acp: $ '+cmd_cp+' '+source_target)
            # call alien_cp !!!
            os.system(cmd_cp+' '+source_target)
            counter += 1
            #subprocess.call([cmd_cp, source_target])
    else: 
        #print(res)
        for subdir in res.split('\n'):
            #print('\t*'*generation + 'subdir = '+subdir)
            if 'no such file or directory' in subdir:
                #print('\t*'*generation+'... continued (no such file or dir')
                continue
            if main_source_dir.endswith(subdir):
                #print('\t*'*generation+'... continued (last file)')
                continue
            if '.' in subdir:
                #print('\t*'*generation+'... continued (just file)')
                if subdir.endswith('.zip'): 
                    print("\n\tacp: inner:::alien_cp alien:///"+main_source_dir+'/'+subdir+"*.zip PATH'")
                    #counter += 1
                continue
            counter = aliencp(main_source_dir+'/'+subdir, generation=generation+1, counter=counter)

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

