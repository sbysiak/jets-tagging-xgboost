import os
import glob


def run_hadd(main_dir='DATA/sim/2017/LHC17f8g/1/255618/001/', filename_suffix='', out_dir='DATA/'):
    histos_lst = []
    jetsTree_lst = []
    for path, subdirs, files in os.walk(main_dir):
        for name in files:
            #if name == 'Kinematics.root':
            if name == 'histos.root':
                #print('\n\n\nRUNNING: aliroot -l -q \'run_recon.C({},\"{}\")\' '.format( n_events, path ))
                histos_lst.append(os.path.join(path, name))
            if name == 'jetsTree.root':
                #print('\n\n\nRUNNING: aliroot -l -q \'run_recon.C({},\"{}\")\' '.format( n_events, path ))
                jetsTree_lst.append(os.path.join(path, name))

    if not filename_suffix:
        suff_lst = main_dir.split('/')
        if suff_lst[-1] == '': suff_lst = suff_lst[:-1]
        filename_suffix = '_'.join(suff_lst[suff_lst.index('2017')+1:])


    #cmd_hadd = 'hadd DATA/histos_hadd_1.root'
    cmd_hadd = 'hadd '+os.path.join(out_dir, 'histos_hadd')+'_'+filename_suffix+'.root'
    for file in histos_lst:
        cmd_hadd += (' '+file)
    print('\n\n\nRUNNING: '+ cmd_hadd)
    os.system(cmd_hadd)

    #cmd_hadd = 'hadd DATA/jetsTree_hadd_1.root'
    cmd_hadd = 'hadd '+os.path.join(out_dir, 'jetsTree_hadd')+'_'+filename_suffix+'.root'
    for file in jetsTree_lst:
        cmd_hadd += (' '+file)
    print('\n\n\nRUNNING: '+ cmd_hadd)
    os.system(cmd_hadd)


if __name__ == '__main__':


    #run_hadd(main_dir='DATA/sim/2017/LHC17f8g/1/255539')
    #run_hadd(main_dir='DATA/sim/2017/LHC17f8g/15/255539')
    #run_hadd(main_dir='DATA/sim/2017/LHC17f8g/20/255539')

    #run_hadd(main_dir='DATA/sim/2017/LHC17f8g/1/255618')
    #run_hadd(main_dir='DATA/sim/2017/LHC17f8g/15/255618')
    #run_hadd(main_dir='DATA/sim/2017/LHC17f8g/20/255618')

    #run_hadd(main_dir='DATA/sim/2017/LHC17f8g/1')
    #run_hadd(main_dir='DATA/sim/2017/LHC17f8g/15')
    #run_hadd(main_dir='DATA/sim/2017/LHC17f8g/20')

    #run_hadd(main_dir='DATA/sim/2017/LHC17f8g')
    run_hadd(main_dir='DATA/sim/2017/LHC17f8g/5/255540')
