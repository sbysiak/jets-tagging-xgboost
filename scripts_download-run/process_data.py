import os
import glob


# run unzip, fastjet and rm unnecessary files (e.g. all but FJ results)


def run_unzip(path):
	#unzip DATA/sim/2017/LHC17f8g/1/255618/010/root_archive.zip -d DATA/sim/2017/LHC17f8g/1/255618/010/
	# -n = dont overwrite, -d = extract to external dir
	cmd_unzip = 'unzip -n '+os.path.join(path, 'root_archive.zip') + ' -d '+path
	print('\nEXEC: pd: '+cmd_unzip)
	os.system(cmd_unzip)


def run_recon(path, n_events):
    cmd_recon = 'aliroot -l -q \'run_recon.C({},\"{}\")\' '.format( n_events, path)
    print('\nEXEC: pd: '+cmd_recon)
    os.system(cmd_recon)

def run_rm(path):
    # REMOVE
    #rm $(ls -I "Kinematics.root" -I "galice.root" -I "jetsTree.root" -I "histos.root" -I "*.zip") # -I file_to_ignore

    #save_lst = ['jetsTree.root', 'histos.root', '*.zip', '*.C']
    #save_lst = ['jetsTree.root', 'histos.root', '*.C', 'Kinematics.root', 'galice.root']
    save_lst = ['jetsTree.root', 'histos.root', '*.C']
    # create empty file, cause without it LS doesnt ignore files after -I
    cmd_rm = '(cd '+path+'; touch file-torm; rm $(ls '
    for file in save_lst: cmd_rm += ' -I "{}"'.format(file)
    cmd_rm += ') )'
    print('\nEXEC: pd: '+cmd_rm)
    os.system(cmd_rm)


def process_data(main_dir='DATA/sim/2017/LHC17f8g/1/255618/001', n_events=200):
    for path, subdirs, files in os.walk(main_dir):
        if 'jetsTree.root' in files and 'histos.root' in files: 
            print('\tpd: jetsTree and histos already in dir')
            run_rm(path)
            continue

        if 'root_archive.zip' in files:
            run_unzip(path)
            run_recon(path, n_events)
            run_rm(path)

        elif 'aod_archive.root' in files or 'QA_archive.root' in files: 
            print('\tpd: only aod_archive and QA_archive -- removing') 
            run_rm(path)

        elif 'Kinematics.root' in files and 'galice.root' in files:
            run_recon(path, n_events)
            run_rm(path)
                

if __name__ == '__main__':


    process_data('DATA/sim/2017/LHC17f8g/15/')

            
