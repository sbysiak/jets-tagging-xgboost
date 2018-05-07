from __future__ import division
import os
import glob
import subprocess


# run unzip, fastjet and rm unnecessary files (e.g. all but FJ results)


def run_unzip(path, zip_name):
	#unzip DATA/sim/2017/LHC17f8g/1/255618/010/root_archive.zip -d DATA/sim/2017/LHC17f8g/1/255618/010/
	# -n = dont overwrite, -d = extract to external dir
	cmd_unzip = 'unzip -n '+os.path.join(path, zip_name) + ' -d '+path
	print('\nEXEC: pd: '+cmd_unzip)
	os.system(cmd_unzip)


def run_fastjet_kinematics(path, n_events):
    cmd_recon = 'aliroot -l -q \'run_recon.C({},\"{}\")\' '.format( n_events, path)
    print('\nEXEC: fj_kine: '+cmd_recon)
    os.system(cmd_recon)


def run_fastjet_aod(lhc_dataset, output_suffix):
	# making list of files to process + running train + renaming it
	# # # making list of files
	# scenario 1: files are downloaded between consecutive executions
	# scenario 2: there are a lot of files downloaded in advance
	#
	# alternative: processing one file per execution -- drawbacks - copying all files N_files times
	# # #
	dest = '../DATA/trains/train_{}'.format(output_suffix)
	cmd_recon = 'cd collect-data; ./run_JetExtractor.sh {}; rm {} -r; mv train {}; cd -'.format(lhc_dataset, dest, dest)
	print('\nEXEC: fj_aod: '+cmd_recon)
	os.system(cmd_recon)


def run_rm(path, files_to_ignore):
    # REMOVE
    #rm $(ls -I "Kinematics.root" -I "galice.root" -I "jetsTree.root" -I "histos.root" -I "*.zip") # -I file_to_ignore

    #save_lst = ['jetsTree.root', 'histos.root', '*.zip', '*.C']
    #save_lst = ['jetsTree.root', 'histos.root', '*.C', 'Kinematics.root', 'galice.root']
    #save_lst = ['jetsTree.root', 'histos.root', '*.C']

	# create empty file, cause without it 'ls' doesnt ignore files after -I
	# mind () around command executing it in sub-shell
	cmd_rm = '(cd '+path+'; touch file-torm; rm $(ls '
	for file in files_to_ignore: cmd_rm += ' -I "{}"'.format(file)
	cmd_rm += ' -I "*.C" -I "*.py"' # security, no need probably
	cmd_rm += ') )'
	print('\nEXEC: pd: '+cmd_rm)
	os.system(cmd_rm)


def process_data_kinematics(main_dir='DATA/sim/2017/LHC17f8g/1/255618/001', fnames={}, n_events=200):

    for path, subdirs, files in os.walk(main_dir):
        if all( [f in files for f in fnames['output']] ):
            print('\tpd: output files already in dir')
            run_rm(path, fnames['ignore_rm'])

        elif all( [f in files for f in fnames['input']] ):
            run_fastjet_kinematics(path, n_events)
            run_rm(path, fnames['ignore_rm'])

        elif fnames['download'] in files:
            run_unzip(path, fnames['download'])
            run_fastjet_kinematics(path, n_events)
            run_rm(path, fnames['ignore_rm'])


def process_data_aod(main_dir='DATA/sim/2017/LHC17f8g/1/255618/001', fnames={}, limits=(0,99999)):
	limit_low, limit_high = limits
	# unzip files -- rm all but AliAOD.root -- run FJ -- rm AliAOD if not in ignored
	for path, subdirs, files in os.walk(main_dir):
		if fnames['download'] in files:
			run_unzip(path, fnames['download'])
			run_rm(path, fnames['input'])

	# prepare file with list of AOD files
	find_cmd = 'find '+os.path.join(os.getcwd(), main_dir)+' -name "{}" | sort'.format(fnames['input'][0])
	all_aod = subprocess.check_output(find_cmd, shell=True).split('\n')
	print 'searching for AODs: \n' + find_cmd
	input_files = 'collect-data/'+fnames['dataset']+'_root_archive_AliAOD.txt'
	with open(input_files, 'w') as file:
	    print('in file opened scope')
	    print('all_aod = {}'.format(all_aod))
	    for aod_file in all_aod:
	        # e.g. '/home/...../DATA/sim/2017/LHC17f8g/1/255618/056/AliAOD.root'
        	if not aod_file: continue
	        chunk_num = int(aod_file.split(os.sep)[-2])

        	# print('chunk num={}, limit_low={}, limit_high={}'.format(chunk_num,limit_low,limit_high))
	        if chunk_num >= limit_low and chunk_num <= limit_high:
	            print 'writing: ', aod_file, '...'
	            file.write(aod_file+'\n')

	os.system('wc -l '+input_files)
	run_fastjet_aod(fnames['dataset'], fnames['output'][0])

	for path, subdirs, files in os.walk(main_dir):
		if path == main_dir: continue
		#print path
		chunk_num = int(path.split(os.sep)[-1])
		# print('removing: chunk num={}, limit_low={}, limit_high={}'.format(chunk_num,limit_low,limit_high))
		if chunk_num >= limit_low and chunk_num <= limit_high:
			run_rm(path, fnames['ignore_rm'])






if __name__ == '__main__':


    process_data('DATA/sim/2017/LHC17f8g/15/')
