# very dirty script for visualization of jets' consistuents 
# with their origin (mothers)
#
# run with:
# $ python fj2graphviz.py input_file_name
#
# script bases on:
#   print ancestors()
#   print("# JET nr ..")
#   print("==== EVENT .. ====")
# from reconstruction.C
#
# input format:
#============ EVENT 198 ==================
#
#		# JET nr 0
#
# 0.K_L0   <-K0_bar(-311)<-s(3)<-s(3)<-s(3)<-g(21)<-g(21)<-g(21)<-g(21)<-g(21)<-g(21)<-g(21)
#	 X <-(1270)<-(1245)<-(891)<-(575)<-(512)<-(183)<-(112)<-(80)<-(66)<-(53)<-(50)
#	 true first: g (21)
#	 true last: s (3)
#	 EXP (end): light
#
#------- true first=0  true last=3  exp=0
#------c: true first=0  true last=0  exp=0
#------b: true first=0  true last=0  exp=0
#
#
#              # JET nr 1
#    ... etc
#    ... etc
#


import os
import sys

if len(sys.argv) > 1: input_fname = sys.argv[1]
else: sys.exit('input filename not given!') 

with open(input_fname) as f_in, open('fj2graphviz/outputs/'+input_fname.replace('.txt', '.dot'), 'w') as f_out:
    head = 'digraph G {\n\nrankdir=LR;\nsize="20"\nminsep="1"\n\n'
    f_out.write(head)


    i_event, i_jet = -1,-1
    buffer = []    
    for line in f_in:
        if 'EVENT' in line: 
            i_event = int(line.split('EVENT')[1].split('=')[0])
            nodes_lst = []
            links_lst = []
            buffer = []
        if 'JET' in line:
            i_jet = int(line.split('nr')[1])
            jetname = 'jet{}_{}'.format(i_event, i_jet)
            jetlabel = jetname
            f_out.write('node [shape = box] {} [ label="{}" color=yellow style=filled]; \n'.format(jetname, jetlabel))
        if '<-' not in line: continue
        if not buffer: 
            buffer = line.split('<-')
            prevnode=''
            prevlabel=''
            continue
        for pid, pnr in zip(buffer, line.split('<-')):
            nodename = 'part{}_{}'.format(i_event, pnr.replace('(','').replace(')','').strip()).replace('\n','')
            nodelabel = pid.replace('\n','')
            if nodename.endswith('X'): nodename = nodename + '_' + str(i_jet) + '_' + nodelabel.split('.')[0].strip()


            if prevlabel != nodelabel or prevnode in nodes_lst:
                if nodename not in nodes_lst: 
                    color = ('color=lightgrey style=filled' if pid == buffer[-1] else '')
                    style = ('fixedsize=true height=0.1 width=0.1'  if (prevlabel == nodelabel) else '')
                    f_out.write('node [shape = circle] {} [ label="{}" {} {} ]; \n'.format(nodename, nodelabel, color, style))
                    nodes_lst.append(nodename)
                if prevnode and [nodename, prevnode] not in links_lst: 
                    style = ('[minlen=0]'  if (prevlabel == nodelabel) else '')
                    f_out.write('{}->{} {} \n'.format(nodename, prevnode, style))
                    links_lst.append([nodename,prevnode])
                prevnode = nodename
                prevlabel = nodelabel


            if pid == buffer[0]:
                f_out.write('{}->{} \n'.format(nodename, jetname))

        buffer = []
    f_out.write('}')


os.system('dot -Tpdf fj2graphviz/outputs/{} -ofj2graphviz/outputs/{}'.format( input_fname.replace('.txt', '.dot'), input_fname.replace('.txt', '.pdf') ))
os.system('dot -Tpng fj2graphviz/outputs/{} -ofj2graphviz/outputs/{}'.format( input_fname.replace('.txt', '.dot'), input_fname.replace('.txt', '.png') ))
