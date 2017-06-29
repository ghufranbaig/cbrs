import math


#Log functions

def fairness_index(shares):
# Fairness index calculation
	vals = list(shares.values())
	num_vals = len(vals)
	if 0 in vals:
		log_vals = []
	else:
		log_vals = [math.log10(i) for i in vals]
	sq_vals = [i**2 for i in vals]
	jains_index = 0
	# jains_index = float((sum(vals)**2))/(num_vals*(sum(sq_vals)))
	# jains_index = eval("%.2f" % jains_index)
	sum_ = eval("%.2f" % sum(vals))
	if len(log_vals)>0 :
		l_sum = eval("%.2f" % sum(log_vals))
	else:
		l_sum = -1
	
	return (sum_,l_sum,jains_index)


def writeInfo(scheme,assign,share_eNB,share_UE,info):
		info.write(scheme + '\n')

		info.write('Assignment per ENB\n')
		info.write(str(assign))
		info.write('\n')
		info.write('Share per ENB: '+str(fairness_index(share_eNB))+'\n')
		info.write(str(share_eNB))
		info.write('\n')
		info.write('Share per UE: '+str(fairness_index(share_UE))+'\n')
		info.write(str(share_UE))
		info.write('\n')
		info.write('\n')

def writegeneralInfo(info,enb_coord,ue_list_to_print,UEs,load,i_map,i_map_dir,ue_interfering,edges,G,u_m):
	info.write('ENBs\n')
	info.write(str(enb_coord))
	info.write('\n')
	info.write('UE LIST\n')
	info.write(str(ue_list_to_print))
	info.write('\n')
	tmp={}
	for a in UEs:
		tmp[a]=[]
		for b in UEs[a]:
			tmp[a].append(u_m[b])
	info.write('UEs Association\n')
	info.write(str(tmp))
	info.write('\n')


	info.write('\n')

        info.write('LOAD\n')
        info.write(str(load))
        info.write('\n')
        info.write('I_MAP\n')
        info.write(str(i_map))
        info.write('\n')
        info.write('I_MAP_DIR\n')
        info.write(str(i_map_dir))
        info.write('\n')
        info.write('UE_INTERFERING\n')
        info.write(str(ue_interfering))
        info.write('\n')
        tmp = {}
        for i in edges:
                tmp[i] = []
                for j in edges[i]:
                        tmp[i].append(u_m[j])
	info.write(str(tmp))
	info.write('\n')
	info.write('GRAPHS\n')
	info.write(str(G))
	info.write('\n')
	info.write('\n')


# Plot functions
def plot_graph(outputDir,filename,i_map,enb_coord,u_m,UEs,l,w):
	#plot topology using Matlab/Octave
	pf = open(outputDir + filename + '.m','w')

	num = len(enb_coord)
	pf.write('clr = hsv('+str(num)+');\n')
	#pf.write("figure('Position',[-10,-10,"+str(l+10)+","+str(w+10)+"]);\n")
	pf.write("f = figure('visible','off');\n") 
	pf.write('hold on;\n')
	delta = 1

	for i in range(num):
		for n in i_map[i]:
			pf.write('plot(['+str(enb_coord[i][0])+' '+ str(enb_coord[n][0])+'],['+str(enb_coord[i][1])+' '+ str(enb_coord[n][1])+'],\'--\',\'Color\',[0.5 0.5 0.5],\'LineWidth\',1);\n')


	for i in range(num):
		pf.write('plot('+str(enb_coord[i][0])+','+str(enb_coord[i][1])+',\'o\',\'Color\',clr('+str(i+1)+',:),\'LineWidth\',3);\n')
		pf.write('text('+str(enb_coord[i][0]-15)+','+str(enb_coord[i][1])+','+'"'+str(i)+'");\n')
		for j in UEs[i]:
			pf.write('plot('+str(j[0])+','+str(j[1])+',\'*\',\'Color\',clr('+str(i+1)+',:),\'LineWidth\',1);\n')
			pf.write('text('+str(j[0]-10)+','+str(j[1])+','+'"'+str(u_m[j])+'","fontsize",5);\n')


	
	pf.write("axis('Position',[-10 "+str(l+10)+" -10 "+str(w+10)+"]);\n")
	pf.write('print '+outputDir+filename+'.jpg;\n')

	pf.write('hold off;\n') 
	pf.close()


def plot_ue_interference(outputDir,filename,edges,enb_coord,u_m,UE,l,w):
	pf = open(outputDir + filename + '.m','w')

	num = len(enb_coord)
	pf.write('clr = hsv('+str(num)+');\n')
	#pf.write("figure('Position',[-10,-10,"+str(l+10)+","+str(w+10)+"]);\n")
	delta = 1

	pf.write('figure;\n') 
	pf.write('hold on;\n') 
	for i in range(num):
		for ue in edges[i]:
			pf.write('plot(['+str(enb_coord[i][0])+' '+ str(ue[0])+'],['+str(enb_coord[i][1])+' '+ str(ue[1])+'],\'--\',\'Color\',clr('+str(i+1)+',:),\'LineWidth\',0.5);\n')


	for i in range(num):
		pf.write('plot('+str(enb_coord[i][0])+','+str(enb_coord[i][1])+',\'o\',\'Color\',clr('+str(i+1)+',:),\'LineWidth\',3);\n')
		pf.write('text('+str(enb_coord[i][0]-15)+','+str(enb_coord[i][1])+','+'"'+str(i)+'");\n')
		for j in UEs[i]:
			pf.write('plot('+str(j[0])+','+str(j[1])+',\'*\',\'Color\',clr('+str(i+1)+',:),\'LineWidth\',1);\n')
			pf.write('text('+str(j[0]-10)+','+str(j[1])+','+'"'+str(u_m[j])+'","fontsize",5);\n')

	pf.write("axis('Position',[-10 "+str(l+10)+" -10 "+str(w+10)+"]);\n")
	pf.write('print '+outputDir+filename+'.jpg;\n')
	pf.write('hold off;\n') 
	pf.write('figure;\n') 
	pf.close()


