indexes_ferm = {}
indexes_us = {}
indexes_wmmf = {}

jain = {}
thr = {}
log = {}

num_exp = 10

for density in range(5,11):
	infile = open('res/'+str(density)+'/comparison.txt')
	stats = []
	for j in range(num_exp):
		vals = []
		infile.readline() 
		for i in range(5):
			infile.readline() 
			line = infile.readline() 
			line.strip()
			line = line.split()
			enb = (eval(line[3][1:-1]),eval(line[4][:-1]),eval(line[5][:-1]))

			line = infile.readline() 
			line.strip()
			line = line.split()
			ue = (eval(line[3][1:-1]),eval(line[4][:-1]),eval(line[5][:-1]))
			vals.append([enb,ue])
		stats.append(vals)


	thruput_ferm = []
	thruput_us = []
	thruput_wmmf = []

	log_ferm = []
	log_us = []
	log_wmmf = []

	jain_ferm = []
	jain_us = []
	jain_wmmf = []

	for run in range(num_exp):
		gran = 1

		mode = 0
		thruput_ferm.append(stats[run][mode][gran][0])
		log_ferm.append(stats[run][mode][gran][1])
		jain_ferm.append(stats[run][mode][gran][2])

		mode = 2
		thruput_us.append(stats[run][mode][gran][0])
		log_us.append(stats[run][mode][gran][1])
		jain_us.append(stats[run][mode][gran][2])

		mode = 4
		thruput_wmmf.append(stats[run][mode][gran][0])
		log_wmmf.append(stats[run][mode][gran][1])
		jain_wmmf.append(stats[run][mode][gran][2])

	z = 0
	for n in range(len(log_ferm)):
		if (log_ferm[n]==-1):
			log_ferm[n] = 0
			log_us[n] = 0
			log_wmmf[n] = 0
			z-=1
		if (log_us[n]==-1):
			log_ferm[n] = 0
			log_us[n] = 0
			log_wmmf[n] = 0
			z-=1
		if (log_wmmf[n]==-1):
			log_ferm[n] = 0
			log_us[n] = 0
			log_wmmf[n] = 0
			z-=1


	if -1 in log_ferm:
		l = -1
	else:
		l=sum(log_ferm)/(num_exp-z)
	indexes_ferm[density] = [sum(thruput_ferm)/num_exp,l,sum(jain_ferm)/num_exp]

	if -1 in log_us:
		l = -1
	else:
		l=sum(log_us)/(num_exp-z)
	indexes_us[density] = [sum(thruput_us)/num_exp,l,sum(jain_us)/num_exp]

	if -1 in log_wmmf:
		l = -1
	else:
		l=sum(log_wmmf)/(num_exp-z)
	indexes_wmmf[density] = [sum(thruput_wmmf)/num_exp,l,sum(jain_wmmf)/num_exp]

	jain[density] = [indexes_ferm[density][2],indexes_us[density][2],indexes_wmmf[density][2]]
	thr[density] = [indexes_ferm[density][0],indexes_us[density][0],indexes_wmmf[density][0]]
	log[density] = (indexes_ferm[density][1],indexes_us[density][1],indexes_wmmf[density][1])



print('JAIN')
for i in jain:
	print (i,'\t',jain[i][0],'\t',jain[i][1])

print('THRUPUT')
for i in thr:
	print (i,'\t',thr[i][0],'\t',thr[i][1])

print('log')
for i in thr:
	print (i,log[i][0],log[i][1])

