import random
import math
import os
import Queue as Q
from copy import deepcopy
from Fermi import Fermi
from KlogN import KlogN
from KlogN import KlogN_Backoff
from max_min import max_min_share
from max_min import WMMF
from cppScript import generate_topo_file
from cppScript import generate_topo_file_udp_flows
from cppScript import generate_topo_file_wifi_flows
from cppScript import generate_topo_file_wifiac_flows
from cppScript import generate_topo_file_udp_dynflows
from cppScript import generate_topo_file_wifi_dynflows
from Fermi import Fermi_RZ
from Fermi import Fermi_UE
from Fermi import Fermi_Pre_Alloc


outputDir = "res/"

SNR = 3
rad = 100
bw = 100
rbgsize = 4
spectrumBW = {6:1.4e6, 15:3.0e6, 25:5.0e6, 50:10.0e6, 75:15.0e6, 100:20.0e6}
SpectralEfficiencyForCqi = [
  0.0,
  0.15, 0.23, 0.38, 0.6, 0.88, 1.18,
  1.48, 1.91, 2.41,
  2.73, 3.32, 3.9, 4.52, 5.12, 5.55
]


# Taken from Rahman's model in ns-3
Noise = 1.65959e-19
# Noise = 3.98107e-21

# Radius of an eNodeB cell, so all users are at distance [rad_i, rad_o], in meters
#rad_o = 100
#rad_i = 30

# Thresholds for calculating interference (dist_thresh is currently not used)
dist_thresh = 180
cqi_thresh = 3



def dist(a,b):
	x = a[0] - b[0]
	y = a[1] - b[1]
	return math.hypot(x,y)

def angle(a,b):
	x = a[0] - b[0]
	y = a[1] - b[1]
	return math.atan2(y,x)

def dBmtoW(dbm):
	return (10**((dbm-30)/10.0))

def WtodBm(W):
	return (10*math.log10(W) + 30)

def RxPowertoDistance(P_r,P_t):
	n = 3.5
	L_0 = 25.686
	d_0 = 1
	return( d_0 * 10**( math.log10(P_t/P_r)/n - L_0/(n*10.0) ) )

def RxPwrtoTxPower (P_r,d):
	n = 3.5
	L_0 = 25.686
	d_0 = 1
	return ( P_r*(10**((L_0 + 10*n*math.log10(d/d_0))/10.0)) )

def getOverHearingRange(enb,ue):
	d = dist(enb,ue)
	TxP = RxPwrtoTxPower (dBmtoW(-90.0),d)
	m_range = RxPowertoDistance(dBmtoW(-100.0),TxP)
	return m_range

def rcvd_pows(bs_coords,ue):
        # Taken from Rahman's model in ns-3
	n = 3.5
	L_0 = 25.686
	d_0 = 1
        P_t = 1
	pows = {}
	dists = {}
	for i in range(len(bs_coords)):
		bs = bs_coords[i]
		d = dist(bs,ue)
		P_r = P_t/(10**((L_0 + 10*n*math.log10(d/d_0))/10))
		pows[i] = P_r
		dists[i] = d
	return (pows,dists)

def gen_eNbs_coord(n,l,w):
	coord = []
	for i in range (n):
		x=random.randint(0,l)
		y=random.randint(0,w)
		coord.append((x,y))
	return coord

def gen_ue_coord(enb,l,w,rad_i,rad_o):
	while(True):
		r = random.uniform(rad_i, rad_o)
		q = random.uniform(0, 2*math.pi)
		x = enb[0] + r*math.cos(q)
		y = enb[1] + r*math.sin(q)
		if ((x>0 and x<l) and (y>0 and y<w)):
			return(x,y)



def max_interferer(interferers,powers):
	max_i = interferers[0]
	for i in interferers:
		if (powers[i] > powers[max_i]):
			max_i = i
	return max_i



def GenerateGraphInfo(enb_coord,UEs):
	n = len(enb_coord)
	Distance = {}
	RxPower = {}
	CQIVals = {}

	for i in range(n):
		rcv_pow = []
		dists = []
		CQIs = []
		Actual_CQI = []

		class1 = [True for a in range(len(UEs[i]))]
		class1_S_E = [True for a in range(len(UEs[i]))]		


		ue_num = 0
		for u in UEs[i]:
			(powers,distances) = rcvd_pows(enb_coord,u)
			rcv_pow.append(powers)
			dists.append(distances)

			snr = {}
			cqis = {}
			for k in range(n):
				if(i==k):
					snr[k] = powers[i]/(Noise*spectrumBW[bw])
				else :
					snr[k] = powers[i]/(powers[k] + Noise*spectrumBW[bw])
				cqis[k] = getCQI(getSpectralEfficiency(snr[k]))
			CQIs.append(cqis)

		Distance[i] = dists
		RxPower[i] = rcv_pow
		CQIVals[i] = CQIs
	
	return (Distance,RxPower,CQIVals)

def GenerateGraph_CQI(Distance,RxPower,CQIVals,enb_coord,UEs,u_m):
	n = len(enb_coord)
	G = {} # intereference graph based on CQI difference
	G_Ue = {} # interference graph based on interference between UEs	
	Class1Clients = {}
	ReuseZoneSize = {}


	ues_intefering = {}
	intf_edges = {}

	for enb in UEs:
		#print (UEs[enb])
		for user in UEs[enb]:
			G_Ue[u_m[user]] = set()
			for user2 in UEs[enb]:
				if (u_m[user] == u_m[user2]):
					continue
				G_Ue[u_m[user]].add(u_m[user2])
	
	for i in range(n):
		G[i] = set()
		intf_edges[i]=set()
		ues_intefering[i] = len(UEs[i])
		Class1Clients[i] = []


	for i in range(n):
		rcv_pow = RxPower[i]
		dists = Distance[i]
		CQIs = CQIVals[i]

		class1 = [True for a in range(len(UEs[i]))]

		#######################################
		# interference map generation based on reported CQI values
		# iterates over all eNBs
		for j in range(n):
			if (i==j):
				continue
			#iterates over all the UEs for eNB j
			for z in range(len(rcv_pow)):
				u_pow = rcv_pow[z]
				u_dist = dists[z]
				u_cqi = CQIs[z]
				
				if (u_cqi[i] - u_cqi[j] > cqi_thresh):
				#if (u_pow[j]>thresh):
				#if (u_dist[j]<dist_thresh):

					class1[z] = False # if it interfering with some othe eNB, its not class 1
					
					G[i].add(j)
					G[j].add(i)

					ues_intefering[j] += 1
					intf_edges[j].add(UEs[i][z])

					for user in UEs[j]:
						G_Ue[u_m[UEs[i][z]]].add(u_m[user])
						G_Ue[u_m[user]].add(u_m[UEs[i][z]])
				
		#######################################

		# class 1 clients for eNB i
		numClass1 = 0
		for c in range(len(class1)):
			if (class1[c]):
				numClass1 += 1
				Class1Clients[i].append(UEs[i][c])
		ReuseZoneSize[i] = numClass1
								
	for i in G:
		G[i] = list(G[i])
	for i in G_Ue:
		G_Ue[i] = list(G_Ue[i])

	return (G,G_Ue,ues_intefering,intf_edges,ReuseZoneSize,Class1Clients)

# Find interference conflicts between eNodeBs and UEs
# Similar to IQhoppingShares - see comments there

def GenerateDirGraph_CQI(Distance,RxPower,CQIVals,enb_coord):

	n = len(enb_coord)
	
	G_dir = {} # directed graph to run sequential assignments

	for i in range(n):
		G_dir[i] = set()
	
	for i in range(n):
		rcv_pow = RxPower[i]
		dists = Distance[i]
		CQIs = CQIVals[i]
		Actual_CQI = []

		#######################################
		# interference map generation based on reported CQI values
		# iterates over all eNBs
		for j in range(n):
			if (i==j):
				continue
			#iterates over all the UEs for eNB j
			for z in range(len(rcv_pow)):
				u_pow = rcv_pow[z]
				u_dist = dists[z]
				u_cqi = CQIs[z]
				
				if (u_cqi[i] - u_cqi[j] > cqi_thresh):
				#if (u_pow[j]>thresh):
				#if (u_dist[j]<dist_thresh):

					G_dir[i].add(j)
				
		#######################################

	for i in G_dir:
		G_dir[i] = list(G_dir[i])

	return (G_dir)

def GenerateGraph_SE(Distance,RxPower,CQIVals,enb_coord,UEs,u_m):
	n = len(enb_coord)
	G_S_E = {} # inteference graph based on spectral efficiency difference
	ReuseZoneSize_S_E = {}
	Class1Clients = {}
	
	ues_intefering = {}
	intf_edges = {}
	
	cqi_debug = open('CQIDebug.txt','w')

	for i in range(n):
		G_S_E[i] = set()

		intf_edges[i]=set()
		ues_intefering[i] = len(UEs[i])
		Class1Clients[i] = []

	for i in range(n):
		rcv_pow = RxPower[i]
		dists = Distance[i]
		CQIs = CQIVals[i]
		Actual_CQI = []

		class1_S_E = [True for a in range(len(UEs[i]))]		


		ue_num = 0
		for u in UEs[i]:
			(powers,distances) = rcvd_pows(enb_coord,u)
			rcv_pow.append(powers)
			dists.append(distances)

			snr = {}
			cqis = {}
			for k in range(n):
				if(i==k):
					snr[k] = powers[i]/(Noise*spectrumBW[bw])
				else :
					snr[k] = powers[i]/(powers[k] + Noise*spectrumBW[bw])
				cqis[k] = getCQI(getSpectralEfficiency(snr[k]))
			CQIs.append(cqis)

			#############################
			# generation of interference map according to the FERMI paper 
			# i.e. percentage differe throughput per block must be less than 25% 
			intf_a = sum(powers.values())-powers[i]+Noise*spectrumBW[bw]

			s_a = powers[i]/ intf_a
			cqi_a = getCQI(getSpectralEfficiency(s_a))

			if ((SpectralEfficiencyForCqi[cqis[i]] - SpectralEfficiencyForCqi[cqi_a]) / SpectralEfficiencyForCqi[cqis[i]] > 0.25) :
				class1_S_E[ue_num] = False; # if UE interfering with some other eNB its class 2

			Actual_CQI.append((cqis[i],cqi_a))
			if (not class1_S_E[ue_num]):
				interferers = []
				for m in range(n):
					if not (m == i):
						interferers.append(m)
				while ((SpectralEfficiencyForCqi[cqis[i]] - SpectralEfficiencyForCqi[cqi_a]) / SpectralEfficiencyForCqi[cqis[i]] > 0.25):
					strongest = max_interferer(interferers,powers)
					intf_a -= powers[strongest]
					interferers.remove(strongest)

					s_a = powers[i]/ intf_a
					cqi_a = getCQI(getSpectralEfficiency(s_a))

					G_S_E[i].add(strongest)
					G_S_E[strongest].add(i)

					ues_intefering[strongest] += 1
					intf_edges[strongest].add(UEs[i][ue_num])


			ue_num += 1
			################################


		# class 1 clients for eNB i

		numClass1 = 0
		for c in range(len(class1_S_E)):
			if (class1_S_E[c]):
				numClass1 += 1
				Class1Clients[i].append(UEs[i][c])
		ReuseZoneSize_S_E[i] = numClass1

		
		cqi_debug.write(str(i)+':')
		cqi_debug.write(str(Actual_CQI) + '\n')
		for cqi in CQIs:
			cqi_debug.write(str(cqi))
			cqi_debug.write('\n')
								
	for i in G_S_E:
		G_S_E[i] = list(G_S_E[i])

	cqi_debug.close()


	return (G_S_E,ues_intefering,intf_edges,ReuseZoneSize,Class1Clients)

# Find connected components of a graph
def connected_graphs(G):
	conn_G = []
	v = list(G.keys())
	q = Q.Queue()
	while (len(v)!=0):
		q.put(v[0])
		g = []
		while q.empty() == False:
			curr = q.get()
			if curr in v:
				g.append(curr)
				v.remove(curr)
			for n in G[curr]:
				if n in v:
					q.put(n)
					v.remove(n)
					g.append(n)
		G_dash = {}
		for i in g:
			G_dash[i] = deepcopy(G[i])
		conn_G.append(G_dash)
	return conn_G



# Find contiguous blocks of channels and make their represenation compact
def make_blks(A):
	if (len(A) == 0):
		return([(-1,-1)])
	tmp = deepcopy(A)
	lo = min(tmp)
	hi = lo
	tmp.remove(lo)
	A_ = []
	while (len(tmp)>0):
		new_hi = min(tmp)
		if (new_hi != hi+1):
			A_.append((lo,hi))
			lo = new_hi
			hi = lo
		else:
			hi = new_hi
		tmp.remove(hi)
	A_.append((lo,hi))
	return A_



def getSpectralEfficiency(snr):
	BER = 0.00005
	return math.log(1+(snr/(-math.log(5*BER)/1.5)),2)

def getCQI(s):
	cqi = 0
	for i in range(len(SpectralEfficiencyForCqi)):
		if(s>SpectralEfficiencyForCqi[i]):
			cqi+=1
		else:
			break
	return (cqi-1)


# Find expected sub-band share per eNodeB for IQ-Hopping
def IQhoppingShares(Distance,RxPower,CQIVals,N,enb_coord,UEs):
	n = len(enb_coord)
	ues_intefering = {}

	interferingENB = {}

	ReuzeClients = {}

	intf_edges = {}

	ctrlIntfrnc = {}
	for i in range(n):
		dists = Distance[i]
		k = 0
		for u in UEs[i]:
			u_dist = dists[k]
			minforue = 1000
			for j in range(n):
				if (i==j):
					continue
				if (u_dist[j] < minforue ):
					minforue = u_dist[j]
			if float(u_dist[i])/minforue < 1.5 :
				ctrlIntfrnc[u] = 0.85
			elif float(u_dist[i])/minforue < 4 :
				ctrlIntfrnc[u] = 0.40
			else :
				#print u_dist
				#print u_dist[i] , minforue
				ctrlIntfrnc[u] = 0.0
			k+=1
				



	for i in range(n):
		ues_intefering[i] = len(UEs[i])
		interferingENB[i] = set() 
		ReuzeClients[i] = set()
		intf_edges[i]=set()

	for i in range(n):
		rcv_pow = RxPower[i]
		dists = Distance[i]
		CQIs = CQIVals[i]

                # Count number of UEs that eNodeB <j> interfers with (ues_interfering)
		for j in range(n):
			if (i==j):
				continue
			#for u_pow in rcv_pow:
			for z in range(len(rcv_pow)):

				prachRange = getOverHearingRange(enb_coord[i],UEs[i][z])
				#print(prachRange)

				u_pow = rcv_pow[z]
				u_dist = dists[z]
				u_cqi = CQIs[z]
				#if (u_pow[j]>thresh):
				if (u_dist[j]<prachRange):
				#if (u_cqi[i] - u_cqi[j] > cqi_thresh):
					ues_intefering[j] += 1
					interferingENB[j].add(i)
					intf_edges[j].add(UEs[i][z])

				ActualSNR = u_pow[i]/(Noise*spectrumBW[bw] + (sum(u_pow.values()) - u_pow[i]) )
				ActualCQI = getCQI(getSpectralEfficiency(ActualSNR))
				#print ActualCQI
				if (u_cqi[i] - ActualSNR <= cqi_thresh):
					ReuzeClients[i].add(z)

	#print ReuzeClients
	

        # Compute the share = #my_UEs / (#my_UEs + #interfering_UEs)
	W = {}
	W_1 = {}
	W_2 = {}
	W_Agressive = {}
	for enb in ues_intefering:
		W[enb] = int(math.floor((N*float(len(UEs[enb]))/ues_intefering[enb]) + 0.5))


		R = len(ReuzeClients[enb])
		I = len(UEs[enb]) - R
		O = ues_intefering[enb] - R - I

		Rz = int (N*(float(R)/(R+I))+0.5 )
		if (I+O != 0):
			W_1[enb] = Rz + int((N-Rz)*(float(I)/(I+O)) + 0.5 ) 
		else :
			W_1[enb] = Rz

		Rz = int (N*(float(R)/(R+max(I,O)))+0.5 )
		if (I+O != 0):
			W_2[enb] = Rz + int((N-Rz)*(float(I)/(I+O)) + 0.5 ) 
		else :
			W_2[enb] = Rz


		other_UEs = float(ues_intefering[enb] - len(UEs[enb]))
		if (len(interferingENB[enb])>0):
			scaled_other = other_UEs/len(interferingENB[enb])
		else:
			scaled_other = other_UEs

		W_Agressive[enb] = int(math.floor(N*len(UEs[enb])/(scaled_other+len(UEs[enb])) + 0.5))	

	#print W
	#print W_1
	#print W_2

	return (W,W_1,W_2,W_Agressive,ctrlIntfrnc,intf_edges)



def fairness_index(shares):
# Fairness index calculation
	vals = list(shares.values())
	num_vals = len(vals)
	if 0 in vals:
		log_vals = []
	else:
		log_vals = [math.log10(i) for i in vals]
	sq_vals = [i**2 for i in vals]
	jains_index = float((sum(vals)**2))/(num_vals*(sum(sq_vals)))
	jains_index = eval("%.2f" % jains_index)
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

def remove_node(e,c_map):
	if (e in c_map):
		for n in c_map:
			if e in c_map[n]:
				c_map[n].remove(e)
		del c_map[e]

def FermiAllocations(UEs,u_m,G,load,N,info,comp):
		assign = []
		alloc = []
		assign_us = []


                # Main calculation
		for c_map in G:
			(channel_assignment,allocated_share) = Fermi(c_map,load,N)
			assign.append(channel_assignment)
			alloc.append(allocated_share)
			#assign_us.append(KlogN(c_map,W,N))

		Assign = {}
		for a in assign:
			Assign.update(a)
		Allocated_share = {}
		for a in Assign:
			share = 0
			for interval in Assign[a]:
				share += interval[1]-interval[0]+1
			Allocated_share[a] = share
		
		share_UE={}
		for a in Allocated_share:
			share = eval("%.2f" % (float(Allocated_share[a])/load[a]))
			for u in UEs[a]:
				share_UE[u_m[u]]=share

		writeInfo('FERMI',Assign,Allocated_share,share_UE,info)

		return (Assign,Allocated_share)

def FermiAllocationsDir(UEs,u_m,G,load,N,info,comp,OpReuseSizes,Class1Clients):
		#print ("DIR")
		assign = []
		alloc = []
		assign_us = []
		class1ClientWeight = {}

		for enb in UEs:
			DesiredSize = int(math.floor(len(Class1Clients[enb])/float(load[enb]) * N + 0.5))

			LeftToAllot = OpReuseSizes[enb]
			for u in Class1Clients[enb]:
				DesiredPerUe = float(DesiredSize)/len(Class1Clients[enb])
				#print (enb,DesiredPerUe,LeftToAllot)
				if (LeftToAllot >= DesiredPerUe):
					class1ClientWeight[u] = 0
					LeftToAllot -= DesiredPerUe
				elif (LeftToAllot > 0):
					class1ClientWeight[u] = 1 - LeftToAllot/DesiredPerUe
					LeftToAllot = 0
				else :
					class1ClientWeight[u] = 1.0
					
				#class1ClientWeight[u] = float(DesiredSize - OpReuseSizes[enb])/(len(Class1Clients[enb]))
			

		load={}
		IsClass1 = {}
		RZSizes = {}

		k=0
		for enb in UEs:
			for user in UEs[enb]:
				if (user in Class1Clients[enb]):
					load[k] = class1ClientWeight[user]
					IsClass1[k] = True
				else:
					load[k]= 1.0
					IsClass1[k] = False
				RZSizes[k] = OpReuseSizes[enb]
				k+=1
	
		#print load
				

		#print load
                # Main calculation
		for c_map in G:
			(channel_assignment,allocated_share) = Fermi_UE(c_map,load,N,IsClass1,RZSizes)
			assign.append(channel_assignment)
			alloc.append(allocated_share)
			#assign_us.append(KlogN(c_map,W,N))

		Assign = {}
		for a in assign:
			Assign.update(a)

		Assign

		Assign_enb = {}
		u = 0

		share_UE={}

		for enb in UEs:
			Assign_enb[enb] = []
			for user in UEs[enb]:
				share = 0
				for interval in Assign[u]:
					if interval not in Assign_enb[enb]:
						Assign_enb[enb].append(interval)
						share += interval[1]-interval[0]+1
				share_UE[u] = share
				u+=1
		#print (Assign_enb)
			
		Allocated_share = {}
		for a in Assign_enb:
			share = 0
			for interval in Assign_enb[a]:
				share += interval[1]-interval[0]+1
			Allocated_share[a] = share
		
		#print(Allocated_share)
		'''		
		share_UE={}
		for a in Allocated_share:
			share = eval("%.2f" % (float(Allocated_share[a])/load[a]))
			for u in UEs[a]:
				share_UE[u_m[u]]=share
		'''

		writeInfo('FERMI-Dir',Assign_enb,Allocated_share,share_UE,info)
		#print ("DIR END")
		return (Assign_enb,Allocated_share,share_UE)

def FermiAllocationsRZ(UEs,u_m,G,load,N,info,comp,OpReuseSize):
		assign = []
		alloc = []
		assign_us = []

		Allocated_share = {}
		Assign = {}

		#print OpReuseSize

		for e in load:
			if (load[e] <= 0):
				#print e
				for c_map in G:
					remove_node(e,c_map)
				Allocated_share[e] = 0 + OpReuseSize[e]
				Assign[e] = [(0,OpReuseSize[e]-1)]
				#Allocated_share[e] = min(0 + OpReuseSize[e],N)

                # Main calculation
		for c_map in G:
			if (len(c_map) == 0):
				continue
			(channel_assignment,allocated_share) = Fermi_RZ(c_map,load,N,OpReuseSize)
			assign.append(channel_assignment)
			alloc.append(allocated_share)
			#assign_us.append(KlogN(c_map,W,N))


		for a in assign:
			Assign.update(a)

		for a in Assign:
			share = 0
			for interval in Assign[a]:
				share += interval[1]-interval[0]+1
			Allocated_share[a] = share
			#Allocated_share[a] = min(share + OpReuseSize[a],N)
		
		share_UE={}
		#for a in Allocated_share:
		#	share = eval("%.2f" % (float(Allocated_share[a])/load[a]))
		for u in UEs[a]:
			share_UE[u_m[u]]=1

		writeInfo('FERMI Reuse Zone',Assign,Allocated_share,share_UE,info)

		return (Assign,Allocated_share)

def FermiAssignIQAlloc(UEs,u_m,G,load,N,info,comp,OpReuseSizes,IQshare,Fermishare,FermishareUE,Class1Clients):


		#print ("DIR")
		assign = []
		alloc = []
		assign_us = []
		class1ClientWeight = {}

		for enb in UEs:
			DesiredSize = int(math.floor(len(Class1Clients[enb])/float(load[enb]) * N + 0.5))
			for u in Class1Clients[enb]:
				class1ClientWeight[u] = float(DesiredSize - OpReuseSizes[enb])/(N*len(Class1Clients[enb]))
			


		Shares = {}
		ShareAvailable = {}
		for enb in UEs:
			ShareAvailable[enb] = min(Fermishare[enb],IQshare[enb])
			OpReuseSizes[enb] = min(OpReuseSizes[enb],ShareAvailable[enb])
			if (Fermishare[enb] <= IQshare[enb]):
				for u in UEs[enb]:
					Shares[u_m[u]] = FermishareUE[u_m[u]]
			else :
				for u in UEs[enb]:
					if u not in Class1Clients[enb]:
						Shares[u_m[u]] = min(ShareAvailable[enb],FermishareUE[u_m[u]])
						ShareAvailable[enb] -= Shares[u_m[u]]
						OpReuseSizes[enb] -= Shares[u_m[u]]

				for u in UEs[enb]:
					if u in Class1Clients[enb]:
						Shares[u_m[u]] = ShareAvailable[enb]
			OpReuseSizes[enb] = max(0,OpReuseSizes[enb])


		load={}
		IsClass1 = {}
		RZSizes = {}

		k=0
		for enb in UEs:
			for user in UEs[enb]:
				if (user in Class1Clients[enb]):
					load[k] = class1ClientWeight[user]
					IsClass1[k] = True
				else:
					load[k]= 1.0
					IsClass1[k] = False
				RZSizes[k] = OpReuseSizes[enb]
				k+=1

		#print load
				

		#print load
                # Main calculation
		for c_map in G:
			(channel_assignment,allocated_share) = Fermi_Pre_Alloc(c_map,load,N,IsClass1,RZSizes,Shares)
			assign.append(channel_assignment)
			alloc.append(allocated_share)
			#assign_us.append(KlogN(c_map,W,N))

		Assign = {}
		for a in assign:
			Assign.update(a)

		Assign

		Assign_enb = {}
		u = 0

		share_UE={}

		for enb in UEs:
			Assign_enb[enb] = []
			for user in UEs[enb]:
				share = 0
				for interval in Assign[u]:
					if interval not in Assign_enb[enb]:
						Assign_enb[enb].append(interval)
						share += interval[1]-interval[0]+1
				share_UE[u] = share
				u+=1
		#print (Assign_enb)
			
		Allocated_share = {}
		for a in Assign_enb:
			share = 0
			for interval in Assign_enb[a]:
				share += interval[1]-interval[0]+1
			Allocated_share[a] = share
		
		#print(Allocated_share)
		'''		
		share_UE={}
		for a in Allocated_share:
			share = eval("%.2f" % (float(Allocated_share[a])/load[a]))
			for u in UEs[a]:
				share_UE[u_m[u]]=share
		'''

		writeInfo('FERMI-IQ',Assign_enb,Allocated_share,share_UE,info)
		#print ("DIR END")
		return (Assign_enb,Allocated_share)


def AdjustLoads(i_map,load,reuseClients,N):
	OpReuseSize = {}
	IsolatedLoads = {}
	for eNB in reuseClients:
		size = float(reuseClients[eNB])/load[eNB] * N
		desired = float(math.floor(size+0.5))
		for eNB_other in i_map[eNB]:
			size_other = float(reuseClients[eNB_other])/load[eNB_other] * N
			if (size_other < size):
				size = size_other
		OpReuseSize[eNB] = int(math.floor(size+0.5))
		if (OpReuseSize[eNB] > 0):
			IsolatedLoads[eNB] = int(math.floor(load[eNB] - reuseClients[eNB] + (1-OpReuseSize[eNB]/desired)*reuseClients[eNB] + 0.5))
		else:
			IsolatedLoads[eNB] = load[eNB]
	return (OpReuseSize,IsolatedLoads)
		



def SequentialAllocations(UEs,u_m,ue_interfering,edges,G,i_map_dir,load,N,info,comp):
		W = {}
		for enb in ue_interfering:
			W[enb] = int(math.floor((N*float(len(UEs[enb]))/ue_interfering[enb]) + 0.5))

                # Main calculation
		(Assign_us,ue_Assign_us) = KlogN(i_map_dir,UEs,W,edges,N)

		Assign_mis = {}
		for a in Assign_us:
			Assign_mis[a]=make_blks(Assign_us[a])

		tmp={}
		for a in ue_Assign_us:
			tmp[u_m[a]]=ue_Assign_us[a]
		
		Allocated_share_us = {}
		for a in Assign_us:
			Allocated_share_us[a]=len(Assign_us[a])

		share_UE = {}
		for a in ue_Assign_us:
			share_UE[u_m[a]]=len(ue_Assign_us[a])

		writeInfo('SEQUENTIAL',Assign_mis,Allocated_share_us,share_UE,info)

		#info.write('Assignment per UE\n')
		#info.write(str(tmp))
		#info.write('\n')
		#info.write('DESIRED SHARE\n')
		#info.write(str(W))
		#info.write('\n')

		return (Assign_mis,Allocated_share_us,W)

def SequentialAllocationsBackOff(UEs,u_m,ue_interfering,edges,G,i_map_dir,load,N,info,comp):
		W = {}
		for enb in ue_interfering:
			W[enb] = int(math.floor((N*float(len(UEs[enb]))/ue_interfering[enb]) + 0.5))
		#(Assign_us,ue_Assign_us) = KlogN(i_map_dir,UEs,W,edges,N)		

                # Main calculation
		(Assign_us2,ue_Assign_us2) = KlogN_Backoff(i_map_dir,UEs,W,edges,N)

		Assign_mis_bo = {}
		for a in Assign_us2:
			Assign_mis_bo[a]=make_blks(Assign_us2[a])

		tmp={}
		for a in ue_Assign_us2:
			tmp[u_m[a]]=ue_Assign_us2[a]

		Allocated_share_us = {}
		for a in Assign_us2:
			Allocated_share_us[a]=len(Assign_us2[a])

		share_UE = {}
		for a in ue_Assign_us2:
			share_UE[u_m[a]]=len(ue_Assign_us2[a])

		writeInfo('SEQUENTIAL-BACKOFF',Assign_mis_bo,Allocated_share_us,share_UE,info)
		#info.write('Assignment per UE\n')
		#info.write(str(tmp))
		#info.write('\n')
		#info.write('DESIRED SHARE\n')
		#info.write(str(W))
		#info.write('\n')
		
		return (Assign_mis_bo,Allocated_share_us)

def SequentialAllocationsFermiShare(UEs,u_m,ue_interfering,edges,G,i_map_dir,load,N,info,comp,FERMIshare):
		W = FERMIshare	

                # Main calculation
		(Assign_us2,ue_Assign_us2) = KlogN_Backoff(i_map_dir,UEs,W,edges,N)

		Assign_mis_bo = {}
		for a in Assign_us2:
			Assign_mis_bo[a]=make_blks(Assign_us2[a])

		tmp={}
		for a in ue_Assign_us2:
			tmp[u_m[a]]=ue_Assign_us2[a]

		Allocated_share_us = {}
		for a in Assign_us2:
			Allocated_share_us[a]=len(Assign_us2[a])

		share_UE = {}
		for a in ue_Assign_us2:
			share_UE[u_m[a]]=len(ue_Assign_us2[a])

		writeInfo('SEQUENTIAL-BACKOFF-FERMI-SHARE',Assign_mis_bo,Allocated_share_us,share_UE,info)
		#info.write('Assignment per UE\n')
		#info.write(str(tmp))
		#info.write('\n')
		#info.write('DESIRED SHARE\n')
		#info.write(str(W))
		#info.write('\n')
		
		return (Assign_mis_bo)

def IQhopping(Distance,RxPower,CQIVals,N,enb_coord,UEs,info):

	(IQshares,IQshares_1,IQshares_2,IQAggressiveShares,ctrlIntfrnc,intf_edges) = IQhoppingShares(Distance,RxPower,CQIVals,N,enb_coord,UEs)

	info.write('IQHopping:\n')
	info.write('Share per ENB: '+str(fairness_index(IQshares))+'\n')
	info.write(str(IQshares))
	info.write('\n')
	info.write('\n')

	info.write('IQHopping_1:\n')
	info.write('Share per ENB: '+str(fairness_index(IQshares))+'\n')
	info.write(str(IQshares_1))
	info.write('\n')
	info.write('\n')

	info.write('IQHopping_2:\n')
	info.write('Share per ENB: '+str(fairness_index(IQshares))+'\n')
	info.write(str(IQshares_2))
	info.write('\n')
	info.write('\n')

	info.write('IQHopping Aggressive:\n')
	info.write('Share per ENB: '+str(fairness_index(IQAggressiveShares))+'\n')
	info.write(str(IQAggressiveShares))
	info.write('\n')
	info.write('\n')


	info.write('Ctrl Interference:\n')
	info.write(str(ctrlIntfrnc))
	info.write('\n')
	info.write('\n')

	return (IQshares,IQshares_1,IQshares_2,IQAggressiveShares,ctrlIntfrnc,intf_edges)

def plot_graph(filename,i_map,enb_coord,u_m,UEs):
	#plot topology using Matlab/Octave
	pf = open(outputDir + filename + '.m','w')

	num = len(enb_coord)
	pf.write('clr = hsv('+str(num)+');\n')
	#pf.write("figure('Position',[-10,-10,"+str(l+10)+","+str(w+10)+"]);\n")
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
	pf.write('print '+filename+'.jpg;\n')

	pf.write('hold off;\n') 
	pf.close()


def plot_ue_interference(filename,edges,enb_coord,u_m,UEs):
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
	pf.write('print '+filename+'.jpg;\n')
	pf.write('hold off;\n') 
	pf.write('figure;\n') 
	pf.close()


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


def operatorNetwork(num,l,w):
	center = (l/2,w/2)
	ENB = []
	UEs = {}
	load = {}
	n = 0
	for k in range(num):
		UEAssoc = {}
		enb_coord = []
		origin = (random.randint(center[0]-50,center[0]+50),random.randint(center[1]-50,center[1]+50))
		r = 300
		q = random.uniform(0, 2*math.pi)
		for i in range(4):
			x = origin[0] + r*math.cos(q)
			y = origin[1] + r*math.sin(q)
			enb_coord.append ((x,y))
			UEAssoc[(x,y)] = []
			for j in range(6):
				ue = gen_ue_coord((x,y),l,w,20,150)
				UEAssoc[(x,y)].append(ue)
			q += math.pi/2

		'''
		for j in range(21):
			x = random.randint(50,l-50)
			y = random.randint(50,w-50)

			ue gen_ue_coord(enb_coord[enb],l,w,20,150)
			(t, dists) = rcvd_pows(enb_coord,(x,y))
			min_dist = min(dists.values())
			for d in dists:
				if dists[d] == min_dist:
					UEAssoc[enb_coord[d]].append((x,y))
		'''
		for e in enb_coord:
			ENB.append(e)
			UEs[n] = UEAssoc[e]
			load[n] = len(UEAssoc[e])
			n+=1
	print(load)

	return (ENB,UEs,load)
	


def main(density,l,w,toytop):

        # number of sub-channels
	N = int(math.ceil(bw / rbgsize))

        # size of the grid (in meters)
	#l = 500
	#w = 500

        # number of eNodeBs
	n = int(math.floor((l*w)*density + 0.5))
	print(n)
	enb_coord = gen_eNbs_coord(n,l,w)

        # array [n] with number of UEs per eNodeB
	load = {}
	for i in range(n):
		load[i] = 3

	no_class1_clients_in_topo = 3

	UEs = {}
	ue_list_to_print={}

	for enb in range(len(load)):
		ues = []
		for j in range(load[enb]):
			#if (j==0 and no_class1_clients_in_topo > 0):
			#	ue = gen_ue_coord(enb_coord[enb],l,w,10,15)
			#	no_class1_clients_in_topo -= 1
			if (j < 1):
				ue = gen_ue_coord(enb_coord[enb],l,w,10,20)
			else:
				ue = gen_ue_coord(enb_coord[enb],l,w,50,150)
			ues.append(ue)

		UEs[enb] = ues

	On = [0.1]
	Off = [0.1]

	#(enb_coord,UEs,load) = operatorNetwork(3,800,800)

	'''
	# single client topologies
	if (toytop == 1):
		#topo 1
		enb_coord = [(100,100),(400,100),(700,100)]
		UEs={0: [(250,100)], 1: [(550,110)], 2: [(550, 90)]}
		load={0: 1,1: 1, 2: 1}
		On = {0: [0.2, 0.1],1: [0.2, 0.1],2: [0.2, 0.1]}
		Off = {0: [0, 0.1],1: [0, 0.1],2: [0, 0.1]}
	elif (toytop == 2):
		#topo 2
		enb_coord = [(100,100),(100,400),(400,400),(400,100)]
		UEs={0: [(90,200)], 1: [(250,410)], 2: [(410, 250)], 3: [(250, 90)]}
		load={0: 1,1: 1, 2: 1, 3: 1}
		On = {0: [0.2, 0.1],1: [0.2, 0.1],2: [0.2, 0.1],3: [0.2, 0.1]}
		Off = {0: [0, 0.1],1: [0, 0.1],2: [0, 0.1],3: [0, 0.1]}
	elif (toytop == 3):
		#topo 3
		enb_coord = [(100,100),(100,400),(400,400),(400,100)]
		UEs={0: [(90,200)], 1: [(250,410)], 2: [(410, 250)], 3: [(250, 250)]}
		load={0: 1,1: 1, 2: 1, 3: 1}
		On = {0: [0.2, 0.1],1: [0.2, 0.1],2: [0.2, 0.1],3: [0.2, 0.1]}
		Off = {0: [0, 0.1],1: [0, 0.1],2: [0, 0.1],3: [0, 0.1]}
	elif (toytop == 4):
		#topo 4
		enb_coord = [(100,100),(400,100),(700,100),(1000,100)]
		UEs={0: [(250,100)], 1: [(550,110)], 2: [(550, 90)], 3: [(850, 90)]}
		load={0: 1,1: 1, 2: 1, 3: 1}
		On = {0: [0.2, 0.1],1: [0.2, 0.1],2: [0.2, 0.1],3: [0.2, 0.1]}
		Off = {0: [0, 0.1],1: [0, 0.1],2: [0, 0.1],3: [0, 0.1]}
	elif (toytop == 5):
		#topo 5
		enb_coord = [(100,400),(400,400),(625,600),(625,200)]
		UEs={0: [(250,400)], 1: [(550,550)], 2: [(625, 400)], 3: [(550, 300)]}
		load={0: 1,1: 1, 2: 1, 3: 1}
		On = {0: [0.2, 0.1],1: [0.2, 0.1],2: [0.2, 0.1],3: [0.2, 0.1]}
		Off = {0: [0, 0.1],1: [0, 0.1],2: [0, 0.1],3: [0, 0.1]}

	'''
	#Topologies in FERMI paper
	if (toytop == 1):
		#topo 1
		#enb_coord = [(100,100),(400,100),(700,100)]
		#enb_coord = enb_coord = [(100,100),(400,100),(700,100)]
		#UEs={0: [(100,150) , (250,100)], 1: [(400, 150), (550,110)], 2: [(700,150), (550, 90)]}
		#UEs={0: [(250,100)], 1: [(550,110)], 2: [(550, 90)]}
		#load={0: 2,1: 2, 2: 2}
		enb_coord = [(100,100),(300,100),(200,0)]
		UEs={0: [(190,95),(190,105),(190,100)], 1: [(210, 100)], 2: [(200,90)]}
		load={0: 3,1: 1,2:1}
		#load={0: 1,1: 1, 2: 1}
		On = {0: [0.2, 0.1],1: [0.2, 0.1],2: [0.2, 0.1]}
		Off = {0: [0, 0.1],1: [0, 0.1],2: [0, 0.1]}
	elif (toytop == 2):
		#topo 2
		'''
		enb_coord = [(100,100),(100,400),(400,400),(400,100)]
		UEs={0: [(100,50) , (90,200)], 1: [(100, 450), (250,410)], 2: [(400,450), (410, 250)], 3: [(400,50), (250, 90)]}
		load={0: 2,1: 2, 2: 2, 3: 2}
		'''
		enb_coord = [(100,100),(500,500)]
		UEs={0: [(150,100), (100, 150)],1: [(550,500)]}
		load={0: 2, 1: 1}
		On = {0: [0.2, 0.1],1: [0.2, 0.1],2: [0.2, 0.1],3: [0.2, 0.1]}
		Off = {0: [0, 0.1],1: [0, 0.1],2: [0, 0.1],3: [0, 0.1]}
	elif (toytop == 3):
		#topo 3
		enb_coord = [(100,100),(100,400),(400,400),(400,100)]
		UEs={0: [(100,50) , (90,200)], 1: [(100, 450), (250,410)], 2: [(400,450), (410, 250)], 3: [(400,50), (250, 250)]}
		load={0: 2,1: 2, 2: 2, 3: 2}
		On = {0: [0.2, 0.1],1: [0.2, 0.1],2: [0.2, 0.1],3: [0.2, 0.1]}
		Off = {0: [0, 0.1],1: [0, 0.1],2: [0, 0.1],3: [0, 0.1]}
	elif (toytop == 4):
		#topo 4
		enb_coord = [(100,100),(400,100),(700,100),(1000,100)]
		UEs={0: [(100,150) , (250,100)], 1: [(400, 150), (550,110)], 2: [(700,150), (550, 90)], 3: [(1000,150), (850, 90)]}
		load={0: 2,1: 2, 2: 2, 3: 2}
		On = {0: [0.2, 0.1],1: [0.2, 0.1],2: [0.2, 0.1],3: [0.2, 0.1]}
		Off = {0: [0, 0.1],1: [0, 0.1],2: [0, 0.1],3: [0, 0.1]}
	elif (toytop == 5):
		#topo 5
		#enb_coord = [(100,400),(400,400),(625,600),(625,200)]
		#UEs={0: [(50,400) , (250,400)], 1: [(400, 350), (550,550)], 2: [(675,600), (625, 400)], 3: [(675,200), (550, 300)]}
		#load={0: 2,1: 2, 2: 2, 3: 2}
		enb_coord = [(400,600),(600,600),(590,500),(200,600)]
		UEs={0: [(500,610) , (300,600)], 1: [(500, 590)], 2: [(690,500)], 3: [(100,600)]}
		load={0: 2,1: 1, 2: 1, 3: 1}
		On = {0: [0.2, 0.1],1: [0.2, 0.1],2: [0.2, 0.1],3: [0.2, 0.1]}
		Off = {0: [0, 0.1],1: [0, 0.1],2: [0, 0.1],3: [0, 0.1]}

	'''
	enb_coord = [(100,0),(200,0)]
	UEs={0: [(120,0)], 1:[(210,5),(210,10),(210,15)]}
	load={0: 1, 1: 3 }
	'''
	'''
	enb_coord = [(556, 352), (672, 711), (589, 664), (509, 602), (387, 249), (298, 244), (434, 703), (779, 669), (106, 583), (629, 157), (136, 785), (132, 736), (636, 192), (555, 169)]
	UEs = {0: [(567.0058472214033, 481.37156639551495), (520.9393371950939, 336.86220828473955), (520.593464670609, 332.63185262736033), (554.7527864306933, 203.63429084930777), (634.4303212817989, 331.52904976965453), (594.428533191026, 407.5636954496506)], 1: [(708.4235710744688, 743.5369868822962), (690.6201581038623, 726.515464488296), (614.3753113477153, 712.4760570490323), (701.2608276038981, 698.3962456590011), (655.8688712626627, 676.179781766384), (584.8381482193378, 646.5102308051949)], 2: [(589.5722183090157, 738.1049804052884), (543.9149425503529, 713.3270596676714), (693.2870551604027, 583.7343019506682), (465.3832888092391, 646.1084488155283), (480.0187353788513, 652.988256578613), (679.7494099411015, 748.5614991076931)], 3: [(447.10620190618073, 482.6872052661287), (529.1804453183474, 558.413437647292), (394.6592368830069, 530.6702869524715), (610.025306452538, 574.1487532650397), (535.1659621401599, 532.928999223494), (470.96760380902487, 686.3610135623658)], 4: [(378.28137664262124, 270.07651249987816), (412.36716347878166, 171.44242817210647), (527.6812946936609, 298.38640315834874), (445.3106685200999, 149.42747454628466), (351.7095039344982, 223.43433796005056), (396.8181571858327, 310.21199717299385)], 5: [(185.84553757060212, 159.07376068482603), (191.06242737028208, 167.74489900647205), (363.36289935904426, 206.16385359543352), (283.5656395933625, 294.28682289922364), (425.7432505032069, 248.4035092364313), (416.0614908809715, 219.4541648457248)], 6: [(422.1582673346475, 752.5807742240354), (474.3681972194192, 685.2581403682973), (360.15169249173147, 754.4604322089085), (342.59111178190875, 673.8057257536167), (363.2118323965816, 753.4265392806343), (383.6145883776824, 647.6811768561874)], 7: [(713.9112127877702, 682.1134729013031), (710.8920056854697, 696.202234763561), (779.1576352498665, 691.0704934778063), (681.4853772200358, 754.886980894657), (779.1108503080255, 635.8161327900546), (643.9666376550286, 660.9111144740939)], 8: [(9.349143626202604, 498.8950883705866), (186.705047921371, 582.6124819975195), (140.0612071807089, 492.2121758593375), (251.61144196736305, 600.8990305750943), (93.43157457372118, 627.4025445408697), (83.23438775321748, 649.7529521537581)], 9: [(637.6510950602419, 180.16462252068612), (619.5702242665637, 229.6559066599279), (619.74125447848, 180.38081910957723), (521.83398477942, 127.63228647833219), (714.5820485896965, 175.88018327469092), (659.5363237138677, 178.62186755919325)], 10: [(63.42897713185157, 734.6891911178855), (141.36091545809208, 714.7244048684037), (203.43774803572194, 786.1404285736185), (83.86082901449473, 661.0231478694357), (41.120593068394584, 796.5386780750434), (223.5268888536554, 738.1905945290034)], 11: [(202.0713386602543, 645.5985404137346), (47.43927277037018, 736.2776999288467), (145.41751130444385, 755.982636903769), (74.19634519478006, 687.1275720983392), (240.05587643119247, 757.1982378583383), (45.88631923206498, 653.2343253945794)], 12: [(555.6905672213386, 214.8642591280897), (579.3884033961027, 227.9559483736894), (528.8087590114715, 104.61391863839873), (586.8450674371186, 176.6617585148029), (702.6398821386954, 254.1053129606898), (597.1585580860813, 194.26510767868598)], 13: [(593.4692102053252, 154.39210995575903), (584.1569002223373, 221.29464726264305), (444.61225024522184, 170.24553115814362), (604.6257226839488, 303.65154646675353), (679.5956235475355, 91.23070094376554), (638.9799997062987, 185.704327268724)]}
	load = {0: 6, 1: 6, 2: 6, 3: 6, 4: 6, 5: 6, 6: 6, 7: 6, 8: 6, 9: 6, 10: 6, 11: 6, 12: 6, 13: 6}
	'''

	#####################################################
	# u_m: maps UE coordinates to ue number
	i=0
	for n in UEs:
		for u in UEs[n]:
			ue_list_to_print[i]=u
			#print(rcvd_pows(enb_coord,u))
			i+=1
	#reversing
	u_m = {}
	for u in ue_list_to_print:
		u_m[ue_list_to_print[u]] = u
	######################################################

	#write info to file
	info = open(outputDir + 'info.txt','w')
	comp = open(outputDir + 'comparison.txt','a')

	#writing variables to file to recreate experiment
	info.write('\n')
	info.write('\n')
	info.write('\n')
	info.write('enb_coord = ')
	info.write(str(enb_coord))
	info.write('\n')
	info.write('UEs = ')
	info.write(str(UEs))
	info.write('\n')
	info.write('load = ')
	info.write(str(load))
	info.write('\n')
	




        # Generating conflict graph
        # Note that <thr> is not used anymore and will be removed in future

	(Distance,RxPower,CQIVals) = GenerateGraphInfo(enb_coord,UEs)
	(i_map,i_map_ue,ue_interfering,edges,reuseClientsNum,Class1Clients) = GenerateGraph_CQI(Distance,RxPower,CQIVals,enb_coord,UEs,u_m)
	i_map_dir = GenerateDirGraph_CQI(Distance,RxPower,CQIVals,enb_coord)
        #(i_map,i_map_dir,i_map_ue,ue_interfering,edges,reuseClientsNum,Class1Clients) = gen_conflict_graph(enb_coord,UEs,u_m)
	info.write('Reuse Zone Clients\n')
	info.write(str(reuseClientsNum)+'\n')
	#print (reuseClients)
	(OpReuseSizes,IsolationLoad) = AdjustLoads(i_map,load,reuseClientsNum,N)
	#print (load)
	#print i_map_ue
	#print(OpReuseSizes)
	info.write('Operational ReuzeZoneSizes\n')
	info.write(str(OpReuseSizes)+'\n')


	#print(OpReuseSizes.values())
	#print(IsolationLoad)

        # Finds connected components withing i_map directional graph
        G = connected_graphs(i_map)
        G_dir = connected_graphs(i_map_dir)
	G_Ue = connected_graphs(i_map_ue)
        
	N = int(math.ceil(bw / rbgsize))


	writegeneralInfo(info,enb_coord,ue_list_to_print,UEs,load,i_map,i_map_dir,ue_interfering,edges,G,u_m)

        # Computes shares and assignments for every scheme
	#(Assign_FERMI_RZ,FERMIshare_RZ) = FermiAllocationsRZ(UEs,u_m,deepcopy(G),IsolationLoad,N,info,comp,OpReuseSizes)
	#print (FERMIshare_RZ)

	(Assign_FERMI,FERMIshare) = FermiAllocations(UEs,u_m,deepcopy(G),load,N,info,comp)
	Assign_FERMI_RZ = Assign_FERMI; # hack
	(Assign_FERMI_Dir,FERMIshare_Dir,FERMIshare_UE) = FermiAllocationsDir(UEs,u_m,deepcopy(G_Ue),load,N,info,comp,OpReuseSizes,Class1Clients)

	#(Assign_FERMI2,FERMIshare2) = FermiAllocations(UEs,u_m,deepcopy(G_dir),load,N,info,comp)

	(Assign_SEQ,Allocated_share_us,W) = SequentialAllocations(UEs,u_m,ue_interfering,edges,deepcopy(G),i_map_dir,load,N,info,comp)

	(Assign_SEQ_bo,Allocated_share_us) = SequentialAllocationsBackOff(UEs,u_m,ue_interfering,edges,deepcopy(G),i_map_dir,load,N,info,comp)

	(Assign_SEQ_FERM) = SequentialAllocationsFermiShare(UEs,u_m,ue_interfering,edges,deepcopy(G),i_map_dir,load,N,info,comp,FERMIshare)

	(IQshares,IQshares_1,IQshares_2,IQAggressiveShares,ctrlIntfrnc,IQedges) = IQhopping(Distance,RxPower,CQIVals,N,enb_coord,UEs,info)

	#(Assign_IQ_FERMI,IQ_FERMIshare) = FermiAssignIQAlloc(UEs,u_m,deepcopy(G_Ue),load,N,info,comp,OpReuseSizes,IQshares,FERMIshare_Dir,FERMIshare_UE,Class1Clients)
	Assign_IQ_FERMI = Assign_FERMI
		

        # Creates ns-3 scripts for different simulations
	#generate_topo_file(outputDir, bw,rbgsize,enb_coord,UEs, Assign_FERMI, Assign_FERMI_RZ, Assign_FERMI_Dir, Assign_IQ_FERMI, IQshares, IQshares_1, FERMIshare_Dir)
	generate_topo_file_udp_flows(outputDir, bw,rbgsize,enb_coord,UEs, Assign_FERMI, Assign_FERMI_RZ, Assign_FERMI_Dir, Assign_IQ_FERMI, IQshares, FERMIshare, IQshares_2,On,Off,IQedges,ctrlIntfrnc)
	generate_topo_file_wifi_flows(outputDir,enb_coord,UEs)
	generate_topo_file_wifiac_flows(outputDir,enb_coord,deepcopy(UEs))

	#generate_topo_file_udp_dynflows(outputDir, bw,rbgsize,enb_coord,UEs, Assign_FERMI, Assign_FERMI_RZ, Assign_FERMI_Dir, Assign_IQ_FERMI, IQshares, IQshares_1, IQshares_2,On,Off,edges,ctrlIntfrnc)
	#generate_topo_file_wifi_dynflows(outputDir,enb_coord,UEs)
	#FERMI 1
	#FERMI_RZ 2
	#FERMI_DIr 3
	#FERMI with IQ share 4
	#IQ DISTRIBUTED 5
	#IQ_1 6 
	#IQ_2 7

	#os.system('mv topo.cc generated-ns3-scripts/Fermi-'+str(thr)+'.cc')
	#os.system('cp generated-ns3-scripts/Fermi-'+str(thr)+'.cc /home/ghufran/workspace/ns-3-dev/scratch/Fermi-'+str(thr)+'.cc')


	info.close()
	comp.close()


	plot_graph('interferencemap', i_map, enb_coord, u_m, UEs)
	os.system('octave ' + outputDir + 'interferencemap.m')
	plot_ue_interference('UEinterferencemap', edges, enb_coord, u_m, UEs)
	os.system('octave ' + outputDir + 'UEinterferencemap.m')

	'''
	drawing = {}
	for i in Assign_mis:
		drawing[i]=[]
		for blk in Assign_mis[i]:
			drawing[i].append('['+str(blk[0])+','+str(blk[0])+','+str(blk[1]+1)+','+str(blk[1]+1)+'],[0,1,1,0]')

	for i in Assign_mis:
		pf.write('subplot('+str(len(UEs)/2)+',2,'+str(i+1)+');\n')
		pf.write('hold on\n')
		for blk in drawing[i]:
			pf.write('fill('+blk+',clr('+str(i+1)+',:));\n')

		k = 2
		for n in i_map[i]:
			if i not in i_map_dir[n]:
				for blk in drawing[n]:
					pf.write("plot("+blk+",'Color',clr("+str(n+1)+",:),'LineWidth',"+str(k)+");\n")
			k -= 0.3
		pf.write("axis('Position',[0 "+str(bw/rbgsize + 1)+" 0 1]);\n")
		pf.write("set(gca,'YTickLabel',[]);\n")
		pf.write('hold off;\n') 
	pf.write('print Blocks.jpg;\n')
	pf.close()
	'''

# Body, generating scripts
#os.system('mkdir ' + outputDir)
for z in range(7,8):
	l = 700
	w = 700
	os.system('mkdir -p res')
	os.system('mkdir -p generated-ns3-scripts')
	dens = float(8)/(l*w)
	for a in range(0,1):
	#for a in [1]:
		#os.system('echo '+str(a)+' >> comparison.txt')
		main(dens,l,w,3)
		#main(dens,l,w,a)

		os.system('mkdir res/'+str(a))
		os.system('mv res/info.txt res/'+str(a)+'/'+'info.txt')
		os.system('mv interferencemap.jpg res/'+str(a)+'/'+'interferencemap.jpg')
		os.system('mv UEinterferencemap.jpg res/'+str(a)+'/'+'UEinterferencemap.jpg')
		os.system('cp res/topo.cc generated-ns3-scripts/Topo_'+str(a)+'.cc')
		os.system('cp res/WiFitopo.cc generated-ns3-scripts/WiFi_'+str(a)+'.cc')
		os.system('cp res/WiFiactopo.cc generated-ns3-scripts/WiFiac_'+str(a)+'.cc')
		os.system('mv res/topo.cc res/'+str(a)+'/Topo.cc')
		os.system('mv res/WiFitopo.cc res/'+str(a)+'/WiFi.cc')

	#os.system('mv comparison.txt res/'+str(z)+'/comparison.txt')
