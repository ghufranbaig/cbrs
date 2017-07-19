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
from Fermi import Fermi_RZ
from Fermi import Fermi_UE
from Fermi import Fermi_Pre_Alloc
from Fermi import getCliques
from Fermi import FermiPreCompute
from helper_misc import *
import numpy as np

timesteps = 1

outputDir = "res/"

SNR = 3
rad = 100
bw = 100
rbgsize = 1
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
def prop_model(d,P_t):
	f = 3.600;
	d = d/1000.0;
	h_b = 30.0;
	h_m = 5.0;

	if (d == 0):
		d = 1

	A_f = 92.4+20*math.log10(d)+20*math.log10(f);
	A_bm = 20.41+9.83*math.log10(d)+7.894*math.log10(f)+9.56*(math.log10(f))**2;
	G_t = math.log10(h_b/200)*(13.958+5.78*(math.log10(d))**2);
	G_r = (42.57+13.7*math.log10(f))*(math.log10(h_m)-0.585);

	PL = A_f + A_bm - G_t - G_r;
	P_r = P_t/(10**(PL/10))

	# Taken from Rahman's model in ns-3	
	# n = 3.5
	# L_0 = 25.686
	# d_0 = 1
	#P_r = P_t/(10**((L_0 + 10*n*math.log10(d/d_0))/10))

	#print (P_t,P_r)
	return P_r

def rcvd_pows(bs_coords,ue):
	# Taken from Rahman's model in ns-3	
	P_t = 1
	pows = {}
	dists = {}
	for i in range(len(bs_coords)):
		bs = bs_coords[i]
		d = dist(bs,ue)
		P_r = prop_model(d,P_t)
		pows[i] = P_r
		dists[i] = d
	return (pows,dists)

def rcvd_pows_bs(bs_coords,enb):
	P_t = 1
	pows = {}
	dists = {}
	curr_bs = bs_coords[enb]
	for i in range(len(bs_coords)):
		bs = bs_coords[i]
		d = dist(bs,curr_bs)
		if d == 0:
			P_r = -1
		else:
			P_r = prop_model(d,P_t)
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

def gen_ue_coord(n,l,w,Op):
	coord = []
	ues = 0
	while(ues < n):
		x=random.randint(0,l)
		y=random.randint(0,w)
		for e in Op.eNBs:
			if (dist(e.coords,(x,y)) < 50):
				coord.append((x,y))
				ues += 1
				break
	return coord

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

	Distance_BS = {}
	RxPower_BS = {}

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


		(powers_bs,dist_bs) = rcvd_pows_bs(enb_coord,i)

		Distance_BS[i] = dist_bs
		RxPower_BS[i] = powers_bs

	
	return (Distance,RxPower,CQIVals,Distance_BS,RxPower_BS)

def GenerateGraph_BS(Distance_BS,RxPower_BS,enb_coord,UEs):
	n = len(enb_coord)
	G = {} # intereference graph based on CQI difference

	rx_power_thresh  = dBmtoW(-80)
	for i in range(n):
		G[i] = set()

	for i in range(n):
		
		#######################################
		# interference map generation based on reported CQI values
		# iterates over all eNBs
		for j in range(n):
			if (i==j):
				continue

			if (RxPower_BS[i][j] > rx_power_thresh):
				G[i].add(j)
				G[j].add(i)
								
	for i in G:
		G[i] = list(G[i])

	return (G)

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


def assign_ues_to_enbs(enbs,ues):
	assign_final = {}
	for u in ues:
		enb_id = 0
		min_d = dist(u,enbs[0])
		for i in range(len(enbs)):
			d=dist(u,enbs[i])
			if (d < min_d):
				enb_id = i
				min_d = d
		assign_final[u] = enb_id
	return assign_final


def remove_node(e,c_map):
	if (e in c_map):
		for n in c_map:
			if e in c_map[n]:
				c_map[n].remove(e)
		del c_map[e]

def FermiAllocations(UEs,u_m,G,load,N,info,comp,logging=False):
		assign = []
		alloc = []
		assign_grid = {}

		for enb in UEs:
			assign_grid[enb] = [0]*N


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

				for j in range(interval[0],interval[1]+1):
					assign_grid[a][j] = 1

			Allocated_share[a] = share
		
		share_UE={}
		for a in Allocated_share:
			if load[a] == 0:
				share = 0
			else:
				share = eval("%.2f" % (float(Allocated_share[a])/load[a]))
			for u in UEs[a]:
				share_UE[u_m[u]]=share

		writeInfo('FERMI',Assign,Allocated_share,share_UE,info)
		#print assign_grid
		if logging:
			print 'share: ', Allocated_share

		return (Assign,Allocated_share,assign_grid)

def FermiAllocationsSimple(UEs,FermiIntfMap,load,N):
		assign = []
		alloc = []
		assign_grid = {}

		for enb in UEs:
			assign_grid[enb] = [0]*N

                # Main calculation
		for c_map in FermiIntfMap:
			(channel_assignment,allocated_share) = FermiPreCompute(c_map[0],load,N,c_map[1],c_map[2],c_map[3])
			assign.append(channel_assignment)
			alloc.append(allocated_share)
			#assign_us.append(KlogN(c_map,W,N))

		Assign = {}
		for a in assign:
			Assign.update(a)



		Alloc_actual = {}
		for a in alloc:
			Alloc_actual.update(a)




		Allocated_share = {}
		for a in Assign:
			share = 0
			for interval in Assign[a]:
				share += interval[1]-interval[0]+1

				for j in range(interval[0],interval[1]+1):
					assign_grid[a][j] = 1

			Allocated_share[a] = share
			#if (load[a] != 0 and Allocated_share[a] == 0):
			#	print a, 'Something Wrong', Alloc_actual[a]
		#info2.write(str(load))
		#info2.write('\n')
		#writeInfo('FERMI',Assign,Allocated_share,Allocated_share,info2)

		return (Assign,Allocated_share,assign_grid)
 

def getSNR(N,assignment,Rx_power,enb,j):
	SNR = []
	for i in range(N):
		if assignment[enb][i] == 0:
			snr = -1
		else:
			I = 0;
			for bs in Rx_power:
				#print assignment[bs][i], bs, enb
				if (assignment[bs][i] == 1):
					I += Rx_power[enb][j][bs]
			snr = getSpectralEfficiency(Rx_power[enb][j][enb]/((Noise*spectrumBW[bw])+I-Rx_power[enb][j][enb])) #/ getSpectralEfficiency(Rx_power[enb][j][enb]/(Noise*spectrumBW[bw]))
			#print snr
		SNR.append(snr)
	return SNR


def process_activity(activity):
	processed = {}
	processed_ue = {}
	processed_bs = {}
	for ac in activity:
		for enb in ac:
			if enb not in processed:
				processed[enb] = {}
				processed_ue[enb] = {}
				processed_bs[enb] = 0
			for ue in ac[enb]:
				if ue not in processed[enb]:
					processed[enb][ue] = []
					processed_ue[enb][ue] = 0
				effi = 0
				for ef in ac[enb][ue]:
					if (ef != -1):
						effi += ef
				processed[enb][ue].append(effi)
				processed_ue[enb][ue] += effi
				processed_bs[enb] += effi




	#print processed_ue	
	return processed_ue	


class Operator:
	def __init__(self, op_id):
		self.ID = op_id
		self.UEs = []
		self.eNBs =[]

	def add_eNBs (self, eNB_coords, offset):
		i = 0
		for e in eNB_coords:
			self.eNBs.append(eNodeB(offset + i, e, self.ID))
			i += 1

	def add_UEs (self, UE_coords, offset):
		i = 0
		for u in UE_coords:
			self.UEs.append(UserE(offset + i, u, self.ID))
			#activity_for_ue = np.random.randint(2, size=(timesteps)).tolist()
			activity_for_ue = np.ones(timesteps).tolist()
			self.UEs[i].activity = activity_for_ue
			i += 1

	def assignUestoeNBs (self):
		ue_coords = []
		enb_coords = []
		for u in self.UEs:
			ue_coords.append(u.coords)
		for e in self.eNBs:
			enb_coords.append(e.coords)

		temp_assign = assign_ues_to_enbs(enb_coords,ue_coords)

		i = 0
		for u in ue_coords:
			assert (cmp(self.UEs[i].coords, u) == 0),"User coords don't match"

			self.UEs[i].eNB = self.eNBs[temp_assign[u]].ID
			self.eNBs[temp_assign[u]].UEs.append(self.UEs[i])
			i += 1

	def getData(self):
		operator_enbs = {}
		operator_ues = {}
		enb_coord = []
		UEs = {}
		UE_activity = {}
		#UE_activity = {}
		for e in self.eNBs:
			enb_coord.append(e.coords)
			UEs[e.ID] = []
			UE_activity[e.ID] = []

		operator_enbs[self.ID] = enb_coord
		operator_ues[self.ID] = []
		for u in self.UEs:
			operator_ues[self.ID].append(u.coords)
			UEs[u.eNB].append(u.coords)
			UE_activity[u.eNB].append(u.activity)
		return (operator_enbs,operator_ues,UEs,enb_coord,UE_activity)



class eNodeB:

    def __init__(self, enb_id, coord, operator):
    	self.ID = enb_id
        self.coords = coord
        self.operator = operator    # creates a new empty list for each dog
        self.UEs = []

    def add_ue(self, ue):
        self.UEs.append(ue)

class UserE:

    def __init__(self, ue_id, coord, operator):
    	self.ID = ue_id
        self.coords = coord
        self.operator = operator
        self.eNB = -1
        self.activity = []

	def assign_eNB(self, enb_id):
		self.eNB = enb_id


def run_ideal (UEs,u_m,G,N,UE_activity,info,comp,timesteps,Rx_power):
	activity = []
	for i in range(timesteps):
		load = {}

		for enb in UEs:
			active_ue = 0
			for j in range(len(UEs[enb])):
				if (UE_activity[enb][j][i] == 1):
					active_ue += 1
			load[enb] = active_ue
		#print load

		(Assign_FERMI,FERMIshare,assign_grid) = FermiAllocations(UEs,u_m,deepcopy(G),load,N,info,comp, True)
		activity.append(run_single (UEs,N,UE_activity,Rx_power,assign_grid,i))

	return process_activity(activity)

	#print activity

def run_baseline1 (UEs,u_m,G,N,UE_activity,info,comp,timesteps,Rx_power):
	activity = []
	for i in range(timesteps):
		load = {}

		for enb in UE_activity:
			load[enb] = 1
		#print load

		(Assign_FERMI,FERMIshare,assign_grid) = FermiAllocations(UEs,u_m,deepcopy(G),load,N,info,comp)
		activity.append(run_single (UEs,N,UE_activity,Rx_power,assign_grid,i))

	return process_activity(activity)
	#print activity

def run_baseline2 (UEs,u_m,G,N,UE_activity,info,comp,timesteps,Rx_power):
	activity = []
	assign_grid = {}
	for i in range(timesteps):
		for enb in UEs:
			assign_grid[enb] = [1]*N

		activity.append(run_single (UEs,N,UE_activity,Rx_power,assign_grid,i))

	#print activity[0]
	return process_activity(activity)

def run_creditBased1a (UEs,u_m,G,N,UE_activity,info,comp,timesteps,Rx_power,Operators):
	activity = []

	credit = {}
	for op in Operators:
		credit[op] = 100.0
	load = {}
	op_active_users = {}
	for i in range(timesteps):
		for op in Operators:
			for enb in op.eNBs:
				load[enb.ID] = credit[op]/len(op.eNBs)
		print load

		(Assign_FERMI,FERMIshare,assign_grid) = FermiAllocations(UEs,u_m,deepcopy(G),load,N,info,comp)
		activity.append(run_single (UEs,N,UE_activity,Rx_power,assign_grid,i))

	return process_activity(activity)

def run_creditBased1b (UEs,u_m,G,N,UE_activity,info,comp,timesteps,Rx_power,Operators):
	activity = []

	credit = {}
	for op in Operators:
		credit[op] = float(len(op.UEs))
	load = {}
	op_active_users = {}
	for i in range(timesteps):
		for op in Operators:
			for enb in op.eNBs:
				load[enb.ID] = credit[op]/len(op.eNBs)
		print load

		(Assign_FERMI,FERMIshare,assign_grid) = FermiAllocations(UEs,u_m,deepcopy(G),load,N,info,comp)
		activity.append(run_single (UEs,N,UE_activity,Rx_power,assign_grid,i))

	return process_activity(activity)

def run_creditBased2a (UEs,u_m,G,N,UE_activity,info,comp,timesteps,Rx_power,Operators):
	activity = []

	credit = {}
	for op in Operators:
		credit[op] = 100.0
	load = {}
	op_active_users = {}
	for i in range(timesteps):
		for op in Operators:
			tot_active = 0
			for enb in op.eNBs:
				active_ue = 0
				for j in range(len(enb.UEs)):
					if (enb.UEs[j].activity[i] == 1):
						active_ue += 1
				load[enb.ID] = active_ue
				tot_active += active_ue
			op_active_users[op] = tot_active
		for op in Operators:
			tot_active = 0
			for enb in op.eNBs:
				load[enb.ID] = credit[op]*load[enb.ID]/op_active_users[op]
		print load

		(Assign_FERMI,FERMIshare,assign_grid) = FermiAllocations(UEs,u_m,deepcopy(G),load,N,info,comp)
		activity.append(run_single (UEs,N,UE_activity,Rx_power,assign_grid,i))

	return process_activity(activity)

def run_creditBased2b (UEs,u_m,G,N,UE_activity,info,comp,timesteps,Rx_power,Operators):
	activity = []

	credit = {}
	for op in Operators:
		credit[op] = float(len(op.UEs))


	load = {}
	op_active_users = {}
	for i in range(timesteps):
		for op in Operators:
			tot_active = 0
			for enb in op.eNBs:
				active_ue = 0
				for j in range(len(enb.UEs)):
					if (enb.UEs[j].activity[i] == 1):
						active_ue += 1
				load[enb.ID] = active_ue
				tot_active += active_ue
			op_active_users[op] = tot_active
		for op in Operators:
			tot_active = 0
			for enb in op.eNBs:
				load[enb.ID] = credit[op]*load[enb.ID]/op_active_users[op]
		print load

		(Assign_FERMI,FERMIshare,assign_grid) = FermiAllocations(UEs,u_m,deepcopy(G),load,N,info,comp)
		activity.append(run_single (UEs,N,UE_activity,Rx_power,assign_grid,i))

	return process_activity(activity)

def getUtility(UE_activity,operators,ac):
	utils = {}
	zeroThroughput = {}
	for op in operators:
		utility = 0	
		zeroThroughput[op.ID] = 0
		for e in op.eNBs:
			enb = e.ID
			active_ue = 0
			for ue in range(len(e.UEs)):
				effi = 0
				for ef in ac[enb][ue]:
					if (ef != -1):
						effi += ef
				if (effi == 0 and UE_activity[enb][ue][0] == 1):
					utility += -100
					zeroThroughput[op.ID] += 1
				elif(effi != 0):
					utility += math.log(effi)
		utils[op.ID] = utility
	return (utils,zeroThroughput)


def run_single (UEs,N,UE_activity,Rx_power,assignment,step):
	curr_rates = {}
	for enb in UEs:
		curr_rates[enb] = {}
		for j in range(len(UEs[enb])):
			if (UE_activity[enb][j][step] == 0):
				curr_rates[enb][j] = []
			else:
				curr_rates[enb][j] = getSNR(N,assignment,Rx_power,enb,j)
	#print curr_rates
	return curr_rates

def getUtilityForOp(UEs,FermiIntfMap,load,N,UE_activity,Rx_power,op,t):
	(Assign_FERMI,FERMIshare,assign_grid) = FermiAllocationsSimple(UEs,FermiIntfMap,deepcopy(load),N)
	ac = run_single (UEs,N,UE_activity,Rx_power,assign_grid,t)
	(util,zeroThroughput) = getUtility(UE_activity,[op],ac)
	return util[op.ID]

def getAllUtil(UEs,FermiIntfMap,load,N,UE_activity,Rx_power,op,t):
	(Assign_FERMI,FERMIshare,assign_grid) = FermiAllocationsSimple(UEs,FermiIntfMap,deepcopy(load),N)
	ac = run_single (UEs,N,UE_activity,Rx_power,assign_grid,t)
	(util,zeroThroughput) = getUtility(UE_activity,op,ac)
	return (util,zeroThroughput)

def find_max_grad(UEs,FermiIntfMap,load,N,UE_activity,Rx_power,op, utility,t):
	#print 'Finding Max grad direction'
	#print 'Utility to beat:', utility
	#print load
	delta = 1
	E1 = 0
	E2 = 0
	util = utility
	optimalLoad = deepcopy(load)
	for e1 in op.eNBs:
		for e2 in op.eNBs:
			if (e1.ID == e2.ID):
				continue
			if (load[e1.ID] <= delta):
				continue
			if (load[e2.ID] == 0):
				continue

			newLoad = deepcopy(load)
			newLoad[e1.ID] -= delta
			newLoad[e2.ID] += delta
			#toprint = {x:newLoad[x] for x in range(0,10)}
			newutility = getUtilityForOp(UEs,FermiIntfMap,newLoad,N,UE_activity,Rx_power,op,t)
			#print e1.ID, '->', e2.ID, ':', newutility
			#print toprint
			#print ''

			if (newutility > util):
				E1 = e1.ID
				E2 = e2.ID
				util = newutility
				optimalLoad = deepcopy(newLoad)
				print 'new max', util
	return (E1,E2,optimalLoad,util)


def maximizeUtility(UEs,FermiIntfMap,load,N,UE_activity,Rx_power,op,t):
	utility = getUtilityForOp(UEs,FermiIntfMap,load,N,UE_activity,Rx_power,op,t)
	print utility
	delta = 1
	changed_utility = 0
	
	it = 0

	while (it < 2):
		(E1,E2,optimalLoad,utility) = find_max_grad(UEs,FermiIntfMap,deepcopy(load),N,UE_activity,Rx_power,op, utility,t)
		if (E1 == E2):
			break
		changed_utility = 1
		load = optimalLoad
		#print 'Iterating on the max gradient', E1 ,'->', E2
		while (load[E1] > delta):
			nload = deepcopy(load)
			nload[E1] -= delta
			nload[E2] += delta
			newutility = getUtilityForOp(UEs,FermiIntfMap,nload,N,UE_activity,Rx_power,op,t)
			if(newutility <= utility):
				break
			utility = newutility
			print utility
			load = nload
		print '\n\n'
		it += 1
	
	print 'Final Utility = ', utility
	return (load,changed_utility)

def run_creditBased2bWith (UEs,u_m,G,N,UE_activity,info,comp,timesteps,Rx_power,Operators):
	activity = []
	util_thresh = 1e-5

	All_util = []
	credit = {}
	for op in Operators:
		credit[op] = float(len(op.UEs))

	FermiIntfMap = []
	for m in G:
		FermiIntfMap.append(getCliques(m))


	load = {}
	op_active_users = {}
	for i in range(timesteps):
		for op in Operators:
			tot_active = 0
			for enb in op.eNBs:
				active_ue = 0
				for j in range(len(enb.UEs)):
					if (enb.UEs[j].activity[i] == 1):
						active_ue += 1
				load[enb.ID] = active_ue
				tot_active += active_ue
			op_active_users[op] = tot_active
		for op in Operators:
			tot_active = 0
			for enb in op.eNBs:
				load[enb.ID] = int(credit[op]*load[enb.ID]/op_active_users[op])
		#print load			
		(origUtil,zeroUe) = getAllUtil(UEs,FermiIntfMap,load,N,UE_activity,Rx_power,Operators,i)
		#info2.write('Orig: \n')
		#info2.write(str(load)+'\n')
		info2.write(str(origUtil)+'\n')
		#info2.write(str(zeroUe)+'\n')
		All_util.append(origUtil)
		finalLoad = {}
		prevUtil = origUtil
		for it in range(100):
			#anyChanged = 0
			print str(it)+'th iteration'
			for op in Operators:
				print op.ID, 'Optimizing'
				(nload,changed) = maximizeUtility(UEs,FermiIntfMap,load,N,UE_activity,Rx_power,op,i)
				(newUtil,zeroUe) = getAllUtil(UEs,FermiIntfMap,nload,N,UE_activity,Rx_power,Operators,i)
				for e in op.eNBs:
					finalLoad[e.ID] = nload[e.ID]
				load = nload
				All_util.append(newUtil)
				#anyChanged += changed
			#(Assign_FERMI,FERMIshare,assign_grid) = FermiAllocations(UEs,u_m,deepcopy(G),load,N,info,comp)
			#run_single (UEs,N,UE_activity,Rx_power,assign_grid,i)

				#print 'Orig: ', origUtil
				#print 'New:',op.ID, newUtil
				#info2.write(str(op.ID) + ' Optimizes \n')
				#info2.write(str(nload) + '\n')
				info2.write(str(newUtil) + '\n')
				#info2.write(str(zeroUe) + '\n')
			#(finalUtil,finalzeroUe) = getAllUtil(UEs,FermiIntfMap,finalLoad,N,UE_activity,Rx_power,Operators,i)

			anyChanged = 0
			for o in prevUtil:
				if (abs(prevUtil[o] - newUtil[o]) > util_thresh):
					anyChanged = 1
			if (anyChanged == 0):
				break
			prevUtil = newUtil
			#info2.write('Everyone Optimizes \n')
			#info2.write(str(finalLoad) + '\n')
			#info2.write(str(finalUtil)+'\n')
			#info2.write(str(finalzeroUe)+'\n')
			#info2.write('\n\n')
	return All_util

def main(density,l,w,toytop):

        # number of sub-channels
	N = 100 #int(math.ceil(bw / rbgsize))
	print N
	operators = 3
	operator_enbs = {}
	operator_ues = {}
	enb_coord = []
	UEs = {}
	UE_activity = {}
	Activity = {}

	Operators = []

	usersPerOperator = {0:50,1:50,2:50}



    # size of the grid (in meters)
	#l = 500
	#w = 500

        # number of eNodeBs
	#n = int(math.floor((l*w)*density + 0.5))
	npo = [10,10,10]
	print(npo)
	j = 0
	k = 0;
	for i in range(operators):
		Op = Operator(i)
		temp_coords = gen_eNbs_coord(npo[i],l,w)
		Op.add_eNBs(temp_coords,j)
		j += len(temp_coords)

		temp_coords = gen_ue_coord(usersPerOperator[i],l,w,Op)
		Op.add_UEs(temp_coords,k)
		k += len(temp_coords)
		Op.assignUestoeNBs()
		Operators.append(Op)


	for op in Operators:
		(a,b,c,d,e) = op.getData()
		operator_enbs.update(a)
		operator_ues.update(b)
		UEs.update(c)
		enb_coord += d
		UE_activity.update(e)



	load = {}
	for i in range(len(enb_coord)):
		load[i] = len(UEs[i])

	print load

	
	ue_list_to_print={}

	'''
	enb_coord = [(100,0),(150,0)]
	UEs={0: [(120,0)], 1:[(110,0)]}
	load={0: 1, 1: 1 }
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

	(Distance,RxPower,CQIVals,Distance_BS,RxPower_BS) = GenerateGraphInfo(enb_coord,UEs)
	i_map = GenerateGraph_BS(Distance_BS,RxPower_BS,enb_coord,UEs)

        # Finds connected components withing i_map directional graph
	G = connected_graphs(i_map)
	#N = int(math.ceil(bw / rbgsize))


	writegeneralInfo(info,enb_coord,ue_list_to_print,UEs,load,i_map,i_map,{},{},G,u_m)

	results = {}
	#(Assign_FERMI,FERMIshare,assign_grid) = FermiAllocations(UEs,u_m,deepcopy(G),load,N,info,comp)
	#run_single (UEs,N,UE_activity,RxPower,assign_grid)
	# res = run_ideal (UEs,u_m,G,N,UE_activity,info,comp,timesteps,RxPower)
	# i = 0
	# for e in res:
	# 	for u in res[e]:
	# 		if i not in results:
	# 			results[i] = []
	# 		results[i].append(res[e][u])
	# 		i += 1
	# res = run_baseline1 (UEs,u_m,G,N,UE_activity,info,comp,timesteps,RxPower)
	# i = 0
	# for e in res:
	# 	for u in res[e]:
	# 		if i not in results:
	# 			results[i] = []
	# 		results[i].append(res[e][u])
	# 		i += 1
	# res = run_baseline2 (UEs,u_m,G,N,UE_activity,info,comp,timesteps,RxPower)
	# i = 0
	# for e in res:
	# 	for u in res[e]:
	# 		if i not in results:
	# 			results[i] = []
	# 		results[i].append(res[e][u])
	# 		i += 1
	# res = run_creditBased1a (UEs,u_m,G,N,UE_activity,info,comp,timesteps,RxPower,Operators)
	# i = 0
	# for e in res:
	# 	for u in res[e]:
	# 		if i not in results:
	# 			results[i] = []
	# 		results[i].append(res[e][u])
	# 		i += 1
	# res = run_creditBased1b (UEs,u_m,G,N,UE_activity,info,comp,timesteps,RxPower,Operators)
	# i = 0
	# for e in res:
	# 	for u in res[e]:
	# 		if i not in results:
	# 			results[i] = []
	# 		results[i].append(res[e][u])
	# 		i += 1

	# res = run_creditBased2a (UEs,u_m,G,N,UE_activity,info,comp,timesteps,RxPower,Operators)
	# i = 0
	# for e in res:
	# 	for u in res[e]:
	# 		if i not in results:
	# 			results[i] = []
	# 		results[i].append(res[e][u])
	# 		i += 1
	# res = run_creditBased2b (UEs,u_m,G,N,UE_activity,info,comp,timesteps,RxPower,Operators)
	# i = 0
	# for e in res:
	# 	for u in res[e]:
	# 		if i not in results:
	# 			results[i] = []
	# 		results[i].append(res[e][u])
	# 		i += 1

	#res = run_ideal (UEs,u_m,G,N,UE_activity,info,comp,timesteps,RxPower)

	
	All_util = run_creditBased2bWith (UEs,u_m,G,N,UE_activity,info,comp,timesteps,RxPower,Operators)
	f = open('results.dat','w')
	#print All_util
	for i in range(len(All_util)):
		for j in range(len(Operators)):
			f.write(str(All_util[i][j]))
			if (j == len(Operators)-1):
				f.write('\n')
			else:
				f.write(',')
	f.close()
	os.system('octave plotUtils.m')
	








	info.close()
	comp.close()


	plot_graph(outputDir,'interferencemap', i_map, enb_coord, u_m, UEs,l,w,npo[0])
	os.system('octave ' + outputDir + 'interferencemap.m')
	#plot_ue_interference(outputDir,'UEinterferencemap', edges, enb_coord, u_m, UEs,l, w)
	#os.system('octave ' + outputDir + 'UEinterferencemap.m')

info2 = open('convergence.txt','w')


# Body, generating scripts
#os.system('mkdir ' + outputDir)
for z in range(10):
	l = 300
	w = 300
	info2.write(str(z)+'\n')
	main(1,l,w,3)
	os.system('mv res/utils.jpg res/utils_'+str(z)+'.jpg')
	os.system('mv res/interferencemap.jpg res/interferencemap_'+str(z)+'.jpg')
	os.system('mv convergence.txt res/convergence_'+str(z)+'.txt')

	#for a in [1]:
		#os.system('echo '+str(a)+' >> comparison.txt')
		#main(dens,l,w,3)
		#main(dens,l,w,a)

		#os.system('mkdir res/'+str(a))
		#os.system('mv res/info.txt res/'+str(a)+'/'+'info.txt')
		#os.system('mv interferencemap.jpg res/'+str(a)+'/'+'interferencemap.jpg')
		#os.system('mv UEinterferencemap.jpg res/'+str(a)+'/'+'UEinterferencemap.jpg')
		#os.system('cp res/topo.cc generated-ns3-scripts/Topo_'+str(a)+'.cc')
		#os.system('cp res/WiFitopo.cc generated-ns3-scripts/WiFi_'+str(a)+'.cc')
		#os.system('cp res/WiFiactopo.cc generated-ns3-scripts/WiFiac_'+str(a)+'.cc')
		#os.system('mv res/topo.cc res/'+str(a)+'/Topo.cc')
		#os.system('mv res/WiFitopo.cc res/'+str(a)+'/WiFi.cc')

	#os.system('mv comparison.txt res/'+str(z)+'/comparison.txt')
info2.close()
