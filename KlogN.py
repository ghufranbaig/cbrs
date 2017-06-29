import random
import math
from copy import deepcopy

steps = 0
def lubys (G):
	global steps
	steps += 1
	N = len(G)
	if N == 0 :
		return set()
	else :
		prio_ = [i for i in range(N)]
		random.shuffle(prio_)
		V = list(G.keys())
		W = list(G.keys())

		prio = {}
		for i in range(N):
			prio[V[i]] = prio_[i]
		for i in V:
			for j in G[i]:
				if prio[i] < prio[j] :
					W.remove(i)
					break
		G_dash = deepcopy(G)
		minus = set()
		for v in W:
			minus.add(v)
			for n in G[v]:
				minus.add(n)

		for m in minus:
			del G_dash[m]
			for v in G_dash :
				if m in G_dash[v]:
					G_dash[v].remove(m)
		return set(W)|lubys(G_dash)
					
		
def expand_graph(G,W):
	N = 0
	mapping = {}
	# Generate mapping from 1 node to w nodes
	for v in G:
		mapping[v] = [i for i in range(N,N+W[v])]
		N += W[v]
	G_dash = {}
	
	#for every node
	for v in mapping:
		#for every newly generated node
		for i in mapping[v]:
			#each new node in the same clique must be its neighbor
			G_dash[i] = deepcopy(mapping[v])
			G_dash[i].remove(i)
			for n in G[v]:
				G_dash[i] += mapping[n]
	return (G_dash,mapping)


def expand_graph2(G,W,UEs,edges):
	#N = 0
	#mapping = {}
	# Generate mapping from 1 node to w nodes
	#for v in G:
	#	mapping[v] = [i for i in range(N,N+W[v])]
	#	N += W[v]

	#print(mapping)
	#print(W)
	N = 0
	mapping = {}
	ue_to_v = {}
	for v in UEs:
		mapping[v] = []
		j=0
		n=len(UEs[v]) 
		x=0
		for u in UEs[v]:
			if (j+1==n):
				w=W[v]-x
			elif(j%n==0):
				w=int(math.ceil(float(W[v])/n))
			else:
				w=int(math.floor(float(W[v])/n))

			x+=w
			j+=1
			m = [i for i in range(N,N+w)]
			N += w
			mapping[v] += m
			ue_to_v[u] = m
	#print(mapping)
			
		
	G_dash = {}
	
	#for every node
	for v in mapping:
		#for every newly generated node
		for i in mapping[v]:
			#each new node in the same clique must be its neighbor
			G_dash[i] = deepcopy(mapping[v])
			G_dash[i].remove(i)

	#add edges for interference

	for v in edges:
		for e in edges[v]:
			#for i in mapping[v]:
				#G_dash[i] += ue_to_v[e]
			for i in ue_to_v[e]:
				G_dash[i] += mapping[v]

	#for e in edges[4]:
	#	print(ue_to_v[e])
	#print(3,ue_to_v[UEs[3][0]])
	#for e in UEs[4]:
	#	print(4,ue_to_v[e])

	#print(G_dash)

	return (G_dash,mapping,ue_to_v)			
		


def KlogN(G,UEs,W,edges,N):


	#(e_G,mapping) = expand_graph(G,W)
	(e_G,mapping,ue_to_v) = expand_graph2(G,W,UEs,edges)
	Assign = {}
	for i in range(N):
		S = lubys(e_G)
		'''
		remove_from_s = set()
		for s in S:
			for n in e_G[s]:
				if n in S:
					remove_from_s.add(s)
		for s in remove_from_s:
			S.remove(s)
		'''
		for v in S:
			Assign[v] = [i]
			# minus v from Graph
			del e_G[v]
			for v_ in e_G :
				if v in e_G[v_]:
					e_G[v_].remove(v)


	subchannels = {}
	subchannels_ue = {}
	#set unallocated nodes to empty
	for v in e_G:
		Assign[v] = []
	#combine the allocations for each vertex
	for v in G:
		subchannels[v] = []
		for n in mapping[v]:
			subchannels[v]+=Assign[n]
	for u in ue_to_v:
		subchannels_ue[u] = []
		for n in ue_to_v[u]:
			subchannels_ue[u]+=Assign[n]
	return (subchannels,subchannels_ue)
			
def KlogN_Backoff(G,UEs,W,edges,N):


	#(e_G,mapping) = expand_graph(G,W)
	(e_G,mapping,ue_to_v) = expand_graph2(G,W,UEs,edges)
	Assign = {}
	for i in range(N):
		global steps
		steps = 0
		S = lubys(e_G)
		#print (steps)

		remove_from_s = set()
		for s in S:
			for n in e_G[s]:
				if n in S:
					remove_from_s.add(s)
		for s in remove_from_s:
			S.remove(s)

		for v in S:
			Assign[v] = [i]
			# minus v from Graph
			del e_G[v]
			for v_ in e_G :
				if v in e_G[v_]:
					e_G[v_].remove(v)


	subchannels = {}
	subchannels_ue = {}
	#set unallocated nodes to empty
	for v in e_G:
		Assign[v] = []
	#combine the allocations for each vertex
	for v in G:
		subchannels[v] = []
		for n in mapping[v]:
			subchannels[v]+=Assign[n]
	for u in ue_to_v:
		subchannels_ue[u] = []
		for n in ue_to_v[u]:
			subchannels_ue[u]+=Assign[n]
	return (subchannels,subchannels_ue)
			
'''	
	
i_map = {0: [1,3],
	1: [0,2],
	2: [1,3,4,5],
	3: [0,2,6],
	4: [2,5],
	5: [2,4],
	6: [3]
	}
load = {0:1/3,1:1/4,2:2/8,3:1/7,4:2/5,5:1/5,6:3/4}
for i in load:
	load[i] = math.floor(20*load[i]+0.5)
KlogN(i_map,load,20)
#print(lubys(i_map))
'''
