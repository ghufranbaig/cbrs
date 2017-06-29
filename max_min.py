import random
import math
import os
from copy import deepcopy

def get_btlneck(v):
	n = 0
	min_ = None
	for i in v:
		if n==0 :
			min_ = i
			n = 1
		else :
			if (v[i]<v[min_]):
				min_ = i
	return min_

def max_min_share(UEs,interfering_ue,N):
	n=0
	share = {}
	ue_list = {}
	constraints = {}
	R = {}
	for enb in UEs:
		R[enb] = N
		c = []
		for u in UEs[enb]:
			share[n]=-1
			ue_list[u] = n
			c.append(n)
			n+=1
		constraints[enb] = c
	v = {}
	for enb in interfering_ue:
		for u in interfering_ue[enb]:
			constraints[enb].append(ue_list[u])
		v[enb]=float(R[enb])/len(constraints[enb])

	num = len(constraints)

	for i in range(num):
		if (len(v)==0):
			break

		b=get_btlneck(v)
		share_ = v[b]
		mark_to_del = []
		for u in constraints[b]:
			share[u] = share_
			mark_to_del.append(u)
			for c in constraints:
				if u in constraints[c]:
					#constraints[c].remove(u)
					R[c] -= share_
		for u in mark_to_del:
			for c in constraints:
				if u in constraints[c]:
					constraints[c].remove(u) 
		del v[b]
		mark_to_del = []
		for x in v:
			if (len(constraints[x])>0):
				v[x]=float(R[x])/len(constraints[x])
			else :
				mark_to_del.append(x)
		for n in mark_to_del:
			del v[n]

	#for i in range(len(share)):
	#	share[i] = int(math.floor(share[i]+0.5)) 
	
	enb_share = {}
	for enb in UEs:
		enb_share[enb]=0
		for u in UEs[enb]:
			enb_share[enb] += share[ue_list[u]]
		enb_share[enb] = int(math.floor(enb_share[enb]+0.5))

	for i in share:
		share[i] = eval("%.2f" % share[i])
	return (share,enb_share)



def get_shares(W,c,R):
	shares={}
	L=0
	for n in c:
		L+=W[n]
	for n in c:
		shares[n]=R*float(W[n])/L
	return shares


def get_min1(shares):
	b,s=-1,-1
	for c in shares:
		for i in shares[c]:
			if b==-1:
				b=i
				s=shares[c][i]
			elif shares[c][i] < s :
				b=i
				s=shares[c][i]
	if (b==-1):
		print('error in get min for WMMF')
	return(b,s)

def get_min(shares,b):
	s=-1
	for c in shares:
		if (b in shares[c]):
			if s==-1:
				s=shares[c][b]
			elif shares[c][b] < s :
				s=shares[c][b]
	if (s==-1):
		print(b,shares,'error in get min for WMMF')
	return(s)
def better_alloc(share,prev):

	a = sorted(share.values())
	b = sorted(prev.values())
		
	return (a>b)
	
	sum1,sum2 = 0,0
	for n in share:
		sum1 += share[n]
		sum2 += prev[n] 
	
	if (sum2>sum1):
		return False
	else:
		return True
	for n in share:
		if (share[n]>prev[n]):
			for k in prev:
				if (n==k):
					continue
				if ((prev[k]<=prev[n]) and (share[k]<prev[k])):
					return False
	

	return True

import itertools
def WMMF (i_map,W,N):

	Resources={}
	Constraints = {}
	v={}
	share = {}
	for enb in i_map:
		Resources[enb] = N
		Constraints[enb] = [enb]
		w=W[enb]
		for x in i_map[enb]:
			Constraints[enb].append(x)
			w+=W[x]
		v[enb]=Resources[enb]/w

	k = 0
	Resources={}
	Constraints = {}
	share={}
	for enb in i_map:
		Resources[k] = N
		Constraints[k] = [enb]
		k+=1
		for x in i_map[enb]:
			Resources[k] = N
			Constraints[k] = [enb, x]
			k+=1
	
	constraints = deepcopy(Constraints)
	R = deepcopy(Resources)
	Shares = {}
	for c in constraints:
		Shares[c] = get_shares(W,constraints[c],R[c])
	num = len(i_map)

	
	for i in range(num):
		(b,s)=get_min1(Shares)
		share[b]=s

		for c in constraints:
			if b in constraints[c]:
				constraints[c].remove(b)
				R[c]-=s

		Shares = {}
		for c in constraints:
			Shares[c] = get_shares(W,constraints[c],R[c])
	for s in share:
		share[s]=int(math.floor(share[s]+0.5))
	return(share)
	
def WMMF2 (i_map,W,N):

	Resources={}
	Constraints = {}
	v={}
	share = {}
	for enb in i_map:
		Resources[enb] = N
		Constraints[enb] = [enb]
		w=W[enb]
		for x in i_map[enb]:
			Constraints[enb].append(x)
			w+=W[x]
		v[enb]=Resources[enb]/w

	nodes = list(Constraints.keys())
	all_perm = list(itertools.permutations(nodes))
	#print (all_perm)
	final_share = {}

	#all_perm = [[1,4,2,3,6,5,7,0],[5,3,7,6,2,1,0,4]]
	for perm in all_perm:
		share={}
		constraints = deepcopy(Constraints)
		R = deepcopy(Resources)
		for i in perm:		
			Shares = {}
			for c in constraints:
				Shares[c] = get_shares(W,constraints[c],R[c])

			s=get_min(Shares,i)
			share[i]=s
			#print constraints
			#print Shares
			#print i,s 
			#print ('\n')
			for c in constraints:
				if i in constraints[c]:
					constraints[c].remove(i)
					R[c]-=s
		if (len(final_share)==0):
			final_share = share
			print share
		elif (better_alloc(share,final_share)):
			final_share = share
		#print(perm,share)

	for s in final_share:
		final_share[s]=int(math.floor(final_share[s]+0.5))
	return (final_share)

'''
i_map = {0: [1], 1: [0,2], 2: [1]}
W = {0: 3, 1: 4, 2: 7}
N = 11
print(WMMF(i_map,W,N))
i_map = {0: [1], 1: [0, 7], 2: [3, 5, 6], 3: [2, 4, 5], 4: [3, 5], 5: [2, 3, 4, 6], 6: [2, 5, 7], 7: [1, 6]}
W = {0: 2, 1: 3, 2: 2, 3: 1, 4: 2, 5: 3, 6: 2, 7: 1}
N = 25
'''



