import math
import networkx as nx
import Queue as Q
from copy import deepcopy

#Tree class
class Node(object):
    def __init__(self, data):
        self.data = data
        self.children = []

    def add_child(self, obj):
        self.children.append(obj)
    def print_tree(self):
	print(self.data)
	print('--------')
	for c in self.children:
		c.print_tree()


#Triangulation
def BFS(graph,start,end,q):
	#print('BFS')
        paths = []
        temp_path = [start]	
        q.put(temp_path)
        while q.empty() == False:
                tmp_path = q.get()
                last_node = tmp_path[len(tmp_path)-1]
                if last_node == end:
                    paths.append(tmp_path)
                    #print ("VALID_PATH : ",tmp_path)
                for link_node in graph[last_node]:
                    if link_node not in tmp_path:
                        new_path = []
                        new_path = tmp_path + [link_node]
                        q.put(new_path)
        return paths


def max_u(u,w):
    m = u[0]
    for i in range(1,len(u)):
        if (w[u[i]]>=w[m]):
            m = u[i]
    return m


# The algorithm to transform the inital graph into chordal graph
# as mentioned in the paper in Algorithm 1, step 1 (and describe in another paper)
# 
# Input: 
#   G : interference graph - map with keys = nodes and value = [adjacent nodes]
# Output:
#   G : chordal interf graph in the same format
#   fill_in : edges added to make the graph chordal, to be removed later

def remove_node(e,c_map):
	if (e in c_map):
		for n in c_map:
			if e in c_map[n]:
				c_map[n].remove(e)
		del c_map[e]
	return c_map

def dfs(graph, start):
    if len(graph) == 0:
       return set()
    visited, stack = set(), [start]
    #print '**********************************'
    #print graph
    #print start
    while stack:
        vertex = stack.pop()
        if vertex not in visited:
            visited.add(vertex)
            #print vertex
            stack.extend(graph[vertex] - visited)
    #print visited
    #print '**********************************'
    return visited	

def Triangulation (G):

    F = set()
    w = {}
    alpha = {}
    un_numbered = []
    for v in G:
        w[v] = 0
        un_numbered.append(v)
    n = len(G)
    #print n

    G_unnumbered = deepcopy(G)
    for node in G_unnumbered:
        G_unnumbered[node] = set(G_unnumbered[node])

    for i in range (n,0,-1):
        #print (un_numbered)
        z = max_u(un_numbered,w)
        alpha[z] = i
        un_numbered.remove(z)
        wz_minus = deepcopy(w)
        for y in un_numbered:
           if z in G[y]:
		w[y] += 1
                #F.add((y,z))
                continue
           G_dash = deepcopy(G_unnumbered)
           for node in G_unnumbered:
               if (node != z and node != y and wz_minus[node] >= wz_minus[y]):
			G_dash = remove_node(node,G_dash)
	   
           p = dfs(G_dash,z)
           #print p
           if y in p :
                w[y] += 1
                F.add((y,z))
        G_unnumbered = remove_node(z,G_unnumbered)
        #print(z,w)


    for a in G:
        G[a]=set(G[a])

    fill_in = []
    for edge in F:
	if edge[1] not in G[edge[0]]:
        	G[edge[0]].add(edge[1])
        	G[edge[1]].add(edge[0])
		fill_in.append(edge)
    for a in G:
        G[a]=list(G[a])
    return (G,fill_in)


def Triangulation2 (G):
    path_queue = Q.Queue();
    F = set()
    w = {}
    alpha = {}
    un_numbered = []
    for v in G:
        w[v] = 0
        un_numbered.append(v)
    n = len(G)
    for i in range (n,0,-1):
        #print (w)
        z = max_u(un_numbered,w)
        #print(z)
        alpha[z] = i
        un_numbered.remove(z)
        wz_minus = deepcopy(w)
        for y in un_numbered:
           p = BFS(G,y,z,path_queue)

           
           inc = False
           for x in p:
               flag = True  
               invalid_path = False             
               if (len(x) > 2):
                    for j in range (1,len(x)-1):
                       if (x[j] not in un_numbered):
                           invalid_path = True
                           flag = False
                           break
                       if (wz_minus[x[j]] >= wz_minus [y]):
                           flag = False
                           break
               if invalid_path :
                       continue
               if flag :
                       inc = True
                       break
           if inc :
                w[y] += 1
                F.add((y,z))
        #print(z,w)


    for a in G:
        G[a]=set(G[a])

    fill_in = []
    for edge in F:
	if edge[1] not in G[edge[0]]:
        	G[edge[0]].add(edge[1])
        	G[edge[1]].add(edge[0])
		fill_in.append(edge)
    for a in G:
        G[a]=list(G[a])
    return (G,fill_in)



# Allocation
def s_tuple(v_i,C,L):
        t_i = 0
        l_i = -1
        for j in range (len(C)):
                if v_i in C[j]:
                        t_i += 1
                        if L[j] > l_i:
                                l_i = L[j]
        return (l_i,t_i)

def init_alloc (v_i,C,l,R):

        assigned = False
        for j in range (len(C)):
                l_sum = 0
                if v_i in C[j]:
                        for k in C[j]:
                                l_sum += l[k]

			if (l_sum > 0): # this check added to incorporate nodes with share = 0
				new_val = math.floor(l[v_i]*R[j]/l_sum + 0.5)
			else:
                        	new_val = 0
                        if  not assigned:
                                A = new_val
                                assigned = True
                        elif new_val < A:
                                A = new_val
        return A

def max_rank(s):
        max_r = 0
        for i in range (1,len(s)):
                if (s[i][0]>s[max_r][0]):
                        max_r = i
                elif (s[i][0]==s[max_r][0]):
                        if (s[i][1]>s[max_r][1]):
                                max_r = i
        return max_r
        
      
def cliques(graph):
	G=nx.Graph()
	for n in graph:
		G.add_node(n)
	for n in graph:
		for e in graph[n]:
			G.add_edge(n,e)


	C = list(nx.chordal_graph_cliques(G))
	for i in range (len(C)):
		C[i] = list(C[i])
	return C



# Calculate Fermi allocation as described in the paper as Algorithm 2
# and also the first part of step 2 of Algorithm 1 in the paper
# Input:
#   G - chordal graph -
#   load - map (node, #users)
#   N - number of sub-channels
#   C - set of maximal cliques, array of cliques, each being an array of belonging nodes
# Output:
#   Alloc - share per eNodeB, map (eNodeB_ID, share)
def Allocate (G,load,N,C):
        V = list(G.keys())
        U = V[:]
        A = {}
        L = []
        R = []
        for j in range(len(C)):
                sum_load = 0
                for i in C[j]:
                        sum_load += load[i]
                L.append(sum_load)
                R.append(N)
        s =[]
        for i in U:
                s.append(s_tuple(i,C,L))
                
        for i in U:
                A[i]=init_alloc (i,C,load,R)
                
        Alloc = {}
        while len(U) > 0:
                v_0 = U[max_rank(s)]
		#print (A)
                Alloc[v_0]=int(A[v_0])
                U.remove(v_0)
                for j in range (len(C)):
                        if v_0 in C[j]:
                                R[j] -= A[v_0]
                                C[j].remove(v_0)
                                L[j] -= load[v_0]
                s =[]
                for i in U:
                        s.append(s_tuple(i,C,L))
                
                for i in U:
                        A[i]=init_alloc (i,C,load,R)
                
        return Alloc


# Assignment & Restoration

def find_unalloc_block(sub_channels,i):
	for n in range (i,len(sub_channels)):
		if sub_channels[n] == 1:
			return (i,n-1)
	return (i,len(sub_channels)-1)

def find_unalloc_block2(sub_channels,i,j):
	for n in range (i,j+1):
		if sub_channels[n] == 1:
			return (i,n-1)
	return (i,j)

def find_unalloc(sub_channels,i):
	for n in range (i,len(sub_channels)):
		if sub_channels[n] == 0:
			return n
	return -1


def get_max_blk(blks_avail):
	blk = blks_avail[0]
	for block in blks_avail:
		if (block[1]-block[0]) > (blk[1]-blk[0]) :
			blk = block
	return blk
	

def allocate_sub_chan(v,sub_channels,alloc):

	start = find_unalloc(sub_channels,0)
	blks_avail = []
	while (start != -1) :

		blk = find_unalloc_block(sub_channels,start)

		if (blk[1]-blk[0]+1) >= alloc :
			return [(blk[0],blk[0]+alloc-1)]
		blks_avail.append(blk)
		start = find_unalloc(sub_channels,blk[1]+1)
	rem = alloc
	assigned = []
	while not (rem == 0 or len(blks_avail)==0):
		blk = get_max_blk(blks_avail)
		assigned.append(blk)
		blks_avail.remove(blk)
		rem -= blk[1]-blk[0]+1
	return assigned
		

def color_channel(assigned,subchannels):
	for block in assigned:
		for i in range(block[0],block[1]+1):
			subchannels[i] = 1
	return subchannels

def makeTree(root,C):
	children = []
	for c in C:
		if len(set(c) & set(root.data)) > 0:
			children.append(c)
	for ch in children:
		C.remove(ch)
	for ch in children:
		root.add_child(makeTree(Node(ch),C))
	return root



# Second part of step 2 of Algorithm 1 in Fermi, summarized in the text
# Assignes the actual channels based on previously calculated shares
# Input:
#   Allocation - from the previous step
#   N - number of channels
#   C - cliques
# Output:
#   Assign - map (nodeID, [channels])
def Assignment(Alloc,N,C):
	
	root = Node(C[-1])
	C.remove(C[-1])
	tree = makeTree(root,C)	
	
	Assign = {}
	subchannels = [0 for i in range(N)]
	q = Q.Queue()
	q.put(tree)
	U = list(Alloc.keys())
	while (q.empty() == False):
		curr_node = q.get()
		subchannels = [0 for i in range(N)]
		for ch in curr_node.children:
			q.put(ch)
		for v in curr_node.data:
			if v in Assign:
				subchannels = color_channel(Assign[v],subchannels)
		for v in curr_node.data:
			if v in U:
				Assign[v]=allocate_sub_chan(v,subchannels,Alloc[v]) 
				U.remove(v)
				subchannels = color_channel(Assign[v],subchannels)
	return Assign


# Second part of step 2 of Algorithm 1 in Fermi, summarized in the text
# Assignes the actual channels based on previously calculated shares
# Input:
#   Allocation - from the previous step
#   N - number of channels
#   C - cliques
#   ReuzeZoneSizes - operational ReuseZoneSizes for each eNB
# Output:
#   Assign - map (nodeID, [channels])
def Assignment_RZ(Alloc,N,C,ReuzeZoneSizes):
	root = Node(C[-1])
	C.remove(C[-1])
	tree = makeTree(root,C)	
	
	Assign = {}

	for e in ReuzeZoneSizes:
		if (ReuzeZoneSizes[e] > 0):
			Assign[e] = [(0,ReuzeZoneSizes[e] - 1)]
	
	subchannels = [0 for i in range(N)]
	q = Q.Queue()
	q.put(tree)
	U = list(Alloc.keys())
	while (q.empty() == False):
		curr_node = q.get()
		#print ("CLIQUE:",curr_node.data)

		# initially all channels are available
		subchannels = [0 for i in range(N)]

		# put the children of all current nodes in queue
		for ch in curr_node.children:
			q.put(ch)

		for v in curr_node.data:
			# color channels that are already assigned in this clique
			if v in Assign:
				subchannels = color_channel(Assign[v],subchannels)

		for v in curr_node.data:
			# check if current not hasn't been allocated its share
			if v in U:
				RZ_exists = False
				# if in U and also in Assign means its an RZ client
				if v in Assign:
					RZ = Assign[v][0]
					RZ_exists = True
				if (Alloc[v] > 0): # check added to incorporate nodes with 0 load
					Assign[v]=allocate_sub_chan(v,subchannels,Alloc[v]) 
					if (RZ_exists):
						Assign[v].append(RZ)
				U.remove(v)
				# if assigned anything to v color it, else assign it an empty list
				if v in Assign:
					subchannels = color_channel(Assign[v],subchannels)
				else :
					Assign[v] = []

	return Assign


def restore(v,free,subchannels):
	rec_blks = []
	for blk in free:
		start = find_unalloc(subchannels,blk[0])
		while(start <= blk[1] and start > -1):
			b = find_unalloc_block2(subchannels,start,blk[1])
			rec_blks.append(b)
			color_channel([b],subchannels)
			start = find_unalloc(subchannels,b[1]+1)
	return rec_blks
		

# Removes edges added to make the graph chordal
# Step 3 in Algorithm 1, and explained in text
# Input:
#   Assign - assignment
#   fill_in - added edges
#   i_map - original interference map
#   N - number of sub-channels
# Output:
#   Assign - updated assignment
def Restoration(Assign,fill_in,i_map,N):
	for e in fill_in:
		# for e[0]
		subchannels = [0 for i in range(N)]
		subchannels = color_channel(Assign[e[0]],subchannels)
		for v in i_map[e[0]]:
			subchannels = color_channel(Assign[v],subchannels)
		restored_bw = restore(e[0],Assign[e[1]],subchannels) 
		if len(restored_bw) != 0:
			for b in restored_bw:
				Assign[e[0]].append(b)
		# for e[1]
		subchannels = [0 for i in range(N)]
		subchannels = color_channel(Assign[e[1]],subchannels)
		for v in i_map[e[1]]:
			subchannels = color_channel(Assign[v],subchannels)
		restored_bw = restore(e[1],Assign[e[0]],subchannels) 
		if len(restored_bw) != 0:
			for b in restored_bw:
				Assign[e[1]].append(b)
	return Assign


def CanTakeChannel(ch,candidate,G,curr):
	#print G[candidate]
	#print curr
	for n in G[candidate]:
			for interval in curr[n]:
				if (ch>=interval[0] and ch<=interval[1]):
					return False
	return True

def freeSC(subchannels):
	free = []
	for i in range(len(subchannels)):
		if (subchannels[i] == 0):
			free.append(i)
	return free
# reclaim sub channels that are unused in a clique
# this can happen in interference map of UEs with class 1 clients
def ReclaimSC(Assign,N,C,IsClass1,i_map):
	root = Node(C[-1])
	C.remove(C[-1])
	tree = makeTree(root,C)	
	
	Reclaim = deepcopy(Assign)
	
	subchannels = [0 for i in range(N)]
	q = Q.Queue()
	q.put(tree)
	#U = list(Alloc.keys())
	while (q.empty() == False):
		curr_node = q.get()
		subchannels = [0 for i in range(N)]
		for ch in curr_node.children:
			q.put(ch)
		for v in curr_node.data:
			if v in Reclaim:
				subchannels = color_channel(Reclaim[v],subchannels)

		class1ForthisClique = []
		for v in curr_node.data:
			if (IsClass1[v]):
				class1ForthisClique.append(v)
		if (sum(subchannels) > 0 and len(class1ForthisClique) > 0): # some sub channels are free in this clique
			availableSC = freeSC(subchannels)
			num = 0
			for ch in availableSC:
				for count in range (len(class1ForthisClique)):
					if (CanTakeChannel(ch,class1ForthisClique[num],i_map,Reclaim)):
						Reclaim[class1ForthisClique[num]].append((ch,ch))
						num += 1
						num %= len(class1ForthisClique)
						break

	return Reclaim


def getCliques(i_map):
	(i_map_,fill_in) = Triangulation (deepcopy(i_map))
	C = cliques(i_map_)
	for c in C:
		c.sort()
	C.sort()
	return (i_map,i_map_,fill_in,C)
# Main function

def FermiPreCompute(i_map,load,N,i_map_,fill_in,C):
	#print (C)
	Alloc = Allocate(i_map_,load,N,deepcopy(C))

	#print(Alloc)
	Assign = Assignment(Alloc,N,deepcopy(C))

	Res = Restoration(Assign,fill_in,i_map,N)

	return (Res,Alloc)

def Fermi(i_map,load,N):

	(i_map_,fill_in) = Triangulation (deepcopy(i_map))
	C = cliques(i_map_)
	for c in C:
		c.sort()
	C.sort()

	#print (C)
	Alloc = Allocate(i_map_,load,N,deepcopy(C))

	#print(Alloc)
	Assign = Assignment(Alloc,N,deepcopy(C))

	Res = Restoration(Assign,fill_in,i_map,N)

	return (Res,Alloc)

# Main function allocations with reuze zone
def Fermi_RZ(i_map,load,N,ReuzeZoneSizes):


	(i_map_,fill_in) = Triangulation (deepcopy(i_map))
	C = cliques(i_map_)
	Alloc = Allocate(i_map_,load,N,deepcopy(C))


	#Adjusting Allocations to incorporate Reuse Zone
	for e in Alloc:
		Alloc[e] = int(math.floor(float(Alloc[e])/N * (N - ReuzeZoneSizes[e]) + 0.5))

	RZ_Sizes = {}
	for e in ReuzeZoneSizes:
		if e in i_map:
			RZ_Sizes[e] = ReuzeZoneSizes[e]

	Assign = Assignment_RZ(Alloc,N,deepcopy(C),RZ_Sizes)


	Res = Restoration(Assign,fill_in,i_map,N)

	return (Res,Alloc)

# Main function allocations where interference map is between UEs instead of eNBs
def Fermi_UE(i_map,load,N,IsClass1,RZSizes):

	#print (i_map)

	(i_map_,fill_in) = Triangulation (deepcopy(i_map))
	C = cliques(i_map_)
	Alloc = Allocate(i_map_,load,N,deepcopy(C))

	#print (Alloc)

	#Adjusting Allocations to incorporate Reuse Zone
	for e in Alloc:
		Alloc[e] = int(math.floor(float(Alloc[e])/N * (N - RZSizes[e]) + 0.5))

	#print ("Load: ",load)
	#print ("Alloc",Alloc)
	#for u in IsClass1:
	#	if not IsClass1[u]:
	#		RZSizes[u] = 0
			

	RZ_Sizes = {}
	for e in RZSizes:
		if e in i_map:
			RZ_Sizes[e] = RZSizes[e]
	#print ("RZSizes",RZSizes)
	Assign = Assignment_RZ(Alloc,N,deepcopy(C),RZ_Sizes)

	#print("Fill in:",fill_in)
	#print("Clique Tree:",C)
	#print("Assignment:",Assign)

	Reclaimed = ReclaimSC(Assign,N,deepcopy(C),IsClass1,i_map)
	#print(Reclaimed)
	Res = Restoration(Reclaimed,fill_in,i_map,N)
	
	for u in IsClass1:
		if not IsClass1[u] and u in i_map:
			if (0,RZSizes[u]-1) in Res[u]:
				Res[u].remove((0,RZSizes[u]-1))

	return (Res,Alloc)


def Fermi_Pre_Alloc(i_map,load,N,IsClass1,RZSizes,Alloc):

	#print (i_map)

	(i_map_,fill_in) = Triangulation (deepcopy(i_map))
	C = cliques(i_map_)
	#Alloc = Allocate(i_map_,load,N,deepcopy(C))

	#print (Alloc)

	#Adjusting Allocations to incorporate Reuse Zone
	#for e in Alloc:
	#	Alloc[e] = int(math.floor(float(Alloc[e])/N * (N - RZSizes[e]) + 0.5))

	#for u in IsClass1:
	#	if not IsClass1[u]:
	#		RZSizes[u] = 0
			

	RZ_Sizes = {}
	for e in RZSizes:
		if e in i_map:
			RZ_Sizes[e] = RZSizes[e]
			if (IsClass1[e]):
				Alloc[e] -= RZ_Sizes[e]

	Assign = Assignment_RZ(Alloc,N,deepcopy(C),RZ_Sizes)

	#print("Fill in:",fill_in)
	#print("Clique Tree:",C)
	#print("Assignment:",Assign)

	Reclaimed = Assign;#ReclaimSC(Assign,N,deepcopy(C),IsClass1,i_map)
	#print(Reclaimed)
	Res = Assign;#Restoration(Reclaimed,fill_in,i_map,N)
	
	for u in IsClass1:
		if not IsClass1[u] and u in i_map:
			if (0,RZSizes[u]-1) in Res[u]:
				Res[u].remove((0,RZSizes[u]-1))

	return (Res,Alloc)

def Fermi_Pre_Alloc2(i_map,load,N,Alloc,ReuzeZoneSizes):


	(i_map_,fill_in) = Triangulation (deepcopy(i_map))
	C = cliques(i_map_)
	#Alloc = Allocate(i_map_,load,N,deepcopy(C))


	#Adjusting Allocations to incorporate Reuse Zone
	for e in Alloc:
		#Alloc[e] = int(math.floor(float(Alloc[e])/N * (N - ReuzeZoneSizes[e]) + 0.5))
		Alloc[e] = max(0,Alloc[e]- ReuzeZoneSizes[e])

	RZ_Sizes = {}
	for e in ReuzeZoneSizes:
		if e in i_map:
			RZ_Sizes[e] = ReuzeZoneSizes[e]

	Assign = Assignment_RZ(Alloc,N,deepcopy(C),RZ_Sizes)


	Res = Assign;#Restoration(Assign,fill_in,i_map,N)

	return (Res,Alloc)

# G = {0:[1,3],1:[2,0],2:[3,1],3:[0,2]}
# (g_,f) = Triangulation2(G)
# C = cliques(g_)
# print g_

