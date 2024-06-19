
import ast
import sys
import math
import numpy as np #v 1.23.3
import networkx as nx
import argparse
import itertools
import pandas as pd
import os
import scipy
from numpy import roots
import random
import pickle
#os.chdir('/home/ll16598/Documents/ARCHITECTURAL_IMMUNITY/ANALYSIS')

import scipy
from numpy import roots
import random
import pickle

import sys

list_paths_skl=[]
list_paths_vol=[]
usr='cluster'
mode='real'#'test'
if usr=='cluster':
    directory="/user/work/ll16598"
    directory_G=directory
    array_id = int(sys.argv[1])-1
    iteration = int(sys.argv[2])
    timestamp = int(sys.argv[3])
    print('Array=',array_id)
else:
    directory='/home/ll16598/Documents/ARCHITECTURAL_IMMUNITY'
    directory_G='/home/ll16598/Documents/ARCHITECTURAL_IMMUNITY/CURRENT_Gs/01-05-23_2/GRAPHS'
    array_id = 1
    iteration = 1
    timestamp = 1
cham_en_key=pd.read_csv(directory+'/abm_cham_en_key.csv')
TREATS = pd.read_csv(directory+'/COLONY_INFO.csv')
dir_G=directory_G+"/BASE"
dir_G_widths_o=directory_G+'/WIDTH_O'
nurse_k = 0
forag_k = 0
nurse_ne_k = 0
forag_ne_k = 0



#TREATS = pd.read_csv('/media/cf19810/One Touch/CT_ANALYSIS/COLONY_INFO.csv')
analysis_df=TREATS
name_list=[]
for i in TREATS['name']:
    name_list.append(i)


#For network creation:
def process_interactions(interactions_list):
    contacts_dict = {}
    for interaction in interactions_list:
        # Sort the tuple to treat (A, B) the same as (B, A)
        sorted_interaction = tuple(sorted(interaction))
        
        # Update the contacts count
        if sorted_interaction in contacts_dict:
            contacts_dict[sorted_interaction] += 1
        else:
            contacts_dict[sorted_interaction] = 1
            
    return [(edge, weight) for edge, weight in contacts_dict.items()]


# In[ ]:


def convert_edge_format(edges):
    new_edges = []
    for edge in edges:
        new_edges.append(tuple([edge[0][0], edge[0][1], edge[1]]))
    return new_edges
    

files_width = [ f.path for f in os.scandir(dir_G_widths_o)]
files = [ f.path for f in os.scandir(dir_G)]


# In[9]:


def calculate_distance (node1, node2):
    np_node1=np.array(node1)
    np_node2=np.array(node2)
    squared_dist = np.sum((np_node1-np_node2)**2, axis=0) #this calculates distance between two points
    dist = np.sqrt(squared_dist)
    return dist




#Graph set-up
G_list=[]
all_names=[]
l=0
while l<len(files) and len(G_list)<100:

    for file in range(0,len(files)):
        day=os.path.basename(files[file])
        filename = "_".join(day.split("_")[:2])  # Split by underscore, take first two parts, join with underscore
        #print(file)
        if l==len(files):
            break
        if filename == name_list[l]:
            G=nx.read_graphml(files[file])
            Gs = sorted(nx.connected_components(G), key=len, reverse=True)
            Gmax = G.subgraph(Gs[0])
            G_list.append(Gmax)
            print(filename)
            all_names.append(filename)
            l+=1
        else:
            l+=0
            
G_list_width=[]
all_names=[]
l=0
while l<len(files) and len(G_list_width)<100:
    for file in range(0,len(files)):
        day=os.path.basename(files_width[file])
        filename = "_".join(day.split("_")[:2])  # Split by underscore, take first two parts, join with underscore
        #print(file)
        if l==len(files):
            break
        if filename == name_list[l]:
            G=nx.read_graphml(files_width[file])
            Gs = sorted(nx.connected_components(G), key=len, reverse=True)
            Gmax = G.subgraph(Gs[0])
            G_list_width.append(Gmax)
            print(filename)
            all_names.append(filename)
            l+=1
        else:
            l+=0



#Ant class
class Ant:
    def __init__(self, spores, caste, identity, edge, t, prev_node, next_node,identity2):
        #self.orientation=orientation
        self.spores=spores
        self.caste=caste
        self.identity=identity
        self.edge=edge#just basic edge
        self.t=t
        self.prev_node=prev_node
        self.next_node=next_node
        self.identity2=identity2
        
#Transmission functions
def draw_from_poisson(lambda_poisson):
    return np.random.poisson(lambda_poisson)
# In[17]:

def get_free_ant_encounters(N_agents):
    lambda_rate_per_agent = 0.0208  # Average interaction rate per agent per second
    #total_pairs = N_agents * (N_agents - 1) / 2
    Lambda_total = lambda_rate_per_agent * N_agents
    expected_interactions = Lambda_total
    return draw_from_poisson(Lambda_total)



def transmit(n, lst, T):
    
    max_attempts = 200  # You can adjust this based on your needs
    attempts = 0

    usage_count = {item: 0 for item in lst}

    # Dictionary to keep track of how many times each item has been used
    usage_count = {item.identity: 0 for item in lst}

    # Set to keep track of pairs that have already been formed
    formed_pairs = set()

    pairs = []

    while len(pairs) < n and len(pairs)<len(lst) and attempts < max_attempts:

        # Shuffle the list
        random.shuffle(lst)
        attempts+=1
        if attempts>=max_attempts:
            print('attemptmax')
            return lst, pairs
        # Try forming a pair
        ant, ant2 = lst[0], lst[1]

        # Ensure we don't have duplicate items in a pair and form pairs with distinct items
        if ant != ant2:
            # Convert pair to a frozenset so order doesn't matter, e.g., (1, 2) is same as (2, 1)
            current_pair = frozenset([ant.identity, ant2.identity])
            #print(current_pair)
            # Check if pair hasn't been formed before and items haven't been used more than twice
            if current_pair not in formed_pairs and usage_count[ant.identity] < 1 and usage_count[ant2.identity] < 1:
                pairs.append(tuple([ant.identity, ant2.identity, ant.caste, ant2.caste, T]))
                formed_pairs.add(current_pair)
                usage_count[ant.identity] += 1
                usage_count[ant2.identity] += 1
                #this is where social net would be added
                if random.random()>B:
                    continue
                #print('transmitted')
                if ant.spores>ant2.spores and ant.spores>sigma:
                    load_difference=ant.spores-ant2.spores
                    if load_difference<=0:
                        continue
                    transmitted=transmission_rate*load_difference*dt
                    ant2.spores+=transmitted
                    ant.spores-=transmitted
                elif ant2.spores>ant.spores and ant2.spores>sigma:
                    load_difference=ant2.spores-ant.spores
                    if load_difference<=0:
                        continue
                    transmitted=transmission_rate*load_difference*dt
                    ant.spores+=transmitted
                    ant2.spores-=transmitted
    return lst, pairs

caste=['Nurse', 'Forager', 'Inoculated_Forager', 'Inoculated_Nurse']
orientation=np.nan
A=0#identity
edge=np.nan
prev_node=np.nan
next_node=np.nan
nurse=Ant(0, caste[0], A, edge, 0, prev_node, next_node,0)
forager=Ant(0, caste[1], A, edge, 0, prev_node, next_node,0)
inoculated=Ant(1, caste[2], A, edge, 0, prev_node, next_node,0)




#Agents
def create_agents(nurse_number, forager_number, inoculated_forager_number, inoculated_nurse_number):
    nurses=[]
    foragers=[]
    inoculated_foragers=[]
    inoculated_nurses=[]
    A=0
    edge=np.nan
    prev_node=np.nan
    next_node=np.nan
    for a in range(0, nurse_number):
        nurse=Ant(0, caste[0], A, edge, 0, prev_node, next_node,0)
        nurses.append(nurse)
        A+=1

    for a in range(0, forager_number):
        forager=Ant(0, caste[1], A, edge, 0, prev_node, next_node,0)
        foragers.append(forager)
        A+=1

    for a in range(0, inoculated_forager_number):
        inoculated=Ant(1, caste[2], A, edge, 0, prev_node, next_node,0)
        inoculated_foragers.append(inoculated)
        A+=1
        # no inoculated nurses
    for a in range(0, inoculated_nurse_number):
        inoculated=Ant(1, caste[2], A, edge, 0, prev_node, next_node,0)
        inoculated_nurses.append(inoculated)
        A+=1
        #we re
    return nurses, foragers, inoculated_foragers, inoculated_nurses




nurse_number=180
forager_number = 180 - nurse_number
inoculated_nurse_number= 0
inoculated_forager_number = 20 - inoculated_nurse_number

transmission_rate=0.00138
ant_movement_speed=0.2#cm/s half a body length per second
ants_cm=0.4
ant_length=0.4
random_entrance=True
junction_capacity1=0.9333236#0.3374093
junction_capacity2=0.9022186#0.2759298
junction_num=2
max_interact=2
capacity_j=junction_num
capacity_e=1#/ants_cm
chamber_capacity1=1.513027
chamber_capacity2=1.450666
B= np.mean(list([0.39,0.36, 0.23])) #(probability of contact 
epidemic_start=21600
dt=1#timestep
Tmax=86400-21600
sigma=0.00024#transmission threshold
rando_stasis=False
if mode=='test':
    Tmax=200
    epidemic_start=100



XMID= 5.564687499999999
YMID= 5.579222039473684
INITIAL_NE=tuple([XMID,YMID,0])


# In[21]:


wed_seq = [96 - i*5 for i in range(20)]
wed_seq.reverse()
mon_seq = [99 - i*5 for i in range(20)]
mon_seq.reverse()
WED_MON_Gs=list(wed_seq+mon_seq)


# In[23]:


all_nest_nodes = []
all_chamber_nodes = []
all_initial_nes=[]
all_junction_nodes = []
all_end_nodes = []
all_edges=[]
all_nodes=[]
for g in range(0,len(G_list)):
#for g in range(0,len(wed_seq)):
    G=G_list[g]
    G2=G_list_width[g]
    

        
    #list_ant=[]
    attributes = nx.get_node_attributes(G, 'TYPE')
    nest_nodes = [node for node, type_value in attributes.items() if type_value == 'NEST EN']
    min_dist=10
    nest_nodes2=[]
    initial_nes=[]
    for ne in nest_nodes:
        x,y,z=ast.literal_eval(G.nodes[ne]['coord'])
        
	
        xy=tuple([x,y,0])
        ne2=tuple([ne,xy])
        nest_nodes2.append(ne2)
        dist=calculate_distance(xy,INITIAL_NE)
        if dist<min_dist:
            initial_ne=tuple([ne,xy])
            min_dist=dist
    initial_nes.append(initial_ne)
    nest_nodes2=list(set(nest_nodes2)-set(initial_nes))
    chamber_nodes = [node for node, type_value in attributes.items() if type_value == 'CHAM']
    junction_nodes = [node for node, type_value in attributes.items() if type_value == 'JUNC']
    end_nodes = [node for node, type_value in attributes.items() if type_value == 'END']
    if len(chamber_nodes)==0:
        chamber_nodes=np.nan


#     nest_nodes = [(node, list_ant) for node in nest_nodes]
#     chamber_nodes = [(node, list_ant) for node in chamber_nodes]
#     junction_nodes = [(node, list_ant) for node in junction_nodes]
#     end_nodes = [(node, list_ant) for node in end_nodes]

    all_initial_nes.append(initial_nes)
    all_nest_nodes.append(nest_nodes2)
    all_chamber_nodes.append(chamber_nodes)
    all_junction_nodes.append(junction_nodes)
    all_end_nodes.append(end_nodes)
    try:
        node_list = chamber_nodes + initial_nes+nest_nodes+ end_nodes+ junction_nodes
    except TypeError:
        node_list = initial_nes+nest_nodes+ end_nodes+ junction_nodes

    all_nodes.append(node_list)
    edge_weights = []
    for u, v, data in G.edges(data=True):
        length = data['weight']
        width = G2.get_edge_data(u, v)['weight']
        if u in nest_nodes:
            x,y,z=ast.literal_eval(G.nodes[u]['coord'])
            xy=tuple([x,y,0])
            u=tuple([u,xy])
        if v in nest_nodes:
            x,y,z=ast.literal_eval(G.nodes[v]['coord'])
            xy=tuple([x,y,0])
            v=tuple([v,xy])   
        list_coords=[]
        #This is subject to change depending on if cells or real coords are used
        l=0
        length_mm=round(length*10)
        for i in range(0, length_mm):
            list_coords.append(l)
            l+=1
        edge_weights.append((tuple([u, v]), length, width))
       # print(length,width)
    all_edges.append(edge_weights)
    



# In[ ]:


DIR_PICK_CHTIME=directory+'/CHTIME_PICKLES'
DIR_PICK_EDGE=directory+'/EDGE_PICKLES'
DIR_PICK_NODE=directory+'/NODE_PICKLES'
DIR_PICK_INT=directory+'/INTERACTION_PICKLES'
DIR_PICK_GRAPHS=directory+'/GRAPH_OUTPUT'
DIR_PICK_INT2=directory+'/INTERACTION_ANT_SPACE_PICKLES'
# In[ ]:

strits=str(epidemic_start)
strarr=str(Tmax)
stramp=str(timestamp+round(random.random()*10))

import copy
ALL_ITS_NODE_ANTS=[]
ALL_ITS_EDGE_ANTS=[]

min_cham_time_list=[]
ALL_NODE_ANTS=[]
ALL_EDGE_ANTS=[]
ALL_NODE_ANTS=[]
ALL_EDGE_ANTS=[]

ALL_INTERACTIONS=[]
ALL_INTERACTIONS_ANT_SPACE=[]
t_save1=[i for i in range(0,epidemic_start,10000)]
t_save2=[i for i in range(epidemic_start, epidemic_start+5000,40)]
t_save3=[i for i in range(epidemic_start+5000, Tmax+100,2000)]
t_save=t_save1+t_save2+t_save3
t_save=t_save1+t_save2+t_save3
#an=int(sys.argv[2])#160#80#20#0#chamber_attractiveness
#ps=int(sys.argv[1])#16#4#pheromone deposited
#py=50#staying factor
#decay_rate=0.05
#max_pheromone=int(sys.argv[6])#550#225
def leaving_probability(pn):
    leaving=False
    try:
        leaving_probability=0.5/(pn+1)
    except ZeroDivisionError:
        leaving_probability=1
    if random.random()<leaving_probability or leaving_probability<0:
        leaving=True
    return leaving
    
g=array_id
kk=WED_MON_Gs[g]
for iterr in range(0, 10):
#wednesday or monday chamber volume
    if kk in wed_seq:
        chamber_capacity=chamber_capacity1
        junction_capacity=junction_capacity1
        max_out_prop=1#0.1

    else:
        chamber_capacity=chamber_capacity2
        junction_capacity=junction_capacity2
        max_out_prop=1#0.1
    chamber_list=all_chamber_nodes[kk]
    if chamber_list is np.nan:
        chamber_list=[]

	#Set-up simulation

    print('SIMULATION ON = ',name_list[kk])

    min_cham_time=Tmax
    ne_list=all_nest_nodes[kk]
    initial_ne_list=all_initial_nes[kk]
    end_list=all_end_nodes[kk]
    junction_list=all_junction_nodes[kk]
    edge_list=all_edges[kk]
    name=name_list[kk]
    G=G_list[kk]
    nurses, foragers, inoculated_foragers, inoculated_nurses=create_agents(nurse_number, forager_number, \
                                                                           inoculated_forager_number, inoculated_nurse_number)
    inoculateds=list(inoculated_foragers+ inoculated_nurses)        
    original_ants=list(nurses+foragers)
    max_outside=180*max_out_prop
    G_NODE_ANTS=[]
    G_EDGE_ANTS=[]
    G_INTERACTIONS=[]
    G_INTERACTIONS_ANT_SPACE=[]
    g_transmission_chamber=[]
    g_transmission_tunnel=[]
    g_transmission_junction=[]
    reached_chamber=False
    #SET-UP THE SIMULATION BY ALLOWING ANTS TO SETTLE IN THE NEST AND RUN FOR A TIME
    if isinstance(chamber_list, float) and np.isnan(chamber_list):
        chamber_list = []
    elif not chamber_list:
        chamber_list = []

    node_list = chamber_list + ne_list + end_list + junction_list+initial_ne_list
    nech_list = chamber_list + ne_list+initial_ne_list
    inodes=ne_list+initial_ne_list
    edge_ants={edge: [] for edge in edge_list}
    node_ants = {node: [] for node in node_list} # create an empty list for each node
    pheromone_dict = {node: 0 for node in chamber_list} 
    d1_list_ne=[node for node in ne_list if len(list(G.neighbors(node[0])))==1]
    d1_list_ine=[node for node in initial_ne_list if len(list(G.neighbors(node[0])))==1]
    d1_list=end_list+d1_list_ne+d1_list_ine
    #break



 #Chamber assignment
    num_nes=len(inodes)
    num_juncs=len(junction_list)
    num_nes=len(inodes)
    num_chams=len(chamber_list)
    used_ants=[]
    if num_chams>0:
        num_nodes=len(node_list)
        prop_ch=num_chams/num_nodes
        num_juncs=len(junction_list)
        num_pc_ch=180*prop_ch
        number_across_nest=180-num_pc_ch
        num_per_ch=round(num_pc_ch/num_chams)
        if num_per_ch<1:
            num_per_ch=1
        if num_pc_ch<(num_per_ch*num_chams):
            num_pc_ch=(num_per_ch*num_chams)

        #num_pj=math.ceil((number_across_nest)/num_juncs)

        num_in_cham=0
        
        while num_in_cham<num_pc_ch:

            for chamber in chamber_list:
                num_in_this_cham=0
                while num_in_this_cham<num_per_ch:
                    for ant in original_ants:
                        if ant.identity in used_ants:
                            continue
                        else:
                            node_ants[chamber].append(ant)
                            used_ants.append(ant.identity)
                            num_in_cham+=1
                            num_in_this_cham+=1
                            no=random.choice(list(G.neighbors(chamber)))
                            for nod in node_list:
                                if nod[0]==no:
                                        ant.next_node=nod
                                        break

                            if num_in_cham>=num_pc_ch:
                                break
                            if num_in_this_cham>=num_per_ch:
                                break
                if num_in_cham>=num_pc_ch:
                    break
        for node in chamber_list:
            print('n in cham:',len(node_ants[node]),num_per_ch)
    node_list2=[node for node in node_list if node not in chamber_list]   

    for ant in original_ants:
        if ant.identity in used_ants:
            continue
        chosen_ne = random.choice(node_list2)
        ant.prev_node=chosen_ne
        if chosen_ne in inodes:
            ant.next_node=random.choice(list(G.neighbors(chosen_ne[0])))
        else:
            no=random.choice(list(G.neighbors(chosen_ne)))
            for nod in node_list:
                if nod[0]==no:
                        ant.next_node=nod
                        break
        node_ants[chosen_ne].append(ant) 

	#Main simulation
    T=0
    while T<Tmax:
        #Placing inoculated ants in entrances at the start time
        if T==epidemic_start:
            id2=0
            for node in chamber_list:
                ants=node_ants[node]
                for ant in ants:
                    ant.identity2=id2
                    id2+=1        
            node_list_no_cham=[node for node in node_list if node not in chamber_list]
            for node in node_list_no_cham:
                for ant in node_ants[node]:
                    ant.identity2=id2
                    id2+=1
            for edge in edge_list:
                ants=edge_ants[edge]
                for ant in ants:
                    ant.identity2=id2
                    id2+=1
            max_outside=200*max_out_prop
            chosen_ne = random.choice(inodes)
            tss=0
            for ant in inoculateds:
                ant.t=T+tss
                tss+=4 #delay next ant by 4s
                ant.prev_node=chosen_ne
                ant.next_node=random.choice(list(G.neighbors(chosen_ne[0])))
                node_ants[chosen_ne].append(ant) # add the ant to the list for the chosen node
        ##FIRST STEP ADVANCING ANTS FROM NODES
        random.shuffle(node_list)
        for node in node_list:
            ants=node_ants[node]
            #print(node_ants)
            random.shuffle(ants)
            choices=[]
            for ant in ants:

                ant.prev_node=node
                leaving=False
                if ant.t>=T: #only ants not 'moving'
                    continue                
                elif node in ne_list and ant.t<T: 
                    leaving=True
                elif node in initial_ne_list and ant.t<T: 
                    leaving=True
                elif node in chamber_list and ant.t<T:
                    leaving=leaving_probability(len(ants))
                elif node in junction_list and ant.t<T:
                    leaving=True
                elif node in d1_list and ant.t<T:
                    leaving=True
                    for edge in edge_list:
                        if node==edge[0][1]:
                            ant.prev_node=edge[0][1]
                            ant.next_node=edge[0][0]
                            ant.t=(T+(edge[1]/ant_movement_speed))
                            chosen_edge=edge
                            break
                        elif node==edge[0][0]:
                            ant.prev_node=edge[0][0]
                            ant.next_node=edge[0][1]
                            ant.t=(T+(edge[1]/ant_movement_speed))
                            chosen_edge=edge
                            break
                    edge_ants[chosen_edge].append(ant)
                    node_ants[node] = [a for a in node_ants[node] if a.identity != ant.identity]  # remove the ant by ID
                    continue
                if leaving==False:
                    continue
                   #Ant is leaving node (if edges are available)
                if leaving==True: 
                #edge choice
                    choices=[]
                    edge_choices=[e for e in edge_list if e[0][1]==ant.prev_node or e[0][0]==ant.prev_node]
                    for edge in edge_choices:  
                        capacity=max((edge[2] * edge[1]) / ants_cm, 1)
                        e_ants=edge_ants[edge]
                        if len(e_ants)>=capacity:
                            continue
                        else:
                            choices.append(edge)
                    if len(choices)==0: #is tunnel full?
                        continue
                    chosen_edge=random.choice(choices)

                    ant.prev_node=node
                    if chosen_edge[0][1] == node:
                        ant.next_node = chosen_edge[0][0]
                    else:
                        ant.next_node = chosen_edge[0][1]

                    ant.t=(T+(chosen_edge[1]/ant_movement_speed))
                    edge_ants[chosen_edge].append(ant) #ant moves t edge
                    node_ants[node] = [a for a in node_ants[node] if a.identity != ant.identity]

        random.shuffle(edge_list)
        for edge in edge_list:
            ants=edge_ants[edge]
            outside_num=np.sum([len(node_ants[node]) for node in inodes]) #len([ant for ne_node in node_list for \
                     #ant in node_ants[ne_node]]) 

            for ant in ants:
                if ant.t>=T:
                   # print('stopped')
                    continue                
                else:
                    num_ants=len(node_ants[ant.next_node])
                    n=ant.next_node
                    if  num_ants>=capacity_j and n in junction_list:
                        if random.uniform(0,1)<0.5:
                            next_no=ant.prev_node
                            prev_no=ant.next_node
                            ant.next_node=next_no
                            ant.prev_no=prev_no
                            ant.t=(T+(edge[1]/ant_movement_speed))
                            continue
                        else:
                            continue
                    elif num_ants>=capacity_e and n in end_list:
                        if random.uniform(0,1)<0.5:
                            next_no=ant.prev_node
                            prev_no=ant.next_node
                            ant.next_node=next_no
                            ant.prev_no=prev_no
                            ant.t=(T+(edge[1]/ant_movement_speed))
                            continue
                        else:
                            continue                            #turn if full

                    elif n in ne_list:
                        n1=n
                        n=random.choice(ne_list)
                        ant.next_node=random.choice(list(G.neighbors(n[0])))
                        dist=calculate_distance(n1[1],n[1])
                        ant.t=(T+(dist/ant_movement_speed))
                    elif n in initial_ne_list:
                        ant.prev_node=n
                        ant.next_node=random.choice(list(G.neighbors(n[0])))
                        dist=random.uniform(0, 8.6) #uniform distribution to edge of petri and back
                        ant.t=(T+(dist/ant_movement_speed))
                    ant.prev_node=n
                    node_ants[n].append(ant)
                    edge_ants[edge] = [a for a in edge_ants[edge] if a.identity != ant.identity]  # remove the ant by ID

        #TRANSMISSION
        if T<epidemic_start and T in t_save:
            node_ants2=copy.deepcopy(node_ants)
            edge_ants2=copy.deepcopy(edge_ants)
            G_NODE_ANTS.append(node_ants2)
            G_EDGE_ANTS.append(edge_ants2)
        if T<epidemic_start:
            T+=dt
            continue
        for node in node_list:
            #transmission cannot take place in entrances or ends
            if node in ne_list or node in initial_ne_list or node in end_list:
                continue
            elif node in chamber_list:
                capacity=chamber_capacity
            ants=node_ants[node]
            ants_in_node=len(ants)
            n_encounters=get_free_ant_encounters(ants_in_node)
            if node in junction_list:
                n_encounters=1#=junction_capacity#change to 1
            if ants_in_node<=1:
                continue
            ants, interacting_pairs=transmit(n_encounters, ants, T)
            node_ants[node]=ants
            G_INTERACTIONS.extend(interacting_pairs)
        #edge transmission            
        for edge in edge_list:
            capacity=edge[2]*edge[1]
            ants=edge_ants[edge]
            ants_in_edge=len(ants)
            #only interact if edge is at capacity
            if ants_in_edge<=1:
                continue
            n_encounters=get_free_ant_encounters(ants_in_edge)
            #print(T, 'transmission in edge')
            ants, interacting_pairs=transmit(n_encounters, ants, T)
            edge_ants[edge]=ants
            G_INTERACTIONS.extend(interacting_pairs)

        T+=dt  
        if T in t_save:
            node_ants2=copy.deepcopy(node_ants)
            edge_ants2=copy.deepcopy(edge_ants)
            G_NODE_ANTS.append(node_ants2)
            G_EDGE_ANTS.append(edge_ants2)






    randooo=str(random.uniform(0,1000))
    OUT_NODE=DIR_PICK_NODE+'/_IT_'+name+'_'+strits+'_'+stramp+randooo+'.pickle'
    with open(OUT_NODE, "wb") as f:
           pickle.dump(G_NODE_ANTS, f)

    OUT_EDGE=DIR_PICK_EDGE+'/_IT_'+name+'_'+strits+'_'+stramp+randooo+'.pickle'
    with open(OUT_EDGE, "wb") as f:
           pickle.dump(G_EDGE_ANTS, f) 

print('ITERATION COMPLETED = ', stramp)

print('ITERATION COMPLETED = ', stramp)






# In[ ]:




