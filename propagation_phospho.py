from igraph import *
from scipy import sparse
import numpy as np
from scipy.stats import *
import os,fishertest,sys,scipy,math
from sklearn.neighbors import KernelDensity
from scipy.spatial.distance import jensenshannon

def denoise_square(G):
	weight=G.strength(G.vs, mode='all', loops=False, weights='weight')
	for i in G.es():
		node_A=i.tuple[0]
		node_B=i.tuple[1]
		den=math.sqrt(weight[node_A]*weight[node_B])
		num=i["weight"]
		G.es[i.index]["weight"]=num/den
	return (G)

def load_gene_names():
	uniprot_to_gene={}
	gene_to_uniprot={}
	localization={}
	f1=open("/nfs/research/petsalaki/shared_folder/diffusion/data/uniprot_to_gene.tab","r")
	seq=f1.readline()
	while(seq!=""):
		seq=seq.strip().split("\t")
		uniprot_to_gene[seq[0]]=seq[1].split(";")[0].strip()
		seq=f1.readline()
	f1=open("/nfs/research/petsalaki/shared_folder/diffusion/data/gene_to_uniprot.tab","r")
	seq=f1.readline()
	while(seq!=""):
		seq=seq.strip().split("\t")
		gene_to_uniprot[seq[0]]=seq[1:]
		seq=f1.readline()
	return uniprot_to_gene,gene_to_uniprot
uniprot_to_gene,gene_to_uniprot=load_gene_names()

def calc_kde(vector):
	obs = len(vector)
	sigma = np.std(vector, ddof=1)
	IQR = (np.percentile(vector, q=75) - np.percentile(vector, q=25)) / 1.3489795003921634
	sigma = min(sigma, IQR)
	bw=1.0
	if sigma > 0:
		bw= sigma * (obs * 3 / 4.0) ** (-1 / 5)
	else:
		IQR = (np.percentile(vector, q=99) - np.percentile(vector, q=1)) / 4.6526957480816815
		if IQR > 0:
			bw = IQR * (obs * 3 / 4.0) ** (-1 / 5)
	band=bw

	kde = KernelDensity(kernel = 'gaussian', bandwidth=bw).fit(vector.reshape(-1,1))
	grid = np.linspace(min(vector)-1,max(vector)+1,len(vector)*100).reshape(-1,1)
	log_dens = kde.score_samples(grid)
	pdf=np.exp(log_dens)
	grid=grid.ravel()
	normalization=sum(pdf.ravel())
	cdf=np.cumsum(pdf)/normalization
	return cdf,grid

def kinase_classification():
	tyr_kinase=[]
	st_kinase=[]
	f1=open("/nfs/research/petsalaki/shared_folder/diffusion/data/pfam_domains.txt")
	seq=f1.readline()
	seq=seq.strip().split("\t")
	tyr_kinase.extend(seq[1:])
	seq=f1.readline()
	seq=seq.strip().split("\t")
	st_kinase.extend(seq[1:])
	return tyr_kinase,st_kinase

def load_seeds(test,graph_nodes,sem_sim):
	mapping={"SimNTO":1,"TCSS":2,"Cosine":3,"KORBEL":4,"gic":5,"resnikbma":6,"resnikmax":7}

	f1=open("/nfs/research/petsalaki/shared_folder/diffusion/data/ssim_mean_std_all.txt")
	seq=f1.readline()
	seq=f1.readline()
	zscores_global={}
	while(seq!=""):
		seq=seq.strip().split("\t")
		zscores_global[seq[0]]=np.array(seq[mapping[sem_sim]].split("|"),dtype=float)
		seq=f1.readline()
	#change the datapath
	data_path="/nfs/research/petsalaki/users/sambor/pilot_project_mdd_bp_scz_networks-main/"
	f1=open(data_path+test)
	seq=f1.readline()
	ssim={}
	seeds_pos={}
	seeds_neg={}
	while(seq!=""):
		seq= seq.strip().split("\t")
		seq[0]=float(seq[0])
		if seq[0]>0.0:
			seeds_pos[seq[1]]=seq[0]
		if seq[0]<-0.0:
			seeds_neg[seq[1]]=-seq[0]
		seq=f1.readline()
	tyr_kinase,st_kinase=kinase_classification()
	tyr_pos=list((set(tyr_kinase).intersection(seeds_pos.keys())).intersection(set(graph_nodes)))
	tyr_neg=list((set(tyr_kinase).intersection(seeds_neg.keys())).intersection(set(graph_nodes)))
	st_pos= list((set(st_kinase).intersection(seeds_pos.keys())).intersection(set(graph_nodes)))
	st_neg=list((set(st_kinase).intersection(seeds_neg.keys())).intersection(set(graph_nodes)))
	all_the_rest_pos=list((set(seeds_pos.keys()).difference(set(tyr_pos+st_pos))).intersection(set(graph_nodes)))
	all_the_rest_neg=list((set(seeds_neg.keys()).difference(set(tyr_neg+st_neg))).intersection(set(graph_nodes)))

	mapping={"SimNTO":2,"TCSS":3,"Cosine":4,"KORBEL":5,"gic":6,"resnikbma":7,"resnikmax":8}
	seeds=[tyr_pos,st_pos,all_the_rest_pos,tyr_neg,st_neg,all_the_rest_neg]
	ssim={}
	for i in list(set(tyr_pos+tyr_neg+st_pos+st_neg+all_the_rest_pos+all_the_rest_neg)):
		ssim[i]={}
		fsim=open("/nfs/research/petsalaki/shared_folder/diffusion/data/chunk_results_all/"+i+"_all.txt")
		seq_sim=fsim.readline()
		seq_sim=fsim.readline()
		while(seq_sim!=""):
			seq_sim=seq_sim.strip().split("\t")
			if float(seq_sim[mapping[sem_sim]])>0:
				ssim[i][seq_sim[1]]=float(seq_sim[mapping[sem_sim]])
			else:
				ssim[i][seq_sim[1]]=0.00001
			seq_sim=fsim.readline()


	return seeds_pos,seeds_neg,seeds,zscores_global,ssim

def gen_val(network):
	number_of_nodes=network.vcount()
	rwr_random_values={}
	empirical_rwr=np.zeros((6,number_of_nodes),dtype=float)

	empirical_values={}
	pvalues={}
	for i in network.vs:
		pvalues[i["name"]]=np.zeros(6,dtype=int)
		empirical_values[i["name"]]=np.zeros(6,dtype=float)

	pos=False
	if ini_pos[0]!="False":
		pos=True
		to_delete=[]
		for j in network.vs.select(name_in=ini_pos):
			to_delete.append(j.index)

	neg=False
	if ini_neg[0]!="False":
		neg=True
		to_delete=[]
		for j in network.vs.select(name_in=ini_neg):
			to_delete.append(j.index)
	flag_pos=0
	flag_neg=0
	for i in enumerate(seeds):
		if len(i[1])>0:
			if pos==True and flag_pos==0 and i[0]<=2:
				network.delete_vertices(to_delete)
				flag_pos=1
			if neg==True and flag_neg==0 and i[0]>2:
				network.delete_vertices(to_delete)
				flag_neg=1

			number_of_nodes=network.vcount()
			reset_vertex=np.zeros(number_of_nodes)

			for j in network.vs.select(name_in=i[1]):
				if i[0]<=2:
					reset_vertex[j.index]=seeds_pos[j["name"]]
				else:
					reset_vertex[j.index]=seeds_neg[j["name"]]

			pagerank=np.array(network.personalized_pagerank(reset=reset_vertex,directed=False, damping=damping, weights='weight'))

			if flag_pos==1 and i[0]<=2:
				for j in to_delete:
					pagerank=np.concatenate((pagerank[:j], [0], pagerank[j:]))

			if flag_neg==1 and i[0]>2:
				for j in to_delete:
					pagerank=np.concatenate((pagerank[:j], [0], pagerank[j:]))

			empirical_rwr[i[0]]=pagerank

		if i[0]==2 and pos==True:
			network = Graph.Read_Ncol("/nfs/research/petsalaki/shared_folder/diffusion/networks/old_net/square/"+sim_type+".txt", weights=True, directed=False)
			pos=False
			flag_pos=0
			to_delete=[]


	for i in enumerate(empirical_rwr.T):
		empirical_values[graph_nodes[i[0]]]=i[1]

	for ii in range(1000):

		network_random = Graph.Read_Ncol("/nfs/research/petsalaki/shared_folder/diffusion/networks/old_net/square/"+sim_type+"_random/"+str(ii)+".txt", weights=True, directed=False)
		nodes=network_random.vs["name"]
		number_of_nodes=network_random.vcount()
		random_rwr=np.zeros((6,number_of_nodes),dtype=float)
		pos=False
		if ini_pos[0]!="False":
			to_delete=[]
			for j in network_random.vs.select(name_in=ini_pos):
				to_delete.append(j.index)
			pos=True
		neg=False
		if ini_neg[0]!="False":
			to_delete=[]
			for j in network_random.vs.select(name_in=ini_neg):
				to_delete.append(j.index)
			neg=True
		flag_pos=0
		flag_neg=0
		for jj in enumerate(seeds):
			if len(jj[1])>0:
				if pos==True and flag_pos==0 and jj[0]<=2:
					network_random.delete_vertices(to_delete)
					flag_pos=1
				if neg==True and flag_neg==0 and jj[0]>2:
					network_random.delete_vertices(to_delete)
					flag_neg=1

				number_of_nodes=network_random.vcount()
				reset_vertex=np.zeros(number_of_nodes)

				for j in network_random.vs.select(name_in=jj[1]):
					if jj[0]<=2:
						reset_vertex[j.index]=seeds_pos[j["name"]]
					else:
						reset_vertex[j.index]=seeds_neg[j["name"]]
				prandom=np.array(network_random.personalized_pagerank(reset=reset_vertex,directed=False, damping=damping, weights='weight'))

				if flag_pos==1 and jj[0]<=2:
					for j in to_delete:
						prandom=np.concatenate((prandom[:j], [0], prandom[j:]))

				if flag_neg==1 and jj[0]>2:
					for j in to_delete:
						prandom=np.concatenate((prandom[:j], [0], prandom[j:]))

				if jj[0]==2 and pos==True:
					network_random = Graph.Read_Ncol("/nfs/research/petsalaki/shared_folder/diffusion/networks/old_net/square/"+sim_type+"_random/"+str(ii)+".txt", weights=True, directed=False)
					flag_pos=0
					pos=False
					to_delete=[]

				random_rwr[jj[0]]=prandom
			else:
				if jj[0]==2:
					network_random = Graph.Read_Ncol("/nfs/research/petsalaki/shared_folder/diffusion/networks/old_net/square/"+sim_type+"_random/"+str(ii)+".txt", weights=True, directed=False)
					flag_pos=0
					pos=False
					to_delete=[]

		random_rwr=random_rwr.T
		for i in enumerate(nodes):
			pvalues[i[1]]=np.add(pvalues[i[1]],np.greater(empirical_values[i[1]],random_rwr[i[0]]))


	f1=open("/nfs/research/petsalaki/users/sambor/pilot_project_mdd_bp_scz_networks-main/"+res_folder+"/"+sim_type+"/"+str(damping)+"/"+test+"/pagerank.txt","w")
	for i in empirical_values:
		f1.write(i+"\t"+"\t".join(map(str,empirical_values[i]))+"\n")
	f1.close()

	f1=open("/nfs/research/petsalaki/users/sambor/pilot_project_mdd_bp_scz_networks-main/"+res_folder+"/"+sim_type+"/"+str(damping)+"/"+test+"/pvalues.txt","w")
	for i in pvalues:
		f1.write(i+"\t"+"\t".join(map(str,pvalues[i]))+"\n")
	f1.close()

	f1=open("/nfs/research/petsalaki/users/sambor/pilot_project_mdd_bp_scz_networks-main/"+res_folder+"/"+sim_type+"/"+str(damping)+"/"+test+"/start_seeds.txt","w")
	for i in seeds[0]+seeds[1]+seeds[2]:
		f1.write(i+"\t"+str(seeds_pos[i])+"\n")

	for i in seeds[3]+seeds[4]+seeds[5]:
		f1.write(i+"\t"+str(-seeds_neg[i])+"\n")
	f1.close()
	return empirical_values,pvalues

def write_results(nodes_kde,seed_nodes,rwr,direction):

	#exp=test.split("/")[1]
	exp=test #.split("/")
	folder="/nfs/research/petsalaki/users/sambor/pilot_project_mdd_bp_scz_networks-main/"+res_folder+"/"+sim_type+"/"+str(damping)+"/"+test
	nodes=seed_nodes
	for j in [0.75,0.8,0.85,0.9,0.95]:
		cluster_folder=folder+"/"+direction+"/cluster/"+str(j)
		fisher_folder_proteome=folder+"/"+direction+"/fisher/proteome/"+str(j)+"/"
		f2=open(cluster_folder+"/"+exp,"w")
		fisher_proteins=[]
		for k in nodes_kde[j]:
			temp=[]
			f2.write(k+"\t"+"\t".join(nodes_kde[j][k])+"\n")
			if len(nodes_kde[j][k])>=2:
				fisher_proteins=fisher_proteins+nodes_kde[j][k]
				for i in nodes_kde[j][k]:
					temp.append(uniprot_to_gene.get(i))
			fishertest.load(list(set(fisher_proteins)),0.05, ["C","F","D","P","R","K","RT","B"],fisher_folder_proteome,nodes,rwr,uniprot_to_gene,"proteome")

def ego_friends(subnet,position_nodes,seed_nodes,sim,zscores):
	zscore_threshold={"SimNTO":0.84, "TCSS":1.645,"resnikbma":1.645,"resnikmax":1.645,"KORBEL":0.84, "gic":1.64}
	nodes_kde={}
	for j in [0.75,0.8,0.85,0.9,0.95]:
		nodes_kde[j]={}

	for i in enumerate(sorted(seed_nodes)):
		ego=subnet.induced_subgraph(subnet.neighborhood(vertices=position_nodes[i[1]], order=2, mode='all', mindist=0),implementation='create_from_scratch')
		nodes=ego.vs["name"]
		n_nodes=len(nodes)
		if n_nodes>5:
			position_ego={k: v for v, k in enumerate(nodes)}
			functional_nodes=[position_ego[i[1]]]

			for j in nodes:
				num=sim[i[1]][j]-float(zscores[j][0])
				den=float(zscores[j][1])
				if den>0 and (num/den)>1.645:
					functional_nodes.append(position_ego[j])

			ego=ego.induced_subgraph(functional_nodes,implementation='copy_and_delete')
			ego=ego.induced_subgraph(ego.neighborhood(vertices=ego.vs.find(i[1]), order=2, mode='all', mindist=0),implementation='copy_and_delete')

			ego.vs.select(_degree=0).delete()
			nodes=ego.vs["name"]
			n_nodes=len(nodes)
			position_ego={k: v for v, k in enumerate(nodes)}

			if len(ego.clusters(mode='strong'))>1 and i[1] in position_ego:

				ego=ego.induced_subgraph(ego.subcomponent(position_ego[i[1]], mode='all'),implementation='copy_and_delete')
				nodes=ego.vs["name"]
				n_nodes=len(nodes)
				position_ego={k: v for v, k in enumerate(nodes)}

			if i[1] in position_ego and n_nodes>5:
				position_ego={k: v for v, k in enumerate(nodes)}
				first_shell=ego.neighbors(position_ego[i[1]], mode='all')

				distances=dict.fromkeys(first_shell,1)
				nodes_id=ego.vs.indices
				second_shell=list(set(nodes_id).difference(set(first_shell+[position_ego[i[1]]])))
				for j in second_shell:
					distances[j]=2
				distances[position_ego[i[1]]]=0
				for j in ego.es:
					node_A=j.tuple[0]
					node_B=j.tuple[1]
					dist=[distances[node_A],distances[node_B]]
					if (dist[0]==1 and dist[1]==1) or (dist[0]==2 and dist[1]==2):
						ego_sim= (sim[i[1]][nodes[node_A]]+sim[i[1]][nodes[node_B]])/2.0
						ego.es[j.index]["weight"]=ego_sim
					elif (dist[0]==1 and dist[1]==2):
						ego_sim=sim[i[1]][nodes[node_B]]
						ego.es[j.index]["weight"]=ego_sim
					elif (dist[0]==2 and dist[1]==1):
						ego_sim=sim[i[1]][nodes[node_A]]
						ego.es[j.index]["weight"]=ego_sim

				ego=denoise_square(ego)

				reset_vertex=np.zeros(n_nodes)
				reset_vertex[position_ego[i[1]]]=1.0
				ego_rwr=np.array(ego.personalized_pagerank(reset=reset_vertex,directed=False, damping=0.85, weights='weight'))
				ssim=[]
				dist=[]

				for j in ego.vs:
					reset_vertex=np.zeros(n_nodes)
					reset_vertex[j.index]=1.0
					ego_node=np.array(ego.personalized_pagerank(reset=reset_vertex,directed=False, damping=0.85, weights='weight'))
					ssim.append(sim[i[1]][j["name"]])
					dist.append(1-jensenshannon(ego_node,ego_rwr))

				ssim=np.array(ssim)
				dist=np.array(dist)

				ssim=1000*np.log2(1+ssim)
				dist=1000*np.log2(1+dist)

				nodes=np.array(nodes)
				cdf_dist,grid_dist=calc_kde(dist)
				cdf_ssim,grid_ssim=calc_kde(ssim)

				for j in [0.75,0.8,0.85,0.9,0.95]:
					if j not in nodes_kde:
						nodes_kde[j]={}

					position=np.searchsorted(cdf_dist,j)
					interpolation=np.argmin([abs(j-cdf_dist[position]),abs(j-cdf_dist[position-1])])
					grid_threshold_dist=grid_dist[position-interpolation]

					position=np.searchsorted(cdf_ssim,j)
					interpolation=np.argmin([abs(j-cdf_ssim[position]),abs(j-cdf_ssim[position-1])])
					grid_threshold_ssim=grid_ssim[position-interpolation]

					nodes_sim=[i[1]]
					nodes_dist=[i[1]]
					for k in enumerate(nodes):
						if ssim[k[0]]>grid_threshold_ssim:
							nodes_sim.append(k[1])
						if dist[k[0]]>grid_threshold_dist:
							nodes_dist.append(k[1])
					nodes_kde[j][i[1]]=list(set(nodes_sim).intersection(nodes_dist))
	return nodes_kde

def ego_filtering(network,rwr,pval,seeds,sim,zscores_global,direction):
	zscore_val={0.5:1.645,1.0:1.282,2:1.036}
	mapping={"SimNTO":1,"TCSS":2,"Cosine":3,"KORBEL":4,"gic":5,"resnikbma":6,"resnikmax":7}

	subnet = network.induced_subgraph(list(network.vs.select(name_in=pval)),implementation='create_from_scratch')
	subnet.vs.select(_degree=0).delete()

	nodes=subnet.vs["name"]
	n_nodes=network.vcount()
	position_nodes={k: v for v, k in enumerate(nodes)}
	seed_nodes=list(set(seeds).intersection(set(nodes)))

	nodes_kde=ego_friends(subnet,position_nodes,seed_nodes,sim,zscores_global)
	write_results(nodes_kde,seed_nodes,rwr,direction)


if __name__ == '__main__':
	#parameters
	sim_type=sys.argv[1]
	damping=float(sys.argv[2])
	test=sys.argv[3]

	ini_pos=sys.argv[4].split(",")
	ini_neg=sys.argv[5].split(",")
	res_folder=sys.argv[6]
	data_folder=sys.argv[7]
	#load network
	network = Graph.Read_Ncol("/nfs/research/petsalaki/shared_folder/diffusion/networks/old_net/square/"+sim_type+".txt", weights=True, directed=False)
	graph_nodes=network.vs["name"]
	number_of_nodes=network.vcount()

	seeds_pos,seeds_neg,seeds,zscores_global,ssim=load_seeds(test,graph_nodes,sim_type)
	empirical_values,pvalues=gen_val(network)
"""
	pvalues_pos=[]
	rwr_pos={}

	pvalues_neg=[]
	rwr_neg={}

	for i in pvalues:
		rwr_pos[i]=np.mean(empirical_values[i][0:3])
		rwr_neg[i]=np.mean(empirical_values[i][3:])
		if np.max(pvalues[i][0:3])>950:
			pvalues_pos.append(i)
		if np.max(pvalues[i][3:])>950:
			pvalues_neg.append(i)

	pvalues_pos=pvalues_pos+seeds[0]+seeds[1]+seeds[2]
	pvalues_neg=pvalues_neg+seeds[3]+seeds[4]+seeds[5]

	#common=list(set(pvalues_pos).intersection(pvalues_neg))

	#remove_from_pos=[]
	#remove_from_neg=[]
	#for i in common:
	#	if rwr_pos[i]<rwr_neg[i]:
	#		remove_from_pos.append(i)
	#	if rwr_neg[i]<rwr_pos[i]:
	#		remove_from_neg.append(i)

	#pvalues_pos=list(set(pvalues_pos).difference(remove_from_pos))
	#pvalues_neg=list(set(pvalues_neg).difference(remove_from_neg))

	pvalues_pos=list(set(pvalues_pos).intersection(graph_nodes))
	pvalues_neg=list(set(pvalues_neg).intersection(graph_nodes))

	ego_filtering(network,rwr_pos,pvalues_pos,seeds_pos,ssim,zscores_global,"upregulated")
	ego_filtering(network,rwr_neg,pvalues_neg,seeds_neg,ssim,zscores_global,"downregulated")
"""
