import pandas as pd
import networkx as nx
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, QED
from rdkit.ML.Descriptors import MoleculeDescriptors
import copy
import matplotlib.pyplot as plt
import multiprocessing
import os
from sklearn.decomposition import PCA


source_path='./source/valid2'
sim_path='./similarity/valid2'
para_path='./connected/para'
photo_path='./connected/photo'



# 构建网络
# G是全连接，G0是非全连接

def create_net(i,threshold):
    
    sim=pd.read_csv(sim_path+'/'+i)
    G = nx.Graph()
    edges = [(sim.iloc[j, 0], sim.iloc[j, 1], sim.iloc[j, 2]) for j in range(sim.shape[0])]
    G.add_weighted_edges_from(edges)
    G0=G.copy()
    
    ceiling=max([weight for e1, e2, weight in G.edges(data='weight')])*threshold
    
    for u, v, weight in G.edges(data='weight'):
        if weight < ceiling:
            G[u][v]['weight'] = 0
            G0.remove_edge(u, v)
    
    return(G,G0)



# 自建参数
def weighted_degree(i,G):
    w_all = sum(wt for u, v, wt in G.edges(i, data='weight'))   
    return(w_all)


# 辉瑞参数
def sorted_degree(G0):   
    sorted_dict = {node: degree for node, degree in sorted(G0.degree(), key=lambda x: x[1], reverse=True)}
    return(sorted_dict)

def pfizer(G0, sorted_dict):
    G00=G0.copy()
    pfizer_predict=[]
    
    if len(sorted_dict) > 5:
        while len(pfizer_predict) < 5:
            current_edges_num=list(sorted_dict.values())[0]
            a=list(sorted_dict.values()).count(current_edges_num)
            new=list(sorted_dict.keys())[:a]
            pfizer_predict.extend(new)
            for node in new:
                G00.remove_node(node)
            sorted_dict=sorted_degree(G00)
   
    else:
        pfizer_predict = list(sorted_dict.keys())
                                  
    return(pfizer_predict)


# 计算参数
def compute_para(i,G,G0,drug,pfizer_predict):
    
    nodes=G.nodes()
    num=len(nodes)
    
    
    clustering=nx.clustering(G, weight='weight')
    
    eccentricity={}
    for j in sorted(nx.connected_components(G0), key=len, reverse=True):
        ecc=nx.eccentricity(G0.subgraph(j), weight='weight')
        eccentricity.update(ecc)
    
    distance_dict = {(e1, e2): 1 / weight for e1, e2, weight in G0.edges(data='weight')}
    nx.set_edge_attributes(G0, distance_dict, 'distance')
    closeness_centrality=nx.closeness_centrality(G0, distance='distance')
    
    degree_centrality=nx.degree_centrality(G0)
    betweenness_centrality=nx.betweenness_centrality(G0, normalized=True, weight='weight')    
    eigenvector_centrality=nx.eigenvector_centrality(G, max_iter=1000, weight='weight')

    
    
    para = pd.DataFrame({
        'pfizer': [int(j in pfizer_predict) for j in nodes],
        'node_num': num,
        'cc_num': [len(nx.node_connected_component(G0, j)) for j in nodes],
        'clustering': [clustering[j] for j in nodes],
        'eccentricity': [eccentricity[j] for j in nodes],
        'closeness_centrality': [closeness_centrality[j] for j in nodes],
        'degree_centrality': [degree_centrality[j] for j in nodes],
        'betweenness_centrality': [betweenness_centrality[j] for j in nodes],
        'eigenvector_centrality': [eigenvector_centrality[j] for j in nodes],
        'weighted_degree': [weighted_degree(j, G) for j in nodes]
    })
        
        
    drug = pd.concat([drug, para], axis=1)
#     drug = drug[(drug['cc_num'] > 1)] 
#     降噪
    
    drug['Mol'] = [Chem.MolFromSmiles(x) for x in drug['Canonical SMILES']]
    for xx,xxx in Descriptors.descList:
        df = pd.DataFrame()
        df[xx] = drug.Mol.map(xxx)
        drug=pd.concat((drug,df),axis=1)
    
    
    df = pd.DataFrame()
    df['qed detail'] = drug['Mol'].map(QED.properties)
    
    
    alerts = [x[7] for x in df['qed detail']]
    drug = pd.concat([drug, pd.DataFrame(alerts,columns=['alerts'])],axis=1)
    
    drug = drug.loc[:, drug.columns != 'Mol']
    
    train_para=list(drug)
    train_para.remove('CAS Registry Number')
    train_para.remove('Canonical SMILES')
    train_para.remove('pfizer')
    train_para.remove('node_num')

    
    if 'key' in train_para:
        train_para.remove('key')

    
    df = pd.DataFrame()

    df = pd.concat([df, drug[train_para].apply(lambda x: (x - x.mean()) / x.std())])

    drug = pd.concat([drug, df.add_prefix('z_')], axis=1)
    
    features1 = df.iloc[:,:8].values
    features2 = df.iloc[:,8:].values
    features2 = features2[:, ~np.isnan(features2).any(axis=0)]
    
    
    pca = PCA(n_components=1)
    principal_components1 = pca.fit_transform(features1)
    principal_components2 = pca.fit_transform(features2)
    
    pca_df = pd.DataFrame(drug['CAS Registry Number'])
    pca_df['PCA features1'] = principal_components1
    pca_df['PCA features2'] = principal_components2
    
    features1_dict = dict(zip(pca_df['CAS Registry Number'], pca_df['PCA features1']))
    features2_dict = dict(zip(pca_df['CAS Registry Number'], pca_df['PCA features2']))
    
    para_name= 'para.' + i.split(".", 1)[1]
    drug.to_csv(para_path +'/'+ para_name,index=False)
    
    return(nodes,features1_dict,features2_dict)






# 画图
def draw(i,G0,nodes,hit,features1_dict,features2_dict):
    
    
    pos = {node: (features1_dict[node], features2_dict[node]) for node in nodes}
    
    fig, ax = plt.subplots()
    
    # nodes
    nx.draw_networkx_nodes(G0, pos, nodelist=nodes, node_size=3,alpha=0.2,node_color='k')
    nx.draw_networkx_nodes(G0, pos, nodelist=hit, node_size=5,node_color='r')
    
    # edges
    # w = [ G0[e[0]][e[1]]['weight'] for e in G0.edges() ]
    hit=str(hit[0])
    edges = G0.edges()
    edge_colors = ['r' if hit in edge else 'orange' for edge in edges]
    alpha_values = [0.5 if hit in edge else 0.1 for edge in edges]
    edge_draw=nx.draw_networkx_edges(G0, pos, width=0.3, alpha=alpha_values, edge_color=edge_colors)
        
    ax.spines['left'].set_color('k')
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_color('k')
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left') 
    
    # 设置坐标轴范围
    plt.xlim(min(features1_dict.values()) - 0.1,max(features1_dict.values()) + 0.1)
    plt.ylim(min(features2_dict.values()) - 0.1,max(features2_dict.values()) + 0.1)
    
    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)
    plt.xlabel('topology parameters', fontsize=8)
    plt.ylabel('molecular descriptors', fontsize=8)
    
    ax=plt.axes()
    
    plt.title(i.split(".", 2)[1]+': '+str(len(nodes)), loc='left')
    plt.axis('off')
    
    fig_name=photo_path+'/'+i.split(".", 2)[1] + '.png'
    plt.savefig(fig_name,dpi=600)
    plt.close()
    
    
 # 创建网络、计算参数、画图 
    
def compute_net(i):
    G,G0=create_net(i,threshold=0.6)
    drug=pd.read_csv(source_path+'/'+i.split('.',1)[1])
    
    if 'key' in drug.columns:
        hit=drug.loc[drug['key'] == 1, 'CAS Registry Number'].tolist()
    else:
        hit=[]
    
    pfizer_predict=pfizer(G0, sorted_dict=sorted_degree(G0))
    
    nodes,features1_dict,features2_dict=compute_para(i,G,G0,drug,pfizer_predict)
    
    draw(i,G0,nodes,hit,features1_dict,features2_dict)
    
    print(i.split('.',2)[1],'para done!')
    
    
def main():
    
    file_list = os.listdir(sim_path)
    pool = multiprocessing.Pool()
    pool.map(compute_net, file_list)
    pool.close()
    pool.join()
    