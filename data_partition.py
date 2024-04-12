import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import networkx as nx
from Levenshtein import distance as lev_distance
import os
import matplotlib.pyplot as plt
import random 
from tqdm import tqdm
import argparse
import sys

def pair2single(data,pair):
  
  # create a dictionary for sequence to dataset
  sequence_to_fold_seq1 = pd.Series(pair.dataset.values, index=pair.seq1).to_dict()
  sequence_to_fold_seq2 = pd.Series(pair.dataset.values, index=pair.seq2).to_dict()

  # initialization dataset
  data['dataset'] = None
  
  # adsign for for each sequence
  for idx, row in data.iterrows():
      sequence = row['Sequence']
      if sequence in sequence_to_fold_seq1:
        data.at[idx, 'dataset'] = sequence_to_fold_seq1[sequence]
      elif sequence in sequence_to_fold_seq2:
        data.at[idx, 'dataset'] = sequence_to_fold_seq2[sequence]
      else:
        data.at[idx, 'dataset'] = 'train'
  
  single = data.dropna(subset=['dataset'])
  return single



def plot_subgraph(subG,label,name,savedir,pos=None):

    os.makedirs(savedir,exist_ok=True)
    plt.figure(figsize=(12,12))
    

    # draw subgraph
    plt.clf()  # clean recent picture
    if pos is None:
      pos = nx.spring_layout(subG)  
    nx.draw(subG, pos, with_labels=False, node_color='green', edge_color='black',node_size=50)
    if label:
      labels = {node: node for node in subG.nodes()}
      nx.draw_networkx_labels(subG, pos, labels=labels, font_size=12)
    plt.title(f'Cluster {name}')

    plt.savefig(os.path.join(savedir,f'{name}-subgraph.png'))



def plot_graphs(pairs, condition, diff=None):
  '''
  
  just draw ac pairs for visualization
  
  '''
  
  if diff is None:
    diff = 5
  print(f'current condition:{condition}')
  print(f'current diff:{diff}')
  threshold = None
  
  if condition == 'levenstein aligned': #1
  
    threshold = 1
           
    filter_info = (np.abs(pairs['dmic'])>np.log10(diff)) & (pairs['levenstein aligned']==threshold)
    origin_pairs = pairs[filter_info].copy()
    
  if condition == 'blosum62 average': #2.58
  
    threshold = 0.9
           
    filter_info = (np.abs(pairs['dmic'])>np.log10(diff)) & (pairs['blosum62 average']>=threshold)
    origin_pairs = pairs[filter_info].copy()
    
  if condition == 'blosum62 bit score': #2.58
  
    threshold = pairs['bit score'].max()*0.9
           
    filter_info = (np.abs(pairs['dmic'])>np.log10(diff)) & (pairs['bit score']>=threshold)
    origin_pairs = pairs[filter_info].copy()
  
  if condition == 'levenstein': # 2
  
    threshold = 1
    
    filter_info = (np.abs(pairs['dmic'])>np.log10(diff)) & (pairs['lev']<=threshold)
    origin_pairs = pairs[filter_info].copy()
    
  if condition == 'sequence identity': # 0.9
  
    threshold = 0.9
     
    filter_info = (np.abs(pairs['dmic'])>np.log10(diff)) & (pairs['seq id']>=threshold)
    origin_pairs = pairs[filter_info].copy()
    
  if condition == 'tanimoto average': # 0.9
  
    threshold = 0.9
    
    filter_info = (np.abs(pairs['dmic'])>np.log10(diff)) & (pairs['tanimoto simularity']>=threshold)
    origin_pairs = pairs[filter_info].copy()
  
  print(f'there are {len(origin_pairs)}/{len(pairs)} ac pairs')
  
  G = nx.Graph()
  for index, row in origin_pairs.iterrows():
      seq1, seq2 = row['seq1'], row['seq2']
      
      G.add_edge(seq1, seq2)  # add edges
  
  savedir = f"./cluster/{condition}/diff-{diff}"
  plot_subgraph(G,label=False,name='the whole graph',savedir=savedir)
  
  # '''
  # recognize connection subgraphs
  clusters = list(nx.connected_components(G))  
  print('there are {len(clusters)} clusters')
  
  cluster_ac = {}
  for i, cluster in enumerate(clusters):
      subG = G.subgraph(cluster).copy()
      positive_edges = list(subG.edges())
      num_positive_edges = len(positive_edges)
      cluster_ac[i+1] = num_positive_edges

  cluster_ac_df = pd.DataFrame.from_dict(cluster_ac,orient='index', columns=['NumEdges']).sort_values(by='NumEdges')
  cluster_ac_df.to_csv(os.path.join(savedir,'cluster_info.csv'))
      
  for i, cluster in enumerate(clusters):
    print(f'the {i}th cluster begain')
        
    # create subgraph
    subG = G.subgraph(cluster).copy()
    plot_subgraph(subG,label=True,name=f'{i+1}',savedir=savedir)
  # '''
  
    



def normal_split(data,pairs,condition,diff=None):
  '''
  
  just put ac into testset, and the rest of the seuqnce is trainset, randomly split 8:2
  
  '''
  
  if diff is None:
    diff = 5
  print(f'current condition:{condition}')
  print(f'current diff:{diff}')
  threshold = None
  
  if condition == 'levenstein aligned': #1
  
    threshold = 1
           
    filter_info = (np.abs(pairs['dmic'])>np.log10(diff)) & (pairs['levenstein aligned']==threshold)
    origin_pairs = pairs[filter_info].copy()
    
  if condition == 'blosum62 average': #2.58
  
    threshold = 0.9
           
    filter_info = (np.abs(pairs['dmic'])>np.log10(diff)) & (pairs['blosum62 average']>=threshold)
    origin_pairs = pairs[filter_info].copy()
    
  if condition == 'blosum62 bit score': #2.58
  
    threshold = pairs['bit score'].max()*0.9
           
    filter_info = (np.abs(pairs['dmic'])>np.log10(diff)) & (pairs['bit score']>=threshold)
    origin_pairs = pairs[filter_info].copy()
  
  if condition == 'levenstein': # 2
  
    threshold = 1
    
    filter_info = (np.abs(pairs['dmic'])>np.log10(diff)) & (pairs['lev']<=threshold)
    origin_pairs = pairs[filter_info].copy()
    
  if condition == 'sequence identity': # 0.9
  
    threshold = 0.9
     
    filter_info = (np.abs(pairs['dmic'])>np.log10(diff)) & (pairs['seq id']>=threshold)
    origin_pairs = pairs[filter_info].copy()
    
  if condition == 'tanimoto average': # 0.9
  
    threshold = 0.9
    
    filter_info = (np.abs(pairs['dmic'])>np.log10(diff)) & (pairs['tanimoto simularity']>=threshold)
    origin_pairs = pairs[filter_info].copy()
  
  print(f'there are {len(origin_pairs)}/{len(pairs)} ac pairs')
  
  origin_pairs['dataset'] = 'test'
  
  dataset = pair2single(data,origin_pairs)
  train_dataset = dataset[dataset['dataset']=='train'].copy()
  test_dataset = dataset[dataset['dataset']=='test'].copy()
  
  
  num_bins = 10  
  train_dataset['label_binned'] = pd.qcut(train_dataset['Activity'], q=num_bins, labels=False, duplicates='drop')
  
  X_train, X_val, y_train, y_val = train_test_split(train_dataset, train_dataset['Activity'], test_size=0.2, stratify=train_dataset['label_binned'], random_state=42)
  
  print(f'the number of data in train:valid:test={X_train.shape[0]}:{X_val.shape[0]}:{test_dataset.shape[0]}')
  
  
  savedir = os.path.join('./data', condition,f"diff_{diff}" )
  
  os.makedirs(savedir, exist_ok=True)
  
  train_dataset.to_csv(os.path.join(savedir,org_file_name.split('/')[-1].replace('.csv','-all-train.csv')))
  X_train.to_csv(os.path.join(savedir,org_file_name.split('/')[-1].replace('.csv','-train.csv')))
  X_val.to_csv(os.path.join(savedir,org_file_name.split('/')[-1].replace('.csv','-valid.csv')))
  test_dataset.to_csv(os.path.join(savedir,org_file_name.split('/')[-1].replace('.csv','-test.csv')))

def parse_args():
    parser = argparse.ArgumentParser(description="Data Partition, generating train/valid/test set.")
    parser.add_argument('--data','-d', type=str, default='./data/grampa_s_aureus_7_25.csv', help='Path to the data CSV file')
    parser.add_argument('--cliff_pairs','-p', type=str, default='./data/grampa_s_aureus_7_25_acpairs_blast.csv', help='the generated AMP-Cliff pairs by "generate_cliffs.py"')
    
    return parser.parse_args()
      
if __name__=='__main__':
    args = parse_args()
    org_file_name = args.data
    pair_file_name = args.cliff_pairs
    pairs = pd.read_csv(pair_file_name)
    data = pd.read_csv(org_file_name)
    
    
    allowed_conditions = ['tanimoto average','blosum62 average','all']:
    
    condition = pair_file_name.split('_')[-2]
    
    if condition not in allowed_conditions:
    
      print(f"Error: The AMP-Cliff pairs file name is illegal, the file name must contain one of {allowed_conditions}, please checck the file.")
      sys.exit(1)

    if condition == 'all':
      for condition in allowed_conditions:
        for diff in [2,3,4,5]:
          normal_split(data,pairs,condition,diff)
    else:
      
      diff = int(pair_file_name.split('_')[-1].split('-')[0])
      normal_split(data,pairs,condition,diff)