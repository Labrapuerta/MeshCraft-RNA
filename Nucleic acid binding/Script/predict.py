import numpy as np
import os, random, pickle
import datetime
from tqdm import tqdm
from sklearn.metrics import auc, roc_auc_score, precision_recall_curve
from sklearn.model_selection import KFold
from sklearn import metrics
import torch
import torch.nn as nn
from torch.utils.data import RandomSampler
import torch_geometric
from torch_geometric.loader import DataLoader
from data import *
from model import *
import argparse
import numpy as np
import os


parser = argparse.ArgumentParser()
parser.add_argument("--dataset_path", type=str, default='../Example/structure_data/')
parser.add_argument("--feature_path", type=str, default='../Example/prottrans/')
parser.add_argument("--num_workers", type=int, default=8)
parser.add_argument("--ligand", type=str, default='DNA')
#parser.add_argument('--fasta_file', type=str, default='../Example/input.fasta')
parser.add_argument('--output_prottrans', type=str, default='../Example/prottrans/')
parser.add_argument('--output_esmfold', type=str, default='../Example/structure_data/')
parser.add_argument('--output_dssp', type=str, default='../Example/structure_data/')
args = parser.parse_args()

test_data = {}
with open(args.fasta_file) as r1:
    fasta_ori = r1.readlines()
for i in range(len(fasta_ori)):
    if fasta_ori[i][0] == ">":
        name = fasta_ori[i].split('>')[1].replace('\n','')
        seq = fasta_ori[i+1].replace('/n','')
        test_data[name] = seq.replace('\n','')


test_dataset = ProteinGraphDataset(test_data, range(len(test_data)), args, task_list=["DNA", "RNA"])
test_dataloader = DataLoader(test_dataset, batch_size = 1, shuffle=False, drop_last=False, num_workers=args.num_workers, prefetch_factor=2)

models = []
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model_class = GLMSite
for fold in range(5):
    state_dict = torch.load(model_path + 'fold%s.ckpt'%fold, device)
    model = model_class(1044, 128, 5, 0.95, 0.1).to(device)
    model.load_state_dict(state_dict)
    model.eval()
    models.append(model)

pred_dict = {}
for data in tqdm(test_dataloader):
    try:
        data = data.to(device)
        h_V = (data.node_s, data.node_v)
        h_E = (data.edge_s, data.edge_v)
        with torch.no_grad():
            outputs = [model(h_V, data.edge_index, h_E, data.seq).sigmoid() for model in models]
            outputs = torch.stack(outputs,0).mean(0) 
            pred_dict[data.name[0]] = outputs
    
    except Exception as e:
        print(e)
    continue

for key in pred_dict.keys():
    with open('../Example/results/{}.txt'.format(key),'w') as w:
        w.writelines('# The results of GLMSite' + '\n')
        w.writelines('Num' + "\t" + 'AA' + "\t" + 'Score' + '\n')
        for i in range(len(pred_dict[key])):
            if args.ligand == 'DNA':
                w.writelines(str(i) + "\t" + test_data[key][i] + "\t" + str(np.round(pred_dict[key][i][0].cpu().numpy(),3)) + '\n')
            else:
                w.writelines(str(i) + "\t" + test_data[key][i] + "\t" + str(np.round(pred_dict[key][i][1].cpu().numpy(),3)) + '\n')
