from Bio.PDB import PDBParser
import subprocess
import os
import re
import point_cloud_utils as pcu
import numpy as np
import pickle as pkl
import pandas as pd
import matplotlib as mpl

class aminoacid_score:
    def __init__(self, pdb_file, aminoacid_score, cmap = 'viridis'):
        self.aminoacid_score = aminoacid_score
        self.pdb_file = pdb_file
        self.aminoacid_afinity = self._aminoacid_afinity()
        self.n_aminoacid_afinity = (self.aminoacid_afinity - np.min(self.aminoacid_afinity)) / (np.max(self.aminoacid_afinity) - np.min(self.aminoacid_afinity))
        self.cmap = mpl.colormaps.get_cmap(cmap)
        self.colors = [self.cmap(score) for score in self.n_aminoacid_afinity]
        self.hex_colors = ["#{:02x}{:02x}{:02x}".format(int(r*255), int(g*255), int(b*255)) for r, g, b, _ in self.colors]
        self.pdb_v = self._read_pdb()
        self.clusters = self._cluster_aminoacids()
    def _aminoacid_afinity(self):
        aminoacid_afinity = []
        with open(self.aminoacid_score, 'r') as aaf:
            for line in aaf:
                part = line.split()
                try:
                    float_value = float(part[2])
                    aminoacid_afinity.append(float_value)
                except ValueError:
                    continue
        return np.array(aminoacid_afinity)

    def _read_pdb(self):
        with open(self.pdb_file, 'r') as f:
            pdb = f.read()
            return pdb
    
    def filter_clusters(clusters):
        clusters.sort(key=lambda x: x[2], reverse=True)

        filtered_clusters = []
        used_indices = set()

        for cluster in clusters:
            indices, _, _ = cluster
            if not any(index in used_indices for index in indices):
                filtered_clusters.append(cluster)
                used_indices.update(indices)
        return filtered_clusters

    def _cluster_aminoacids(self):
        score = pd.read_csv(self.aminoacid_score, delimiter= '\t', index_col=0)
        sorted_scores = score.sort_values(by = 'Score', ascending = False)
        clusters = sorted_scores.iloc[:5]
        aminoacid_groups = []
        scored_function = []
        for i in clusters.index:
            aminoacid_groups.append([i-1,i,i+1])
        for i in aminoacid_groups:
            cluster = score.iloc[i]
            scored_function.append((cluster.index.tolist(), cluster['AA'].tolist(), cluster['Score'].mean()))

        filtered_clusters = filter_clusters(scored_function)
        return filtered_clusters
    

    

class mesher(aminoacid_score):
    def __init__(self, pdb_file,faces = 15000):
        self.msms = '/content/MeshCraft-RNA/msms_i86_64Linux2_2.6.1/msms.x86_64Linux2.2.6.1'
        self.pdb_file = pdb_file
        self.name = self._pdb_name()
        self.faces = str(faces)
        self.output_files = os.path.join(os.getcwd(), 'rna_structures')
        self.pqr_name = self.name + '.pqr'
        self.pqr_file = self._pdb2pqr()
        self.xyzrn_name = self.name + '.xyzr'
        self.vert_name = self.name + '.vert'
        self.face_name = self.name + '.face'
        self.xyzrn_file, self.model_charges = self._pqr2xyzrn()
        self.object_name = self.name + '.obj'
        self._xyzrn2mesh()
        self.vert_file = os.path.join(self.output_files, self.vert_name)
        self.face_file = os.path.join(self.output_files, self.face_name)
        self.object_file = self._toobj()
        self.watertight_object = self._watertight()

    def _pdb_name(self):
        parser = PDBParser()
        structure = parser.get_structure('protein', self.pdb_file)
        if structure.header['name'] != '':
            name = structure.header['name']
            return name
        elif structure.header['idcode'] != '':
            name = structure.header['idcode'] 
            return name
        return 'unknown_protein'
    
    def _pdb2pqr(self):
        try:
            subprocess.run(["pdb2pqr", "-ff=AMBER", f"{self.pdb_file}", f"{os.path.join(self.output_files,self.pqr_name)}"], capture_output= True)
            os.remove(f'{os.path.join(self.output_files, self.name)}.log')
            return os.path.join(self.output_files,self.pqr_name)
        except:
            print('Unable to process RNA forcefields')
            return None
        
    def _pqr2xyzrn(self):
        with open(self.pqr_file) as f:
                selected_columns = []
                charges = []
                for row in f:
                    row = row.split()
                    if len(row) == 10:
                        selected_columns.append(f"  {row[5]}   {row[6]}    {row[7]} {row[9]}")
                        charges.append(float(row[8]))

        with open(os.path.join(self.output_files, self.xyzrn_name), 'w') as f:
                for row in selected_columns:
                    f.write(row + '\n')
        
        with open(os.path.join(self.output_files, self.name + 'charges.pkl'), 'wb') as f:
                pkl.dump(charges, f)
        
        return os.path.join(self.output_files, self.xyzrn_name), os.path.join(self.output_files, self.name + 'charges.pkl')
    
    def _xyzrn2mesh(self):
        try:
            subprocess.run([self.msms, '-if ', self.xyzrn_file , '-of ', os.path.join(self.output_files,self.name), '-density', '1', '-no_header'], capture_output= True)  
        except:
            print('Unable to mesh RNA')
    
    def _toobj(self):
        verts = []
        normals = []
        faces = []
        with open(self.vert_file, 'r') as f:
            for line in f:
                parts = line.split()
                verts.append([float(parts[0]), float(parts[1]), float(parts[2])])
                normals.append([float(parts[3]), float(parts[4]), float(parts[5])])
        
        with open(self.face_file, 'r') as f:
            for line in f:
                parts = line.split()
                faces.append([int(parts[0]), int(parts[1]), int(parts[2])])

        with open(os.path.join(self.output_files,self.object_name), 'w') as f:
            for row in verts:
                f.write(f'v {row[0]} {row[1]} {row[2]}\n')
            for row in normals:
                f.write(f'vn {row[0]} {row[1]} {row[2]}\n')
            for row in faces:
                f.write(f'f {row[0]}//{row[0]} {row[1]}//{row[1]} {row[2]}//{row[2]}\n')
            
        return os.path.join(self.output_files,self.object_name)
    
    def _watertight(self):
        v,f, vn = pcu.load_mesh_vfn(self.object_file)
        vw, fw = pcu.make_mesh_watertight(v,f, 20000)
        pcu.save_mesh_vf(os.path.join(self.output_files,self.name+'wt.obj'),vw, fw,dtype=np.float32)
        return os.path.join(self.output_files,self.name+'wt.obj')
    