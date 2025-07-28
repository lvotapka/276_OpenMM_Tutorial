"""
cluster.py

Given an MD trajectory, this script will cluster the trajectory by RMSD
and write the cluster centroids as PDB files.
"""
import mdtraj as md
import numpy as np
import scipy.cluster.hierarchy
from scipy.spatial.distance import squareform

# Define important variables first
trajectory_file_name = "tryp_benz.pdb"
align_by_selection = "name CA" # alpha carbons
clustering_linkage_method = "average" # may also be 'single', 'complete', etc.
rmsd_cluster_cutoff = 3.75 # in Angstroms

# Load the trajectory file
traj = md.load(trajectory_file_name)

# Align the trajectory by alpha carbons
align_indices = traj.topology.select(align_by_selection)
traj.superpose(traj, atom_indices=align_indices, ref_atom_indices=align_indices)

# Compute all pairwise RMSDs between conformations
distances = np.empty((traj.n_frames, traj.n_frames))
for i in range(traj.n_frames):
    distances[i] = md.rmsd(traj, traj, i)
print('Max pairwise rmsd: %f nm' % np.max(distances))

# Clustering only accepts reduced form. Squareform's checks are too stringent
assert np.all(distances - distances.T < 1e-6)
reduced_distances = squareform(distances, checks=False)

# Implement the average-link clustering
linkage = scipy.cluster.hierarchy.linkage(
    reduced_distances, method=clustering_linkage_method)

# Make lists of cluster members
flattened_clusters = scipy.cluster.hierarchy.fcluster(
    linkage, rmsd_cluster_cutoff, criterion="distance")
num_clusters = np.max(flattened_clusters)
cluster_frames = []
for cluster_index in range(1, num_clusters+1):
    cluster_frame_list = []
    for i, cluster_index2 in enumerate(flattened_clusters):
        if cluster_index2 == cluster_index:
            cluster_frame_list.append(i)
    cluster_frames.append(cluster_frame_list)

# Find the cluster centroids
for i, cluster in enumerate(cluster_frames):
    cluster_distances = np.empty((len(cluster), len(cluster)))
    for j, distance_index in enumerate(cluster):
        for k, distance_index2 in enumerate(cluster):
            cluster_distances[j,k] = distances[distance_index, distance_index2]
            
    beta = 1
    index = np.exp(-beta*cluster_distances / cluster_distances.std()).sum(axis=1).argmax()
    
    print("cluster:", i, "centroid frame:", cluster[index])
    
    # Save centroid structures
    centroid_name = "tryp_ben_centroid_{}.pdb".format(i)
    traj[cluster[index]].save(centroid_name)

