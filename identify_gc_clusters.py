"""
Script to identify distinct GC-content subpopulations in mRNA data
and add cluster assignments to mrna_metrics.tsv
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import pandas as pd

# Read the existing TSV file
input_file = "results/mrna_metrics.tsv"
df = pd.read_csv(input_file, sep="\t")

print(f"Loaded {len(df)} mRNA sequences")
print(f"GC content range: {df['gc_content'].min():.4f} to {df['gc_content'].max():.4f}")
print(f"GC content mean: {df['gc_content'].mean():.4f}")
print(f"GC content median: {df['gc_content'].median():.4f}")
print()

# Prepare data for clustering
gc_values = df['gc_content'].values.reshape(-1, 1)

# Test different numbers of clusters (2 to 6)
silhouette_scores = []
inertias = []
k_range = range(2, 7)

for k in k_range:
    kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
    labels = kmeans.fit_predict(gc_values)
    silhouette_scores.append(silhouette_score(gc_values, labels))
    inertias.append(kmeans.inertia_)

# Find optimal k using silhouette score
optimal_k = k_range[np.argmax(silhouette_scores)]
print(f"Silhouette scores for different k values:")
for k, score in zip(k_range, silhouette_scores):
    marker = " <- OPTIMAL" if k == optimal_k else ""
    print(f"  k={k}: {score:.4f}{marker}")
print()

# Fit final model with optimal k
kmeans_final = KMeans(n_clusters=optimal_k, random_state=42, n_init=10)
df['gc_cluster'] = kmeans_final.fit_predict(gc_values)

# Get cluster centers and sort clusters by GC content
cluster_centers = kmeans_final.cluster_centers_.flatten()
sorted_cluster_indices = np.argsort(cluster_centers)

# Remap cluster labels so they are ordered by GC content (0=lowest, etc.)
cluster_mapping = {old_idx: new_idx for new_idx, old_idx in enumerate(sorted_cluster_indices)}
df['gc_cluster'] = df['gc_cluster'].map(cluster_mapping)

# Recalculate centers with new labels
cluster_stats = df.groupby('gc_cluster')['gc_content'].agg(['count', 'mean', 'std', 'min', 'max'])
print(f"Identified {optimal_k} distinct GC-content subpopulations:\n")
print(cluster_stats)
print()

# Save updated file
output_file = "results/mrna_metrics.tsv"
df.to_csv(output_file, sep="\t", index=False, float_format="%.4f")
print(f"Updated file saved to: {output_file}")

# Create visualization
fig, axes = plt.subplots(2, 1, figsize=(10, 8))

# Plot 1: Histogram with cluster overlay
ax1 = axes[0]
ax1.hist(df['gc_content'], bins=50, alpha=0.5, color='gray', edgecolor='black', density=True)
colors = plt.cm.Set2(np.linspace(0, 1, optimal_k))
for cluster_id in range(optimal_k):
    cluster_data = df[df['gc_cluster'] == cluster_id]['gc_content']
    ax1.hist(cluster_data, bins=50, alpha=0.6, label=f'Cluster {cluster_id}', 
             color=colors[cluster_id], density=True)
ax1.set_xlabel('GC Content', fontsize=12)
ax1.set_ylabel('Density', fontsize=12)
ax1.set_title('GC Content Distribution with Clusters', fontsize=14, fontweight='bold')
ax1.legend()
ax1.grid(alpha=0.3)

# Plot 2: Cluster assignment plot
ax2 = axes[1]
for cluster_id in range(optimal_k):
    cluster_data = df[df['gc_cluster'] == cluster_id]['gc_content']
    cluster_mean = cluster_stats.loc[cluster_id, 'mean']
    ax2.scatter(cluster_data, [cluster_id]*len(cluster_data), 
                alpha=0.4, s=20, color=colors[cluster_id], label=f'Cluster {cluster_id}')
    ax2.axvline(cluster_mean, color=colors[cluster_id], linestyle='--', linewidth=2, alpha=0.8)

ax2.set_xlabel('GC Content', fontsize=12)
ax2.set_ylabel('Cluster', fontsize=12)
ax2.set_title('Cluster Assignments', fontsize=14, fontweight='bold')
ax2.set_yticks(range(optimal_k))
ax2.grid(alpha=0.3)
ax2.legend(loc='upper right')

plt.tight_layout()
plt.savefig('results/gc_clusters.png', dpi=200, bbox_inches='tight')
print("Cluster visualization saved to: results/gc_clusters.png")
plt.close()

# Plot elbow and silhouette plots
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax1 = axes[0]
ax1.plot(k_range, inertias, 'bo-', linewidth=2, markersize=8)
ax1.axvline(optimal_k, color='red', linestyle='--', alpha=0.7, label=f'Optimal k={optimal_k}')
ax1.set_xlabel('Number of Clusters (k)', fontsize=12)
ax1.set_ylabel('Inertia', fontsize=12)
ax1.set_title('Elbow Method', fontsize=14, fontweight='bold')
ax1.grid(alpha=0.3)
ax1.legend()

ax2 = axes[1]
ax2.plot(k_range, silhouette_scores, 'go-', linewidth=2, markersize=8)
ax2.axvline(optimal_k, color='red', linestyle='--', alpha=0.7, label=f'Optimal k={optimal_k}')
ax2.set_xlabel('Number of Clusters (k)', fontsize=12)
ax2.set_ylabel('Silhouette Score', fontsize=12)
ax2.set_title('Silhouette Score', fontsize=14, fontweight='bold')
ax2.grid(alpha=0.3)
ax2.legend()

plt.tight_layout()
plt.savefig('results/cluster_optimization.png', dpi=200, bbox_inches='tight')
print("Cluster optimization plots saved to: results/cluster_optimization.png")
