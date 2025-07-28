#!/usr/bin/env python3
"""
Simple STRING DB Pathway Explorer
=================================

Quick start script for exploring protein pathways with STRING DB API
"""

import requests
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

def quick_pathway_analysis(start_protein="IL10", max_levels=3, score_threshold=400):
    """
    Quick pathway analysis starting from a protein
    
    Parameters:
    -----------
    start_protein : str
        Starting protein (e.g., "IL10")
    max_levels : int
        Number of levels to explore (default: 3)
    score_threshold : int
        Minimum interaction confidence (0-1000, default: 400)
    """
    
    # STRING DB configuration
    base_url = "https://string-db.org/api"
    species_id = 9606  # Homo sapiens
    
    print(f"Starting pathway analysis for {start_protein}")
    print(f"Exploring {max_levels} levels with score threshold {score_threshold}")
    print("=" * 60)
    
    # Step 1: Get starting protein info
    def get_protein_info(protein_name):
        url = f"{base_url}/json/resolve"
        params = {
            'identifier': protein_name,
            'species': species_id,
            'limit': 1
        }
        
        response = requests.get(url, params=params)
        data = response.json()
        
        if data:
            return {
                'string_id': data[0]['stringId'],
                'name': data[0]['preferredName'],
                'description': data[0]['description']
            }
        return None
    
    # Step 2: Get interactions for a set of proteins
    def get_interactions(protein_ids, threshold=400):
        if isinstance(protein_ids, str):
            protein_ids = [protein_ids]
            
        url = f"{base_url}/json/network"
        params = {
            'identifiers': '%0d'.join(protein_ids),
            'species': species_id,
            'required_score': threshold,
            'add_white_nodes': 5
        }
        
        response = requests.get(url, params=params)
        return response.json()
    
    # Initialize tracking
    G = nx.Graph()
    protein_data = {}
    level_proteins = defaultdict(set)
    
    # Get starting protein
    start_info = get_protein_info(start_protein)
    if not start_info:
        print(f"Could not find protein: {start_protein}")
        return None
    
    current_proteins = [start_info['string_id']]
    explored = {start_info['string_id']}
    protein_data[start_info['name']] = start_info
    level_proteins[0].add(start_info['name'])
    
    print(f"Starting protein: {start_info['name']} - {start_info['description']}")
    print()
    
    # Explore each level
    for level in range(1, max_levels + 1):
        print(f"Level {level}: Exploring {len(current_proteins)} proteins...")
        
        # Get interactions
        interactions = get_interactions(current_proteins, score_threshold)
        
        if not interactions:
            print(f"No interactions found at level {level}")
            break
        
        # Process interactions
        next_level = set()
        interactions_added = 0
        
        for interaction in interactions:
            protein1 = interaction['preferredName_A']
            protein2 = interaction['preferredName_B']
            string_id1 = interaction['stringId_A']
            string_id2 = interaction['stringId_B']
            score = interaction['score']
            
            # Add to graph
            G.add_edge(protein1, protein2, score=score, level=level)
            
            # Track protein info
            if protein1 not in protein_data:
                protein_data[protein1] = {
                    'string_id': string_id1,
                    'name': protein1,
                    'level': level
                }
                level_proteins[level].add(protein1)
            
            if protein2 not in protein_data:
                protein_data[protein2] = {
                    'string_id': string_id2,
                    'name': protein2,
                    'level': level
                }
                level_proteins[level].add(protein2)
            
            # Collect for next level
            for string_id in [string_id1, string_id2]:
                if string_id not in explored:
                    next_level.add(string_id)
                    explored.add(string_id)
            
            interactions_added += 1
        
        print(f"  Added {interactions_added} interactions")
        print(f"  Found {len(next_level)} new proteins")
        
        # Limit proteins for next level (prevent explosion)
        if len(next_level) > 20:
            next_level = list(next_level)[:20]
            print(f"  Limited to 20 proteins for next level")
        
        current_proteins = list(next_level)
        
        if not current_proteins:
            print(f"No more proteins to explore after level {level}")
            break
        
        print()
    
    # Results summary
    print("Analysis Summary:")
    print(f"Total proteins: {G.number_of_nodes()}")
    print(f"Total interactions: {G.number_of_edges()}")
    print()
    
    # Level breakdown
    print("Proteins by level:")
    for level, proteins in level_proteins.items():
        if proteins:
            print(f"  Level {level}: {len(proteins)} proteins")
            if level == 0:
                print(f"    {', '.join(list(proteins))}")
            else:
                sample = list(proteins)[:5]
                print(f"    {', '.join(sample)}{' ...' if len(proteins) > 5 else ''}")
    print()
    
    # Top hub proteins
    degrees = dict(G.degree())
    top_hubs = sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:10]
    print("Top hub proteins (most connections):")
    for protein, degree in top_hubs:
        print(f"  {protein}: {degree} connections")
    print()
    
    # Visualize network
    plt.figure(figsize=(14, 10))
    
    # Node colors by level
    node_colors = []
    for node in G.nodes():
        if node in level_proteins[0]:
            node_colors.append(0)  # Starting protein
        else:
            level = protein_data.get(node, {}).get('level', 1)
            node_colors.append(level)
    
    # Node sizes by degree
    node_sizes = [degrees[node] * 50 + 100 for node in G.nodes()]
    
    # Layout
    pos = nx.spring_layout(G, k=2, iterations=50)
    
    # Draw network
    nx.draw_networkx_nodes(G, pos, 
                          node_color=node_colors, 
                          node_size=node_sizes,
                          cmap=plt.cm.viridis, 
                          alpha=0.8)
    
    # Edge thickness by score
    edge_scores = [G[u][v]['score'] for u, v in G.edges()]
    max_score = max(edge_scores) if edge_scores else 1
    edge_widths = [score / max_score * 3 + 0.5 for score in edge_scores]
    
    nx.draw_networkx_edges(G, pos, 
                          width=edge_widths, 
                          alpha=0.6, 
                          edge_color='gray')
    
    # Labels for important nodes only
    important_nodes = [node for node, degree in top_hubs[:15]]
    labels = {node: node for node in important_nodes}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8)
    
    plt.title(f"Protein Interaction Network: {start_protein}\n"
             f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}, "
             f"Levels: {max_levels}")
    
    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, 
                              norm=plt.Normalize(vmin=0, vmax=max_levels))
    sm.set_array([])
    cbar = plt.colorbar(sm)
    cbar.set_label('Discovery Level', rotation=270, labelpad=20)
    
    plt.axis('off')
    plt.tight_layout()
    plt.show()
    
    return G, protein_data, level_proteins

# Run the analysis
if __name__ == "__main__":
    # Example 1: IL10 pathway
    print("Example 1: IL10 Pathway Analysis")
    G, proteins, levels = quick_pathway_analysis("IL10", max_levels=3, score_threshold=400)
    
    print("\n" + "=" * 60)
    
    # Example 2: TNF pathway with different parameters
    print("\nExample 2: TNF Pathway Analysis (stricter threshold)")
    G2, proteins2, levels2 = quick_pathway_analysis("TNF", max_levels=4, score_threshold=600)
    
    print("\n" + "=" * 60)
    print("Analysis complete!")
    print("\nYou can customize the analysis by changing:")
    print("- start_protein: Any protein name (e.g., 'IL6', 'IFNG', 'TP53')")
    print("- max_levels: Number of interaction levels to explore")
    print("- score_threshold: Confidence threshold (0-1000)")
    print("  * 150: low confidence")
    print("  * 400: medium confidence")
    print("  * 700: high confidence")
    print("  * 900: highest confidence")