/*
An Experimental Study on Hub Labeling based Shortest Path Algorithms [Experiments and Analyses]

Authors: Ye Li, Leong Hou U, Man Lung Yiu, Ngai Meng Kou
Contact: yb47438@umac.mo
Affiliation: University of Macau

The MIT License (MIT)

Copyright (c) 2016 University of Macau

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#pragma once
#ifndef COVER_ORDERING_H
#define COVER_ORDERING_H

#include<algorithm>
#include<unordered_set>
//#include<unordered_map>
#include<time.h>
#include "graph.h"
#include "graph_search.h"
#include "labels.h"
#include "time_util.h" 
#include "heap.h"
#include "paras.h"
#ifdef _WIN32
	#include<google/sparsehash/sparseconfig.h>
#endif
    #include<google/dense_hash_map>
	#include<google/dense_hash_set>

#define numOfVertices SP_Constants::numOfVertices
#define numOfEdges SP_Constants::numOfEdges
#define WEIGHTED_FLAG SP_Constants::WEIGHTED_FLAG  
#define DIRECTED_FLAG SP_Constants::DIRECTED_FLAG
#define CNT 16
#define cover_value_type long long
#include<unordered_map>  
#include<cmath>

using namespace time_util;
  
class COrdering {
public:
	vector<NodeID> inv; // Fetch the original vertex id by a given ranking.
	vector<NodeID> rank; // Fetch the ranking of a given vertex id.
	
	void Relabel(Graph& graph) {
		
		for (NodeID v = 0; v < numOfVertices; ++v) rank[inv[v]] = v;
		
		// Array Representation
		vector<EdgeID> new_vertices(numOfVertices + 1);
		vector<NodeID> new_edges;
		new_edges.reserve(graph.edges.size());
		for (NodeID ranking = 0; ranking < numOfVertices; ++ranking) {
			NodeID originalVertex = inv[ranking];
			for (EdgeID eid = graph.vertices[originalVertex]; eid < graph.vertices[originalVertex + 1]; ++eid)
				new_edges.push_back(rank[graph.edges[eid]]);
			new_vertices[ranking + 1] = new_edges.size();
		}
		graph.vertices.swap(new_vertices);
		graph.edges.swap(new_edges);

		if (DIRECTED_FLAG == true) {
			vector<EdgeID> r_new_vertices(numOfVertices + 1);
			vector<NodeID> r_new_edges;
			r_new_edges.reserve(graph.r_edges.size());
			for (NodeID ranking = 0; ranking < numOfVertices; ++ranking) {
				NodeID originalVertex = inv[ranking];
				for (EdgeID eid = graph.r_vertices[originalVertex]; eid < graph.r_vertices[originalVertex + 1]; ++eid)
					r_new_edges.push_back(rank[graph.r_edges[eid]]);
				r_new_vertices[ranking + 1] = r_new_edges.size();
			}
			graph.r_vertices.swap(r_new_vertices);
			graph.r_edges.swap(r_new_edges);
		}

	}

	void Relabel(WGraph& wgraph) {

		for (NodeID v = 0; v < numOfVertices; ++v) rank[inv[v]] = v;

		/*vector<vector<NodeID> > new_adj(numOfVertices);
		vector<vector<EdgeWeight> > new_adj_weight(numOfVertices);
		vector<vector<NodeID> > new_r_adj(numOfVertices);
		vector<vector<EdgeWeight> > new_r_adj_weight(numOfVertices);
		for (NodeID v = 0; v < numOfVertices; ++v) {
			for (NodeID i = 0; i < wgraph.adj[v].size(); ++i) {
				new_adj[rank[v]].push_back(rank[wgraph.adj[v][i]]);
				new_adj_weight[rank[v]].push_back(wgraph.adj_weight[v][i]);
			}
			if (DIRECTED_FLAG == true) {
				for (NodeID i = 0; i < wgraph.r_adj[v].size(); ++i) {
					new_r_adj[rank[v]].push_back(rank[wgraph.r_adj[v][i]]);
					new_r_adj_weight[rank[v]].push_back(wgraph.r_adj_weight[v][i]);
				}
			}
		}
		wgraph.adj.swap(new_adj);
		wgraph.adj_weight.swap(new_adj_weight);
		if (DIRECTED_FLAG == true) {
			wgraph.r_adj.swap(new_r_adj);
			wgraph.r_adj_weight.swap(new_r_adj_weight);
		}*/

		// Array Representation
		vector<EdgeID> new_vertices(numOfVertices + 1);
		vector<NodeEdgeWeightPair> new_edges;
		new_edges.reserve(wgraph.edges.size());
		for (NodeID ranking = 0; ranking < numOfVertices; ++ranking) {
			NodeID originalVertex = inv[ranking];
			for (EdgeID eid = wgraph.vertices[originalVertex]; eid < wgraph.vertices[originalVertex + 1]; ++eid)
				new_edges.push_back(make_pair(rank[wgraph.edges[eid].first],wgraph.edges[eid].second));
			new_vertices[ranking + 1] = new_edges.size();
		}
		wgraph.vertices.swap(new_vertices);
		wgraph.edges.swap(new_edges);

		if (DIRECTED_FLAG == true) {
			vector<EdgeID> r_new_vertices(numOfVertices + 1);
			vector<NodeEdgeWeightPair> r_new_edges;
			r_new_edges.reserve(wgraph.r_edges.size());
			for (NodeID ranking = 0; ranking < numOfVertices; ++ranking) {
				NodeID originalVertex = inv[ranking];
				for (EdgeID eid = wgraph.r_vertices[originalVertex]; eid < wgraph.r_vertices[originalVertex + 1]; ++eid)
					r_new_edges.push_back(make_pair(rank[wgraph.r_edges[eid].first],wgraph.r_edges[eid].second));
				r_new_vertices[ranking + 1] = r_new_edges.size();
			}
			wgraph.r_vertices.swap(r_new_vertices);
			wgraph.r_edges.swap(r_new_edges);
		}
	} 

	void ReswapLabel(Graph& graph) {

		vector<vector<NodeID> > new_adj(numOfVertices);
		for (NodeID v = 0; v < numOfVertices; ++v) {
			for (NodeID i = 0; i < graph.adj[v].size(); ++i) {
				new_adj[v].push_back(inv[graph.adj[rank[v]][i]]);
			}
		}
		graph.adj.swap(new_adj);
	}

	void save_rank(const char* order_file) {
		ofstream ofs(order_file);
		for (int i = 0; i < numOfVertices; ++i) {
			ofs << inv[i] << endl;
		}
		ofs.close();
	}

	~COrdering() {
		inv.clear();
		rank.clear();
	}

};



class Coverage_Ordering : public COrdering {
	typedef	vector<NodeID> tree;
	public:
		Label labels;
		DLabel dlabels;
		NodeID last_available;
		
		//tree parent_tree: while a vertex v is not in the tree, parent_tree[v] = numOfVertices
	
	NodeID labeling_source_bfs(NodeID source, tree& parent_tree, vector<NodeID>& coverage, vector<NodeID>& descendants, vector<NodeID>& root_hop, vector<NodeID>& last_hop, vector<NodeID>& que, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<bool>& usd, NodeID ranking, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, Graph& graph) { 
			descendants.clear();
			NodeID visited_arcs = 0;
			
			int max_hop = -1;
			
			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];
					
			for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}
			
			NodeID que_t0 = 0, que_t1 = 0, que_h = 0;

			que[que_h++] = source;
			vis[source] = true;
			que_t1 = que_h;
			parent_tree[source] = numOfVertices;
			root_hop[source] = 0;
			last_hop[source] = 0;

			for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					
					
					NodeID v = que[que_i];
					if (usd[v]) continue; 

					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
					for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
						NodeID w = tmp_idx_v.first[i];
						EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
						if (td <= d) {
							goto pruned_forward;
						}
					}
 
				//	if(source == 28674)
				//		cout << v << " ," << tmp_idx_v.first.size()  << "," << descendants.size() << endl;

					tmp_idx_v.first.back() = ranking;
					tmp_idx_v.second.back() = d;

					tmp_idx_v.first.push_back(numOfVertices);
					tmp_idx_v.second.push_back(INF_WEIGHT);

					descendants.push_back(v);
					visited_arcs++;
					 
				//	if(source == 28674)
				//		cout << "2" << endl;
					for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid){
	
						NodeID w = graph.edges[eid];
						if (!vis[w]) {											
							que[que_h++] = w;
							parent_tree[w] = v;
							root_hop[w] = root_hop[v] + 1;
							last_hop[w] = root_hop[w];
							if(max_hop < root_hop[w])
								max_hop = root_hop[w];
							vis[w] = true;
						}
					}
				//	if(source == 28674)
				//		cout << "3" << endl;
				pruned_forward:
					{}
				}
				que_t0 = que_t1;
				que_t1 = que_h;
			}
			for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
			for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
			}
			
			usd[source] = true;
			return visited_arcs;
		}
	
	NodeID labeling_source_dij(NodeID source, tree& parent_tree, vector<NodeID>& coverage, vector<NodeID>& descendants, vector<NodeID>& root_hop, vector<NodeID>& last_hop, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, queue<NodeID>& visited_que, vector<EdgeWeight>& distances, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<bool>& usd, NodeID ranking, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, WGraph& wgraph){
			
		descendants.clear();
		NodeID visited_arcs = 0;
		
		const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];

		for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
		}

		pqueue.update(source, 0);
		distances[source] = 0;
		parent_tree[source] = numOfVertices;
		root_hop[source] = 0;
		
		while (!pqueue.empty()) {

			NodeID v;
			EdgeWeight v_d;
			pqueue.extract_min(v, v_d);
			pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
			vis[v] = true;
			visited_que.push(v);

			if (usd[v]) continue;
			for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
				NodeID w = tmp_idx_v.first[i];
				EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
				if (td <= v_d) {
					goto pruned;
				}
			}

			// Traverse
			tmp_idx_v.first.back() = ranking;
			tmp_idx_v.second.back() = v_d;
			tmp_idx_v.first.push_back(numOfVertices);
			tmp_idx_v.second.push_back(INF_WEIGHT);
			descendants.push_back(v);
			visited_arcs++;

			for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid) {
				NodeID w = wgraph.edges[eid].first;
				EdgeWeight w_d = wgraph.edges[eid].second + v_d;
				if (!vis[w]) {
					if( distances[w] > w_d ){
						pqueue.update(w, w_d);
						distances[w] = w_d;
						parent_tree[w] = v;
						root_hop[w] = root_hop[v] + 1;
					}
				}
			}
			pruned: 
				{}
			}

			while (!visited_que.empty()) {
				NodeID vis_v = visited_que.front();
				visited_que.pop();
				vis[vis_v] = false;
				distances[vis_v] = INF_WEIGHT;
				pqueue.clear(vis_v);
			}

			pqueue.clear_n();

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
			usd[source] = true;
			
			return visited_arcs;
	}
	
	NodeID labeling_source_bfs_directed(NodeID source, tree& parent_tree, tree& r_parent_tree, vector<NodeID>& coverage, vector<NodeID>& r_coverage, vector<NodeID>& descendants, vector<NodeID>& r_descendants,vector<NodeID>& root_hop, vector<NodeID>& r_root_hop, vector<NodeID>& que, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<EdgeWeight>& r_dst_r, vector<bool>& usd, NodeID ranking, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& r_tmp_idx,Graph& graph) { 
			descendants.clear();
			r_descendants.clear();
			NodeID visited_arcs = 0;
			
			// Forward search.
			// Initialize forward labels of r.
			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}

			NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
			que[que_h++] = source;
			vis[source] = true;
			que_t1 = que_h;
			parent_tree[source] = numOfVertices;
			root_hop[source] = 0;

			for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID v = que[que_i];
					pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];
					//index_t &idx_v = index_[inv[v]];

					if (usd[v]) continue;

					// Pruned by the forward labels of r and backward labels of v in the forward search from r when reaching v.
					for (size_t i = 0; i < r_tmp_idx_v.first.size(); ++i) {
						NodeID w = r_tmp_idx_v.first[i];
						EdgeWeight td = r_tmp_idx_v.second[i] + dst_r[w];
						if (td <= d) {
							goto pruned_forward;
						}
					}

					// Traverse
					r_tmp_idx_v.first.back() = ranking;
					r_tmp_idx_v.second.back() = d;
					r_tmp_idx_v.first.push_back(numOfVertices);
					r_tmp_idx_v.second.push_back(INF_WEIGHT);
					descendants.push_back(v);
					
					for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid){
						NodeID w = graph.edges[eid];
						if (!vis[w]) {
							que[que_h++] = w;
							vis[w] = true;
							parent_tree[w] = v;
							root_hop[w] = root_hop[v] + 1;
						}
					}
				pruned_forward:
					{}
				}
				que_t0 = que_t1;
				que_t1 = que_h;
			}
			for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;


			// Backward search.
			// Initialize backward labels of r.
			const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[source];

			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i) {
				r_dst_r[r_tmp_idx_r.first[i]] = r_tmp_idx_r.second[i];
			}


			que_t0 = 0, que_t1 = 0, que_h = 0;
			que[que_h++] = source;
			vis[source] = true;
			que_t1 = que_h;
			r_parent_tree[source] = numOfVertices;
			r_root_hop[source] = 0;

			for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID v = que[que_i];
					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
					//index_t &idx_v = index_[inv[v]];

					if (usd[v]) continue;

					// Pruned by the backward labels of r and forward labels of v in the backward search from r when reaching v (v->r path).
					for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
						NodeID w = tmp_idx_v.first[i];
						EdgeWeight td = tmp_idx_v.second[i] + r_dst_r[w];
						if (td <= d) {
							goto pruned_backward;
						}
					}

					// Traverse
					tmp_idx_v.first.back() = ranking;
					tmp_idx_v.second.back() = d;
					tmp_idx_v.first.push_back(numOfVertices);
					tmp_idx_v.second.push_back(INF_WEIGHT);
					r_descendants.push_back(v);

					// Array Representation
					for (EdgeID eid = graph.r_vertices[v]; eid < graph.r_vertices[v + 1]; ++eid) {
						NodeID w = graph.r_edges[eid];

						if (!vis[w]) {
							que[que_h++] = w;
							vis[w] = true;
							r_parent_tree[w] = v;
							r_root_hop[w] = r_root_hop[v] + 1;
						}
					}
				pruned_backward:
					{}
				}
				que_t0 = que_t1;
				que_t1 = que_h;
			}
			for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i)
				r_dst_r[r_tmp_idx_r.first[i]] = INF_WEIGHT;
 
			usd[source] = true;
			return visited_arcs;
		}
	
	NodeID labeling_source_dij_directed(NodeID source, tree& parent_tree, tree& r_parent_tree, vector<NodeID>& coverage, vector<NodeID>& r_coverage, vector<NodeID>& descendants, vector<NodeID>& r_descendants, vector<NodeID>& root_hop, vector<NodeID>& r_root_hop, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, queue<NodeID>& visited_que, vector<EdgeWeight>& distances, vector<bool>& vis, vector<EdgeWeight>& dst_r,vector<EdgeWeight>& r_dst_r, vector<bool>& usd, NodeID ranking, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& r_tmp_idx,WGraph& wgraph){
			
		descendants.clear();
		r_descendants.clear();
		NodeID visited_arcs = 0;
		
		const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];

		for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
		}

		pqueue.update(source, 0);
		distances[source] = 0;		
		parent_tree[source] = numOfVertices;
		root_hop[source] = 0;
		
		while (!pqueue.empty()) {

			NodeID v;
			EdgeWeight v_d;
			pqueue.extract_min(v, v_d);
			pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];
			vis[v] = true;
			visited_que.push(v);
 
			if (usd[v]) continue;
			for (size_t i = 0; i < r_tmp_idx_v.first.size(); ++i) {
				NodeID w = r_tmp_idx_v.first[i];
				EdgeWeight td = r_tmp_idx_v.second[i] + dst_r[w];
				if (td <= v_d) {
					goto pruned_forward;
				}
			}

			// Traverse
			r_tmp_idx_v.first.back() = ranking;
			r_tmp_idx_v.second.back() = v_d;
			r_tmp_idx_v.first.push_back(numOfVertices);
			r_tmp_idx_v.second.push_back(INF_WEIGHT);
			descendants.push_back(v);
			visited_arcs++;

			for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid) {
				NodeID w = wgraph.edges[eid].first;
				EdgeWeight w_d = wgraph.edges[eid].second + v_d;
				if (!vis[w]) {
					if( distances[w] > w_d ){
						pqueue.update(w, w_d);
						distances[w] = w_d;
						parent_tree[w] = v;
						root_hop[w] = root_hop[v] + 1;
					}
				}
			}
			pruned_forward: 
				{}
			}

			while (!visited_que.empty()) {
				NodeID vis_v = visited_que.front();
				visited_que.pop();
				vis[vis_v] = false;
				distances[vis_v] = INF_WEIGHT;
				pqueue.clear(vis_v); 
			}

			pqueue.clear_n();

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
			
			
			const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[source];

			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i) {
				dst_r[r_tmp_idx_r.first[i]] = r_tmp_idx_r.second[i];
			}

			pqueue.update(source, 0);
			distances[source] = 0;		
			
			
			// reverse search
			while (!pqueue.empty()) {

				NodeID v;
				EdgeWeight v_d;
				pqueue.extract_min(v, v_d);
				pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
				vis[v] = true;
				visited_que.push(v);

				if (usd[v]) continue;
				for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
					NodeID w = tmp_idx_v.first[i];
					EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
					if (td <= v_d) {
						goto pruned_backward;
					}
				}

				// Traverse
				tmp_idx_v.first.back() = ranking;
				tmp_idx_v.second.back() = v_d;
				tmp_idx_v.first.push_back(numOfVertices);
				tmp_idx_v.second.push_back(INF_WEIGHT);
				r_descendants.push_back(v);
				visited_arcs++;

				for (EdgeID eid = wgraph.r_vertices[v]; eid < wgraph.r_vertices[v + 1]; ++eid) {
					NodeID w = wgraph.r_edges[eid].first;
					EdgeWeight w_d = wgraph.r_edges[eid].second + v_d;
					if (!vis[w]) {
						if( distances[w] > w_d ){
							pqueue.update(w, w_d);
							distances[w] = w_d;
							r_parent_tree[w] = v;
							r_root_hop[w] = root_hop[v] + 1;
						}
					}
				}
				pruned_backward: 
					{}
			}

			while (!visited_que.empty()) {
				NodeID vis_v = visited_que.front();
				visited_que.pop();
				vis[vis_v] = false;
				distances[vis_v] = INF_WEIGHT;
				pqueue.clear(vis_v);
			}

			pqueue.clear_n();

			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i)
				dst_r[r_tmp_idx_r.first[i]] = INF_WEIGHT;
			
			usd[source] = true;
			
			return visited_arcs;
	}
	
	
	NodeID bfs_walk(NodeID source, tree& parent_tree, vector<NodeID>& coverage,  vector<NodeID>& root_hop, vector<NodeID>& last_hop, vector<NodeID>& que, vector<NodeID>& que2,vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<bool>& usd, Graph& graph) { 
			
			NodeID visited_arcs = 0;
			
			NodeID que_t0 = 0, que_t1 = 0, que_h = 0, que_h2 = 0;

			que2[que_h2++] = source;
			vis[source] = true;
			 
			for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID v = que2[que_i];
					
					for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid){
						NodeID w = graph.edges[eid];
						if ( parent_tree[w] == v){
							que2[que_h++] = w;
							vis[w] == true;
						}
					}	
				}					
				que_t0 = que_t1;
				que_t1 = que_h;
			}

			que_t0 = 0, que_t1 = 0, que_h = 0;
			que_t1 = que_h;
			que[que_h++] = source;

			for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID v = que[que_i];
					if (usd[v]) continue; 

					EdgeWeight td = root_hop[source] + root_hop[v];
					if (td <= d) {
						goto pruned_forward;
					}
					
					visited_arcs++;
					 
					for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid){

						NodeID w = graph.edges[eid];

						if (!vis[w]) {
							que[que_h++] = w;
							vis[w] = true;
						}
					}
				pruned_forward:
					{}
				} 
				que_t0 = que_t1;
				que_t1 = que_h;
			}

			for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
			for (size_t i = 0; i < que_h2; ++i) vis[que2[i]] = false;
			
			return visited_arcs;
		}
	
	
	NodeID bfs_one_hop_walk(NodeID source, tree& parent_tree, vector<NodeID>& coverage,  vector<NodeID>& root_hop, vector<NodeID>& last_hop, vector<NodeID>& que, vector<NodeID>& que2,vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<bool>& usd, Graph& graph) { 
			
			NodeID visited_arcs = 0;
			
			NodeID que_t0 = 0, que_t1 = 0, que_h = 0, que_h2 = 0;

			que2[que_h2++] = source;
			vis[source] = true;
			 
			for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID v = que2[que_i];
					
					for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid){
						NodeID w = graph.edges[eid];
						if ( parent_tree[w] == v){
							que2[que_h++] = w;
							vis[w] == true;
						}
					}	
				}					
				que_t0 = que_t1;
				que_t1 = que_h;
			}

			que_t0 = 0, que_t1 = 0, que_h = 0;
			que_t1 = que_h;
			que[que_h++] = source;

			for (EdgeWeight d = 0; que_t0 < que_h && d < 4; d = d + 1) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID v = que[que_i];
					if (usd[v]) continue; 

					EdgeWeight td = root_hop[source] + root_hop[v];
					if (td <= d) {
						goto pruned_forward;
					}
					
					visited_arcs++;
					 
					for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid){

						NodeID w = graph.edges[eid];

						if (!vis[w]) {
							que[que_h++] = w;
							vis[w] = true;
						}
					}
				pruned_forward:
					{}
				}
				que_t0 = que_t1;
				que_t1 = que_h;
			}

			for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
			for (size_t i = 0; i < que_h2; ++i) vis[que2[i]] = false;
			
			return visited_arcs;
		}
		
	void calcover(vector<NodeID>& descendants, tree& parent_tree, vector<NodeID>& coverage, vector<NodeID>& root_hop, vector<NodeID>& last_alive, vector<NodeID>& depth){
		for (NodeID di = descendants.size() - 1; di > -1; --di) {
			NodeID dv = descendants[di];
			// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
			coverage[dv]++;
			//acc[dv]++;
			if(parent_tree[dv]!= numOfVertices){
				coverage[parent_tree[dv]] += coverage[dv];	
				//acc[parent_tree[dv]]++;	
				if(depth[dv] < root_hop[dv])
					depth[dv] = root_hop[dv];
				if(depth[dv] > root_hop[parent_tree[dv]]){
					depth[parent_tree[dv]] = depth[dv];
				}
				last_alive[dv] = coverage[dv];
			}
		}
		return;
	}
	
	void calcover(vector<NodeID>& descendants, tree& parent_tree, vector<NodeID>& coverage, vector<NodeID>& root_hop){
		for (NodeID di = descendants.size() - 1; di > -1; --di) {
			NodeID dv = descendants[di];
			coverage[dv]++;
			if(parent_tree[dv]!= numOfVertices){
				coverage[parent_tree[dv]] += coverage[dv];	
			}
		}
		return;
	}
	
	
	void clear_tmp(vector<NodeID>& descendants, vector<NodeID>& coverage, tree& parent_tree, vector<NodeID>& root_hop, vector<NodeID>& depth){
		for (NodeID di = descendants.size() - 1; di > -1; --di) {
			NodeID dv = descendants[di];
			// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
			coverage[dv] = 0;
			if(parent_tree[dv]!= numOfVertices)
				coverage[parent_tree[dv]] = 0;	
			parent_tree[dv] = numOfVertices;
			root_hop[dv] = 0;
			depth[dv] = 0;
		}
		descendants.clear();
		return;
	}
	
	NodeID findComputeAncestor(NodeID v, vector<NodeID>& coverage, tree& parent_tree, vector<NodeID>& root_hop){
		NodeID anv = parent_tree[v];
		
		for(NodeID hop = 1; hop < root_hop[v] / 2; ++hop){
			anv = parent_tree[anv];
		}						
		return anv;
	}

	double normalize(NodeID mina, NodeID maxa, NodeID minb, NodeID maxb, NodeID a, NodeID b){
		double a_score = (double) a * (double)mina / (double) maxa;
		double b_score = (double) b * (double)minb / (double) maxb;
		return a_score + b_score;
	}
	
	vector<NodeID> findSigPath(NodeID source, vector<NodeID>& coverage, tree& parent_tree, vector<bool>& usd, vector<NodeID>& upwardsPow, Graph& graph){
		vector<NodeID> sigpath;		
		upwardsPow.clear();
		NodeID v = source;
		while(true){			
			NodeID max_v = v;
			NodeID max_cover = -1;			
			bool non_leaf_flag = false;
			NodeID childrenNeighbors = 0;
			for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid){
					NodeID w = graph.edges[eid];
					if (parent_tree[w] == v && usd[w] == false) {
						non_leaf_flag = true;
						childrenNeighbors++;
						if(max_cover < coverage[w]){
							max_cover = coverage[w];
							max_v = w;
						}/*
						if(max_cover < graph.vertices[w + 1] - graph.vertices[w]){
							max_cover = graph.vertices[w + 1] - graph.vertices[w];
							max_v = w;
						}*/
					}					
			}
			/*if( v != source){
				if(childrenNeighbors == 0)
					upwardsPow.push_back(0);
				else 
					upwardsPow.push_back( (coverage[v] / childrenNeighbors ) * ( graph.vertices[v + 1 ] - graph.vertices[v] - childrenNeighbors)) ;
			}*/
			if(non_leaf_flag == false) break;
			v = max_v;
			sigpath.push_back(v);
		}		
		//cout << "sig path len:" << sigpath.size() << endl;
		return sigpath;
	}
	
	vector<NodeID> findRevSigPath(NodeID source, vector<NodeID>& coverage, vector<NodeID>& r_coverage, tree& parent_tree,tree& r_parent_tree, vector<bool>& usd, Graph& graph){
		vector<NodeID> sigpath;		
		NodeID v = source;
		NodeID max_degree = -1;
		while(true){			
			NodeID max_v = v;
			NodeID max_cover = -1;
			bool non_leaf_flag = false;
			if(max_degree < (graph.vertices[v + 1] - graph.vertices[v]) * (graph.r_vertices[v + 1] - graph.r_vertices[v]))
				max_degree = (graph.vertices[v + 1] - graph.vertices[v]) * (graph.r_vertices[v + 1] - graph.r_vertices[v]);
			for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid){
					NodeID w = graph.edges[eid];
					if (parent_tree[w] == v && usd[w] == false) {
						non_leaf_flag = true;
						if(max_cover < coverage[w] + r_coverage[w]){
							max_cover = coverage[w] + r_coverage[w];
							max_v = w;
						}
					}					
			}
			if(non_leaf_flag == false) break;
			v = max_v;
			sigpath.push_back(v);
		}
		
		vector<NodeID> rev_sigpath;		
		v = source;
		NodeID max_degree_rev = -1;
		while(true){			
			NodeID max_v = v;
			NodeID max_cover = -1;
			bool non_leaf_flag = false;
			if(max_degree_rev < (graph.r_vertices[v + 1] - graph.r_vertices[v]) * (graph.vertices[v + 1] - graph.vertices[v]))
				max_degree_rev = (graph.r_vertices[v + 1] - graph.r_vertices[v]) * (graph.vertices[v + 1] - graph.vertices[v]);
			for (EdgeID eid = graph.r_vertices[v]; eid < graph.r_vertices[v + 1]; ++eid){
					NodeID w = graph.r_edges[eid];
					if (r_parent_tree[w] == v && usd[w] == false) {
						non_leaf_flag = true;
						if(max_cover < r_coverage[w] + coverage[w]){
							max_cover =  r_coverage[w] + coverage[w];
							max_v = w;
						}
					}					
			} 
			if(non_leaf_flag == false) break;
			v = max_v;
			rev_sigpath.push_back(v);
		}
		/*
		cout << "sig path len:" << sigpath.size() << endl;
		NodeID maxd = -1;
		NodeID maxcd = -1;
		NodeID maxd_choice = -1;
		NodeID maxcd_choice = -1;
		if(sigpath.size()!=0){
			for(NodeID i = 0; i < sigpath.size(); ++i){
			//	cout << i << ":" << (graph.vertices[sigpath[i]+1] - graph.vertices[sigpath[i]]) *(graph.r_vertices[sigpath[i]+1] - graph.r_vertices[sigpath[i]]) << " ";
				if(maxd<(graph.vertices[sigpath[i]+1] - graph.vertices[sigpath[i]]) *(graph.r_vertices[sigpath[i]+1] - graph.r_vertices[sigpath[i]])){
					maxd_choice = i;
					maxd = (graph.vertices[sigpath[i]+1] - graph.vertices[sigpath[i]]) *(graph.r_vertices[sigpath[i]+1] - graph.r_vertices[sigpath[i]]);
				}
			}
			cout << endl;				
			for(NodeID i = 0; i < sigpath.size() - 1; ++i){
				//cout << i << ":" << (graph.vertices[sigpath[i]+1] - graph.vertices[sigpath[i]]) *(graph.r_vertices[sigpath[i]+1] - graph.r_vertices[sigpath[i]]) * abs(coverage[sigpath[i+1]] - coverage[sigpath[i]]) << " ";
				if(maxcd < (graph.vertices[sigpath[i]+1] - graph.vertices[sigpath[i]]) *(graph.r_vertices[sigpath[i]+1] - graph.r_vertices[sigpath[i]]) * abs(coverage[sigpath[i+1]] - coverage[sigpath[i]])){
					maxcd = (graph.vertices[sigpath[i]+1] - graph.vertices[sigpath[i]]) *(graph.r_vertices[sigpath[i]+1] - graph.r_vertices[sigpath[i]]) * abs(coverage[sigpath[i+1]] - coverage[sigpath[i]]);
					maxcd_choice = i;
				}
			}
		//	cout << endl;
			//cout << maxd_choice << " vs. " << maxcd_choice << endl;

			//for(NodeID i = 0; i < sigpath.size() - 1; ++i)
			//	cout << i << ":" << (graph.vertices[sigpath[i+1]+1] - graph.vertices[sigpath[i+1]]) *(graph.r_vertices[sigpath[i+1]+1] - graph.r_vertices[sigpath[i+1]]) * abs(coverage[sigpath[i+1]] - coverage[sigpath[i]]) << " ";
			//cout << endl;
		}
		
		cout << "rev sig path len:" << rev_sigpath.size() << endl;
		maxd_choice = -1;
		maxcd_choice = -1;
		NodeID rev_maxd = -1;
		NodeID rev_maxcd = -1;
		if(rev_sigpath.size()!=0){			
			for(NodeID i = 0; i < rev_sigpath.size(); ++i){
				//cout << i << ":" << (graph.vertices[rev_sigpath[i]+1] - graph.vertices[rev_sigpath[i]]) * (graph.r_vertices[rev_sigpath[i]+1] - graph.r_vertices[rev_sigpath[i]]) << " ";
				if(rev_maxd <  (graph.vertices[rev_sigpath[i]+1] - graph.vertices[rev_sigpath[i]]) * (graph.r_vertices[rev_sigpath[i]+1] - graph.r_vertices[rev_sigpath[i]])){
					rev_maxd =  (graph.vertices[rev_sigpath[i]+1] - graph.vertices[rev_sigpath[i]]) * (graph.r_vertices[rev_sigpath[i]+1] - graph.r_vertices[rev_sigpath[i]]);
					maxd_choice = i;
				}
			}
			cout << endl;			
			//for(NodeID i = 0; i < rev_sigpath.size() - 1; ++i)
			//	cout << i << ":" << abs(r_coverage[rev_sigpath[i+1]] - r_coverage[rev_sigpath[i]]) << " ";
			//cout << endl;		
			for(NodeID i = 0; i < rev_sigpath.size() - 1; ++i){
				//cout << i << ":" << (graph.vertices[rev_sigpath[i]+1] - graph.vertices[rev_sigpath[i]]) * (graph.r_vertices[rev_sigpath[i]+1] - graph.r_vertices[rev_sigpath[i]]) * abs(r_coverage[rev_sigpath[i+1]] - r_coverage[rev_sigpath[i]]) << " ";
				if(rev_maxcd <(graph.vertices[rev_sigpath[i]+1] - graph.vertices[rev_sigpath[i]]) * (graph.r_vertices[rev_sigpath[i]+1] - graph.r_vertices[rev_sigpath[i]]) * abs(r_coverage[rev_sigpath[i+1]] - r_coverage[rev_sigpath[i]]) ){
					rev_maxcd = (graph.vertices[rev_sigpath[i]+1] - graph.vertices[rev_sigpath[i]]) * (graph.r_vertices[rev_sigpath[i]+1] - graph.r_vertices[rev_sigpath[i]]) * abs(r_coverage[rev_sigpath[i+1]] - r_coverage[rev_sigpath[i]]);
					maxcd_choice = i;
				}
			}
		//	cout << endl;
		//	cout << maxd_choice << " vs. " << maxcd_choice << endl;
		}
		*/
		
		
		if(max_degree < max_degree_rev){
			reverse(rev_sigpath.begin(), rev_sigpath.end());
			sigpath = rev_sigpath;
		} 
		
		return sigpath;
	}
	
	vector<NodeID> findRevSigPath(NodeID source, vector<NodeID>& coverage, vector<NodeID>& r_coverage, tree& parent_tree,tree& r_parent_tree, vector<bool>& usd, WGraph& wgraph, NodeID& maxdegree, bool& reverse_flag){
		vector<NodeID> sigpath;		
		NodeID v = source;
		NodeID max_degree = -1;
		while(true){			
			NodeID max_v = v;
			NodeID max_cover = -1;
			bool non_leaf_flag = false;
			if(max_degree < (wgraph.vertices[v + 1] - wgraph.vertices[v]) * (wgraph.r_vertices[v + 1] - wgraph.r_vertices[v]))
				max_degree = (wgraph.vertices[v + 1] - wgraph.vertices[v]) * (wgraph.r_vertices[v + 1] - wgraph.r_vertices[v]);
			for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid){
					NodeID w = wgraph.edges[eid].first;
					if (parent_tree[w] == v && usd[w] == false) {
						non_leaf_flag = true;
						NodeID c_coverage = coverage[w] * (wgraph.r_vertices[v + 1] - wgraph.r_vertices[v]);// + r_coverage[w];//(wgraph.r_vertices[v + 1] - wgraph.r_vertices[v]); //coverage[w] + r_coverage[w]						
						if(max_cover < c_coverage){
							max_cover = c_coverage;
							max_v = w;
						}
					}					
			}
			if(non_leaf_flag == false) break;
			v = max_v;
			sigpath.push_back(v);
		}
		
		vector<NodeID> rev_sigpath;		 
		v = source;
		NodeID max_degree_rev = -1;
		while(true){			
			NodeID max_v = v;
			NodeID max_cover = -1;
			bool non_leaf_flag = false;
			if(max_degree_rev < (wgraph.r_vertices[v + 1] - wgraph.r_vertices[v]) * (wgraph.vertices[v + 1] - wgraph.vertices[v]))
				max_degree_rev = (wgraph.r_vertices[v + 1] - wgraph.r_vertices[v]) * (wgraph.vertices[v + 1] - wgraph.vertices[v]);
			for (EdgeID eid = wgraph.r_vertices[v]; eid < wgraph.r_vertices[v + 1]; ++eid){
					NodeID w = wgraph.r_edges[eid].first;
					if (r_parent_tree[w] == v && usd[w] == false) {
						non_leaf_flag = true;
						NodeID c_coverage = r_coverage[w]* (wgraph.vertices[v + 1] - wgraph.vertices[v]);// + coverage[w];// (wgraph.vertices[v + 1] - wgraph.vertices[v]); //r_coverage[w] + coverage[w];
						if(max_cover < c_coverage){
							max_cover =  c_coverage;
							max_v = w;
						}
					}					
			} 
			if(non_leaf_flag == false) break;
			v = max_v;
			rev_sigpath.push_back(v);
		}
		
		maxdegree = max_degree;
		reverse_flag = false;
		
		if(max_degree < max_degree_rev){
			reverse(rev_sigpath.begin(), rev_sigpath.end());
			sigpath = rev_sigpath;
			maxdegree = max_degree_rev;
			reverse_flag = true;
		} /*
		if(coverage[0] < r_coverage[0]){
			reverse_flag = true;
		}*/
		
		return sigpath;
	}
	
	
	vector<NodeID> findSigPath(NodeID source, vector<NodeID>& coverage, tree& parent_tree, vector<bool>& usd,vector<NodeID>& upwardsPow, WGraph& wgraph, NodeID& maxdegree){
		vector<NodeID> sigpath;		
		NodeID v = source;
		upwardsPow.clear();
		maxdegree = -1;
		while(true){			
			NodeID max_v = v;
			NodeID max_cover = -1;
			NodeID max_v_degree;
			bool non_leaf_flag = false;
			for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid) {
				NodeID w = wgraph.edges[eid].first;
				if (parent_tree[w] == v && usd[w] == false) {
					non_leaf_flag = true;
					if(max_cover < coverage[w]){
						max_cover = coverage[w];
						max_v = w;
						max_v_degree = wgraph.vertices[w+1] - wgraph.vertices[w];
					}
				/*	if(max_cover < abs(coverage[v] - coverage[w]) * (wgraph.vertices[w + 1] - wgraph.vertices[w]) ){
						max_cover = abs(coverage[v] - coverage[w]) * (wgraph.vertices[w + 1] - wgraph.vertices[w]);
						max_v = w;
					}*/
				}					
			}
			if(non_leaf_flag == false) break;
			v = max_v;
			if(maxdegree < max_v_degree)
				maxdegree = max_v_degree;
			sigpath.push_back(v);
		}		
		//cout << "sig path len:" << sigpath.size() << endl;
		return sigpath;
	}
	
	/*
	vector<NodeID> findSigPath(NodeID source, vector<NodeID>& coverage, tree& parent_tree, vector<bool>& usd,vector<NodeID>& upwardsPow, WGraph& wgraph){
		vector<NodeID> sigpath;		
		NodeID v = source;
		upwardsPow.clear();
		while(true){			
			NodeID max_v = v;
			NodeID max_cover = -1;
			bool non_leaf_flag = false;
			NodeID childrenNeighbors = 0;			
			for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid) {
				NodeID w = wgraph.edges[eid].first;
				if (parent_tree[w] == v && usd[w] == false) {
					non_leaf_flag = true;
					childrenNeighbors++;
					if(max_cover < coverage[w]){
						max_cover = coverage[w];
						max_v = w;
					}
				}					
			}
			if( v != source){
				if(childrenNeighbors == 0)
					upwardsPow.push_back(0);
				else
					upwardsPow.push_back( coverage[v] / childrenNeighbors * ( wgraph.vertices[v + 1 ] - wgraph.vertices[v] - childrenNeighbors)) ;
			}
			if(non_leaf_flag == false) break;
			v = max_v;
			sigpath.push_back(v);
		}		
		//cout << "sig path len:" << sigpath.size() << endl;
		return sigpath;
	}
	*/
	
	NodeID calupwardFanout(vector<NodeID>& descendants, vector<NodeID>& upwardcoverage, vector<NodeID>& coverage, tree& parent_tree, vector<NodeID>& root_hop, Graph& graph){
		for(NodeID i = 0; i < descendants.size(); ++i){
			NodeID v = descendants[i];
			if(root_hop[v] == 0) continue;
			for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid){
				NodeID w = graph.edges[eid];
				if (root_hop[v] == 1 && root_hop[w] != 0 ) {
					upwardcoverage[v]++;
				}

				if(root_hop[v] > 1){
					if(root_hop[w] <= root_hop[v])
						upwardcoverage[v]++;
				}
				
				/*if (root_hop[v] > root_hop[w] && root_hop[w] != 0)
					upwardcoverage[w] += upwardcoverage[v];	
				*/
			}				
		}
	}
	
	void clear_upwardtmp(vector<NodeID>& descendants, vector<NodeID>& upwardcoverage){
		for(NodeID i = 0; i < descendants.size(); ++i){
			upwardcoverage[descendants[i]] = 0;
		}
	}
	
	NodeID calupwardTree(NodeID v, tree& parent_tree, vector<NodeID>& root_hop, Graph& graph){
		NodeID calup = 0;
		for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid){
			NodeID w = graph.edges[eid];
			if (parent_tree[w] != v) {
				calup += root_hop[w] + 1;
			}
		}
		return calup;
	}
	
	void candidate_gen(vector<NodeID>& candidate,vector<NodeID>& inv, vector<bool>& usd){
		NodeID candidate_size = 300;
		candidate.clear();
		for(NodeID i = last_available; i < numOfVertices; ++i){
			NodeID v = inv[i];
			// Candidate
			if(usd[v] == false)
				candidate.push_back(v);
			if(candidate.size() > candidate_size - 1)
				break;
		}
		
		while(usd[inv[last_available]] == true)
			last_available++;
		
		return;
	}
	
	NodeID topcan(vector<NodeID>& inv, vector<bool>& usd, NodeID& last_available){
		while(true){
			NodeID v = inv[last_available];
			if(usd[v] == false)
				return v;
			else
				last_available++;
		}
	}
	
	/*NodeID pick(Graph& graph, vector<bool>& usd, NodeID iter_num, vector<NodeID>& candidate, vector<NodeID>& coverage, vector<NodeID>& acc_count, tree& parent_tree, vector<NodeID>& root_hop, vector<NodeID>& last_hop, vector<NodeID>& last_alive, vector<NodeID>& depth){
		
		NodeID chosen = numOfVertices;
		
		priority_queue<pair<double, NodeID> > combine;
		
		double alpha = 1;
		double beta = 0;
		
		NodeID maxhop = -1;
		NodeID minhop = numOfVertices;
		for(NodeID i = 0; i < candidate.size(); ++i){
			NodeID v = candidate[i];
			if(root_hop[v] != 0 ){
				if(maxhop < root_hop[v])
					maxhop = root_hop[v];
				if(minhop > root_hop[v])
					minhop = root_hop[v];
			}else{
				if(maxhop < last_hop[v])
					maxhop = last_hop[v];
				if(minhop > last_hop[v])
					minhop = last_hop[v];
			}
		}
		
		srand(time(NULL));
		
		bool sum_flag = false;
		int roll;
		double pvalue;
			
		for(NodeID i = 0; i < candidate.size(); ++i){
			NodeID v = candidate[i];
			roll = rand()%(maxhop-minhop)+minhop;
			
			if(root_hop[v] != 0 ){
				if(roll <= root_hop[v])
					sum_flag = true;
				else
					sum_flag = false;
			}else{
				if(roll <= last_hop[v])
					sum_flag = true;
				else
					sum_flag = false;
			}
			
			if(sum_flag == true)
				pvalue = sum_neighbors(v, graph, usd,  coverage, acc_count, parent_tree, root_hop);
			else 
				pvalue = coverage[v];
			
			combine.push(make_pair(pvalue, v));
		}
		
		
		for(NodeID i = 0; i < candidate.size(); ++i){
			NodeID v = candidate[i];
			
			combine.push(make_pair(coverage[v] * root_hop[v], v));
		}
		chosen = combine.top().second;
		return chosen;
		
	}*/
/*
	NodeID pick(Graph& graph, vector<bool>& usd, NodeID iter_num, vector<NodeID>& candidate, vector<NodeID>& coverage, vector<NodeID>& second_coverage){
		
		NodeID chosen = numOfVertices;
		
		priority_queue<pair<double, NodeID> > combine;
		
		for(NodeID i = 0; i < candidate.size(); ++i){
			NodeID v = candidate[i];			
			combine.push(make_pair(coverage[v] + second_coverage[v], v));
		}
		chosen = combine.top().second;
		return chosen;
		
	}
	*/
	
	NodeID pick(Graph& graph, vector<bool>& usd, NodeID iter_num,  vector<NodeID>& candidate, vector<NodeID>& coverage, vector<NodeID>& acc_count, tree& parent_tree, vector<NodeID>& root_hop, vector<NodeID>& last_hop, vector<NodeID>& last_alive, vector<NodeID>& depth){
		
		NodeID chosen = numOfVertices;
		
		priority_queue<pair<double, NodeID> > combine;
		
		for(NodeID i = 0; i < candidate.size(); ++i){
			NodeID v = candidate[i];	
			NodeID pv = findComputeAncestor(v, coverage, parent_tree, root_hop);
			NodeID cover_v = sum_neighbors(v, pv, graph, usd, coverage, acc_count, parent_tree, root_hop, last_alive, last_hop);
			
			/*
			NodeID cover_v = coverage[pv];
			if(root_hop[pv] == 0)
				cover_v = coverage[v];
			*/
			combine.push(make_pair(cover_v, v));
		}
		chosen = combine.top().second;
		coverage[chosen] = combine.top().first;
		return chosen;
		
	}
	
	NodeID simplepick(Graph& graph, vector<bool>& usd, NodeID iter_num, vector<NodeID>& candidate, vector<NodeID>& coverage, vector<NodeID>& root_hop){
		
		NodeID chosen = numOfVertices;
		
		priority_queue<pair<double, NodeID> > combine;
		
		for(NodeID i = 0; i < candidate.size(); ++i){
			NodeID v = candidate[i];			
			combine.push(make_pair(coverage[v] * root_hop[v], v));
		}
		chosen = combine.top().second;
		return chosen;
		
	}
	
	NodeID sigpathpick(Graph& graph, vector<bool>& usd, NodeID iter_num, vector<NodeID>& candidate, vector<NodeID>& coverage, vector<NodeID>& root_hop){
		
		NodeID chosen = numOfVertices;
		
		priority_queue<pair<double, NodeID> > combine;
		
		for(NodeID i = 0; i < candidate.size(); ++i){
			NodeID v = candidate[i];			
			combine.push(make_pair(coverage[v] * root_hop[v], v));
		}
		chosen = combine.top().second;
		return chosen;
		
	}
	
	/*
	NodeID sum_neighbors(NodeID v, Graph& graph, vector<bool>& usd,  vector<NodeID>& coverage, vector<NodeID>& acc_count, tree& parent_tree, vector<NodeID>& root_hop, vector<NodeID>& last_alive){
		priority_queue<pair<NodeID, NodeID> > process_order;
		unordered_set<NodeID> processed;
		processed.insert(v);
		acc_count[v] = coverage[v];
		
		for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid){			
			NodeID w = graph.edges[eid];
			if(usd[w]) continue;
			// process the ones in higher positions in the tree 
			process_order.push(make_pair(-root_hop[w], w));
		}
		
		while(process_order.empty() != true){
			NodeID w = process_order.top().second;
			NodeID hop = process_order.top().first;
			process_order.pop();
			if(hop != 0){
				if(processed.find(parent_tree[w]) != processed.end()) continue;
				processed.insert(w);
				acc_count[v] += coverage[w];
			} else
				acc_count[v] += last_alive[w];
		}
		processed.clear();
		NodeID pvalue = acc_count[v];
		acc_count[v] = 0;
		return pvalue;
	}
	*/
	NodeID sum_neighbors(NodeID v, NodeID pv, Graph& graph, vector<bool>& usd,  vector<NodeID>& coverage, vector<NodeID>& acc_count, tree& parent_tree, vector<NodeID>& root_hop, vector<NodeID>& last_alive, vector<NodeID>& last_hop){
		priority_queue<pair<NodeID, NodeID> > process_order;
		unordered_set<NodeID> processed;
		processed.insert(v);
		if(root_hop[pv] == 0)
			acc_count[v] = coverage[v];
		else
			acc_count[v] = coverage[pv];
		
		for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid){			
			NodeID w = graph.edges[eid];
			if(usd[w]) continue;
			process_order.push(make_pair(-root_hop[w], w));
		}
		
		while(process_order.empty() != true){
			NodeID w = process_order.top().second;
			NodeID hop = process_order.top().first;
			process_order.pop();
			if(hop != 0){
				if(processed.find(parent_tree[w]) != processed.end()){
					processed.insert(w);
					continue;
				}
				processed.insert(w);
				if(w == pv) continue;
				acc_count[v] += coverage[w] * root_hop[w];
			}
			else
				acc_count[v] += last_alive[w] * last_hop[w];
		}
		processed.clear();
		NodeID pvalue = acc_count[v];
		acc_count[v] = 0;
		return pvalue;
	}
	
	int build_tmp_tree(NodeID source, tree& parent_tree, vector<NodeID>& coverage, vector<NodeID>& descendants, vector<NodeID>& root_hop, vector<NodeID>& last_hop, vector<NodeID>& que, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<bool>& usd, NodeID ranking, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, Graph& graph) { 
			descendants.clear();
			
			NodeID visited_arcs = 0;
			
			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];
					
			for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}
			
			NodeID que_t0 = 0, que_t1 = 0, que_h = 0;

			que[que_h++] = source;
			vis[source] = true;
			que_t1 = que_h;
			parent_tree[source] = numOfVertices;
			root_hop[source] = 0;
			last_hop[source] = 0;

			for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID v = que[que_i];
					if (usd[v]) continue; 

					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
					for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
						NodeID w = tmp_idx_v.first[i];
						EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
						if (td <= d) {
							goto pruned_forward;
						}
					}

					descendants.push_back(v);
					 
					for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid){

						NodeID w = graph.edges[eid];

						if (!vis[w]) {
							que[que_h++] = w;
							parent_tree[w] = v;
							root_hop[w] = root_hop[v] + 1;
							last_hop[w] = root_hop[w];
							vis[w] = true;
						}
					}
				pruned_forward:
					{}
				}
				que_t0 = que_t1;
				que_t1 = que_h;
			}

			for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
			for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
			}
			
			return 0;
		}
	
	void calcover_secondtree(vector<NodeID>& descendants, tree& parent_tree, tree& first_parent_tree,vector<NodeID>& coverage, vector<NodeID>& root_hop, vector<NodeID>& last_alive, vector<NodeID>& depth){
		for (NodeID di = descendants.size() - 1; di > -1; --di) {
			NodeID dv = descendants[di];
			// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
			coverage[dv]++;
			//acc[dv]++;
			if(parent_tree[dv]!= numOfVertices){
				if(first_parent_tree[dv] != parent_tree[dv])
					coverage[parent_tree[dv]] += coverage[dv];	
				//acc[parent_tree[dv]]++;	
				if(depth[dv] < root_hop[dv])
					depth[dv] = root_hop[dv];
				if(depth[dv] > root_hop[parent_tree[dv]]){
					depth[parent_tree[dv]] = depth[dv];
				}
				last_alive[dv] = coverage[dv];
			}
		}
		return;
	}
	
	void walk_stats(Graph& graph, vector<NodeID> border){
		cout << "Building Coverage Ordering Based Labels" << endl;
		inv.resize(numOfVertices);
		rank.resize(numOfVertices);
		
		vector<NodeID> deg_inv;
		vector<NodeID> deg_rank;
		deg_inv.resize(numOfVertices);
		deg_rank.resize(numOfVertices);
				
		vector<pair<float, NodeID> > deg(numOfVertices);
		for (size_t v = 0; v < numOfVertices; ++v) {
				deg[v] = make_pair((graph.vertices[v + 1] - graph.vertices[v]) + float(rand()) / RAND_MAX, v);
		}
		sort(deg.rbegin(), deg.rend());
		for (size_t v = 0; v < numOfVertices; ++v) deg_inv[v] = deg[v].second;
		
		last_available = 0;
		
		vector<NodeID> parent_tree(numOfVertices, numOfVertices);
		vector<NodeID> root_hop(numOfVertices, 0);	
		vector<NodeID> coverage(numOfVertices, 0);
		vector<NodeID> depth(numOfVertices, 0);
		vector<NodeID> last_alive(numOfVertices, 0);		
		vector<NodeID> last_hop(numOfVertices, 0);
		vector<NodeID> descendants;
		descendants.reserve(numOfVertices);
		
		vector<NodeID> second_parent_tree(numOfVertices, numOfVertices);
		vector<NodeID> second_root_hop(numOfVertices, 0);	
		vector<NodeID> second_coverage(numOfVertices, 0);
		vector<NodeID> second_depth(numOfVertices, 0);
		vector<NodeID> second_last_alive(numOfVertices, 0);		
		vector<NodeID> second_last_hop(numOfVertices, 0);
		vector<NodeID> second_descendants;
		second_descendants.reserve(numOfVertices);
				
		vector<NodeID> max_hop_set;
		
		vector<NodeID> candidate;
		candidate.reserve(numOfVertices);
				
		vector<NodeID> acc_count;
		acc_count.resize(numOfVertices);
		
		vector<NodeID> que(numOfVertices);
		vector<NodeID> que2(numOfVertices);
		vector<bool> vis(numOfVertices);
		vector<bool> usd(numOfVertices, false);
// Preparing basic structure for pl algorithm.
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));
				 
		NodeID chosen = deg_inv[0];
		
		NodeID taken = numOfVertices; // /10;
		NodeID last_iter_cover = 0;
		NodeID last_iter_hop = 0;
		unordered_map<int, double> actual_stats;
		unordered_map<int, double> over_stats;
		
		vector<NodeID> upwardcoverage(numOfVertices, 0);
		
		long long currentsum = 0;
		
		for(NodeID i = 0; i < taken; ++i){
			inv[i] = chosen; 
			rank[chosen] = i;
			if(usd[chosen] == true) cout << "shit" << endl;
			NodeID actual_cover = labeling_source_bfs(chosen, parent_tree, coverage, descendants, root_hop, last_hop, que, vis, dst_r, usd, i, tmp_idx, graph);
			currentsum += actual_cover;
			/*
			if(last_iter_hop > 0){
				actual_stats[last_iter_hop] += actual_cover;
				over_stats[last_iter_hop] += last_iter_cover;
			}
			cout << last_iter_cover << " vs. " << actual_cover << " on " << last_iter_hop << " hops" << endl;
			*/
			if(i == numOfVertices - 1 ) break;  
			
			/*for(NodeID j = descendants.size() - 1; j >= 0; --j){
				if(c_max_hop == root_hop[descendants[j]])
					if(usd[descendants[j]] == false)
						max_hop_set.push_back(descendants[j]);
			}*/
			
		//	NodeID second_source = max_hop_set[rand()%max_hop_set.size()];
			
		//	build_tmp_tree(second_source, second_parent_tree, second_coverage, second_descendants, second_root_hop, second_last_hop, que, vis, dst_r, usd, i, tmp_idx, graph);
			calcover(descendants, parent_tree, coverage, root_hop, last_alive, depth);
			
			
			
		//	calcover_secondtree(second_descendants, second_parent_tree, parent_tree, second_coverage, second_root_hop, second_last_alive, second_depth);
			 
		//	NodeID topcan = candidate_gen(candidate, deg_inv, usd); 
			
	/*		vector<NodeID> level_howmany(50, 0);
			vector<NodeID> level_counts(50, 0);
		
			for(NodeID j = descendants.size() - 1; j >= 0; --j){
				NodeID tv = descendants[j];
				level_howmany[root_hop[tv]]++;
			}
		*/	
			int iterround = 30;
			//for(int j = 0; j < candidate.size(); ++j){
			//	NodeID tv = candidate[j];
			/*
			calupwardFanout(descendants, upwardcoverage, coverage, parent_tree,root_hop, graph);
			
			for(int j = descendants.size() - 1; j >= 0; --j){
				NodeID tv = descendants[j];
				if(root_hop[tv] == 0 ) continue;
				if(level_counts[root_hop[tv]] >= iterround) continue;
				int rollnum = rand()%level_howmany[root_hop[tv]];
				if(rollnum < iterround){
					//NodeID walk = bfs_walk(tv,  parent_tree, coverage,  root_hop, last_hop, que,que2,vis,  dst_r, usd,  graph);
					//NodeID walk = bfs_one_hop_walk(tv,  parent_tree, coverage,  root_hop, last_hop, que,que2,vis,  dst_r, usd,  graph);
					//NodeID walk = calupwardTree(tv, parent_tree, root_hop, graph);	
					NodeID walk = upwardcoverage[tv];
					level_counts[root_hop[tv]]++; 
					cout <<  tv << " " << root_hop[tv] << " " << coverage[tv] << " " << walk << " " << graph.vertices[tv + 1] - graph.vertices[tv] << " ";
				}
			}  
			*/
			vector<NodeID> upwardsPow;
			vector<NodeID> sigpath = findSigPath(chosen, coverage, parent_tree, usd, upwardsPow, graph);
		
		
		if(sigpath.size()==0){
			chosen = topcan(deg_inv, usd, last_available);
		}
			
		NodeID maxDegree = -1;
		NodeID chosenhop = -1;
	
		for(NodeID j = 0; j < sigpath.size(); ++j){
			NodeID curDegree = graph.vertices[sigpath[j] +1] - graph.vertices[sigpath[j]];
			if(curDegree > maxDegree){
				maxDegree = curDegree;
				chosen = sigpath[j];
				chosenhop = j;
			}
			//cout << bfs_walk(sigpath[j], parent_tree, coverage,  root_hop, last_hop, que,que2,vis, dst_r, usd, graph) + coverage[sigpath[j]] << "," << curDegree << "," << coverage[sigpath[j]] << " " ;
		} 
		//cout << endl;
		
		/*
		if(sigpath.size()>2){
			NodeID mintrend = INF_WEIGHT;
			NodeID maxtrend = -1;
			NodeID mindegree = INF_WEIGHT;
			NodeID maxdegree = -1;
			NodeID chosentrend = 0;
			//cout << "Coverage Trends:";
			for(NodeID j = 0; j < sigpath.size() - 1; ++j){
				NodeID curtrend = abs(coverage[sigpath[j+1]] - coverage[sigpath[j]]);
				NodeID curdegree = (graph.vertices[sigpath[j] +1] - graph.vertices[sigpath[j]]);
				if(maxtrend <  curtrend){
					 maxtrend = curtrend;
					//chosentrend = j;
					//chosen = sigpath[j];
					//chosenhop = j;
				}
				if(mintrend > curtrend){
					mintrend = curtrend;
				}
				if(maxdegree < curdegree) maxdegree = curdegree;
				if(mindegree > curdegree) mindegree = curdegree;
			} 
			
			double max_score = -1;
			//cout << "Coverage Trends:";
			for(NodeID j = 0; j < sigpath.size() - 1; ++j){
				NodeID curtrend = abs(coverage[sigpath[j+1]] - coverage[sigpath[j]]);
				NodeID curdegree = (graph.vertices[sigpath[j] +1] - graph.vertices[sigpath[j]]);
				double curscore =  normalize(mindegree, maxdegree, mintrend, maxtrend, curdegree, curtrend);
				if(max_score < curscore){
					curscore = max_score;
					chosen = sigpath[j];
					chosenhop = j;
				}
			} 
		}*/
		
	/*	if(sigpath.size()<=1){
			chosen = topcan(deg_inv, usd);
		}
					
	
		if(sigpath.size()>1){
			NodeID maxtrend = -1;
			NodeID chosentrend = 0;
			cout << "Coverage Trends:";
			for(NodeID j = 1; j < sigpath.size() - 1; ++j){
				//NodeID curDegree = wgraph.vertices[sigpath[j] +1] - wgraph.vertices[sigpath[j]];
				//cout << " " << curDegree;
				//NodeID curDegree = wgraph.vertices[sigpath[j] +1] - wgraph.vertices[sigpath[j]];
				cout << " " <<  j << ":" << abs(coverage[sigpath[j+1]] - coverage[sigpath[j]]);
				if(maxtrend <  abs(coverage[sigpath[j+1]] - coverage[sigpath[j]])){
					 maxtrend = abs(coverage[sigpath[j+1]] - coverage[sigpath[j]]);
					chosentrend = j;
				}
			}
		}
		cout << endl;
		*/
	/*	for(NodeID j = 0; j < numOfVertices; ++j){
			if(usd[deg_inv[j]] == false ){
				if(root_hop[deg_inv[j]] == chosenhop){
					chosen = deg_inv[j];
				}
			}
		}*/
		
	//		cout << "iteration:" << i << " - " << currentsum << " - " << usd[chosen] << "-" << chosenhop << " - " << sigpath.size() << endl;
			 
	/*	if(usd[chosen] == true){
			cout << inv[i] << ", shit," << chosen << endl;
			for(NodeID j = 1; j < sigpath.size(); ++j){
				cout << " " << sigpath[j];
			}
			cout << endl;
		} 
		*/
			
	/*		
			for(NodeID j = 1; j < hop_counts.size(); ++j){
				if(hop_counts[j] != 0 )
					if(j < sigpath.size())
						cout <<  sigpath[j] << " " << j << " " << 0 << " " << hop_walks[j] / hop_counts[j] << " " << graph.vertices[sigpath[j] + 1] - graph.vertices[sigpath[j]] << " ";
					else
						cout <<  j << " " << j << " " << 0 << " " << hop_walks[j] / hop_counts[j] << " " << j << " ";

			}
			
			cout << endl; 
			*/
			//clear_upwardtmp(descendants, upwardcoverage);
			
			/*cout << "candidate" << endl;
			for(int j = 0; j < candidate.size(); ++j) 
				cout << candidate[j] << " ";
			cout << endl;*/
			
			 
			//chosen = simplepick(graph, usd, i, candidate,  coverage, root_hop);
		//	cout << "5:" << chosen << "," << graph.vertices[chosen+1] - graph.vertices[chosen] << "," << graph.vertices[candidate[0]+1] - graph.vertices[candidate[0]] << endl;
		//	cout << i << "," << last_available << "," << candidate.size() << endl;
			//chosen = pick(graph, usd, i, candidate, coverage, acc_count, parent_tree, root_hop, last_hop, last_alive, depth);
	//		chosen = pick(graph, usd, i, candidate,  coverage, second_coverage);
			
		/*	
			if(usd[border[i+1]] == false)
				if(acc_count[chosen] - 1 > 0 && acc_count[border[i+1]] - 1 > 0 )
				cout << root_hop[chosen] << "," << coverage[chosen] << "," << (acc[chosen] - coverage[chosen])/(acc_count[chosen] - 1) << "\t" << root_hop[border[i+1]] << "," << coverage[border[i+1]] << "," << (acc[border[i+1]] - coverage[border[i+1]])/(acc_count[border[i+1]]-1) << endl;
			
			chosen = border[i+1];
			*/
			
		//	last_iter_cover = coverage[chosen];
		//	last_iter_hop = root_hop[chosen];
			//if( i == taken )
			//	cout << i << "," << root_hop[chosen] << "," << last_hop[chosen] << endl;
			
			clear_tmp(descendants, coverage, parent_tree, root_hop, depth);
		//	clear_tmp(second_descendants, second_coverage, second_parent_tree, second_root_hop, second_depth);
		}
		
	/*	for(unordered_map<int, double>::iterator it = over_stats.begin(); it != over_stats.end(); ++it){
			cout << (*it).first << ":" << (double)((*it).second / actual_stats[(*it).first]) << endl;
		}
		*/
		vector<NodeID> remaining;
		remaining.reserve(numOfVertices - taken);
		
		
		for(NodeID i = 0; i < numOfVertices; ++i){
			if(usd[deg_inv[i]] == false)
				remaining.push_back(deg_inv[i]);
		}
		
		for(NodeID i = taken; i < numOfVertices; ++i){
			chosen = remaining[i-taken];
		if(usd[chosen] == true) cout << "shit" << endl;
			inv[i] = chosen;
			rank[chosen] = i;
			labeling_source_bfs(chosen, parent_tree, coverage, descendants, root_hop, last_hop, que, vis, dst_r, usd, i, tmp_idx, graph);
			clear_tmp(descendants, coverage, parent_tree, root_hop, depth);
		}
		
		
		labels.index_.resize(numOfVertices);
		
		double amount = 0; 
		for (NodeID v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			
			amount += k ;
			
			labels.index_[v].spt_v.resize(k);
			labels.index_[v].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_d[i] = tmp_idx[v].second[i];
			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();

			tmp_idx[v].first.shrink_to_fit();
			tmp_idx[v].second.shrink_to_fit();
		}
		
		cout << "avg label size:" << (double)(amount)/(double)numOfVertices - 1 << endl;

	}
	
	void undirected_weighted_sigpath(WGraph& wgraph){ 
		cout << "Building Coverage Ordering Based Labels" << endl;
		inv.resize(numOfVertices);
		rank.resize(numOfVertices);
		
		vector<NodeID> deg_inv;
		vector<NodeID> deg_rank;
		deg_inv.resize(numOfVertices);
		deg_rank.resize(numOfVertices);
		
		vector<pair<float, NodeID> > deg(numOfVertices);
		for (size_t v = 0; v < numOfVertices; ++v) {
				deg[v] = make_pair((wgraph.vertices[v + 1] - wgraph.vertices[v]) + float(rand()) / RAND_MAX, v);
		}
		sort(deg.rbegin(), deg.rend());
		for (size_t v = 0; v < numOfVertices; ++v) deg_inv[v] = deg[v].second;
		
		last_available = 0;
		
		vector<NodeID> parent_tree(numOfVertices, numOfVertices);
		vector<NodeID> root_hop(numOfVertices, 0);	
		vector<NodeID> coverage(numOfVertices, 0);
		vector<NodeID> depth(numOfVertices, 0);
		vector<NodeID> last_alive(numOfVertices, 0);		
		vector<NodeID> last_hop(numOfVertices, 0);
		vector<NodeID> descendants;
		descendants.reserve(numOfVertices);
		
		vector<NodeID> candidate;
		candidate.reserve(numOfVertices);
				
		vector<NodeID> acc_count;
		acc_count.resize(numOfVertices);
		
		vector<bool> vis(numOfVertices);
		queue<NodeID> visited_que;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);
		vector<bool> usd(numOfVertices, false);

		// Preparing basic structure for pl algorithm.
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));
				 
		NodeID chosen = deg_inv[0];
		 
		NodeID taken = numOfVertices; // /5;
		NodeID last_iter_cover = 0;
		NodeID last_iter_hop = 0;
		unordered_map<int, double> actual_stats;
		unordered_map<int, double> over_stats;
		
		vector<NodeID> upwardcoverage(numOfVertices, 0);
		
		long long currentsum = 0;
		
		for(NodeID i = 0; i < taken; ++i){
			inv[i] = chosen; 
			rank[chosen] = i;
			if(usd[chosen] == true) cout << "shit" << endl;
			NodeID actual_cover = labeling_source_dij(chosen, parent_tree, coverage, descendants, root_hop, last_hop, pqueue, visited_que, distances, vis, dst_r, usd, i, tmp_idx, wgraph);
		//	NodeID actual_cover = labeling_source_bfs(chosen, parent_tree, coverage, descendants, root_hop, last_hop, que, vis, dst_r, usd, i, tmp_idx, graph);
			currentsum += actual_cover;
		
			if(i == numOfVertices - 1 ) break;  
			
			calcover(descendants, parent_tree, coverage, root_hop, last_alive, depth);
			
			vector<NodeID> upwardsPow;
			NodeID max_degree = -1;
			vector<NodeID> sigpath = findSigPath(chosen, coverage, parent_tree, usd,upwardsPow, wgraph, max_degree);
		
		 
			if(sigpath.size()<=2){
					chosen = topcan(deg_inv, usd, last_available);
			}
			
			
			if(sigpath.size() > 2){
				NodeID maxDegree = -1;
				NodeID chosenhop = -1;
				if( (double)max_degree / (double)sigpath.size() < 2 && sigpath.size() > 3){			
					NodeID max_degree_v = -1;
					max_degree = -1;
					//NodeID min_degree = INF_WEIGHT;
					if(sigpath.size()>2){
						NodeID maxtrend = -1;
						NodeID chosentrend = 0;
						//cout << "Coverage Trends:";
						for(NodeID j = 0; j < sigpath.size() - 1; ++j){
							//NodeID curDegree = wgraph.vertices[sigpath[j] +1] - wgraph.vertices[sigpath[j]];
							//cout << " " << curDegree;
							//NodeID curDegree = wgraph.vertices[sigpath[j] +1] - wgraph.vertices[sigpath[j]];
							//cout << " " <<  j << ":" << abs(coverage[sigpath[j+1]] - coverage[sigpath[j]]);
							NodeID cd = wgraph.vertices[sigpath[j] +1] - wgraph.vertices[sigpath[j]];
							if(max_degree < cd){
								max_degree = cd; 
								max_degree_v = sigpath[j];
							}
							//if(min_degree > cd)
							//	min_degree = cd;
							if(maxtrend <  cd * abs(coverage[sigpath[j+1]] - coverage[sigpath[j]])){
								 maxtrend = cd * abs(coverage[sigpath[j+1]] - coverage[sigpath[j]]);
								chosentrend = j;
								chosen = sigpath[j];
							}
						}
					}
				}else{
					for(NodeID j = 0; j < sigpath.size(); ++j){
						NodeID curDegree = wgraph.vertices[sigpath[j] +1] - wgraph.vertices[sigpath[j]];
						if(curDegree > maxDegree){
							maxDegree = curDegree;
							chosen = sigpath[j];
							chosenhop = j;
						}
					}
				}
			}
			
		
		
	//	if(sigpath.size() < 10 && sigpath.size() >2 &&  max_degree / min_degree > 10)
		//	chosen = max_degree_v;  
		//cout << endl;
		//cout << chosentrend << endl;
		 
		//cout << "iteration:" << i << " - " << currentsum << " - " << usd[chosen] << "-" << chosenhop << " - " << sigpath.size() << endl;
			
			clear_tmp(descendants, coverage, parent_tree, root_hop, depth);
		}
		
		vector<NodeID> remaining;
		remaining.reserve(numOfVertices - taken);
		
		
		for(NodeID i = 0; i < numOfVertices; ++i){
			if(usd[deg_inv[i]] == false)
				remaining.push_back(deg_inv[i]);
		}
		
		for(NodeID i = taken; i < numOfVertices; ++i){
			chosen = remaining[i-taken];
		if(usd[chosen] == true) cout << "shit" << endl;
			inv[i] = chosen;
			rank[chosen] = i;
			labeling_source_dij(chosen, parent_tree, coverage, descendants, root_hop, last_hop, pqueue, visited_que, distances, vis, dst_r, usd, i, tmp_idx, wgraph);
			clear_tmp(descendants, coverage, parent_tree, root_hop, depth);
		}
		
		
		labels.index_.resize(numOfVertices);
		
		double amount = 0; 
		for (NodeID v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			
			amount += k ;
			
			labels.index_[v].spt_v.resize(k);
			labels.index_[v].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_d[i] = tmp_idx[v].second[i];
			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();

			tmp_idx[v].first.shrink_to_fit();
			tmp_idx[v].second.shrink_to_fit();
		}
		cout << "avg label size:" << (double)(amount)/(double)numOfVertices - 1 << endl;
	}
	
	void directed_unweighted_sigpath(Graph& graph){
		cout << "Building Coverage Ordering Based Labels" << endl;
		inv.resize(numOfVertices);
		rank.resize(numOfVertices);
		
		vector<NodeID> deg_inv;
		vector<NodeID> deg_rank;
		deg_inv.resize(numOfVertices);
		deg_rank.resize(numOfVertices);
		
		vector<pair<float, NodeID> > deg(numOfVertices);
		for (size_t v = 0; v < numOfVertices; ++v) {
				deg[v] = make_pair((graph.vertices[v + 1] - graph.vertices[v]) * (graph.r_vertices[v + 1] - graph.r_vertices[v])+ float(rand()) / RAND_MAX, v);
		}
		sort(deg.rbegin(), deg.rend());
		for (size_t v = 0; v < numOfVertices; ++v) deg_inv[v] = deg[v].second;
		
		last_available = 0;
		
		vector<NodeID> parent_tree(numOfVertices, numOfVertices);
		vector<NodeID> r_parent_tree(numOfVertices, numOfVertices);
		vector<NodeID> root_hop(numOfVertices, 0);	
		vector<NodeID> r_root_hop(numOfVertices, 0);	
		vector<NodeID> coverage(numOfVertices, 0);	
		vector<NodeID> r_coverage(numOfVertices, 0);
		vector<NodeID> depth(numOfVertices, 0);
		vector<NodeID> last_alive(numOfVertices, 0);		
		vector<NodeID> last_hop(numOfVertices, 0);
		vector<NodeID> descendants;
		descendants.reserve(numOfVertices);
		vector<NodeID> r_descendants;
		r_descendants.reserve(numOfVertices);
		
		vector<NodeID> candidate;
		candidate.reserve(numOfVertices);
				
		vector<NodeID> acc_count;
		acc_count.resize(numOfVertices);
		
		vector<bool> vis(numOfVertices);
		vector<bool> usd(numOfVertices, false);
		vector<NodeID> que(numOfVertices);
		// Preparing basic structure for pl algorithm.
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<EdgeWeight> r_dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));
				 
		NodeID chosen = deg_inv[0];
		 
		NodeID taken = numOfVertices;// / 10;
		NodeID last_iter_cover = 0;
		NodeID last_iter_hop = 0;
		unordered_map<int, double> actual_stats;
		unordered_map<int, double> over_stats;
		
		vector<NodeID> upwardcoverage(numOfVertices, 0);
		
		long long currentsum = 0;
		
		for(NodeID i = 0; i < taken; ++i){
			inv[i] = chosen; 
			rank[chosen] = i;
			if(usd[chosen] == true) cout << "shit" << endl;
			NodeID actual_cover =	labeling_source_bfs_directed(chosen, parent_tree, r_parent_tree, coverage, r_coverage, descendants, r_descendants, root_hop,  r_root_hop,  que,  vis,  dst_r, r_dst_r,  usd, i,  tmp_idx, r_tmp_idx,graph);
			currentsum += actual_cover;
		
			if(i == numOfVertices - 1 ) break;  
			
			calcover(descendants, parent_tree, coverage, root_hop);
			calcover(r_descendants, r_parent_tree, r_coverage, r_root_hop);
			

			vector<NodeID> sigpath = findRevSigPath(chosen, coverage, r_coverage, parent_tree, r_parent_tree, usd, graph);
			//vector<NodeID> sigrevpath = findSigRevPath(chosen, r_coverage, r_parent_tree, usd, graph);
		
		 
		if(sigpath.size()<=2){
			chosen = topcan(deg_inv, usd, last_available);
		}
			
		NodeID maxDegree = -1;
		NodeID minDegree = 	9999999999;
		NodeID chosenhop = -1;
		
		for(NodeID j = 0; j < sigpath.size(); ++j){
			NodeID fcd = graph.vertices[sigpath[j] +1] - graph.vertices[sigpath[j]];
			NodeID bcd = graph.r_vertices[sigpath[j] +1] - graph.r_vertices[sigpath[j]];
			NodeID curDegree = (graph.vertices[sigpath[j] +1] - graph.vertices[sigpath[j]]) * (graph.r_vertices[sigpath[j] +1] - graph.r_vertices[sigpath[j]]); 
			//NodeID curDegree = coverage[sigpath[j]] * r_coverage[sigpath[j]];
			//cout << j << ":" << curDegree << " ";
			if(curDegree > maxDegree){
				maxDegree = curDegree; 
				chosen = sigpath[j];
				chosenhop = j;
			}
			if(curDegree < minDegree)
				minDegree = curDegree;
		}
	//	if(maxDegree != -1)
	//		cout << maxDegree / minDegree << "," << sigpath.size() << endl;
		
		//cout << endl;	
		
		/*if(sigpath.size()>2){
			NodeID maxtrend = -1;
			NodeID chosentrend = 0; 
			//cout << "Coverage Trends:";
			for(NodeID j = 1; j < sigpath.size() - 1; ++j){
				//NodeID curDegree = wgraph.vertices[sigpath[j] +1] - wgraph.vertices[sigpath[j]];
				//cout << " " << curDegree;
				//NodeID curDegree = wgraph.vertices[sigpath[j] +1] - wgraph.vertices[sigpath[j]];
				//cout << " " <<  j << ":" << abs(coverage[sigpath[j+1]] - coverage[sigpath[j]]);
				if(maxtrend <  abs(coverage[sigpath[j+1]]   - coverage[sigpath[j]] )){
					 maxtrend = abs(coverage[sigpath[j+1]]   - coverage[sigpath[j]]);
					chosentrend = j;
					chosen = sigpath[j];
					chosenhop = j;
				}
			} 
		} */
		
		if(usd[chosen] == true){
			for(NodeID j = 0; j < sigpath.size() ; ++j)
				cout << usd[sigpath[j]] << ",";
			cout << endl;
		}
		
	/*	for(NodeID j = 0; j < sigpath.size(); ++j){
			cout << j << ":" << coverage[sigpath[j]] << " ";
		}
		cout << endl;
		
		for(NodeID j = 0; j < sigpath.size(); ++j){
			cout << j << ":" << r_coverage[sigpath[j]] << " ";
		}		
		cout << endl;
		
		for(NodeID j = 0; j < sigpath.size() - 1; ++j){
			cout << j << ":" << abs(coverage[sigpath[j]] - coverage[sigpath[j+1]]) << " ";
		}
		cout << endl;
		
		for(NodeID j = 0; j < sigpath.size() - 1; ++j){
			cout << j << ":" << abs(r_coverage[sigpath[j]] - r_coverage[sigpath[j+1]]) << " ";
		}
		
		cout << endl;
		//cout << chosentrend << endl;
		*/
			//cout << "iteration:" << i << " - " << currentsum << " - " << usd[chosen] << "-" << chosenhop << " - " << sigpath.size() << endl;
			
			clear_tmp(descendants, coverage, parent_tree, root_hop, depth);
			clear_tmp(r_descendants, r_coverage, r_parent_tree, r_root_hop, depth);
		}
		
		vector<NodeID> remaining;
		remaining.reserve(numOfVertices - taken);
		
		
		for(NodeID i = 0; i < numOfVertices; ++i){
			if(usd[deg_inv[i]] == false)
				remaining.push_back(deg_inv[i]);
		}
		
		for(NodeID i = taken; i < numOfVertices; ++i){
			chosen = remaining[i-taken];
		if(usd[chosen] == true) cout << "shit" << endl;
			inv[i] = chosen;
			rank[chosen] = i;
			labeling_source_bfs_directed(chosen, parent_tree, r_parent_tree, coverage, r_coverage, descendants, r_descendants, root_hop,  r_root_hop,  que,  vis,  dst_r, r_dst_r,  usd, i,  tmp_idx, r_tmp_idx,graph);
			clear_tmp(descendants, coverage, parent_tree, root_hop, depth);
			clear_tmp(r_descendants, r_coverage, r_parent_tree, r_root_hop, depth);
		}
		
		
		dlabels.index_.resize(numOfVertices);
		dlabels.bindex_.resize(numOfVertices);
		
		double amount = 0; 		
		for (NodeID v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			
			amount += k ;
			
			dlabels.index_[v].spt_v.resize(k);
			dlabels.index_[v].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_d[i] = tmp_idx[v].second[i];
			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();

			tmp_idx[v].first.shrink_to_fit();
			tmp_idx[v].second.shrink_to_fit();
   
			k = r_tmp_idx[v].first.size();
			
			amount += k ;
			
			dlabels.bindex_[v].spt_v.resize(k);
			dlabels.bindex_[v].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_v[i] = r_tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_d[i] = r_tmp_idx[v].second[i];
			r_tmp_idx[v].first.clear();
			r_tmp_idx[v].second.clear();

			r_tmp_idx[v].first.shrink_to_fit();
			r_tmp_idx[v].second.shrink_to_fit();
		}
		amount = (double)(amount) / 2;
		cout << "avg label size:" << (double)(amount)/(double)numOfVertices - 1 << endl;
	}
	
	void directed_weighted_sigpath(WGraph& wgraph){
		cout << "Building Coverage Ordering Based Labels" << endl;
		inv.resize(numOfVertices);
		rank.resize(numOfVertices);
		
		vector<NodeID> deg_inv;
		vector<NodeID> deg_rank;
		deg_inv.resize(numOfVertices);
		deg_rank.resize(numOfVertices);
		
		vector<pair<float, NodeID> > deg(numOfVertices);
		for (size_t v = 0; v < numOfVertices; ++v) {
				deg[v] = make_pair((wgraph.vertices[v + 1] - wgraph.vertices[v]) * (wgraph.r_vertices[v + 1] - wgraph.r_vertices[v])+ float(rand()) / RAND_MAX, v);
		}
		sort(deg.rbegin(), deg.rend());
		for (size_t v = 0; v < numOfVertices; ++v) deg_inv[v] = deg[v].second;
		
		last_available = 0;
		
		vector<NodeID> parent_tree(numOfVertices, numOfVertices);
		vector<NodeID> r_parent_tree(numOfVertices, numOfVertices);
		vector<NodeID> root_hop(numOfVertices, 0);	
		vector<NodeID> r_root_hop(numOfVertices, 0);	
		vector<NodeID> coverage(numOfVertices, 0);	
		vector<NodeID> r_coverage(numOfVertices, 0);
		vector<NodeID> depth(numOfVertices, 0);
		vector<NodeID> last_alive(numOfVertices, 0);		
		vector<NodeID> last_hop(numOfVertices, 0);
		vector<NodeID> descendants;
		descendants.reserve(numOfVertices);
		vector<NodeID> r_descendants;
		r_descendants.reserve(numOfVertices);
		
		vector<NodeID> candidate;
		candidate.reserve(numOfVertices);
				
		vector<NodeID> acc_count;
		acc_count.resize(numOfVertices);
		
		vector<bool> vis(numOfVertices);
		vector<bool> usd(numOfVertices, false);
		queue<NodeID> visited_que;
		
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);
		
		// Preparing basic structure for pl algorithm.
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<EdgeWeight> r_dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));
				 
		NodeID chosen = deg_inv[0];
		 
		NodeID taken = numOfVertices;// / 10;
		NodeID last_iter_cover = 0;
		NodeID last_iter_hop = 0;
		unordered_map<int, double> actual_stats;
		unordered_map<int, double> over_stats;
		
		vector<NodeID> upwardcoverage(numOfVertices, 0);
		
		long long currentsum = 0;
		
		for(NodeID i = 0; i < taken; ++i){
			inv[i] = chosen; 
			rank[chosen] = i;
			if(usd[chosen] == true) cout << "shit" << endl;
			NodeID actual_cover =	labeling_source_dij_directed(chosen, parent_tree, r_parent_tree, coverage, r_coverage, descendants, r_descendants, root_hop,  r_root_hop,  pqueue,  visited_que, distances, vis,  dst_r, r_dst_r,  usd, i,  tmp_idx, r_tmp_idx, wgraph);
			currentsum += actual_cover;
			
			if(i == numOfVertices - 1 ) break;   
			
			calcover(descendants, parent_tree, coverage, root_hop);
			calcover(r_descendants, r_parent_tree, r_coverage, r_root_hop);
			
			NodeID max_degree;
			bool reverse_flag;
			vector<NodeID> sigpath = findRevSigPath(chosen, coverage, r_coverage, parent_tree, r_parent_tree, usd, wgraph, max_degree, reverse_flag);
			//vector<NodeID> sigrevpath = findSigRevPath(chosen, r_coverage, r_parent_tree, usd, graph);
		
		 
			if(sigpath.size()<=2){
				chosen = topcan(deg_inv, usd, last_available);
			}
				
			NodeID maxDegree = -1;
			NodeID chosenhop = -1; 
			if(sigpath.size() > 2){
				
				if( (double)max_degree / (double)sigpath.size() < 2 && sigpath.size() > 3){			
					NodeID maxtrend = -1;
					NodeID chosentrend = 0;
					//cout << "Coverage Trends:";
					for(NodeID j = 0; j < sigpath.size() - 1; ++j){
						NodeID cd;
						NodeID ccoverage;
						if(reverse_flag == false){
							ccoverage = abs(coverage[sigpath[j+1]] - coverage[sigpath[j]]);
							cd = wgraph.r_vertices[sigpath[j] +1] - wgraph.r_vertices[sigpath[j]];
						}
						else{
							ccoverage = abs(r_coverage[sigpath[j+1]] - r_coverage[sigpath[j]]);
							cd = wgraph.vertices[sigpath[j] +1] - wgraph.vertices[sigpath[j]];
						}
						if(maxtrend <  cd * ccoverage){
							maxtrend = cd * ccoverage;
							chosentrend = j;
							chosen = sigpath[j];
						}
					} 
				}
				else{
					max_degree = -1;
					for(NodeID j = 0; j < sigpath.size(); ++j){
						NodeID curDegree = (wgraph.vertices[sigpath[j] +1] - wgraph.vertices[sigpath[j]]) * (wgraph.r_vertices[sigpath[j] +1] - wgraph.r_vertices[sigpath[j]]); 
						//NodeID curDegree = coverage[sigpath[j]] * r_coverage[sigpath[j]];
						//cout << j << ":" << curDegree << " ";
						if(curDegree > maxDegree){
							maxDegree = curDegree; 
							chosen = sigpath[j];
							chosenhop = j;
						}
					}
				}
				
			}
		
			
		
	//	NodeID max_degree_v = -1;
	//	NodeID max_degree = -1;
		//NodeID min_degree = INF_WEIGHT;
		/*if(sigpath.size()>2){
			NodeID maxtrend = -1;
			NodeID chosentrend = 0;
			//cout << "Coverage Trends:";
			for(NodeID j = 0; j < sigpath.size() - 1; ++j){
				//NodeID curDegree = wgraph.vertices[sigpath[j] +1] - wgraph.vertices[sigpath[j]];
				//cout << " " << curDegree;
				//NodeID curDegree = wgraph.vertices[sigpath[j] +1] - wgraph.vertices[sigpath[j]];
				//cout << " " <<  j << ":" << abs(coverage[sigpath[j+1]] - coverage[sigpath[j]]);
				NodeID cd = wgraph.vertices[sigpath[j] +1] - wgraph.vertices[sigpath[j]];
			//	if(max_degree < cd){
			//		max_degree = cd; 
			//		max_degree_v = sigpath[j];
				//}
				//if(min_degree > cd)
				//	min_degree = cd;
				if(maxtrend <  cd * (abs(coverage[sigpath[j+1]] - coverage[sigpath[j]]) + abs(r_coverage[sigpath[j+1]] - r_coverage[sigpath[j]]))){
					 maxtrend = cd * (abs(coverage[sigpath[j+1]] - coverage[sigpath[j]]) + abs(r_coverage[sigpath[j+1]] - r_coverage[sigpath[j]]));
					chosentrend = j;
					chosen = sigpath[j];
				}
			}
		}*/
		
		if(usd[chosen] == true){
			for(NodeID j = 0; j < sigpath.size() ; ++j)
				cout << usd[sigpath[j]] << ",";
			cout << endl;
		}
		
	/*	for(NodeID j = 0; j < sigpath.size(); ++j){
			cout << j << ":" << coverage[sigpath[j]] << " ";
		}
		cout << endl;
		
		for(NodeID j = 0; j < sigpath.size(); ++j){
			cout << j << ":" << r_coverage[sigpath[j]] << " ";
		}		
		cout << endl;
		
		for(NodeID j = 0; j < sigpath.size() - 1; ++j){
			cout << j << ":" << abs(coverage[sigpath[j]] - coverage[sigpath[j+1]]) << " ";
		}
		cout << endl;
		
		for(NodeID j = 0; j < sigpath.size() - 1; ++j){
			cout << j << ":" << abs(r_coverage[sigpath[j]] - r_coverage[sigpath[j+1]]) << " ";
		}
		
		cout << endl;
		//cout << chosentrend << endl;
		*/
			//cout << "iteration:" << i << " - " << currentsum << " - " << usd[chosen] << "-" << chosenhop << " - " << sigpath.size() << endl;
			
			clear_tmp(descendants, coverage, parent_tree, root_hop, depth);
			clear_tmp(r_descendants, r_coverage, r_parent_tree, r_root_hop, depth);
		}
		
		vector<NodeID> remaining;
		remaining.reserve(numOfVertices - taken);
		
		
		for(NodeID i = 0; i < numOfVertices; ++i){
			if(usd[deg_inv[i]] == false)
				remaining.push_back(deg_inv[i]);
		}
		
		for(NodeID i = taken; i < numOfVertices; ++i){
			chosen = remaining[i-taken];
		if(usd[chosen] == true) cout << "shit" << endl;
			inv[i] = chosen;
			rank[chosen] = i;
			labeling_source_dij_directed(chosen, parent_tree, r_parent_tree, coverage, r_coverage, descendants, r_descendants, root_hop,  r_root_hop,  pqueue,  visited_que, distances, vis,  dst_r, r_dst_r,  usd, i,  tmp_idx, r_tmp_idx, wgraph);
			clear_tmp(descendants, coverage, parent_tree, root_hop, depth);
			clear_tmp(r_descendants, r_coverage, r_parent_tree, r_root_hop, depth);
		}
		
		
		dlabels.index_.resize(numOfVertices);
		dlabels.bindex_.resize(numOfVertices);
		
		double amount = 0; 		
		for (NodeID v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			
			amount += k ;
			
			dlabels.index_[v].spt_v.resize(k);
			dlabels.index_[v].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_d[i] = tmp_idx[v].second[i];
			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();

			tmp_idx[v].first.shrink_to_fit();
			tmp_idx[v].second.shrink_to_fit();

			k = r_tmp_idx[v].first.size();
			
			amount += k ;
			
			dlabels.bindex_[v].spt_v.resize(k);
			dlabels.bindex_[v].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_v[i] = r_tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_d[i] = r_tmp_idx[v].second[i];
			r_tmp_idx[v].first.clear();
			r_tmp_idx[v].second.clear();

			r_tmp_idx[v].first.shrink_to_fit();
			r_tmp_idx[v].second.shrink_to_fit();
		}
		amount = (double)(amount) / 2;
		cout << "avg label size:" << (double)(amount)/(double)numOfVertices - 1 << endl;
	}
	
	/*	
	Coverage_Ordering(Graph& graph){
		cout << "Building Coverage Ordering Based Labels" << endl;
		inv.resize(numOfVertices);
		rank.resize(numOfVertices);
		
		acc.resize(numOfVertices);
		acc_count.resize(numOfVertices);
		
		vector<NodeID> deg_inv;
		vector<NodeID> deg_rank;
		deg_inv.resize(numOfVertices);
		deg_rank.resize(numOfVertices);
		
		vector<pair<float, NodeID> > deg(numOfVertices);
		for (size_t v = 0; v < numOfVertices; ++v) {
			//if (DIRECTED_FLAG == true)
				//deg[v] = make_pair((graph.vertices[v+1] - graph.vertices[v]) * (graph.r_vertices[v + 1] - graph.r_vertices[v]) + float(rand()) / RAND_MAX, v); // In-degree + Out-degree.
			//else
				deg[v] = make_pair((graph.vertices[v + 1] - graph.vertices[v]) + float(rand()) / RAND_MAX, v);
		}
		sort(deg.rbegin(), deg.rend());
		for (size_t v = 0; v < numOfVertices; ++v) deg_inv[v] = deg[v].second;
		
		last_available = 0;
		
		vector<NodeID> parent_tree(numOfVertices, numOfVertices);
		vector<NodeID> root_hop(numOfVertices, 0);	
		vector<NodeID> coverage(numOfVertices, 0);
		vector<NodeID> depth(numOfVertices, 0);
		vector<NodeID> candidate;
		candidate.reserve(numOfVertices);
		

			
		vector<NodeID> descendants;
		descendants.reserve(numOfVertices);
		vector<NodeID> que(numOfVertices);
		vector<bool> vis(numOfVertices);
		vector<bool> usd(numOfVertices, false);
// Preparing basic structure for pl algorithm.
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));
				 
				
		 		
		NodeID chosen = deg_inv[0];
		
		
		
		for(NodeID i = 0; i < numOfVertices; ++i){
			inv[i] = chosen; 
			rank[chosen] = i;
			
			labeling_source_bfs(chosen, parent_tree, coverage, descendants, root_hop, que, vis, dst_r, usd, i, tmp_idx, graph);
			if(i == numOfVertices - 1 ) break; 
			calcover(descendants, parent_tree, coverage, root_hop, depth);
			
			candidate_gen(candidate, deg_inv, usd);
			
		//	cout << i << "," << last_available << "," << candidate.size() << endl;
 
			
			chosen = pick(i, candidate, coverage, acc_count, parent_tree, root_hop, depth);
					
			
			clear_tmp(descendants, coverage, parent_tree, root_hop, depth);
		}
		
		labels.index_.resize(numOfVertices);
		
		double amount = 0;
		for (NodeID v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			
			amount += k ;
			
			labels.index_[v].spt_v.resize(k);
			labels.index_[v].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_d[i] = tmp_idx[v].second[i];
			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();

			tmp_idx[v].first.shrink_to_fit();
			tmp_idx[v].second.shrink_to_fit();
		}
		
		cout << "avg label size:" << (double)(amount)/(double)numOfVertices << endl;
		
		
	} 
	*/ 
	
	Coverage_Ordering(Graph& graph, vector<NodeID> border, bool test){
		cout << "Building Coverage Ordering Based Labels" << endl;
		inv.resize(numOfVertices);
		rank.resize(numOfVertices);
		
		
		vector<NodeID> deg_inv;
		vector<NodeID> deg_rank;
		deg_inv.resize(numOfVertices);
		deg_rank.resize(numOfVertices);
		
		vector<pair<float, NodeID> > deg(numOfVertices);
		for (size_t v = 0; v < numOfVertices; ++v) {
				deg[v] = make_pair((graph.vertices[v + 1] - graph.vertices[v]) + float(rand()) / RAND_MAX, v);
		}
		sort(deg.rbegin(), deg.rend());
		for (size_t v = 0; v < numOfVertices; ++v) deg_inv[v] = deg[v].second;
		
		last_available = 0;
		
		vector<NodeID> parent_tree(numOfVertices, numOfVertices);
		vector<NodeID> root_hop(numOfVertices, 0);	
		vector<NodeID> coverage(numOfVertices, 0);
		vector<NodeID> depth(numOfVertices, 0);
		vector<NodeID> last_alive(numOfVertices, 0);		
		vector<NodeID> last_hop(numOfVertices, 0);
		vector<NodeID> descendants;
		descendants.reserve(numOfVertices);
		
		vector<NodeID> second_parent_tree(numOfVertices, numOfVertices);
		vector<NodeID> second_root_hop(numOfVertices, 0);	
		vector<NodeID> second_coverage(numOfVertices, 0);
		vector<NodeID> second_depth(numOfVertices, 0);
		vector<NodeID> second_last_alive(numOfVertices, 0);		
		vector<NodeID> second_last_hop(numOfVertices, 0);
		vector<NodeID> second_descendants;
		second_descendants.reserve(numOfVertices);
				
				
		vector<NodeID> max_hop_set;
		
		vector<NodeID> candidate;
		candidate.reserve(numOfVertices);
				
		vector<NodeID> acc_count;
		acc_count.resize(numOfVertices);
		
		
		vector<NodeID> que(numOfVertices);
		vector<bool> vis(numOfVertices);
		vector<bool> usd(numOfVertices, false);
// Preparing basic structure for pl algorithm.
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));
				 
		 		
		NodeID chosen = deg_inv[0];
		
		NodeID taken = 2000;
		NodeID last_iter_cover = 0;
		NodeID last_iter_hop = 0;
		unordered_map<int, double> actual_stats;
		unordered_map<int, double> over_stats;
		
		for(NodeID i = 0; i < taken; ++i){
			inv[i] = chosen; 
			rank[chosen] = i;
			
			NodeID actual_cover = labeling_source_bfs(chosen, parent_tree, coverage, descendants, root_hop, last_hop, que, vis, dst_r, usd, i, tmp_idx, graph);
			
			if(last_iter_hop > 0){
				actual_stats[last_iter_hop] += actual_cover;
				over_stats[last_iter_hop] += last_iter_cover;
			}
			cout << last_iter_cover << " vs. " << actual_cover << " on " << last_iter_hop << " hops" << endl;
			
			if(i == numOfVertices - 1 ) break; 
			
			/*for(NodeID j = descendants.size() - 1; j >= 0; --j){
				if(c_max_hop == root_hop[descendants[j]])
					if(usd[descendants[j]] == false)
						max_hop_set.push_back(descendants[j]);
			}*/
			
		//	NodeID second_source = max_hop_set[rand()%max_hop_set.size()];
			
		//	build_tmp_tree(second_source, second_parent_tree, second_coverage, second_descendants, second_root_hop, second_last_hop, que, vis, dst_r, usd, i, tmp_idx, graph);
			
			calcover(descendants, parent_tree, coverage, root_hop, last_alive, depth);
			
		//	calcover_secondtree(second_descendants, second_parent_tree, parent_tree, second_coverage, second_root_hop, second_last_alive, second_depth);
			
			candidate_gen(candidate, deg_inv, usd);
			
		//	cout << i << "," << last_available << "," << candidate.size() << endl;
			chosen = pick(graph, usd, i, candidate, coverage, acc_count, parent_tree, root_hop, last_hop, last_alive, depth);
	//		chosen = pick(graph, usd, i, candidate,  coverage, second_coverage);
	
		/*	
			if(usd[border[i+1]] == false)
				if(acc_count[chosen] - 1 > 0 && acc_count[border[i+1]] - 1 > 0 )
				cout << root_hop[chosen] << "," << coverage[chosen] << "," << (acc[chosen] - coverage[chosen])/(acc_count[chosen] - 1) << "\t" << root_hop[border[i+1]] << "," << coverage[border[i+1]] << "," << (acc[border[i+1]] - coverage[border[i+1]])/(acc_count[border[i+1]]-1) << endl;
			
			chosen = border[i+1];
			*/
			
			last_iter_cover = coverage[chosen];
			last_iter_hop = root_hop[chosen];
			//if( i == taken )
			//	cout << i << "," << root_hop[chosen] << "," << last_hop[chosen] << endl;
			
			clear_tmp(descendants, coverage, parent_tree, root_hop, depth);
		//	clear_tmp(second_descendants, second_coverage, second_parent_tree, second_root_hop, second_depth);
		}
		
		for(unordered_map<int, double>::iterator it = over_stats.begin(); it != over_stats.end(); ++it){
			cout << (*it).first << ":" << (double)((*it).second / actual_stats[(*it).first]) << endl;
		}
		
		vector<NodeID> remaining;
		remaining.reserve(numOfVertices - taken);
		
		
		for(NodeID i = 0; i < numOfVertices; ++i){
			if(usd[deg_inv[i]] == false)
				remaining.push_back(deg_inv[i]);
		}
		
		for(NodeID i = taken; i < numOfVertices; ++i){
			chosen = remaining[i-taken];
			inv[i] = chosen;
			rank[chosen] = i;
			labeling_source_bfs(chosen, parent_tree, coverage, descendants, root_hop, last_hop, que, vis, dst_r, usd, i, tmp_idx, graph);
			clear_tmp(descendants, coverage, parent_tree, root_hop, depth);
		}
		
		
		labels.index_.resize(numOfVertices);
		
		double amount = 0; 
		for (NodeID v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			
			amount += k ;
			
			labels.index_[v].spt_v.resize(k);
			labels.index_[v].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_d[i] = tmp_idx[v].second[i];
			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();

			tmp_idx[v].first.shrink_to_fit();
			tmp_idx[v].second.shrink_to_fit();
		}
		
		cout << "avg label size:" << (double)(amount)/(double)numOfVertices - 1<< endl;
		
		
	} 
	
	Coverage_Ordering(Graph& graph, vector<NodeID> border){
		walk_stats(graph, border);
	}
	
	Coverage_Ordering(WGraph& wgraph){
		undirected_weighted_sigpath(wgraph);
	}
	
	Coverage_Ordering(WGraph& wgraph, bool directed){
		undirected_weighted_sigpath(wgraph);
	}
	
	Coverage_Ordering(Graph& graph, bool D_flags){
		directed_unweighted_sigpath(graph);
	}
	
	Coverage_Ordering(WGraph& wgraph, bool D_flags, bool W_flags){
		directed_weighted_sigpath(wgraph);
	}
};

#endif
