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
#ifndef COVER_ORDERING_COMPRESS_H
#define COVER_ORDERING_COMPRESS_H

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
#include <sstream> 
#include<bitset>

using namespace time_util;
  
class COrdering_C {
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

	~COrdering_C() {
		inv.clear();
		rank.clear();
	}

};



namespace std
{
    class nodeid_vector_hasher {
	public:
		std::size_t operator()(std::pair<NodeID, std::vector<NodeID> > const& pairvec) const {
		  std::size_t seed = pairvec.second.size();
		  seed ^= pairvec.first + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		  for(auto& i : pairvec.second) {
			seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		  }
		  return seed;
		}
};
}

//using boost::hash_combine
template <class T>
inline void hash_combine(std::size_t& seed, T const& v)
{
    seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}


class Coverage_Ordering_Compress : public COrdering_C {
	typedef	vector<NodeID> tree;
	public:
	//	PLabel plabels;
		NodeID last_available;
		CLabel clabels;
		long children_size;
		long r_children_size;
		bool SECOND_LEVEL = false;
		bool DIRECTED = false;
		
		
	NodeID labeling_source_bfs(NodeID source, tree& parent_tree, vector<NodeID>& coverage, vector<NodeID>& descendants, vector<NodeID>& root_hop, vector<NodeID>& last_hop, vector<NodeID>& que, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<bool>& usd, NodeID ranking, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<NodeID> > >&	tmp_idx_token_parents, Graph& graph) { 
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
					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_token_parents_v = tmp_idx_token_parents[v];
			
			
					for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
						NodeID w = tmp_idx_v.first[i];
						EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
						if(tmp_idx_v.second[i] == d + dst_r[w] && tmp_idx_token_parents_v.first[i] == numOfVertices){
							tmp_idx_token_parents_v.first[i] = ranking;//tmp_idx_token_parents_v.first.size();
							tmp_idx_token_parents_v.second[i] = dst_r[w];
						}						
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

					tmp_idx_token_parents_v.first.back() = numOfVertices;
					tmp_idx_token_parents_v.second.back() = INF_WEIGHT;
					tmp_idx_token_parents_v.first.push_back(numOfVertices);
					tmp_idx_token_parents_v.second.push_back(INF_WEIGHT);
					
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
	
		
		//tree parent_tree: while a vertex v is not in the tree, parent_tree[v] = numOfVertices
	NodeID labeling_source_dij(NodeID source, tree& parent_tree, vector<NodeID>& coverage, vector<NodeID>& descendants, vector<NodeID>& root_hop, vector<NodeID>& last_hop, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, queue<NodeID>& visited_que, vector<EdgeWeight>& distances, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<bool>& usd, NodeID ranking, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx,vector<pair<vector<NodeID>, vector<NodeID> > >& tmp_idx_parents, vector<NodeID>& parents, vector<pair<vector<NodeID>, vector<NodeID> > >&	tmp_idx_token_parents, WGraph& wgraph){
			
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
		parents[source] = source; 
		
		while (!pqueue.empty()) {

			NodeID v;
			EdgeWeight v_d;
			pqueue.extract_min(v, v_d);
			pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
			pair<vector<NodeID>, vector<NodeID> > &tmp_idx_parent_v = tmp_idx_parents[v];
			pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_token_parents_v = tmp_idx_token_parents[v];
			
			vis[v] = true;
			visited_que.push(v);

			if (usd[v]) continue;
			for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
				NodeID w = tmp_idx_v.first[i];
				EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
				if(tmp_idx_v.second[i] == v_d + dst_r[w] && tmp_idx_token_parents_v.first[i] == numOfVertices){
					tmp_idx_token_parents_v.first[i] = ranking;//tmp_idx_token_parents_v.first.size();
					tmp_idx_token_parents_v.second[i] = dst_r[w];
				}
				if (td <= v_d) {
					goto pruned;
				}
			}

			// Traverse
			tmp_idx_v.first.back() = ranking;
			tmp_idx_v.second.back() = v_d;
			tmp_idx_v.first.push_back(numOfVertices);
			tmp_idx_v.second.push_back(INF_WEIGHT);
			

			tmp_idx_parent_v.first.back() = ranking;
			tmp_idx_parent_v.second.back() = parents[v];
			tmp_idx_parent_v.first.push_back(numOfVertices);
			tmp_idx_parent_v.second.push_back(numOfVertices);
			
			tmp_idx_token_parents_v.first.back() = numOfVertices;
			tmp_idx_token_parents_v.second.back() = INF_WEIGHT;
			tmp_idx_token_parents_v.first.push_back(numOfVertices);
			tmp_idx_token_parents_v.second.push_back(INF_WEIGHT);
			
			
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
						parents[w] = v;
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
				parents[vis_v] = numOfVertices;
			}

			pqueue.clear_n();

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
			usd[source] = true;
			
			return visited_arcs;
	}
	
	NodeID labeling_source_bfs_directed(NodeID source, tree& parent_tree, tree& r_parent_tree, vector<NodeID>& coverage, vector<NodeID>& r_coverage, vector<NodeID>& descendants, vector<NodeID>& r_descendants,vector<NodeID>& root_hop, vector<NodeID>& r_root_hop, vector<NodeID>& que, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<EdgeWeight>& r_dst_r, vector<bool>& usd, NodeID ranking, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& r_tmp_idx, vector<pair<vector<NodeID>, vector<NodeID> > >&	tmp_idx_token_parents, vector<pair<vector<NodeID>, vector<NodeID> > >&	r_tmp_idx_token_parents,Graph& graph) { 
			descendants.clear();
			r_descendants.clear();
			NodeID visited_arcs = 0;
			
			// Forward search.
			// Initialize forward labels of r.
			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}
			
			
			const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[source];

			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i) {
				r_dst_r[r_tmp_idx_r.first[i]] = r_tmp_idx_r.second[i];
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
					pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_token_parents_v = r_tmp_idx_token_parents[v];

					//index_t &idx_v = index_[inv[v]];

					if (usd[v]) continue;

					// Pruned by the forward labels of r and backward labels of v in the forward search from r when reaching v.
					for (size_t i = 0; i < r_tmp_idx_v.first.size(); ++i) {
						NodeID w = r_tmp_idx_v.first[i];
						EdgeWeight td = r_tmp_idx_v.second[i] + dst_r[w];						
						if(r_tmp_idx_v.second[i] == d + r_dst_r[w] && r_tmp_idx_token_parents_v.first[i] == numOfVertices){
							r_tmp_idx_token_parents_v.first[i] = ranking;//tmp_idx_token_parents_v.first.size();
							r_tmp_idx_token_parents_v.second[i] = r_dst_r[w];
						}	
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
					
					r_tmp_idx_token_parents_v.first.back() = numOfVertices;
					r_tmp_idx_token_parents_v.second.back() = INF_WEIGHT;
					r_tmp_idx_token_parents_v.first.push_back(numOfVertices);
					r_tmp_idx_token_parents_v.second.push_back(INF_WEIGHT);
					
					
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


			// Backward search.
			// Initialize backward labels of r.


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
					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_token_parents_v = tmp_idx_token_parents[v];

					//index_t &idx_v = index_[inv[v]];

					if (usd[v]) continue;
 
					// Pruned by the backward labels of r and forward labels of v in the backward search from r when reaching v (v->r path).
					for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
						NodeID w = tmp_idx_v.first[i];
						EdgeWeight td = tmp_idx_v.second[i] + r_dst_r[w];
						if(tmp_idx_v.second[i] == d + dst_r[w] && tmp_idx_token_parents_v.first[i] == numOfVertices){
							tmp_idx_token_parents_v.first[i] = ranking;//tmp_idx_token_parents_v.first.size();
							tmp_idx_token_parents_v.second[i] = dst_r[w];
						}		
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

					tmp_idx_token_parents_v.first.back() = numOfVertices;
					tmp_idx_token_parents_v.second.back() = INF_WEIGHT;
					tmp_idx_token_parents_v.first.push_back(numOfVertices);
					tmp_idx_token_parents_v.second.push_back(INF_WEIGHT);
					
					
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
			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i)
				r_dst_r[r_tmp_idx_r.first[i]] = INF_WEIGHT;
 
			usd[source] = true;
			return visited_arcs;
		}
	
	
	NodeID labeling_source_dij_directed(NodeID source, tree& parent_tree, tree& r_parent_tree, vector<NodeID>& coverage, vector<NodeID>& r_coverage, vector<NodeID>& descendants, vector<NodeID>& r_descendants, vector<NodeID>& root_hop, vector<NodeID>& r_root_hop, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, queue<NodeID>& visited_que, vector<EdgeWeight>& distances, vector<bool>& vis, vector<EdgeWeight>& dst_r,vector<EdgeWeight>& r_dst_r, vector<bool>& usd, NodeID ranking, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& r_tmp_idx , vector<pair<vector<NodeID>, vector<NodeID> > >&	tmp_idx_token_parents, vector<pair<vector<NodeID>, vector<NodeID> > >&	r_tmp_idx_token_parents, WGraph& wgraph){
			
		descendants.clear();
		r_descendants.clear();
		NodeID visited_arcs = 0;
		
		const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];

		for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
		}
		
		
		const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[source];
		for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i) {
			r_dst_r[r_tmp_idx_r.first[i]] = r_tmp_idx_r.second[i];
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
			pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_token_parents_v = r_tmp_idx_token_parents[v];
			vis[v] = true;
			visited_que.push(v);
 
			if (usd[v]) continue;
			for (size_t i = 0; i < r_tmp_idx_v.first.size(); ++i) {
				NodeID w = r_tmp_idx_v.first[i];
				EdgeWeight td = r_tmp_idx_v.second[i] + dst_r[w];						
				if(r_tmp_idx_v.second[i] == v_d + r_dst_r[w] && r_tmp_idx_token_parents_v.first[i] == numOfVertices){
					r_tmp_idx_token_parents_v.first[i] = ranking;//tmp_idx_token_parents_v.first.size();
					r_tmp_idx_token_parents_v.second[i] = r_dst_r[w];
				}					
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
			
			r_tmp_idx_token_parents_v.first.back() = numOfVertices;
			r_tmp_idx_token_parents_v.second.back() = INF_WEIGHT;
			r_tmp_idx_token_parents_v.first.push_back(numOfVertices);
			r_tmp_idx_token_parents_v.second.push_back(INF_WEIGHT);

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

			
			
		//	const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[source];

	//		for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i) {
	//			dst_r[r_tmp_idx_r.first[i]] = r_tmp_idx_r.second[i];
	//		}

			pqueue.update(source, 0);
			distances[source] = 0;		
			
			
			// reverse search
			while (!pqueue.empty()) {

				NodeID v;
				EdgeWeight v_d;
				pqueue.extract_min(v, v_d);
				pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
				pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_token_parents_v = tmp_idx_token_parents[v];
				vis[v] = true;
				visited_que.push(v);

				if (usd[v]) continue;
				for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
					NodeID w = tmp_idx_v.first[i];
					//NodeID w = tmp_idx_v.first[i];
					EdgeWeight td = tmp_idx_v.second[i] + r_dst_r[w];
					if(tmp_idx_v.second[i] == v_d + dst_r[w] && tmp_idx_token_parents_v.first[i] == numOfVertices){
						tmp_idx_token_parents_v.first[i] = ranking;//tmp_idx_token_parents_v.first.size();
						tmp_idx_token_parents_v.second[i] = dst_r[w];
					}		
					
				//	EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
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
				
				
				tmp_idx_token_parents_v.first.back() = numOfVertices;
				tmp_idx_token_parents_v.second.back() = INF_WEIGHT;
				tmp_idx_token_parents_v.first.push_back(numOfVertices);
				tmp_idx_token_parents_v.second.push_back(INF_WEIGHT);

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

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;			
			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i)
				r_dst_r[r_tmp_idx_r.first[i]] = INF_WEIGHT;
			
			usd[source] = true;
			
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
	
	double normalize(NodeID mina, NodeID maxa, NodeID minb, NodeID maxb, NodeID a, NodeID b){
		double a_score = (double) a * (double)mina / (double) maxa;
		double b_score = (double) b * (double)minb / (double) maxb;
		return a_score + b_score;
	}
	
	string token_string(NodeID r, vector<NodeID>& token){
		stringstream ss;
			ss << r;
		for(NodeID i = 0; i < token.size(); ++i){
			ss << "_";
			ss << token[i];
		}
		return ss.str();
	}
	
	NodeID convertlist(NodeID& token_id, vector<token_t>& tokens_list, vector<vector<NodeID> >& tmp_tokens, vector<vector<EdgeWeight> >& tmp_tokens_distances, unordered_map<pair<NodeID, vector<NodeID> >, NodeID, nodeid_vector_hasher>& token_map, pair<vector<NodeID>, vector<EdgeWeight> > & tmp_idx_v, pair<vector<NodeID>, vector<NodeID> >& tmp_idx_parents_v, pair<vector<NodeID>, vector<NodeID> >&	tmp_idx_token_parents_v){
		
		NodeID lsize = tmp_idx_v.first.size();		
		NodeID anchorid = tmp_idx_v.first[lsize-2];
		
		for(NodeID i = 0; i < lsize - 1; ++i){
			//Add to parent_tree
			NodeID h = tmp_idx_v.first[i];
			NodeID hparent = tmp_idx_token_parents_v.first[i];
			EdgeWeight hparent_dis = tmp_idx_token_parents_v.second[i];
						
			NodeID tid = h;
			// non-trival tokens 
			if(tmp_tokens[h].size() != 0){
				vector<NodeID>& tmp_token_h = tmp_tokens[h];
				//string tstring = token_string(h, tmp_token_h);
				vector<EdgeWeight>& tmp_tokens_distances_h = tmp_tokens_distances[h];
				pair<NodeID, vector<NodeID> > token_key = make_pair(h, tmp_token_h);
				// New token 
				if(token_map.find(token_key) == token_map.end()){
					token_map[token_key] = token_id;
					token_t new_token;
					NodeID csize = tmp_token_h.size();
					new_token.sptc_v = (NodeID*)memalign(64, (csize + 1)* sizeof(NodeID));
					new_token.sptc_d = (EdgeWeight*)memalign(64, (csize + 1) * sizeof(EdgeWeight));			
					new_token.sptc_v[0] = h;
					new_token.sptc_d[0] = csize;
					children_size += (csize + 1);
					for(NodeID j = 0; j < csize; ++j){
						new_token.sptc_v[j+1] = tmp_token_h[j];
						new_token.sptc_d[j+1] = tmp_tokens_distances_h[j];
					}					
					tokens_list.push_back(new_token);
					tid = token_id;
					token_id++;
				}else // Already exist	
					tid = token_map[token_key];
			}
			
			//trival tokens
			if(i == lsize - 2)
				anchorid = tid;
			
			if(hparent!=numOfVertices){
				tmp_tokens[hparent].push_back(tid);
				tmp_tokens_distances[hparent].push_back(hparent_dis);
			}
		}
		
		
		return anchorid;
	}
	
	NodeID convertlist(NodeID& token_id, vector<token_t>& tokens_list, vector<vector<NodeID> >& tmp_tokens, vector<vector<EdgeWeight> >& tmp_tokens_distances, unordered_map<pair<NodeID, vector<NodeID> >, NodeID, nodeid_vector_hasher>& token_map, pair<vector<NodeID>, vector<EdgeWeight> > & tmp_idx_v, pair<vector<NodeID>, vector<NodeID> >&	tmp_idx_token_parents_v, bool REVERSE){
		
		NodeID lsize = tmp_idx_v.first.size();		
		NodeID anchorid = tmp_idx_v.first[lsize-2];
		
		for(NodeID i = 0; i < lsize - 1; ++i){
			//Add to parent_tree
			NodeID h = tmp_idx_v.first[i];
			NodeID hparent = tmp_idx_token_parents_v.first[i];
			EdgeWeight hparent_dis = tmp_idx_token_parents_v.second[i];
						
			NodeID tid = h;
			// non-trival tokens 
			if(tmp_tokens[h].size() != 0){
				vector<NodeID>& tmp_token_h = tmp_tokens[h];
				//string tstring = token_string(h, tmp_token_h);
				vector<EdgeWeight>& tmp_tokens_distances_h = tmp_tokens_distances[h];
				pair<NodeID, vector<NodeID> > token_key = make_pair(h, tmp_token_h);
				// New token 
				if(token_map.find(token_key) == token_map.end()){
					token_map[token_key] = token_id;
					token_t new_token;
					NodeID csize = tmp_token_h.size();
					new_token.sptc_v = (NodeID*)memalign(64, (csize + 1)* sizeof(NodeID));
					new_token.sptc_d = (EdgeWeight*)memalign(64, (csize + 1) * sizeof(EdgeWeight));			
					new_token.sptc_v[0] = h;
					new_token.sptc_d[0] = csize;
					if(REVERSE)
						r_children_size += (csize + 1);
					else
						children_size += (csize + 1);
					for(NodeID j = 0; j < csize; ++j){
						new_token.sptc_v[j+1] = tmp_token_h[j];
						new_token.sptc_d[j+1] = tmp_tokens_distances_h[j];
					}					
					tokens_list.push_back(new_token);
					tid = token_id;
					token_id++;
				}else // Already exist	
					tid = token_map[token_key];
			}
			
			//trival tokens
			if(i == lsize - 2)
				anchorid = tid;
			
			if(hparent!=numOfVertices){
				tmp_tokens[hparent].push_back(tid);
				tmp_tokens_distances[hparent].push_back(hparent_dis);
			}
		}
		
		
		return anchorid;
	}
	
	void validation(pair<vector<NodeID>, vector<EdgeWeight> > & tmp_idx_v, pair<vector<NodeID>, vector<NodeID> >& tmp_idx_token_parents_v, vector<token_t>& tokens_list, NodeID anchorid, vector<NodeID>& que){
		unordered_set<NodeID> cset;
		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;		
		que[que_h++] = anchorid;
		que_t1 = que_h;
		if(anchorid < numOfVertices) return;
		
		for (; que_t0 < que_h;) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID tid = que[que_i];
				
				const token_t& token_v = tokens_list[tid - numOfVertices];
				
				
				NodeID r = token_v.sptc_v[0];
				EdgeWeight csize = token_v.sptc_d[0];
				
				// hashing, can be replaced by 1024 linear probing for efficiency.
				cset.insert(r);
				
				for (EdgeWeight i = 0; i < csize; ++i){
					NodeID w = token_v.sptc_v[i+1];
					if( w < numOfVertices){
						// hashing, can be replaced by 1024 linear probing for efficiency.
						cset.insert(w);
					}else{
						que[que_h++] = w;
					}
				}
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}
		
		for(NodeID i = 0; i < tmp_idx_v.first.size() - 1;++i){
			if(cset.find(tmp_idx_v.first[i]) == cset.end()){
				cout << anchorid << "..." << tmp_idx_v.first[i] << "," << i << "," << tmp_idx_token_parents_v.first[i] << "," << cset.size() << "," << tmp_idx_token_parents_v.first.size() << endl;
			}
		}
		
	}
	
	void validation(pair<vector<NodeID>, vector<EdgeWeight> > & tmp_idx_v, pair<vector<NodeID>, vector<NodeID> >& tmp_idx_token_parents_v, vector<token_t>& tokens_list,vector<token_t>& supertokens, NodeID anchorid, vector<NodeID>& que){
		unordered_set<NodeID> cset;
		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;		
		que[que_h++] = anchorid;
		que_t1 = que_h;
		if(anchorid < numOfVertices) return;
		for (; que_t0 < que_h;) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID tid = que[que_i];				
				const token_t& token_v = tokens_list[tid - numOfVertices];
				
				
				NodeID r = token_v.sptc_v[0];					
				EdgeWeight ssize = token_v.sptc_d[0];
				cset.insert(r);
				token_t& supertoken_r = supertokens[r];
				
				NodeID rcsize = supertoken_r.sptc_v[0];
				
				EdgeWeight fsize = supertoken_r.sptc_d[0];
							
				EdgeWeight spos = 0;
				
				/*
				vector<bool>& one_level_t = one_level[tid - numOfVertices];
				
				for(NodeID i = 0; i < rcsize; ++i){
					
					if(one_level_t[i]){
						NodeID w = supertoken_r.sptc_v[i+1];
						if( w < numOfVertices){// hashing, can be replaced by 1024 linear probing for efficiency.
							cset.insert(w);										
						}else{
							que[que_h++] = w;
						}	
					}
				}*/
				
				
				for(EdgeWeight i = 0; i < fsize; ++i){
					unsigned char fmask = token_v.sptc_fbv[i];						
					bitset<8> fbs(fmask);
					for(NodeID j = 0; j < 8; ++j){							
						if(fbs[7-j]){
							unsigned char smask = token_v.sptc_sbv[spos++];
							bitset<8> sbs(smask);
							for(NodeID k = 0; k < 8; ++k){
								if(sbs[7-k]){
									//if( (i * 8 + j) * 8 + k  +  1 > rcsize) cout << (i * 8 + j) * 8 + k  +  1 << "vs" << rcsize << "vs" << fsize << ":" << fbs<< "," << sbs << " i:" << i << " j:" << j << " fbs[j]:" << fbs[8-j] <<   endl;
									
									NodeID w = supertoken_r.sptc_v[ (i * 8 + j) * 8 + k  +  1]; 
									
								//	if(anchorid == 9751)
									//	cout << w << "," << cset.size() << "," << tmp_idx_v.first.size() << endl; 
									//cout << r << "," << w << endl;
									if( w < numOfVertices){// hashing, can be replaced by 1024 linear probing for efficiency.
										cset.insert(w);										
									}else{
										que[que_h++] = w;
									}										
								}
							}
							//if(spos == ssize) break;
						}							
					}
					//if(spos == ssize) break;						
				}					
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}
		 		
		/*if(cset.size() != tmp_idx_v.first.size())
			cout << cset.size() << " vs. " << tmp_idx_v.first.size() << endl;
		
		for(NodeID i = 0; i < tmp_idx_v.first.size() - 1;++i){
			if(cset.find(tmp_idx_v.first[i]) == cset.end()){
				cout << anchorid << "..." << tmp_idx_v.first[i] << "," << i << "," << tmp_idx_token_parents_v.first[i] << "," << cset.size() << "," << tmp_idx_token_parents_v.first.size() << endl;
			}
		}*/
		
	}
	
	void validation(pair<vector<NodeID>, vector<EdgeWeight> > & tmp_idx_v, pair<vector<NodeID>, vector<NodeID> >& tmp_idx_token_parents_v, vector<token_t>& tokens_list, token_t* supertokens, NodeID anchorid, vector<NodeID>& que){
		unordered_set<NodeID> cset;
		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;		
		que[que_h++] = anchorid;
		que_t1 = que_h;
		if(anchorid < numOfVertices) return;
		for (; que_t0 < que_h;) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID tid = que[que_i];				
				const token_t& token_v = tokens_list[tid - numOfVertices];
				
				
				NodeID r = token_v.sptc_v[0];					
				EdgeWeight ssize = token_v.sptc_d[0];
				cset.insert(r);
				token_t& supertoken_r = supertokens[r];
				
				NodeID rcsize = supertoken_r.sptc_v[0];
				
				EdgeWeight fsize = supertoken_r.sptc_d[0];
							
				EdgeWeight spos = 0;
				
				/*
				vector<bool>& one_level_t = one_level[tid - numOfVertices];
				
				for(NodeID i = 0; i < rcsize; ++i){
					
					if(one_level_t[i]){
						NodeID w = supertoken_r.sptc_v[i+1];
						if( w < numOfVertices){// hashing, can be replaced by 1024 linear probing for efficiency.
							cset.insert(w);										
						}else{
							que[que_h++] = w;
						}	
					}
				}*/
				
				
				for(EdgeWeight i = 0; i < fsize; ++i){
					unsigned char fmask = token_v.sptc_fbv[i];						
					bitset<8> fbs(fmask);
					for(NodeID j = 0; j < 8; ++j){							
						if(fbs[7-j]){
							unsigned char smask = token_v.sptc_sbv[spos++];
							bitset<8> sbs(smask);
							for(NodeID k = 0; k < 8; ++k){
								if(sbs[7-k]){
									//if( (i * 8 + j) * 8 + k  +  1 > rcsize) cout << (i * 8 + j) * 8 + k  +  1 << "vs" << rcsize << "vs" << fsize << ":" << fbs<< "," << sbs << " i:" << i << " j:" << j << " fbs[j]:" << fbs[8-j] <<   endl;
									
									NodeID w = supertoken_r.sptc_v[ (i * 8 + j) * 8 + k  +  1]; 
									
								//	if(anchorid == 9751)
									//	cout << w << "," << cset.size() << "," << tmp_idx_v.first.size() << endl; 
									//cout << r << "," << w << endl;
									if( w < numOfVertices){// hashing, can be replaced by 1024 linear probing for efficiency.
										cset.insert(w);										
									}else{
										que[que_h++] = w;
									}										
								}
							}
							//if(spos == ssize) break;
						}							
					}
					//if(spos == ssize) break;						
				}					
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}
		 		
		if(cset.size() != tmp_idx_v.first.size()-1)
			cout << cset.size() << " vs. " << tmp_idx_v.first.size() << endl;
		/*
		for(NodeID i = 0; i < tmp_idx_v.first.size() - 1;++i){
			if(cset.find(tmp_idx_v.first[i]) == cset.end()){
				cout << anchorid << "..." << tmp_idx_v.first[i] << "," << i << "," << tmp_idx_token_parents_v.first[i] << "," << cset.size() << "," << tmp_idx_token_parents_v.first.size() << endl;
			}
		}*/
		
	}
	
	void validation(pair<vector<NodeID>, vector<EdgeWeight> > & tmp_idx_v, pair<vector<NodeID>, vector<NodeID> >& tmp_idx_token_parents_v, token_t* tokens_list, token_t* supertokens, NodeID anchorid, vector<NodeID>& que){
		unordered_set<NodeID> cset;
		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;		
		que[que_h++] = anchorid;
		que_t1 = que_h;
		if(anchorid < numOfVertices) return;
		for (; que_t0 < que_h;) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID tid = que[que_i];				
				const token_t& token_v = tokens_list[tid - numOfVertices];
				
				
				NodeID r = token_v.sptc_v[0];					
				EdgeWeight ssize = token_v.sptc_d[0];
				cset.insert(r);
				token_t& supertoken_r = supertokens[r];
				
				NodeID rcsize = supertoken_r.sptc_v[0];
				
				EdgeWeight fsize = supertoken_r.sptc_d[0];
							
				EdgeWeight spos = 0;
				
				/*
				vector<bool>& one_level_t = one_level[tid - numOfVertices];
				
				for(NodeID i = 0; i < rcsize; ++i){
					
					if(one_level_t[i]){
						NodeID w = supertoken_r.sptc_v[i+1];
						if( w < numOfVertices){// hashing, can be replaced by 1024 linear probing for efficiency.
							cset.insert(w);										
						}else{
							que[que_h++] = w;
						}	
					}
				}*/
				
				
				for(EdgeWeight i = 0; i < fsize; ++i){
					unsigned char fmask = token_v.sptc_fbv[i];						
					bitset<8> fbs(fmask);
					for(NodeID j = 0; j < 8; ++j){							
						if(fbs[7-j]){
							unsigned char smask = token_v.sptc_sbv[spos++];
							bitset<8> sbs(smask);
							for(NodeID k = 0; k < 8; ++k){
								if(sbs[7-k]){
									//if( (i * 8 + j) * 8 + k  +  1 > rcsize) cout << (i * 8 + j) * 8 + k  +  1 << "vs" << rcsize << "vs" << fsize << ":" << fbs<< "," << sbs << " i:" << i << " j:" << j << " fbs[j]:" << fbs[8-j] <<   endl;
									
									NodeID w = supertoken_r.sptc_v[ (i * 8 + j) * 8 + k  +  1]; 
									
								//	if(anchorid == 9751)
									//	cout << w << "," << cset.size() << "," << tmp_idx_v.first.size() << endl; 
									//cout << r << "," << w << endl;
									if( w < numOfVertices){// hashing, can be replaced by 1024 linear probing for efficiency.
										cset.insert(w);										
									}else{
										que[que_h++] = w;
									}										
								}
							}
							//if(spos == ssize) break;
						}							
					}
					//if(spos == ssize) break;						
				}					
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}
		 		
		if(cset.size() != tmp_idx_v.first.size() - 1)
			cout << cset.size() << " vs. " << tmp_idx_v.first.size() << endl;
		/*
		for(NodeID i = 0; i < tmp_idx_v.first.size() - 1;++i){
			if(cset.find(tmp_idx_v.first[i]) == cset.end()){
				cout << anchorid << "..." << tmp_idx_v.first[i] << "," << i << "," << tmp_idx_token_parents_v.first[i] << "," << cset.size() << "," << tmp_idx_token_parents_v.first.size() << endl;
			}
		}*/
		
	}
	
	
	void converttokens(vector<pair<vector<NodeID>, vector<EdgeWeight> > > & tmp_idx, vector<pair<vector<NodeID>, vector<NodeID> > >& tmp_idx_parents, vector<pair<vector<NodeID>, vector<NodeID> > >&	tmp_idx_token_parents){
		vector<token_t> tokens_list;
		vector<vector<NodeID> > tmp_tokens(numOfVertices);
		vector<vector<EdgeWeight> > tmp_tokens_distances(numOfVertices);
		unordered_map<pair<NodeID, vector<NodeID> >, NodeID, nodeid_vector_hasher> token_map;
		children_size = 0;
		NodeID token_id = numOfVertices;
		clabels.anchor_p = (NodeID*)memalign(64, (numOfVertices)* sizeof(NodeID));
//		vector<NodeID> que(numOfVertices, numOfVertices);
		for(NodeID v = 0; v < numOfVertices; ++v){
			clabels.anchor_p[v] = convertlist(token_id, tokens_list, tmp_tokens, tmp_tokens_distances, token_map, tmp_idx[v], tmp_idx_parents[v],tmp_idx_token_parents[v]);
			
			//validation(tmp_idx[v], tmp_idx_token_parents[v], tokens_list, clabels.anchor_p[v], que);
			
			NodeID lsize = tmp_idx[v].first.size();
			
			for(NodeID i = 0; i < lsize - 1; ++i){			
				NodeID h = tmp_idx[v].first[i];
				NodeID hparent = tmp_idx_token_parents[v].first[i];
				if( hparent != numOfVertices ){
					if(tmp_tokens[hparent].empty() == false){
						tmp_tokens[hparent].clear();
						tmp_tokens_distances[hparent].clear(); 
					}
				}
			}
		}
		
		clabels.numOfTokens = tokens_list.size();
		clabels.tokenindex_p = (token_t*)memalign(64, tokens_list.size() * sizeof(token_t));
		for(NodeID i = 0; i < tokens_list.size(); ++i){
			clabels.tokenindex_p[i] = tokens_list[i];
		}
		
	}
	
	void converttokens(vector<pair<vector<NodeID>, vector<EdgeWeight> > > & tmp_idx, vector<pair<vector<NodeID>, vector<NodeID> > >& tmp_idx_token_parents, vector<pair<vector<NodeID>, vector<EdgeWeight> > > & r_tmp_idx, vector<pair<vector<NodeID>, vector<NodeID> > >& r_tmp_idx_token_parents){
		vector<token_t> tokens_list;
		vector<token_t> r_tokens_list;
		vector<vector<NodeID> > tmp_tokens(numOfVertices);
		vector<vector<EdgeWeight> > tmp_tokens_distances(numOfVertices);
		
		vector<vector<NodeID> > r_tmp_tokens(numOfVertices);
		vector<vector<EdgeWeight> > r_tmp_tokens_distances(numOfVertices);
		
		unordered_map<pair<NodeID, vector<NodeID> >, NodeID, nodeid_vector_hasher> token_map;
		unordered_map<pair<NodeID, vector<NodeID> >, NodeID, nodeid_vector_hasher> r_token_map;
		children_size = 0;
		NodeID token_id = numOfVertices;
		NodeID r_token_id = numOfVertices;
		
		
		clabels.anchor_p = (NodeID*)memalign(64, (numOfVertices)* sizeof(NodeID));
		clabels.r_anchor_p = (NodeID*)memalign(64, (numOfVertices)* sizeof(NodeID));
			
		//vector<NodeID> que(numOfVertices, numOfVertices);
		for(NodeID v = 0; v < numOfVertices; ++v){
			clabels.anchor_p[v] = convertlist(token_id, tokens_list, tmp_tokens, tmp_tokens_distances, token_map, tmp_idx[v], tmp_idx_token_parents[v], false);
			clabels.r_anchor_p[v] = convertlist(r_token_id, r_tokens_list, r_tmp_tokens, r_tmp_tokens_distances, r_token_map, r_tmp_idx[v], r_tmp_idx_token_parents[v], true);
		//if(REVERSE == true)
		//	validation(tmp_idx[v], tmp_idx_token_parents[v], tokens_list, clabels.r_anchor_p[v], que);
		 
			NodeID lsize = tmp_idx[v].first.size();
			 
			for(NodeID i = 0; i < lsize - 1; ++i){			
				NodeID h = tmp_idx[v].first[i];
				NodeID hparent = tmp_idx_token_parents[v].first[i];
				if( hparent != numOfVertices ){
					if(tmp_tokens[hparent].empty() == false){
						tmp_tokens[hparent].clear();
						tmp_tokens_distances[hparent].clear(); 
					}
				}
			}
			
			lsize = r_tmp_idx[v].first.size();				 
			for(NodeID i = 0; i < lsize - 1; ++i){			
				NodeID h = r_tmp_idx[v].first[i];
				NodeID hparent = r_tmp_idx_token_parents[v].first[i];
				if( hparent != numOfVertices ){
					if(r_tmp_tokens[hparent].empty() == false){
						r_tmp_tokens[hparent].clear();
						r_tmp_tokens_distances[hparent].clear(); 
					}
				}
			}		
		}
		
		
		//vector<vector<bool> > one_level(tokens_list.size());
		
		clabels.numOfTokens = tokens_list.size();
		clabels.tokenindex_p = (token_t*)memalign(64, clabels.numOfTokens * sizeof(token_t));
		for(NodeID i = 0; i < clabels.numOfTokens; ++i){
			clabels.tokenindex_p[i] = tokens_list[i];				
		}
		
		clabels.r_numOfTokens = r_tokens_list.size();
		clabels.r_tokenindex_p = (token_t*)memalign(64, clabels.r_numOfTokens * sizeof(token_t));			
		for(NodeID i = 0; i < clabels.r_numOfTokens; ++i){
			clabels.r_tokenindex_p[i] = r_tokens_list[i]; 
		}			
		
					
		if(SECOND_LEVEL){	
			vector<token_t> supertokens(numOfVertices);					
			convertsupertokens(tokens_list, supertokens);
			
			clabels.supertokenindex_p = (token_t*)memalign(64, numOfVertices * sizeof(token_t));				
			for(NodeID i = 0; i < supertokens.size(); ++i){
				clabels.supertokenindex_p[i] = supertokens[i];
			}
			for(NodeID i = 0; i < clabels.numOfTokens; ++i){
				clabels.tokenindex_p[i] = tokens_list[i];				
			}
		
			//convertsupertokens(clabels.tokenindex_p, clabels.supertokenindex_p, clabels.numOfTokens);
			
		//	vector<NodeID> que(numOfVertices, numOfVertices);
		//	for(NodeID v = 0; v < numOfVertices; ++v){ 
				
				//if(REVERSE) { 
					//cout << v << endl;
					//validation(tmp_idx[v], tmp_idx_token_parents[v], tokens_list,supertokens,  clabels.anchor_p[v], que);	
					//validation(tmp_idx[v], tmp_idx_token_parents[v], tokens_list, clabels.supertokenindex_p,  clabels.anchor_p[v], que);	
			//		validation(tmp_idx[v], tmp_idx_token_parents[v], clabels.tokenindex_p, clabels.supertokenindex_p,  clabels.anchor_p[v], que);	

				//} 
			//}
			//convertsupertokens(clabels.r_tokenindex_p, clabels.r_supertokenindex_p, clabels.r_numOfTokens);
			
			vector<token_t> r_supertokens(numOfVertices);	
			convertsupertokens(r_tokens_list, r_supertokens); 	
			
			clabels.r_supertokenindex_p = (token_t*)memalign(64, numOfVertices * sizeof(token_t));
			for(NodeID i = 0; i < r_supertokens.size(); ++i){ 
				clabels.r_supertokenindex_p[i] = r_supertokens[i];
			} 
			for(NodeID i = 0; i < clabels.r_numOfTokens; ++i){
				clabels.r_tokenindex_p[i] = r_tokens_list[i]; 
			}	
			 
			//ector<NodeID> que(numOfVertices, numOfVertices);
			//for(NodeID v = 0; v < numOfVertices; ++v){
				 
				//if(REVERSE) {
					//cout << v << endl;
					//validation(tmp_idx[v], tmp_idx_token_parents[v], r_tokens_list, r_supertokens,  clabels.r_anchor_p[v], que);	
					//validation(r_tmp_idx[v], r_tmp_idx_token_parents[v], r_tokens_list, clabels.r_supertokenindex_p,  clabels.r_anchor_p[v], que);	
			//		validation(r_tmp_idx[v], r_tmp_idx_token_parents[v], clabels.r_tokenindex_p, clabels.r_supertokenindex_p,  clabels.r_anchor_p[v], que);	
				//} 
				//} 
				
			//}
		}
		
		
		
	}
	 
	
	void converttokens(vector<pair<vector<NodeID>, vector<EdgeWeight> > > & tmp_idx, vector<pair<vector<NodeID>, vector<NodeID> > >& tmp_idx_token_parents, bool REVERSE){
		vector<token_t> tokens_list;
		vector<token_t> r_tokens_list;
		vector<vector<NodeID> > tmp_tokens(numOfVertices);
		vector<vector<EdgeWeight> > tmp_tokens_distances(numOfVertices);
		
		vector<vector<NodeID> > r_tmp_tokens(numOfVertices);
		vector<vector<EdgeWeight> > r_tmp_tokens_distances(numOfVertices);
		
		unordered_map<pair<NodeID, vector<NodeID> >, NodeID, nodeid_vector_hasher> token_map;
		unordered_map<pair<NodeID, vector<NodeID> >, NodeID, nodeid_vector_hasher> r_token_map;
		if(REVERSE == false)
			children_size = 0;
		else
			r_children_size = 0;
		NodeID token_id = numOfVertices;
		NodeID r_token_id = numOfVertices;
		
		
		if(REVERSE == false)
			clabels.anchor_p = (NodeID*)memalign(64, (numOfVertices)* sizeof(NodeID));
		else
			clabels.r_anchor_p = (NodeID*)memalign(64, (numOfVertices)* sizeof(NodeID));
			
		//vector<NodeID> que(numOfVertices, numOfVertices);
		for(NodeID v = 0; v < numOfVertices; ++v){
			if(REVERSE == false)
				clabels.anchor_p[v] = convertlist(token_id, tokens_list, tmp_tokens, tmp_tokens_distances, token_map, tmp_idx[v], tmp_idx_token_parents[v], REVERSE);
			else
				clabels.r_anchor_p[v] = convertlist(r_token_id, r_tokens_list, r_tmp_tokens, r_tmp_tokens_distances, r_token_map, tmp_idx[v], tmp_idx_token_parents[v], REVERSE);
			//if(REVERSE == true)
			//	validation(tmp_idx[v], tmp_idx_token_parents[v], tokens_list, clabels.r_anchor_p[v], que);
			 
			 if(REVERSE == false){
				NodeID lsize = tmp_idx[v].first.size();
				 
				for(NodeID i = 0; i < lsize - 1; ++i){			
					NodeID h = tmp_idx[v].first[i];
					NodeID hparent = tmp_idx_token_parents[v].first[i];
					if( hparent != numOfVertices ){
						if(tmp_tokens[hparent].empty() == false){
							tmp_tokens[hparent].clear();
							tmp_tokens_distances[hparent].clear(); 
						}
					}
				}
			 }else{
				NodeID lsize = tmp_idx[v].first.size();				 
				for(NodeID i = 0; i < lsize - 1; ++i){			
					NodeID h = tmp_idx[v].first[i];
					NodeID hparent = tmp_idx_token_parents[v].first[i];
					if( hparent != numOfVertices ){
						if(r_tmp_tokens[hparent].empty() == false){
							r_tmp_tokens[hparent].clear();
							r_tmp_tokens_distances[hparent].clear(); 
						}
					}
				}
				 
			 }
		}
		
		
		//vector<vector<bool> > one_level(tokens_list.size());
		
		if(REVERSE == false){
			clabels.numOfTokens = tokens_list.size();
			clabels.tokenindex_p = (token_t*)memalign(64, clabels.numOfTokens * sizeof(token_t));
			for(NodeID i = 0; i < clabels.numOfTokens; ++i){
				clabels.tokenindex_p[i] = tokens_list[i];				
			}
		}else{
			clabels.r_numOfTokens = r_tokens_list.size();
			clabels.r_tokenindex_p = (token_t*)memalign(64, clabels.r_numOfTokens * sizeof(token_t));			
			for(NodeID i = 0; i < clabels.r_numOfTokens; ++i){
				clabels.r_tokenindex_p[i] = r_tokens_list[i]; 
			}			
		}
					
		if(REVERSE == false){
			if(SECOND_LEVEL){	
				vector<token_t> supertokens(numOfVertices);					
				convertsupertokens(tokens_list, supertokens);
				
				clabels.supertokenindex_p = (token_t*)memalign(64, numOfVertices * sizeof(token_t));				
				for(NodeID i = 0; i < supertokens.size(); ++i){
					clabels.supertokenindex_p[i] = supertokens[i];
				}
				for(NodeID i = 0; i < clabels.numOfTokens; ++i){
					clabels.tokenindex_p[i] = tokens_list[i]; 
				}	
				//convertsupertokens(clabels.tokenindex_p, clabels.supertokenindex_p, clabels.numOfTokens);
				
				//vector<NodeID> que(numOfVertices, numOfVertices);
				//for(NodeID v = 0; v < numOfVertices; ++v){
					
					//if(REVERSE) { 
						//cout << v << endl;
						//validation(tmp_idx[v], tmp_idx_token_parents[v], tokens_list,supertokens,  clabels.anchor_p[v], que);	
						//validation(tmp_idx[v], tmp_idx_token_parents[v], tokens_list, clabels.supertokenindex_p,  clabels.anchor_p[v], que);	
					//} 
				//}
			}  
		}else{ 
			if(SECOND_LEVEL){
				
				//convertsupertokens(clabels.r_tokenindex_p, clabels.r_supertokenindex_p, clabels.r_numOfTokens);
				
				vector<token_t> r_supertokens(numOfVertices);	
				convertsupertokens(r_tokens_list, r_supertokens); 	
				
				clabels.r_supertokenindex_p = (token_t*)memalign(64, numOfVertices * sizeof(token_t));
				for(NodeID i = 0; i < r_supertokens.size(); ++i){ 
					clabels.r_supertokenindex_p[i] = r_supertokens[i];
				} 
				
				for(NodeID i = 0; i < clabels.r_numOfTokens; ++i){
					clabels.r_tokenindex_p[i] = r_tokens_list[i]; 
				}	
				 /*
				//vector<NodeID> que(numOfVertices, numOfVertices);
				for(NodeID v = 0; v < numOfVertices; ++v){
					
					//if(REVERSE) {
						cout << v << endl;
						//validation(tmp_idx[v], tmp_idx_token_parents[v], r_tokens_list, r_supertokens,  clabels.r_anchor_p[v], que);	
						validation(tmp_idx[v], tmp_idx_token_parents[v], r_tokens_list, clabels.r_supertokenindex_p,  clabels.r_anchor_p[v], que);	
					//} 
				}*/
			}
		}
		
		
	}
	
	void convertsupertokens(vector<token_t>& tokens_list, vector<token_t>& supertokens){
		vector<unordered_map<NodeID, NodeID> > sorted_sp(numOfVertices);
		vector<unordered_map<NodeID, EdgeWeight> > dis_sp(numOfVertices);
		
				
		NodeID total_supertoken_children = 0;  
		
		for(NodeID t = 0 ; t < tokens_list.size(); ++t){
			token_t& token = tokens_list[t];
			NodeID r = token.sptc_v[0];
			EdgeWeight csize = token.sptc_d[0];			
			unordered_map<NodeID, NodeID>& rc_map = sorted_sp[r];	
			unordered_map<NodeID, EdgeWeight>& rd_map = dis_sp[r];		
			for(EdgeWeight i = 0; i < csize; ++i){
				NodeID cid = token.sptc_v[i + 1];
				rc_map[cid]++;
				rd_map[cid] = token.sptc_d[i + 1];
			}			
		}		 
		
		//supertokens.resize(numOfVertices);
		// Creating super tokens, sorting the children based on frequency
		for(NodeID v = 0; v < numOfVertices; ++v){
			vector<pair<NodeID, NodeID> > sorted_tmp;
			unordered_map<NodeID, NodeID>& vc_map = sorted_sp[v];
			unordered_map<NodeID, EdgeWeight>& vd_map = dis_sp[v];
			for(unordered_map<NodeID, NodeID>::iterator it = vc_map.begin(); it != vc_map.end(); ++it){
				sorted_tmp.push_back(make_pair((*it).second, (*it).first));
			}
			
			sort(sorted_tmp.rbegin(), sorted_tmp.rend());
			 
			token_t new_token;
			EdgeWeight csize = sorted_tmp.size();
			new_token.sptc_v = (NodeID*)memalign(64, (csize + 1)* sizeof(NodeID));
			new_token.sptc_d = (EdgeWeight*)memalign(64, (csize + 1) * sizeof(EdgeWeight));		
			new_token.sptc_v[0] = csize;
			new_token.sptc_d[0] = ceil((double)ceil((double)csize / (double)8) / (double)8);
			for(NodeID i = 0; i < csize; ++i){
				NodeID cid = sorted_tmp[i].second;
				new_token.sptc_v[i + 1] = cid;
				new_token.sptc_d[i + 1] = vd_map[cid];
			}
			
			supertokens[v] = new_token;
			total_supertoken_children += csize; 
		}
		
		
		// Converting each tokens to supertokens
		vector<bool> isChild(numOfVertices + tokens_list.size(), false);
		
		for(NodeID t = 0 ; t < tokens_list.size(); ++t){
			token_t& token = tokens_list[t];
			NodeID r = token.sptc_v[0];
			EdgeWeight csize = token.sptc_d[0];				
			
			if(csize == 0) continue;
			
		/*	vector<pair<NodeID, NodeID> > sorted_tmp;			
			unordered_map<NodeID, NodeID>& rc_map = sorted_sp[r];	*/
			
			for(EdgeWeight i = 0; i < csize; ++i){
				NodeID cid = token.sptc_v[i + 1];
				isChild[cid] = true;
			}			
			
			//sort(sorted_tmp.rbegin(), sorted_tmp.rend());
			
			const token_t& supertoken_r = supertokens[r];
			//vector<bool>& one_level_t  = one_level[t];
			
			//NodeID second_level_length = ceil((double)supertoken_r.sptc_d[0] / (double)8);
			NodeID first_level_length = ceil((double)supertoken_r.sptc_v[0] / (double)8);
			//cout << first_level_length << "," << second_level_length << endl;
			vector<bool> first_level_bv(first_level_length);
			vector<bool> second_level_bv;
			 
		//	NodeID t1 = 0;
			NodeID t2 = 1;
		//	NodeID tt = sorted_tmp[t1].second;
			NodeID st = supertoken_r.sptc_v[t2];
			
	//		NodeID ctsize = sorted_tmp.size();
			NodeID stsize = supertoken_r.sptc_v[0];
			
			//one_level_t.resize(stsize, false); 
			
			vector<bool> tmp_set(8, false);
			for(NodeID i = 0; i < first_level_length; ++i){
				NodeID in_batch = false;
				
				//if(t1 != ctsize && t2 != stsize)
				fill(tmp_set.begin(), tmp_set.end(), false);
				 
				for(NodeID j = 0; j < 8; ++j){
					
				//	if(t1 == ctsize) break;
					if(t2 == (stsize + 1)) break;
					 
				//	tt = sorted_tmp[t1].second;
					st = supertoken_r.sptc_v[t2];
					
					if(isChild[st]){
						tmp_set[j] = true;
						in_batch = true;
					}					
					t2++;					
				}
				if(in_batch == false) 
					first_level_bv[i] = false;
				else{
					first_level_bv[i] = true;					
					for(NodeID j = 0; j < 8; ++j)
						second_level_bv.push_back(tmp_set[j]);							
				}
			}
			
			 
			for(EdgeWeight i = 0; i < csize; ++i){
				NodeID cid = token.sptc_v[i + 1];
				isChild[cid] = false;
			}	
			
			NodeID first_level_int_length = ceil((double)first_level_length/(double)8); // bytes
			
			NodeID second_level_int_length = ceil((double)second_level_bv.size() / (double)8); // bytes
		//	if(t < 10) cout << "sl:" << r << "," << second_level_int_length << " vs " << second_level_bv.size() << " vs "  << token.sptc_d[0] <<  ";" << stsize << " vs " << first_level_length <<  endl;
			//supertoken_r.sptc_v[0] = stsize; //how many children for this supertoken
			//supertoken_r.sptc_d[0] = first_level_int_length; //how many uchar to store for this token referring to this supertoken = stsize / 8 / 8
			token.sptc_d[0] = second_level_int_length; //how many uchar to store for this token in second level
			
			// convert first_level_bv -> uint8_t* sptc_fbv
			// convert second_level_bv -> uint8_t* sptc_sbv
			// first_level_bv % 8 == 0; 
			// second_level_bv % 8 == 0;
			token.sptc_fbv = (unsigned char*)memalign(64, first_level_int_length * sizeof(unsigned char));
			token.sptc_sbv = (unsigned char*)memalign(64, second_level_int_length * sizeof(unsigned char));
			
		//	cout << "first:" << endl;
			for(NodeID i = 0; i < first_level_int_length; ++i){
				token.sptc_fbv[i] = 0;
				for(NodeID j = 0; j < 8; ++j){
					token.sptc_fbv[i] = token.sptc_fbv[i] << 1;
					if(first_level_bv[i * 8 + j]) ++token.sptc_fbv[i];
				}  
				/*
				bitset<8> x(token.sptc_fbv[i]);
				cout << x << endl;
				for(NodeID j = 0; j < 8; ++j){
					if(first_level_bv[i * 8 + j]) 
						cout << "1";
					else
						cout << "0";
				}			
				cout << endl;	*/	
			}
			//cout << endl;
			
			//cout << "second:" << endl;
			for(NodeID i = 0; i < second_level_int_length; ++i){
				token.sptc_sbv[i] = 0;
				for(NodeID j = 0; j < 8; ++j){
					token.sptc_sbv[i] = token.sptc_sbv[i] << 1;
					if(second_level_bv[i * 8 + j]) ++token.sptc_sbv[i];
				}
				
			//	bitset<8> x(token.sptc_sbv[i]);
			//	cout << x ;
			/*	for(NodeID j = 0; j < 8; ++j){
					if(second_level_bv[i * 8 + j]) 
						cout << "1";
					else
						cout << "0";
				}	*/
			//	cout << ",";
			}

			
			//cout << endl;	
			 
		}	 
		 
		
		 
		cout << " Number of Supertokens: " << supertokens.size() << endl;
		cout << " Average Children of Supertokens: " << (double)total_supertoken_children / (double) supertokens.size() << endl;
		
	}
	
	void convertsupertokens(token_t* tokens_list, token_t* supertokens, NodeID numOfTokens){
		vector<unordered_map<NodeID, NodeID> > sorted_sp(numOfVertices);
		vector<unordered_map<NodeID, EdgeWeight> > dis_sp(numOfVertices);
		
				
		NodeID total_supertoken_children = 0;  
		
		for(NodeID t = 0 ; t < numOfTokens; ++t){
			token_t& token = tokens_list[t];
			NodeID r = token.sptc_v[0];
			EdgeWeight csize = token.sptc_d[0];			
			unordered_map<NodeID, NodeID>& rc_map = sorted_sp[r];	
			unordered_map<NodeID, EdgeWeight>& rd_map = dis_sp[r];		
			for(EdgeWeight i = 0; i < csize; ++i){
				NodeID cid = token.sptc_v[i + 1];
				rc_map[cid]++;
				rd_map[cid] = token.sptc_d[i + 1];
			}			
		}		 
		
		supertokens = (token_t*)memalign(64, numOfVertices * sizeof(token_t));
		// Creating super tokens, sorting the children based on frequency
		for(NodeID v = 0; v < numOfVertices; ++v){
			vector<pair<NodeID, NodeID> > sorted_tmp;
			unordered_map<NodeID, NodeID>& vc_map = sorted_sp[v];
			unordered_map<NodeID, EdgeWeight>& vd_map = dis_sp[v];
			for(unordered_map<NodeID, NodeID>::iterator it = vc_map.begin(); it != vc_map.end(); ++it){
				sorted_tmp.push_back(make_pair((*it).second, (*it).first));
			}
			
			sort(sorted_tmp.rbegin(), sorted_tmp.rend());
			 
			token_t new_token;
			EdgeWeight csize = sorted_tmp.size();
			new_token.sptc_v = (NodeID*)memalign(64, (csize + 1)* sizeof(NodeID));
			new_token.sptc_d = (EdgeWeight*)memalign(64, (csize + 1) * sizeof(EdgeWeight));		
			new_token.sptc_v[0] = csize;
			new_token.sptc_d[0] = ceil((double)ceil((double)csize / (double)8) / (double)8);
			for(NodeID i = 0; i < csize; ++i){
				NodeID cid = sorted_tmp[i].second;
				new_token.sptc_v[i + 1] = cid;
				new_token.sptc_d[i + 1] = vd_map[cid];
			}
			
			supertokens[v] = new_token;
			total_supertoken_children += csize; 
		}
		
		
		// Converting each tokens to supertokens
		vector<bool> isChild(numOfVertices + numOfTokens, false);
		
		for(NodeID t = 0 ; t < numOfTokens; ++t){
			token_t& token = tokens_list[t];
			NodeID r = token.sptc_v[0];
			EdgeWeight csize = token.sptc_d[0];				
			
			if(csize == 0) continue;
			
		/*	vector<pair<NodeID, NodeID> > sorted_tmp;			
			unordered_map<NodeID, NodeID>& rc_map = sorted_sp[r];	*/
			
			for(EdgeWeight i = 0; i < csize; ++i){
				NodeID cid = token.sptc_v[i + 1];
				isChild[cid] = true;
			}			
			
			//sort(sorted_tmp.rbegin(), sorted_tmp.rend());
			
			const token_t& supertoken_r = supertokens[r];
			//vector<bool>& one_level_t  = one_level[t];
			
			//NodeID second_level_length = ceil((double)supertoken_r.sptc_d[0] / (double)8);
			NodeID first_level_length = ceil((double)supertoken_r.sptc_v[0] / (double)8);
			//cout << first_level_length << "," << second_level_length << endl;
			vector<bool> first_level_bv(first_level_length);
			vector<bool> second_level_bv;
			 
		//	NodeID t1 = 0;
			NodeID t2 = 1;
		//	NodeID tt = sorted_tmp[t1].second;
			NodeID st = supertoken_r.sptc_v[t2];
			
	//		NodeID ctsize = sorted_tmp.size();
			NodeID stsize = supertoken_r.sptc_v[0];
			
			//one_level_t.resize(stsize, false); 
			
			vector<bool> tmp_set(8, false);
			for(NodeID i = 0; i < first_level_length; ++i){
				NodeID in_batch = false;
				
				//if(t1 != ctsize && t2 != stsize)
				fill(tmp_set.begin(), tmp_set.end(), false);
				 
				for(NodeID j = 0; j < 8; ++j){
					
				//	if(t1 == ctsize) break;
					if(t2 == (stsize + 1)) break;
					 
				//	tt = sorted_tmp[t1].second;
					st = supertoken_r.sptc_v[t2]; 
					
					if(isChild[st]){
						tmp_set[j] = true;
						in_batch = true;
					}					
					t2++;					
				}
				if(in_batch == false) 
					first_level_bv[i] = false;
				else{
					first_level_bv[i] = true;					
					for(NodeID j = 0; j < 8; ++j)
						second_level_bv.push_back(tmp_set[j]);							
				}
			}
			
			 
			for(EdgeWeight i = 0; i < csize; ++i){
				NodeID cid = token.sptc_v[i + 1];
				isChild[cid] = false;
			}	
			
			NodeID first_level_int_length = ceil((double)first_level_length/(double)8); // bytes
			
			NodeID second_level_int_length = ceil((double)second_level_bv.size() / (double)8); // bytes
		//	if(t < 10) cout << "sl:" << r << "," << second_level_int_length << " vs " << second_level_bv.size() << " vs "  << token.sptc_d[0] <<  ";" << stsize << " vs " << first_level_length <<  endl;
			//supertoken_r.sptc_v[0] = stsize; //how many children for this supertoken
			//supertoken_r.sptc_d[0] = first_level_int_length; //how many uchar to store for this token referring to this supertoken = stsize / 8 / 8
			token.sptc_d[0] = second_level_int_length; //how many uchar to store for this token in second level
			
			// convert first_level_bv -> uint8_t* sptc_fbv
			// convert second_level_bv -> uint8_t* sptc_sbv
			// first_level_bv % 8 == 0; 
			// second_level_bv % 8 == 0;
			token.sptc_fbv = (unsigned char*)memalign(64, first_level_int_length * sizeof(unsigned char));
			token.sptc_sbv = (unsigned char*)memalign(64, second_level_int_length * sizeof(unsigned char));
			
		//	cout << "first:" << endl;
			for(NodeID i = 0; i < first_level_int_length; ++i){
				token.sptc_fbv[i] = 0;
				for(NodeID j = 0; j < 8; ++j){
					token.sptc_fbv[i] = token.sptc_fbv[i] << 1;
					if(first_level_bv[i * 8 + j]) ++token.sptc_fbv[i];
				}  
				/*
				bitset<8> x(token.sptc_fbv[i]);
				cout << x << endl;
				for(NodeID j = 0; j < 8; ++j){
					if(first_level_bv[i * 8 + j]) 
						cout << "1";
					else
						cout << "0";
				}			
				cout << endl;	*/	
			}
			//cout << endl;
			
			//cout << "second:" << endl;
			for(NodeID i = 0; i < second_level_int_length; ++i){
				token.sptc_sbv[i] = 0;
				for(NodeID j = 0; j < 8; ++j){
					token.sptc_sbv[i] = token.sptc_sbv[i] << 1;
					if(second_level_bv[i * 8 + j]) ++token.sptc_sbv[i];
				}
				
			//	bitset<8> x(token.sptc_sbv[i]);
			//	cout << x ;
			/*	for(NodeID j = 0; j < 8; ++j){
					if(second_level_bv[i * 8 + j]) 
						cout << "1";
					else
						cout << "0";
				}	*/
			//	cout << ",";
			}

			
			//cout << endl;	
			 
		}	 
		 
		
		 
		//cout << " Number of Supertokens: " << supertokens.size() << endl;
		cout << " Average Children of Supertokens: " << (double)total_supertoken_children / (double) numOfVertices << endl;
		
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
				 
				 
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx_token_parents(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, numOfVertices)));
				
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx_token_parents(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, numOfVertices)));		 
		
		
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
			NodeID actual_cover =	labeling_source_dij_directed(chosen, parent_tree, r_parent_tree, coverage, r_coverage, descendants, r_descendants, root_hop,  r_root_hop,  pqueue,  visited_que, distances, vis,  dst_r, r_dst_r,  usd, i,  tmp_idx, r_tmp_idx,tmp_idx_token_parents, r_tmp_idx_token_parents, wgraph);
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
			labeling_source_dij_directed(chosen, parent_tree, r_parent_tree, coverage, r_coverage, descendants, r_descendants, root_hop,  r_root_hop,  pqueue,  visited_que, distances, vis,  dst_r, r_dst_r,  usd, i,  tmp_idx, r_tmp_idx, tmp_idx_token_parents, r_tmp_idx_token_parents,wgraph);
			clear_tmp(descendants, coverage, parent_tree, root_hop, depth);
			clear_tmp(r_descendants, r_coverage, r_parent_tree, r_root_hop, depth);
		}
		
		
		converttokens(tmp_idx,  tmp_idx_token_parents, false);
		//converttokens(tmp_idx,  tmp_idx_token_parents, r_tmp_idx, r_tmp_idx_token_parents);
		cout << clabels.numOfTokens << " Tokens in total" << endl;
		cout << (double)children_size / (double) clabels.numOfTokens << " average children number" << endl;
		
		converttokens(r_tmp_idx,  r_tmp_idx_token_parents, true);
		cout << clabels.r_numOfTokens << " Tokens in total" << endl;
		cout << (double)r_children_size / (double) clabels.r_numOfTokens << " average children number" << endl;
		
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
		
		vector<pair<vector<NodeID>, vector<NodeID> > >
			tmp_idx_parents(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<NodeID>(1, numOfVertices)));
		
		
		
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx_token_parents(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, numOfVertices)));
				
		vector<NodeID> parents(numOfVertices, numOfVertices);
		
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
			NodeID actual_cover = labeling_source_dij(chosen, parent_tree, coverage, descendants, root_hop, last_hop, pqueue, visited_que, distances, vis, dst_r, usd, i, tmp_idx,tmp_idx_parents ,parents,tmp_idx_token_parents, wgraph);
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
			labeling_source_dij(chosen, parent_tree, coverage, descendants, root_hop, last_hop, pqueue, visited_que, distances, vis, dst_r, usd, i, tmp_idx,tmp_idx_parents, parents, tmp_idx_token_parents, wgraph);
			clear_tmp(descendants, coverage, parent_tree, root_hop, depth);
		}
		
		//converttokens(tmp_idx,  tmp_idx_parents, tmp_idx_token_parents);
		converttokens(tmp_idx,  tmp_idx_token_parents, false);
		cout << clabels.numOfTokens << " Tokens in total" << endl; 
		cout << (double)children_size / (double) clabels.numOfTokens << " average children number" << endl;
		
	}
		
	
	void walk_stats(Graph& graph){
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
					
			vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx_token_parents(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, numOfVertices)));
					 
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
				NodeID actual_cover = labeling_source_bfs(chosen, parent_tree, coverage, descendants, root_hop, last_hop, que, vis, dst_r, usd, i, tmp_idx, tmp_idx_token_parents, graph);
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
			labeling_source_bfs(chosen, parent_tree, coverage, descendants, root_hop, last_hop, que, vis, dst_r, usd, i, tmp_idx, tmp_idx_token_parents, graph);
			clear_tmp(descendants, coverage, parent_tree, root_hop, depth);
		}
		
		
		converttokens(tmp_idx,  tmp_idx_token_parents, false);
		cout << clabels.numOfTokens << " Tokens in total" << endl;
		cout << (double)children_size / (double) clabels.numOfTokens << " average children number" << endl;
		
		/*
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
		
		cout << "avg label size:" << (double)(amount)/(double)numOfVertices - 1 << endl;*/

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
				 
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx_token_parents(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, numOfVertices)));
				
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx_token_parents(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, numOfVertices)));
		 
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
			NodeID actual_cover =	labeling_source_bfs_directed(chosen, parent_tree, r_parent_tree, coverage, r_coverage, descendants, r_descendants, root_hop,  r_root_hop,  que,  vis,  dst_r, r_dst_r,  usd, i,  tmp_idx, r_tmp_idx, tmp_idx_token_parents, r_tmp_idx_token_parents, graph);
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

		
		if(usd[chosen] == true){
			for(NodeID j = 0; j < sigpath.size() ; ++j)
				cout << usd[sigpath[j]] << ",";
			cout << endl;
		}
			
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
			labeling_source_bfs_directed(chosen, parent_tree, r_parent_tree, coverage, r_coverage, descendants, r_descendants, root_hop,  r_root_hop,  que,  vis,  dst_r, r_dst_r,  usd, i,  tmp_idx, r_tmp_idx,  tmp_idx_token_parents, r_tmp_idx_token_parents, graph);
			clear_tmp(descendants, coverage, parent_tree, root_hop, depth);
			clear_tmp(r_descendants, r_coverage, r_parent_tree, r_root_hop, depth);
		}
		
		converttokens(tmp_idx,  tmp_idx_token_parents, false);
		//converttokens(tmp_idx,  tmp_idx_token_parents, r_tmp_idx, r_tmp_idx_token_parents);
		cout << clabels.numOfTokens << " Tokens in total" << endl;
		cout << (double)children_size / (double) clabels.numOfTokens << " average children number" << endl;
		
		converttokens(r_tmp_idx,  r_tmp_idx_token_parents, true);
		cout << clabels.r_numOfTokens << " Tokens in total" << endl;
		cout << (double)r_children_size / (double) clabels.r_numOfTokens << " average children number" << endl;
		
		//cout << "avg label size:" << (double)(amount)/(double)numOfVertices - 1 << endl;
	}
	
	Coverage_Ordering_Compress(Graph& graph, bool directed_flag = false, bool second_level_flag = false){
		if(second_level_flag)
			SECOND_LEVEL = true; 
		if(directed_flag)
			DIRECTED = true; 
		
		if(DIRECTED == false)
			walk_stats(graph);
		else{
			directed_unweighted_sigpath(graph);
			vector<NodeID> dis_vec(numOfVertices, INF_WEIGHT);
			vector<long> ts_vec(numOfVertices, -1);
			vector<NodeID> que(numOfVertices);
			vector<EdgeWeight> que_d(numOfVertices);
			/*
			if(SECOND_LEVEL)
				cout << clabels.query_p_two_level_d(71697, 40393, 0 , dis_vec, ts_vec, que, que_d) << endl;
			else
				cout << clabels.query_p_d(71697, 40393, 0 , dis_vec, ts_vec, que, que_d) << endl;
			*/
		}
	}

	
	Coverage_Ordering_Compress(WGraph& wgraph, bool directed_flag = false, bool second_level_flag = false){
		
		if(second_level_flag)
			SECOND_LEVEL = true; 
		if(directed_flag)
			DIRECTED = true; 
		if(DIRECTED_FLAG == false)
			undirected_weighted_sigpath(wgraph);
		else
			directed_weighted_sigpath(wgraph);
		
	
	}
	
	
};



#endif
