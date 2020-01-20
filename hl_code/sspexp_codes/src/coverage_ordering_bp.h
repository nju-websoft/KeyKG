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
#ifndef COVER_ORDERING_BP_H
#define COVER_ORDERING_BP_H

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
  
class COrdering_BP {
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

	~COrdering_BP() {
		inv.clear();
		rank.clear();
	}

};


template<int kNumBitParallelRoots = 50>
class Coverage_Ordering_BP : public COrdering_BP {
	typedef	vector<NodeID> tree;
	public:
		BPLabel<kNumBitParallelRoots> bplabels;
		DBPLabel<kNumBitParallelRoots> dbplabels;
		NodeID last_available;
		
		//tree parent_tree: while a vertex v is not in the tree, parent_tree[v] = numOfVertices
	
	NodeID labeling_source_bfs(NodeID source, tree& parent_tree, vector<NodeID>& coverage, vector<NodeID>& descendants, vector<NodeID>& root_hop, vector<NodeID>& last_hop, vector<NodeID>& que, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<bool>& usd, NodeID ranking, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<EdgeWeight>& tmp_d, vector<std::pair<uint64_t, uint64_t> >& tmp_s, vector<std::pair<NodeID, NodeID> >& sibling_es, vector<std::pair<NodeID, NodeID> >& child_es, NodeID& skip_count, NodeID& bp_seed, vector<NodeID>& inv, vector<NodeID>& rank, Graph& graph) { 
			descendants.clear();
			NodeID visited_arcs = 0;			

			index_t_bp<kNumBitParallelRoots>*& index_ = bplabels.index_bp;
			if (bp_seed < kNumBitParallelRoots) {
				
			//	NodeID i_bpspt = ranking - skip_count;
				
				usd[source] = true;
				fill(tmp_d.begin(), tmp_d.end(), INF_WEIGHT);
				fill(tmp_s.begin(), tmp_s.end(), std::make_pair(0, 0));
				
				int que_t0 = 0, que_t1 = 0, que_h = 0;
				que[que_h++] = source;
				tmp_d[source] = 0;
				que_t1 = que_h;
				
				
				int ns = 0;
				vector<pair<NodeID, NodeID> > adj_r(graph.vertices[source + 1] - graph.vertices[source]);
				for (EdgeID eid = graph.vertices[source]; eid < graph.vertices[source + 1]; eid++) {
					NodeID nw = graph.edges[eid];
					NodeID dnw = graph.vertices[nw+1] - graph.vertices[nw];
					adj_r[eid - graph.vertices[source]] = make_pair(dnw, nw);
				}

				sort(adj_r.rbegin(), adj_r.rend());
				
				for (size_t i = 0; i < adj_r.size(); ++i) {
					NodeID v = adj_r[i].second;
					if (!usd[v]) {
						skip_count++;
						usd[v] = true;
						rank[v] = skip_count + ranking;
						inv[skip_count+ranking] = v;
						que[que_h++] = v;
						tmp_d[v] = 1;
						tmp_s[v].first = 1ULL << ns;
						if (++ns == 64) break;
					}
				}
				
				for (EdgeWeight d = 0; que_t0 < que_h; ++d) {
					int num_sibling_es = 0, num_child_es = 0;

					for (int que_i = que_t0; que_i < que_t1; ++que_i) {
						NodeID v = que[que_i]; 

						for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; eid++) {
							NodeID tv = graph.edges[eid];
							EdgeWeight td = d + 1;

							if (d > tmp_d[tv]);
							else if (d == tmp_d[tv]) {
								if (v < tv) { 
									sibling_es[num_sibling_es].first = v;
									sibling_es[num_sibling_es].second = tv;
									++num_sibling_es;
								}
							}
							else {
								if (tmp_d[tv] == INF_WEIGHT) {
									que[que_h++] = tv;
									tmp_d[tv] = td;
								}
								child_es[num_child_es].first = v;
								child_es[num_child_es].second = tv;
								++num_child_es;
							}
						}
					}

					for (int i = 0; i < num_sibling_es; ++i) {
						int v = sibling_es[i].first, w = sibling_es[i].second;
						tmp_s[v].second |= tmp_s[w].first;
						tmp_s[w].second |= tmp_s[v].first;
					}
					for (int i = 0; i < num_child_es; ++i) {
						int v = child_es[i].first, c = child_es[i].second;
						tmp_s[c].first |= tmp_s[v].first;
						tmp_s[c].second |= tmp_s[v].second;
					}

					que_t0 = que_t1;
					que_t1 = que_h;
				}

				for (NodeID v = 0; v < numOfVertices; ++v) {
					index_[v].bpspt_d[bp_seed] = tmp_d[v];
					index_[v].bpspt_s[bp_seed][0] = tmp_s[v].first;
					index_[v].bpspt_s[bp_seed][1] = tmp_s[v].second & ~tmp_s[v].first;
				}
			}
			else{//Building Normal Labels				
				int max_hop = -1;
				
				const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];
				index_t_bp<kNumBitParallelRoots> &idx_r = index_[source];
						
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

						
						index_t_bp<kNumBitParallelRoots> &idx_v = index_[v];
						pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
					
						for (int i = 0; i < kNumBitParallelRoots; ++i) {
							EdgeWeight td = idx_r.bpspt_d[i] + idx_v.bpspt_d[i];
							if (td - 2 <= d) {
								/*td +=
								(idx_r.bpspt_s[i][0] & idx_v.bpspt_s[i][0]) ? -2 :
								((idx_r.bpspt_s[i][0] & idx_v.bpspt_s[i][1]) |
								(idx_r.bpspt_s[i][1] & idx_v.bpspt_s[i][0]))
								? -1 : 0;*/
								td +=
									(idx_r.bpspt_s[i][0] & idx_v.bpspt_s[i][0]) ? -2 :
									((idx_r.bpspt_s[i][0] & idx_v.bpspt_s[i][1]) |
										(idx_r.bpspt_s[i][1] & idx_v.bpspt_s[i][0]))
									? -1 : 0;
									if (td <= d) goto pruned_forward;
								}
						}

						
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
			}
			return visited_arcs;
		}
     
	NodeID labeling_source_bfs_directed(NodeID source, tree& parent_tree, tree& r_parent_tree, vector<NodeID>& coverage, vector<NodeID>& r_coverage, vector<NodeID>& descendants, vector<NodeID>& r_descendants,vector<NodeID>& root_hop, vector<NodeID>& r_root_hop, vector<NodeID>& que, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<EdgeWeight>& r_dst_r, vector<bool>& usd, vector<bool>& r_usd, NodeID ranking, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& r_tmp_idx, vector<EdgeWeight>& tmp_d, vector<EdgeWeight>& r_tmp_d, vector<std::pair<uint64_t, uint64_t> >& tmp_s, vector<std::pair<uint64_t, uint64_t> >& r_tmp_s,vector<std::pair<NodeID, NodeID> >& sibling_es,vector<std::pair<NodeID, NodeID> >& r_sibling_es, vector<std::pair<NodeID, NodeID> >& child_es, vector<std::pair<NodeID, NodeID> >& r_child_es, NodeID& skip_count, NodeID& bp_seed, vector<NodeID>& inv, vector<NodeID>& rank, Graph& graph) { 
			descendants.clear();
			r_descendants.clear();
			NodeID visited_arcs = 0;
			
			index_t_bp<kNumBitParallelRoots>*& index_ = dbplabels.index_bp;
			index_t_bp<kNumBitParallelRoots>*& bindex_ = dbplabels.bindex_bp;

			if (bp_seed < kNumBitParallelRoots) {
				
			//	NodeID i_bpspt = ranking - skip_count;
				
				r_usd[source] = true;
				fill(r_tmp_d.begin(), r_tmp_d.end(), INF_WEIGHT);
				fill(r_tmp_s.begin(), r_tmp_s.end(), std::make_pair(0, 0));
				
				int que_t0 = 0, que_t1 = 0, que_h = 0;
				que[que_h++] = source;
				r_tmp_d[source] = 0;
				que_t1 = que_h; 
				
				
				int ns = 0;
				vector<NodeID> adj_r(graph.vertices[source + 1] - graph.vertices[source]);
				for (EdgeID eid = graph.vertices[source]; eid < graph.vertices[source + 1]; eid++) {
					NodeID nw = graph.edges[eid];
					//NodeID dnw = (graph.vertices[nw+1] - graph.vertices[nw]) * (graph.r_vertices[nw+1] - graph.r_vertices[nw]);
					//adj_r[eid - graph.vertices[source]] = make_pair(dnw, nw);
					adj_r[eid - graph.vertices[source]] =  nw;
				}

				sort(adj_r.begin(), adj_r.end());
								
				
				vector<NodeID> r_adj_r(graph.r_vertices[source + 1] - graph.r_vertices[source]);
				for (EdgeID eid = graph.r_vertices[source]; eid < graph.r_vertices[source + 1]; eid++) {
					NodeID nw = graph.r_edges[eid];
					//NodeID dnw = (graph.vertices[nw+1] - graph.vertices[nw]) * (graph.r_vertices[nw+1] - graph.r_vertices[nw]);
					//r_adj_r[eid - graph.r_vertices[source]] = make_pair(dnw, nw);
					r_adj_r[eid - graph.r_vertices[source]] = nw;

				}

				sort(r_adj_r.begin(), r_adj_r.end());			
				
				
				vector<NodeID> common_adj;
				set_intersection(adj_r.begin(), adj_r.end(), r_adj_r.begin(), r_adj_r.end(), back_inserter(common_adj));
				
				// Adj is sorted by descending  value.
				vector<pair<NodeID, NodeID> > sort_by_on_degree(common_adj.size());		

				for (size_t i = 0; i < common_adj.size(); ++i) {
					NodeID v = common_adj[i];

					NodeID dv = (graph.vertices[v+1] - graph.vertices[v]) * (graph.r_vertices[v+1] - graph.r_vertices[v]);

					sort_by_on_degree[i] = make_pair(dv, v);
					//sort_by_on_sum[i] = make_pair((graph.vertices[v+1] - graph.vertices[v]) * (graph.r_vertices[v + 1] - graph.r_vertices[v]), v);
				}

				sort(sort_by_on_degree.rbegin(), sort_by_on_degree.rend());
				for (size_t i = 0; i < common_adj.size(); ++i)
					common_adj[i] = sort_by_on_degree[i].second;
				
				for (size_t i = 0; i < common_adj.size(); ++i) {
					NodeID v = common_adj[i];
					if (!r_usd[v]) {
						r_usd[v] = true;
						que[que_h++] = v;
						r_tmp_d[v] = 1;
						r_tmp_s[v].first = 1ULL << ns;
						
						skip_count++;
						rank[v] = skip_count + ranking;
						inv[skip_count+ranking] = v;
						if (++ns == 64) break;
					}
				}
				
				for (EdgeWeight d = 0; que_t0 < que_h; ++d) {
					int num_sibling_es = 0, num_child_es = 0;

					for (int que_i = que_t0; que_i < que_t1; ++que_i) {
						NodeID v = que[que_i]; 

						for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; eid++) {
							NodeID tv = graph.edges[eid];
							EdgeWeight td = d + 1;

							if (d > r_tmp_d[tv]);
							else if (d == r_tmp_d[tv]) {
							//	if (v < tv) { 
									r_sibling_es[num_sibling_es].first = v;
									r_sibling_es[num_sibling_es].second = tv;
									++num_sibling_es;
							//	}
							}
							else {
								if (r_tmp_d[tv] == INF_WEIGHT) {
									que[que_h++] = tv;
									r_tmp_d[tv] = td;
								}
								r_child_es[num_child_es].first = v;
								r_child_es[num_child_es].second = tv;
								++num_child_es;
							}
						}
					}

					for (int i = 0; i < num_sibling_es; ++i) {
						int v = r_sibling_es[i].first, w = r_sibling_es[i].second;
						//tmp_s[v].second |= tmp_s[w].first;
						r_tmp_s[w].second |= r_tmp_s[v].first;
					}
					for (int i = 0; i < num_child_es; ++i) {
						int v = r_child_es[i].first, c = r_child_es[i].second;
						r_tmp_s[c].first |= r_tmp_s[v].first;
						r_tmp_s[c].second |= r_tmp_s[v].second;
					}
 
					que_t0 = que_t1;
					que_t1 = que_h;
				}

				for (NodeID v = 0; v < numOfVertices; ++v) {
					bindex_[v].bpspt_d[bp_seed] = r_tmp_d[v];
					bindex_[v].bpspt_s[bp_seed][0] = r_tmp_s[v].first;
					bindex_[v].bpspt_s[bp_seed][1] = r_tmp_s[v].second & ~r_tmp_s[v].first;
				}
				
				
				//backward search	usd[r] = true;
				usd[source] = true;
				
				fill(tmp_d.begin(), tmp_d.end(), INF_WEIGHT);
				fill(tmp_s.begin(), tmp_s.end(), std::make_pair(0, 0));

				que_t0 = 0, que_t1 = 0, que_h = 0;
				que[que_h++] = source;
				tmp_d[source] = 0;
				que_t1 = que_h;

				ns = 0;
				
				for (size_t i = 0; i < common_adj.size(); ++i) {
					NodeID v = common_adj[i];
					if (!usd[v]) {
						usd[v] = true;
						que[que_h++] = v;
						tmp_d[v] = 1;
						tmp_s[v].first = 1ULL << ns;
						if (++ns == 64) break;
					}
				}
				
				for (EdgeWeight d = 0; que_t0 < que_h; ++d) {
					int num_sibling_es = 0, num_child_es = 0;

					for (int que_i = que_t0; que_i < que_t1; ++que_i) {
						NodeID v = que[que_i]; 

						for (EdgeID eid = graph.r_vertices[v]; eid < graph.r_vertices[v + 1]; eid++) {
							NodeID tv = graph.r_edges[eid];
							EdgeWeight td = d + 1;

							if (d > tmp_d[tv]);
							else if (d == tmp_d[tv]) {
								//if (v < tv) { 
									sibling_es[num_sibling_es].first = v;
									sibling_es[num_sibling_es].second = tv;
									++num_sibling_es;
								//}
							}
							else {
								if (tmp_d[tv] == INF_WEIGHT) {
									que[que_h++] = tv;
									tmp_d[tv] = td;
								}
								child_es[num_child_es].first = v;
								child_es[num_child_es].second = tv;
								++num_child_es;
							}
						}
					}

					for (int i = 0; i < num_sibling_es; ++i) {
						int v = sibling_es[i].first, w = sibling_es[i].second;
						//tmp_s[v].second |= tmp_s[w].first;
						tmp_s[w].second |= tmp_s[v].first;
					}
					for (int i = 0; i < num_child_es; ++i) {
						int v = child_es[i].first, c = child_es[i].second;
						tmp_s[c].first |= tmp_s[v].first;
						tmp_s[c].second |= tmp_s[v].second;
					}

					que_t0 = que_t1;
					que_t1 = que_h;
				}

				for (NodeID v = 0; v < numOfVertices; ++v) {
					index_[v].bpspt_d[bp_seed] = tmp_d[v];
					index_[v].bpspt_s[bp_seed][0] = tmp_s[v].first;
					index_[v].bpspt_s[bp_seed][1] = tmp_s[v].second & ~tmp_s[v].first;
				}
			}
			else{//Building Normal Labels	
			
				// Forward search.
				// Initialize forward labels of r.
				const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];
				index_t_bp<kNumBitParallelRoots> &idx_r = index_[source];

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
						index_t_bp<kNumBitParallelRoots> &r_idx_v = bindex_[v];

						//index_t &idx_v = index_[inv[v]];

						if (usd[v]) continue;

						for (int i = 0; i < kNumBitParallelRoots; ++i) {
							EdgeWeight td = idx_r.bpspt_d[i] + r_idx_v.bpspt_d[i];
							if (td - 2 <= d) {
								td +=
									(idx_r.bpspt_s[i][0] & r_idx_v.bpspt_s[i][0]) ? -2 :
									((idx_r.bpspt_s[i][0] & r_idx_v.bpspt_s[i][1]) |
										(idx_r.bpspt_s[i][1] & r_idx_v.bpspt_s[i][0]))
									? -1 : 0;
								if (td <= d) goto pruned_forward;
							}
						}
						
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
				index_t_bp<kNumBitParallelRoots> &r_idx_r = bindex_[source];

				
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
						index_t_bp<kNumBitParallelRoots> &idx_v = index_[v];
						//index_t &idx_v = index_[inv[v]];

						if (usd[v]) continue;

						for (int i = 0; i < kNumBitParallelRoots; ++i) {
							EdgeWeight td = r_idx_r.bpspt_d[i] + idx_v.bpspt_d[i];
							if (td - 2 <= d) {
								/*td +=
									(r_idx_r.bpspt_s[i][0] & idx_v.bpspt_s[i][0]) ? -2 :
									((r_idx_r.bpspt_s[i][0] & idx_v.bpspt_s[i][1]) |
										(r_idx_r.bpspt_s[i][1] & idx_v.bpspt_s[i][0]))
									? -1 : 0;*/
								td +=
									( idx_v.bpspt_s[i][0]) & r_idx_r.bpspt_s[i][0] ? -2 :
									((idx_v.bpspt_s[i][1] & r_idx_r.bpspt_s[i][0]) |
										(idx_v.bpspt_s[i][0] & r_idx_r.bpspt_s[i][1]))
									? -1 : 0;
								if (td <= d) goto pruned_backward;
							}
						}
						
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
				r_usd[source] = true; 
			}
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
				 
				 
		vector<EdgeWeight> tmp_d(numOfVertices);
		vector<std::pair<uint64_t, uint64_t> > tmp_s(numOfVertices);
		vector<std::pair<NodeID, NodeID> > sibling_es(numOfEdges);
		vector<std::pair<NodeID, NodeID> > child_es(numOfEdges);
		
		index_t_bp<kNumBitParallelRoots>*& index_ = bplabels.index_bp;
		index_ = (index_t_bp<kNumBitParallelRoots>*)memalign(64, numOfVertices * sizeof(index_t_bp<kNumBitParallelRoots>));
		
		
		NodeID chosen = deg_inv[0];
		
		NodeID taken = numOfVertices; // /10;
		NodeID last_iter_cover = 0;
		NodeID last_iter_hop = 0;
		unordered_map<int, double> actual_stats;
		unordered_map<int, double> over_stats;
		
		vector<NodeID> upwardcoverage(numOfVertices, 0);
		
		long long currentsum = 0;
		
		NodeID bp_seed = 0;
		
		for(NodeID i = 0; i < numOfVertices; ++i){
			
			inv[i] = chosen; 
			rank[chosen] = i;
			if(usd[chosen] == true) cout << "shit" << endl;
			int skip_count = 0;
			NodeID actual_cover = labeling_source_bfs(chosen, parent_tree, coverage, descendants, root_hop, last_hop, que, vis, dst_r, usd, i, tmp_idx, tmp_d, tmp_s, sibling_es, child_es, skip_count, bp_seed, inv, rank, graph);
			currentsum += actual_cover;
			bp_seed++;
			
			
			
			if(i == numOfVertices - 1 ) break;  
			

			
			calcover(descendants, parent_tree, coverage, root_hop, last_alive, depth);
			
			
			
			int iterround = 30;
			
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
			} 
			
			clear_tmp(descendants, coverage, parent_tree, root_hop, depth);
			
			i += skip_count;
		}
		
		double amount = 0; 
		for (NodeID v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			
			amount += k ; 
			
			
			index_[v].spt_v = (NodeID*)memalign(64, k * sizeof(NodeID));
			index_[v].spt_d = (EdgeWeight*)memalign(64, k * sizeof(EdgeWeight));
			for (NodeID i = 0; i < k; ++i) index_[v].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) index_[v].spt_d[i] = tmp_idx[v].second[i];
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
		vector<bool> r_usd(numOfVertices, false);
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
				 
				 
						 
		vector<EdgeWeight> tmp_d(numOfVertices);
		vector<std::pair<uint64_t, uint64_t> > tmp_s(numOfVertices);
		vector<EdgeWeight> r_tmp_d(numOfVertices);
		vector<std::pair<uint64_t, uint64_t> > r_tmp_s(numOfVertices);
		vector<std::pair<NodeID, NodeID> > sibling_es(numOfEdges);
		vector<std::pair<NodeID, NodeID> > child_es(numOfEdges);
		vector<std::pair<NodeID, NodeID> > r_sibling_es(numOfEdges);
		vector<std::pair<NodeID, NodeID> > r_child_es(numOfEdges);
		
		index_t_bp<kNumBitParallelRoots>*& index_ = dbplabels.index_bp;
		index_ = (index_t_bp<kNumBitParallelRoots>*)memalign(64, numOfVertices * sizeof(index_t_bp<kNumBitParallelRoots>));
		index_t_bp<kNumBitParallelRoots>*& bindex_ = dbplabels.bindex_bp;
		bindex_ = (index_t_bp<kNumBitParallelRoots>*)memalign(64, numOfVertices * sizeof(index_t_bp<kNumBitParallelRoots>));
		
		NodeID bp_seed = 0;
				 
		NodeID chosen = deg_inv[0];
		 
		NodeID taken = numOfVertices;// / 10;
		NodeID last_iter_cover = 0;
		NodeID last_iter_hop = 0;
		unordered_map<int, double> actual_stats;
		unordered_map<int, double> over_stats;
		
		vector<NodeID> upwardcoverage(numOfVertices, 0);
		
		long long currentsum = 0;
		
		for(NodeID i = 0; i < taken; ++i){
			//cout << i << endl;
			inv[i] = chosen; 
			rank[chosen] = i;
			if(usd[chosen] == true) cout << "shit" << endl;
			
			int skip_count = 0;
			NodeID actual_cover =	labeling_source_bfs_directed(chosen, parent_tree, r_parent_tree, coverage, r_coverage, descendants, r_descendants, root_hop,  r_root_hop,  que,  vis,  dst_r, r_dst_r,  usd, r_usd, i,  tmp_idx, r_tmp_idx, tmp_d,  r_tmp_d,  tmp_s,  r_tmp_s, sibling_es, r_sibling_es,  child_es,  r_child_es, skip_count, bp_seed,  inv,  rank, graph);
			currentsum += actual_cover;
			bp_seed++;
		
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
			
			i += skip_count;
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
			/*inv[i] = chosen;
			rank[chosen] = i;
			labeling_source_bfs_directed(chosen, parent_tree, r_parent_tree, coverage, r_coverage, descendants, r_descendants, root_hop,  r_root_hop,  que,  vis,  dst_r, r_dst_r,  usd, i,  tmp_idx, r_tmp_idx, tmp_d,  r_tmp_d,  tmp_s,  r_tmp_s, sibling_es, r_sibling_es,  child_es,  r_child_es, skip_count, bp_seed,  inv,  rank, graph);
			clear_tmp(descendants, coverage, parent_tree, root_hop, depth);
			clear_tmp(r_descendants, r_coverage, r_parent_tree, r_root_hop, depth);*/
		}	
				
		double amount = 0; 
		for (NodeID v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			
			amount += k ; 			
			index_[v].spt_v = (NodeID*)memalign(64, k * sizeof(NodeID));
			index_[v].spt_d = (EdgeWeight*)memalign(64, k * sizeof(EdgeWeight));
			for (NodeID i = 0; i < k; ++i) index_[v].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) index_[v].spt_d[i] = tmp_idx[v].second[i];
			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();

			tmp_idx[v].first.shrink_to_fit(); 
			tmp_idx[v].second.shrink_to_fit();
			 
			k = r_tmp_idx[v].first.size();
			amount += k ; 			
			bindex_[v].spt_v = (NodeID*)memalign(64, k * sizeof(NodeID));
			bindex_[v].spt_d = (EdgeWeight*)memalign(64, k * sizeof(EdgeWeight));
			for (NodeID i = 0; i < k; ++i) bindex_[v].spt_v[i] = r_tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) bindex_[v].spt_d[i] = r_tmp_idx[v].second[i];
			r_tmp_idx[v].first.clear();
			r_tmp_idx[v].second.clear();

			r_tmp_idx[v].first.shrink_to_fit();
			r_tmp_idx[v].second.shrink_to_fit();
			
		}
		
		
		amount = (double)(amount) / 2;
		cout << "avg label size:" << (double)(amount)/(double)numOfVertices - 1 << endl;
	}
	
	Coverage_Ordering_BP(Graph& graph, vector<NodeID> border){
		walk_stats(graph, border);
	}
	
	Coverage_Ordering_BP(Graph& graph, bool D_flags){
		directed_unweighted_sigpath(graph);
	}
		
};

#endif
