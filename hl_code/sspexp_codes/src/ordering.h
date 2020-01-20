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
#ifndef ORDERING_H
#define ORDERING_H

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
#include <unordered_map>
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

using namespace time_util;




class Ordering {
public:
	vector<NodeID> inv; // Fetch the original vertex id by a given ranking.
	vector<NodeID> rank; // Fetch the ranking of a given vertex id.
	
	void Relabel(Graph& graph) {
		
		for (NodeID v = 0; v < numOfVertices; ++v) rank[inv[v]] = v;
		
		/*vector<vector<NodeID> > new_adj(numOfVertices);
		vector<vector<NodeID> > new_r_adj(numOfVertices);

		for (NodeID v = 0; v < numOfVertices; ++v) {
			for (int i = 0; i < graph.adj[v].size(); ++i) {
				new_adj[rank[v]].push_back(rank[graph.adj[v][i]]);
			}
			if (DIRECTED_FLAG == true) {
				for (NodeID i = 0; i < graph.r_adj[v].size(); ++i) {
					new_r_adj[rank[v]].push_back(rank[graph.r_adj[v][i]]);
				}
			}
		}
		graph.adj.swap(new_adj);
		if (DIRECTED_FLAG == true) {
			graph.r_adj.swap(new_r_adj);
		}*/

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

	~Ordering() {
		inv.clear();
		rank.clear();
	}

};

class Degree_Ordering : public Ordering {
public:

	Degree_Ordering(Graph& graph) {
		inv.resize(numOfVertices);
		rank.resize(numOfVertices);

		vector<pair<float, NodeID> > deg(numOfVertices);
	//	srand((unsigned)time(NULL));
		srand(100);
		//for (size_t v = 0; v < numOfVertices; ++v) {
		//	if (DIRECTED_FLAG == true)
		//		deg[v] = make_pair(graph.adj[v].size() * graph.r_adj[v].size() + float(rand()) / RAND_MAX, v); // In-degree + Out-degree.
		//	else
		//		deg[v] = make_pair(graph.adj[v].size() + float(rand()) / RAND_MAX, v);
		//}
		for (size_t v = 0; v < numOfVertices; ++v) {
			if (DIRECTED_FLAG == true)
				deg[v] = make_pair((graph.vertices[v+1] - graph.vertices[v]) * (graph.r_vertices[v + 1] - graph.r_vertices[v]) + float(rand()) / RAND_MAX, v); // In-degree + Out-degree.
			else
				deg[v] = make_pair((graph.vertices[v + 1] - graph.vertices[v]) + float(rand()) / RAND_MAX, v);
		}
		sort(deg.rbegin(), deg.rend());
		for (size_t v = 0; v < numOfVertices; ++v) inv[v] = deg[v].second;
	
		Relabel(graph);	
	}	

	Degree_Ordering(WGraph& wgraph) {
		inv.resize(numOfVertices);
		rank.resize(numOfVertices);

		vector<pair<float, NodeID> > deg(numOfVertices);
		srand(100);
		//for (size_t v = 0; v < numOfVertices; ++v) {
		//	if (DIRECTED_FLAG == true)
		//		deg[v] = make_pair(wgraph.adj[v].size() + wgraph.r_adj[v].size() + float(rand()) / RAND_MAX, v); // In-degree + Out-degree.
		//	else
		//		deg[v] = make_pair(wgraph.adj[v].size() + float(rand()) / RAND_MAX, v);
		//}
		for (size_t v = 0; v < numOfVertices; ++v) {
			if (DIRECTED_FLAG == true)
				deg[v] = make_pair((wgraph.vertices[v + 1] - wgraph.vertices[v]) * (wgraph.r_vertices[v + 1] - wgraph.r_vertices[v]) + float(rand()) / RAND_MAX, v);
			else
				deg[v] = make_pair((wgraph.vertices[v + 1] - wgraph.vertices[v]) + float(rand()) / RAND_MAX, v);
		}
		sort(deg.rbegin(), deg.rend());
		for (size_t v = 0; v < numOfVertices; ++v) inv[v] = deg[v].second;

		Relabel(wgraph);
	}
};

class Degree_1_Ordering : public Ordering{
public:
	Degree_1_Ordering(Graph& graph) {
		inv.resize(numOfVertices, -1);
		rank.resize(numOfVertices, -1);
		vector<pair<float, NodeID> > deg(numOfVertices);

		for (size_t v = 0; v < numOfVertices; ++v) {
			if (DIRECTED_FLAG == true)
				deg[v] = make_pair(graph.adj[v].size() + graph.r_adj.size(), v); // In-degree + Out-degree.
			else
				deg[v] = make_pair(graph.adj[v].size() + float(rand()) / RAND_MAX, v);
		}
		sort(deg.rbegin(), deg.rend());

		unordered_set<NodeID> picked_this_round;
		NodeID picked_rank = 0;

		while (picked_rank != numOfVertices) {
			for (size_t i = 0; i < numOfVertices; ++i) {
				NodeID v = deg[i].second;
				if (inv[v] == -1 && picked_this_round.find(v) == picked_this_round.end()) {
					inv[v] = picked_rank++;
					for (size_t u = 0; u < graph.adj[v].size(); ++u) {
						picked_this_round.insert(graph.adj[v][u]);
					}
				}				
			}
			picked_this_round.clear();
		}

		Relabel(graph);
	}


};

class Closeness_Ordering : public Ordering {
public:
	Closeness_Ordering(NodeID k, Graph& graph) { // Number of sample roots.
		
		if (k <= 0) {
			cout << " Should pick at least 1 seed!" << endl;
			return;
		}
		
		inv.resize(numOfVertices, -1);
		rank.resize(numOfVertices, -1);
		
		vector<NodeID> roots(k);

		srand((unsigned)time(NULL));
		vector<NodeID> vlist(numOfVertices);
		for (size_t v = 0; v < numOfVertices; ++ v) vlist[v] = v;
		random_shuffle(vlist.begin(), vlist.end());
		for (size_t v = 0; v < k; ++v) roots[v] = vlist[v];
		vlist.clear();

		vector<EdgeWeight> closeness(numOfVertices);

		// K iterations to calculate the approximated closeness.
		vector<EdgeWeight> distances;
		vector<bool> vis;
		vector<NodeID> que;
		NodeID que_h = 0;
		graph_search::BFS_init(distances, que, que_h, vis);

		for (size_t r = 0; r < k; ++r) {
			NodeID root = roots[r];
			graph_search::BFS(distances, que, que_h, vis, root, graph);
			for (size_t v = 0; v < numOfVertices; ++v)
				if (distances[v] < INF_WEIGHT - 5) 
					closeness[v] += distances[v];
			graph_search::BFS_clear(distances, que, que_h, vis);
			if (DIRECTED_FLAG == true) {
				graph_search::backward_BFS(distances, que, que_h, vis, root, graph);
				for (size_t v = 0; v < numOfVertices; ++v)
					if (distances[v] < INF_WEIGHT - 5)
						closeness[v] += distances[v];
				graph_search::BFS_clear(distances, que, que_h, vis);
			}
		}
		/*for (size_t v = 0; v < numOfVertices; ++v)
			if (closeness_counter[v] != 0)
				closeness[v] = closeness[v] / closeness_counter[v];
			else
				closeness[v] = INF_WEIGHT;*/

		vector<pair<EdgeWeight, NodeID> > close_rank(numOfVertices);
		for (size_t v = 0; v < numOfVertices; ++v)
			close_rank[v] = make_pair(closeness[v], v);

		sort(close_rank.begin(), close_rank.end());

		for (size_t v = 0; v < numOfVertices; ++v) {
			inv[v] = close_rank[v].second;
		}

		Relabel(graph);
	}
	
	Closeness_Ordering(NodeID k, WGraph& wgraph) {

		if (k <= 0) {
			cout << " Should pick at least 1 seed!" << endl;
			return;
		}

		inv.resize(numOfVertices, -1);
		rank.resize(numOfVertices, -1);

		vector<NodeID> roots(k);
		srand((unsigned)time(NULL));
		vector<NodeID> vlist(numOfVertices);
		for (size_t v = 0; v < numOfVertices; ++v) vlist[v] = v;
		random_shuffle(vlist.begin(), vlist.end());
		for (size_t v = 0; v < k; ++v) roots[v] = vlist[v];
		vlist.clear();

		vector<EdgeWeight> closeness(numOfVertices);

		// K iterations to calculate the approximated closeness.
		vector<EdgeWeight> distances;
		vector<bool> vis;
		queue<NodeID> visited_que;
		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);
		graph_search::Dijkstra_init(distances, pqueue, visited_que, vis);

		for (size_t r = 0; r < k; ++r) {
			NodeID root = roots[r];
			graph_search::Dijkstra(distances, pqueue, visited_que, vis, root, wgraph);
			for (size_t v = 0; v < numOfVertices; ++v)
				if (distances[v] < INF_WEIGHT - 5)
					closeness[v] += distances[v];

			graph_search::Dijkstra_clear(distances, pqueue, visited_que, vis);
			if (DIRECTED_FLAG == true) {
				graph_search::backward_Dijkstra(distances, pqueue, visited_que, vis, root, wgraph);
				for (size_t v = 0; v < numOfVertices; ++v)
					closeness[v] += distances[v];
				graph_search::Dijkstra_clear(distances, pqueue, visited_que, vis);
			}
		}

		//for (size_t v = 0; v < numOfVertices; ++v)  closeness[v] = closeness[v] / (double)k;

		vector<pair<EdgeWeight, NodeID> > close_rank(numOfVertices);
		for (size_t v = 0; v <numOfVertices; ++v)
			close_rank[v] = make_pair(closeness[v], v);

		sort(close_rank.begin(), close_rank.end());

		for (size_t v = 0; v < numOfVertices; ++v) {
			inv[v] = close_rank[v].second;
		}

		Relabel(wgraph);
	}	

};

class Betweenness_Ordering :public Ordering {

	typedef	vector<NodeID> large_tree; // A V-sized parent-pointer array representing sampled trees. -1 for those vertices which do not appear in the tree.
	//typedef	unordered_map<NodeID, NodeID> small_tree; // A smaller hashmap, mapping a tree node vertex to its parent vertex.
	typedef	google::dense_hash_map<NodeID, NodeID> small_tree; // A smaller hashmap, mapping a tree node vertex to its parent vertex.
	
public:
	vector<double> iteration_time;
	Label labels;
	DLabel dlabels;
	PLabel plabels;
	DPLabel dplabels;
	CLabel clabels;
	bool SECOND_LEVEL;
	
	long children_size;
	long r_children_size;

	vector<large_tree> ltrees;
	vector<small_tree> strees;
	vector<NodeID> trees_pointers;
	double init_time;
	double updating_time;
	double labeling_time;
	double adding_time;
	double selecting_time;
	long long bounded_resources;
	int num_of_trees;
	long total_resources;
	long long alive_resources;
	double u1_time;
	double u2_time;
	double u3_time = 0;
	double u4_time = 0;
	double u5_time = 0;
	double u6_time = 0;
	double u7_time = 0;
	double u8_time = 0;
	double u9_time = 0;
	double u10_time = 0;
	double u11_time = 0;
	double u12_time = 0;
	double u81_time = 0;
	double u82_time = 0;
	double u91_time = 0;
	double u92_time = 0;

	small_tree empty_small_tree = small_tree();
	bool BUILD_SMALL_TREE_FLAG = false;
	double get_small_parent_time;
	bool switch_small;
	
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

	//using boost::hash_combine
	template <class T>
	inline void hash_combine(std::size_t& seed, T const& v)
	{
		seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
	}
	
	int build_large_tree_bfs(NodeID source, vector<int>& ltree, vector<NodeID>& descendants, vector<NodeID>& descendants_parents, vector<NodeID>& que, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<bool>& usd, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& r_tmp_idx, /*vector<vector<NodeID> >& adj*/ Graph& graph, int accum_new_arcs, int labeled_arcs_bound) {
		NodeID new_tree_arcs = 0;
		
		descendants.clear();
		ltree.resize(numOfVertices, -1);
		//ltree[source] = -1; // As the root, it has no parent.
		descendants_parents[source] = -1;
		
		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
		const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];
		const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[source];
		
		que[que_h++] = source;
		vis[source] = true;
		que_t1 = que_h;
		++new_tree_arcs;

		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
		}

		for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID v = que[que_i];
						
				if (usd[v] == true) continue;

			//	int adjv_size = adj[v].size();

				if (DIRECTED_FLAG == false) {
					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
					int tmp_idx_v_size = tmp_idx_v.first.size();
					for (int i = 0; i < tmp_idx_v_size; ++i) {
						NodeID w = tmp_idx_v.first[i];
						EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
						if (td <= d) {
							// num of hubs touched so far when building a new tree.
							goto pruned;
						}
					}
				}else{
					pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];

					// Pruned by the forward labels of r and backward labels of v in the forward search from r when reaching v.
					int r_tmp_idx_v_size = r_tmp_idx_v.first.size();
					for (int i = 0; i < r_tmp_idx_v_size; ++i) {
						NodeID w = r_tmp_idx_v.first[i];
						EdgeWeight td = r_tmp_idx_v.second[i] + dst_r[w];
						if (td <= d) {
							// num of hubs touched so far when building a new tree.
							goto pruned;
						}
					}
				}				

				descendants.push_back(v);
				++new_tree_arcs;
				if (accum_new_arcs + new_tree_arcs > labeled_arcs_bound) {
					goto stop_growing;
				}

				// Array Representation
		/*		for (int i = 0; i < adjv_size; ++i) {
					NodeID w = adj[v][i];*/
				for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid){
					NodeID w = graph.edges[eid];

					if (!vis[w]) {
						que[que_h++] = w;
						descendants_parents[w] = v;
						//ltree[w] = v; // Set parent nodes. 
						vis[w] = true;		// num of hubs touched so far when building a new tree.					
					}
				}
			pruned:
				{}
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}

		stop_growing:
		{}

		for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
		}
		for (int i = 0; i < descendants.size(); ++i) {
			ltree[descendants[i]] = descendants_parents[descendants[i]];
		}

		return new_tree_arcs;
	}

	int build_small_tree_bfs(NodeID source, small_tree& stree, vector<NodeID>& descendants, vector<NodeID>& descendants_parents, vector<NodeID>& que, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<bool>& usd, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& r_tmp_idx, /*vector<vector<NodeID> >& adj*/ Graph& graph, int accum_new_arcs, int labeled_arcs_bound) {
		NodeID new_tree_arcs = 0;
		
		descendants.clear();
		stree.clear();
	//	stree.reserve(V / 8);
		descendants_parents[source] = -1;
		//stree.insert(make_pair(source, -1));// As the root, it has no parent.

		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
		const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];
		const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[source];
		
		que[que_h++] = source;
		vis[source] = true;
		que_t1 = que_h;
		++new_tree_arcs;

		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
		}

		for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID v = que[que_i];

				if (usd[v] == true) continue;

				//int adjv_size = adj[v].size();
				if (DIRECTED_FLAG == false) {
					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
					int tmp_idx_v_size = tmp_idx_v.first.size();
					for (int i = 0; i < tmp_idx_v_size; ++i) {
						NodeID w = tmp_idx_v.first[i];
						EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
						if (td <= d) {
							// num of hubs touched so far when building a new tree.
							goto pruned;
						}
					}
				}
				else {
					pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];

					// Pruned by the forward labels of r and backward labels of v in the forward search from r when reaching v.
					int r_tmp_idx_v_size = r_tmp_idx_v.first.size();
					for (int i = 0; i < r_tmp_idx_v_size; ++i) {
						NodeID w = r_tmp_idx_v.first[i];
						EdgeWeight td = r_tmp_idx_v.second[i] + dst_r[w];
						if (td <= d) {
							// num of hubs touched so far when building a new tree.
							goto pruned;
						}
					}
				}

				descendants.push_back(v);
				// num of hubs touched so far when building a new tree.
				++new_tree_arcs;
				if (accum_new_arcs + new_tree_arcs > labeled_arcs_bound) {
					goto stop_growing;
				}

				// Array Representation
				/*for (int i = 0; i < adjv_size; ++i) {
					NodeID w = adj[v][i];*/
				for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid){
					NodeID w = graph.edges[eid];

					if (!vis[w]) {
						que[que_h++] = w;
						descendants_parents[w] = v;
						//stree.insert(make_pair(w, v));// Set parent nodes. 
						vis[w] = true;					
					}
				}
			pruned:
				{}
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}

	stop_growing:
		{}

		for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
		}

		for (int i = 0; i < descendants.size(); ++i) {
			stree.insert(make_pair(descendants[i], descendants_parents[descendants[i]]));
		}

		return new_tree_arcs;
	}
	
	int build_large_tree_dij(NodeID source, large_tree& ltree, vector<NodeID>& descendants, vector<NodeID>& descendants_parents, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, queue<NodeID> &visited_que, vector<EdgeWeight>& distances, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<bool>& usd, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& r_tmp_idx, WGraph& wgraph, int accum_new_arcs, int labeled_arcs_bound) {
		NodeID new_tree_arcs = 0;

		descendants.clear();
		ltree.resize(numOfVertices, -1);
	//	ltree[source] = -1; // As the root, it has no parent.
		descendants_parents[source] = -1;

		const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];
		const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[source];

		/*vector<vector<NodeID> >& adj = wgraph.adj;
		vector<vector<EdgeWeight> >& adj_w = wgraph.adj_weight;*/


		pqueue.update(source, 0);
		//++new_tree_arcs;

		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
		}

		visited_que.push(source);
		while(!pqueue.empty()){
				NodeID v;
				EdgeWeight v_d;
				pqueue.extract_min(v, v_d);
				vis[v] = true;


				if (usd[v] == true) continue;
				//int adj_v_size = adj[v].size();

				if (DIRECTED_FLAG == false) {
					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
					for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
						NodeID w = tmp_idx_v.first[i];
						EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
						if (td <= v_d) {
							goto pruned;
						}
					}
				}
				else {
					pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];

					// Pruned by the forward labels of r and backward labels of v in the forward search from r when reaching v.
					for (size_t i = 0; i < r_tmp_idx_v.first.size(); ++i) {
						NodeID w = r_tmp_idx_v.first[i];
						EdgeWeight td = r_tmp_idx_v.second[i] + dst_r[w];

						if (td <= v_d) {
							goto pruned;
						}
					}
				}

				// num of hubs touched so far when building a new tree.
				++new_tree_arcs;
				descendants.push_back(v);
				
				if (accum_new_arcs + new_tree_arcs > labeled_arcs_bound) {
					break;
				}

				/*for (size_t i = 0; i < adj_v_size; ++i) {
					NodeID w = adj[v][i];
					EdgeWeight w_d = adj_w[v][i] + v_d;		*/

				 //Array Representations.
				for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid){
					
					NodeID w = wgraph.edges[eid].first;
					EdgeWeight w_d = wgraph.edges[eid].second + v_d;

					//if (wgraph.adj[v][eid - wgraph.vertices[v]] != w || wgraph.adj_weight[v][eid - wgraph.vertices[v]] + v_d != w_d) cout << "we got problems build large tree." << endl;

					if (!vis[w]) {
						if (distances[w] > w_d) {
							pqueue.update(w, w_d);
							distances[w] = w_d;
							//ltree[w] = v; // Set parent nodes. 
							descendants_parents[w] = v;
							visited_que.push(w);
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
	/*	while (!pqueue.empty()) {
			NodeID v;
			EdgeWeight v_d;
			pqueue.extract_min(v, v_d);
			vis[v] = false;
			distances[v] = INF_WEIGHT;
		}*/

		pqueue.clear_n();

		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
		}

		for (int i = 0; i < descendants.size(); ++i)
			ltree[descendants[i]] = descendants_parents[descendants[i]];

		return new_tree_arcs;
	}

	int build_small_tree_dij(NodeID source, small_tree& stree, vector<NodeID>& descendants, vector<NodeID>& descendants_parents, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, queue<NodeID> &visited_que, vector<EdgeWeight>& distances, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<bool>& usd, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& r_tmp_idx, WGraph& wgraph, int accum_new_arcs, int labeled_arcs_bound) {
		NodeID new_tree_arcs = 0;

		descendants.clear();
		stree.clear();
		//	ltree[source] = -1; // As the root, it has no parent.
		descendants_parents[source] = -1;

		const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];
		const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[source];

		/*vector<vector<NodeID> >& adj = wgraph.adj;
		vector<vector<EdgeWeight> >& adj_w = wgraph.adj_weight;*/


		pqueue.update(source, 0);
		//++new_tree_arcs;

		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
		}
		visited_que.push(source);
		while (!pqueue.empty()) {
			NodeID v;
			EdgeWeight v_d;
			pqueue.extract_min(v, v_d);
			vis[v] = true;
			//visited_que.push(v);
			


			if (usd[v] == true) continue;

			//int adj_v_size = adj[v].size();

			if (DIRECTED_FLAG == false) {
				pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
				for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
					NodeID w = tmp_idx_v.first[i];
					EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
					if (td <= v_d) {
						goto pruned;
					}
				}
			}
			else {
				pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];

				// Pruned by the forward labels of r and backward labels of v in the forward search from r when reaching v.
				for (size_t i = 0; i < r_tmp_idx_v.first.size(); ++i) {
					NodeID w = r_tmp_idx_v.first[i];
					EdgeWeight td = r_tmp_idx_v.second[i] + dst_r[w];

					if (td <= v_d) {
						goto pruned;
					}
				}
			}

			// num of hubs touched so far when building a new tree.
			++new_tree_arcs;
			descendants.push_back(v);

			if (accum_new_arcs + new_tree_arcs > labeled_arcs_bound) {
				break;
			}
			
		/*	for (size_t i = 0; i < adj_v_size; ++i) {
				NodeID w = adj[v][i];
				EdgeWeight w_d = adj_w[v][i] + v_d;*/

			//Array Representation
			for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid){
				NodeID w = wgraph.edges[eid].first;
				EdgeWeight w_d = wgraph.edges[eid].second + v_d;


				//if (wgraph.adj[v][eid - wgraph.vertices[v]] != w || wgraph.adj_weight[v][eid - wgraph.vertices[v]] + v_d != w_d ) cout << "we got problems build small tree." << endl;

				if (!vis[w]) {
					if (distances[w] > w_d) {
						pqueue.update(w, w_d);
						distances[w] = w_d;
						//ltree[w] = v; // Set parent nodes. 
						descendants_parents[w] = v;
						visited_que.push(w);
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
	
		//while (!pqueue.empty()) {
		//	NodeID v;
		//	EdgeWeight v_d;
		//	pqueue.extract_min(v, v_d);
		//	//vis[v] = false;
		//	//distances[v] = INF_WEIGHT;
		//}
		
		pqueue.clear_n();
	//	pqueue.clear();

		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
		}

		for (int i = 0; i < descendants.size(); ++i)
			stree.insert(make_pair(descendants[i], descendants_parents[descendants[i]]));

		return new_tree_arcs;
	}
	
	void get_descendants(int tree_id, NodeID picked_v, vector<NodeID>& descendants, vector<NodeID>& que, vector<bool>& vis, vector<vector<NodeID> >& adj) {
		descendants.clear();

		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;

		que[que_h++] = picked_v;
		vis[picked_v] = true;
		que_t1 = que_h;


		bool isLarge = false;
		int stree_idx = trees_pointers[tree_id];
		if (stree_idx == -1) 
			isLarge = true;
		
		small_tree& stree_tid = (stree_idx  != -1 ) ? strees[stree_idx] : empty_small_tree;
		large_tree& ltree_tid = ltrees[tree_id];


		for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID v = que[que_i];
				descendants.push_back(v);
				vector<EdgeWeight>& adj_v = adj[v];
				NodeID adj_v_size = adj_v.size();
				for (size_t i = 0; i <adj_v_size; ++i) {
					NodeID w = adj_v[i];

					NodeID parentID = 0;
					if (isLarge)
						parentID = get_large_parent(w, ltree_tid);
					else
						parentID = get_small_parent(w, stree_tid);

					if (!vis[w] && parentID == v) {
						que[que_h++] = w;
						vis[w] = true;
					}
				}
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}

		for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
	}

	void get_descendants(int tree_id, NodeID picked_v, vector<NodeID>& descendants, vector<NodeID>& descendants_parents, vector<NodeID>& que, vector<bool>& vis, vector<vector<NodeID> >& adj, const bool isLarge, small_tree& stree_tid, large_tree& ltree_tid) {
	
		descendants.clear();
	//	descendants.reserve(V / 8);

		descendants.push_back(picked_v);
		descendants_parents[picked_v] = -1;

		vis[picked_v] = true;

	/*	bool isLarge = false;
		int stree_idx = trees_pointers[tree_id];
		if (stree_idx == -1)
			isLarge = true;

		small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
		large_tree& ltree_tid = ltrees[tree_id];*/


		for (int i = 0; i < descendants.size(); ++i) {
			NodeID v = descendants[i];
			vector<EdgeWeight>& adj_v = adj[v];
			NodeID adj_v_size = adj_v.size();
			for (int j = 0; j < adj_v_size; ++j) {
				NodeID w = adj_v[j];

				NodeID parentID = 0;


				if (isLarge)
					parentID = get_large_parent(w, ltree_tid);
				else {
				//	double get_small_parent_incre_time = GetCurrentTimeSec();
					parentID = get_small_parent(w, stree_tid);
					//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
				}


				if (!vis[w] && parentID == v) {
					descendants.push_back(w);
					descendants_parents[w] = v;
					vis[w] = true;
				}
			}
		}
		int dsize = descendants.size();
		for (int i = 0; i <dsize; ++i) vis[descendants[i]] = false;
	}
	
	void get_descendants(int tree_id, NodeID picked_v, vector<NodeID>& descendants, vector<NodeID>& descendants_parents, vector<NodeID>& que, vector<bool>& vis, Graph& graph, const bool isLarge, small_tree& stree_tid, large_tree& ltree_tid) {

		descendants.clear();
		//	descendants.reserve(V / 8);

		descendants.push_back(picked_v);
		descendants_parents[picked_v] = -1;

		vis[picked_v] = true;

		/*	bool isLarge = false;
		int stree_idx = trees_pointers[tree_id];
		if (stree_idx == -1)
		isLarge = true;

		small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
		large_tree& ltree_tid = ltrees[tree_id];*/


		for (int i = 0; i < descendants.size(); ++i) {
			NodeID v = descendants[i];

			//vector<EdgeWeight>& adj_v = adj[v];
			//NodeID adj_v_size = adj_v.size();

			//for (int j = 0; j < adj_v_size; ++j) {
			//NodeID w = adj_v[j];

			//Array Representation

			for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid) {

				NodeID w = graph.edges[eid];

				NodeID parentID = 0;


				if (isLarge)
					parentID = get_large_parent(w, ltree_tid);
				else {
					//	double get_small_parent_incre_time = GetCurrentTimeSec();
					parentID = get_small_parent(w, stree_tid);
					//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
				}


				if (!vis[w] && parentID == v) {
					descendants.push_back(w);
					descendants_parents[w] = v;
					vis[w] = true;
				}
			}
		}
		int dsize = descendants.size();
		for (int i = 0; i <dsize; ++i) vis[descendants[i]] = false;
	}


	void get_descendants(int tree_id, NodeID picked_v, vector<NodeID>& descendants, vector<NodeID>& descendants_parents, vector<NodeID>& que, vector<bool>& vis, WGraph& wgraph, const bool isLarge, small_tree& stree_tid, large_tree& ltree_tid) {

		descendants.clear();
		//	descendants.reserve(V / 8);

		descendants.push_back(picked_v);
		descendants_parents[picked_v] = -1;

		vis[picked_v] = true;

		/*	bool isLarge = false;
		int stree_idx = trees_pointers[tree_id];
		if (stree_idx == -1)
		isLarge = true;

		small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
		large_tree& ltree_tid = ltrees[tree_id];*/


		for (int i = 0; i < descendants.size(); ++i) {
			NodeID v = descendants[i];

		/*	vector<EdgeWeight>& adj_v = wgraph.adj[v];
			NodeID adj_v_size = adj_v.size();
			
			for (int j = 0; j < adj_v_size; ++j) {
				NodeID w = adj_v[j];*/

			//Array Representation

			for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1];++eid){
				
				NodeID w = wgraph.edges[eid].first;

				//if(wgraph.adj[v][eid - wgraph.vertices[v]] != w) cout << "we got problems. get descent" << endl;

				NodeID parentID = 0;


				if (isLarge)
					parentID = get_large_parent(w, ltree_tid);
				else {
					//	double get_small_parent_incre_time = GetCurrentTimeSec();
					parentID = get_small_parent(w, stree_tid);
					//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
				}


				if (!vis[w] && parentID == v) {
					descendants.push_back(w);
					descendants_parents[w] = v;
					vis[w] = true;
				}
			}
		}
		int dsize = descendants.size();
		for (int i = 0; i <dsize; ++i) vis[descendants[i]] = false;
	}

	void get_ascendants(int tree_id, NodeID picked_v, vector<NodeID>& ascendants, vector<NodeID>& que, vector<bool>& vis, vector<vector<NodeID> >& adj, vector<vector<NodeID> >& r_adj) {
		ascendants.clear();

		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
		//	vector<vector<NodeID> >& adj = graph.adj;
		//vector<vector<NodeID> >& r_adj = graph.r_adj;

		que[que_h++] = picked_v;
		vis[picked_v] = true;
		que_t1 = que_h;


		bool isLarge = false;
		int stree_idx = trees_pointers[tree_id];
		if (stree_idx == -1)
			isLarge = true;

		small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
		large_tree& ltree_tid = ltrees[tree_id];

		for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID v = que[que_i];
				ascendants.push_back(v);

				NodeID parentID = 0;
				if (isLarge)
					parentID = get_large_parent(v, ltree_tid);
				else
					parentID = get_small_parent(v, stree_tid);

				if (DIRECTED_FLAG == false) {
					vector<EdgeWeight>& adj_v = adj[v];
					int adj_v_size = adj_v.size();
					for (int i = 0; i < adj_v_size; ++i) {
						NodeID w = adj_v[i];


						//if (!vis[w] && get_parent(tree_id, v) == w) {
						if (!vis[w] && parentID == w) {
							que[que_h++] = w;
							vis[w] = true;
						}
					}
				}
				else {
					vector<EdgeWeight>& r_adj_v = r_adj[v];
					int r_adj_v_size = r_adj_v.size();


					for (int i = 0; i < r_adj_v_size; ++i) {
						NodeID w = r_adj_v[i];

						//	if (!vis[w] && get_parent(tree_id, v) == w) {
						if (!vis[w] && parentID == w) {
							que[que_h++] = w;
							vis[w] = true;
						}
					}
				}
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}

		for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
	}

	void get_ascendants(int tree_id, NodeID picked_v, vector<NodeID>& ascendants, vector<NodeID>& que, vector<bool>& vis, vector<vector<NodeID> >& adj, vector<vector<NodeID> >& r_adj, const bool isLarge, small_tree& stree_tid, large_tree& ltree_tid) {
		ascendants.clear();

		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
	//	vector<vector<NodeID> >& adj = graph.adj;
		//vector<vector<NodeID> >& r_adj = graph.r_adj;

		que[que_h++] = picked_v;
		vis[picked_v] = true;
		que_t1 = que_h;


		//bool isLarge = false;
		//int stree_idx = trees_pointers[tree_id];
		//if (stree_idx == -1) 
		//	isLarge = true;
		//
		//small_tree& stree_tid = (stree_idx != -1 )? strees[stree_idx] : empty_small_tree;
		//large_tree& ltree_tid = ltrees[tree_id];

		for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID v = que[que_i];
				ascendants.push_back(v);

				NodeID parentID = 0;
				if (isLarge)
					parentID = get_large_parent(v, ltree_tid);
				else {
					//double get_small_parent_incre_time = GetCurrentTimeSec();
					parentID = get_small_parent(v, stree_tid);
					//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
				}
				if (DIRECTED_FLAG == false) {
					vector<EdgeWeight>& adj_v = adj[v];
					int adj_v_size = adj_v.size();
					for (int i = 0; i < adj_v_size; ++i) {
						NodeID w = adj_v[i];

						//if (!vis[w] && get_parent(tree_id, v) == w) {
						if (!vis[w] &&parentID == w) {
							que[que_h++] = w;
							vis[w] = true;
						}
					}
				}else {
					vector<EdgeWeight>& r_adj_v = r_adj[v];
					int r_adj_v_size = r_adj_v.size();			

					for (int i = 0; i < r_adj_v_size; ++i) {
						NodeID w = r_adj_v[i];
						
						if (!vis[w] && parentID == w) {
							que[que_h++] = w;
							vis[w] = true;
						}
					}
				}
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}

		for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
	}

	void get_ascendants(int tree_id, NodeID picked_v, vector<NodeID>& ascendants, vector<NodeID>& que, vector<bool>& vis, Graph& graph, const bool isLarge, small_tree& stree_tid, large_tree& ltree_tid) {
		ascendants.clear();

		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
		//	vector<vector<NodeID> >& adj = graph.adj;
		//vector<vector<NodeID> >& r_adj = graph.r_adj;

		que[que_h++] = picked_v;
		vis[picked_v] = true;
		que_t1 = que_h;


		//bool isLarge = false;
		//int stree_idx = trees_pointers[tree_id];
		//if (stree_idx == -1) 
		//	isLarge = true;
		//
		//small_tree& stree_tid = (stree_idx != -1 )? strees[stree_idx] : empty_small_tree;
		//large_tree& ltree_tid = ltrees[tree_id];

		for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID v = que[que_i];
				ascendants.push_back(v);

				NodeID parentID = 0;
				if (isLarge)
					parentID = get_large_parent(v, ltree_tid);
				else {
					//double get_small_parent_incre_time = GetCurrentTimeSec();
					parentID = get_small_parent(v, stree_tid);
					//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
				}
				if (DIRECTED_FLAG == false) {
					//vector<EdgeWeight>& adj_v = adj[v];
					//	int adj_v_size = adj_v.size();
					//	for (int i = 0; i < adj_v_size; ++i) {
					//		NodeID w = adj_v[i];

					//Array Representation
					for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid) {

						NodeID w = graph.edges[eid];

						//if (!vis[w] && get_parent(tree_id, v) == w) {
						if (!vis[w] && parentID == w) {
							que[que_h++] = w;
							vis[w] = true;
						}
					}
				}
				else {
					//vector<EdgeWeight>& r_adj_v = r_adj[v];
					//	int r_adj_v_size = r_adj_v.size();
					//for (int i = 0; i < r_adj_v_size; ++i) {
					//	NodeID w = r_adj_v[i];

					//Array Representation
					for (EdgeID eid = graph.r_vertices[v]; eid < graph.r_vertices[v + 1]; ++eid) {

						NodeID w = graph.r_edges[eid];

						if (!vis[w] && parentID == w) {
							que[que_h++] = w;
							vis[w] = true;
						}
					}
				}
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}

		for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
	}

	void get_ascendants(int tree_id, NodeID picked_v, vector<NodeID>& ascendants, vector<NodeID>& que, vector<bool>& vis, WGraph& wgraph, const bool isLarge, small_tree& stree_tid, large_tree& ltree_tid) {
		ascendants.clear();

		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
		//	vector<vector<NodeID> >& adj = graph.adj;
		//vector<vector<NodeID> >& r_adj = graph.r_adj;

		que[que_h++] = picked_v;
		vis[picked_v] = true;
		que_t1 = que_h;


		//bool isLarge = false;
		//int stree_idx = trees_pointers[tree_id];
		//if (stree_idx == -1) 
		//	isLarge = true;
		//
		//small_tree& stree_tid = (stree_idx != -1 )? strees[stree_idx] : empty_small_tree;
		//large_tree& ltree_tid = ltrees[tree_id];

		for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID v = que[que_i];
				ascendants.push_back(v);

				NodeID parentID = 0;
				if (isLarge)
					parentID = get_large_parent(v, ltree_tid);
				else {
					//double get_small_parent_incre_time = GetCurrentTimeSec();
					parentID = get_small_parent(v, stree_tid);
					//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
				}
				if (DIRECTED_FLAG == false) {
				/*	vector<EdgeWeight>& adj_v = wgraph.adj[v];
					int adj_v_size = adj_v.size();
					for (int i = 0; i < adj_v_size; ++i) {
						NodeID w = adj_v[i];*/

					//Array Representation
					for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid){

						NodeID w = wgraph.edges[eid].first;

						//if (!vis[w] && get_parent(tree_id, v) == w) {
						if (!vis[w] && parentID == w) {
							que[que_h++] = w;
							vis[w] = true;
						}
					}
				}
				else {
				/*	vector<EdgeWeight>& r_adj_v = wgraph.r_adj[v];
					int r_adj_v_size = r_adj_v.size();
					for (int i = 0; i < r_adj_v_size; ++i) {
						NodeID w = r_adj_v[i];*/

					//Array Representation
					for (EdgeID eid = wgraph.r_vertices[v]; eid < wgraph.r_vertices[v + 1]; ++eid){

						NodeID w = wgraph.r_edges[eid].first;

					//	if (wgraph.r_adj[v][eid - wgraph.r_vertices[v]] != w) cout << "we got problems get ascendants." << endl;

						if (!vis[w] && parentID == w) {
							que[que_h++] = w;
							vis[w] = true;
						}
					}
				}
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}

		for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
	}

	int labeling_source_bfs(NodeID source, vector<NodeID>& que, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<EdgeWeight>& r_dst_r, vector<bool>& usd, NodeID ranking, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& r_tmp_idx, vector<pair<vector<NodeID>, vector<NodeID> > >& tmp_idx_token_parents, vector<pair<vector<NodeID>, vector<NodeID> > >& r_tmp_idx_token_parents, Graph& graph) {
		NodeID labeled_arcs = 0;
		NodeID visited_arcs = 0;

		//vector<vector<NodeID> >& adj = graph.adj;
		const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];
		const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[source];
			
		if(DIRECTED_FLAG == true){
			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i) {
				r_dst_r[r_tmp_idx_r.first[i]] = r_tmp_idx_r.second[i];
			}
		}
			
		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
		}
		
		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;

		que[que_h++] = source;
		vis[source] = true;
		que_t1 = que_h;

		for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID v = que[que_i];
				// hubs touched so far when building this labels
				if (usd[v]) continue; 


				if (DIRECTED_FLAG == false) {
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

					labeled_arcs++;
					tmp_idx_v.first.back() = ranking;
					tmp_idx_v.second.back() = d;
					tmp_idx_v.first.push_back(numOfVertices);
					tmp_idx_v.second.push_back(INF_WEIGHT);
					
					tmp_idx_token_parents_v.first.back() = numOfVertices;
					tmp_idx_token_parents_v.second.back() = INF_WEIGHT;
					tmp_idx_token_parents_v.first.push_back(numOfVertices);
					tmp_idx_token_parents_v.second.push_back(INF_WEIGHT);

					
				}
				else {
					pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];
					pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_token_parents_v = r_tmp_idx_token_parents[v];

					// Pruned by the forward labels of r and backward labels of v in the forward search from r when reaching v.
					for (size_t i = 0; i < r_tmp_idx_v.first.size(); ++i) {
						NodeID w = r_tmp_idx_v.first[i];
						EdgeWeight td = r_tmp_idx_v.second[i] + dst_r[w];
						if(r_tmp_idx_v.second[i] == d + r_dst_r[w] && r_tmp_idx_token_parents_v.first[i] == numOfVertices){
							r_tmp_idx_token_parents_v.first[i] = ranking;//tmp_idx_token_parents_v.first.size();
							r_tmp_idx_token_parents_v.second[i] = r_dst_r[w];
						}
						if (td <= d) {
							//hubs touched so far when building this labels
							goto pruned_forward;
						}
					}
					

					// arcs touched so far when building this labels
					labeled_arcs++;
					r_tmp_idx_v.first.back() = ranking;
					r_tmp_idx_v.second.back() = d;
					r_tmp_idx_v.first.push_back(numOfVertices);
					r_tmp_idx_v.second.push_back(INF_WEIGHT);
					
					r_tmp_idx_token_parents_v.first.back() = numOfVertices;
					r_tmp_idx_token_parents_v.second.back() = INF_WEIGHT;
					r_tmp_idx_token_parents_v.first.push_back(numOfVertices);
					r_tmp_idx_token_parents_v.second.push_back(INF_WEIGHT);
					
				}

				// Array Representation
				/*for (size_t i = 0; i < adj[v].size(); ++i) {
					NodeID w = adj[v][i];*/
					
				visited_arcs += (graph.vertices[v + 1] - graph.vertices[v]);
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


		// Backward search.
		if (DIRECTED_FLAG == true) {


			//vector<vector<NodeID> >& r_adj = graph.r_adj;

			que_t0 = 0, que_t1 = 0, que_h = 0;

			que[que_h++] = source;
			vis[source] = true;
			que_t1 = que_h;

			for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID v = que[que_i];

					if (usd[v]) continue;

					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_token_parents_v = tmp_idx_token_parents[v];

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

					labeled_arcs++;
					tmp_idx_v.first.back() = ranking;
					tmp_idx_v.second.back() = d;
					tmp_idx_v.first.push_back(numOfVertices);
					tmp_idx_v.second.push_back(INF_WEIGHT);
					
					tmp_idx_token_parents_v.first.back() = numOfVertices;
					tmp_idx_token_parents_v.second.back() = INF_WEIGHT;
					tmp_idx_token_parents_v.first.push_back(numOfVertices);
					tmp_idx_token_parents_v.second.push_back(INF_WEIGHT);
					
				//	for (size_t i = 0; i < r_adj[v].size(); ++i) {
					//	NodeID w = r_adj[v][i];
					
					visited_arcs += (graph.r_vertices[v + 1] - graph.r_vertices[v]);
					for (EdgeID eid = graph.r_vertices[v]; eid < graph.r_vertices[v + 1]; ++eid){
						NodeID w = graph.r_edges[eid];
						if (!vis[w]) {				
							// hubs touched so far when building this labels							
							que[que_h++] = w;
							vis[w] = true;
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
		}
		
		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) 
			dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;		
		//cout << visited_arcs << endl;
		
		usd[source] = true;
		return labeled_arcs;
	}

	int labeling_source_bfs(NodeID source, vector<NodeID>& que, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<bool>& usd, NodeID ranking, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& r_tmp_idx, Graph& graph) {
		NodeID labeled_arcs = 0;
		NodeID visited_arcs = 0;

		//vector<vector<NodeID> >& adj = graph.adj;
		const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];
		const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[source];
				
		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
		}
		
		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;

		que[que_h++] = source;
		vis[source] = true;
		que_t1 = que_h;

		for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID v = que[que_i];
				// hubs touched so far when building this labels
				if (usd[v]) continue; 


				if (DIRECTED_FLAG == false) {
					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
					for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
						NodeID w = tmp_idx_v.first[i];
						EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
						if (td <= d) {
							goto pruned_forward;
						}
					}

					labeled_arcs++;
					tmp_idx_v.first.back() = ranking;
					tmp_idx_v.second.back() = d;
					tmp_idx_v.first.push_back(numOfVertices);
					tmp_idx_v.second.push_back(INF_WEIGHT);
				}
				else {
					pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];

					// Pruned by the forward labels of r and backward labels of v in the forward search from r when reaching v.
					for (size_t i = 0; i < r_tmp_idx_v.first.size(); ++i) {
						NodeID w = r_tmp_idx_v.first[i];
						EdgeWeight td = r_tmp_idx_v.second[i] + dst_r[w];
						if (td <= d) {
							//hubs touched so far when building this labels
							goto pruned_forward;
						}
					}
					

					// arcs touched so far when building this labels
					labeled_arcs++;
					r_tmp_idx_v.first.back() = ranking;
					r_tmp_idx_v.second.back() = d;
					r_tmp_idx_v.first.push_back(numOfVertices);
					r_tmp_idx_v.second.push_back(INF_WEIGHT);
				}

				// Array Representation
				/*for (size_t i = 0; i < adj[v].size(); ++i) {
					NodeID w = adj[v][i];*/
					
				visited_arcs += (graph.vertices[v + 1] - graph.vertices[v]);
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
		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
		}
		

		// Backward search.
		if (DIRECTED_FLAG == true) {

			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i) {
				dst_r[r_tmp_idx_r.first[i]] = r_tmp_idx_r.second[i];
			}

			//vector<vector<NodeID> >& r_adj = graph.r_adj;

			que_t0 = 0, que_t1 = 0, que_h = 0;

			que[que_h++] = source;
			vis[source] = true;
			que_t1 = que_h;

			for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID v = que[que_i];

					if (usd[v]) continue;

					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];

					for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
						NodeID w = tmp_idx_v.first[i];
						EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
						if (td <= d) {
							goto pruned_backward;
						}
					}

					labeled_arcs++;
					tmp_idx_v.first.back() = ranking;
					tmp_idx_v.second.back() = d;
					tmp_idx_v.first.push_back(numOfVertices);
					tmp_idx_v.second.push_back(INF_WEIGHT);

				//	for (size_t i = 0; i < r_adj[v].size(); ++i) {
					//	NodeID w = r_adj[v][i];
					
					visited_arcs += (graph.r_vertices[v + 1] - graph.r_vertices[v]);
					for (EdgeID eid = graph.r_vertices[v]; eid < graph.r_vertices[v + 1]; ++eid){
						NodeID w = graph.r_edges[eid];
						if (!vis[w]) {				
							// hubs touched so far when building this labels							
							que[que_h++] = w;
							vis[w] = true;
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
				dst_r[r_tmp_idx_r.first[i]] = INF_WEIGHT;
		}
		
		//cout << visited_arcs << endl;
		
		usd[source] = true;
		return labeled_arcs;
	}

	int labeling_source_dij(NodeID source, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, queue<NodeID> &visited_que, vector<EdgeWeight>& distances, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<EdgeWeight>& r_dst_r, vector<bool>& usd, NodeID ranking, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& r_tmp_idx, vector<pair<vector<NodeID>, vector<NodeID> > >& tmp_idx_token_parents, vector<pair<vector<NodeID>, vector<NodeID> > >& r_tmp_idx_token_parents, WGraph& wgraph) {
		NodeID labeled_arcs = 0;

		NodeID visited_arcs = 0;
	/*	vector<vector<NodeID> >& adj = wgraph.adj;
		vector<vector<EdgeWeight> >& adj_w = wgraph.adj_weight;*/

		const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];
		const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[source];

		
		if(DIRECTED_FLAG == true){
			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i) {
				r_dst_r[r_tmp_idx_r.first[i]] = r_tmp_idx_r.second[i];
			}
		}
		
		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
		}

		pqueue.update(source, 0);


		while (!pqueue.empty()) {
			NodeID v;
			EdgeWeight v_d;
			pqueue.extract_min(v, v_d);
			vis[v] = true;
			visited_que.push(v);

			if (usd[v]) continue;


			if (DIRECTED_FLAG == false) {
				pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
				pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_token_parents_v = tmp_idx_token_parents[v];

				
				for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
					NodeID w = tmp_idx_v.first[i];
					EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
					if(tmp_idx_v.second[i] == v_d + dst_r[w] && tmp_idx_token_parents_v.first[i] == numOfVertices){
							tmp_idx_token_parents_v.first[i] = ranking;//tmp_idx_token_parents_v.first.size();
							tmp_idx_token_parents_v.second[i] = dst_r[w];
						}
					if (td <= v_d) {
						goto pruned_forward;
					}
				}

				// arcs touched so far when building this labels
				labeled_arcs++;

				tmp_idx_v.first.back() = ranking;
				tmp_idx_v.second.back() = v_d;
				tmp_idx_v.first.push_back(numOfVertices);
				tmp_idx_v.second.push_back(INF_WEIGHT);
				
				tmp_idx_token_parents_v.first.back() = numOfVertices;
				tmp_idx_token_parents_v.second.back() = INF_WEIGHT;
				tmp_idx_token_parents_v.first.push_back(numOfVertices);
				tmp_idx_token_parents_v.second.push_back(INF_WEIGHT);

			}
			else {
				pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];
				pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_token_parents_v = r_tmp_idx_token_parents[v];

				// Pruned by the forward labels of r and backward labels of v in the forward search from r when reaching v.
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

				// arcs touched so far when building this labels
				labeled_arcs++;
				
				r_tmp_idx_v.first.back() = ranking;
				r_tmp_idx_v.second.back() = v_d;
				r_tmp_idx_v.first.push_back(numOfVertices);
				r_tmp_idx_v.second.push_back(INF_WEIGHT);
				
				r_tmp_idx_token_parents_v.first.back() = numOfVertices;
				r_tmp_idx_token_parents_v.second.back() = INF_WEIGHT;
				r_tmp_idx_token_parents_v.first.push_back(numOfVertices);
				r_tmp_idx_token_parents_v.second.push_back(INF_WEIGHT);
				
			}

	/*		for (size_t i = 0; i < adj[v].size(); ++i) {
			
				NodeID w = adj[v][i];
				EdgeWeight w_d = adj_w[v][i] + v_d;*/

			//Array Representation
			
			visited_arcs += (wgraph.vertices[v + 1] - wgraph.vertices[v]);
			
			for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid){

				NodeID w = wgraph.edges[eid].first;
				EdgeWeight w_d = wgraph.edges[eid].second + v_d;

		//		if (wgraph.adj[v][eid - wgraph.vertices[v]] != w || wgraph.adj_weight[v][eid - wgraph.vertices[v]] + v_d != w_d) cout << "we got problems labeling tree." << endl;
			
				if (!vis[w]) {
					if (distances[w] > w_d) {
						pqueue.update(w, w_d);
						distances[w] = w_d;
					}
				}
			}
		pruned_forward:
			{}
		}
		
		// Clear foward search structures.
		while (!visited_que.empty()) {
			NodeID vis_v = visited_que.front();
			visited_que.pop();
			vis[vis_v] = false;
			distances[vis_v] = INF_WEIGHT;
		}


		// Backward search.
		if (DIRECTED_FLAG == true) {

	/*		vector<vector<NodeID> >& r_adj = wgraph.r_adj;
			vector<vector<EdgeWeight> >& r_adj_w = wgraph.r_adj_weight;*/			


			pqueue.update(source, 0);
			while (!pqueue.empty()) {
				NodeID v;
				EdgeWeight v_d;
				pqueue.extract_min(v, v_d);
				vis[v] = true;
				visited_que.push(v);


				if (usd[v]) continue;


				pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
				pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_token_parents_v = tmp_idx_token_parents[v];
				
				for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
					NodeID w = tmp_idx_v.first[i];
					EdgeWeight td = tmp_idx_v.second[i] + r_dst_r[w];
					if(tmp_idx_v.second[i] == v_d + dst_r[w] && tmp_idx_token_parents_v.first[i] == numOfVertices){
							tmp_idx_token_parents_v.first[i] = ranking;//tmp_idx_token_parents_v.first.size();
							tmp_idx_token_parents_v.second[i] = dst_r[w];
					}	
					if (td <= v_d) {
						goto pruned_backward;
					}
				}

				// arcs touched so far when building this labels
				labeled_arcs++;
				tmp_idx_v.first.back() = ranking;
				tmp_idx_v.second.back() = v_d;
				tmp_idx_v.first.push_back(numOfVertices);
				tmp_idx_v.second.push_back(INF_WEIGHT);

				tmp_idx_token_parents_v.first.back() = numOfVertices;
				tmp_idx_token_parents_v.second.back() = INF_WEIGHT;
				tmp_idx_token_parents_v.first.push_back(numOfVertices);
				tmp_idx_token_parents_v.second.push_back(INF_WEIGHT);
				
				
		/*		for (size_t i = 0; i < r_adj[v].size(); ++i) {
					NodeID w = r_adj[v][i];
					EdgeWeight w_d = r_adj_w[v][i] + v_d;*/
				
				////Array Representation
				
				visited_arcs += (wgraph.r_vertices[v + 1] - wgraph.r_vertices[v]);
				
				for (EdgeID eid = wgraph.r_vertices[v]; eid < wgraph.r_vertices[v + 1]; ++eid){

					NodeID w = wgraph.r_edges[eid].first;
					EdgeWeight w_d = wgraph.r_edges[eid].second + v_d ;

					if (!vis[w]) {
						if (distances[w] > w_d) {
							pqueue.update(w, w_d);
							distances[w] = w_d;
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
				r_dst_r[r_tmp_idx_r.first[i]] = INF_WEIGHT;		
		}
	
		for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
			dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
		
		//cout << visited_arcs << endl;
	
		usd[source] = true;
		return labeled_arcs;
	}

	
	int labeling_source_dij(NodeID source, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, queue<NodeID> &visited_que, vector<EdgeWeight>& distances, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<bool>& usd, NodeID ranking, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& r_tmp_idx, WGraph& wgraph) {
		NodeID labeled_arcs = 0;

		NodeID visited_arcs = 0;
	/*	vector<vector<NodeID> >& adj = wgraph.adj;
		vector<vector<EdgeWeight> >& adj_w = wgraph.adj_weight;*/

		const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];
		const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[source];

		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
		}

		pqueue.update(source, 0);


		while (!pqueue.empty()) {
			NodeID v;
			EdgeWeight v_d;
			pqueue.extract_min(v, v_d);
			vis[v] = true;
			visited_que.push(v);

			if (usd[v]) continue;


			if (DIRECTED_FLAG == false) {
				pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
				for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
					NodeID w = tmp_idx_v.first[i];
					EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
					if (td <= v_d) {
						goto pruned_forward;
					}
				}

				// arcs touched so far when building this labels
				labeled_arcs++;

				tmp_idx_v.first.back() = ranking;
				tmp_idx_v.second.back() = v_d;
				tmp_idx_v.first.push_back(numOfVertices);
				tmp_idx_v.second.push_back(INF_WEIGHT);
			}
			else {
				pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];

				// Pruned by the forward labels of r and backward labels of v in the forward search from r when reaching v.
				for (size_t i = 0; i < r_tmp_idx_v.first.size(); ++i) {
					NodeID w = r_tmp_idx_v.first[i];
					EdgeWeight td = r_tmp_idx_v.second[i] + dst_r[w];
					if (td <= v_d) {
						goto pruned_forward;
					}
				}

				// arcs touched so far when building this labels
				labeled_arcs++;
				
				r_tmp_idx_v.first.back() = ranking;
				r_tmp_idx_v.second.back() = v_d;
				r_tmp_idx_v.first.push_back(numOfVertices);
				r_tmp_idx_v.second.push_back(INF_WEIGHT);
			}

	/*		for (size_t i = 0; i < adj[v].size(); ++i) {
			
				NodeID w = adj[v][i];
				EdgeWeight w_d = adj_w[v][i] + v_d;*/

			//Array Representation
			
			visited_arcs += (wgraph.vertices[v + 1] - wgraph.vertices[v]);
			
			for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid){

				NodeID w = wgraph.edges[eid].first;
				EdgeWeight w_d = wgraph.edges[eid].second + v_d;

		//		if (wgraph.adj[v][eid - wgraph.vertices[v]] != w || wgraph.adj_weight[v][eid - wgraph.vertices[v]] + v_d != w_d) cout << "we got problems labeling tree." << endl;
			
				if (!vis[w]) {
					if (distances[w] > w_d) {
						pqueue.update(w, w_d);
						distances[w] = w_d;
					}
				}
			}
		pruned_forward:
			{}
		}
		
		// Clear foward search structures.
		while (!visited_que.empty()) {
			NodeID vis_v = visited_que.front();
			visited_que.pop();
			vis[vis_v] = false;
			distances[vis_v] = INF_WEIGHT;
		}

		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
		}


		// Backward search.
		if (DIRECTED_FLAG == true) {

	/*		vector<vector<NodeID> >& r_adj = wgraph.r_adj;
			vector<vector<EdgeWeight> >& r_adj_w = wgraph.r_adj_weight;*/

			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i) {
				dst_r[r_tmp_idx_r.first[i]] = r_tmp_idx_r.second[i];
			}


			pqueue.update(source, 0);
			while (!pqueue.empty()) {
				NodeID v;
				EdgeWeight v_d;
				pqueue.extract_min(v, v_d);
				vis[v] = true;
				visited_que.push(v);


				if (usd[v]) continue;


				pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];

				for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
					NodeID w = tmp_idx_v.first[i];
					EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
					if (td <= v_d) {
						goto pruned_backward;
					}
				}

				// arcs touched so far when building this labels
				labeled_arcs++;
				tmp_idx_v.first.back() = ranking;
				tmp_idx_v.second.back() = v_d;
				tmp_idx_v.first.push_back(numOfVertices);
				tmp_idx_v.second.push_back(INF_WEIGHT);

		/*		for (size_t i = 0; i < r_adj[v].size(); ++i) {
					NodeID w = r_adj[v][i];
					EdgeWeight w_d = r_adj_w[v][i] + v_d;*/
				
				////Array Representation
				
				visited_arcs += (wgraph.r_vertices[v + 1] - wgraph.r_vertices[v]);
				
				for (EdgeID eid = wgraph.r_vertices[v]; eid < wgraph.r_vertices[v + 1]; ++eid){

					NodeID w = wgraph.r_edges[eid].first;
					EdgeWeight w_d = wgraph.r_edges[eid].second + v_d ;

					if (!vis[w]) {
						if (distances[w] > w_d) {
							pqueue.update(w, w_d);
							distances[w] = w_d;
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
		}
	
		
		//cout << visited_arcs << endl;
	
		usd[source] = true;
		return labeled_arcs;
	}
		
	int labeling_source_bfs_path(NodeID source, vector<NodeID>& que, vector<bool>& vis, vector<NodeID>& parents, vector<EdgeWeight>& dst_r, vector<bool>& usd, NodeID ranking, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& r_tmp_idx, vector<pair<vector<NodeID>, vector<NodeID> > >& tmp_idx_parent, vector<pair<vector<NodeID>, vector<NodeID> > >& r_tmp_idx_parent, Graph& graph) {
		NodeID labeled_arcs = 0;

		//vector<vector<NodeID> >& adj = graph.adj;
		const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];
		const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[source];

		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
		}

		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;

		que[que_h++] = source;
		vis[source] = true;
		que_t1 = que_h;
		parents[source] = source;

		for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID v = que[que_i];
				// hubs touched so far when building this labels
				if (usd[v]) continue;


				if (DIRECTED_FLAG == false) {
					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
					pair<vector<NodeID>, vector<NodeID> > &tmp_idx_parent_v = tmp_idx_parent[v];

					for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
						NodeID w = tmp_idx_v.first[i];
						EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
						if (td <= d) {
							goto pruned_forward;
						}
					}

					labeled_arcs++;
					tmp_idx_v.first.back() = ranking;
					tmp_idx_v.second.back() = d;
					tmp_idx_v.first.push_back(numOfVertices);
					tmp_idx_v.second.push_back(INF_WEIGHT);


					tmp_idx_parent_v.first.back() = ranking;
					tmp_idx_parent_v.second.back() = parents[v];
					tmp_idx_parent_v.first.push_back(numOfVertices);
					tmp_idx_parent_v.second.push_back(numOfVertices);
				}


				// Array Representation
				/*for (size_t i = 0; i < adj[v].size(); ++i) {
				NodeID w = adj[v][i];*/
				for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid) {

					NodeID w = graph.edges[eid];

					if (!vis[w]) {
						que[que_h++] = w;
						vis[w] = true;
						parents[w] = v;
					}
				}
			pruned_forward:
				{}
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}

		for (size_t i = 0; i < que_h; ++i) {
			vis[que[i]] = false;
			parents[que[i]] = numOfVertices;
		}
		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
		}

		parents[source] = source;
		

		usd[source] = true;
		return labeled_arcs;
	}
 
	
	int labeling_source_dij_path(NodeID source, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, queue<NodeID> &visited_que, vector<EdgeWeight>& distances, vector<bool>& vis, vector<NodeID>& parents, vector<EdgeWeight>& dst_r, vector<bool>& usd, NodeID ranking, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& r_tmp_idx, vector<pair<vector<NodeID>, vector<NodeID> > >& tmp_idx_parent, vector<pair<vector<NodeID>, vector<NodeID> > >& r_tmp_idx_parent, WGraph& wgraph) {
		NodeID labeled_arcs = 0;

		/*	vector<vector<NodeID> >& adj = wgraph.adj;
		vector<vector<EdgeWeight> >& adj_w = wgraph.adj_weight;*/

		const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];
		const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[source];

		//const pair<vector<NodeID>, vector<NodeID> > &tmp_idx_parent_r = tmp_idx_parent[source];
		//const pair<vector<NodeID>, vector<NodeID> > &r_tmp_idx_parent_r = r_tmp_idx_parent[source];

		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
		}

		pqueue.update(source, 0);
		parents[source] = source;

		while (!pqueue.empty()) {
			NodeID v;
			EdgeWeight v_d;
			pqueue.extract_min(v, v_d);
			vis[v] = true;
			visited_que.push(v);

			if (usd[v]) continue;


			if (DIRECTED_FLAG == false) {
				pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];

				pair<vector<NodeID>, vector<NodeID> > &tmp_idx_parent_v = tmp_idx_parent[v];

				for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
					NodeID w = tmp_idx_v.first[i];
					EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
					if (td <= v_d) {
						goto pruned_forward;
					}
				}

				// arcs touched so far when building this labels
				labeled_arcs++;

				tmp_idx_v.first.back() = ranking;
				tmp_idx_v.second.back() = v_d;
				tmp_idx_v.first.push_back(numOfVertices);
				tmp_idx_v.second.push_back(INF_WEIGHT);

				tmp_idx_parent_v.first.back() = ranking;
				tmp_idx_parent_v.second.back() = parents[v];
				tmp_idx_parent_v.first.push_back(numOfVertices);
				tmp_idx_parent_v.second.push_back(numOfVertices);

				
			}
			else {
				pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];

				pair<vector<NodeID>, vector<NodeID> > &r_tmp_idx_parent_v = r_tmp_idx_parent[v];

				// Pruned by the forward labels of r and backward labels of v in the forward search from r when reaching v.
				for (size_t i = 0; i < r_tmp_idx_v.first.size(); ++i) {
					NodeID w = r_tmp_idx_v.first[i];
					EdgeWeight td = r_tmp_idx_v.second[i] + dst_r[w];
					if (td <= v_d) {
						goto pruned_forward;
					}
				}

				// arcs touched so far when building this labels
				labeled_arcs++;

				r_tmp_idx_v.first.back() = ranking;
				r_tmp_idx_v.second.back() = v_d;
				r_tmp_idx_v.first.push_back(numOfVertices);
				r_tmp_idx_v.second.push_back(INF_WEIGHT);

				r_tmp_idx_parent_v.first.back() = ranking;
				r_tmp_idx_parent_v.second.back() = parents[v];
				r_tmp_idx_parent_v.first.push_back(numOfVertices);
				r_tmp_idx_parent_v.second.push_back(numOfVertices);

			}

			/*		for (size_t i = 0; i < adj[v].size(); ++i) {

			NodeID w = adj[v][i];
			EdgeWeight w_d = adj_w[v][i] + v_d;*/

			//Array Representation
			for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid) {

				NodeID w = wgraph.edges[eid].first;
				EdgeWeight w_d = wgraph.edges[eid].second + v_d;

				//		if (wgraph.adj[v][eid - wgraph.vertices[v]] != w || wgraph.adj_weight[v][eid - wgraph.vertices[v]] + v_d != w_d) cout << "we got problems labeling tree." << endl;

				if (!vis[w]) {
					if (distances[w] > w_d) {
						pqueue.update(w, w_d);
						distances[w] = w_d;
						parents[w] = v;
					}
				}
			}
		pruned_forward:
			{}
		}

		// Clear foward search structures.
		while (!visited_que.empty()) {
			NodeID vis_v = visited_que.front();
			visited_que.pop();
			vis[vis_v] = false;
			distances[vis_v] = INF_WEIGHT;
			parents[vis_v] = numOfVertices;
		}

		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
		}

		parents[source] = source;

		// Backward search.
		if (DIRECTED_FLAG == true) {

			/*		vector<vector<NodeID> >& r_adj = wgraph.r_adj;
			vector<vector<EdgeWeight> >& r_adj_w = wgraph.r_adj_weight;*/

			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i) {
				dst_r[r_tmp_idx_r.first[i]] = r_tmp_idx_r.second[i];
			}


			pqueue.update(source, 0);
			while (!pqueue.empty()) {
				NodeID v;
				EdgeWeight v_d;
				pqueue.extract_min(v, v_d);
				vis[v] = true;
				visited_que.push(v);


				if (usd[v]) continue;


				pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
				pair<vector<NodeID>, vector<NodeID> > &tmp_idx_parent_v = tmp_idx_parent[v];

				for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
					NodeID w = tmp_idx_v.first[i];
					EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
					if (td <= v_d) {
						goto pruned_backward;
					}
				}

				// arcs touched so far when building this labels
				labeled_arcs++;
				tmp_idx_v.first.back() = ranking;
				tmp_idx_v.second.back() = v_d;
				tmp_idx_v.first.push_back(numOfVertices);
				tmp_idx_v.second.push_back(INF_WEIGHT);

				tmp_idx_parent_v.first.back() = ranking;
				tmp_idx_parent_v.second.back() = parents[v];
				tmp_idx_parent_v.first.push_back(numOfVertices);
				tmp_idx_parent_v.second.push_back(numOfVertices);

				/*		for (size_t i = 0; i < r_adj[v].size(); ++i) {
				NodeID w = r_adj[v][i];
				EdgeWeight w_d = r_adj_w[v][i] + v_d;*/

				////Array Representation
				for (EdgeID eid = wgraph.r_vertices[v]; eid < wgraph.r_vertices[v + 1]; ++eid) {

					NodeID w = wgraph.r_edges[eid].first;
					EdgeWeight w_d = wgraph.r_edges[eid].second + v_d;

					if (!vis[w]) {
						if (distances[w] > w_d) {
							pqueue.update(w, w_d);
							distances[w] = w_d;
							parents[w] = v;
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
				parents[vis_v] = numOfVertices;
			}
			parents[source] = numOfVertices;
			pqueue.clear_n();

			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i)
				dst_r[r_tmp_idx_r.first[i]] = INF_WEIGHT;
		}

		usd[source] = true;
		return labeled_arcs;
	}

	void shrink_tree(int tid, NodeID root, large_tree& ltree_tid) {
		small_tree new_stree;
		//GOOGLE DENSE MAP
		new_stree.set_empty_key(numOfVertices+1);
		new_stree.set_deleted_key(numOfVertices+2);
		for (NodeID v = 0; v < ltree_tid.size(); ++v) {
			int pt = ltree_tid[v];
			if (pt != -1)
				new_stree.insert(make_pair(v, pt));
		}
		new_stree.insert(make_pair(root, -1));
		strees.push_back(new_stree);
		trees_pointers[tid] = strees.size() - 1;
		alive_resources += (long long)new_stree.size();

		//ltrees[tid].clear();
		// FIX
		ltree_tid.clear();
		ltree_tid.shrink_to_fit();
		alive_resources -= (long long)numOfVertices;
	}

	void enlarge_tree(int tid, NodeID root, small_tree& stree_tid, large_tree& ltree_tid) {
		ltree_tid.resize(numOfVertices,-1);
		
		for (small_tree::iterator it = stree_tid.begin(); it != stree_tid.end(); ++it) {			
			NodeID v = (*it).first;
			NodeID pt = (*it).second;
			ltree_tid[v] = pt;
		}

		ltree_tid[root] = -1;

		ltrees.push_back(ltree_tid);

		trees_pointers[tid] = -1;
		alive_resources += (long long)numOfVertices;

		alive_resources -= (long long)stree_tid.size();
		stree_tid.clear();
		stree_tid.resize(0);
	}

	int get_large_parent(NodeID v, large_tree& ltree) {
		return ltree[v];
	}

	int get_small_parent_withouttest(NodeID v, small_tree& stree) {

		return stree[v];

	}

	int get_small_parent(NodeID v, small_tree& stree) {
		small_tree::iterator it = stree.find(v);

		if (it != stree.end())
			return stree[v];
		else
			return -1;
	}


	void remove_large_node(NodeID v, large_tree& ltree) {
		ltree[v] = -1;
	}
	void remove_small_node(NodeID v, small_tree& stree) {
		stree.erase(v);
	//	stree.resize(0);
		//stree[v] = -1;
	}

	//GOOGLE
	void cnt16_2max_queue(vector<cover_value_type>& picked_value, vector<cover_value_type>& sum, vector<bool>& has_picked, /*unordered_set<NodeID>& affected_v*/ google::dense_hash_set<NodeID>& affected_v, vector<bool>& affected_list, const vector<vector<cover_value_type> >& cnt16, bool& switch_small, bool& switch_set, benchmark::heap<2, cover_value_type, NodeID>& value_heap) {

		// First Round, build the entire cnt16 one time.
		if (picked_value.size() == 0) {
			picked_value.resize(numOfVertices, 0);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				NodeID max1 = 0;
				NodeID max2 = 0;
				const vector<cover_value_type>& cnt16_v = cnt16[v];
				cover_value_type &picked_value_v = picked_value[v];
				cover_value_type& sum_v = sum[v];
				sum_v = 0;
				for (int i = 0; i < CNT; ++i) {
					sum_v += cnt16_v[i];
					if (cnt16_v[i] > max2)
						max2 = cnt16_v[i];
					if (max2 > max1) {
						NodeID tmp = max1;
						max1 = max2;
						max2 = tmp;
					}
				}
				picked_value_v = sum_v - max1 - max2;
				if (switch_small == true) {
					value_heap.update(v, -picked_value_v);
				}
			}
		}

		//for (unordered_set<NodeID>::iterator it = affected_v.begin(); it != affected_v.end(); ++it) {
		if (!switch_set) {
			for (int v = 0; v < numOfVertices; ++v) {
				//NodeID v = *it;
				if (affected_list[v] == false) continue;
				affected_list[v] = false;
				if (has_picked[v] == true) continue;
				//if (v == best_vid) continue;
				cover_value_type& sum_v = sum[v];

				/*	if (sum_v < best_so_far) {
				picked_value[v] = sum_v;
				if (switch_small == true)
				value_heap.update(v, sum_v);
				lazy_updates[v] = true;
				continue;
				}*/
				NodeID max1 = 0;
				NodeID max2 = 0;
				const vector<cover_value_type>& cnt16_v = cnt16[v];
				cover_value_type& picked_value_v = picked_value[v];


				sum_v = 0;
				for (int i = 0; i < CNT; ++i) {
						sum_v += cnt16_v[i];
					if (cnt16_v[i] > max2)
						max2 = cnt16_v[i];
					if (max2 > max1) {
						int tmp = max1;
						max1 = max2;
						max2 = tmp;
					}
				}
				picked_value_v = sum_v - max1 - max2;
				//if (picked_value_v > best_so_far) best_so_far = picked_value_v;
				if (switch_small == true) {
					value_heap.update(v, -picked_value_v);
				}
			}
		}
		else {
			//GOOGLE
			
			//for (unordered_set<NodeID>::iterator it = affected_v.begin(); it != affected_v.end(); ++it) {
			for (google::dense_hash_set<NodeID>::iterator it = affected_v.begin(); it != affected_v.end(); ++it) {

				NodeID v = *it;
				affected_list[v] = false;
				if (has_picked[v] == true) continue;
				//if (v == best_vid) continue;
				cover_value_type& sum_v = sum[v];

				//if (sum_v < best_so_far) {
				//	picked_value[v] = sum_v;
				//	if (switch_small == true)
				//		value_heap.update(v, sum_v);
				//	lazy_updates[v] = true;
				//	continue;
				//}
				NodeID max1 = 0;
				NodeID max2 = 0;
				const vector<cover_value_type>& cnt16_v = cnt16[v];
				cover_value_type& picked_value_v = picked_value[v];


				sum_v = 0;
				for (int i = 0; i < CNT; ++i) {
						sum_v += cnt16_v[i];
					if (cnt16_v[i] > max2)
						max2 = cnt16_v[i];
					if (max2 > max1) {
						int tmp = max1;
						max1 = max2;
						max2 = tmp;
					}
				}
				picked_value_v = sum_v - max1 - max2;
				//if (picked_value_v > best_so_far) best_so_far = picked_value_v;
				if (switch_small == true) {
					value_heap.update(v, -picked_value_v);
				}
			}
			affected_v.clear();
		}
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
		
		clabels.total_children = children_size;
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

	
	Betweenness_Ordering(NodeID k, double beta, Graph& graph, NodeID stopv, long long bound_times, long long trans_times, bool slevel, bool directed_flags) {
		
		double start_time = GetCurrentTimeSec();
		DIRECTED_FLAG = directed_flags;
		SECOND_LEVEL = slevel;
		
		if (k < CNT) {
			cout << " Should pick at least 16 seed!" << endl;
			return;
		}
		iteration_time.resize(numOfVertices);

		updating_time = 0;
		labeling_time = 0;
		adding_time = 0;
		selecting_time = 0;
		init_time = 0;
		u1_time = 0;
		u2_time = 0;
		u3_time = 0;
		num_of_trees = 0;
		total_resources = 0;
		get_small_parent_time = 0;

		inv.resize(numOfVertices, -1);
		rank.resize(numOfVertices, -1);

		if (DIRECTED_FLAG == false) {
			labels.index_.resize(numOfVertices);
		}
		else { //if (DIRECTED_FLAG == true) {
			dlabels.index_.resize(numOfVertices);
			dlabels.bindex_.resize(numOfVertices);
		}

		// Resource Bounding.
		NodeID new_tree_arcs = 0;
		NodeID labeled_arcs = 0;
		alive_resources = 0;
		bounded_resources = (long long)10 * (long long)k * (long long)numOfVertices;

		bool stop_flag = false;

		vector<NodeID> st_roots; // The root vertices of the sampled trees.
		vector<NodeID> st_sizes; // The size of the sampled trees.

		NodeID size_thr = numOfVertices < 8 ? numOfVertices : numOfVertices / 8; // When a tree shrink to a size smaller than the size_thr, it should be converted into a small_tree representation.
																				 //	vector<vector<int> > cover;


		init_time = GetCurrentTimeSec();
		// Building the initial k trees.
	//	srand((unsigned)time(NULL));
		srand(0x728f9d0);
		

		vector<NodeID> descendants_parents(numOfVertices, -1);

		ltrees.reserve(numOfVertices);
		strees.reserve(numOfVertices);
		trees_pointers.reserve(numOfVertices);
		st_roots.reserve(numOfVertices);
		st_sizes.reserve(numOfVertices);
		st_roots.resize(k);
		st_sizes.resize(k);

		ltrees.resize(k);
		strees.clear();
		trees_pointers.resize(k, -1); // Point the sampling trees to its small tree representations.
									  //	cover.resize(k, vector<int>(V, 0));

		vector<NodeID> ascendants;
		ascendants.reserve(numOfVertices);
		vector<NodeID> descendants;
		descendants.reserve(numOfVertices);
		vector<NodeID> que(numOfVertices);
		vector<bool> vis(numOfVertices);
		vector<bool> has_picked(numOfVertices);
		vector<vector<cover_value_type> > cnt16(numOfVertices, vector<cover_value_type>(CNT));
		//		vector<bool> affected(V, true);
		vector<bool> picked_as_root(numOfVertices, false);
		vector<bool> trivial_root(numOfVertices, false);
		int picked_count = 0;
		int used_as_cover = 0;
		vector<cover_value_type> sum(numOfVertices);
		vector<bool> usd(numOfVertices, false);
		vector<cover_value_type> cover_tid(numOfVertices, 0);

		//unordered_set<NodeID> affected_v; //GOOGLE
		google::dense_hash_set<NodeID> affected_v;
		affected_v.set_empty_key(numOfVertices + 1);
		affected_v.set_deleted_key(numOfVertices + 2);

		vector<bool> affected_list(numOfVertices, false);
		int affected_count = 0;
		bool switch_set = false;

		//vector<vector<int> > cover_record(V, vector<int>(V, 0));

		//vector<unordered_set<int> > reverse_index(0); //GOOGLE
//		vector<unordered_set<int> > reverse_index(0);
		vector<google::dense_hash_set<NodeID> > reverse_index;
		reverse_index.reserve(numOfVertices);

		// Preparing basic structure for pl algorithm.
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<EdgeWeight> r_dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT))); // Backward labels.

		
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx_token_parents(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, numOfVertices)));
				
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx_token_parents(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, numOfVertices)));

													 //	int best_sum = 0;
													 //	NodeID best_sum_vid = 0;
													 //vector<bool> lazy_updates(V, false);

													 // Building the initial k sampled trees. 
		for (int i = 0; i < k; ++i) {
			NodeID r = rand() % numOfVertices;
			if (picked_as_root[r] == true) { --i; continue; }

			st_roots[i] = r;
			large_tree &ltree_i = ltrees[i];
			picked_as_root[r] = true;
			picked_count++;
			used_as_cover++;

			// Array Representation
			//build_large_tree_bfs(r, ltree_i, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph.adj, 0, numOfVertices);
			build_large_tree_bfs(r, ltree_i, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph, 0, numOfVertices);

			st_sizes[i] = descendants.size();
			alive_resources += (long long)numOfVertices;

			// Turn large tree to small tree.
			// 1. Build new small tree with ltree_i and add it to strees;
			// 2. Set tree_pointers[i] to the position of the strees;
			// 3. Clear the large tree ltree_i by setting its size to 0.
			if (st_sizes[i] < size_thr)
				shrink_tree(i, st_roots[i], ltrees[i]);


			// Calculate the subtree size of each nodes			.
			//	vector<int>& cover_i = cover[i];

			bool isLarge = false;
			int stree_idx = trees_pointers[i];
			if (stree_idx == -1)
				isLarge = true;

			small_tree& stree_i = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;

			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];

				++cover_tid[q]; // Increment because it hitsthe r-q path.
				++cnt16[q][i%CNT];
				//++sum[q];
				//	if (sum[q] > best_sum) {
				//		best_sum = sum[q];
				//		best_sum_vid = q;
				//	}
				////	cover_record[i][q]++;

				int parent_node;
				if (isLarge)
					parent_node = get_large_parent(q, ltree_i);
				else
					parent_node = get_small_parent(q, stree_i);


				if (parent_node > -1) {
					cover_tid[parent_node] += cover_tid[q];
					cnt16[parent_node][i%CNT] += cover_tid[q];
					//sum[parent_node] += cover_i[q];
					//	cover_record[i][parent_node] += cover_record[i][q];
				}
			}	
			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];
				cover_tid[q] = 0;
			}
		}


		//cnt16_aggreate(cnt16, has_picked, cover);

		vector<cover_value_type> picked_value;
		benchmark::heap<2, cover_value_type, NodeID> value_heap(numOfVertices);
		bool switch_small = false;

		//cnt16_2max(picked_value, sum, has_picked,affected, cnt16, switch_small, value_heap);
		//		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap, best_sum_vid, lazy_updates);
		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);

		init_time = GetCurrentTimeSec() - init_time;


		cout << "Finish initial sampling " << init_time << "secs." << endl;

		//for (NodeID v = 0; v < V; ++v)	value_heap.update(v, -picked_value[v]); // Min-heap.

		double max_value = -1;


		// Perform V iterations to rank all vertices.
		NodeID ranking = 0;
		for (; ranking < numOfVertices; ++ranking) {
		/*	if (ranking % (numOfVertices / 10) == 0)
				cout << "picking the " << ranking << "-th orders." << endl;*/

			double time_this_iteration = GetCurrentTimeSec();

			//best_sum = 0;

			NodeID pv = 0;
			cover_value_type pv_cover = -1;

			if (ranking == stopv) {
				stop_flag = true;
				cout << "topv labeling time:" << labeling_time << endl;
				break;
			}
			double selecting_acc_time = GetCurrentTimeSec();
			// In the beginning, pick the best cover vertex by a linear scan.
			if (switch_small == false) {
				double u1_acc_time = GetCurrentTimeSec();
				//while (true) {
				for (NodeID v = 0; v < numOfVertices; ++v) {
					if (has_picked[v] == true) continue;
					if (pv_cover == picked_value[v]) {
						if (sum[pv] < sum[v]) {
							pv_cover = picked_value[v];
							pv = v;
						}
					}
					else if (pv_cover < picked_value[v]) {
						pv_cover = picked_value[v];
						pv = v;
					}
				}
				/*	if (lazy_updates[pv] == true) {
				NodeID max1 = 0;
				NodeID max2 = 0;
				const vector<int>& cnt16_pv = cnt16[pv];
				int& picked_value_v = picked_value[pv];

				int sum_v = sum[pv];
				for (int i = 0; i < CNT; ++i) {
				if (cnt16_pv[i] > max2)
				max2 = cnt16_pv[i];
				if (max2 > max1) {
				int tmp = max1;
				max1 = max2;
				max2 = tmp;
				}
				}
				picked_value_v = sum_v - max1 - max2;

				lazy_updates[pv] = false;
				}
				else
				break;*/
				//}

				if (sum[pv] < numOfVertices / 8 && switch_small == false) {
					cout << ranking << endl;
					cout << pv_cover << endl;
					cout << sum[pv] << endl;
					cout << picked_count << endl;
					cout << "Switch to max-heap representations." << endl;
					switch_small = true;
				}
				u1_time += (GetCurrentTimeSec() - u1_acc_time);
			}
			else {// Later, utilizing the max-heap directly
				double u2_acc_time = GetCurrentTimeSec();
				//while (true) {
				value_heap.extract_min(pv, pv_cover);
				pv_cover = -pv_cover;
				/*	if (lazy_updates[pv] == true) {
				NodeID max1 = 0;
				NodeID max2 = 0;
				const vector<int>& cnt16_pv = cnt16[pv];
				int& picked_value_v = picked_value[pv];

				int sum_v = sum[pv];
				for (int i = 0; i < CNT; ++i) {
				if (cnt16_pv[i] > max2)
				max2 = cnt16_pv[i];
				if (max2 > max1) {
				int tmp = max1;
				max1 = max2;
				max2 = tmp;
				}
				}
				picked_value_v = sum_v - max1 - max2;
				value_heap.update(pv, -picked_value_v);
				lazy_updates[pv] = false;
				}
				else
				break;
				}*/
				u2_time += (GetCurrentTimeSec() - u2_acc_time);

			}

			selecting_time += (GetCurrentTimeSec() - selecting_acc_time);

			//// FIX
			//if (switch_small == true && reverse_index.size() > 0 && sum[pv] <= reverse_index[pv].size()) {
			//	cout << "I think we can stop here " << GetCurrentTimeSec() - start_time << " ranking:" << ranking << endl;
			//	break;
			//}

			has_picked[pv] = true;
			double labeling_acc_time = GetCurrentTimeSec();
			labeled_arcs = labeling_source_bfs(pv, que, vis, dst_r, r_dst_r, usd, ranking, tmp_idx, r_tmp_idx, tmp_idx_token_parents, r_tmp_idx_token_parents, graph);
			labeling_time += (GetCurrentTimeSec() - labeling_acc_time);
			inv[ranking] = pv;
			rank[pv] = ranking;

			//cout << ranking << "\t" << labeled_arcs << endl;


			int current_trees_num = st_roots.size();

			// 4 steps in updating:
			double updating_acc_time = GetCurrentTimeSec();
			if (switch_small == false) { // We still iterate all alive trees to calculate the cover.	
				for (int tid = 0; tid < current_trees_num; ++tid) {

					// This tree has already be taken because its root has been selected.
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					// pv has already been covered in tree tid.
					int pvpid = 0;
					if (isLarge)
						pvpid = get_large_parent(pv, ltree_tid);
					else
						pvpid = get_small_parent(pv, stree_tid);
					if (pvpid == -1 && st_roots[tid] != pv) // pv has been covered in tree tid 
						continue;


					// Get the descendants and ascendants of pv in tree tid.

					// Array Representation
					//// 1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;				 
					//get_descendants(tid, pv, descendants, descendants_parents, que, vis, graph.adj, isLarge, stree_tid, ltree_tid);

					//get_ascendants(tid, pv, ascendants, que, vis, graph.adj, graph.r_adj, isLarge, stree_tid, ltree_tid);
				
					get_descendants(tid, pv, descendants, descendants_parents, que, vis, graph, isLarge, stree_tid, ltree_tid);

					get_ascendants(tid, pv, ascendants, que, vis, graph, isLarge, stree_tid, ltree_tid);

					NodeID subtract_size = descendants.size();
					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
												//cover[tid][av] -= subtract_size;
						cnt16[av][tid%CNT] -= subtract_size;
						//sum[av] -= subtract_size;
						//	cover_record[tid][av] -= subtract_size;
						//	if (cover_record[tid][av] < 0)
						//	cout << "1: " << cover_record[tid][av] << endl;
						if (affected_list[av] == false) {
							affected_list[av] = true;
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
						}
					}


					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
				//	vector<int> cover_tid(numOfVertices, 0);
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];

						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//sum[dv] -= cover_tid[dv];
						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];

					/*	if (cnt16[dv][tid%CNT] < 0)
							cout << "a6 " << cnt16[dv][tid%CNT] << endl;*/

						//	cover_record[tid][dv] -= cover_tid[dv];
						//	if (cover_record[tid][dv] < 0) {
						//	cout << "2: " << cover_record[tid][dv]  << " " << cover_tid[dv]  << " " << tid << " " << dv <<  endl;
						//	}

						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else 
							remove_small_node(dv, stree_tid);

						
						//	remove_node(tid, dv);
						//affected[dv] = true;
						if (affected_list[dv] == false) {
							affected_count++;
							affected_list[dv] = true;
							if (switch_set)
								affected_v.insert(dv);
						}
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;
					}


					// Drop the size of subtree.
					st_sizes[tid] -= subtract_size;
					if (trees_pointers[tid] != -1) alive_resources -= (cover_value_type)subtract_size;
				//	if (st_sizes[tid] < 0) {
				//		cout << "a5" << st_sizes[tid] << endl;
				//	}
					if (st_sizes[tid] > 0) {
						// 3. Convert shrink trees to small trees representataions.				
						if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
							shrink_tree(tid, st_roots[tid], ltrees[tid]);
						}
					}

					if (st_roots[tid] == pv)
						if (isLarge) {
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
							//ltrees.shrink_to_fit();
						}
						else {
							stree_tid.clear();
							stree_tid.resize(0);
							//strees.shrink_to_fit();
						}

				}
			}
			else { // We use the reverse index to update trees.

				   // Build the reverse_index first (one time).
				if (reverse_index.size() == 0) {
					reverse_index.resize(numOfVertices);
					for (NodeID v = 0; v < numOfVertices; ++v) {
						reverse_index[v].set_empty_key(numOfVertices + 1);
						reverse_index[v].set_deleted_key(numOfVertices + 2);
					}
					for (int tid = 0; tid < current_trees_num; ++tid) {

						if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue; //This whole tree has already been covered. 

						bool isLarge = false;
						int stree_idx = trees_pointers[tid];
						if (stree_idx == -1)
							isLarge = true;

						small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
						large_tree& ltree_tid = ltrees[tid];

						for (int v = 0; v < numOfVertices; ++v) {

							// Get the descendants and ascendants of pv in tree tid.

							//	if (cover[tid].find(pv) == cover[tid].end()) continue; // This pv has been covered in this tree before.
							//	if (cover[tid][pv] <= 0) continue; // modified-5-14
							int vpid = 0;
							if (isLarge)
								vpid = get_large_parent(v, ltree_tid);
							else
								vpid = get_small_parent(v, stree_tid);

							// To testify whether v is still uncovered in tree tid
							// v has parent vpid != -1, but if v is the root of tree tid, it still can be uncovered.
							if ((vpid != -1 || st_roots[tid] == v))
								reverse_index[v].insert(tid);
							//if (cover[tid].find(v) != cover[tid].end())
							//	reverse_index[v].insert(tid);
						}
					}
				}

				double u4_acc_time = GetCurrentTimeSec();

				//	double u4_acc_time = GetCurrentTimeSec();
				// Access the alive trees directly by the reverse index.
				//				cout << reverse_index[pv].size() << " " << sum[pv] << endl;

				
//				for (unordered_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) { // GOOGLE
				for (google::dense_hash_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) {


					int tid = *it;
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					// Get the descendants and ascendants of pv in tree tid.

					//	if (cover[tid].find(pv) == cover[tid].end()) continue; // This pv has been covered in this tree before.
					//	if (cover[tid][pv] <= 0) continue; // modified-5-14

					int pvpid = 0;
					if (isLarge)
						pvpid = get_large_parent(pv, ltree_tid);
					else {
						//	double get_small_parent_incre_time = GetCurrentTimeSec();
						pvpid = get_small_parent(pv, stree_tid);
						//	get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}

					if (pvpid == -1 && st_roots[tid] != pv)
						continue;

					// Array Representation
					//// 1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;
					//get_descendants(tid, pv, descendants, descendants_parents, que, vis, graph.adj, isLarge, stree_tid, ltree_tid);

					//get_ascendants(tid, pv, ascendants, que, vis, graph.adj, graph.r_adj, isLarge, stree_tid, ltree_tid);

					get_descendants(tid, pv, descendants, descendants_parents, que, vis, graph, isLarge, stree_tid, ltree_tid);

					get_ascendants(tid, pv, ascendants, que, vis, graph, isLarge, stree_tid, ltree_tid);

					NodeID subtract_size = descendants.size();

					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
												//						cover[tid][av] -= subtract_size;

						cnt16[av][tid%CNT] -= subtract_size;
						//sum[av] -= subtract_size;

						//	cover_record[tid][av] -= subtract_size;
						//	if (cover_record[tid][av] < 0) {
						//	cout << "4: " << cover_record[tid][av] << " " << subtract_size << " " << tid << " " << endl;
						//	}

						//affected[av] = true;
						/*if (cnt16[av][tid%CNT] < 0)
						cout << " a2 " << cnt16[av][tid%CNT]  << " " << av << " " << st_roots[tid] << " " << tid << endl;*/
						if (affected_list[av] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
							affected_list[av] = true;
						}
					}

					//double u2_acc_time = GetCurrentTimeSec();
					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					//vector<int> cover_tid(numOfVertices, 0);
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
					//	sum[dv] -= cover_tid[dv];
						//cover[tid].erase(dv);
						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];

						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else 
							remove_small_node(dv, stree_tid);
							
						//	remove_node(tid, dv);
						//affected[dv] = true;

						//	cover_record[tid][dv] -= cover_tid[dv];
						//if (cover_record[tid][dv] < 0)
						//	cout << "3: " << cover_record[tid][dv] << " " << cover_tid[dv] << " " << di << " " << descendants.size() << " " << picked_as_root[st_roots[tid]] << endl;

						/*if (cnt16[dv][tid%CNT]  < 0)
						cout << " a3 " << cnt16[dv][tid%CNT] << " " << dv << " " << st_roots[tid] << " " << tid << endl;
						*/
						if (dv != pv)
							reverse_index[dv].erase(tid);
						//in_tree_count[dv]--;


						if (affected_list[dv] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(dv);
							affected_list[dv] = true;
						}
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;
					}

					if (st_roots[tid] == pv)
						if (isLarge) {
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
						}
						else {
							stree_tid.clear();
							stree_tid.resize(0);
						}

					st_sizes[tid] -= subtract_size;
				//	if (st_sizes[tid]  < 0)
				//		cout << " a4 " << st_sizes[tid] << endl;

					if (trees_pointers[tid] != -1) alive_resources -= subtract_size;


					// 3. Convert shrink trees to small trees representataions.	
					if (st_sizes[tid] > 0) {
						if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
							shrink_tree(tid, st_roots[tid], ltrees[tid]);
						}
					}

				}
				reverse_index[pv].clear();
				reverse_index[pv].resize(0);
				u4_time += (GetCurrentTimeSec() - u4_acc_time);
			}

			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			// 4. Sample new trees based on.
			// If we sample a new tree, we need to dynamically increment the following data structures:
			// 1) ltrees
			// 2) trees_pointers (-1 or pointer to strees)
			// 3) cover
			// 4) st_roots
			// 5) st_sizes 
			// 6) strees (maybe)
			// two constraints: cl > ct and 10kn total vertices
			double adding_acc_time = GetCurrentTimeSec();
			new_tree_arcs = 0;

			//cout << "ranking:" << ranking << endl;
			//cout << "pick cover:" << pv_cover << endl;
			//cout << "pick count:" << picked_count << endl;
			//cout << "labeled_arcs:" << labeled_arcs << endl;
			//cout << "new_tree_arcs:" << new_tree_arcs << endl;
			//cout << "alive_resources:" << alive_resources << endl;
			//cout << "bounded_resources:" << bounded_resources << endl;
			//cout << "#tree:" << cover.size() << endl;
			int how_many = 0;
			//cout << "beta*labeld_arcs: " << beta * labeled_arcs << " new tree _arcs: " << new_tree_arcs << " alive " << alive_resources << " bounded " << bounded_resources << endl;
			NodeID sam_this_round = 0;
			while (beta * labeled_arcs > new_tree_arcs && alive_resources < bounded_resources) {

				if (picked_count == numOfVertices) break;



				NodeID r = rand() % numOfVertices;
				if (picked_as_root[r] == true) continue;
				if (trivial_root[r] == true) continue;
				int tmp_tree_arcs = 0;

				large_tree new_ltree;
				small_tree new_stree;
				//GOOGLE DENSE MAP
				new_stree.set_empty_key(numOfVertices + 1);
				new_stree.set_deleted_key(numOfVertices + 2);

				//	BUILD_SMALL_TREE_FLAG = false;
				bool isLarge = false;
				if (!BUILD_SMALL_TREE_FLAG) {

					double u5_acc_time = GetCurrentTimeSec();
					// Array Representation
					//tmp_tree_arcs = build_large_tree_bfs(r, new_ltree, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph.adj, new_tree_arcs, beta * labeled_arcs);
					tmp_tree_arcs = build_large_tree_bfs(r, new_ltree, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph, new_tree_arcs, beta * labeled_arcs);
					u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;

					if (descendants.size() <size_thr)
						BUILD_SMALL_TREE_FLAG = true;

					if (descendants.size() < 2) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					isLarge = true;
					how_many++;
					used_as_cover++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					st_roots.push_back(r);
					trees_pointers.push_back(-1);
					st_sizes.push_back(descendants.size());

					//new_tree_arcs += descendants.size();
					alive_resources += numOfVertices;

					// Turn large tree to small tree.
					// 1. Build new small tree with ltree_i and add it to strees;
					// 2. Set tree_pointers[i] to the position of the strees;
					// 3. Clear the large tree ltree_i by setting its size to 0.
					int i = st_roots.size() - 1;
					if (st_sizes[i] < size_thr) {
						shrink_tree(i, st_roots[i], new_ltree);
						isLarge = false;
					}
					ltrees.push_back(new_ltree);
				}
				else {//Build Small Tree Directly
					double u5_acc_time = GetCurrentTimeSec();
					// Array Representation
					//	tmp_tree_arcs = build_small_tree_bfs(r, new_stree, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph.adj, new_tree_arcs, beta * labeled_arcs);
					tmp_tree_arcs = build_small_tree_bfs(r, new_stree, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph, new_tree_arcs, beta * labeled_arcs);
					u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;

					if (descendants.size() >= size_thr)
						BUILD_SMALL_TREE_FLAG = false;

					if (descendants.size() < 2) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					isLarge = false;

					how_many++;
					used_as_cover++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					st_roots.push_back(r);
					trees_pointers.push_back(strees.size()); //we assume it is inserted.
					st_sizes.push_back(descendants.size());

					alive_resources += new_stree.size();
					//strees.push_back(new_stree);
					//// new_ltree is empty.
					//ltrees.push_back(new_ltree);

					int i = st_roots.size() - 1;
					if (st_sizes[i] >= size_thr) {
						enlarge_tree(i, st_roots[i], new_stree, new_ltree);
						new_stree.clear();
						new_stree.resize(0);
						isLarge = true;
					}
					else {
						strees.push_back(new_stree);
						ltrees.push_back(new_ltree);
					}
				}

				int i = st_roots.size() - 1;

				int numOfStrees = strees.size();
				small_tree &new_stree_ref = numOfStrees != 0 ? strees[strees.size() - 1] : empty_small_tree;

				// Calculate the subtree size of each nodes.
				//				unordered_map<int, int> cover_new;
			//	vector<cover_value_type> cover_new(numOfVertices, 0); // modified-5-14

				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];

					++cover_tid[q]; // Increment because it hitsthe r-q path. 
					++cnt16[q][i%CNT];
					//++sum[q];
					/*	if (sum[q] > best_sum) {
					best_sum = sum[q];
					best_sum_vid = q;
					}*/
					//affected[q] = true;
					//	cover_record[i][q]++;
					if (affected_list[q] == false) {
						affected_count++;
						if (switch_set)
							affected_v.insert(q);
						affected_list[q] = true;
					}

					if (switch_small == true) {
						reverse_index[q].insert(i);
					}

					int parent_node = -1;

					if (isLarge)
						parent_node = get_large_parent(q, new_ltree);
					else
						parent_node = get_small_parent(q, new_stree_ref);

					if (parent_node > -1) {
						cover_tid[parent_node] += cover_tid[q];
						cnt16[parent_node][i%CNT] += cover_tid[q];
					//	sum[parent_node] += cover_new[q];
						//	cover_record[i][parent_node] += cover_record[i][q];
					}
				}
				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];
					cover_tid[q] = 0;
				}

				sam_this_round++;
			}

			//cout << ranking << "\t" << sam_this_round << endl;


			adding_time += (GetCurrentTimeSec() - adding_acc_time);
			//cout << ranking << ":  beta:" << beta << " beta * labeled_arcs " << beta * labeled_arcs << " new_tree_arcs " << new_tree_arcs << " how many " << how_many << endl;

			//double u3_acc_time = GetCurrentTimeSec();

			updating_acc_time = GetCurrentTimeSec();

			// Build the cover value heap for one time.
			if (switch_small == true && value_heap.empty())
				for (int nid = 0; nid < numOfVertices; ++nid)
					if (has_picked[nid] == false)
						value_heap.update(nid, -picked_value[nid]);

			//	cnt16_2max(picked_value, sum, has_picked, affected, cnt16, switch_small, value_heap);		
			//double u3_acc_time = GetCurrentTimeSec();

			//			cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap, best_sum_vid, lazy_updates);
			cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);
			//u3_time += (GetCurrentTimeSec() - u3_acc_time);

			if (affected_count < numOfVertices / 8)
				switch_set = true;
			else
				switch_set = false;
			affected_count = 0;

			//iteration_time[ranking] = (GetCurrentTimeSec() - updating_acc_time);
			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			time_this_iteration = GetCurrentTimeSec() - time_this_iteration;
			iteration_time[ranking] = time_this_iteration;
		}


		cout << "picked count:" << picked_count << endl;
		cout << "picked as tree covers:" << used_as_cover << endl;
		/*for (int v = 0; v < V; ++v)
		if (has_picked[v] == false) {
		cout << "wrong 2" << endl;
		cout << picked_value[v] << endl;
		}*/
		
		//// FIX
		//for (int v = 0; v < numOfVertices; ++v){
		//	for (int i = 0; i < CNT; ++i) {
		//		if (cnt16[v][i] != 0)
		//			cout << "mi?" << endl;
		//	}
		//}

		if (stop_flag == true) return;

		double labeling_acc_time = GetCurrentTimeSec();
		// Order the remaining vertices with degree, if any.
		if (ranking != numOfVertices) {
			vector<NodeID> degree(numOfVertices);
			vector<NodeID> r_degree(numOfVertices);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				/*vector<NodeID>& adj_v = graph.adj[v];
				for (NodeID i = 0; i < adj_v.size(); ++i) {
				if (has_picked[adj_v[i]] != true)
				degree[v]++;
				}*/

				for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid) {
					NodeID w = graph.edges[eid];
					if (has_picked[w] != true)
						degree[v]++;
				}

				if (DIRECTED_FLAG == true) {
					/*vector<NodeID>& r_adj_v = graph.r_adj[v];
					for (NodeID i = 0; i < r_adj_v.size(); ++i) {
					if (has_picked[r_adj_v[i]] != true)
					r_degree[v]++;
					}*/
					for (EdgeID eid = graph.r_vertices[v]; eid < graph.r_vertices[v + 1]; ++eid) {
						NodeID w = graph.r_edges[eid];
						if (has_picked[w] != true)
							r_degree[v]++;
					}
				}
			}

			vector<pair<NodeID, NodeID> > degree_v;
			degree_v.reserve(numOfVertices - ranking);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				if (DIRECTED_FLAG == true)
					degree_v.push_back(make_pair(degree[v] * r_degree[v], v));
				else
					degree_v.push_back(make_pair(degree[v], v));
			}
			sort(degree_v.rbegin(), degree_v.rend());

			int count = 0;
			for (; ranking < numOfVertices; ++ranking) {
				NodeID tpv = degree_v[count++].second;
				inv[ranking] = tpv;
				rank[tpv] = ranking;
				labeling_source_bfs(tpv, que, vis, dst_r, r_dst_r, usd, ranking, tmp_idx, r_tmp_idx, tmp_idx_token_parents, r_tmp_idx_token_parents, graph);
			}
		}
		labeling_time += (GetCurrentTimeSec() - labeling_acc_time);

		 
		converttokens(tmp_idx,  tmp_idx_token_parents, false);
		//converttokens(tmp_idx,  tmp_idx_token_parents, r_tmp_idx, r_tmp_idx_token_parents);
		cout << clabels.numOfTokens << " Tokens in total" << endl;
		cout << (double)children_size / (double) clabels.numOfTokens << " average children number" << endl;
		
		if (DIRECTED_FLAG == true) {
			converttokens(r_tmp_idx,  r_tmp_idx_token_parents, true);
			cout << clabels.r_numOfTokens << " Tokens in total" << endl;
			cout << (double)r_children_size / (double) clabels.r_numOfTokens << " average children number" << endl;
		}
		
		
		// Finish the constructions.
		/*
		for (NodeID v = 0; v < numOfVertices; ++v) {
			if (DIRECTED_FLAG == false) {
				
				NodeID k = tmp_idx[v].first.size();
				labels.index_[v].spt_v.resize(k);
				labels.index_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();

				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();

			}
			else {
				NodeID k = tmp_idx[v].first.size();
				dlabels.index_[v].spt_v.resize(k);
				dlabels.index_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();

				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();

				k = r_tmp_idx[v].first.size();
				dlabels.bindex_[v].spt_v.resize(k);
				dlabels.bindex_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_v[i] = r_tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_d[i] = r_tmp_idx[v].second[i];
				r_tmp_idx[v].first.clear();
				r_tmp_idx[v].second.clear();

				r_tmp_idx[v].first.shrink_to_fit();
				r_tmp_idx[v].second.shrink_to_fit();
			}
		}*/

		num_of_trees = st_roots.size();

		double small_tree_num = 0;

		for (int i = 0; i < num_of_trees; ++i) {
			if (trees_pointers[i] != -1) {
				small_tree_num++;
				//small_tree_num += strees[trees_pointers[i]].size();
			}
		}

		total_resources = small_tree_num;
		/*for (int i = 0; i < cnt16.size(); ++i){
		for (int j = 0; j < cnt16[i].size(); ++j)
		if (cnt16[i][j] != 0)
		cout << "aaaa " << cnt16[i][j] << endl;
		}*/

	}
	
	Betweenness_Ordering(NodeID k, double beta, Graph& graph, NodeID stopv, long long bound_times = 10, long long trans_times = 8) {

		double start_time = GetCurrentTimeSec();

		if (k < CNT) {
			cout << " Should pick at least 16 seed!" << endl;
			return;
		}
		iteration_time.resize(numOfVertices);

		updating_time = 0;
		labeling_time = 0;
		adding_time = 0;
		selecting_time = 0;
		init_time = 0;
		u1_time = 0;
		u2_time = 0;
		u3_time = 0;
		num_of_trees = 0;
		total_resources = 0;
		get_small_parent_time = 0;

		inv.resize(numOfVertices, -1);
		rank.resize(numOfVertices, -1);

		if (DIRECTED_FLAG == false) {
			labels.index_.resize(numOfVertices);
		}
		else { //if (DIRECTED_FLAG == true) {
			dlabels.index_.resize(numOfVertices);
			dlabels.bindex_.resize(numOfVertices);
		}

		// Resource Bounding.
		NodeID new_tree_arcs = 0;
		NodeID labeled_arcs = 0;
		alive_resources = 0;
		bounded_resources = (long long)10 * (long long)k * (long long)numOfVertices;

		bool stop_flag = false;

		vector<NodeID> st_roots; // The root vertices of the sampled trees.
		vector<NodeID> st_sizes; // The size of the sampled trees.

		NodeID size_thr = numOfVertices < 8 ? numOfVertices : numOfVertices / 8; // When a tree shrink to a size smaller than the size_thr, it should be converted into a small_tree representation.
																				 //	vector<vector<int> > cover;


		init_time = GetCurrentTimeSec();
		// Building the initial k trees.
	//	srand((unsigned)time(NULL));
		srand(0x728f9d0);
		

		vector<NodeID> descendants_parents(numOfVertices, -1);

		ltrees.reserve(numOfVertices);
		strees.reserve(numOfVertices);
		trees_pointers.reserve(numOfVertices);
		st_roots.reserve(numOfVertices);
		st_sizes.reserve(numOfVertices);
		st_roots.resize(k);
		st_sizes.resize(k);

		ltrees.resize(k);
		strees.clear();
		trees_pointers.resize(k, -1); // Point the sampling trees to its small tree representations.
									  //	cover.resize(k, vector<int>(V, 0));

		vector<NodeID> ascendants;
		ascendants.reserve(numOfVertices);
		vector<NodeID> descendants;
		descendants.reserve(numOfVertices);
		vector<NodeID> que(numOfVertices);
		vector<bool> vis(numOfVertices);
		vector<bool> has_picked(numOfVertices);
		vector<vector<cover_value_type> > cnt16(numOfVertices, vector<cover_value_type>(CNT));
		//		vector<bool> affected(V, true);
		vector<bool> picked_as_root(numOfVertices, false);
		vector<bool> trivial_root(numOfVertices, false);
		int picked_count = 0;
		int used_as_cover = 0;
		vector<cover_value_type> sum(numOfVertices);
		vector<bool> usd(numOfVertices, false);
		vector<cover_value_type> cover_tid(numOfVertices, 0);

		//unordered_set<NodeID> affected_v; //GOOGLE
		google::dense_hash_set<NodeID> affected_v;
		affected_v.set_empty_key(numOfVertices + 1);
		affected_v.set_deleted_key(numOfVertices + 2);

		vector<bool> affected_list(numOfVertices, false);
		int affected_count = 0;
		bool switch_set = false;

		//vector<vector<int> > cover_record(V, vector<int>(V, 0));

		//vector<unordered_set<int> > reverse_index(0); //GOOGLE
//		vector<unordered_set<int> > reverse_index(0);
		vector<google::dense_hash_set<NodeID> > reverse_index;
		reverse_index.reserve(numOfVertices);

		// Preparing basic structure for pl algorithm.
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT))); // Backward labels.

													 //	int best_sum = 0;
													 //	NodeID best_sum_vid = 0;
													 //vector<bool> lazy_updates(V, false);

													 // Building the initial k sampled trees. 
		for (int i = 0; i < k; ++i) {
			NodeID r = rand() % numOfVertices;
			if (picked_as_root[r] == true) { --i; continue; }

			st_roots[i] = r;
			large_tree &ltree_i = ltrees[i];
			picked_as_root[r] = true;
			picked_count++;
			used_as_cover++;

			// Array Representation
			//build_large_tree_bfs(r, ltree_i, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph.adj, 0, numOfVertices);
			build_large_tree_bfs(r, ltree_i, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph, 0, numOfVertices);

			st_sizes[i] = descendants.size();
			alive_resources += (long long)numOfVertices;

			// Turn large tree to small tree.
			// 1. Build new small tree with ltree_i and add it to strees;
			// 2. Set tree_pointers[i] to the position of the strees;
			// 3. Clear the large tree ltree_i by setting its size to 0.
			if (st_sizes[i] < size_thr)
				shrink_tree(i, st_roots[i], ltrees[i]);


			// Calculate the subtree size of each nodes			.
			//	vector<int>& cover_i = cover[i];

			bool isLarge = false;
			int stree_idx = trees_pointers[i];
			if (stree_idx == -1)
				isLarge = true;

			small_tree& stree_i = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;

			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];

				++cover_tid[q]; // Increment because it hitsthe r-q path.
				++cnt16[q][i%CNT];
				//++sum[q];
				//	if (sum[q] > best_sum) {
				//		best_sum = sum[q];
				//		best_sum_vid = q;
				//	}
				////	cover_record[i][q]++;

				int parent_node;
				if (isLarge)
					parent_node = get_large_parent(q, ltree_i);
				else
					parent_node = get_small_parent(q, stree_i);


				if (parent_node > -1) {
					cover_tid[parent_node] += cover_tid[q];
					cnt16[parent_node][i%CNT] += cover_tid[q];
					//sum[parent_node] += cover_i[q];
					//	cover_record[i][parent_node] += cover_record[i][q];
				}
			}	
			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];
				cover_tid[q] = 0;
			}
		}


		//cnt16_aggreate(cnt16, has_picked, cover);

		vector<cover_value_type> picked_value;
		benchmark::heap<2, cover_value_type, NodeID> value_heap(numOfVertices);
		bool switch_small = false;

		//cnt16_2max(picked_value, sum, has_picked,affected, cnt16, switch_small, value_heap);
		//		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap, best_sum_vid, lazy_updates);
		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);

		init_time = GetCurrentTimeSec() - init_time;


		cout << "Finish initial sampling " << init_time << "secs." << endl;

		//for (NodeID v = 0; v < V; ++v)	value_heap.update(v, -picked_value[v]); // Min-heap.

		double max_value = -1;


		// Perform V iterations to rank all vertices.
		NodeID ranking = 0;
		for (; ranking < numOfVertices; ++ranking) {
		/*	if (ranking % (numOfVertices / 10) == 0)
				cout << "picking the " << ranking << "-th orders." << endl;*/

			double time_this_iteration = GetCurrentTimeSec();

			//best_sum = 0;

			NodeID pv = 0;
			cover_value_type pv_cover = -1;

			if (ranking == stopv) {
				stop_flag = true;
				cout << "topv labeling time:" << labeling_time << endl;
				break;
			}
			double selecting_acc_time = GetCurrentTimeSec();
			// In the beginning, pick the best cover vertex by a linear scan.
			if (switch_small == false) {
				double u1_acc_time = GetCurrentTimeSec();
				//while (true) {
				for (NodeID v = 0; v < numOfVertices; ++v) {
					if (has_picked[v] == true) continue;
					if (pv_cover == picked_value[v]) {
						if (sum[pv] < sum[v]) {
							pv_cover = picked_value[v];
							pv = v;
						}
					}
					else if (pv_cover < picked_value[v]) {
						pv_cover = picked_value[v];
						pv = v;
					}
				}
				/*	if (lazy_updates[pv] == true) {
				NodeID max1 = 0;
				NodeID max2 = 0;
				const vector<int>& cnt16_pv = cnt16[pv];
				int& picked_value_v = picked_value[pv];

				int sum_v = sum[pv];
				for (int i = 0; i < CNT; ++i) {
				if (cnt16_pv[i] > max2)
				max2 = cnt16_pv[i];
				if (max2 > max1) {
				int tmp = max1;
				max1 = max2;
				max2 = tmp;
				}
				}
				picked_value_v = sum_v - max1 - max2;

				lazy_updates[pv] = false;
				}
				else
				break;*/
				//}

				if (sum[pv] < numOfVertices / 8 && switch_small == false) {
					cout << ranking << endl;
					cout << pv_cover << endl;
					cout << sum[pv] << endl;
					cout << picked_count << endl;
					cout << "Switch to max-heap representations." << endl;
					switch_small = true;
				}
				u1_time += (GetCurrentTimeSec() - u1_acc_time);
			}
			else {// Later, utilizing the max-heap directly
				double u2_acc_time = GetCurrentTimeSec();
				//while (true) {
				value_heap.extract_min(pv, pv_cover);
				pv_cover = -pv_cover;
				/*	if (lazy_updates[pv] == true) {
				NodeID max1 = 0;
				NodeID max2 = 0;
				const vector<int>& cnt16_pv = cnt16[pv];
				int& picked_value_v = picked_value[pv];

				int sum_v = sum[pv];
				for (int i = 0; i < CNT; ++i) {
				if (cnt16_pv[i] > max2)
				max2 = cnt16_pv[i];
				if (max2 > max1) {
				int tmp = max1;
				max1 = max2;
				max2 = tmp;
				}
				}
				picked_value_v = sum_v - max1 - max2;
				value_heap.update(pv, -picked_value_v);
				lazy_updates[pv] = false;
				}
				else
				break;
				}*/
				u2_time += (GetCurrentTimeSec() - u2_acc_time);

			}

			selecting_time += (GetCurrentTimeSec() - selecting_acc_time);

			//// FIX
			//if (switch_small == true && reverse_index.size() > 0 && sum[pv] <= reverse_index[pv].size()) {
			//	cout << "I think we can stop here " << GetCurrentTimeSec() - start_time << " ranking:" << ranking << endl;
			//	break;
			//}

			has_picked[pv] = true;
			double labeling_acc_time = GetCurrentTimeSec();
			labeled_arcs = labeling_source_bfs(pv, que, vis, dst_r, usd, ranking, tmp_idx, r_tmp_idx, graph);
			labeling_time += (GetCurrentTimeSec() - labeling_acc_time);
			inv[ranking] = pv;
			rank[pv] = ranking;

			//cout << ranking << "\t" << labeled_arcs << endl;


			int current_trees_num = st_roots.size();

			// 4 steps in updating:
			double updating_acc_time = GetCurrentTimeSec();
			if (switch_small == false) { // We still iterate all alive trees to calculate the cover.	
				for (int tid = 0; tid < current_trees_num; ++tid) {

					// This tree has already be taken because its root has been selected.
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					// pv has already been covered in tree tid.
					int pvpid = 0;
					if (isLarge)
						pvpid = get_large_parent(pv, ltree_tid);
					else
						pvpid = get_small_parent(pv, stree_tid);
					if (pvpid == -1 && st_roots[tid] != pv) // pv has been covered in tree tid 
						continue;


					// Get the descendants and ascendants of pv in tree tid.

					// Array Representation
					//// 1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;				 
					//get_descendants(tid, pv, descendants, descendants_parents, que, vis, graph.adj, isLarge, stree_tid, ltree_tid);

					//get_ascendants(tid, pv, ascendants, que, vis, graph.adj, graph.r_adj, isLarge, stree_tid, ltree_tid);
				
					get_descendants(tid, pv, descendants, descendants_parents, que, vis, graph, isLarge, stree_tid, ltree_tid);

					get_ascendants(tid, pv, ascendants, que, vis, graph, isLarge, stree_tid, ltree_tid);

					NodeID subtract_size = descendants.size();
					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
												//cover[tid][av] -= subtract_size;
						cnt16[av][tid%CNT] -= subtract_size;
						//sum[av] -= subtract_size;
						//	cover_record[tid][av] -= subtract_size;
						//	if (cover_record[tid][av] < 0)
						//	cout << "1: " << cover_record[tid][av] << endl;
						if (affected_list[av] == false) {
							affected_list[av] = true;
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
						}
					}


					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
				//	vector<int> cover_tid(numOfVertices, 0);
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];

						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//sum[dv] -= cover_tid[dv];
						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];

					/*	if (cnt16[dv][tid%CNT] < 0)
							cout << "a6 " << cnt16[dv][tid%CNT] << endl;*/

						//	cover_record[tid][dv] -= cover_tid[dv];
						//	if (cover_record[tid][dv] < 0) {
						//	cout << "2: " << cover_record[tid][dv]  << " " << cover_tid[dv]  << " " << tid << " " << dv <<  endl;
						//	}

						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else 
							remove_small_node(dv, stree_tid);

						
						//	remove_node(tid, dv);
						//affected[dv] = true;
						if (affected_list[dv] == false) {
							affected_count++;
							affected_list[dv] = true;
							if (switch_set)
								affected_v.insert(dv);
						}
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;
					}


					// Drop the size of subtree.
					st_sizes[tid] -= subtract_size;
					if (trees_pointers[tid] != -1) alive_resources -= (cover_value_type)subtract_size;
				//	if (st_sizes[tid] < 0) {
				//		cout << "a5" << st_sizes[tid] << endl;
				//	}
					if (st_sizes[tid] > 0) {
						// 3. Convert shrink trees to small trees representataions.				
						if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
							shrink_tree(tid, st_roots[tid], ltrees[tid]);
						}
					}

					if (st_roots[tid] == pv)
						if (isLarge) {
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
							//ltrees.shrink_to_fit();
						}
						else {
							stree_tid.clear();
							stree_tid.resize(0);
							//strees.shrink_to_fit();
						}

				}
			}
			else { // We use the reverse index to update trees.

				   // Build the reverse_index first (one time).
				if (reverse_index.size() == 0) {
					reverse_index.resize(numOfVertices);
					for (NodeID v = 0; v < numOfVertices; ++v) {
						reverse_index[v].set_empty_key(numOfVertices + 1);
						reverse_index[v].set_deleted_key(numOfVertices + 2);
					}
					for (int tid = 0; tid < current_trees_num; ++tid) {

						if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue; //This whole tree has already been covered. 

						bool isLarge = false;
						int stree_idx = trees_pointers[tid];
						if (stree_idx == -1)
							isLarge = true;

						small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
						large_tree& ltree_tid = ltrees[tid];

						for (int v = 0; v < numOfVertices; ++v) {

							// Get the descendants and ascendants of pv in tree tid.

							//	if (cover[tid].find(pv) == cover[tid].end()) continue; // This pv has been covered in this tree before.
							//	if (cover[tid][pv] <= 0) continue; // modified-5-14
							int vpid = 0;
							if (isLarge)
								vpid = get_large_parent(v, ltree_tid);
							else
								vpid = get_small_parent(v, stree_tid);

							// To testify whether v is still uncovered in tree tid
							// v has parent vpid != -1, but if v is the root of tree tid, it still can be uncovered.
							if ((vpid != -1 || st_roots[tid] == v))
								reverse_index[v].insert(tid);
							//if (cover[tid].find(v) != cover[tid].end())
							//	reverse_index[v].insert(tid);
						}
					}
				}

				double u4_acc_time = GetCurrentTimeSec();

				//	double u4_acc_time = GetCurrentTimeSec();
				// Access the alive trees directly by the reverse index.
				//				cout << reverse_index[pv].size() << " " << sum[pv] << endl;

				
//				for (unordered_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) { // GOOGLE
				for (google::dense_hash_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) {


					int tid = *it;
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					// Get the descendants and ascendants of pv in tree tid.

					//	if (cover[tid].find(pv) == cover[tid].end()) continue; // This pv has been covered in this tree before.
					//	if (cover[tid][pv] <= 0) continue; // modified-5-14

					int pvpid = 0;
					if (isLarge)
						pvpid = get_large_parent(pv, ltree_tid);
					else {
						//	double get_small_parent_incre_time = GetCurrentTimeSec();
						pvpid = get_small_parent(pv, stree_tid);
						//	get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}

					if (pvpid == -1 && st_roots[tid] != pv)
						continue;

					// Array Representation
					//// 1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;
					//get_descendants(tid, pv, descendants, descendants_parents, que, vis, graph.adj, isLarge, stree_tid, ltree_tid);

					//get_ascendants(tid, pv, ascendants, que, vis, graph.adj, graph.r_adj, isLarge, stree_tid, ltree_tid);

					get_descendants(tid, pv, descendants, descendants_parents, que, vis, graph, isLarge, stree_tid, ltree_tid);

					get_ascendants(tid, pv, ascendants, que, vis, graph, isLarge, stree_tid, ltree_tid);

					NodeID subtract_size = descendants.size();

					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
												//						cover[tid][av] -= subtract_size;

						cnt16[av][tid%CNT] -= subtract_size;
						//sum[av] -= subtract_size;

						//	cover_record[tid][av] -= subtract_size;
						//	if (cover_record[tid][av] < 0) {
						//	cout << "4: " << cover_record[tid][av] << " " << subtract_size << " " << tid << " " << endl;
						//	}

						//affected[av] = true;
						/*if (cnt16[av][tid%CNT] < 0)
						cout << " a2 " << cnt16[av][tid%CNT]  << " " << av << " " << st_roots[tid] << " " << tid << endl;*/
						if (affected_list[av] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
							affected_list[av] = true;
						}
					}

					//double u2_acc_time = GetCurrentTimeSec();
					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					//vector<int> cover_tid(numOfVertices, 0);
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
					//	sum[dv] -= cover_tid[dv];
						//cover[tid].erase(dv);
						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];

						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else 
							remove_small_node(dv, stree_tid);
							
						//	remove_node(tid, dv);
						//affected[dv] = true;

						//	cover_record[tid][dv] -= cover_tid[dv];
						//if (cover_record[tid][dv] < 0)
						//	cout << "3: " << cover_record[tid][dv] << " " << cover_tid[dv] << " " << di << " " << descendants.size() << " " << picked_as_root[st_roots[tid]] << endl;

						/*if (cnt16[dv][tid%CNT]  < 0)
						cout << " a3 " << cnt16[dv][tid%CNT] << " " << dv << " " << st_roots[tid] << " " << tid << endl;
						*/
						if (dv != pv)
							reverse_index[dv].erase(tid);
						//in_tree_count[dv]--;


						if (affected_list[dv] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(dv);
							affected_list[dv] = true;
						}
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;
					}

					if (st_roots[tid] == pv)
						if (isLarge) {
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
						}
						else {
							stree_tid.clear();
							stree_tid.resize(0);
						}

					st_sizes[tid] -= subtract_size;
				//	if (st_sizes[tid]  < 0)
				//		cout << " a4 " << st_sizes[tid] << endl;

					if (trees_pointers[tid] != -1) alive_resources -= subtract_size;


					// 3. Convert shrink trees to small trees representataions.	
					if (st_sizes[tid] > 0) {
						if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
							shrink_tree(tid, st_roots[tid], ltrees[tid]);
						}
					}

				}
				reverse_index[pv].clear();
				reverse_index[pv].resize(0);
				u4_time += (GetCurrentTimeSec() - u4_acc_time);
			}

			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			// 4. Sample new trees based on.
			// If we sample a new tree, we need to dynamically increment the following data structures:
			// 1) ltrees
			// 2) trees_pointers (-1 or pointer to strees)
			// 3) cover
			// 4) st_roots
			// 5) st_sizes 
			// 6) strees (maybe)
			// two constraints: cl > ct and 10kn total vertices
			double adding_acc_time = GetCurrentTimeSec();
			new_tree_arcs = 0;

			//cout << "ranking:" << ranking << endl;
			//cout << "pick cover:" << pv_cover << endl;
			//cout << "pick count:" << picked_count << endl;
			//cout << "labeled_arcs:" << labeled_arcs << endl;
			//cout << "new_tree_arcs:" << new_tree_arcs << endl;
			//cout << "alive_resources:" << alive_resources << endl;
			//cout << "bounded_resources:" << bounded_resources << endl;
			//cout << "#tree:" << cover.size() << endl;
			int how_many = 0;
			//cout << "beta*labeld_arcs: " << beta * labeled_arcs << " new tree _arcs: " << new_tree_arcs << " alive " << alive_resources << " bounded " << bounded_resources << endl;
			NodeID sam_this_round = 0;
			while (beta * labeled_arcs > new_tree_arcs && alive_resources < bounded_resources) {

				if (picked_count == numOfVertices) break;



				NodeID r = rand() % numOfVertices;
				if (picked_as_root[r] == true) continue;
				if (trivial_root[r] == true) continue;
				int tmp_tree_arcs = 0;

				large_tree new_ltree;
				small_tree new_stree;
				//GOOGLE DENSE MAP
				new_stree.set_empty_key(numOfVertices + 1);
				new_stree.set_deleted_key(numOfVertices + 2);

				//	BUILD_SMALL_TREE_FLAG = false;
				bool isLarge = false;
				if (!BUILD_SMALL_TREE_FLAG) {

					double u5_acc_time = GetCurrentTimeSec();
					// Array Representation
					//tmp_tree_arcs = build_large_tree_bfs(r, new_ltree, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph.adj, new_tree_arcs, beta * labeled_arcs);
					tmp_tree_arcs = build_large_tree_bfs(r, new_ltree, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph, new_tree_arcs, beta * labeled_arcs);
					u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;

					if (descendants.size() <size_thr)
						BUILD_SMALL_TREE_FLAG = true;

					if (descendants.size() < 2) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					isLarge = true;
					how_many++;
					used_as_cover++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					st_roots.push_back(r);
					trees_pointers.push_back(-1);
					st_sizes.push_back(descendants.size());

					//new_tree_arcs += descendants.size();
					alive_resources += numOfVertices;

					// Turn large tree to small tree.
					// 1. Build new small tree with ltree_i and add it to strees;
					// 2. Set tree_pointers[i] to the position of the strees;
					// 3. Clear the large tree ltree_i by setting its size to 0.
					int i = st_roots.size() - 1;
					if (st_sizes[i] < size_thr) {
						shrink_tree(i, st_roots[i], new_ltree);
						isLarge = false;
					}
					ltrees.push_back(new_ltree);
				}
				else {//Build Small Tree Directly
					double u5_acc_time = GetCurrentTimeSec();
					// Array Representation
					//	tmp_tree_arcs = build_small_tree_bfs(r, new_stree, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph.adj, new_tree_arcs, beta * labeled_arcs);
					tmp_tree_arcs = build_small_tree_bfs(r, new_stree, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph, new_tree_arcs, beta * labeled_arcs);
					u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;

					if (descendants.size() >= size_thr)
						BUILD_SMALL_TREE_FLAG = false;

					if (descendants.size() < 2) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					isLarge = false;

					how_many++;
					used_as_cover++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					st_roots.push_back(r);
					trees_pointers.push_back(strees.size()); //we assume it is inserted.
					st_sizes.push_back(descendants.size());

					alive_resources += new_stree.size();
					//strees.push_back(new_stree);
					//// new_ltree is empty.
					//ltrees.push_back(new_ltree);

					int i = st_roots.size() - 1;
					if (st_sizes[i] >= size_thr) {
						enlarge_tree(i, st_roots[i], new_stree, new_ltree);
						new_stree.clear();
						new_stree.resize(0);
						isLarge = true;
					}
					else {
						strees.push_back(new_stree);
						ltrees.push_back(new_ltree);
					}
				}

				int i = st_roots.size() - 1;

				int numOfStrees = strees.size();
				small_tree &new_stree_ref = numOfStrees != 0 ? strees[strees.size() - 1] : empty_small_tree;

				// Calculate the subtree size of each nodes.
				//				unordered_map<int, int> cover_new;
			//	vector<cover_value_type> cover_new(numOfVertices, 0); // modified-5-14

				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];

					++cover_tid[q]; // Increment because it hitsthe r-q path. 
					++cnt16[q][i%CNT];
					//++sum[q];
					/*	if (sum[q] > best_sum) {
					best_sum = sum[q];
					best_sum_vid = q;
					}*/
					//affected[q] = true;
					//	cover_record[i][q]++;
					if (affected_list[q] == false) {
						affected_count++;
						if (switch_set)
							affected_v.insert(q);
						affected_list[q] = true;
					}

					if (switch_small == true) {
						reverse_index[q].insert(i);
					}

					int parent_node = -1;

					if (isLarge)
						parent_node = get_large_parent(q, new_ltree);
					else
						parent_node = get_small_parent(q, new_stree_ref);

					if (parent_node > -1) {
						cover_tid[parent_node] += cover_tid[q];
						cnt16[parent_node][i%CNT] += cover_tid[q];
					//	sum[parent_node] += cover_new[q];
						//	cover_record[i][parent_node] += cover_record[i][q];
					}
				}
				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];
					cover_tid[q] = 0;
				}

				sam_this_round++;
			}

			//cout << ranking << "\t" << sam_this_round << endl;


			adding_time += (GetCurrentTimeSec() - adding_acc_time);
			//cout << ranking << ":  beta:" << beta << " beta * labeled_arcs " << beta * labeled_arcs << " new_tree_arcs " << new_tree_arcs << " how many " << how_many << endl;

			//double u3_acc_time = GetCurrentTimeSec();

			updating_acc_time = GetCurrentTimeSec();

			// Build the cover value heap for one time.
			if (switch_small == true && value_heap.empty())
				for (int nid = 0; nid < numOfVertices; ++nid)
					if (has_picked[nid] == false)
						value_heap.update(nid, -picked_value[nid]);

			//	cnt16_2max(picked_value, sum, has_picked, affected, cnt16, switch_small, value_heap);		
			//double u3_acc_time = GetCurrentTimeSec();

			//			cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap, best_sum_vid, lazy_updates);
			cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);
			//u3_time += (GetCurrentTimeSec() - u3_acc_time);

			if (affected_count < numOfVertices / 8)
				switch_set = true;
			else
				switch_set = false;
			affected_count = 0;

			//iteration_time[ranking] = (GetCurrentTimeSec() - updating_acc_time);
			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			time_this_iteration = GetCurrentTimeSec() - time_this_iteration;
			iteration_time[ranking] = time_this_iteration;
		}


		cout << "picked count:" << picked_count << endl;
		cout << "picked as tree covers:" << used_as_cover << endl;
		/*for (int v = 0; v < V; ++v)
		if (has_picked[v] == false) {
		cout << "wrong 2" << endl;
		cout << picked_value[v] << endl;
		}*/
		
		//// FIX
		//for (int v = 0; v < numOfVertices; ++v){
		//	for (int i = 0; i < CNT; ++i) {
		//		if (cnt16[v][i] != 0)
		//			cout << "mi?" << endl;
		//	}
		//}

		if (stop_flag == true) return;

		double labeling_acc_time = GetCurrentTimeSec();
		// Order the remaining vertices with degree, if any.
		if (ranking != numOfVertices) {
			vector<NodeID> degree(numOfVertices);
			vector<NodeID> r_degree(numOfVertices);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				/*vector<NodeID>& adj_v = graph.adj[v];
				for (NodeID i = 0; i < adj_v.size(); ++i) {
				if (has_picked[adj_v[i]] != true)
				degree[v]++;
				}*/

				for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid) {
					NodeID w = graph.edges[eid];
					if (has_picked[w] != true)
						degree[v]++;
				}

				if (DIRECTED_FLAG == true) {
					/*vector<NodeID>& r_adj_v = graph.r_adj[v];
					for (NodeID i = 0; i < r_adj_v.size(); ++i) {
					if (has_picked[r_adj_v[i]] != true)
					r_degree[v]++;
					}*/
					for (EdgeID eid = graph.r_vertices[v]; eid < graph.r_vertices[v + 1]; ++eid) {
						NodeID w = graph.r_edges[eid];
						if (has_picked[w] != true)
							r_degree[v]++;
					}
				}
			}

			vector<pair<NodeID, NodeID> > degree_v;
			degree_v.reserve(numOfVertices - ranking);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				if (DIRECTED_FLAG == true)
					degree_v.push_back(make_pair(degree[v] * r_degree[v], v));
				else
					degree_v.push_back(make_pair(degree[v], v));
			}
			sort(degree_v.rbegin(), degree_v.rend());

			int count = 0;
			for (; ranking < numOfVertices; ++ranking) {
				NodeID tpv = degree_v[count++].second;
				inv[ranking] = tpv;
				rank[tpv] = ranking;
				labeling_source_bfs(tpv, que, vis, dst_r, usd, ranking, tmp_idx, r_tmp_idx, graph);
			}
		}
		labeling_time += (GetCurrentTimeSec() - labeling_acc_time);

		// Finish the constructions.
		for (NodeID v = 0; v < numOfVertices; ++v) {
			if (DIRECTED_FLAG == false) {
				NodeID k = tmp_idx[v].first.size();
				labels.index_[v].spt_v.resize(k);
				labels.index_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();

				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();

			}
			else {
				NodeID k = tmp_idx[v].first.size();
				dlabels.index_[v].spt_v.resize(k);
				dlabels.index_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();

				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();

				k = r_tmp_idx[v].first.size();
				dlabels.bindex_[v].spt_v.resize(k);
				dlabels.bindex_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_v[i] = r_tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_d[i] = r_tmp_idx[v].second[i];
				r_tmp_idx[v].first.clear();
				r_tmp_idx[v].second.clear();

				r_tmp_idx[v].first.shrink_to_fit();
				r_tmp_idx[v].second.shrink_to_fit();
			}
		}

		num_of_trees = st_roots.size();

		double small_tree_num = 0;

		for (int i = 0; i < num_of_trees; ++i) {
			if (trees_pointers[i] != -1) {
				small_tree_num++;
				//small_tree_num += strees[trees_pointers[i]].size();
			}
		}

		total_resources = small_tree_num;
		/*for (int i = 0; i < cnt16.size(); ++i){
		for (int j = 0; j < cnt16[i].size(); ++j)
		if (cnt16[i][j] != 0)
		cout << "aaaa " << cnt16[i][j] << endl;
		}*/

	}

	Betweenness_Ordering(NodeID k, double beta, Graph& graph, NodeID stopv, long long bound_times, long long trans_times, double mem_bound) {

		double start_time = GetCurrentTimeSec();

		if (k < CNT) {
			cout << " Should pick at least 16 seed!" << endl;
			return;
		}
		iteration_time.resize(numOfVertices);

		bound_times = 10;
		trans_times = 8;
		updating_time = 0;
		labeling_time = 0;
		adding_time = 0;
		selecting_time = 0;
		init_time = 0;
		u1_time = 0;
		u2_time = 0;
		u3_time = 0;
		num_of_trees = 0;
		total_resources = 0;
		get_small_parent_time = 0;
		
		inv.resize(numOfVertices, -1);
		rank.resize(numOfVertices, -1);
		 
		if (DIRECTED_FLAG == false) {
			labels.index_.resize(numOfVertices);
		}
		else { //if (DIRECTED_FLAG == true) {
			dlabels.index_.resize(numOfVertices);
			dlabels.bindex_.resize(numOfVertices);
		}

		// Resource Bounding.
		NodeID new_tree_arcs = 0;
		NodeID labeled_arcs = 0;
		long long alive_memory = 0;

		bounded_resources = mem_bound * 1024 * 1024 * 1024;

		bool stop_flag = false;

		vector<NodeID> st_roots; // The root vertices of the sampled trees.
		vector<NodeID> st_sizes; // The size of the sampled trees.
		st_roots.reserve(numOfVertices);
		st_sizes.reserve(numOfVertices);
		ltrees.reserve(numOfVertices);
		strees.reserve(numOfVertices);
		trees_pointers.reserve(numOfVertices);
		alive_memory += numOfVertices * sizeof(NodeID) * 2;
		alive_memory += numOfVertices * sizeof(large_tree);
		alive_memory += numOfVertices * sizeof(small_tree);
		alive_memory += numOfVertices * sizeof(NodeID);

		NodeID size_thr = numOfVertices < 8 ? numOfVertices : numOfVertices / 8; // When a tree shrink to a size smaller than the size_thr, it should be converted into a small_tree representation.
																				 //	vector<vector<int> > cover;


		init_time = GetCurrentTimeSec();
		// Building the initial k trees.
		//srand((unsigned)time(NULL));
		srand(0x728f9d0);

		vector<NodeID> descendants_parents(numOfVertices, -1);
		alive_memory += numOfVertices * sizeof(NodeID);


		st_roots.resize(k);
		st_sizes.resize(k);

		ltrees.resize(k);
		strees.clear();
		trees_pointers.resize(k, -1); // Point the sampling trees to its small tree representations.
									  //	cover.resize(k, vector<int>(V, 0));
		

		vector<NodeID> que(numOfVertices);
		vector<bool> vis(numOfVertices);

		vector<NodeID> ascendants;
		ascendants.reserve(numOfVertices);
		vector<NodeID> descendants;
		descendants.reserve(numOfVertices);
		vector<bool> has_picked(numOfVertices);
		vector<vector<cover_value_type> > cnt16(numOfVertices, vector<cover_value_type>(CNT));

		alive_memory += numOfVertices * sizeof(NodeID);
		alive_memory += numOfVertices * sizeof(NodeID);
		alive_memory += numOfVertices * sizeof(bool);
		alive_memory += numOfVertices * CNT * sizeof(cover_value_type);

		//		vector<bool> affected(V, true);
		vector<bool> picked_as_root(numOfVertices, false);
		vector<bool> trivial_root(numOfVertices, false);
		int picked_count = 0;
		int used_as_cover = 0;
		vector<cover_value_type> sum(numOfVertices);
		vector<bool> usd(numOfVertices, false);
		vector<cover_value_type> cover_tid(numOfVertices, 0);
		alive_memory += numOfVertices * sizeof(bool);
		alive_memory += numOfVertices * sizeof(bool);
		alive_memory += numOfVertices * sizeof(cover_value_type);
		alive_memory += numOfVertices * sizeof(bool);
		alive_memory += numOfVertices * sizeof(cover_value_type);



		//unordered_set<NodeID> affected_v; //GOOGLE
		google::dense_hash_set<NodeID> affected_v;
		affected_v.set_empty_key(numOfVertices + 1);
		affected_v.set_deleted_key(numOfVertices + 2);

		// google dense hash set... how to measure its memory usage?


		vector<bool> affected_list(numOfVertices, false);
		int affected_count = 0;
		bool switch_set = false;
		alive_memory += numOfVertices * sizeof(bool);

		//vector<vector<int> > cover_record(V, vector<int>(V, 0));

		//vector<unordered_set<int> > reverse_index(0); //GOOGLE
		//		vector<unordered_set<int> > reverse_index(0);
		vector<google::dense_hash_set<NodeID> > reverse_index;
		reverse_index.reserve(numOfVertices);
		alive_memory += numOfVertices * sizeof(google::dense_hash_set<NodeID>);

		// google dense hash set... how to measure its memory usage?


		// Preparing basic structure for pl algorithm.
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT))); // Backward labels.

													 //	int best_sum = 0;
													 //	NodeID best_sum_vid = 0;
													 //vector<bool> lazy_updates(V, false);

		vector<cover_value_type> dependencies(numOfVertices, 0);

		cout << "Auxilariy Data Memory Usage " << alive_memory << endl;
		// Building the initial k sampled trees. 
		for (int i = 0; i < k; ++i) {
			NodeID r = rand() % numOfVertices;
			if (picked_as_root[r] == true) { --i; continue; }

			st_roots[i] = r;
			large_tree &ltree_i = ltrees[i];
			picked_as_root[r] = true;
			picked_count++;
			used_as_cover++;

			// Array Representation
			//build_large_tree_bfs(r, ltree_i, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph.adj, 0, numOfVertices);
			build_large_tree_bfs(r, ltree_i, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph, 0, numOfVertices);

			st_sizes[i] = descendants.size();
			alive_memory += numOfVertices * sizeof(NodeID);

			// Turn large tree to small tree.
			// 1. Build new small tree with ltree_i and add it to strees;
			// 2. Set tree_pointers[i] to the position of the strees;
			// 3. Clear the large tree ltree_i by setting its size to 0.
			if (st_sizes[i] < size_thr) {
				shrink_tree(i, st_roots[i], ltrees[i]);
				alive_memory -= numOfVertices * sizeof(NodeID);
				alive_memory += strees[trees_pointers[i]].bucket_count()  *  sizeof(NodeID) * 2;

				/*if (strees[trees_pointers[i]].bucket_count() * sizeof(NodeID) * 2 > numOfVertices * sizeof(NodeID)) {
					cout << "Are you kidding me?" << endl;
					cout << strees[trees_pointers[i]].size() << " : " << strees[trees_pointers[i]].bucket_count() * 2 << endl;
				}*/
			}


			// Calculate the subtree size of each nodes			.
			//	vector<int>& cover_i = cover[i];


			bool isLarge = false;
			int stree_idx = trees_pointers[i];
			if (stree_idx == -1)
				isLarge = true;


			small_tree& stree_i = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;

			alive_memory += descendants.size() * sizeof(NodeID); // The size occupied by reverse index.

			//vector<int> cover_i(numOfVertices, 0);
			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];

				++cover_tid[q]; // Increment because it hitsthe r-q path.
				++cnt16[q][i%CNT];
				++dependencies[q];

				int parent_node;
				if (isLarge)
					parent_node = get_large_parent(q, ltree_i);
				else
					parent_node = get_small_parent(q, stree_i);


				if (parent_node > -1) {
					cover_tid[parent_node] += cover_tid[q];
					cnt16[parent_node][i%CNT] += cover_tid[q];

					dependencies[parent_node] += cover_tid[q];
					//sum[parent_node] += cover_i[q];
					//	cover_record[i][parent_node] += cover_record[i][q];
				}
			}
			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];
				cover_tid[q] = 0;
			}
	
		}


		//cnt16_aggreate(cnt16, has_picked, cover);

		vector<cover_value_type> picked_value;
		benchmark::heap<2, cover_value_type, NodeID> value_heap(numOfVertices);
		bool switch_small = false;

		// Size of value heap.
		alive_memory += numOfVertices * sizeof(NodeID);

		//cnt16_2max(picked_value, sum, has_picked,affected, cnt16, switch_small, value_heap);
		//		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap, best_sum_vid, lazy_updates);
		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);

		init_time = GetCurrentTimeSec() - init_time;


		cout << "Finish initial sampling " << init_time << "secs." << endl;
		cout << "Init Trees Memory Usage " << alive_memory << endl;
		//for (NodeID v = 0; v < V; ++v)	value_heap.update(v, -picked_value[v]); // Min-heap.

		double max_value = -1;


		// Perform V iterations to rank all vertices.
		NodeID ranking = 0;
		for (; ranking < numOfVertices; ++ranking) {
			if (ranking % (numOfVertices / 10) == 0)
				cout << "picking the " << ranking << "-th orders." << endl;

			double time_this_iteration = GetCurrentTimeSec();

			//best_sum = 0;

			NodeID pv = 0;
			cover_value_type pv_cover = -1;
			cover_value_type second_pv_cover = -1;
			NodeID spv = 0;

			if (ranking == stopv) {
				stop_flag = true;
				cout << "topv labeling time:" << labeling_time << endl;	
				break;
			}
			double selecting_acc_time = GetCurrentTimeSec();
			// In the beginning, pick the best cover vertex by a linear scan.
			if (switch_small == false) {
				double u1_acc_time = GetCurrentTimeSec();
				//while (true) {
				for (NodeID v = 0; v < numOfVertices; ++v) {
					if (has_picked[v] == true) continue;

					//if (pv_cover == picked_value[v]) {
					//	if (sum[pv] < sum[v]) {
					//		pv_cover = picked_value[v];
					//		pv = v;
					//	}
					//}
					//else if (pv_cover < picked_value[v]) {
					//	pv_cover = picked_value[v];
					//	pv = v;
					//}
					u1_time += (GetCurrentTimeSec() - u1_acc_time);

					if (second_pv_cover < picked_value[v]) {
						second_pv_cover = picked_value[v];
						spv = v;
						if (pv_cover < second_pv_cover) {
							cover_value_type tmp = second_pv_cover;
							second_pv_cover = pv_cover;
							pv_cover = tmp;
							NodeID tmpv = spv;
							spv = pv;
							pv = tmpv;
						}
						else if (pv_cover == second_pv_cover) {
							if (sum[pv] < sum[spv]) {
								cover_value_type tmp = second_pv_cover;
								second_pv_cover = pv_cover;
								pv_cover = tmp;
								NodeID tmpv = spv;
								spv = pv;
								pv = tmpv;
							}
						}
					}
					else if (second_pv_cover == picked_value[v]) {
						if (sum[spv] < sum[v]) {
							second_pv_cover = picked_value[v];
							spv = v;
						}
						if (sum[pv] < sum[spv]) {
							cover_value_type tmp = second_pv_cover;
							second_pv_cover = pv_cover;
							pv_cover = tmp;
							NodeID tmpv = spv;
							spv = pv;
							pv = tmpv;
						}
					}

				}			

				if (sum[pv] < numOfVertices / 8 && switch_small == false) {
					cout << ranking << endl;
					cout << pv_cover << endl;
					cout << sum[pv] << endl;
					cout << picked_count << endl;
					cout << "Switch to max-heap representations." << endl;
					switch_small = true;
				}
			}
			else {// Later, utilizing the max-heap directly
				double u2_acc_time = GetCurrentTimeSec();
				//while (true) {
				value_heap.extract_min(pv, pv_cover);
				pv_cover = -pv_cover;				
				u2_time += (GetCurrentTimeSec() - u2_acc_time);

				second_pv_cover = -value_heap.top();
				spv = value_heap.top_value();

			}

			cout << dependencies[pv] << "\t" << dependencies[spv] << "\t" << pv_cover << "\t" << second_pv_cover << endl;

			selecting_time += (GetCurrentTimeSec() - selecting_acc_time);

			//// Memcon
			if (switch_small == true && reverse_index.size() > 0 && sum[pv] <= reverse_index[pv].size()) {
				cout << "I think we can stop here " << GetCurrentTimeSec() - start_time << " ranking:" << ranking << endl;
				break;
			}

			// memcon
			/*if (switch_small&& reverse_index.size() > 0)
				cout << ranking << "," << alive_memory << "," << bounded_resources << "," << picked_count << "," << sum[pv] << "," << reverse_index[pv].size() << endl;*/

			/*if (pv_cover == 0 && picked_count == numOfVertices) {
				cout << "time to stop" << endl;
				break;
			}*/

			has_picked[pv] = true;
			double labeling_acc_time = GetCurrentTimeSec();
			labeled_arcs = labeling_source_bfs(pv, que, vis, dst_r, usd, ranking, tmp_idx, r_tmp_idx, graph);
			labeling_time += (GetCurrentTimeSec() - labeling_acc_time);
			inv[ranking] = pv;
			rank[pv] = ranking;

			int current_trees_num = st_roots.size();

			// 4 steps in updating:
			double updating_acc_time = GetCurrentTimeSec();
			//long long released_memory = 0;

			if (switch_small == false) { // We still iterate all alive trees to calculate the cover.	
				for (int tid = 0; tid < current_trees_num; ++tid) {

					// This tree has already be taken because its root has been selected.
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					// pv has already been covered in tree tid.
					int pvpid = 0;
					if (isLarge)
						pvpid = get_large_parent(pv, ltree_tid);
					else
						pvpid = get_small_parent(pv, stree_tid);
					if (pvpid == -1 && st_roots[tid] != pv) // pv has been covered in tree tid 
						continue;


					// Get the descendants and ascendants of pv in tree tid.
					get_descendants(tid, pv, descendants, descendants_parents, que, vis, graph, isLarge, stree_tid, ltree_tid);

					get_ascendants(tid, pv, ascendants, que, vis, graph, isLarge, stree_tid, ltree_tid);

					NodeID subtract_size = descendants.size();
					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
												//cover[tid][av] -= subtract_size;
						cnt16[av][tid%CNT] -= subtract_size;
						if (affected_list[av] == false) {
							affected_list[av] = true;
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
						}
					}


					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					//	vector<int> cover_tid(numOfVertices, 0);
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];

						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//sum[dv] -= cover_tid[dv];
						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];


						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else {
							alive_memory -= stree_tid.bucket_count() * sizeof(NodeID) * 2;
							remove_small_node(dv, stree_tid);
							stree_tid.resize(0);
							alive_memory += stree_tid.bucket_count() * sizeof(NodeID) * 2;
						}

						//	remove_node(tid, dv);
						//affected[dv] = true;
						if (affected_list[dv] == false) {
							affected_count++;
							affected_list[dv] = true;
							if (switch_set)
								affected_v.insert(dv);
						}
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;
					}


					alive_memory -= descendants.size() * sizeof(NodeID); // The size released by reverse index.
				//	released_memory += descendants.size() * sizeof(NodeID); // Released memory.

					// Drop the size of subtree.
					st_sizes[tid] -= subtract_size;
					//if (trees_pointers[tid] != -1) alive_memory -= subtract_size * sizeof(NodeID) * 2;

					//	if (st_sizes[tid] < 0) {
					//		cout << "a5" << st_sizes[tid] << endl;
					//	}
					if (st_sizes[tid] > 0) {
						// 3. Convert shrink trees to small trees representataions.				
						if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
							shrink_tree(tid, st_roots[tid], ltrees[tid]);
							alive_memory -= numOfVertices * sizeof(NodeID);
							alive_memory += strees[trees_pointers[tid]].bucket_count() * sizeof(NodeID) * 2;

						//	released_memory += (numOfVertices * sizeof(NodeID)) - strees[trees_pointers[tid]].bucket_count() * sizeof(NodeID) * 2; // Released Memory.
							
						/*	if (strees[trees_pointers[tid]].bucket_count() * sizeof(NodeID) * 2 > numOfVertices * sizeof(NodeID)) {
								cout << "Are you kidding me?" << endl;
								cout << strees[trees_pointers[tid]].size() << " : " << strees[trees_pointers[tid]].bucket_count() * 2 << endl;
							}*/
						}
					}
					
					// Remove trees
					if (st_roots[tid] == pv)
						if (isLarge) {
							alive_memory -= numOfVertices * sizeof(NodeID);
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
							//ltrees.shrink_to_fit();
						
							//released_memory += numOfVertices * sizeof(NodeID); // Released memory.
						
						}
						else {
							alive_memory -= stree_tid.bucket_count() * sizeof(NodeID) * 2;
							stree_tid.clear();
							stree_tid.resize(0);
							//strees.shrink_to_fit();
							
							//released_memory += stree_tid.bucket_count() * sizeof(NodeID) * 2; // Released memory.
						}				
				}
			}
			else { // We use the reverse index to update trees.

				   // Build the reverse_index first (one time).
				if (reverse_index.size() == 0) {
					reverse_index.resize(numOfVertices);
					for (NodeID v = 0; v < numOfVertices; ++v) {
						reverse_index[v].set_empty_key(numOfVertices + 1);
						reverse_index[v].set_deleted_key(numOfVertices + 2);
					}
					for (int tid = 0; tid < current_trees_num; ++tid) {

						if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue; //This whole tree has already been covered. 

						bool isLarge = false;
						int stree_idx = trees_pointers[tid];
						if (stree_idx == -1)
							isLarge = true;

						small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
						large_tree& ltree_tid = ltrees[tid];

						for (int v = 0; v < numOfVertices; ++v) {
							// Get the descendants and ascendants of pv in tree tid.

							//	if (cover[tid].find(pv) == cover[tid].end()) continue; // This pv has been covered in this tree before.
							//	if (cover[tid][pv] <= 0) continue; // modified-5-14
							int vpid = 0;
							if (isLarge)
								vpid = get_large_parent(v, ltree_tid);
							else
								vpid = get_small_parent(v, stree_tid);

							// To testify whether v is still uncovered in tree tid
							// v has parent vpid != -1, but if v is the root of tree tid, it still can be uncovered.
							if ((vpid != -1 || st_roots[tid] == v))
								reverse_index[v].insert(tid);
							//if (cover[tid].find(v) != cover[tid].end())
							//	reverse_index[v].insert(tid);
						}
					}
				}

				double u4_acc_time = GetCurrentTimeSec();

				//	double u4_acc_time = GetCurrentTimeSec();
				// Access the alive trees directly by the reverse index.
				//				cout << reverse_index[pv].size() << " " << sum[pv] << endl;


				//				for (unordered_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) { // GOOGLE

				for (google::dense_hash_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) {
					
					int tid = *it;
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;
					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					// Get the descendants and ascendants of pv in tree tid.

					//	if (cover[tid].find(pv) == cover[tid].end()) continue; // This pv has been covered in this tree before.
					//	if (cover[tid][pv] <= 0) continue; // modified-5-14

					int pvpid = 0;
					if (isLarge)
						pvpid = get_large_parent(pv, ltree_tid);
					else {
						//	double get_small_parent_incre_time = GetCurrentTimeSec();
						pvpid = get_small_parent(pv, stree_tid);
						//	get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}

					if (pvpid == -1 && st_roots[tid] != pv)
						continue;

					// Array Representation
					//// 1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;
					//get_descendants(tid, pv, descendants, descendants_parents, que, vis, graph.adj, isLarge, stree_tid, ltree_tid);

					//get_ascendants(tid, pv, ascendants, que, vis, graph.adj, graph.r_adj, isLarge, stree_tid, ltree_tid);

					get_descendants(tid, pv, descendants, descendants_parents, que, vis, graph, isLarge, stree_tid, ltree_tid);

					get_ascendants(tid, pv, ascendants, que, vis, graph, isLarge, stree_tid, ltree_tid);

					NodeID subtract_size = descendants.size();

					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
												//						cover[tid][av] -= subtract_size;

						cnt16[av][tid%CNT] -= subtract_size;
						//sum[av] -= subtract_size;

						//	cover_record[tid][av] -= subtract_size;
						//	if (cover_record[tid][av] < 0) {
						//	cout << "4: " << cover_record[tid][av] << " " << subtract_size << " " << tid << " " << endl;
						//	}

						//affected[av] = true;
						/*if (cnt16[av][tid%CNT] < 0)
						cout << " a2 " << cnt16[av][tid%CNT]  << " " << av << " " << st_roots[tid] << " " << tid << endl;*/
						if (affected_list[av] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
							affected_list[av] = true;
						}
					}

					//double u2_acc_time = GetCurrentTimeSec();
					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					//vector<int> cover_tid(numOfVertices, 0);
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//	sum[dv] -= cover_tid[dv];
						//cover[tid].erase(dv);
						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];

						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else {
							alive_memory -= stree_tid.bucket_count() * sizeof(NodeID) * 2;
							remove_small_node(dv, stree_tid);
							stree_tid.resize(0);
							alive_memory += stree_tid.bucket_count() * sizeof(NodeID) * 2;
						}
						//	remove_node(tid, dv);
						//affected[dv] = true;

						//	cover_record[tid][dv] -= cover_tid[dv];
						//if (cover_record[tid][dv] < 0)
						//	cout << "3: " << cover_record[tid][dv] << " " << cover_tid[dv] << " " << di << " " << descendants.size() << " " << picked_as_root[st_roots[tid]] << endl;

						/*if (cnt16[dv][tid%CNT]  < 0)
						cout << " a3 " << cnt16[dv][tid%CNT] << " " << dv << " " << st_roots[tid] << " " << tid << endl;
						*/
						if (dv != pv)
							reverse_index[dv].erase(tid);
						//in_tree_count[dv]--;


						if (affected_list[dv] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(dv);
							affected_list[dv] = true;
						}
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;
					}



					st_sizes[tid] -= subtract_size;
					//if (trees_pointers[tid] != -1) alive_memory -= subtract_size * sizeof(NodeID) * 2;

					//alive_memory -= subtract_size * sizeof(NodeID); // The size released by reverse index.

					//	if (st_sizes[tid]  < 0)
					//		cout << " a4 " << st_sizes[tid] << endl;
					// Remove trees
					if (st_roots[tid] == pv)
						if (isLarge) {
							alive_memory -= numOfVertices * sizeof(NodeID);
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
							//ltrees.shrink_to_fit();

							//cout << "release large" << endl;
							//released_memory += numOfVertices * sizeof(NodeID); // Released Memory.
						}
						else {
							alive_memory -= stree_tid.bucket_count() * sizeof(NodeID) * 2;
							stree_tid.clear();
							stree_tid.resize(0);
							//strees.shrink_to_fit();

							//cout << "release small" << endl;
							//released_memory += strees[trees_pointers[tid]].bucket_count() * sizeof(NodeID) * 2; // Released Memory.
						}

					// 3. Convert shrink trees to small trees representataions.	
					if (st_sizes[tid] > 0) {
						if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
							shrink_tree(tid, st_roots[tid], ltrees[tid]);
							alive_memory -= numOfVertices * sizeof(NodeID);
							alive_memory += strees[trees_pointers[tid]].bucket_count() * sizeof(NodeID) * 2;


							//released_memory += ((numOfVertices * sizeof(NodeID)) - strees[trees_pointers[tid]].bucket_count() * sizeof(NodeID) * 2); // Released Memory.
							
							/*if (strees[trees_pointers[tid]].bucket_count() * sizeof(NodeID) * 2 > numOfVertices * sizeof(NodeID)) {
								cout << "Are you kidding me?" << endl;
								cout << strees[trees_pointers[tid]].size() << " : " << strees[trees_pointers[tid]].bucket_count() * 2 << endl;
							}*/
						}
					}					

				}
			
				alive_memory -= reverse_index[pv].bucket_count()  * sizeof(NodeID);

				reverse_index[pv].clear();
				reverse_index[pv].resize(0);
				u4_time += (GetCurrentTimeSec() - u4_acc_time);
			}

			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			// 4. Sample new trees based on.
			// If we sample a new tree, we need to dynamically increment the following data structures:
			// 1) ltrees
			// 2) trees_pointers (-1 or pointer to strees)
			// 3) cover
			// 4) st_roots
			// 5) st_sizes 
			// 6) strees (maybe)
			// two constraints: cl > ct and 10kn total vertices
			double adding_acc_time = GetCurrentTimeSec();
			new_tree_arcs = 0;

			//cout << "ranking:" << ranking << endl;
			//cout << "pick cover:" << pv_cover << endl;
			//cout << "pick count:" << picked_count << endl;
			//cout << "labeled_arcs:" << labeled_arcs << endl;
			//cout << "new_tree_arcs:" << new_tree_arcs << endl;
			//cout << "alive_resources:" << alive_resources << endl;
			//cout << "bounded_resources:" << bounded_resources << endl;
			//cout << "#tree:" << cover.size() << endl;
			int how_many = 0;
			//cout << "beta*labeld_arcs: " << beta * labeled_arcs << " new tree _arcs: " << new_tree_arcs << " alive " << alive_resources << " bounded " << bounded_resources << endl;
			
			//cout << alive_memory << ":" << picked_count << endl;
		//	long long added_memory = 0;
			//cout << alive_memory << " vs " << bounded_resources << endl;


			while (beta * labeled_arcs > new_tree_arcs&& alive_memory < bounded_resources) {
					
				if (picked_count == numOfVertices) break;
				//cout << ranking << ":" << picked_count << ":" << alive_memory << ":" << bounded_resources << endl;

				NodeID r = rand() % numOfVertices;
				if (picked_as_root[r] == true) continue;
				if (trivial_root[r] == true) continue;
				int tmp_tree_arcs = 0;

				large_tree new_ltree;
				small_tree new_stree;
				//GOOGLE DENSE MAP
				new_stree.set_empty_key(numOfVertices + 1);
				new_stree.set_deleted_key(numOfVertices + 2);

				//	BUILD_SMALL_TREE_FLAG = false;
				bool isLarge = false;
				if (!BUILD_SMALL_TREE_FLAG) {
					
					double u5_acc_time = GetCurrentTimeSec();
					// Array Representation
					//tmp_tree_arcs = build_large_tree_bfs(r, new_ltree, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph.adj, new_tree_arcs, beta * labeled_arcs);
					tmp_tree_arcs = build_large_tree_bfs(r, new_ltree, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph, 0, numOfVertices);
					u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;

				//	cout << "large tree tree arcs:" << tmp_tree_arcs << endl;

					if (descendants.size() <size_thr)
						BUILD_SMALL_TREE_FLAG = true;

					if (descendants.size() < 2) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

				//	cout << "large tree d size:" << descendants.size() << endl;

					isLarge = true;
					how_many++;
					used_as_cover++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					st_roots.push_back(r);
					trees_pointers.push_back(-1);
					st_sizes.push_back(descendants.size());

					//new_tree_arcs += descendants.size();
				//	alive_resources += numOfVertices;

					alive_memory += numOfVertices * sizeof(NodeID);
					//added_memory += numOfVertices * sizeof(NodeID); // Added Memory.
					//cout << new_ltree.capacity() << endl;

					// Turn large tree to small tree.
					// 1. Build new small tree with ltree_i and add it to strees;
					// 2. Set tree_pointers[i] to the position of the strees;
					// 3. Clear the large tree ltree_i by setting its size to 0.
					int i = st_roots.size() - 1;
					if (st_sizes[i] < size_thr) {
						shrink_tree(i, st_roots[i], new_ltree);
						isLarge = false;
						alive_memory -= numOfVertices * sizeof(NodeID);
						alive_memory += strees[trees_pointers[i]].bucket_count() * sizeof(NodeID) * 2;
						
						//added_memory += (strees[trees_pointers[i]].bucket_count() * sizeof(NodeID) * 2) - numOfVertices * sizeof(NodeID); // Added Memory.				
					}
					ltrees.push_back(new_ltree);
				}
				else {//Build Small Tree Directly
					double u5_acc_time = GetCurrentTimeSec();
					// Array Representation
					//	tmp_tree_arcs = build_small_tree_bfs(r, new_stree, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph.adj, new_tree_arcs, beta * labeled_arcs);
					//tmp_tree_arcs = build_small_tree_bfs(r, new_stree, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph, alive_memory / (2*sizeof(NodeID)), bounded_resources / (2*sizeof(NodeID)));
					
					// Gradually Sampling
					tmp_tree_arcs = build_small_tree_bfs(r, new_stree, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph, new_tree_arcs, beta * labeled_arcs);
					u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;

					if (descendants.size() >= size_thr)
						BUILD_SMALL_TREE_FLAG = false;

				//	cout << "small tree tree arcs:" << tmp_tree_arcs << endl;

					if (descendants.size() < 2) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

				//	cout << "small tree d size:" << descendants.size() << endl;

					isLarge = false;

					how_many++;
					used_as_cover++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					st_roots.push_back(r);
					trees_pointers.push_back(strees.size()); //we assume it is inserted.
					st_sizes.push_back(descendants.size());

					//cout << new_stree.size() << endl;
				//	alive_memory += new_stree.max_size() * sizeof(NodeID) * 2;
					//strees.push_back(new_stree);
					//// new_ltree is empty.
					//ltrees.push_back(new_ltree);
					alive_memory += new_stree.bucket_count() * sizeof(NodeID) * 2;
					//added_memory += new_stree.bucket_count() * sizeof(NodeID) * 2; // Added Memory.

					int i = st_roots.size() - 1;
					if (st_sizes[i] >= size_thr) {
						alive_memory -= new_stree.bucket_count() * sizeof(NodeID) * 2;
						alive_memory += numOfVertices * sizeof(NodeID);

						//added_memory = (numOfVertices * sizeof(NodeID)) - new_stree.bucket_count() * sizeof(NodeID) * 2; // Added Memory.

						enlarge_tree(i, st_roots[i], new_stree, new_ltree);
						new_stree.clear();
						new_stree.resize(0);
						isLarge = true;
					}
					else {
						strees.push_back(new_stree);
						ltrees.push_back(new_ltree);
					}
				}

				int i = st_roots.size() - 1;

				int numOfStrees = strees.size();
				small_tree &new_stree_ref = numOfStrees != 0 ? strees[strees.size() - 1] : empty_small_tree;

				// Calculate the subtree size of each nodes.
				//				unordered_map<int, int> cover_new;
				//	vector<cover_value_type> cover_new(numOfVertices, 0); // modified-5-14

				cover_value_type cover_seen = 0;

				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];

					++cover_tid[q]; // Increment because it hitsthe r-q path. 
					++cnt16[q][i%CNT];

					++dependencies[q];

					//++sum[q];
					/*	if (sum[q] > best_sum) {
					best_sum = sum[q];
					best_sum_vid = q;
					}*/
					//affected[q] = true;
					//	cover_record[i][q]++;
					if (affected_list[q] == false) {
						affected_count++;
						if (switch_set)
							affected_v.insert(q);
						affected_list[q] = true;
					}

					if (switch_small == true) {
						reverse_index[q].insert(i);
					}

					int parent_node = -1;

					if (isLarge)
						parent_node = get_large_parent(q, new_ltree);
					else
						parent_node = get_small_parent(q, new_stree_ref);

					if (parent_node > -1) {
						cover_tid[parent_node] += cover_tid[q];
						cnt16[parent_node][i%CNT] += cover_tid[q];

						dependencies[parent_node] += cover_tid[q];
						//	sum[parent_node] += cover_new[q];
						//	cover_record[i][parent_node] += cover_record[i][q];
					}

					if (parent_node == r)
						cover_seen += cover_tid[q];
				}


				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];
					cover_tid[q] = 0;
				}
				
			
				alive_memory += descendants.size() * sizeof(NodeID); // Size occupied by reverse index.

				//added_memory += descendants.size() * sizeof(NodeID); // Size occupied by reverse index.
			}

			//cout << "Round " << ranking << " add - released " << added_memory - released_memory << " released:" << released_memory << endl;

			//// memcon
			//if (switch_small)
			//	cout << ranking << "," << alive_memory << "," << bounded_resources << "," << picked_count << "," << sum[pv] << "," << reverse_index[pv].size() << endl;			
			adding_time += (GetCurrentTimeSec() - adding_acc_time);
			//cout << ranking << ":  beta:" << beta << " beta * labeled_arcs " << beta * labeled_arcs << " new_tree_arcs " << new_tree_arcs << " how many " << how_many << endl;

			//double u3_acc_time = GetCurrentTimeSec();

			updating_acc_time = GetCurrentTimeSec();

			// Build the cover value heap for one time.
			if (switch_small == true && value_heap.empty())
				for (int nid = 0; nid < numOfVertices; ++nid)
					if (has_picked[nid] == false)
						value_heap.update(nid, -picked_value[nid]);

			//	cnt16_2max(picked_value, sum, has_picked, affected, cnt16, switch_small, value_heap);		
			//double u3_acc_time = GetCurrentTimeSec();

			//			cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap, best_sum_vid, lazy_updates);
			cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);
			//u3_time += (GetCurrentTimeSec() - u3_acc_time);

			if (affected_count < numOfVertices / 8)
				switch_set = true;
			else
				switch_set = false;
			affected_count = 0;

			//iteration_time[ranking] = (GetCurrentTimeSec() - updating_acc_time);
			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			time_this_iteration = GetCurrentTimeSec() - time_this_iteration;
			iteration_time[ranking] = time_this_iteration;
		}


		cout << "picked as tree covers:" << used_as_cover << endl;
		/*for (int v = 0; v < V; ++v)
		if (has_picked[v] == false) {
		cout << "wrong 2" << endl;
		cout << picked_value[v] << endl;
		}*/


		if (stop_flag == true) return;

		double labeling_acc_time = GetCurrentTimeSec();
		// Order the remaining vertices with degree, if any.
		if (ranking != numOfVertices) {
			vector<NodeID> degree(numOfVertices);
			vector<NodeID> r_degree(numOfVertices);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				/*vector<NodeID>& adj_v = graph.adj[v];
				for (NodeID i = 0; i < adj_v.size(); ++i) {
					if (has_picked[adj_v[i]] != true)
						degree[v]++;
				}*/

				for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid) {
					NodeID w = graph.edges[eid];
					if (has_picked[w] != true)
						degree[v]++;
				}

				if (DIRECTED_FLAG == true) {
					/*vector<NodeID>& r_adj_v = graph.r_adj[v];
					for (NodeID i = 0; i < r_adj_v.size(); ++i) {
						if (has_picked[r_adj_v[i]] != true)
							r_degree[v]++;
					}*/
					for (EdgeID eid = graph.r_vertices[v]; eid < graph.r_vertices[v + 1]; ++eid) {
						NodeID w = graph.r_edges[eid];
						if (has_picked[w] != true)
							r_degree[v]++;
					}
				}
			}

			vector<pair<NodeID, NodeID> > degree_v;
			degree_v.reserve(numOfVertices - ranking);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				if (DIRECTED_FLAG == true)
					degree_v.push_back(make_pair(degree[v] * r_degree[v], v));
				else
					degree_v.push_back(make_pair(degree[v], v));
			}
			sort(degree_v.rbegin(), degree_v.rend());

			int count = 0;
			for (; ranking < numOfVertices; ++ranking) {
				NodeID tpv = degree_v[count++].second;
				inv[ranking] = tpv;
				rank[tpv] = ranking;
				labeling_source_bfs(tpv, que, vis, dst_r, usd, ranking, tmp_idx, r_tmp_idx, graph);
			}
		}
		labeling_time += (GetCurrentTimeSec() - labeling_acc_time);

		// Finish the constructions.
		for (NodeID v = 0; v < numOfVertices; ++v) {
			if (DIRECTED_FLAG == false) {
				NodeID k = tmp_idx[v].first.size();
				labels.index_[v].spt_v.resize(k);
				labels.index_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();

				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();

			}
			else {
				NodeID k = tmp_idx[v].first.size();
				dlabels.index_[v].spt_v.resize(k);
				dlabels.index_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();

				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();

				k = r_tmp_idx[v].first.size();
				dlabels.bindex_[v].spt_v.resize(k);
				dlabels.bindex_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_v[i] = r_tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_d[i] = r_tmp_idx[v].second[i];
				r_tmp_idx[v].first.clear();
				r_tmp_idx[v].second.clear();

				r_tmp_idx[v].first.shrink_to_fit();
				r_tmp_idx[v].second.shrink_to_fit();
			}
		}

		num_of_trees = st_roots.size();

		double small_tree_num = 0;

		for (int i = 0; i < num_of_trees; ++i) {
			if (trees_pointers[i] != -1) {
				small_tree_num++;
				//small_tree_num += strees[trees_pointers[i]].size();
			}
		}

		total_resources = small_tree_num;
		/*for (int i = 0; i < cnt16.size(); ++i){
		for (int j = 0; j < cnt16[i].size(); ++j)
		if (cnt16[i][j] != 0)
		cout << "aaaa " << cnt16[i][j] << endl;
		}*/

	}



		

	Betweenness_Ordering(NodeID k, double beta, WGraph& wgraph, NodeID stopv, long long bound_times = 10, long long trans_times = 8) {

		double start_time = GetCurrentTimeSec();

		if (k < CNT) {
			cout << " Should pick at least 16 seed!" << endl;
			return;
		}
		iteration_time.resize(numOfVertices);

		updating_time = 0;
		labeling_time = 0;
		adding_time = 0;
		selecting_time = 0;
		init_time = 0;
		u1_time = 0;
		u2_time = 0;
		u3_time = 0;
		num_of_trees = 0;
		total_resources = 0;


		inv.resize(numOfVertices, -1);
		rank.resize(numOfVertices, -1);

		if (DIRECTED_FLAG == false) {
			labels.index_.resize(numOfVertices);
		}
		else { //if (DIRECTED_FLAG == true) {
			dlabels.index_.resize(numOfVertices);
			dlabels.bindex_.resize(numOfVertices);
		}

		// Resource Bounding.
		NodeID new_tree_arcs = 0;
		NodeID labeled_arcs = 0;
		alive_resources = 0;
		bounded_resources = bound_times * (long long)k * (long long)numOfVertices;

		vector<NodeID> st_roots; // The root vertices of the sampled trees.
		vector<NodeID> st_sizes; // The size of the sampled trees.

		NodeID size_thr = numOfVertices < 8 ? numOfVertices : numOfVertices / trans_times; // When a tree shrink to a size smaller than the size_thr, it should be converted into a small_tree representation.

		init_time = GetCurrentTimeSec();
		// Building the initial k trees.
	//	srand((unsigned)time(NULL));
		srand(0x728f9d0);


		ltrees.reserve(numOfVertices);
		strees.reserve(numOfVertices);
		trees_pointers.reserve(numOfVertices);
		st_roots.reserve(numOfVertices);
		st_sizes.reserve(numOfVertices);
		st_roots.resize(k);
		st_sizes.resize(k);

		ltrees.resize(k);
		//strees.clear();
		trees_pointers.resize(k, -1); // Point the sampling trees to its tree representations. large tree [0, V) small tree [V, 2V).
		//vector<unordered_map<int, int> > cover(k, unordered_map<int, int>());

		vector<NodeID> ascendants;
		ascendants.reserve(numOfVertices);
		vector<NodeID> descendants;
		descendants.reserve(numOfVertices);
		vector<NodeID> que(numOfVertices);
		vector<bool> vis(numOfVertices);
		vector<bool> has_picked(numOfVertices);
		vector<vector<cover_value_type> > cnt16(numOfVertices, vector<cover_value_type>(CNT));
		//		vector<bool> affected(V, true);
		vector<bool> picked_as_root(numOfVertices, false);
		vector<bool> trivial_root(numOfVertices, false);
		int picked_count = 0;
		int used_as_cover = 0;
		vector<cover_value_type> sum(numOfVertices);
		bool stop_flag = false;


		// tmp structure used in the dijkstra searches. 

		// Binary Heap
		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);

		queue<NodeID> visited_que;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
		vector<bool> usd(numOfVertices, false);

		//GOOGLE
//		unordered_set<NodeID> affected_v;
		google::dense_hash_set<NodeID> affected_v;
		affected_v.set_empty_key(numOfVertices + 1);
		affected_v.set_deleted_key(numOfVertices + 2);

		vector<cover_value_type> cover_tid(numOfVertices, 0);
		//vector<cover_value_type> cover_new(numOfVertices, 0); // modified-5-14


		vector<bool> affected_list(numOfVertices, false);
		int affected_count = 0;
		bool switch_set = false;

		//vector<unordered_set<int> > reverse_index(0);
		//GOOGLE DENSE SET
		vector<google::dense_hash_set<NodeID> > reverse_index;
		//vector<NodeID> in_tree_count(numOfVertices);

		// Preparing basic structure for pl algorithm.
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT))); // Backward labels.

		vector<NodeID> descendants_parents(numOfVertices, -1);

		// Building the initial k sampled trees. 
		for (int i = 0; i < k; ++i) {
			NodeID r = rand() % numOfVertices;
			if (picked_as_root[r] == true) { --i; continue; }

			st_roots[i] = r;
			large_tree &ltree_i = ltrees[i];
			picked_as_root[r] = true;
			picked_count++;
			used_as_cover++;

			build_large_tree_dij(r, ltree_i, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, 0, numOfVertices);
			st_sizes[i] = descendants.size();
			alive_resources += (long long)numOfVertices;

			// Turn large tree to small tree.
			// 1. Build new small tree with ltree_i and add it to strees;
			// 2. Set tree_pointers[i] to the position of the strees;
			// 3. Clear the large tree ltree_i by setting its size to 0.
			if (st_sizes[i] < size_thr)
				shrink_tree(i, st_roots[i], ltrees[i]);

			bool isLarge = false;
			int stree_idx = trees_pointers[i];
			if (stree_idx == -1)
				isLarge = true;

			small_tree& stree_i = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;

			// Calculate the subtree size of each nodes.
		//	unordered_map<int, int>& cover_i = cover[i];
			//vector<int> cover_i(numOfVertices, 0);

			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];

				++cover_tid[q]; // Increment because it hits the r-q path.
				++cnt16[q][i%CNT];
				//++sum[q];

				int parent_node;

				if (isLarge)
					parent_node = get_large_parent(q, ltree_i);
				else
					parent_node = get_small_parent(q, stree_i);

				if (parent_node > -1) {
					cover_tid[parent_node] += cover_tid[q];
					cnt16[parent_node][i%CNT] += cover_tid[q];
					//	sum[parent_node] += cover_i[q];
				}
			}
			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];
				cover_tid[q] = 0;
			}
		}



		//cnt16_aggreate(cnt16, has_picked, cover);

		vector<cover_value_type> picked_value;
		benchmark::heap<2, cover_value_type, NodeID> value_heap(numOfVertices);
		//bool switch_small = false;
		switch_small = false;

		//cnt16_2max(picked_value, sum, has_picked,affected, cnt16, switch_small, value_heap);
//		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, cnt16, switch_small, value_heap);
		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);

		init_time = GetCurrentTimeSec() - init_time;

		cout << "Finish initial sampling " << init_time << "secs." << endl;

		//for (NodeID v = 0; v < V; ++v)	value_heap.update(v, -picked_value[v]); // Min-heap.

		double max_value = -1;

		// Perform V iterations to rank all vertices.
		NodeID ranking = 0;
		for (; ranking < numOfVertices; ++ranking) {
			if (ranking % (numOfVertices / 10) == 0)
				cout << "picking the " << ranking << "-th orders." << endl;

			double time_this_iteration = GetCurrentTimeSec();

			NodeID pv = 0;
			long long pv_cover = -1;

			if (ranking == stopv) {
				stop_flag = true;
				cout << "topv labeling time:" << labeling_time << endl;
				break;
			}

			double selecting_acc_time = GetCurrentTimeSec();
			// In the beginning, pick the best cover vertex by a linear scan.
			if (switch_small == false) {
				double u1_acc_time = GetCurrentTimeSec();

				for (NodeID v = 0; v < numOfVertices; ++v) {
					if (has_picked[v] == true) continue;
					if (pv_cover == picked_value[v]) {
						if (sum[pv] < sum[v]) {
							pv_cover = picked_value[v];
							pv = v;
						}
					}
					else if (pv_cover < picked_value[v]) {
						pv_cover = picked_value[v];
						pv = v;
					}
				}

				if (sum[pv] < numOfVertices / 8 && switch_small == false) {
					//cout << "u3 time:" << u3_time << endl;
					cout << ranking << endl;
					cout << pv_cover << endl;
					cout << sum[pv] << endl;
					cout << picked_count << endl;
					cout << "Switch to max-heap representations." << endl;
					switch_small = true;
				}
				u1_time += (GetCurrentTimeSec() - u1_acc_time);

			}
			else {// Later, utilizing the max-heap directly
				double u2_acc_time = GetCurrentTimeSec();

				value_heap.extract_min(pv, pv_cover);
				pv_cover = -pv_cover;
				u2_time += (GetCurrentTimeSec() - u2_acc_time);
				/*
				int max1 = 0;
				int max2 = 0;
				int max1v = -1;
				int max2v = -1;
				for(int  i = 0; i < CNT; ++i){
					if(max2 <= cnt16[pv][i]){
						max2 = cnt16[pv][i];
						max2v = i;
					}
					if(max1 <= max2){
						int tmp = max1;
						max1 = max2;
						max2 = tmp;
						int tmpv = max1v;
						max1v = max2v;
						max2v = max1v;
					}
				}
				long long sumstd = 0;
				long long summean = 0;
				for(int  i = 0; i < CNT; ++i){
					if(i != max1v && i != max2v){
						sumstd += (long long)cnt16[pv][i] * (long long)cnt16[pv][i];
						summean += (long long)cnt16[pv][i];
					}
				}
				sumstd = sumstd/(long long)(CNT-2);
				summean = summean/(long long)(CNT-2);
				cout << sqrt(sumstd - summean * summean) << endl;*/
			}

			selecting_time += (GetCurrentTimeSec() - selecting_acc_time);

			//// FIX
			//if (switch_small == true && reverse_index.size() > 0 && sum[pv] <= reverse_index[pv].size() * 10) {//sum[pv] <= in_tree_count[pv] * 10){//
			//	cout << "I think we can stop here " << GetCurrentTimeSec() - start_time << " ranking:" << ranking << endl;
			//	break;
			//}

			//if (has_picked[pv] == true)
			//	cout << "wrong1" << endl;

			has_picked[pv] = true;
			double labeling_acc_time = GetCurrentTimeSec();
			labeled_arcs = labeling_source_dij(pv, pqueue, visited_que, distances, vis, dst_r, usd, ranking, tmp_idx, r_tmp_idx, wgraph);
			labeling_time += (GetCurrentTimeSec() - labeling_acc_time);
			inv[ranking] = pv;
			rank[pv] = ranking;

			int current_trees_num = st_roots.size();

			// 4 steps in updating:
			double updating_acc_time = GetCurrentTimeSec();
			if (switch_small == false) { // We still iterate all alive trees to calculate the cover.	
			//	double u3_acc_time = GetCurrentTimeSec();
				for (int tid = 0; tid < current_trees_num; ++tid) {

					// Get the descendants and ascendants of pv in tree tid.
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					// pv has already been covered in tree tid.
					int pvpid = 0;
					if (isLarge)
						pvpid = get_large_parent(pv, ltree_tid);
					else {
						//double get_small_parent_incre_time = GetCurrentTimeSec();
						pvpid = get_small_parent(pv, stree_tid);
						//	get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}
					if (pvpid == -1 && st_roots[tid] != pv) // pv has been covered in tree tid 
						continue;


					double u3_acc_time = GetCurrentTimeSec();

					// //1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;				 
					//get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph.adj, isLarge, stree_tid, ltree_tid);
					//get_ascendants(tid, pv, ascendants, que, vis, wgraph.adj, wgraph.r_adj, isLarge, stree_tid, ltree_tid);

					//Array Representation
					get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph, isLarge, stree_tid, ltree_tid);
					get_ascendants(tid, pv, ascendants, que, vis, wgraph, isLarge, stree_tid, ltree_tid);

					u3_time += (GetCurrentTimeSec() - u3_acc_time);


					double u4_acc_time = GetCurrentTimeSec();
					NodeID subtract_size = descendants.size();
					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
					//	cover[tid][av] -= subtract_size;
						cnt16[av][tid%CNT] -= subtract_size;
						//sum[av] -= subtract_size;
						//affected[av] = true;

						if (affected_list[av] == false) {
							affected_list[av] = true;
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
						}
					}
					u4_time += (GetCurrentTimeSec() - u4_acc_time);


					double u5_acc_time = GetCurrentTimeSec();
					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					//vector<int> cover_tid(numOfVertices, 0);
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//sum[dv] -= cover_tid[dv];

						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];
						if (cnt16[dv][tid%CNT] < 0)
							cout << "a6 " << cnt16[dv][tid%CNT] << endl;

						/*if (isLarge)
							remove_large_node(dv, ltree_tid);
						else
							remove_small_node(dv, stree_tid);*/

							//						remove_node(tid, dv);
													//affected[dv] = true;

						if (affected_list[dv] == false) {
							affected_count++;
							affected_list[dv] = true;
							if (switch_set)
								affected_v.insert(dv);
						}
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;

						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else 
							remove_small_node(dv, stree_tid);
					}
					u5_time += (GetCurrentTimeSec() - u5_acc_time);

					// Drop the size of subtree.
					st_sizes[tid] -= subtract_size;
					if (trees_pointers[tid] != -1) alive_resources -= (long long)subtract_size;

					// 3. Convert shrink trees to small trees representataions.	
					if (st_sizes[tid] > 0) {
						if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
							shrink_tree(tid, st_roots[tid], ltrees[tid]);
						}
					}

					if (st_roots[tid] == pv)
						if (isLarge) {
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
							//ltrees.shrink_to_fit();
						}
						else {
							stree_tid.clear();
							stree_tid.resize(0);
							//strees.shrink_to_fit();
						}

				}
				//	u3_time += (GetCurrentTimeSec() - u3_acc_time);
			}
			else { // We use the reverse index to update trees.
			   // Build the reverse_index first (one time).
				if (reverse_index.size() == 0) {
					reverse_index.resize(numOfVertices);
					for (NodeID v = 0; v < numOfVertices; ++v) {
						reverse_index[v].set_empty_key(numOfVertices + 1);
						reverse_index[v].set_deleted_key(numOfVertices + 2);
					}
					for (int tid = 0; tid < current_trees_num; ++tid) {
						if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue; //This whole tree has already been covered. 

						bool isLarge = false;
						int stree_idx = trees_pointers[tid];
						if (stree_idx == -1)
							isLarge = true;

						small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
						large_tree& ltree_tid = ltrees[tid];

						for (int v = 0; v < numOfVertices; ++v) {

							int vpid = 0;
							if (isLarge)
								vpid = get_large_parent(v, ltree_tid);
							else {
								//double get_small_parent_incre_time = GetCurrentTimeSec();
								vpid = get_small_parent(v, stree_tid);
								//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
							}
							if ((vpid != -1 || st_roots[tid] == v)) {
								reverse_index[v].insert(tid);
								//in_tree_count[v]++;
							}
						}
					}
				}

				//cout << reverse_index[pv].size() << "," << sum[pv] << "," <<  (double)sum[pv] / (double)reverse_index[pv].size() << endl;
				// Access the alive trees directly by the reverse index.
				//for (unordered_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) {
				for (google::dense_hash_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) {
					int tid = *it;
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					double u6_acc_time = GetCurrentTimeSec();
					/*
							int parent_node;
							if (isLarge)
								parent_node = get_large_parent(pv, ltree_tid);
							else {
								//double get_small_parent_incre_time = GetCurrentTimeSec();
								parent_node = get_small_parent(pv, stree_tid);
								//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
							}
							if(parent_node == -1 && pv != st_roots[tid])
								continue;
						*/
					double u7_acc_time = GetCurrentTimeSec();

					//	// 1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;
					//	//double u4_acc_time = getcurrenttimesec();
					//	get_descendants(tid, pv, descendants, descendants_parents,que, vis, wgraph.adj, islarge, stree_tid, ltree_tid);
					////	u4_time += (getcurrenttimesec() - u4_acc_time);

					//	get_ascendants(tid, pv, ascendants, que, vis, wgraph.adj, wgraph.r_adj, isLarge, stree_tid, ltree_tid);

						// Array Representation
					get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph, isLarge, stree_tid, ltree_tid);
					get_ascendants(tid, pv, ascendants, que, vis, wgraph, isLarge, stree_tid, ltree_tid);

					NodeID subtract_size = descendants.size();
					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
						//cover[tid][av] -= subtract_size;
						cnt16[av][tid%CNT] -= subtract_size;
						//affected[av] = true;
					//	sum[av] -= subtract_size;

						if (affected_list[av] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
							affected_list[av] = true;
						}
					}
					u7_time += (GetCurrentTimeSec() - u7_acc_time);

					//double u2_acc_time = GetCurrentTimeSec();
					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					double u8_acc_time = GetCurrentTimeSec();
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//	sum[dv] -= cover_tid[dv];

						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];
						//	if (cnt16[dv][tid%CNT] < 0)
						//		cout << "a6 " << cnt16[dv][tid%CNT] << endl;

							//if(dv == pv){
								/*if (isLarge)
									remove_large_node(dv, ltree_tid);
								else
									remove_small_node(dv, stree_tid);*/
									//}
									//affected[dv] = true;

									//double u82_acc_time = GetCurrentTimeSec();
						if (dv != pv)
							reverse_index[dv].erase(tid);
						//if (dv != pv)
						//	in_tree_count[dv]--;

					//	double u81_acc_time = GetCurrentTimeSec();
						if (affected_list[dv] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(dv);
							affected_list[dv] = true;
						}

						//u81_time += (GetCurrentTimeSec() - u81_acc_time);

					//	u82_time += (GetCurrentTimeSec() - u82_acc_time);
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;
						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else 
							remove_small_node(dv, stree_tid);
					}
					u8_time += (GetCurrentTimeSec() - u8_acc_time);

					if (st_roots[tid] == pv)
						if (isLarge) {
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
						}
						else{
							stree_tid.clear();
							stree_tid.resize(0);
						}


					st_sizes[tid] -= subtract_size;
					if (trees_pointers[tid] != -1) alive_resources -= (long long)subtract_size;
					//	if (st_sizes[tid] < 0) {
					//		cout << "a5" << st_sizes[tid] << endl;
					//	}

						// 3. Convert shrink trees to small trees representataions.				
					if (st_sizes[tid] > 0) {
						if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
							shrink_tree(tid, st_roots[tid], ltrees[tid]);
						}
					}

					//u2_time += (GetCurrentTimeSec() - u2_acc_time);
				}
				reverse_index[pv].clear();
				reverse_index[pv].resize(0);

			}
			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			// 4. Sample new trees based on.
			// If we sample a new tree, we need to dynamically increment the following data structures:
			// 1) ltrees
			// 2) trees_pointers (-1 or pointer to strees)
			// 3) cover
			// 4) st_roots
			// 5) st_sizes 
			// 6) strees (maybe)
			// two constraints: cl > ct and 10kn total vertices
			double adding_acc_time = GetCurrentTimeSec();


			/*cout << "ranking:" << ranking << endl;
			cout << "pick cover:" << pv_cover << endl;
			cout << "pick count:" << picked_count << endl;
			cout << "labeled_arcs:" << labeled_arcs << endl;
			cout << "new_tree_arcs:" << new_tree_arcs << endl;
			cout << "alive_resources:" << alive_resources << endl;
			cout << "bounded_resources:" << bounded_resources << endl;
			cout << "#tree:" << cover.size() << endl;
			cout << "u5/adding" << u5_time / adding_time << endl;
*/
			new_tree_arcs = 0;
			int how_many = 0;
			while (beta * labeled_arcs > new_tree_arcs && alive_resources < bounded_resources) {

				if (picked_count == numOfVertices) break;

				NodeID r = rand() % numOfVertices;
				if (picked_as_root[r] == true) continue;
				if (trivial_root[r] == true) continue;

				large_tree new_ltree;
				small_tree new_stree;
				//GOOGLE DENSE MAP
				new_stree.set_empty_key(numOfVertices + 1);
				new_stree.set_deleted_key(numOfVertices + 2);
				bool isLarge = false;
				int tmp_tree_arcs = 0;

				if (!BUILD_SMALL_TREE_FLAG) {

					//double u5_acc_time = GetCurrentTimeSec();
					double u9_acc_time = GetCurrentTimeSec();
					tmp_tree_arcs = build_large_tree_dij(r, new_ltree, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, new_tree_arcs, beta * labeled_arcs);
					u9_time += (GetCurrentTimeSec() - u9_acc_time);
					//u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;

					if (descendants.size() < size_thr)
						BUILD_SMALL_TREE_FLAG = true;

					if (descendants.size() < 100) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					//	double u91_acc_time = GetCurrentTimeSec();

					isLarge = true;
					how_many++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					used_as_cover++;
					st_roots.push_back(r);
					trees_pointers.push_back(-1);
					st_sizes.push_back(descendants.size());

					//new_tree_arcs += descendants.size();
					alive_resources += (long long)numOfVertices;

					//	u91_time += (GetCurrentTimeSec() - u91_acc_time);
						// Turn large tree to small tree.
						// 1. Build new small tree with ltree_i and add it to strees;
						// 2. Set tree_pointers[i] to the position of the strees;
						// 3. Clear the large tree ltree_i by setting its size to 0.

					//	double u92_acc_time = GetCurrentTimeSec();
					int i = st_roots.size() - 1;
					if (st_sizes[i] < size_thr) {
						shrink_tree(i, st_roots[i], new_ltree);
						// FIX
						isLarge = false;
						/*for (int di = descendants.size() - 1; di > -1; --di) {
							NodeID q = descendants[di];
							NodeID v_parent = get_large_parent(q, new_ltree);
							NodeID s_parent = get_small_parent(q, strees[strees.size() - 1]);
							if (v_parent != s_parent) cout << v_parent << " vs " << s_parent << endl;
						}*/
					}
					ltrees.push_back(new_ltree);
					//	u92_time += (GetCurrentTimeSec() - u92_acc_time);
				}
				else {//Build Small Tree Directly
					//double u5_acc_time = GetCurrentTimeSec();
					double u9_acc_time = GetCurrentTimeSec();
					tmp_tree_arcs = build_small_tree_dij(r, new_stree, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, new_tree_arcs, beta * labeled_arcs);
					//u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;
					u9_time += (GetCurrentTimeSec() - u9_acc_time);

					if (descendants.size() >= size_thr)
						BUILD_SMALL_TREE_FLAG = false;

					if (descendants.size() < 100) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					isLarge = false;

					how_many++;
					used_as_cover++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					st_roots.push_back(r);
					trees_pointers.push_back(strees.size()); //we assume it is inserted.
					st_sizes.push_back(descendants.size());

					alive_resources += (long long)new_stree.size();
					//strees.push_back(new_stree);
					//// new_ltree is empty.
					//ltrees.push_back(new_ltree);

					int i = st_roots.size() - 1;
					if (st_sizes[i] >= size_thr) {
						enlarge_tree(i, st_roots[i], new_stree, new_ltree);
						new_stree.clear();
						new_stree.resize(0);
						isLarge = true;
					}
					else {
						strees.push_back(new_stree);
						ltrees.push_back(new_ltree);
					}

				}

				int i = st_roots.size() - 1;

				int numOfStrees = strees.size();
				small_tree &new_stree_ref = numOfStrees != 0 ? strees[strees.size() - 1] : empty_small_tree;

				// Calculate the subtree size of each nodes.

				double u10_acc_time = GetCurrentTimeSec();
				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];

					++cover_tid[q]; // Increment because it hitsthe r-q path.
					++cnt16[q][i%CNT];
					//++sum[q];
					//affected[q] = true;

					if (affected_list[q] == false) {
						affected_count++;
						if (switch_set)
							affected_v.insert(q);
						affected_list[q] = true;
					}

					if (switch_small == true) {
						reverse_index[q].insert(i);
						//in_tree_count[q]++;
					}


					int parent_node;
					if (isLarge)
						parent_node = get_large_parent(q, new_ltree);
					else {
						//double get_small_parent_incre_time = GetCurrentTimeSec();
						parent_node = get_small_parent(q, new_stree_ref);
						//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}
					if (parent_node > -1) {
						cover_tid[parent_node] += cover_tid[q];
						cnt16[parent_node][i%CNT] += cover_tid[q];
						//sum[parent_node] += cover_new[q];
					}
				}
				u10_time += (GetCurrentTimeSec() - u10_acc_time);
				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];
					cover_tid[q] = 0;
				}
				//cout << "beta*labeld_arcs: " << beta * labeled_arcs << " new tree _arcs: " << new_tree_arcs << " alive " << alive_resources << " bounded " << bounded_resources << endl;

			}
			adding_time += (GetCurrentTimeSec() - adding_acc_time);


			//double u3_acc_time = GetCurrentTimeSec();

			updating_acc_time = GetCurrentTimeSec();

			// Build the cover value heap for one time.
			if (switch_small == true && value_heap.empty())
				for (int nid = 0; nid < numOfVertices; ++nid)
					if (has_picked[nid] == false)
						value_heap.update(nid, -picked_value[nid]);

			//	cnt16_2max(picked_value, sum, has_picked, affected, cnt16, switch_small, value_heap);
		//c	cnt16_2max_queue(picked_value, sum, has_picked, affected_v, cnt16, switch_small, value_heap);

			double u11_acc_time = GetCurrentTimeSec();
			cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);
			u11_time += (GetCurrentTimeSec() - u11_acc_time);
			//	iteration_time[ranking] = (GetCurrentTimeSec() - updating_acc_time);
			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			if (affected_count < numOfVertices / 8)
				switch_set = true;
			else
				switch_set = false;
			affected_count = 0;

			time_this_iteration = GetCurrentTimeSec() - time_this_iteration;
			iteration_time[ranking] = time_this_iteration;

		}

		cout << "picked as tree covers:" << used_as_cover << endl;
		cout << "picked count:" << picked_count << endl;
		
		//// FIX
		//for (int v = 0; v < numOfVertices; ++v){
		//	for (int i = 0; i < CNT; ++i) {
		//		if (cnt16[v][i] != 0)
		//			cout << "mi?" << endl;
		//	}
		//}
		if (stop_flag == true) return;

		double labeling_acc_time = GetCurrentTimeSec();
		// Order the remaining vertices with degree, if any.
		if (ranking != numOfVertices) {
			vector<NodeID> degree(numOfVertices);
			vector<NodeID> r_degree(numOfVertices);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				/*vector<NodeID>& adj_v = graph.adj[v];
				for (NodeID i = 0; i < adj_v.size(); ++i) {
				if (has_picked[adj_v[i]] != true)
				degree[v]++;
				}*/

				for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid) {
					NodeID w = wgraph.edges[eid].first;
					if (has_picked[w] != true)
						degree[v]++;
				}

				if (DIRECTED_FLAG == true) {
					/*vector<NodeID>& r_adj_v = graph.r_adj[v];
					for (NodeID i = 0; i < r_adj_v.size(); ++i) {
					if (has_picked[r_adj_v[i]] != true)
					r_degree[v]++;
					}*/
					for (EdgeID eid = wgraph.r_vertices[v]; eid < wgraph.r_vertices[v + 1]; ++eid) {
						NodeID w = wgraph.r_edges[eid].first;
						if (has_picked[w] != true)
							r_degree[v]++;
					}
				}
			}
			
			vector<pair<NodeID, NodeID> > degree_v;
			degree_v.reserve(numOfVertices - ranking);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				if (DIRECTED_FLAG == true)
					degree_v.push_back(make_pair(degree[v] * r_degree[v], v));
				else
					degree_v.push_back(make_pair(degree[v], v));
			}
			sort(degree_v.rbegin(), degree_v.rend());

			int count = 0;
			for (; ranking < numOfVertices; ++ranking) {
				NodeID tpv = degree_v[count++].second;
				inv[ranking] = tpv;
				rank[tpv] = ranking;
				labeling_source_dij(tpv, pqueue, visited_que, distances, vis, dst_r, usd, ranking, tmp_idx, r_tmp_idx, wgraph);
			}
		}
		labeling_time += (GetCurrentTimeSec() - labeling_acc_time);


		// Finish the constructions.
		for (NodeID v = 0; v < numOfVertices; ++v) {
			if (DIRECTED_FLAG == false) {
				NodeID k = tmp_idx[v].first.size();
				labels.index_[v].spt_v.resize(k);
				labels.index_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();
				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();

			}
			else {
				NodeID k = tmp_idx[v].first.size();
				dlabels.index_[v].spt_v.resize(k);
				dlabels.index_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();
				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();

				k = r_tmp_idx[v].first.size();
				dlabels.bindex_[v].spt_v.resize(k);
				dlabels.bindex_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_v[i] = r_tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_d[i] = r_tmp_idx[v].second[i];
				r_tmp_idx[v].first.clear();
				r_tmp_idx[v].second.clear();
				r_tmp_idx[v].first.shrink_to_fit();
				r_tmp_idx[v].second.shrink_to_fit();
			}
		}

		num_of_trees = st_roots.size();

		double small_tree_num = 0;

		for (int i = 0; i < num_of_trees; ++i) {
			if (trees_pointers[i] != -1) {
				small_tree_num++;
				//small_tree_num += strees[trees_pointers[i]].size();
			}
		}
		
		total_resources = small_tree_num;
		
	}
	
	Betweenness_Ordering(NodeID k, double beta, WGraph& wgraph, NodeID stopv, long long bound_times, long long trans_times , bool slevel, bool directed_flags) {
		DIRECTED_FLAG = directed_flags;
		SECOND_LEVEL = slevel;
		
		double start_time = GetCurrentTimeSec();

		if (k < CNT) {
			cout << " Should pick at least 16 seed!" << endl;
			return;
		}
		iteration_time.resize(numOfVertices);

		updating_time = 0;
		labeling_time = 0;
		adding_time = 0;
		selecting_time = 0;
		init_time = 0;
		u1_time = 0;
		u2_time = 0;
		u3_time = 0;
		num_of_trees = 0;
		total_resources = 0;


		inv.resize(numOfVertices, -1);
		rank.resize(numOfVertices, -1);

		/*
		if (DIRECTED_FLAG == false) {
			labels.index_.resize(numOfVertices);
		}
		else { //if (DIRECTED_FLAG == true) {
			dlabels.index_.resize(numOfVertices);
			dlabels.bindex_.resize(numOfVertices);
		}*/

		// Resource Bounding.
		NodeID new_tree_arcs = 0;
		NodeID labeled_arcs = 0;
		alive_resources = 0;
		bounded_resources = bound_times * (long long)k * (long long)numOfVertices;

		vector<NodeID> st_roots; // The root vertices of the sampled trees.
		vector<NodeID> st_sizes; // The size of the sampled trees.

		NodeID size_thr = numOfVertices < 8 ? numOfVertices : numOfVertices / trans_times; // When a tree shrink to a size smaller than the size_thr, it should be converted into a small_tree representation.

		init_time = GetCurrentTimeSec();
		// Building the initial k trees.
	//	srand((unsigned)time(NULL));
		srand(0x728f9d0);


		ltrees.reserve(numOfVertices);
		strees.reserve(numOfVertices);
		trees_pointers.reserve(numOfVertices);
		st_roots.reserve(numOfVertices);
		st_sizes.reserve(numOfVertices);
		st_roots.resize(k);
		st_sizes.resize(k);

		ltrees.resize(k);
		//strees.clear();
		trees_pointers.resize(k, -1); // Point the sampling trees to its tree representations. large tree [0, V) small tree [V, 2V).
		//vector<unordered_map<int, int> > cover(k, unordered_map<int, int>());

		vector<NodeID> ascendants;
		ascendants.reserve(numOfVertices);
		vector<NodeID> descendants;
		descendants.reserve(numOfVertices);
		vector<NodeID> que(numOfVertices);
		vector<bool> vis(numOfVertices);
		vector<bool> has_picked(numOfVertices);
		vector<vector<cover_value_type> > cnt16(numOfVertices, vector<cover_value_type>(CNT));
		//		vector<bool> affected(V, true);
		vector<bool> picked_as_root(numOfVertices, false);
		vector<bool> trivial_root(numOfVertices, false);
		int picked_count = 0;
		int used_as_cover = 0;
		vector<cover_value_type> sum(numOfVertices);
		bool stop_flag = false;

		

		// tmp structure used in the dijkstra searches. 

		// Binary Heap
		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);

		queue<NodeID> visited_que;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
		vector<bool> usd(numOfVertices, false);

		//GOOGLE
//		unordered_set<NodeID> affected_v;
		google::dense_hash_set<NodeID> affected_v;
		affected_v.set_empty_key(numOfVertices + 1);
		affected_v.set_deleted_key(numOfVertices + 2);

		vector<cover_value_type> cover_tid(numOfVertices, 0);
		//vector<cover_value_type> cover_new(numOfVertices, 0); // modified-5-14


		vector<bool> affected_list(numOfVertices, false);
		int affected_count = 0;
		bool switch_set = false;

		//vector<unordered_set<int> > reverse_index(0);
		//GOOGLE DENSE SET
		vector<google::dense_hash_set<NodeID> > reverse_index;
		//vector<NodeID> in_tree_count(numOfVertices);

		// Preparing basic structure for pl algorithm.
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<EdgeWeight> r_dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT))); // Backward labels.

		
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
		tmp_idx_token_parents(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
			vector<EdgeWeight>(1, numOfVertices)));
			
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx_token_parents(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, numOfVertices)));

				
		vector<NodeID> descendants_parents(numOfVertices, -1);

		// Building the initial k sampled trees. 
		for (int i = 0; i < k; ++i) {
			NodeID r = rand() % numOfVertices;
			if (picked_as_root[r] == true) { --i; continue; }

			st_roots[i] = r;
			large_tree &ltree_i = ltrees[i];
			picked_as_root[r] = true;
			picked_count++;
			used_as_cover++;

			build_large_tree_dij(r, ltree_i, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, 0, numOfVertices);
			st_sizes[i] = descendants.size();
			alive_resources += (long long)numOfVertices;

			// Turn large tree to small tree.
			// 1. Build new small tree with ltree_i and add it to strees;
			// 2. Set tree_pointers[i] to the position of the strees;
			// 3. Clear the large tree ltree_i by setting its size to 0.
			if (st_sizes[i] < size_thr)
				shrink_tree(i, st_roots[i], ltrees[i]);

			bool isLarge = false;
			int stree_idx = trees_pointers[i];
			if (stree_idx == -1)
				isLarge = true;

			small_tree& stree_i = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;

			// Calculate the subtree size of each nodes.
		//	unordered_map<int, int>& cover_i = cover[i];
			//vector<int> cover_i(numOfVertices, 0);

			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];

				++cover_tid[q]; // Increment because it hits the r-q path.
				++cnt16[q][i%CNT];
				//++sum[q];

				int parent_node;

				if (isLarge)
					parent_node = get_large_parent(q, ltree_i);
				else
					parent_node = get_small_parent(q, stree_i);

				if (parent_node > -1) {
					cover_tid[parent_node] += cover_tid[q];
					cnt16[parent_node][i%CNT] += cover_tid[q];
					//	sum[parent_node] += cover_i[q];
				}
			}
			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];
				cover_tid[q] = 0;
			}
		}



		//cnt16_aggreate(cnt16, has_picked, cover);

		vector<cover_value_type> picked_value;
		benchmark::heap<2, cover_value_type, NodeID> value_heap(numOfVertices);
		//bool switch_small = false;
		switch_small = false;

		//cnt16_2max(picked_value, sum, has_picked,affected, cnt16, switch_small, value_heap);
//		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, cnt16, switch_small, value_heap);
		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);

		init_time = GetCurrentTimeSec() - init_time;

		cout << "Finish initial sampling " << init_time << "secs." << endl;

		//for (NodeID v = 0; v < V; ++v)	value_heap.update(v, -picked_value[v]); // Min-heap.

		double max_value = -1;

		// Perform V iterations to rank all vertices.
		NodeID ranking = 0;
		for (; ranking < numOfVertices; ++ranking) {
			if (ranking % (numOfVertices / 10) == 0)
				cout << "picking the " << ranking << "-th orders." << endl;

			double time_this_iteration = GetCurrentTimeSec();

			NodeID pv = 0;
			long long pv_cover = -1;

			if (ranking == stopv) {
				stop_flag = true;
				cout << "topv labeling time:" << labeling_time << endl;
				break;
			}

			double selecting_acc_time = GetCurrentTimeSec();
			// In the beginning, pick the best cover vertex by a linear scan.
			if (switch_small == false) {
				double u1_acc_time = GetCurrentTimeSec();

				for (NodeID v = 0; v < numOfVertices; ++v) {
					if (has_picked[v] == true) continue;
					if (pv_cover == picked_value[v]) {
						if (sum[pv] < sum[v]) {
							pv_cover = picked_value[v];
							pv = v;
						}
					}
					else if (pv_cover < picked_value[v]) {
						pv_cover = picked_value[v];
						pv = v;
					}
				}

				if (sum[pv] < numOfVertices / 8 && switch_small == false) {
					//cout << "u3 time:" << u3_time << endl;
					cout << ranking << endl;
					cout << pv_cover << endl;
					cout << sum[pv] << endl;
					cout << picked_count << endl;
					cout << "Switch to max-heap representations." << endl;
					switch_small = true;
				}
				u1_time += (GetCurrentTimeSec() - u1_acc_time);

			}
			else {// Later, utilizing the max-heap directly
				double u2_acc_time = GetCurrentTimeSec();

				value_heap.extract_min(pv, pv_cover);
				pv_cover = -pv_cover;
				u2_time += (GetCurrentTimeSec() - u2_acc_time);
				/*
				int max1 = 0;
				int max2 = 0;
				int max1v = -1;
				int max2v = -1;
				for(int  i = 0; i < CNT; ++i){
					if(max2 <= cnt16[pv][i]){
						max2 = cnt16[pv][i];
						max2v = i;
					}
					if(max1 <= max2){
						int tmp = max1;
						max1 = max2;
						max2 = tmp;
						int tmpv = max1v;
						max1v = max2v;
						max2v = max1v;
					}
				}
				long long sumstd = 0;
				long long summean = 0;
				for(int  i = 0; i < CNT; ++i){
					if(i != max1v && i != max2v){
						sumstd += (long long)cnt16[pv][i] * (long long)cnt16[pv][i];
						summean += (long long)cnt16[pv][i];
					}
				}
				sumstd = sumstd/(long long)(CNT-2);
				summean = summean/(long long)(CNT-2);
				cout << sqrt(sumstd - summean * summean) << endl;*/
			}

			selecting_time += (GetCurrentTimeSec() - selecting_acc_time);

			//// FIX
			//if (switch_small == true && reverse_index.size() > 0 && sum[pv] <= reverse_index[pv].size() * 10) {//sum[pv] <= in_tree_count[pv] * 10){//
			//	cout << "I think we can stop here " << GetCurrentTimeSec() - start_time << " ranking:" << ranking << endl;
			//	break;
			//}

			//if (has_picked[pv] == true)
			//	cout << "wrong1" << endl;

			has_picked[pv] = true;
			double labeling_acc_time = GetCurrentTimeSec();
			labeled_arcs = labeling_source_dij(pv, pqueue, visited_que, distances, vis, dst_r, r_dst_r, usd, ranking, tmp_idx, r_tmp_idx, tmp_idx_token_parents, r_tmp_idx_token_parents, wgraph);
			labeling_time += (GetCurrentTimeSec() - labeling_acc_time);
			inv[ranking] = pv;
			rank[pv] = ranking;

			int current_trees_num = st_roots.size();

			// 4 steps in updating:
			double updating_acc_time = GetCurrentTimeSec();
			if (switch_small == false) { // We still iterate all alive trees to calculate the cover.	
			//	double u3_acc_time = GetCurrentTimeSec();
				for (int tid = 0; tid < current_trees_num; ++tid) {

					// Get the descendants and ascendants of pv in tree tid.
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					// pv has already been covered in tree tid.
					int pvpid = 0;
					if (isLarge)
						pvpid = get_large_parent(pv, ltree_tid);
					else {
						//double get_small_parent_incre_time = GetCurrentTimeSec();
						pvpid = get_small_parent(pv, stree_tid);
						//	get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}
					if (pvpid == -1 && st_roots[tid] != pv) // pv has been covered in tree tid 
						continue;


					double u3_acc_time = GetCurrentTimeSec();

					// //1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;				 
					//get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph.adj, isLarge, stree_tid, ltree_tid);
					//get_ascendants(tid, pv, ascendants, que, vis, wgraph.adj, wgraph.r_adj, isLarge, stree_tid, ltree_tid);

					//Array Representation
					get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph, isLarge, stree_tid, ltree_tid);
					get_ascendants(tid, pv, ascendants, que, vis, wgraph, isLarge, stree_tid, ltree_tid);

					u3_time += (GetCurrentTimeSec() - u3_acc_time);


					double u4_acc_time = GetCurrentTimeSec();
					NodeID subtract_size = descendants.size();
					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
					//	cover[tid][av] -= subtract_size;
						cnt16[av][tid%CNT] -= subtract_size;
						//sum[av] -= subtract_size;
						//affected[av] = true;

						if (affected_list[av] == false) {
							affected_list[av] = true;
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
						}
					}
					u4_time += (GetCurrentTimeSec() - u4_acc_time);


					double u5_acc_time = GetCurrentTimeSec();
					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					//vector<int> cover_tid(numOfVertices, 0);
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//sum[dv] -= cover_tid[dv];

						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];
						if (cnt16[dv][tid%CNT] < 0)
							cout << "a6 " << cnt16[dv][tid%CNT] << endl;

						/*if (isLarge)
							remove_large_node(dv, ltree_tid);
						else
							remove_small_node(dv, stree_tid);*/

							//						remove_node(tid, dv);
													//affected[dv] = true;

						if (affected_list[dv] == false) {
							affected_count++;
							affected_list[dv] = true;
							if (switch_set)
								affected_v.insert(dv);
						}
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;

						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else 
							remove_small_node(dv, stree_tid);
					}
					u5_time += (GetCurrentTimeSec() - u5_acc_time);

					// Drop the size of subtree.
					st_sizes[tid] -= subtract_size;
					if (trees_pointers[tid] != -1) alive_resources -= (long long)subtract_size;

					// 3. Convert shrink trees to small trees representataions.	
					if (st_sizes[tid] > 0) {
						if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
							shrink_tree(tid, st_roots[tid], ltrees[tid]);
						}
					}

					if (st_roots[tid] == pv)
						if (isLarge) {
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
							//ltrees.shrink_to_fit();
						}
						else {
							stree_tid.clear();
							stree_tid.resize(0);
							//strees.shrink_to_fit();
						}

				}
				//	u3_time += (GetCurrentTimeSec() - u3_acc_time);
			}
			else { // We use the reverse index to update trees.
			   // Build the reverse_index first (one time).
				if (reverse_index.size() == 0) {
					reverse_index.resize(numOfVertices);
					for (NodeID v = 0; v < numOfVertices; ++v) {
						reverse_index[v].set_empty_key(numOfVertices + 1);
						reverse_index[v].set_deleted_key(numOfVertices + 2);
					}
					for (int tid = 0; tid < current_trees_num; ++tid) {
						if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue; //This whole tree has already been covered. 

						bool isLarge = false;
						int stree_idx = trees_pointers[tid];
						if (stree_idx == -1)
							isLarge = true;

						small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
						large_tree& ltree_tid = ltrees[tid];

						for (int v = 0; v < numOfVertices; ++v) {

							int vpid = 0;
							if (isLarge)
								vpid = get_large_parent(v, ltree_tid);
							else {
								//double get_small_parent_incre_time = GetCurrentTimeSec();
								vpid = get_small_parent(v, stree_tid);
								//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
							}
							if ((vpid != -1 || st_roots[tid] == v)) {
								reverse_index[v].insert(tid);
								//in_tree_count[v]++;
							}
						}
					}
				}

				//cout << reverse_index[pv].size() << "," << sum[pv] << "," <<  (double)sum[pv] / (double)reverse_index[pv].size() << endl;
				// Access the alive trees directly by the reverse index.
				//for (unordered_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) {
				for (google::dense_hash_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) {
					int tid = *it;
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					double u6_acc_time = GetCurrentTimeSec();
					/*
							int parent_node;
							if (isLarge)
								parent_node = get_large_parent(pv, ltree_tid);
							else {
								//double get_small_parent_incre_time = GetCurrentTimeSec();
								parent_node = get_small_parent(pv, stree_tid);
								//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
							}
							if(parent_node == -1 && pv != st_roots[tid])
								continue;
						*/
					double u7_acc_time = GetCurrentTimeSec();

					//	// 1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;
					//	//double u4_acc_time = getcurrenttimesec();
					//	get_descendants(tid, pv, descendants, descendants_parents,que, vis, wgraph.adj, islarge, stree_tid, ltree_tid);
					////	u4_time += (getcurrenttimesec() - u4_acc_time);

					//	get_ascendants(tid, pv, ascendants, que, vis, wgraph.adj, wgraph.r_adj, isLarge, stree_tid, ltree_tid);

						// Array Representation
					get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph, isLarge, stree_tid, ltree_tid);
					get_ascendants(tid, pv, ascendants, que, vis, wgraph, isLarge, stree_tid, ltree_tid);

					NodeID subtract_size = descendants.size();
					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
						//cover[tid][av] -= subtract_size;
						cnt16[av][tid%CNT] -= subtract_size;
						//affected[av] = true;
					//	sum[av] -= subtract_size;

						if (affected_list[av] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
							affected_list[av] = true;
						}
					}
					u7_time += (GetCurrentTimeSec() - u7_acc_time);

					//double u2_acc_time = GetCurrentTimeSec();
					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					double u8_acc_time = GetCurrentTimeSec();
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//	sum[dv] -= cover_tid[dv];

						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];
						//	if (cnt16[dv][tid%CNT] < 0)
						//		cout << "a6 " << cnt16[dv][tid%CNT] << endl;

							//if(dv == pv){
								/*if (isLarge)
									remove_large_node(dv, ltree_tid);
								else
									remove_small_node(dv, stree_tid);*/
									//}
									//affected[dv] = true;

									//double u82_acc_time = GetCurrentTimeSec();
						if (dv != pv)
							reverse_index[dv].erase(tid);
						//if (dv != pv)
						//	in_tree_count[dv]--;

					//	double u81_acc_time = GetCurrentTimeSec();
						if (affected_list[dv] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(dv);
							affected_list[dv] = true;
						}

						//u81_time += (GetCurrentTimeSec() - u81_acc_time);

					//	u82_time += (GetCurrentTimeSec() - u82_acc_time);
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;
						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else 
							remove_small_node(dv, stree_tid);
					}
					u8_time += (GetCurrentTimeSec() - u8_acc_time);

					if (st_roots[tid] == pv)
						if (isLarge) {
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
						}
						else{
							stree_tid.clear();
							stree_tid.resize(0);
						}


					st_sizes[tid] -= subtract_size;
					if (trees_pointers[tid] != -1) alive_resources -= (long long)subtract_size;
					//	if (st_sizes[tid] < 0) {
					//		cout << "a5" << st_sizes[tid] << endl;
					//	}

						// 3. Convert shrink trees to small trees representataions.				
					if (st_sizes[tid] > 0) {
						if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
							shrink_tree(tid, st_roots[tid], ltrees[tid]);
						}
					}

					//u2_time += (GetCurrentTimeSec() - u2_acc_time);
				}
				reverse_index[pv].clear();
				reverse_index[pv].resize(0);

			}
			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			// 4. Sample new trees based on.
			// If we sample a new tree, we need to dynamically increment the following data structures:
			// 1) ltrees
			// 2) trees_pointers (-1 or pointer to strees)
			// 3) cover
			// 4) st_roots
			// 5) st_sizes 
			// 6) strees (maybe)
			// two constraints: cl > ct and 10kn total vertices
			double adding_acc_time = GetCurrentTimeSec();


			/*cout << "ranking:" << ranking << endl;
			cout << "pick cover:" << pv_cover << endl;
			cout << "pick count:" << picked_count << endl;
			cout << "labeled_arcs:" << labeled_arcs << endl;
			cout << "new_tree_arcs:" << new_tree_arcs << endl;
			cout << "alive_resources:" << alive_resources << endl;
			cout << "bounded_resources:" << bounded_resources << endl;
			cout << "#tree:" << cover.size() << endl;
			cout << "u5/adding" << u5_time / adding_time << endl;
*/
			new_tree_arcs = 0;
			int how_many = 0;
			while (beta * labeled_arcs > new_tree_arcs && alive_resources < bounded_resources) {

				if (picked_count == numOfVertices) break;

				NodeID r = rand() % numOfVertices;
				if (picked_as_root[r] == true) continue;
				if (trivial_root[r] == true) continue;

				large_tree new_ltree;
				small_tree new_stree;
				//GOOGLE DENSE MAP
				new_stree.set_empty_key(numOfVertices + 1);
				new_stree.set_deleted_key(numOfVertices + 2);
				bool isLarge = false;
				int tmp_tree_arcs = 0;

				if (!BUILD_SMALL_TREE_FLAG) {

					//double u5_acc_time = GetCurrentTimeSec();
					double u9_acc_time = GetCurrentTimeSec();
					tmp_tree_arcs = build_large_tree_dij(r, new_ltree, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, new_tree_arcs, beta * labeled_arcs);
					u9_time += (GetCurrentTimeSec() - u9_acc_time);
					//u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;

					if (descendants.size() < size_thr)
						BUILD_SMALL_TREE_FLAG = true;

					if (descendants.size() < 100) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					//	double u91_acc_time = GetCurrentTimeSec();

					isLarge = true;
					how_many++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					used_as_cover++;
					st_roots.push_back(r);
					trees_pointers.push_back(-1);
					st_sizes.push_back(descendants.size());

					//new_tree_arcs += descendants.size();
					alive_resources += (long long)numOfVertices;

					//	u91_time += (GetCurrentTimeSec() - u91_acc_time);
						// Turn large tree to small tree.
						// 1. Build new small tree with ltree_i and add it to strees;
						// 2. Set tree_pointers[i] to the position of the strees;
						// 3. Clear the large tree ltree_i by setting its size to 0.

					//	double u92_acc_time = GetCurrentTimeSec();
					int i = st_roots.size() - 1;
					if (st_sizes[i] < size_thr) {
						shrink_tree(i, st_roots[i], new_ltree);
						// FIX
						isLarge = false;
						/*for (int di = descendants.size() - 1; di > -1; --di) {
							NodeID q = descendants[di];
							NodeID v_parent = get_large_parent(q, new_ltree);
							NodeID s_parent = get_small_parent(q, strees[strees.size() - 1]);
							if (v_parent != s_parent) cout << v_parent << " vs " << s_parent << endl;
						}*/
					}
					ltrees.push_back(new_ltree);
					//	u92_time += (GetCurrentTimeSec() - u92_acc_time);
				}
				else {//Build Small Tree Directly
					//double u5_acc_time = GetCurrentTimeSec();
					double u9_acc_time = GetCurrentTimeSec();
					tmp_tree_arcs = build_small_tree_dij(r, new_stree, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, new_tree_arcs, beta * labeled_arcs);
					//u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;
					u9_time += (GetCurrentTimeSec() - u9_acc_time);

					if (descendants.size() >= size_thr)
						BUILD_SMALL_TREE_FLAG = false;

					if (descendants.size() < 100) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					isLarge = false;

					how_many++;
					used_as_cover++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					st_roots.push_back(r);
					trees_pointers.push_back(strees.size()); //we assume it is inserted.
					st_sizes.push_back(descendants.size());

					alive_resources += (long long)new_stree.size();
					//strees.push_back(new_stree);
					//// new_ltree is empty.
					//ltrees.push_back(new_ltree);

					int i = st_roots.size() - 1;
					if (st_sizes[i] >= size_thr) {
						enlarge_tree(i, st_roots[i], new_stree, new_ltree);
						new_stree.clear();
						new_stree.resize(0);
						isLarge = true;
					}
					else {
						strees.push_back(new_stree);
						ltrees.push_back(new_ltree);
					}

				}

				int i = st_roots.size() - 1;

				int numOfStrees = strees.size();
				small_tree &new_stree_ref = numOfStrees != 0 ? strees[strees.size() - 1] : empty_small_tree;

				// Calculate the subtree size of each nodes.

				double u10_acc_time = GetCurrentTimeSec();
				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];

					++cover_tid[q]; // Increment because it hitsthe r-q path.
					++cnt16[q][i%CNT];
					//++sum[q];
					//affected[q] = true;

					if (affected_list[q] == false) {
						affected_count++;
						if (switch_set)
							affected_v.insert(q);
						affected_list[q] = true;
					}

					if (switch_small == true) {
						reverse_index[q].insert(i);
						//in_tree_count[q]++;
					}


					int parent_node;
					if (isLarge)
						parent_node = get_large_parent(q, new_ltree);
					else {
						//double get_small_parent_incre_time = GetCurrentTimeSec();
						parent_node = get_small_parent(q, new_stree_ref);
						//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}
					if (parent_node > -1) {
						cover_tid[parent_node] += cover_tid[q];
						cnt16[parent_node][i%CNT] += cover_tid[q];
						//sum[parent_node] += cover_new[q];
					}
				}
				u10_time += (GetCurrentTimeSec() - u10_acc_time);
				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];
					cover_tid[q] = 0;
				}
				//cout << "beta*labeld_arcs: " << beta * labeled_arcs << " new tree _arcs: " << new_tree_arcs << " alive " << alive_resources << " bounded " << bounded_resources << endl;

			}
			adding_time += (GetCurrentTimeSec() - adding_acc_time);


			//double u3_acc_time = GetCurrentTimeSec();

			updating_acc_time = GetCurrentTimeSec();

			// Build the cover value heap for one time.
			if (switch_small == true && value_heap.empty())
				for (int nid = 0; nid < numOfVertices; ++nid)
					if (has_picked[nid] == false)
						value_heap.update(nid, -picked_value[nid]);

			//	cnt16_2max(picked_value, sum, has_picked, affected, cnt16, switch_small, value_heap);
		//c	cnt16_2max_queue(picked_value, sum, has_picked, affected_v, cnt16, switch_small, value_heap);

			double u11_acc_time = GetCurrentTimeSec();
			cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);
			u11_time += (GetCurrentTimeSec() - u11_acc_time);
			//	iteration_time[ranking] = (GetCurrentTimeSec() - updating_acc_time);
			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			if (affected_count < numOfVertices / 8)
				switch_set = true;
			else
				switch_set = false;
			affected_count = 0;

			time_this_iteration = GetCurrentTimeSec() - time_this_iteration;
			iteration_time[ranking] = time_this_iteration;

		}

		cout << "picked as tree covers:" << used_as_cover << endl;
		cout << "picked count:" << picked_count << endl;
		
		//// FIX
		//for (int v = 0; v < numOfVertices; ++v){
		//	for (int i = 0; i < CNT; ++i) {
		//		if (cnt16[v][i] != 0)
		//			cout << "mi?" << endl;
		//	}
		//}
		if (stop_flag == true) return;

		double labeling_acc_time = GetCurrentTimeSec();
		// Order the remaining vertices with degree, if any.
		if (ranking != numOfVertices) {
			vector<NodeID> degree(numOfVertices);
			vector<NodeID> r_degree(numOfVertices);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				/*vector<NodeID>& adj_v = graph.adj[v];
				for (NodeID i = 0; i < adj_v.size(); ++i) {
				if (has_picked[adj_v[i]] != true)
				degree[v]++;
				}*/

				for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid) {
					NodeID w = wgraph.edges[eid].first;
					if (has_picked[w] != true)
						degree[v]++;
				}

				if (DIRECTED_FLAG == true) {
					/*vector<NodeID>& r_adj_v = graph.r_adj[v];
					for (NodeID i = 0; i < r_adj_v.size(); ++i) {
					if (has_picked[r_adj_v[i]] != true)
					r_degree[v]++;
					}*/
					for (EdgeID eid = wgraph.r_vertices[v]; eid < wgraph.r_vertices[v + 1]; ++eid) {
						NodeID w = wgraph.r_edges[eid].first;
						if (has_picked[w] != true)
							r_degree[v]++;
					}
				}
			}
			
			vector<pair<NodeID, NodeID> > degree_v;
			degree_v.reserve(numOfVertices - ranking);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				if (DIRECTED_FLAG == true)
					degree_v.push_back(make_pair(degree[v] * r_degree[v], v));
				else
					degree_v.push_back(make_pair(degree[v], v));
			}
			sort(degree_v.rbegin(), degree_v.rend());

			int count = 0;
			for (; ranking < numOfVertices; ++ranking) {
				NodeID tpv = degree_v[count++].second;
				inv[ranking] = tpv;
				rank[tpv] = ranking;
				labeling_source_dij(tpv, pqueue, visited_que, distances, vis, dst_r, r_dst_r, usd, ranking, tmp_idx, r_tmp_idx, tmp_idx_token_parents, r_tmp_idx_token_parents,  wgraph);
			}
		}
		labeling_time += (GetCurrentTimeSec() - labeling_acc_time);


		// Finish the constructions.
		/*
		for (NodeID v = 0; v < numOfVertices; ++v) {
			if (DIRECTED_FLAG == false) {
				NodeID k = tmp_idx[v].first.size();
				labels.index_[v].spt_v.resize(k);
				labels.index_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();
				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();

			}
			else {
				NodeID k = tmp_idx[v].first.size();
				dlabels.index_[v].spt_v.resize(k);
				dlabels.index_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();
				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();

				k = r_tmp_idx[v].first.size();
				dlabels.bindex_[v].spt_v.resize(k);
				dlabels.bindex_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_v[i] = r_tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_d[i] = r_tmp_idx[v].second[i];
				r_tmp_idx[v].first.clear();
				r_tmp_idx[v].second.clear();
				r_tmp_idx[v].first.shrink_to_fit();
				r_tmp_idx[v].second.shrink_to_fit();
			}
		}*/

		converttokens(tmp_idx,  tmp_idx_token_parents, false);
		//converttokens(tmp_idx,  tmp_idx_token_parents, r_tmp_idx, r_tmp_idx_token_parents);
		cout << clabels.numOfTokens << " Tokens in total" << endl;
		cout << (double)children_size / (double) clabels.numOfTokens << " average children number" << endl;
		 
		if (DIRECTED_FLAG == true) {
			converttokens(r_tmp_idx,  r_tmp_idx_token_parents, true);
			cout << clabels.r_numOfTokens << " Tokens in total" << endl;
			cout << (double)r_children_size / (double) clabels.r_numOfTokens << " average children number" << endl;
		}
		
		num_of_trees = st_roots.size();

		double small_tree_num = 0;

		for (int i = 0; i < num_of_trees; ++i) {
			if (trees_pointers[i] != -1) {
				small_tree_num++;
				//small_tree_num += strees[trees_pointers[i]].size();
			}
		}
		
		total_resources = small_tree_num;
		
	}
	
	Betweenness_Ordering(NodeID k, double beta, WGraph& wgraph, NodeID stopv, long long bound_times, long long trans_times, double mem_bound) {

		double start_time = GetCurrentTimeSec();

		if (k < CNT) {
			cout << " Should pick at least 16 seed!" << endl;
			return;
		}
		iteration_time.resize(numOfVertices);

		updating_time = 0;
		labeling_time = 0;
		adding_time = 0;
		selecting_time = 0;
		init_time = 0;
		u1_time = 0;
		u2_time = 0;
		u3_time = 0;
		num_of_trees = 0;
		total_resources = 0;

		inv.resize(numOfVertices, -1);
		rank.resize(numOfVertices, -1);

		if (DIRECTED_FLAG == false) {
			labels.index_.resize(numOfVertices);
		}
		else { //if (DIRECTED_FLAG == true) {
			dlabels.index_.resize(numOfVertices);
			dlabels.bindex_.resize(numOfVertices);
		}

		// Resource Bounding.
		NodeID new_tree_arcs = 0;
		NodeID labeled_arcs = 0;
		alive_resources = 0;
		long long alive_memory = 0;
		bounded_resources = mem_bound * 1024 * 1024 * 1024;

		vector<NodeID> st_roots; // The root vertices of the sampled trees.
		vector<NodeID> st_sizes; // The size of the sampled trees.
		st_roots.reserve(numOfVertices);
		st_sizes.reserve(numOfVertices);
		alive_memory += numOfVertices * sizeof(NodeID) * 2;

		NodeID size_thr = numOfVertices < 8 ? numOfVertices : numOfVertices / trans_times; // When a tree shrink to a size smaller than the size_thr, it should be converted into a small_tree representation.

		init_time = GetCurrentTimeSec();
		// Building the initial k trees.
		//srand((unsigned)time(NULL));
		srand(0x728f9d0);

		st_roots.resize(k);
		st_sizes.resize(k);

		ltrees.reserve(numOfVertices);
		strees.reserve(numOfVertices);
		trees_pointers.reserve(numOfVertices);

		alive_memory += numOfVertices * sizeof(large_tree);
		alive_memory += numOfVertices * sizeof(small_tree);
		alive_memory += numOfVertices * sizeof(NodeID);

		ltrees.resize(k);
		strees.clear();
		trees_pointers.resize(k, -1); // Point the sampling trees to its tree representations. large tree [0, V) small tree [V, 2V).
									  //vector<unordered_map<int, int> > cover(k, unordered_map<int, int>());

		vector<NodeID> que(numOfVertices);
		vector<bool> vis(numOfVertices);

		// Additional Auxilarity Data.
		vector<NodeID> ascendants;
		ascendants.reserve(numOfVertices);
		vector<NodeID> descendants;
		descendants.reserve(numOfVertices);
		vector<bool> has_picked(numOfVertices);
		vector<vector<cover_value_type> > cnt16(numOfVertices, vector<cover_value_type>(CNT));
		//		vector<bool> affected(V, true);
		vector<bool> picked_as_root(numOfVertices, false);
		vector<bool> trivial_root(numOfVertices, false);
		int picked_count = 0;
		int used_as_cover = 0;
		vector<cover_value_type> sum(numOfVertices);
		bool stop_flag = false;

		alive_memory += numOfVertices * sizeof(NodeID) * 2;
		alive_memory += numOfVertices * sizeof(bool) * 3;
		alive_memory += numOfVertices * CNT * sizeof(cover_value_type);
		alive_memory += numOfVertices * sizeof(cover_value_type);

		// tmp structure used in the dijkstra searches. 

		// Binary Heap
		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);

		queue<NodeID> visited_que;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
		vector<bool> usd(numOfVertices, false);

		//GOOGLE
		//		unordered_set<NodeID> affected_v;
		google::dense_hash_set<NodeID> affected_v;
		affected_v.set_empty_key(numOfVertices + 1);
		affected_v.set_deleted_key(numOfVertices + 2);

		//alive memory for affected set.


		vector<cover_value_type> cover_tid(numOfVertices, 0);
		//vector<cover_value_type> cover_new(numOfVertices, 0); // modified-5-14
		vector<bool> affected_list(numOfVertices, false);
		int affected_count = 0;
		bool switch_set = false;

		alive_memory += numOfVertices * sizeof(cover_value_type);
		alive_memory += numOfVertices * sizeof(bool);


		//vector<unordered_set<int> > reverse_index(0);
		//GOOGLE DENSE SET
		vector<google::dense_hash_set<NodeID> > reverse_index;
		//vector<NodeID> in_tree_count(numOfVertices);

		// Preparing basic structure for pl algorithm.
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT))); // Backward labels.

		vector<NodeID> descendants_parents(numOfVertices, -1);

		alive_memory += numOfVertices * sizeof(NodeID);

		// Building the initial k sampled trees. 
		for (int i = 0; i < k; ++i) {
			NodeID r = rand() % numOfVertices;
			if (picked_as_root[r] == true) { --i; continue; }

			st_roots[i] = r;
			large_tree &ltree_i = ltrees[i];
			picked_as_root[r] = true;
			picked_count++;
			used_as_cover++;

			build_large_tree_dij(r, ltree_i, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, 0, numOfVertices);
			st_sizes[i] = descendants.size();
			alive_memory += numOfVertices * sizeof(NodeID);

			// Turn large tree to small tree.
			// 1. Build new small tree with ltree_i and add it to strees;
			// 2. Set tree_pointers[i] to the position of the strees;
			// 3. Clear the large tree ltree_i by setting its size to 0.
			if (st_sizes[i] < size_thr) {
				shrink_tree(i, st_roots[i], ltrees[i]);
				alive_memory -= numOfVertices * sizeof(NodeID);
				alive_memory += strees[trees_pointers[i]].bucket_count() * sizeof(NodeID) * 2;			
			}

			bool isLarge = false;
			int stree_idx = trees_pointers[i];
			if (stree_idx == -1)
				isLarge = true;

			small_tree& stree_i = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;

			// Calculate the subtree size of each nodes.
			//	unordered_map<int, int>& cover_i = cover[i];
			//vector<int> cover_i(numOfVertices, 0);

			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];

				++cover_tid[q]; // Increment because it hits the r-q path.
				++cnt16[q][i%CNT];
				//++sum[q];

				int parent_node;

				if (isLarge)
					parent_node = get_large_parent(q, ltree_i);
				else
					parent_node = get_small_parent(q, stree_i);

				if (parent_node > -1) {
					cover_tid[parent_node] += cover_tid[q];
					cnt16[parent_node][i%CNT] += cover_tid[q];
					//	sum[parent_node] += cover_i[q];
				}
			}
			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];
				cover_tid[q] = 0;
			}
			
			alive_memory += descendants.size() * sizeof(NodeID); // The size occupied by reverse index.
		}



		//cnt16_aggreate(cnt16, has_picked, cover);

		vector<cover_value_type> picked_value;
		benchmark::heap<2, cover_value_type, NodeID> value_heap(numOfVertices);
		//bool switch_small = false;
		switch_small = false;

		//cnt16_2max(picked_value, sum, has_picked,affected, cnt16, switch_small, value_heap);
		//		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, cnt16, switch_small, value_heap);
		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);

		init_time = GetCurrentTimeSec() - init_time;

		cout << "Finish initial sampling " << init_time << "secs." << endl;

		//for (NodeID v = 0; v < V; ++v)	value_heap.update(v, -picked_value[v]); // Min-heap.

		double max_value = -1;

		// Perform V iterations to rank all vertices.
		NodeID ranking = 0;
		for (; ranking < numOfVertices; ++ranking) {
			if (ranking % (numOfVertices / 10) == 0)
				cout << "picking the " << ranking << "-th orders." << endl;

			double time_this_iteration = GetCurrentTimeSec();

			NodeID pv = 0;
			long long pv_cover = -1;

			if (ranking == stopv) {
				stop_flag = true;
				cout << "topv labeling time:" << labeling_time << endl;
				break;
			}

			double selecting_acc_time = GetCurrentTimeSec();
			// In the beginning, pick the best cover vertex by a linear scan.
			if (switch_small == false) {
				double u1_acc_time = GetCurrentTimeSec();

				for (NodeID v = 0; v < numOfVertices; ++v) {
					if (has_picked[v] == true) continue;
					if (pv_cover == picked_value[v]) {
						if (sum[pv] < sum[v]) {
							pv_cover = picked_value[v];
							pv = v;
						}
					}
					else if (pv_cover < picked_value[v]) {
						pv_cover = picked_value[v];
						pv = v;
					}
				}

				/*
				int max1 = 0;
				int max2 = 0;
				int max1v = -1;
				int max2v = -1;
				for(int  i = 0; i < CNT; ++i){
				if(max2 <= cnt16[pv][i]){
				max2 = cnt16[pv][i];
				max2v = i;
				}
				if(max1 <= max2){
				int tmp = max1;
				max1 = max2;
				max2 = tmp;
				int tmpv = max1v;
				max1v = max2v;
				max2v = max1v;
				}
				}
				long long sumstd = 0;
				long long summean = 0;
				for(int  i = 0; i < CNT; ++i){
				if(i != max1v && i != max2v){
				sumstd += (long long)cnt16[pv][i] * (long long)cnt16[pv][i];
				summean += (long long)cnt16[pv][i];
				}
				}
				sumstd = sumstd/(long long)(CNT-2);
				summean = summean/(long long)(CNT-2);
				cout << sqrt(sumstd - summean * summean) << endl;*/

				if (sum[pv] < numOfVertices / 8 && switch_small == false) {
					//cout << "u3 time:" << u3_time << endl;
					cout << ranking << endl;
					cout << pv_cover << endl;
					cout << sum[pv] << endl;
					cout << picked_count << endl;
					cout << "Switch to max-heap representations." << endl;
					switch_small = true;
				}
				u1_time += (GetCurrentTimeSec() - u1_acc_time);

			}
			else {// Later, utilizing the max-heap directly
				double u2_acc_time = GetCurrentTimeSec();

				value_heap.extract_min(pv, pv_cover);
				pv_cover = -pv_cover;
				u2_time += (GetCurrentTimeSec() - u2_acc_time);
				/*
				int max1 = 0;
				int max2 = 0;
				int max1v = -1;
				int max2v = -1;
				for(int  i = 0; i < CNT; ++i){
				if(max2 <= cnt16[pv][i]){
				max2 = cnt16[pv][i];
				max2v = i;
				}
				if(max1 <= max2){
				int tmp = max1;
				max1 = max2;
				max2 = tmp;
				int tmpv = max1v;
				max1v = max2v;
				max2v = max1v;
				}
				}
				long long sumstd = 0;
				long long summean = 0;
				for(int  i = 0; i < CNT; ++i){
				if(i != max1v && i != max2v){
				sumstd += (long long)cnt16[pv][i] * (long long)cnt16[pv][i];
				summean += (long long)cnt16[pv][i];
				}
				}
				sumstd = sumstd/(long long)(CNT-2);
				summean = summean/(long long)(CNT-2);
				cout << sqrt(sumstd - summean * summean) << endl;*/
			}

			selecting_time += (GetCurrentTimeSec() - selecting_acc_time);

			// Memcon
			if (switch_small == true && reverse_index.size() > 0 && sum[pv] <= reverse_index[pv].size() * 10) {//sum[pv] <= in_tree_count[pv] * 10){//
				cout << "I think we can stop here " << GetCurrentTimeSec() - start_time << " ranking:" << ranking << endl;
				break;
			}
			// memcon
			/*if (switch_small&& reverse_index.size() > 0)
				cout << ranking << "," << alive_memory << "," << bounded_resources << "," << picked_count << "," << sum[pv] << "," << reverse_index[pv].size() << endl;*/


			/*if (pv_cover == 0 && picked_count == numOfVertices) {
				cout << "time to stop" << endl;
				break;
			}*/

			//if (has_picked[pv] == true)
			//	cout << "wrong1" << endl;

			has_picked[pv] = true;
			double labeling_acc_time = GetCurrentTimeSec();
			labeled_arcs = labeling_source_dij(pv, pqueue, visited_que, distances, vis, dst_r, usd, ranking, tmp_idx, r_tmp_idx, wgraph);
			labeling_time += (GetCurrentTimeSec() - labeling_acc_time);
			inv[ranking] = pv;
			rank[pv] = ranking;
			int current_trees_num = st_roots.size();

			// 4 steps in updating:
			double updating_acc_time = GetCurrentTimeSec();
			if (switch_small == false) { // We still iterate all alive trees to calculate the cover.	
										 //	double u3_acc_time = GetCurrentTimeSec();
				for (int tid = 0; tid < current_trees_num; ++tid) {

					// Get the descendants and ascendants of pv in tree tid.
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					// pv has already been covered in tree tid.
					int pvpid = 0;
					if (isLarge)
						pvpid = get_large_parent(pv, ltree_tid);
					else {
						//double get_small_parent_incre_time = GetCurrentTimeSec();
						pvpid = get_small_parent(pv, stree_tid);
						//	get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}
					if (pvpid == -1 && st_roots[tid] != pv) // pv has been covered in tree tid 
						continue;


					double u3_acc_time = GetCurrentTimeSec();

					// //1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;				 
					//get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph.adj, isLarge, stree_tid, ltree_tid);
					//get_ascendants(tid, pv, ascendants, que, vis, wgraph.adj, wgraph.r_adj, isLarge, stree_tid, ltree_tid);

					//Array Representation
					get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph, isLarge, stree_tid, ltree_tid);
					get_ascendants(tid, pv, ascendants, que, vis, wgraph, isLarge, stree_tid, ltree_tid);

					u3_time += (GetCurrentTimeSec() - u3_acc_time);


					double u4_acc_time = GetCurrentTimeSec();
					NodeID subtract_size = descendants.size();
					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
												//	cover[tid][av] -= subtract_size;
						cnt16[av][tid%CNT] -= subtract_size;
						//sum[av] -= subtract_size;
						//affected[av] = true;

						if (affected_list[av] == false) {
							affected_list[av] = true;
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
						}
					}
					u4_time += (GetCurrentTimeSec() - u4_acc_time);


					double u5_acc_time = GetCurrentTimeSec();
					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					//vector<int> cover_tid(numOfVertices, 0);
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//sum[dv] -= cover_tid[dv];

						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];
						/*if (cnt16[dv][tid%CNT] < 0)
							cout << "a6 " << cnt16[dv][tid%CNT] << endl;*/

						/*if (isLarge)
						remove_large_node(dv, ltree_tid);
						else
						remove_small_node(dv, stree_tid);*/

						//						remove_node(tid, dv);
						//affected[dv] = true;

						if (affected_list[dv] == false) {
							affected_count++;
							affected_list[dv] = true;
							if (switch_set)
								affected_v.insert(dv);
						}
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;

						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else {
							alive_memory -= stree_tid.bucket_count() * sizeof(NodeID) * 2;
							remove_small_node(dv, stree_tid);
							stree_tid.resize(0);
							alive_memory += stree_tid.bucket_count() * sizeof(NodeID) * 2;
						}
					}
					u5_time += (GetCurrentTimeSec() - u5_acc_time);


					alive_memory -= descendants.size() * sizeof(NodeID); // The size released by reverse index.

					// Drop the size of subtree.
					st_sizes[tid] -= subtract_size;
					//if (trees_pointers[tid] != -1) alive_memory -= subtract_size * sizeof(NodeID) * 2;

					// 3. Convert shrink trees to small trees representataions.	
					if (st_sizes[tid] > 0) {
						if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
							shrink_tree(tid, st_roots[tid], ltrees[tid]);
							alive_memory -= numOfVertices * sizeof(NodeID);
							alive_memory += strees[trees_pointers[tid]].bucket_count() * sizeof(NodeID) * 2;
						}
					}

					if (st_roots[tid] == pv)
						if (isLarge) {
							alive_memory -= numOfVertices * sizeof(NodeID);
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
							//ltrees.shrink_to_fit();
						}
						else {
							alive_memory -= stree_tid.bucket_count() * sizeof(NodeID) * 2;
							stree_tid.clear();
							stree_tid.resize(0);
							//strees.shrink_to_fit();
						}

				}
				//	u3_time += (GetCurrentTimeSec() - u3_acc_time);
			}
			else { // We use the reverse index to update trees.
				   // Build the reverse_index first (one time).
				if (reverse_index.size() == 0) {
					reverse_index.resize(numOfVertices);
					for (NodeID v = 0; v < numOfVertices; ++v) {
						reverse_index[v].set_empty_key(numOfVertices + 1);
						reverse_index[v].set_deleted_key(numOfVertices + 2);
					}
					for (int tid = 0; tid < current_trees_num; ++tid) {
						if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue; //This whole tree has already been covered. 

						bool isLarge = false;
						int stree_idx = trees_pointers[tid];
						if (stree_idx == -1)
							isLarge = true;

						small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
						large_tree& ltree_tid = ltrees[tid];

						for (int v = 0; v < numOfVertices; ++v) {

							int vpid = 0;
							if (isLarge)
								vpid = get_large_parent(v, ltree_tid);
							else {
								//double get_small_parent_incre_time = GetCurrentTimeSec();
								vpid = get_small_parent(v, stree_tid);
								//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
							}
							if ((vpid != -1 || st_roots[tid] == v)) {
								reverse_index[v].insert(tid);
								//in_tree_count[v]++;
							}
						}
					}
				}

				//cout << reverse_index[pv].size() << "," << sum[pv] << "," <<  (double)sum[pv] / (double)reverse_index[pv].size() << endl;
				// Access the alive trees directly by the reverse index.
				//for (unordered_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) {
				for (google::dense_hash_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) {
					int tid = *it;
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					double u6_acc_time = GetCurrentTimeSec();
					/*
					int parent_node;
					if (isLarge)
					parent_node = get_large_parent(pv, ltree_tid);
					else {
					//double get_small_parent_incre_time = GetCurrentTimeSec();
					parent_node = get_small_parent(pv, stree_tid);
					//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}
					if(parent_node == -1 && pv != st_roots[tid])
					continue;
					*/
					double u7_acc_time = GetCurrentTimeSec();

					//	// 1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;
					//	//double u4_acc_time = getcurrenttimesec();
					//	get_descendants(tid, pv, descendants, descendants_parents,que, vis, wgraph.adj, islarge, stree_tid, ltree_tid);
					////	u4_time += (getcurrenttimesec() - u4_acc_time);

					//	get_ascendants(tid, pv, ascendants, que, vis, wgraph.adj, wgraph.r_adj, isLarge, stree_tid, ltree_tid);

					// Array Representation
					get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph, isLarge, stree_tid, ltree_tid);
					get_ascendants(tid, pv, ascendants, que, vis, wgraph, isLarge, stree_tid, ltree_tid);

					NodeID subtract_size = descendants.size();
					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
												//cover[tid][av] -= subtract_size;
						cnt16[av][tid%CNT] -= subtract_size;
						//affected[av] = true;
						//	sum[av] -= subtract_size;

						if (affected_list[av] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
							affected_list[av] = true;
						}
					}
					u7_time += (GetCurrentTimeSec() - u7_acc_time);

					//double u2_acc_time = GetCurrentTimeSec();
					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					double u8_acc_time = GetCurrentTimeSec();
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//	sum[dv] -= cover_tid[dv];

						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];
						//	if (cnt16[dv][tid%CNT] < 0)
						//		cout << "a6 " << cnt16[dv][tid%CNT] << endl;

						//if(dv == pv){
						/*if (isLarge)
						remove_large_node(dv, ltree_tid);
						else
						remove_small_node(dv, stree_tid);*/
						//}
						//affected[dv] = true;

						//double u82_acc_time = GetCurrentTimeSec();
						if (dv != pv)
							reverse_index[dv].erase(tid);
						//if (dv != pv)
						//	in_tree_count[dv]--;

						//	double u81_acc_time = GetCurrentTimeSec();
						if (affected_list[dv] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(dv);
							affected_list[dv] = true;
						}

						//u81_time += (GetCurrentTimeSec() - u81_acc_time);

						//	u82_time += (GetCurrentTimeSec() - u82_acc_time);
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;
						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else{
							alive_memory -= stree_tid.bucket_count() * sizeof(NodeID) * 2;
							remove_small_node(dv, stree_tid);
							stree_tid.resize(0);
							alive_memory += stree_tid.bucket_count() * sizeof(NodeID) * 2;
						}
					}
					u8_time += (GetCurrentTimeSec() - u8_acc_time);

					if (st_roots[tid] == pv)
						if (isLarge) {
							alive_memory -= ltree_tid.size() * sizeof(NodeID);
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
						}
						else {
							alive_memory -= stree_tid.bucket_count() * sizeof(NodeID) * 2;
							stree_tid.clear();
							stree_tid.resize(0);
						}


					//alive_memory -= descendants.size() * sizeof(NodeID); // The size released by reverse index.
					st_sizes[tid] -= subtract_size;
					//if (trees_pointers[tid] != -1) alive_memory -= subtract_size * sizeof(NodeID) * 2;
					//	if (st_sizes[tid] < 0) {
					//		cout << "a5" << st_sizes[tid] << endl;
					//	}

					// 3. Convert shrink trees to small trees representataions.				
					if (st_sizes[tid] > 0) {
						if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
							shrink_tree(tid, st_roots[tid], ltrees[tid]);
							alive_memory -= numOfVertices * sizeof(NodeID);
							alive_memory += strees[trees_pointers[tid]].size() * sizeof(NodeID) * 2;
						}
					}

					//u2_time += (GetCurrentTimeSec() - u2_acc_time);
				}
				alive_memory -= reverse_index[pv].bucket_count() * sizeof(NodeID);
				reverse_index[pv].clear();
				reverse_index[pv].resize(0);

			}
			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			// 4. Sample new trees based on.
			// If we sample a new tree, we need to dynamically increment the following data structures:
			// 1) ltrees
			// 2) trees_pointers (-1 or pointer to strees)
			// 3) cover
			// 4) st_roots
			// 5) st_sizes 
			// 6) strees (maybe)
			// two constraints: cl > ct and 10kn total vertices
			double adding_acc_time = GetCurrentTimeSec();


			/*cout << "ranking:" << ranking << endl;
			cout << "pick cover:" << pv_cover << endl;
			cout << "pick count:" << picked_count << endl;
			cout << "labeled_arcs:" << labeled_arcs << endl;
			cout << "new_tree_arcs:" << new_tree_arcs << endl;
			cout << "alive_resources:" << alive_resources << endl;
			cout << "bounded_resources:" << bounded_resources << endl;
			cout << "#tree:" << cover.size() << endl;
			cout << "u5/adding" << u5_time / adding_time << endl;
			*/
			new_tree_arcs = 0;
			int how_many = 0;
			while (alive_memory < bounded_resources) {

				if (picked_count == numOfVertices) break;

				NodeID r = rand() % numOfVertices;
				if (picked_as_root[r] == true) continue;
				if (trivial_root[r] == true) continue;

				large_tree new_ltree;
				small_tree new_stree;
				//GOOGLE DENSE MAP
				new_stree.set_empty_key(numOfVertices + 1);
				new_stree.set_deleted_key(numOfVertices + 2);
				bool isLarge = false;
				int tmp_tree_arcs = 0;

				if (!BUILD_SMALL_TREE_FLAG) {

					//double u5_acc_time = GetCurrentTimeSec();
					double u9_acc_time = GetCurrentTimeSec();
					tmp_tree_arcs = build_large_tree_dij(r, new_ltree, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, 0, numOfVertices);
					u9_time += (GetCurrentTimeSec() - u9_acc_time);
					//u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;

					if (descendants.size() <size_thr)
						BUILD_SMALL_TREE_FLAG = true;

					if (descendants.size() < 100) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					//	double u91_acc_time = GetCurrentTimeSec();

					isLarge = true;
					how_many++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					used_as_cover++;
					st_roots.push_back(r);
					trees_pointers.push_back(-1);
					st_sizes.push_back(descendants.size());

					//new_tree_arcs += descendants.size();
					//alive_resources += (long long)numOfVertices;


					alive_memory += numOfVertices * sizeof(NodeID);

					//	u91_time += (GetCurrentTimeSec() - u91_acc_time);
					// Turn large tree to small tree.
					// 1. Build new small tree with ltree_i and add it to strees;
					// 2. Set tree_pointers[i] to the position of the strees;
					// 3. Clear the large tree ltree_i by setting its size to 0.

					//	double u92_acc_time = GetCurrentTimeSec();
					int i = st_roots.size() - 1;
					if (st_sizes[i] < size_thr) {
						shrink_tree(i, st_roots[i], new_ltree);
						isLarge = false;
						alive_memory -= numOfVertices * sizeof(NodeID);
						alive_memory += strees[trees_pointers[i]].bucket_count() * sizeof(NodeID) * 2;
					
					}
					ltrees.push_back(new_ltree);
					//	u92_time += (GetCurrentTimeSec() - u92_acc_time);


				}
				else {//Build Small Tree Directly
					  //double u5_acc_time = GetCurrentTimeSec();
					double u9_acc_time = GetCurrentTimeSec();
					tmp_tree_arcs = build_small_tree_dij(r, new_stree, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, alive_memory/(2*sizeof(NodeID)), bounded_resources/(2*sizeof(NodeID)) );
					//u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;
					u9_time += (GetCurrentTimeSec() - u9_acc_time);

					if (descendants.size() >= size_thr)
						BUILD_SMALL_TREE_FLAG = false;

					if (descendants.size() < 100) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					isLarge = false;

					how_many++;
					used_as_cover++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					st_roots.push_back(r);
					trees_pointers.push_back(strees.size()); //we assume it is inserted.
					st_sizes.push_back(descendants.size());

					alive_memory += new_stree.bucket_count() * sizeof(NodeID) * 2;

					//alive_resources += (long long)new_stree.size();
					//strees.push_back(new_stree);
					//// new_ltree is empty.
					//ltrees.push_back(new_ltree);

					int i = st_roots.size() - 1;
					if (st_sizes[i] >= size_thr) {
						enlarge_tree(i, st_roots[i], new_stree, new_ltree);
						alive_memory -= new_stree.bucket_count() * sizeof(NodeID) * 2;
						alive_memory += numOfVertices * sizeof(NodeID);
						new_stree.clear();
						new_stree.resize(0);
						isLarge = true;
					}
					else {
						strees.push_back(new_stree);
						ltrees.push_back(new_ltree);
					}

				}

				int i = st_roots.size() - 1;

				int numOfStrees = strees.size();
				small_tree &new_stree_ref = numOfStrees != 0 ? strees[strees.size() - 1] : empty_small_tree;

				// Calculate the subtree size of each nodes.

				double u10_acc_time = GetCurrentTimeSec();
				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];

					++cover_tid[q]; // Increment because it hitsthe r-q path.
					++cnt16[q][i%CNT];
					//++sum[q];
					//affected[q] = true;

					if (affected_list[q] == false) {
						affected_count++;
						if (switch_set)
							affected_v.insert(q);
						affected_list[q] = true;
					}

					if (switch_small == true) {
						reverse_index[q].insert(i);
						//in_tree_count[q]++;
					}


					int parent_node;
					if (isLarge)
						parent_node = get_large_parent(q, new_ltree);
					else {
						//double get_small_parent_incre_time = GetCurrentTimeSec();
						parent_node = get_small_parent(q, new_stree_ref);
						//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}
					if (parent_node > -1) {
						cover_tid[parent_node] += cover_tid[q];
						cnt16[parent_node][i%CNT] += cover_tid[q];
						//sum[parent_node] += cover_new[q];
					}
				}
				u10_time += (GetCurrentTimeSec() - u10_acc_time);
				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];
					cover_tid[q] = 0;
				}


				alive_memory += descendants.size() * sizeof(NodeID); // The size occupied by reverse index.

				//cout << "beta*labeld_arcs: " << beta * labeled_arcs << " new tree _arcs: " << new_tree_arcs << " alive " << alive_resources << " bounded " << bounded_resources << endl;

			}
			//// memcon
			//if(switch_small)
			//	cout << ranking << "," << alive_memory << "," << bounded_resources << "," << picked_count << "," << sum[pv] << "," << reverse_index[pv].size() << endl;
			adding_time += (GetCurrentTimeSec() - adding_acc_time);


			//double u3_acc_time = GetCurrentTimeSec();

			updating_acc_time = GetCurrentTimeSec();

			// Build the cover value heap for one time.
			if (switch_small == true && value_heap.empty())
				for (int nid = 0; nid < numOfVertices; ++nid)
					if (has_picked[nid] == false)
						value_heap.update(nid, -picked_value[nid]);

			//	cnt16_2max(picked_value, sum, has_picked, affected, cnt16, switch_small, value_heap);
			//c	cnt16_2max_queue(picked_value, sum, has_picked, affected_v, cnt16, switch_small, value_heap);

			double u11_acc_time = GetCurrentTimeSec();
			cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);
			u11_time += (GetCurrentTimeSec() - u11_acc_time);
			//	iteration_time[ranking] = (GetCurrentTimeSec() - updating_acc_time);
			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			if (affected_count < numOfVertices / 8)
				switch_set = true;
			else
				switch_set = false;
			affected_count = 0;

			time_this_iteration = GetCurrentTimeSec() - time_this_iteration;
			iteration_time[ranking] = time_this_iteration;

		}

		cout << "picked as tree covers:" << used_as_cover << endl;
		cout << "picked count:" << picked_count << endl;
		//for (int v = 0; v < V; ++v)
		//	if (has_picked[v] == false) {
		//		cout << "wrong 2" << endl;
		//		cout << picked_value[v] << endl;
		//	}

		if (stop_flag == true) return;

		double labeling_acc_time = GetCurrentTimeSec();
		// Order the remaining vertices with degree, if any.
		if (ranking != numOfVertices) {
			vector<NodeID> degree(numOfVertices);
			vector<NodeID> r_degree(numOfVertices);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				/*vector<NodeID>& adj_v = graph.adj[v];
				for (NodeID i = 0; i < adj_v.size(); ++i) {
				if (has_picked[adj_v[i]] != true)
				degree[v]++;
				}*/

				for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid) {
					NodeID w = wgraph.edges[eid].first;
					if (has_picked[w] != true)
						degree[v]++;
				}

				if (DIRECTED_FLAG == true) {
					/*vector<NodeID>& r_adj_v = graph.r_adj[v];
					for (NodeID i = 0; i < r_adj_v.size(); ++i) {
					if (has_picked[r_adj_v[i]] != true)
					r_degree[v]++;
					}*/
					for (EdgeID eid = wgraph.r_vertices[v]; eid < wgraph.r_vertices[v + 1]; ++eid) {
						NodeID w = wgraph.r_edges[eid].first;
						if (has_picked[w] != true)
							r_degree[v]++;
					}
				}
			}

			vector<pair<NodeID, NodeID> > degree_v;
			degree_v.reserve(numOfVertices - ranking);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				if (DIRECTED_FLAG == true)
					degree_v.push_back(make_pair(degree[v] * r_degree[v], v));
				else
					degree_v.push_back(make_pair(degree[v], v));
			}
			sort(degree_v.rbegin(), degree_v.rend());

			int count = 0;
			for (; ranking < numOfVertices; ++ranking) {
				NodeID tpv = degree_v[count++].second;
				inv[ranking] = tpv;
				rank[tpv] = ranking;
				labeling_source_dij(tpv, pqueue, visited_que, distances, vis, dst_r, usd, ranking, tmp_idx, r_tmp_idx, wgraph);
			}
		}
		labeling_time += (GetCurrentTimeSec() - labeling_acc_time);


		// Finish the constructions.
		for (NodeID v = 0; v < numOfVertices; ++v) {
			if (DIRECTED_FLAG == false) {
				NodeID k = tmp_idx[v].first.size();
				labels.index_[v].spt_v.resize(k);
				labels.index_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();
				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();

			}
			else {
				NodeID k = tmp_idx[v].first.size();
				dlabels.index_[v].spt_v.resize(k);
				dlabels.index_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();
				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();

				k = r_tmp_idx[v].first.size();
				dlabels.bindex_[v].spt_v.resize(k);
				dlabels.bindex_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_v[i] = r_tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_d[i] = r_tmp_idx[v].second[i];
				r_tmp_idx[v].first.clear();
				r_tmp_idx[v].second.clear();
				r_tmp_idx[v].first.shrink_to_fit();
				r_tmp_idx[v].second.shrink_to_fit();
			}
		}

		num_of_trees = st_roots.size();

		double small_tree_num = 0;

		for (int i = 0; i < num_of_trees; ++i) {
			if (trees_pointers[i] != -1) {
				small_tree_num++;
				//small_tree_num += strees[trees_pointers[i]].size();
			}
		}

		total_resources = small_tree_num;

	}

	Betweenness_Ordering(NodeID k, double beta, WGraph& wgraph, NodeID stopv, long long bound_times, long long trans_times, double mem_bound, double c_paras) {

		double start_time = GetCurrentTimeSec();

		if (k < CNT) {
			cout << " Should pick at least 16 seed!" << endl;
			return;
		}
		iteration_time.resize(numOfVertices);

		updating_time = 0;
		labeling_time = 0;
		adding_time = 0;
		selecting_time = 0;
		init_time = 0;
		u1_time = 0;
		u2_time = 0;
		u3_time = 0;
		num_of_trees = 0;
		total_resources = 0;

		inv.resize(numOfVertices, -1);
		rank.resize(numOfVertices, -1);

		if (DIRECTED_FLAG == false) {
			labels.index_.resize(numOfVertices);
		}
		else { //if (DIRECTED_FLAG == true) {
			dlabels.index_.resize(numOfVertices);
			dlabels.bindex_.resize(numOfVertices);
		}

		vector<pair<NodeID, NodeID> > nodeid_degree(numOfVertices);
		vector<NodeID> d_order(numOfVertices);
		double refactor_by_jaccard = 1;
		for (NodeID v = 0; v < numOfVertices; ++v) {
			if (DIRECTED_FLAG == false)
				nodeid_degree[v] = make_pair(wgraph.vertices[v + 1] - wgraph.vertices[v], v);
			else
				if (wgraph.vertices[v + 1] - wgraph.vertices[v] != 0 && wgraph.vertices[v + 1] - wgraph.r_vertices[v] != 0)
					nodeid_degree[v] = make_pair((wgraph.vertices[v + 1] - wgraph.vertices[v]) * (wgraph.r_vertices[v + 1] - wgraph.r_vertices[v]), v);
				else if (wgraph.vertices[v + 1] - wgraph.vertices[v] != 0)
					nodeid_degree[v] = make_pair(wgraph.vertices[v + 1] - wgraph.vertices[v], v);
				else
					nodeid_degree[v] = make_pair(wgraph.r_vertices[v + 1] - wgraph.r_vertices[v], v);
		}
		sort(nodeid_degree.rbegin(), nodeid_degree.rend());
		for (NodeID v = 0; v < numOfVertices; ++v) {
			d_order[v] = nodeid_degree[v].second;
		}
		unordered_set<NodeID> top_degree, top_bo, top_union;
		NodeID intersection_count = 0;

		// Resource Bounding.
		NodeID new_tree_arcs = 0;
		NodeID labeled_arcs = 0;
		alive_resources = 0;
		bounded_resources = bound_times * (long long)k * (long long)numOfVertices;

		vector<NodeID> st_roots; // The root vertices of the sampled trees.
		vector<NodeID> st_sizes; // The size of the sampled trees.

		NodeID size_thr = numOfVertices < 8 ? numOfVertices : numOfVertices / trans_times; // When a tree shrink to a size smaller than the size_thr, it should be converted into a small_tree representation.

		init_time = GetCurrentTimeSec();
		// Building the initial k trees.
		//srand((unsigned)time(NULL));

		srand(0x728f9d0);

		// Reserve pointers of vectors.
		ltrees.reserve(numOfVertices);
		strees.reserve(numOfVertices);
		trees_pointers.reserve(numOfVertices);
		st_roots.reserve(numOfVertices);
		st_sizes.reserve(numOfVertices);

		st_roots.resize(k);
		st_sizes.resize(k);

		ltrees.resize(k);
		strees.clear();
		trees_pointers.resize(k, -1); // Point the sampling trees to its tree representations. large tree [0, V) small tree [V, 2V).
									  //vector<unordered_map<int, int> > cover(k, unordered_map<int, int>());

		vector<NodeID> ascendants;
		ascendants.reserve(numOfVertices);
		vector<NodeID> descendants;
		descendants.reserve(numOfVertices);
		vector<NodeID> que(numOfVertices);
		vector<bool> vis(numOfVertices);
		vector<bool> has_picked(numOfVertices);
		vector<vector<cover_value_type> > cnt16(numOfVertices, vector<cover_value_type>(CNT));
		//		vector<bool> affected(V, true);
		vector<bool> picked_as_root(numOfVertices, false);
		vector<bool> trivial_root(numOfVertices, false);
		int picked_count = 0;
		int used_as_cover = 0;
		vector<cover_value_type> sum(numOfVertices);
		bool stop_flag = false;


		// tmp structure used in the dijkstra searches. 

		// Binary Heap
		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);

		queue<NodeID> visited_que;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
		vector<bool> usd(numOfVertices, false);

		//GOOGLE
		//		unordered_set<NodeID> affected_v;
		google::dense_hash_set<NodeID> affected_v;
		affected_v.set_empty_key(numOfVertices + 1);
		affected_v.set_deleted_key(numOfVertices + 2);

		vector<cover_value_type> cover_tid(numOfVertices, 0);
		//vector<cover_value_type> cover_new(numOfVertices, 0); // modified-5-14


		vector<bool> affected_list(numOfVertices, false);
		int affected_count = 0;
		bool switch_set = false;

		//vector<unordered_set<int> > reverse_index(0);
		//GOOGLE DENSE SET
		vector<google::dense_hash_set<NodeID> > reverse_index;
		reverse_index.reserve(numOfVertices);
		//vector<NodeID> in_tree_count(numOfVertices);

		// Preparing basic structure for pl algorithm.
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT))); // Backward labels.

		vector<NodeID> descendants_parents(numOfVertices, -1);


		cover_value_type current_universe = 0;
		cover_value_type covered_universe = 0;

		// Building the initial k sampled trees. 
		for (int i = 0; i < k; ++i) {
			NodeID r = rand() % numOfVertices;
			if (picked_as_root[r] == true) { --i; continue; }

			st_roots[i] = r;
			large_tree &ltree_i = ltrees[i];
			picked_as_root[r] = true;
			picked_count++;
			used_as_cover++;

			build_large_tree_dij(r, ltree_i, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, 0, numOfVertices);
			st_sizes[i] = descendants.size();
			alive_resources += (long long)numOfVertices;

			// Turn large tree to small tree.
			// 1. Build new small tree with ltree_i and add it to strees;
			// 2. Set tree_pointers[i] to the position of the strees;
			// 3. Clear the large tree ltree_i by setting its size to 0.
			if (st_sizes[i] < size_thr)
				shrink_tree(i, st_roots[i], ltrees[i]);

			bool isLarge = false;
			int stree_idx = trees_pointers[i];
			if (stree_idx == -1)
				isLarge = true;

			small_tree& stree_i = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;

			// Calculate the subtree size of each nodes.
			//	unordered_map<int, int>& cover_i = cover[i];
			//vector<int> cover_i(numOfVertices, 0);

			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];

				++cover_tid[q]; // Increment because it hits the r-q path.
				++cnt16[q][i%CNT];
				//++sum[q];

				int parent_node;

				if (isLarge)
					parent_node = get_large_parent(q, ltree_i);
				else
					parent_node = get_small_parent(q, stree_i);

				if (parent_node > -1) {
					cover_tid[parent_node] += cover_tid[q];
					cnt16[parent_node][i%CNT] += cover_tid[q];
					//	sum[parent_node] += cover_i[q];
				}
			}
			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];
				cover_tid[q] = 0;
			}

			current_universe += descendants.size();
		}


		//cnt16_aggreate(cnt16, has_picked, cover);

		vector<cover_value_type> picked_value;
		benchmark::heap<2, cover_value_type, NodeID> value_heap(numOfVertices);
		//bool switch_small = false;
		switch_small = false;

		//cnt16_2max(picked_value, sum, has_picked,affected, cnt16, switch_small, value_heap);
		//		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, cnt16, switch_small, value_heap);
		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);

		init_time = GetCurrentTimeSec() - init_time;

		cout << "Finish initial sampling " << init_time << "secs." << endl;

		//for (NodeID v = 0; v < V; ++v)	value_heap.update(v, -picked_value[v]); // Min-heap.

		double max_value = -1;

		// Perform V iterations to rank all vertices.
		NodeID ranking = 0;

		NodeID avg_size = current_universe / k;
		//double original_c_paras = c_paras;
		//c_paras = 99999 * c_paras;



	//	NodeID avg_size = current_universe / k;
	//	c_paras = c_paras * (double)numOfVertices / (double)avg_size;
		double original_beta = beta;
		double beta_lower = 0.1;
		double beta_upper = 1.25;
		//double firstit_start_time = GetCurrentTimeSec();
		for (; ranking < numOfVertices; ++ranking) {


			//double current_time = GetCurrentTimeSec() - firstit_start_time;
			if (ranking % (numOfVertices / 10) == 0)
				cout << "picking the " << ranking << "-th orders." << endl;

			double time_this_iteration = GetCurrentTimeSec();

			NodeID pv = 0;
			cover_value_type pv_cover = -1;
			cover_value_type second_pv_cover = -1;
			NodeID spv = 0;

			/*
			if (ranking == stopv) {
				stop_flag = true;
				cout << "topv labeling time:" << labeling_time << endl;
				break;
			}*/

			double selecting_acc_time = GetCurrentTimeSec();
			// In the beginning, pick the best cover vertex by a linear scan.
			if (switch_small == false) {
				double u1_acc_time = GetCurrentTimeSec();
				//while (true) {
				cover_value_type sum_pv = sum[pv];
				cover_value_type sum_spv = sum[spv];
				for (NodeID v = 0; v < numOfVertices; ++v) {
					if (has_picked[v] == true) continue;

					if (second_pv_cover < picked_value[v]) {
						second_pv_cover = picked_value[v];
						spv = v;
						if (pv_cover < second_pv_cover) {
							cover_value_type tmp = second_pv_cover;
							second_pv_cover = pv_cover;
							pv_cover = tmp;
							NodeID tmpv = spv;
							spv = pv;
							pv = tmpv;
						}
						else if (pv_cover == second_pv_cover) {
							if (sum_pv < sum_spv) {
								cover_value_type tmp = second_pv_cover;
								second_pv_cover = pv_cover;
								pv_cover = tmp;
								NodeID tmpv = spv;
								spv = pv;
								pv = tmpv;
							}
						}
					}
					else if (second_pv_cover == picked_value[v]) {
						if (sum_spv < sum[v]) {
							second_pv_cover = picked_value[v];
							spv = v;
							if (sum_pv < sum_spv) {
								cover_value_type tmp = second_pv_cover;
								second_pv_cover = pv_cover;
								pv_cover = tmp;
								NodeID tmpv = spv;
								spv = pv;
								pv = tmpv;
							}
						}
					}
				}

				if (sum[pv] < numOfVertices / 8 && switch_small == false) {
					//cout << "u3 time:" << u3_time << endl;
					cout << ranking << endl;
					cout << pv_cover << endl;
					cout << sum[pv] << endl;
					cout << picked_count << endl;
					cout << "Switch to max-heap representations." << endl;
					switch_small = true;
				}
				u1_time += (GetCurrentTimeSec() - u1_acc_time);

			}
			else {// Later, utilizing the max-heap directly
				double u2_acc_time = GetCurrentTimeSec();

				value_heap.extract_min(pv, pv_cover);
				pv_cover = -pv_cover;

				if (ranking == numOfVertices - 1) continue;
				second_pv_cover = -value_heap.top();
				spv = value_heap.top_value();

				u2_time += (GetCurrentTimeSec() - u2_acc_time);

			}

			selecting_time += (GetCurrentTimeSec() - selecting_acc_time);

			//// FIX
			//if (switch_small == true && reverse_index.size() > 0 && sum[pv] <= reverse_index[pv].size() * 10) {//sum[pv] <= in_tree_count[pv] * 10){//
			//	cout << "I think we can stop here " << GetCurrentTimeSec() - start_time << " ranking:" << ranking << endl;
			//	break;
			//}

			//if (has_picked[pv] == true)
			//	cout << "wrong1" << endl;


			has_picked[pv] = true;
			double labeling_acc_time = GetCurrentTimeSec();
			labeled_arcs = labeling_source_dij(pv, pqueue, visited_que, distances, vis, dst_r, usd, ranking, tmp_idx, r_tmp_idx, wgraph);
			labeling_time += (GetCurrentTimeSec() - labeling_acc_time);
			inv[ranking] = pv;
			rank[pv] = ranking;

			top_degree.insert(d_order[ranking]);
			top_bo.insert(pv);
			top_union.insert(d_order[ranking]);
			top_union.insert(pv);

			if (ranking < 100) {
				top_bo.insert(pv);
				top_union.insert(d_order[ranking]);
				top_union.insert(pv);

				if (top_degree.find(pv) != top_degree.end()) {
					intersection_count++;
				}

				/*
				if (ranking != pv) {
				if (top_degree.find(pv) != top_degree.end())
				intersection_count++;

				if (top_bo.find(d_order[ranking]) != top_bo.end())
				intersection_count++;
				}
				else
				intersection_count++;

				if (intersection_count == 0 || top_union.size() == 0)
				refactor_by_jaccard = 1000;
				else
				refactor_by_jaccard = (double)top_union.size() / (double)intersection_count;
				*/
			}


			/*if (ranking != pv) {
				if (top_degree.find(pv) != top_degree.end())
					intersection_count++;

				if (top_bo.find(d_order[ranking]) != top_bo.end())
					intersection_count++;
			}
			else
				intersection_count++;*/
			if (ranking == 100) {
				if (intersection_count == 0 || top_union.size() == 0)
					refactor_by_jaccard = 0;
				else
					refactor_by_jaccard = (double)intersection_count / (double)top_union.size();
			}
			//cout << (double)top_union.size() / (double)intersection_count << endl;

			//cout << ranking << "\t" << current_time << "\t" << labeled_arcs << "\t" << pv_cover << "\t" << sum[pv] << "\t" << second_pv_cover << "\t" << sum[spv] << endl;
			//if (ranking < k * 10) {
			//	if (DIRECTED_FLAG == false)
			//		covered_universe += labeled_arcs;
			//	else
			//		covered_universe += labeled_arcs / 2;
			//}

			//if (ranking == numOfVertices * 0.1) {
			//	//c_paras = original_c_paras * (double)numOfVertices / (covered_universe / k);
			//	//c_paras = 10000;
			//	cout << avg_size << "\t" << covered_universe << "\t" << current_universe << endl;
			//}

			//cout << (double)labeled_arcs / (double)numOfVertices << endl;

			if (ranking == numOfVertices - 1) {
				ranking++;
				break;
			}


			int current_trees_num = st_roots.size();

			// 4 steps in updating:
			double updating_acc_time = GetCurrentTimeSec();
			if (switch_small == false) { // We still iterate all alive trees to calculate the cover.	
										 //	double u3_acc_time = GetCurrentTimeSec();
				for (int tid = 0; tid < current_trees_num; ++tid) {

					// Get the descendants and ascendants of pv in tree tid.
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					// pv has already been covered in tree tid.
					int pvpid = 0;
					if (isLarge)
						pvpid = get_large_parent(pv, ltree_tid);
					else {
						//double get_small_parent_incre_time = GetCurrentTimeSec();
						pvpid = get_small_parent(pv, stree_tid);
						//	get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}
					if (pvpid == -1 && st_roots[tid] != pv) // pv has been covered in tree tid 
						continue;


					double u3_acc_time = GetCurrentTimeSec();

					// //1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;				 
					//get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph.adj, isLarge, stree_tid, ltree_tid);
					//get_ascendants(tid, pv, ascendants, que, vis, wgraph.adj, wgraph.r_adj, isLarge, stree_tid, ltree_tid);

					//Array Representation
					get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph, isLarge, stree_tid, ltree_tid);
					get_ascendants(tid, pv, ascendants, que, vis, wgraph, isLarge, stree_tid, ltree_tid);

					u3_time += (GetCurrentTimeSec() - u3_acc_time);


					double u4_acc_time = GetCurrentTimeSec();
					NodeID subtract_size = descendants.size();
					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
												//	cover[tid][av] -= subtract_size;
						cnt16[av][tid%CNT] -= subtract_size;
						//sum[av] -= subtract_size;
						//affected[av] = true;

						if (affected_list[av] == false) {
							affected_list[av] = true;
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
						}
					}
					u4_time += (GetCurrentTimeSec() - u4_acc_time);


					double u5_acc_time = GetCurrentTimeSec();
					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					//vector<int> cover_tid(numOfVertices, 0);
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//sum[dv] -= cover_tid[dv];

						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];
						if (cnt16[dv][tid%CNT] < 0)
							cout << "a6 " << cnt16[dv][tid%CNT] << endl;

						/*if (isLarge)
						remove_large_node(dv, ltree_tid);
						else
						remove_small_node(dv, stree_tid);*/

						//						remove_node(tid, dv);
						//affected[dv] = true;

						if (affected_list[dv] == false) {
							affected_count++;
							affected_list[dv] = true;
							if (switch_set)
								affected_v.insert(dv);
						}
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;

						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else
							remove_small_node(dv, stree_tid);
					}
					
					current_universe -= descendants.size();
					
					u5_time += (GetCurrentTimeSec() - u5_acc_time);

					// Drop the size of subtree.
					st_sizes[tid] -= subtract_size;
					if (trees_pointers[tid] != -1) alive_resources -= (long long)subtract_size;

					// 3. Convert shrink trees to small trees representataions.	
					if (st_sizes[tid] > 0) {
						if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
							shrink_tree(tid, st_roots[tid], ltrees[tid]);
						}
					}

					if (st_roots[tid] == pv)
						if (isLarge) {
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
							//ltrees.shrink_to_fit();
						}
						else {
							stree_tid.clear();
							stree_tid.resize(0);
							//strees.shrink_to_fit();
						}


				}
				//	u3_time += (GetCurrentTimeSec() - u3_acc_time);
			}
			else { // We use the reverse index to update trees.
				   // Build the reverse_index first (one time).
				if (reverse_index.size() == 0) {
					reverse_index.resize(numOfVertices);
					for (NodeID v = 0; v < numOfVertices; ++v) {
						reverse_index[v].set_empty_key(numOfVertices + 1);
						reverse_index[v].set_deleted_key(numOfVertices + 2);
					}
					for (int tid = 0; tid < current_trees_num; ++tid) {
						if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue; //This whole tree has already been covered. 

						bool isLarge = false;
						int stree_idx = trees_pointers[tid];
						if (stree_idx == -1)
							isLarge = true;

						small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
						large_tree& ltree_tid = ltrees[tid];

						for (int v = 0; v < numOfVertices; ++v) {

							int vpid = 0;
							if (isLarge)
								vpid = get_large_parent(v, ltree_tid);
							else {
								//double get_small_parent_incre_time = GetCurrentTimeSec();
								vpid = get_small_parent(v, stree_tid);
								//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
							}
							if ((vpid != -1 || st_roots[tid] == v)) {
								reverse_index[v].insert(tid);
								//in_tree_count[v]++;
							}
						}
					}
				}

				//cout << reverse_index[pv].size() << "," << sum[pv] << "," <<  (double)sum[pv] / (double)reverse_index[pv].size() << endl;
				// Access the alive trees directly by the reverse index.
				//for (unordered_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) {
				for (google::dense_hash_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) {
					int tid = *it;
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					double u6_acc_time = GetCurrentTimeSec();
					/*
					int parent_node;
					if (isLarge)
					parent_node = get_large_parent(pv, ltree_tid);
					else {
					//double get_small_parent_incre_time = GetCurrentTimeSec();
					parent_node = get_small_parent(pv, stree_tid);
					//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}
					if(parent_node == -1 && pv != st_roots[tid])
					continue;
					*/
					double u7_acc_time = GetCurrentTimeSec();

					//	// 1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;
					//	//double u4_acc_time = getcurrenttimesec();
					//	get_descendants(tid, pv, descendants, descendants_parents,que, vis, wgraph.adj, islarge, stree_tid, ltree_tid);
					////	u4_time += (getcurrenttimesec() - u4_acc_time);

					//	get_ascendants(tid, pv, ascendants, que, vis, wgraph.adj, wgraph.r_adj, isLarge, stree_tid, ltree_tid);

					// Array Representation
					get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph, isLarge, stree_tid, ltree_tid);
					get_ascendants(tid, pv, ascendants, que, vis, wgraph, isLarge, stree_tid, ltree_tid);

					NodeID subtract_size = descendants.size();
					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
												//cover[tid][av] -= subtract_size;
						cnt16[av][tid%CNT] -= subtract_size;
						//affected[av] = true;
						//	sum[av] -= subtract_size;

						if (affected_list[av] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
							affected_list[av] = true;
						}
					}
					u7_time += (GetCurrentTimeSec() - u7_acc_time);

					//double u2_acc_time = GetCurrentTimeSec();
					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					double u8_acc_time = GetCurrentTimeSec();
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//	sum[dv] -= cover_tid[dv];

						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];
						//	if (cnt16[dv][tid%CNT] < 0)
						//		cout << "a6 " << cnt16[dv][tid%CNT] << endl;

						//if(dv == pv){
						/*if (isLarge)
						remove_large_node(dv, ltree_tid);
						else
						remove_small_node(dv, stree_tid);*/
						//}
						//affected[dv] = true;

						//double u82_acc_time = GetCurrentTimeSec();
						if (dv != pv)
							reverse_index[dv].erase(tid);
						//if (dv != pv)
						//	in_tree_count[dv]--;

						//	double u81_acc_time = GetCurrentTimeSec();
						if (affected_list[dv] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(dv);
							affected_list[dv] = true;
						}

						//u81_time += (GetCurrentTimeSec() - u81_acc_time);

						//	u82_time += (GetCurrentTimeSec() - u82_acc_time);
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;
						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else
							remove_small_node(dv, stree_tid);
					}
					u8_time += (GetCurrentTimeSec() - u8_acc_time);


					current_universe -= descendants.size();

					if (st_roots[tid] == pv)
						if (isLarge) {
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
						}
						else {
							stree_tid.clear();
							stree_tid.resize(0);
						}


						st_sizes[tid] -= subtract_size;
						if (trees_pointers[tid] != -1) alive_resources -= (long long)subtract_size;
						//	if (st_sizes[tid] < 0) {
						//		cout << "a5" << st_sizes[tid] << endl;
						//	}

						// 3. Convert shrink trees to small trees representataions.				
						if (st_sizes[tid] > 0) {
							if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
								shrink_tree(tid, st_roots[tid], ltrees[tid]);
							}
						}

						//u2_time += (GetCurrentTimeSec() - u2_acc_time);

				}
				reverse_index[pv].clear();
				reverse_index[pv].resize(0);

			}
			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			// 4. Sample new trees based on.
			// If we sample a new tree, we need to dynamically increment the following data structures:
			// 1) ltrees
			// 2) trees_pointers (-1 or pointer to strees)
			// 3) cover
			// 4) st_roots
			// 5) st_sizes 
			// 6) strees (maybe)
			// two constraints: cl > ct and 10kn total vertices
			double adding_acc_time = GetCurrentTimeSec();


			/*cout << "ranking:" << ranking << endl;
			cout << "pick cover:" << pv_cover << endl;
			cout << "pick count:" << picked_count << endl;
			cout << "labeled_arcs:" << labeled_arcs << endl;
			cout << "new_tree_arcs:" << new_tree_arcs << endl;
			cout << "alive_resources:" << alive_resources << endl;
			cout << "bounded_resources:" << bounded_resources << endl;
			cout << "#tree:" << cover.size() << endl;
			cout << "u5/adding" << u5_time / adding_time << endl;
			*/
			new_tree_arcs = 0;
			int how_many = 0;
			

			cover_value_type sample_constraint = 10000000 * current_universe / st_roots.size() * refactor_by_jaccard;
	//		cout << refactor_by_jaccard << "\t" << pv << "\t" << d_order[ranking] <<  endl;
			cover_value_type sum_spv = sum[spv];

			if (ranking == 100) {
				if (refactor_by_jaccard < 0.1)
					beta = 1;
				else if (refactor_by_jaccard >= 0.1 && refactor_by_jaccard < 0.2)
					beta = 0.75;
				else if (refactor_by_jaccard >= 0.2 && refactor_by_jaccard < 0.5)
					beta = 0.5;
				else if (refactor_by_jaccard >= 0.5)
					beta = 0.25;
			}

			//if (ranking <= 0.1 * numOfVertices) {
			//	
			///*	beta = original_beta * refactor_by_jaccard / c_paras;
			//	if (beta > original_beta)
			//		beta = original_beta;*/


			//	double tratio = (double)(1 - (double)1 / (double)refactor_by_jaccard);
			//	beta = (beta_upper - beta_lower) * tratio + beta_lower;
			//		
			//}  
			//else
			//	beta = original_beta;
			/*if (ranking <= 0.1 * numOfVertices) {
				cout << refactor_by_jaccard << "\t" << beta << endl;
			}*/
		/*	if (ranking <= 0.1 * numOfVertices) {
				cout << refactor_by_jaccard << "\t" << beta << endl;
			}*/

			while (beta * labeled_arcs > new_tree_arcs && alive_resources < bounded_resources) { // && sum_spv < sample_constraint) {

				if (picked_count == numOfVertices) break;

				NodeID r = rand() % numOfVertices;
				if (picked_as_root[r] == true) continue;
				if (trivial_root[r] == true) continue;

				large_tree new_ltree;
				small_tree new_stree;
				//GOOGLE DENSE MAP
				new_stree.set_empty_key(numOfVertices + 1);
				new_stree.set_deleted_key(numOfVertices + 2);
				bool isLarge = false;
				int tmp_tree_arcs = 0;

				if (!BUILD_SMALL_TREE_FLAG) {

					//double u5_acc_time = GetCurrentTimeSec();
					double u9_acc_time = GetCurrentTimeSec();
					tmp_tree_arcs = build_large_tree_dij(r, new_ltree, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, new_tree_arcs, beta * labeled_arcs);
					u9_time += (GetCurrentTimeSec() - u9_acc_time);
					//u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;

					if (descendants.size() < size_thr)
						BUILD_SMALL_TREE_FLAG = true;

					if (descendants.size() < 100) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					//	double u91_acc_time = GetCurrentTimeSec();

					isLarge = true;
					how_many++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					used_as_cover++;
					st_roots.push_back(r);
					trees_pointers.push_back(-1);
					st_sizes.push_back(descendants.size());

					//new_tree_arcs += descendants.size();
					alive_resources += (long long)numOfVertices;

					//	u91_time += (GetCurrentTimeSec() - u91_acc_time);
					// Turn large tree to small tree.
					// 1. Build new small tree with ltree_i and add it to strees;
					// 2. Set tree_pointers[i] to the position of the strees;
					// 3. Clear the large tree ltree_i by setting its size to 0.

					//	double u92_acc_time = GetCurrentTimeSec();
					int i = st_roots.size() - 1;
					if (st_sizes[i] < size_thr) {
						shrink_tree(i, st_roots[i], new_ltree);
						// FIX
						isLarge = false;
						/*for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID q = descendants[di];
						NodeID v_parent = get_large_parent(q, new_ltree);
						NodeID s_parent = get_small_parent(q, strees[strees.size() - 1]);
						if (v_parent != s_parent) cout << v_parent << " vs " << s_parent << endl;
						}*/
					}
					ltrees.push_back(new_ltree);
					//	u92_time += (GetCurrentTimeSec() - u92_acc_time);
				}
				else {//Build Small Tree Directly
					  //double u5_acc_time = GetCurrentTimeSec();
					double u9_acc_time = GetCurrentTimeSec();
					tmp_tree_arcs = build_small_tree_dij(r, new_stree, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, new_tree_arcs, beta * labeled_arcs);
					//u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;
					u9_time += (GetCurrentTimeSec() - u9_acc_time);

					if (descendants.size() >= size_thr)
						BUILD_SMALL_TREE_FLAG = false;

					if (descendants.size() < 100) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					isLarge = false;

					how_many++;
					used_as_cover++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					st_roots.push_back(r);
					trees_pointers.push_back(strees.size()); //we assume it is inserted.
					st_sizes.push_back(descendants.size());

					alive_resources += (long long)new_stree.size();
					//strees.push_back(new_stree);
					//// new_ltree is empty.
					//ltrees.push_back(new_ltree);

					int i = st_roots.size() - 1;
					if (st_sizes[i] >= size_thr) {
						enlarge_tree(i, st_roots[i], new_stree, new_ltree);
						new_stree.clear();
						new_stree.resize(0);
						isLarge = true;
					}
					else {
						strees.push_back(new_stree);
						ltrees.push_back(new_ltree);
					}

				}

				int i = st_roots.size() - 1;

				int numOfStrees = strees.size();
				small_tree &new_stree_ref = numOfStrees != 0 ? strees[strees.size() - 1] : empty_small_tree;

				// Calculate the subtree size of each nodes.

				double u10_acc_time = GetCurrentTimeSec();
				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];

					++cover_tid[q]; // Increment because it hitsthe r-q path.
					++cnt16[q][i%CNT];
					//++sum[q];
					//affected[q] = true;

					if (affected_list[q] == false) {
						affected_count++;
						if (switch_set)
							affected_v.insert(q);
						affected_list[q] = true;
					}

					if (switch_small == true) {
						reverse_index[q].insert(i);
						//in_tree_count[q]++;
					}


					int parent_node;
					if (isLarge)
						parent_node = get_large_parent(q, new_ltree);
					else {
						//double get_small_parent_incre_time = GetCurrentTimeSec();
						parent_node = get_small_parent(q, new_stree_ref);
						//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}
					if (parent_node > -1) {
						cover_tid[parent_node] += cover_tid[q];
						cnt16[parent_node][i%CNT] += cover_tid[q];
						//sum[parent_node] += cover_new[q];
					}
				}
				u10_time += (GetCurrentTimeSec() - u10_acc_time);
				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];
					cover_tid[q] = 0;
				}


				current_universe += descendants.size(); 
				sample_constraint = c_paras * current_universe / st_roots.size();


				//cout << "beta*labeld_arcs: " << beta * labeled_arcs << " new tree _arcs: " << new_tree_arcs << " alive " << alive_resources << " bounded " << bounded_resources << endl;

			}
			adding_time += (GetCurrentTimeSec() - adding_acc_time);


			//double u3_acc_time = GetCurrentTimeSec();

			updating_acc_time = GetCurrentTimeSec();

			// Build the cover value heap for one time.
			if (switch_small == true && value_heap.empty())
				for (int nid = 0; nid < numOfVertices; ++nid)
					if (has_picked[nid] == false)
						value_heap.update(nid, -picked_value[nid]);

			//	cnt16_2max(picked_value, sum, has_picked, affected, cnt16, switch_small, value_heap);
			//c	cnt16_2max_queue(picked_value, sum, has_picked, affected_v, cnt16, switch_small, value_heap);

			double u11_acc_time = GetCurrentTimeSec();
			cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);
			u11_time += (GetCurrentTimeSec() - u11_acc_time);
			//	iteration_time[ranking] = (GetCurrentTimeSec() - updating_acc_time);
			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			if (affected_count < numOfVertices / 8)
				switch_set = true;
			else
				switch_set = false;
			affected_count = 0;

			time_this_iteration = GetCurrentTimeSec() - time_this_iteration;
			iteration_time[ranking] = time_this_iteration;

		}

		cout << "avg size:" << avg_size << endl;
		cout << "picked as tree covers:" << used_as_cover << endl;
		cout << "picked count:" << picked_count << endl;

		//// FIX
		//for (int v = 0; v < numOfVertices; ++v){
		//	for (int i = 0; i < CNT; ++i) {
		//		if (cnt16[v][i] != 0)
		//			cout << "mi?" << endl;
		//	}
		//}
		if (stop_flag == true) return;

		double labeling_acc_time = GetCurrentTimeSec();
		// Order the remaining vertices with degree, if any.
		if (ranking != numOfVertices) {
			vector<NodeID> degree(numOfVertices);
			vector<NodeID> r_degree(numOfVertices);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				/*vector<NodeID>& adj_v = graph.adj[v];
				for (NodeID i = 0; i < adj_v.size(); ++i) {
				if (has_picked[adj_v[i]] != true)
				degree[v]++;
				}*/

				for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid) {
					NodeID w = wgraph.edges[eid].first;
					if (has_picked[w] != true)
						degree[v]++;
				}

				if (DIRECTED_FLAG == true) {
					/*vector<NodeID>& r_adj_v = graph.r_adj[v];
					for (NodeID i = 0; i < r_adj_v.size(); ++i) {
					if (has_picked[r_adj_v[i]] != true)
					r_degree[v]++;
					}*/
					for (EdgeID eid = wgraph.r_vertices[v]; eid < wgraph.r_vertices[v + 1]; ++eid) {
						NodeID w = wgraph.r_edges[eid].first;
						if (has_picked[w] != true)
							r_degree[v]++;
					}
				}
			}

			vector<pair<NodeID, NodeID> > degree_v;
			degree_v.reserve(numOfVertices - ranking);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				if (DIRECTED_FLAG == true)
					degree_v.push_back(make_pair(degree[v] * r_degree[v], v));
				else
					degree_v.push_back(make_pair(degree[v], v));
			}
			sort(degree_v.rbegin(), degree_v.rend());

			int count = 0;
			for (; ranking < numOfVertices; ++ranking) {
				NodeID tpv = degree_v[count++].second;
				inv[ranking] = tpv;
				rank[tpv] = ranking;
				labeling_source_dij(tpv, pqueue, visited_que, distances, vis, dst_r, usd, ranking, tmp_idx, r_tmp_idx, wgraph);
			}
		}
		labeling_time += (GetCurrentTimeSec() - labeling_acc_time);


		// Finish the constructions.
		for (NodeID v = 0; v < numOfVertices; ++v) {
			if (DIRECTED_FLAG == false) {
				NodeID k = tmp_idx[v].first.size();
				labels.index_[v].spt_v.resize(k);
				labels.index_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();
				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();

			}
			else {
				NodeID k = tmp_idx[v].first.size();
				dlabels.index_[v].spt_v.resize(k);
				dlabels.index_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();
				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();

				k = r_tmp_idx[v].first.size();
				dlabels.bindex_[v].spt_v.resize(k);
				dlabels.bindex_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_v[i] = r_tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_d[i] = r_tmp_idx[v].second[i];
				r_tmp_idx[v].first.clear();
				r_tmp_idx[v].second.clear();
				r_tmp_idx[v].first.shrink_to_fit();
				r_tmp_idx[v].second.shrink_to_fit();
			}
		}

		num_of_trees = st_roots.size();

		double small_tree_num = 0;

		for (int i = 0; i < num_of_trees; ++i) {
			if (trees_pointers[i] != -1) {
				small_tree_num++;
				//small_tree_num += strees[trees_pointers[i]].size();
			}
		}

		total_resources = small_tree_num;

	}

	Betweenness_Ordering(NodeID k, double beta, WGraph& wgraph, NodeID stopv, long long bound_times, long long trans_times, double mem_bound, double c_paras, double es_paras) {

		double start_time = GetCurrentTimeSec();

		if (k < CNT) {
			cout << " Should pick at least 16 seed!" << endl;
			return;
		}
		iteration_time.resize(numOfVertices);

		updating_time = 0;
		labeling_time = 0;
		adding_time = 0;
		selecting_time = 0;
		init_time = 0;
		u1_time = 0;
		u2_time = 0;
		u3_time = 0;
		num_of_trees = 0;
		total_resources = 0;

		inv.resize(numOfVertices, -1);
		rank.resize(numOfVertices, -1);

		if (DIRECTED_FLAG == false) {
			labels.index_.resize(numOfVertices);
		}
		else { //if (DIRECTED_FLAG == true) {
			dlabels.index_.resize(numOfVertices);
			dlabels.bindex_.resize(numOfVertices);
		}

		// Resource Bounding.
		NodeID new_tree_arcs = 0;
		NodeID labeled_arcs = 0;
		alive_resources = 0;
		bounded_resources = bound_times * (long long)k * (long long)numOfVertices;

		vector<NodeID> st_roots; // The root vertices of the sampled trees.
		vector<NodeID> st_sizes; // The size of the sampled trees.

		NodeID size_thr = numOfVertices < 8 ? numOfVertices : numOfVertices / trans_times; // When a tree shrink to a size smaller than the size_thr, it should be converted into a small_tree representation.

		init_time = GetCurrentTimeSec();
		// Building the initial k trees.
		//srand((unsigned)time(NULL));

		srand(0x728f9d0);

		// Reserve pointers of vectors.
		ltrees.reserve(numOfVertices);
		strees.reserve(numOfVertices);
		trees_pointers.reserve(numOfVertices);
		st_roots.reserve(numOfVertices);
		st_sizes.reserve(numOfVertices);

		st_roots.resize(k);
		st_sizes.resize(k);

		ltrees.resize(k);
		strees.clear();
		trees_pointers.resize(k, -1); // Point the sampling trees to its tree representations. large tree [0, V) small tree [V, 2V).
									  //vector<unordered_map<int, int> > cover(k, unordered_map<int, int>());

		vector<NodeID> ascendants;
		ascendants.reserve(numOfVertices);
		vector<NodeID> descendants;
		descendants.reserve(numOfVertices);
		vector<NodeID> que(numOfVertices);
		vector<bool> vis(numOfVertices);
		vector<bool> has_picked(numOfVertices);
		vector<vector<cover_value_type> > cnt16(numOfVertices, vector<cover_value_type>(CNT));
		//		vector<bool> affected(V, true);
		vector<bool> picked_as_root(numOfVertices, false);
		vector<bool> trivial_root(numOfVertices, false);
		int picked_count = 0;
		int used_as_cover = 0;
		vector<cover_value_type> sum(numOfVertices);
		bool stop_flag = false;


		// tmp structure used in the dijkstra searches. 

		// Binary Heap
		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);

		queue<NodeID> visited_que;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
		vector<bool> usd(numOfVertices, false);

		//GOOGLE
		//		unordered_set<NodeID> affected_v;
		google::dense_hash_set<NodeID> affected_v;
		affected_v.set_empty_key(numOfVertices + 1);
		affected_v.set_deleted_key(numOfVertices + 2);

		vector<cover_value_type> cover_tid(numOfVertices, 0);
		//vector<cover_value_type> cover_new(numOfVertices, 0); // modified-5-14


		vector<bool> affected_list(numOfVertices, false);
		int affected_count = 0;
		bool switch_set = false;

		//vector<unordered_set<int> > reverse_index(0);
		//GOOGLE DENSE SET
		vector<google::dense_hash_set<NodeID> > reverse_index;
		reverse_index.reserve(numOfVertices);
		//vector<NodeID> in_tree_count(numOfVertices);

		// Preparing basic structure for pl algorithm.
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT))); // Backward labels.

		vector<NodeID> descendants_parents(numOfVertices, -1);


		cover_value_type current_universe = 0;

		// Building the initial k sampled trees. 
		for (int i = 0; i < k; ++i) {
			NodeID r = rand() % numOfVertices;
			if (picked_as_root[r] == true) { --i; continue; }

			st_roots[i] = r;
			large_tree &ltree_i = ltrees[i];
			picked_as_root[r] = true;
			picked_count++;
			used_as_cover++;

			build_large_tree_dij(r, ltree_i, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, 0, numOfVertices);
			st_sizes[i] = descendants.size();
			alive_resources += (long long)numOfVertices;

			// Turn large tree to small tree.
			// 1. Build new small tree with ltree_i and add it to strees;
			// 2. Set tree_pointers[i] to the position of the strees;
			// 3. Clear the large tree ltree_i by setting its size to 0.
			if (st_sizes[i] < size_thr)
				shrink_tree(i, st_roots[i], ltrees[i]);

			bool isLarge = false;
			int stree_idx = trees_pointers[i];
			if (stree_idx == -1)
				isLarge = true;

			small_tree& stree_i = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;

			// Calculate the subtree size of each nodes.
			//	unordered_map<int, int>& cover_i = cover[i];
			//vector<int> cover_i(numOfVertices, 0);

			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];

				++cover_tid[q]; // Increment because it hits the r-q path.
				++cnt16[q][i%CNT];
				//++sum[q];

				int parent_node;

				if (isLarge)
					parent_node = get_large_parent(q, ltree_i);
				else
					parent_node = get_small_parent(q, stree_i);

				if (parent_node > -1) {
					cover_tid[parent_node] += cover_tid[q];
					cnt16[parent_node][i%CNT] += cover_tid[q];
					//	sum[parent_node] += cover_i[q];
				}
			}
			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];
				cover_tid[q] = 0;
			}

			current_universe += descendants.size();
		}


		//cnt16_aggreate(cnt16, has_picked, cover);

		vector<cover_value_type> picked_value;
		benchmark::heap<2, cover_value_type, NodeID> value_heap(numOfVertices);
		//bool switch_small = false;
		switch_small = false;

		//cnt16_2max(picked_value, sum, has_picked,affected, cnt16, switch_small, value_heap);
		//		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, cnt16, switch_small, value_heap);
		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);

		init_time = GetCurrentTimeSec() - init_time;

		cout << "Finish initial sampling " << init_time << "secs." << endl;

		//for (NodeID v = 0; v < V; ++v)	value_heap.update(v, -picked_value[v]); // Min-heap.

		double max_value = -1;

		// Perform V iterations to rank all vertices.
		NodeID ranking = 0;


		NodeID avg_size = current_universe / k;
		c_paras = c_paras * (double)numOfVertices / (double)avg_size;

		for (; ranking < numOfVertices; ++ranking) {
			if (ranking % (numOfVertices / 10) == 0)
				cout << "picking the " << ranking << "-th orders." << endl;

			double time_this_iteration = GetCurrentTimeSec();

			NodeID pv = 0;
			cover_value_type pv_cover = -1;
			cover_value_type second_pv_cover = -1;
			NodeID spv = 0;

			if (ranking == stopv) {
				stop_flag = true;
				cout << "topv labeling time:" << labeling_time << endl;
				break;
			}

			double selecting_acc_time = GetCurrentTimeSec();
			// In the beginning, pick the best cover vertex by a linear scan.
			if (switch_small == false) {
				double u1_acc_time = GetCurrentTimeSec();
				//while (true) {
				cover_value_type sum_pv = sum[pv];
				cover_value_type sum_spv = sum[spv];
				for (NodeID v = 0; v < numOfVertices; ++v) {
					if (has_picked[v] == true) continue;

					if (second_pv_cover < picked_value[v]) {
						second_pv_cover = picked_value[v];
						spv = v;
						if (pv_cover < second_pv_cover) {
							cover_value_type tmp = second_pv_cover;
							second_pv_cover = pv_cover;
							pv_cover = tmp;
							NodeID tmpv = spv;
							spv = pv;
							pv = tmpv;
						}
						else if (pv_cover == second_pv_cover) {
							if (sum_pv < sum_spv) {
								cover_value_type tmp = second_pv_cover;
								second_pv_cover = pv_cover;
								pv_cover = tmp;
								NodeID tmpv = spv;
								spv = pv;
								pv = tmpv;
							}
						}
					}
					else if (second_pv_cover == picked_value[v]) {
						if (sum_spv < sum[v]) {
							second_pv_cover = picked_value[v];
							spv = v;
							if (sum_pv < sum_spv) {
								cover_value_type tmp = second_pv_cover;
								second_pv_cover = pv_cover;
								pv_cover = tmp;
								NodeID tmpv = spv;
								spv = pv;
								pv = tmpv;
							}
						}
					}
				}

				if (sum[pv] < numOfVertices / 8 && switch_small == false) {
					//cout << "u3 time:" << u3_time << endl;
					cout << ranking << endl;
					cout << pv_cover << endl;
					cout << sum[pv] << endl;
					cout << picked_count << endl;
					cout << "Switch to max-heap representations." << endl;
					switch_small = true;
				}
				u1_time += (GetCurrentTimeSec() - u1_acc_time);

			}
			else {// Later, utilizing the max-heap directly
				double u2_acc_time = GetCurrentTimeSec();

				value_heap.extract_min(pv, pv_cover);
				pv_cover = -pv_cover;

				if (ranking == numOfVertices - 1) continue;
				second_pv_cover = -value_heap.top();
				spv = value_heap.top_value();

				u2_time += (GetCurrentTimeSec() - u2_acc_time);

			}

			selecting_time += (GetCurrentTimeSec() - selecting_acc_time);

			//// FIX
			if (switch_small == true && reverse_index.size() > 0 && sum[pv] <= reverse_index[pv].size() * es_paras) {//sum[pv] <= in_tree_count[pv] * 10){//
				cout << "I think we can stop here " << GetCurrentTimeSec() - start_time << " ranking:" << ranking << endl;
				break;
			}

			//if (has_picked[pv] == true)
			//	cout << "wrong1" << endl;

			has_picked[pv] = true;
			double labeling_acc_time = GetCurrentTimeSec();
			labeled_arcs = labeling_source_dij(pv, pqueue, visited_que, distances, vis, dst_r, usd, ranking, tmp_idx, r_tmp_idx, wgraph);
			labeling_time += (GetCurrentTimeSec() - labeling_acc_time);
			inv[ranking] = pv;
			rank[pv] = ranking;

			if (ranking == numOfVertices - 1) {
				ranking++;
				break;
			}


			int current_trees_num = st_roots.size();

			// 4 steps in updating:
			double updating_acc_time = GetCurrentTimeSec();
			if (switch_small == false) { // We still iterate all alive trees to calculate the cover.	
										 //	double u3_acc_time = GetCurrentTimeSec();
				for (int tid = 0; tid < current_trees_num; ++tid) {

					// Get the descendants and ascendants of pv in tree tid.
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					// pv has already been covered in tree tid.
					int pvpid = 0;
					if (isLarge)
						pvpid = get_large_parent(pv, ltree_tid);
					else {
						//double get_small_parent_incre_time = GetCurrentTimeSec();
						pvpid = get_small_parent(pv, stree_tid);
						//	get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}
					if (pvpid == -1 && st_roots[tid] != pv) // pv has been covered in tree tid 
						continue;


					double u3_acc_time = GetCurrentTimeSec();

					// //1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;				 
					//get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph.adj, isLarge, stree_tid, ltree_tid);
					//get_ascendants(tid, pv, ascendants, que, vis, wgraph.adj, wgraph.r_adj, isLarge, stree_tid, ltree_tid);

					//Array Representation
					get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph, isLarge, stree_tid, ltree_tid);
					get_ascendants(tid, pv, ascendants, que, vis, wgraph, isLarge, stree_tid, ltree_tid);

					u3_time += (GetCurrentTimeSec() - u3_acc_time);


					double u4_acc_time = GetCurrentTimeSec();
					NodeID subtract_size = descendants.size();
					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
												//	cover[tid][av] -= subtract_size;
						cnt16[av][tid%CNT] -= subtract_size;
						//sum[av] -= subtract_size;
						//affected[av] = true;

						if (affected_list[av] == false) {
							affected_list[av] = true;
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
						}
					}
					u4_time += (GetCurrentTimeSec() - u4_acc_time);


					double u5_acc_time = GetCurrentTimeSec();
					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					//vector<int> cover_tid(numOfVertices, 0);
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//sum[dv] -= cover_tid[dv];

						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];
						if (cnt16[dv][tid%CNT] < 0)
							cout << "a6 " << cnt16[dv][tid%CNT] << endl;

						/*if (isLarge)
						remove_large_node(dv, ltree_tid);
						else
						remove_small_node(dv, stree_tid);*/

						//						remove_node(tid, dv);
						//affected[dv] = true;

						if (affected_list[dv] == false) {
							affected_count++;
							affected_list[dv] = true;
							if (switch_set)
								affected_v.insert(dv);
						}
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;

						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else
							remove_small_node(dv, stree_tid);
					}

					current_universe -= descendants.size();

					u5_time += (GetCurrentTimeSec() - u5_acc_time);

					// Drop the size of subtree.
					st_sizes[tid] -= subtract_size;
					if (trees_pointers[tid] != -1) alive_resources -= (long long)subtract_size;

					// 3. Convert shrink trees to small trees representataions.	
					if (st_sizes[tid] > 0) {
						if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
							shrink_tree(tid, st_roots[tid], ltrees[tid]);
						}
					}

					if (st_roots[tid] == pv)
						if (isLarge) {
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
							//ltrees.shrink_to_fit();
						}
						else {
							stree_tid.clear();
							stree_tid.resize(0);
							//strees.shrink_to_fit();
						}

				}
				//	u3_time += (GetCurrentTimeSec() - u3_acc_time);
			}
			else { // We use the reverse index to update trees.
				   // Build the reverse_index first (one time).
				if (reverse_index.size() == 0) {
					reverse_index.resize(numOfVertices);
					for (NodeID v = 0; v < numOfVertices; ++v) {
						reverse_index[v].set_empty_key(numOfVertices + 1);
						reverse_index[v].set_deleted_key(numOfVertices + 2);
					}
					for (int tid = 0; tid < current_trees_num; ++tid) {
						if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue; //This whole tree has already been covered. 

						bool isLarge = false;
						int stree_idx = trees_pointers[tid];
						if (stree_idx == -1)
							isLarge = true;

						small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
						large_tree& ltree_tid = ltrees[tid];

						for (int v = 0; v < numOfVertices; ++v) {

							int vpid = 0;
							if (isLarge)
								vpid = get_large_parent(v, ltree_tid);
							else {
								//double get_small_parent_incre_time = GetCurrentTimeSec();
								vpid = get_small_parent(v, stree_tid);
								//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
							}
							if ((vpid != -1 || st_roots[tid] == v)) {
								reverse_index[v].insert(tid);
								//in_tree_count[v]++;
							}
						}
					}
				}

				//cout << reverse_index[pv].size() << "," << sum[pv] << "," <<  (double)sum[pv] / (double)reverse_index[pv].size() << endl;
				// Access the alive trees directly by the reverse index.
				//for (unordered_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) {
				for (google::dense_hash_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) {
					int tid = *it;
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					double u6_acc_time = GetCurrentTimeSec();
					/*
					int parent_node;
					if (isLarge)
					parent_node = get_large_parent(pv, ltree_tid);
					else {
					//double get_small_parent_incre_time = GetCurrentTimeSec();
					parent_node = get_small_parent(pv, stree_tid);
					//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}
					if(parent_node == -1 && pv != st_roots[tid])
					continue;
					*/
					double u7_acc_time = GetCurrentTimeSec();

					//	// 1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;
					//	//double u4_acc_time = getcurrenttimesec();
					//	get_descendants(tid, pv, descendants, descendants_parents,que, vis, wgraph.adj, islarge, stree_tid, ltree_tid);
					////	u4_time += (getcurrenttimesec() - u4_acc_time);

					//	get_ascendants(tid, pv, ascendants, que, vis, wgraph.adj, wgraph.r_adj, isLarge, stree_tid, ltree_tid);

					// Array Representation
					get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph, isLarge, stree_tid, ltree_tid);
					get_ascendants(tid, pv, ascendants, que, vis, wgraph, isLarge, stree_tid, ltree_tid);

					NodeID subtract_size = descendants.size();
					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
												//cover[tid][av] -= subtract_size;
						cnt16[av][tid%CNT] -= subtract_size;
						//affected[av] = true;
						//	sum[av] -= subtract_size;

						if (affected_list[av] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
							affected_list[av] = true;
						}
					}
					u7_time += (GetCurrentTimeSec() - u7_acc_time);

					//double u2_acc_time = GetCurrentTimeSec();
					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					double u8_acc_time = GetCurrentTimeSec();
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//	sum[dv] -= cover_tid[dv];

						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];
						//	if (cnt16[dv][tid%CNT] < 0)
						//		cout << "a6 " << cnt16[dv][tid%CNT] << endl;

						//if(dv == pv){
						/*if (isLarge)
						remove_large_node(dv, ltree_tid);
						else
						remove_small_node(dv, stree_tid);*/
						//}
						//affected[dv] = true;

						//double u82_acc_time = GetCurrentTimeSec();
						if (dv != pv)
							reverse_index[dv].erase(tid);
						//if (dv != pv)
						//	in_tree_count[dv]--;

						//	double u81_acc_time = GetCurrentTimeSec();
						if (affected_list[dv] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(dv);
							affected_list[dv] = true;
						}

						//u81_time += (GetCurrentTimeSec() - u81_acc_time);

						//	u82_time += (GetCurrentTimeSec() - u82_acc_time);
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;
						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else
							remove_small_node(dv, stree_tid);
					}
					u8_time += (GetCurrentTimeSec() - u8_acc_time);


					current_universe -= descendants.size();

					if (st_roots[tid] == pv)
						if (isLarge) {
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
						}
						else {
							stree_tid.clear();
							stree_tid.resize(0);
						}


						st_sizes[tid] -= subtract_size;
						if (trees_pointers[tid] != -1) alive_resources -= (long long)subtract_size;
						//	if (st_sizes[tid] < 0) {
						//		cout << "a5" << st_sizes[tid] << endl;
						//	}

						// 3. Convert shrink trees to small trees representataions.				
						if (st_sizes[tid] > 0) {
							if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
								shrink_tree(tid, st_roots[tid], ltrees[tid]);
							}
						}

						//u2_time += (GetCurrentTimeSec() - u2_acc_time);
				}
				reverse_index[pv].clear();
				reverse_index[pv].resize(0);

			}
			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			// 4. Sample new trees based on.
			// If we sample a new tree, we need to dynamically increment the following data structures:
			// 1) ltrees
			// 2) trees_pointers (-1 or pointer to strees)
			// 3) cover
			// 4) st_roots
			// 5) st_sizes 
			// 6) strees (maybe)
			// two constraints: cl > ct and 10kn total vertices
			double adding_acc_time = GetCurrentTimeSec();


			/*cout << "ranking:" << ranking << endl;
			cout << "pick cover:" << pv_cover << endl;
			cout << "pick count:" << picked_count << endl;
			cout << "labeled_arcs:" << labeled_arcs << endl;
			cout << "new_tree_arcs:" << new_tree_arcs << endl;
			cout << "alive_resources:" << alive_resources << endl;
			cout << "bounded_resources:" << bounded_resources << endl;
			cout << "#tree:" << cover.size() << endl;
			cout << "u5/adding" << u5_time / adding_time << endl;
			*/
			new_tree_arcs = 0;
			int how_many = 0;


			cover_value_type sample_constraint = c_paras * current_universe / st_roots.size();
			cover_value_type sum_spv = sum[spv];

			while ( beta * labeled_arcs > new_tree_arcs && alive_resources < bounded_resources && sum_spv < sample_constraint) {

				if (picked_count == numOfVertices) break;

				NodeID r = rand() % numOfVertices;
				if (picked_as_root[r] == true) continue;
				if (trivial_root[r] == true) continue;

				large_tree new_ltree;
				small_tree new_stree;
				//GOOGLE DENSE MAP
				new_stree.set_empty_key(numOfVertices + 1);
				new_stree.set_deleted_key(numOfVertices + 2);
				bool isLarge = false;
				int tmp_tree_arcs = 0;

				if (!BUILD_SMALL_TREE_FLAG) {

					//double u5_acc_time = GetCurrentTimeSec();
					double u9_acc_time = GetCurrentTimeSec();
					tmp_tree_arcs = build_large_tree_dij(r, new_ltree, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, new_tree_arcs, beta * labeled_arcs);
					u9_time += (GetCurrentTimeSec() - u9_acc_time);
					//u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;

					if (descendants.size() < size_thr)
						BUILD_SMALL_TREE_FLAG = true;

					if (descendants.size() < 100) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					//	double u91_acc_time = GetCurrentTimeSec();

					isLarge = true;
					how_many++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					used_as_cover++;
					st_roots.push_back(r);
					trees_pointers.push_back(-1);
					st_sizes.push_back(descendants.size());

					//new_tree_arcs += descendants.size();
					alive_resources += (long long)numOfVertices;

					//	u91_time += (GetCurrentTimeSec() - u91_acc_time);
					// Turn large tree to small tree.
					// 1. Build new small tree with ltree_i and add it to strees;
					// 2. Set tree_pointers[i] to the position of the strees;
					// 3. Clear the large tree ltree_i by setting its size to 0.

					//	double u92_acc_time = GetCurrentTimeSec();
					int i = st_roots.size() - 1;
					if (st_sizes[i] < size_thr) {
						shrink_tree(i, st_roots[i], new_ltree);
						// FIX
						isLarge = false;
						/*for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID q = descendants[di];
						NodeID v_parent = get_large_parent(q, new_ltree);
						NodeID s_parent = get_small_parent(q, strees[strees.size() - 1]);
						if (v_parent != s_parent) cout << v_parent << " vs " << s_parent << endl;
						}*/
					}
					ltrees.push_back(new_ltree);
					//	u92_time += (GetCurrentTimeSec() - u92_acc_time);
				}
				else {//Build Small Tree Directly
					  //double u5_acc_time = GetCurrentTimeSec();
					double u9_acc_time = GetCurrentTimeSec();
					tmp_tree_arcs = build_small_tree_dij(r, new_stree, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, new_tree_arcs, beta * labeled_arcs);
					//u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;
					u9_time += (GetCurrentTimeSec() - u9_acc_time);

					if (descendants.size() >= size_thr)
						BUILD_SMALL_TREE_FLAG = false;

					if (descendants.size() < 100) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					isLarge = false;

					how_many++;
					used_as_cover++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					st_roots.push_back(r);
					trees_pointers.push_back(strees.size()); //we assume it is inserted.
					st_sizes.push_back(descendants.size());

					alive_resources += (long long)new_stree.size();
					//strees.push_back(new_stree);
					//// new_ltree is empty.
					//ltrees.push_back(new_ltree);

					int i = st_roots.size() - 1;
					if (st_sizes[i] >= size_thr) {
						enlarge_tree(i, st_roots[i], new_stree, new_ltree);
						new_stree.clear();
						new_stree.resize(0);
						isLarge = true;
					}
					else {
						strees.push_back(new_stree);
						ltrees.push_back(new_ltree);
					}

				}

				int i = st_roots.size() - 1;

				int numOfStrees = strees.size();
				small_tree &new_stree_ref = numOfStrees != 0 ? strees[strees.size() - 1] : empty_small_tree;

				// Calculate the subtree size of each nodes.

				double u10_acc_time = GetCurrentTimeSec();
				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];

					++cover_tid[q]; // Increment because it hitsthe r-q path.
					++cnt16[q][i%CNT];
					//++sum[q];
					//affected[q] = true;

					if (affected_list[q] == false) {
						affected_count++;
						if (switch_set)
							affected_v.insert(q);
						affected_list[q] = true;
					}

					if (switch_small == true) {
						reverse_index[q].insert(i);
						//in_tree_count[q]++;
					}


					int parent_node;
					if (isLarge)
						parent_node = get_large_parent(q, new_ltree);
					else {
						//double get_small_parent_incre_time = GetCurrentTimeSec();
						parent_node = get_small_parent(q, new_stree_ref);
						//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}
					if (parent_node > -1) {
						cover_tid[parent_node] += cover_tid[q];
						cnt16[parent_node][i%CNT] += cover_tid[q];
						//sum[parent_node] += cover_new[q];
					}
				}
				u10_time += (GetCurrentTimeSec() - u10_acc_time);
				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];
					cover_tid[q] = 0;
				}


				current_universe += descendants.size();
				sample_constraint = c_paras * current_universe / st_roots.size();


				//cout << "beta*labeld_arcs: " << beta * labeled_arcs << " new tree _arcs: " << new_tree_arcs << " alive " << alive_resources << " bounded " << bounded_resources << endl;

			}
			adding_time += (GetCurrentTimeSec() - adding_acc_time);


			//double u3_acc_time = GetCurrentTimeSec();

			updating_acc_time = GetCurrentTimeSec();

			// Build the cover value heap for one time.
			if (switch_small == true && value_heap.empty())
				for (int nid = 0; nid < numOfVertices; ++nid)
					if (has_picked[nid] == false)
						value_heap.update(nid, -picked_value[nid]);

			//	cnt16_2max(picked_value, sum, has_picked, affected, cnt16, switch_small, value_heap);
			//c	cnt16_2max_queue(picked_value, sum, has_picked, affected_v, cnt16, switch_small, value_heap);

			double u11_acc_time = GetCurrentTimeSec();
			cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);
			u11_time += (GetCurrentTimeSec() - u11_acc_time);
			//	iteration_time[ranking] = (GetCurrentTimeSec() - updating_acc_time);
			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			if (affected_count < numOfVertices / 8)
				switch_set = true;
			else
				switch_set = false;
			affected_count = 0;

			time_this_iteration = GetCurrentTimeSec() - time_this_iteration;
			iteration_time[ranking] = time_this_iteration;

		}

		cout << "picked as tree covers:" << used_as_cover << endl;
		cout << "picked count:" << picked_count << endl;

		//// FIX
		//for (int v = 0; v < numOfVertices; ++v){
		//	for (int i = 0; i < CNT; ++i) {
		//		if (cnt16[v][i] != 0)
		//			cout << "mi?" << endl;
		//	}
		//}
		if (stop_flag == true) return;

		double labeling_acc_time = GetCurrentTimeSec();
		// Order the remaining vertices with degree, if any.
		if (ranking != numOfVertices) {
			vector<NodeID> degree(numOfVertices);
			vector<NodeID> r_degree(numOfVertices);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				/*vector<NodeID>& adj_v = graph.adj[v];
				for (NodeID i = 0; i < adj_v.size(); ++i) {
				if (has_picked[adj_v[i]] != true)
				degree[v]++;
				}*/

				for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid) {
					NodeID w = wgraph.edges[eid].first;
					if (has_picked[w] != true)
						degree[v]++;
				}

				if (DIRECTED_FLAG == true) {
					/*vector<NodeID>& r_adj_v = graph.r_adj[v];
					for (NodeID i = 0; i < r_adj_v.size(); ++i) {
					if (has_picked[r_adj_v[i]] != true)
					r_degree[v]++;
					}*/
					for (EdgeID eid = wgraph.r_vertices[v]; eid < wgraph.r_vertices[v + 1]; ++eid) {
						NodeID w = wgraph.r_edges[eid].first;
						if (has_picked[w] != true)
							r_degree[v]++;
					}
				}
			}

			vector<pair<NodeID, NodeID> > degree_v;
			degree_v.reserve(numOfVertices - ranking);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				if (DIRECTED_FLAG == true)
					degree_v.push_back(make_pair(degree[v] * r_degree[v], v));
				else
					degree_v.push_back(make_pair(degree[v], v));
			}
			sort(degree_v.rbegin(), degree_v.rend());

			int count = 0;
			for (; ranking < numOfVertices; ++ranking) {
				NodeID tpv = degree_v[count++].second;
				inv[ranking] = tpv;
				rank[tpv] = ranking;
				labeling_source_dij(tpv, pqueue, visited_que, distances, vis, dst_r, usd, ranking, tmp_idx, r_tmp_idx, wgraph);
			}
		}
		labeling_time += (GetCurrentTimeSec() - labeling_acc_time);


		// Finish the constructions.
		for (NodeID v = 0; v < numOfVertices; ++v) {
			if (DIRECTED_FLAG == false) {
				NodeID k = tmp_idx[v].first.size();
				labels.index_[v].spt_v.resize(k);
				labels.index_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();
				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();

			}
			else {
				NodeID k = tmp_idx[v].first.size();
				dlabels.index_[v].spt_v.resize(k);
				dlabels.index_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();
				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();

				k = r_tmp_idx[v].first.size();
				dlabels.bindex_[v].spt_v.resize(k);
				dlabels.bindex_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_v[i] = r_tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_d[i] = r_tmp_idx[v].second[i];
				r_tmp_idx[v].first.clear();
				r_tmp_idx[v].second.clear();
				r_tmp_idx[v].first.shrink_to_fit();
				r_tmp_idx[v].second.shrink_to_fit();
			}
		}

		num_of_trees = st_roots.size();

		double small_tree_num = 0;

		for (int i = 0; i < num_of_trees; ++i) {
			if (trees_pointers[i] != -1) {
				small_tree_num++;
				//small_tree_num += strees[trees_pointers[i]].size();
			}
		}

		total_resources = small_tree_num;

	}

	Betweenness_Ordering(NodeID k, double beta, WGraph& wgraph, NodeID stopv, bool path_flag) {
		 
		double start_time = GetCurrentTimeSec();

		if (k < CNT) {
			cout << " Should pick at least 16 seed!" << endl;
			return;
		}
		iteration_time.resize(numOfVertices);

		long long bound_times = 10;
		long long trans_times = 8;
		updating_time = 0;
		labeling_time = 0;
		adding_time = 0;
		selecting_time = 0;
		init_time = 0;
		u1_time = 0;
		u2_time = 0;
		u3_time = 0;
		num_of_trees = 0;
		total_resources = 0;


		inv.resize(numOfVertices, -1);
		rank.resize(numOfVertices, -1);

		if (DIRECTED_FLAG == false) {
			plabels.index_.resize(numOfVertices);
		}
		else { //if (DIRECTED_FLAG == true) {
			dplabels.index_.resize(numOfVertices);
			dplabels.bindex_.resize(numOfVertices);
		}

		// Resource Bounding.
		NodeID new_tree_arcs = 0;
		NodeID labeled_arcs = 0;
		alive_resources = 0;
		bounded_resources = bound_times * (long long)k * (long long)numOfVertices;

		vector<NodeID> st_roots; // The root vertices of the sampled trees.
		vector<NodeID> st_sizes; // The size of the sampled trees.

		NodeID size_thr = numOfVertices < 8 ? numOfVertices : numOfVertices / trans_times; // When a tree shrink to a size smaller than the size_thr, it should be converted into a small_tree representation.

		init_time = GetCurrentTimeSec();
		// Building the initial k trees.
		//	srand((unsigned)time(NULL));
		srand(0x728f9d0);


		ltrees.reserve(numOfVertices);
		strees.reserve(numOfVertices);
		trees_pointers.reserve(numOfVertices);
		st_roots.reserve(numOfVertices);
		st_sizes.reserve(numOfVertices);
		st_roots.resize(k);
		st_sizes.resize(k);

		ltrees.resize(k);
		//strees.clear();
		trees_pointers.resize(k, -1); // Point the sampling trees to its tree representations. large tree [0, V) small tree [V, 2V).
									  //vector<unordered_map<int, int> > cover(k, unordered_map<int, int>());

		vector<NodeID> ascendants;
		ascendants.reserve(numOfVertices);
		vector<NodeID> descendants;
		descendants.reserve(numOfVertices);
		vector<NodeID> que(numOfVertices);
		vector<bool> vis(numOfVertices);
		vector<NodeID> parents(numOfVertices, numOfVertices);
		vector<bool> has_picked(numOfVertices);
		vector<vector<cover_value_type> > cnt16(numOfVertices, vector<cover_value_type>(CNT));
		//		vector<bool> affected(V, true);
		vector<bool> picked_as_root(numOfVertices, false);
		vector<bool> trivial_root(numOfVertices, false);
		int picked_count = 0;
		int used_as_cover = 0;
		vector<cover_value_type> sum(numOfVertices);
		bool stop_flag = false;


		// tmp structure used in the dijkstra searches. 

		// Binary Heap
		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);

		queue<NodeID> visited_que;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
		vector<bool> usd(numOfVertices, false);

		//GOOGLE
		//		unordered_set<NodeID> affected_v;
		google::dense_hash_set<NodeID> affected_v;
		affected_v.set_empty_key(numOfVertices + 1);
		affected_v.set_deleted_key(numOfVertices + 2);

		vector<cover_value_type> cover_tid(numOfVertices, 0);
		//vector<cover_value_type> cover_new(numOfVertices, 0); // modified-5-14


		vector<bool> affected_list(numOfVertices, false);
		int affected_count = 0;
		bool switch_set = false;

		//vector<unordered_set<int> > reverse_index(0);
		//GOOGLE DENSE SET
		vector<google::dense_hash_set<NodeID> > reverse_index;
		//vector<NodeID> in_tree_count(numOfVertices);

		// Preparing basic structure for pl algorithm.
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT))); // Backward labels.


		vector<pair<vector<NodeID>, vector<NodeID> > >
			tmp_idx_parent(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<NodeID>(1, numOfVertices)));

		vector<pair<vector<NodeID>, vector<NodeID> > >
			r_tmp_idx_parent(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<NodeID>(1, numOfVertices))); // Backward labels.

		vector<NodeID> descendants_parents(numOfVertices, -1);

		// Building the initial k sampled trees. 
		for (int i = 0; i < k; ++i) {
			NodeID r = rand() % numOfVertices;
			if (picked_as_root[r] == true) { --i; continue; }

			st_roots[i] = r;
			large_tree &ltree_i = ltrees[i];
			picked_as_root[r] = true;
			picked_count++;
			used_as_cover++;

			build_large_tree_dij(r, ltree_i, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, 0, numOfVertices);
			st_sizes[i] = descendants.size();
			alive_resources += (long long)numOfVertices;

			// Turn large tree to small tree.
			// 1. Build new small tree with ltree_i and add it to strees;
			// 2. Set tree_pointers[i] to the position of the strees;
			// 3. Clear the large tree ltree_i by setting its size to 0.
			if (st_sizes[i] < size_thr)
				shrink_tree(i, st_roots[i], ltrees[i]);

			bool isLarge = false;
			int stree_idx = trees_pointers[i];
			if (stree_idx == -1)
				isLarge = true;

			small_tree& stree_i = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;

			// Calculate the subtree size of each nodes.
			//	unordered_map<int, int>& cover_i = cover[i];
			//vector<int> cover_i(numOfVertices, 0);

			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];

				++cover_tid[q]; // Increment because it hits the r-q path.
				++cnt16[q][i%CNT];
				//++sum[q];

				int parent_node;

				if (isLarge)
					parent_node = get_large_parent(q, ltree_i);
				else
					parent_node = get_small_parent(q, stree_i);

				if (parent_node > -1) {
					cover_tid[parent_node] += cover_tid[q];
					cnt16[parent_node][i%CNT] += cover_tid[q];
					//	sum[parent_node] += cover_i[q];
				}
			}
			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];
				cover_tid[q] = 0;
			}
		}



		//cnt16_aggreate(cnt16, has_picked, cover);

		vector<cover_value_type> picked_value;
		benchmark::heap<2, cover_value_type, NodeID> value_heap(numOfVertices);
		//bool switch_small = false;
		switch_small = false;

		//cnt16_2max(picked_value, sum, has_picked,affected, cnt16, switch_small, value_heap);
		//		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, cnt16, switch_small, value_heap);
		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);

		init_time = GetCurrentTimeSec() - init_time;

		cout << "Finish initial sampling " << init_time << "secs." << endl;

		//for (NodeID v = 0; v < V; ++v)	value_heap.update(v, -picked_value[v]); // Min-heap.

		double max_value = -1;

		// Perform V iterations to rank all vertices.
		NodeID ranking = 0;
		for (; ranking < numOfVertices; ++ranking) {
			if (ranking % (numOfVertices / 10) == 0)
				cout << "picking the " << ranking << "-th orders." << endl;

			double time_this_iteration = GetCurrentTimeSec();

			NodeID pv = 0;
			long long pv_cover = -1;

			if (ranking == stopv) {
				stop_flag = true;
				cout << "topv labeling time:" << labeling_time << endl;
				break;
			}

			double selecting_acc_time = GetCurrentTimeSec();
			// In the beginning, pick the best cover vertex by a linear scan.
			if (switch_small == false) {
				double u1_acc_time = GetCurrentTimeSec();

				for (NodeID v = 0; v < numOfVertices; ++v) {
					if (has_picked[v] == true) continue;
					if (pv_cover == picked_value[v]) {
						if (sum[pv] < sum[v]) {
							pv_cover = picked_value[v];
							pv = v;
						}
					}
					else if (pv_cover < picked_value[v]) {
						pv_cover = picked_value[v];
						pv = v;
					}
				}

				if (sum[pv] < numOfVertices / 8 && switch_small == false) {
					//cout << "u3 time:" << u3_time << endl;
					cout << ranking << endl;
					cout << pv_cover << endl;
					cout << sum[pv] << endl;
					cout << picked_count << endl;
					cout << "Switch to max-heap representations." << endl;
					switch_small = true;
				}
				u1_time += (GetCurrentTimeSec() - u1_acc_time);

			}
			else {// Later, utilizing the max-heap directly
				double u2_acc_time = GetCurrentTimeSec();

				value_heap.extract_min(pv, pv_cover);
				pv_cover = -pv_cover;
				u2_time += (GetCurrentTimeSec() - u2_acc_time);
				/*
				int max1 = 0;
				int max2 = 0;
				int max1v = -1;
				int max2v = -1;
				for(int  i = 0; i < CNT; ++i){
				if(max2 <= cnt16[pv][i]){
				max2 = cnt16[pv][i];
				max2v = i;
				}
				if(max1 <= max2){
				int tmp = max1;
				max1 = max2;
				max2 = tmp;
				int tmpv = max1v;
				max1v = max2v;
				max2v = max1v;
				}
				}
				long long sumstd = 0;
				long long summean = 0;
				for(int  i = 0; i < CNT; ++i){
				if(i != max1v && i != max2v){
				sumstd += (long long)cnt16[pv][i] * (long long)cnt16[pv][i];
				summean += (long long)cnt16[pv][i];
				}
				}
				sumstd = sumstd/(long long)(CNT-2);
				summean = summean/(long long)(CNT-2);
				cout << sqrt(sumstd - summean * summean) << endl;*/
			}

			selecting_time += (GetCurrentTimeSec() - selecting_acc_time);

			//// FIX
			//if (switch_small == true && reverse_index.size() > 0 && sum[pv] <= reverse_index[pv].size() * 10) {//sum[pv] <= in_tree_count[pv] * 10){//
			//	cout << "I think we can stop here " << GetCurrentTimeSec() - start_time << " ranking:" << ranking << endl;
			//	break;
			//}

			//if (has_picked[pv] == true)
			//	cout << "wrong1" << endl;

			has_picked[pv] = true;
			double labeling_acc_time = GetCurrentTimeSec();
			labeled_arcs = labeling_source_dij_path(pv, pqueue, visited_que, distances, vis, parents, dst_r, usd, ranking, tmp_idx, r_tmp_idx, tmp_idx_parent, r_tmp_idx_parent, wgraph);
			labeling_time += (GetCurrentTimeSec() - labeling_acc_time);
			inv[ranking] = pv;
			rank[pv] = ranking;

			int current_trees_num = st_roots.size();

			// 4 steps in updating:
			double updating_acc_time = GetCurrentTimeSec();
			if (switch_small == false) { // We still iterate all alive trees to calculate the cover.	
										 //	double u3_acc_time = GetCurrentTimeSec();
				for (int tid = 0; tid < current_trees_num; ++tid) {

					// Get the descendants and ascendants of pv in tree tid.
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					// pv has already been covered in tree tid.
					int pvpid = 0;
					if (isLarge)
						pvpid = get_large_parent(pv, ltree_tid);
					else {
						//double get_small_parent_incre_time = GetCurrentTimeSec();
						pvpid = get_small_parent(pv, stree_tid);
						//	get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}
					if (pvpid == -1 && st_roots[tid] != pv) // pv has been covered in tree tid 
						continue;


					double u3_acc_time = GetCurrentTimeSec();

					// //1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;				 
					//get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph.adj, isLarge, stree_tid, ltree_tid);
					//get_ascendants(tid, pv, ascendants, que, vis, wgraph.adj, wgraph.r_adj, isLarge, stree_tid, ltree_tid);

					//Array Representation
					get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph, isLarge, stree_tid, ltree_tid);
					get_ascendants(tid, pv, ascendants, que, vis, wgraph, isLarge, stree_tid, ltree_tid);

					u3_time += (GetCurrentTimeSec() - u3_acc_time);


					double u4_acc_time = GetCurrentTimeSec();
					NodeID subtract_size = descendants.size();
					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
												//	cover[tid][av] -= subtract_size;
						cnt16[av][tid%CNT] -= subtract_size;
						//sum[av] -= subtract_size;
						//affected[av] = true;

						if (affected_list[av] == false) {
							affected_list[av] = true;
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
						}
					}
					u4_time += (GetCurrentTimeSec() - u4_acc_time);


					double u5_acc_time = GetCurrentTimeSec();
					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					//vector<int> cover_tid(numOfVertices, 0);
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//sum[dv] -= cover_tid[dv];

						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];
						if (cnt16[dv][tid%CNT] < 0)
							cout << "a6 " << cnt16[dv][tid%CNT] << endl;

						/*if (isLarge)
						remove_large_node(dv, ltree_tid);
						else
						remove_small_node(dv, stree_tid);*/

						//						remove_node(tid, dv);
						//affected[dv] = true;

						if (affected_list[dv] == false) {
							affected_count++;
							affected_list[dv] = true;
							if (switch_set)
								affected_v.insert(dv);
						}
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;

						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else
							remove_small_node(dv, stree_tid);
					}
					u5_time += (GetCurrentTimeSec() - u5_acc_time);

					// Drop the size of subtree.
					st_sizes[tid] -= subtract_size;
					if (trees_pointers[tid] != -1) alive_resources -= (long long)subtract_size;

					// 3. Convert shrink trees to small trees representataions.	
					if (st_sizes[tid] > 0) {
						if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
							shrink_tree(tid, st_roots[tid], ltrees[tid]);
						}
					}

					if (st_roots[tid] == pv)
						if (isLarge) {
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
							//ltrees.shrink_to_fit();
						}
						else {
							stree_tid.clear();
							stree_tid.resize(0);
							//strees.shrink_to_fit();
						}

				}
				//	u3_time += (GetCurrentTimeSec() - u3_acc_time);
			}
			else { // We use the reverse index to update trees.
				   // Build the reverse_index first (one time).
				if (reverse_index.size() == 0) {
					reverse_index.resize(numOfVertices);
					for (NodeID v = 0; v < numOfVertices; ++v) {
						reverse_index[v].set_empty_key(numOfVertices + 1);
						reverse_index[v].set_deleted_key(numOfVertices + 2);
					}
					for (int tid = 0; tid < current_trees_num; ++tid) {
						if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue; //This whole tree has already been covered. 

						bool isLarge = false;
						int stree_idx = trees_pointers[tid];
						if (stree_idx == -1)
							isLarge = true;

						small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
						large_tree& ltree_tid = ltrees[tid];

						for (int v = 0; v < numOfVertices; ++v) {

							int vpid = 0;
							if (isLarge)
								vpid = get_large_parent(v, ltree_tid);
							else {
								//double get_small_parent_incre_time = GetCurrentTimeSec();
								vpid = get_small_parent(v, stree_tid);
								//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
							}
							if ((vpid != -1 || st_roots[tid] == v)) {
								reverse_index[v].insert(tid);
								//in_tree_count[v]++;
							}
						}
					}
				}

				//cout << reverse_index[pv].size() << "," << sum[pv] << "," <<  (double)sum[pv] / (double)reverse_index[pv].size() << endl;
				// Access the alive trees directly by the reverse index.
				//for (unordered_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) {
				for (google::dense_hash_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) {
					int tid = *it;
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					double u6_acc_time = GetCurrentTimeSec();
					/*
					int parent_node;
					if (isLarge)
					parent_node = get_large_parent(pv, ltree_tid);
					else {
					//double get_small_parent_incre_time = GetCurrentTimeSec();
					parent_node = get_small_parent(pv, stree_tid);
					//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}
					if(parent_node == -1 && pv != st_roots[tid])
					continue;
					*/
					double u7_acc_time = GetCurrentTimeSec();

					//	// 1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;
					//	//double u4_acc_time = getcurrenttimesec();
					//	get_descendants(tid, pv, descendants, descendants_parents,que, vis, wgraph.adj, islarge, stree_tid, ltree_tid);
					////	u4_time += (getcurrenttimesec() - u4_acc_time);

					//	get_ascendants(tid, pv, ascendants, que, vis, wgraph.adj, wgraph.r_adj, isLarge, stree_tid, ltree_tid);

					// Array Representation
					get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph, isLarge, stree_tid, ltree_tid);
					get_ascendants(tid, pv, ascendants, que, vis, wgraph, isLarge, stree_tid, ltree_tid);

					NodeID subtract_size = descendants.size();
					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
												//cover[tid][av] -= subtract_size;
						cnt16[av][tid%CNT] -= subtract_size;
						//affected[av] = true;
						//	sum[av] -= subtract_size;

						if (affected_list[av] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
							affected_list[av] = true;
						}
					}
					u7_time += (GetCurrentTimeSec() - u7_acc_time);

					//double u2_acc_time = GetCurrentTimeSec();
					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					double u8_acc_time = GetCurrentTimeSec();
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//	sum[dv] -= cover_tid[dv];

						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];
						//	if (cnt16[dv][tid%CNT] < 0)
						//		cout << "a6 " << cnt16[dv][tid%CNT] << endl;

						//if(dv == pv){
						/*if (isLarge)
						remove_large_node(dv, ltree_tid);
						else
						remove_small_node(dv, stree_tid);*/
						//}
						//affected[dv] = true;

						//double u82_acc_time = GetCurrentTimeSec();
						if (dv != pv)
							reverse_index[dv].erase(tid);
						//if (dv != pv)
						//	in_tree_count[dv]--;

						//	double u81_acc_time = GetCurrentTimeSec();
						if (affected_list[dv] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(dv);
							affected_list[dv] = true;
						}

						//u81_time += (GetCurrentTimeSec() - u81_acc_time);

						//	u82_time += (GetCurrentTimeSec() - u82_acc_time);
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;
						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else
							remove_small_node(dv, stree_tid);
					}
					u8_time += (GetCurrentTimeSec() - u8_acc_time);

					if (st_roots[tid] == pv)
						if (isLarge) {
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
						}
						else {
							stree_tid.clear();
							stree_tid.resize(0);
						}


						st_sizes[tid] -= subtract_size;
						if (trees_pointers[tid] != -1) alive_resources -= (long long)subtract_size;
						//	if (st_sizes[tid] < 0) {
						//		cout << "a5" << st_sizes[tid] << endl;
						//	}

						// 3. Convert shrink trees to small trees representataions.				
						if (st_sizes[tid] > 0) {
							if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
								shrink_tree(tid, st_roots[tid], ltrees[tid]);
							}
						}

						//u2_time += (GetCurrentTimeSec() - u2_acc_time);
				}
				reverse_index[pv].clear();
				reverse_index[pv].resize(0);

			}
			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			// 4. Sample new trees based on.
			// If we sample a new tree, we need to dynamically increment the following data structures:
			// 1) ltrees
			// 2) trees_pointers (-1 or pointer to strees)
			// 3) cover
			// 4) st_roots
			// 5) st_sizes 
			// 6) strees (maybe)
			// two constraints: cl > ct and 10kn total vertices
			double adding_acc_time = GetCurrentTimeSec();


			/*cout << "ranking:" << ranking << endl;
			cout << "pick cover:" << pv_cover << endl;
			cout << "pick count:" << picked_count << endl;
			cout << "labeled_arcs:" << labeled_arcs << endl;
			cout << "new_tree_arcs:" << new_tree_arcs << endl;
			cout << "alive_resources:" << alive_resources << endl;
			cout << "bounded_resources:" << bounded_resources << endl;
			cout << "#tree:" << cover.size() << endl;
			cout << "u5/adding" << u5_time / adding_time << endl;
			*/
			new_tree_arcs = 0;
			int how_many = 0;
			while (beta * labeled_arcs > new_tree_arcs && alive_resources < bounded_resources) {

				if (picked_count == numOfVertices) break;

				NodeID r = rand() % numOfVertices;
				if (picked_as_root[r] == true) continue;
				if (trivial_root[r] == true) continue;

				large_tree new_ltree;
				small_tree new_stree;
				//GOOGLE DENSE MAP
				new_stree.set_empty_key(numOfVertices + 1);
				new_stree.set_deleted_key(numOfVertices + 2);
				bool isLarge = false;
				int tmp_tree_arcs = 0;

				if (!BUILD_SMALL_TREE_FLAG) {

					//double u5_acc_time = GetCurrentTimeSec();
					double u9_acc_time = GetCurrentTimeSec();
					tmp_tree_arcs = build_large_tree_dij(r, new_ltree, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, new_tree_arcs, beta * labeled_arcs);
					u9_time += (GetCurrentTimeSec() - u9_acc_time);
					//u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;

					if (descendants.size() < size_thr)
						BUILD_SMALL_TREE_FLAG = true;

					if (descendants.size() < 100) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					//	double u91_acc_time = GetCurrentTimeSec();

					isLarge = true;
					how_many++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					used_as_cover++;
					st_roots.push_back(r);
					trees_pointers.push_back(-1);
					st_sizes.push_back(descendants.size());

					//new_tree_arcs += descendants.size();
					alive_resources += (long long)numOfVertices;

					//	u91_time += (GetCurrentTimeSec() - u91_acc_time);
					// Turn large tree to small tree.
					// 1. Build new small tree with ltree_i and add it to strees;
					// 2. Set tree_pointers[i] to the position of the strees;
					// 3. Clear the large tree ltree_i by setting its size to 0.

					//	double u92_acc_time = GetCurrentTimeSec();
					int i = st_roots.size() - 1;
					if (st_sizes[i] < size_thr) {
						shrink_tree(i, st_roots[i], new_ltree);
						// FIX
						isLarge = false;
						/*for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID q = descendants[di];
						NodeID v_parent = get_large_parent(q, new_ltree);
						NodeID s_parent = get_small_parent(q, strees[strees.size() - 1]);
						if (v_parent != s_parent) cout << v_parent << " vs " << s_parent << endl;
						}*/
					}
					ltrees.push_back(new_ltree);
					//	u92_time += (GetCurrentTimeSec() - u92_acc_time);
				}
				else {//Build Small Tree Directly
					  //double u5_acc_time = GetCurrentTimeSec();
					double u9_acc_time = GetCurrentTimeSec();
					tmp_tree_arcs = build_small_tree_dij(r, new_stree, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, new_tree_arcs, beta * labeled_arcs);
					//u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;
					u9_time += (GetCurrentTimeSec() - u9_acc_time);

					if (descendants.size() >= size_thr)
						BUILD_SMALL_TREE_FLAG = false;

					if (descendants.size() < 100) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					isLarge = false;

					how_many++;
					used_as_cover++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					st_roots.push_back(r);
					trees_pointers.push_back(strees.size()); //we assume it is inserted.
					st_sizes.push_back(descendants.size());

					alive_resources += (long long)new_stree.size();
					//strees.push_back(new_stree);
					//// new_ltree is empty.
					//ltrees.push_back(new_ltree);

					int i = st_roots.size() - 1;
					if (st_sizes[i] >= size_thr) {
						enlarge_tree(i, st_roots[i], new_stree, new_ltree);
						new_stree.clear();
						new_stree.resize(0);
						isLarge = true;
					}
					else {
						strees.push_back(new_stree);
						ltrees.push_back(new_ltree);
					}

				}

				int i = st_roots.size() - 1;

				int numOfStrees = strees.size();
				small_tree &new_stree_ref = numOfStrees != 0 ? strees[strees.size() - 1] : empty_small_tree;

				// Calculate the subtree size of each nodes.

				double u10_acc_time = GetCurrentTimeSec();
				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];

					++cover_tid[q]; // Increment because it hitsthe r-q path.
					++cnt16[q][i%CNT];
					//++sum[q];
					//affected[q] = true;

					if (affected_list[q] == false) {
						affected_count++;
						if (switch_set)
							affected_v.insert(q);
						affected_list[q] = true;
					}

					if (switch_small == true) {
						reverse_index[q].insert(i);
						//in_tree_count[q]++;
					}


					int parent_node;
					if (isLarge)
						parent_node = get_large_parent(q, new_ltree);
					else {
						//double get_small_parent_incre_time = GetCurrentTimeSec();
						parent_node = get_small_parent(q, new_stree_ref);
						//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}
					if (parent_node > -1) {
						cover_tid[parent_node] += cover_tid[q];
						cnt16[parent_node][i%CNT] += cover_tid[q];
						//sum[parent_node] += cover_new[q];
					}
				}
				u10_time += (GetCurrentTimeSec() - u10_acc_time);
				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];
					cover_tid[q] = 0;
				}
				//cout << "beta*labeld_arcs: " << beta * labeled_arcs << " new tree _arcs: " << new_tree_arcs << " alive " << alive_resources << " bounded " << bounded_resources << endl;

			}
			adding_time += (GetCurrentTimeSec() - adding_acc_time);


			//double u3_acc_time = GetCurrentTimeSec();

			updating_acc_time = GetCurrentTimeSec();

			// Build the cover value heap for one time.
			if (switch_small == true && value_heap.empty())
				for (int nid = 0; nid < numOfVertices; ++nid)
					if (has_picked[nid] == false)
						value_heap.update(nid, -picked_value[nid]);

			//	cnt16_2max(picked_value, sum, has_picked, affected, cnt16, switch_small, value_heap);
			//c	cnt16_2max_queue(picked_value, sum, has_picked, affected_v, cnt16, switch_small, value_heap);

			double u11_acc_time = GetCurrentTimeSec();
			cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);
			u11_time += (GetCurrentTimeSec() - u11_acc_time);
			//	iteration_time[ranking] = (GetCurrentTimeSec() - updating_acc_time);
			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			if (affected_count < numOfVertices / 8)
				switch_set = true;
			else
				switch_set = false;
			affected_count = 0;

			time_this_iteration = GetCurrentTimeSec() - time_this_iteration;
			iteration_time[ranking] = time_this_iteration;

		}

		cout << "picked as tree covers:" << used_as_cover << endl;
		cout << "picked count:" << picked_count << endl;

		//// FIX
		//for (int v = 0; v < numOfVertices; ++v){
		//	for (int i = 0; i < CNT; ++i) {
		//		if (cnt16[v][i] != 0)
		//			cout << "mi?" << endl;
		//	}
		//}
		if (stop_flag == true) return;

		double labeling_acc_time = GetCurrentTimeSec();
		// Order the remaining vertices with degree, if any.
		if (ranking != numOfVertices) {
			vector<NodeID> degree(numOfVertices);
			vector<NodeID> r_degree(numOfVertices);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				/*vector<NodeID>& adj_v = graph.adj[v];
				for (NodeID i = 0; i < adj_v.size(); ++i) {
				if (has_picked[adj_v[i]] != true)
				degree[v]++;
				}*/

				for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid) {
					NodeID w = wgraph.edges[eid].first;
					if (has_picked[w] != true)
						degree[v]++;
				}

				if (DIRECTED_FLAG == true) {
					/*vector<NodeID>& r_adj_v = graph.r_adj[v];
					for (NodeID i = 0; i < r_adj_v.size(); ++i) {
					if (has_picked[r_adj_v[i]] != true)
					r_degree[v]++;
					}*/
					for (EdgeID eid = wgraph.r_vertices[v]; eid < wgraph.r_vertices[v + 1]; ++eid) {
						NodeID w = wgraph.r_edges[eid].first;
						if (has_picked[w] != true)
							r_degree[v]++;
					}
				}
			}

			vector<pair<NodeID, NodeID> > degree_v;
			degree_v.reserve(numOfVertices - ranking);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				if (DIRECTED_FLAG == true)
					degree_v.push_back(make_pair(degree[v] * r_degree[v], v));
				else
					degree_v.push_back(make_pair(degree[v], v));
			}
			sort(degree_v.rbegin(), degree_v.rend());

			int count = 0;
			for (; ranking < numOfVertices; ++ranking) {
				NodeID tpv = degree_v[count++].second;
				inv[ranking] = tpv;
				rank[tpv] = ranking;
				labeling_source_dij(tpv, pqueue, visited_que, distances, vis, dst_r, usd, ranking, tmp_idx, r_tmp_idx, wgraph);
			}
		}
		labeling_time += (GetCurrentTimeSec() - labeling_acc_time);


		// Finish the constructions.
		for (NodeID v = 0; v < numOfVertices; ++v) {
			if (DIRECTED_FLAG == false) {
				NodeID k = tmp_idx[v].first.size();
				plabels.index_[v].spt_v.resize(k);
				plabels.index_[v].spt_d.resize(k);
				plabels.index_[v].spt_p.resize(k);
				for (NodeID i = 0; i < k; ++i) plabels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) plabels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				for (NodeID i = 0; i < k; ++i) plabels.index_[v].spt_p[i] = tmp_idx_parent[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();
				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();
				tmp_idx_parent[v].first.clear();
				tmp_idx_parent[v].second.clear();
				tmp_idx_parent[v].first.shrink_to_fit();
				tmp_idx_parent[v].second.shrink_to_fit();
			}
			else {
				NodeID k = tmp_idx[v].first.size();
				dplabels.index_[v].spt_v.resize(k);
				dplabels.index_[v].spt_d.resize(k);
				dplabels.index_[v].spt_p.resize(k);
				for (NodeID i = 0; i < k; ++i) dplabels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dplabels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				for (NodeID i = 0; i < k; ++i) dplabels.index_[v].spt_p[i] = tmp_idx_parent[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();
				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();
				tmp_idx_parent[v].first.clear();
				tmp_idx_parent[v].second.clear();
				tmp_idx_parent[v].first.shrink_to_fit();
				tmp_idx_parent[v].second.shrink_to_fit();

				k = r_tmp_idx[v].first.size();
				dplabels.bindex_[v].spt_v.resize(k);
				dplabels.bindex_[v].spt_p.resize(k);
				dplabels.bindex_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) dplabels.bindex_[v].spt_v[i] = r_tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dplabels.bindex_[v].spt_d[i] = r_tmp_idx[v].second[i];
				for (NodeID i = 0; i < k; ++i) dplabels.bindex_[v].spt_p[i] = r_tmp_idx_parent[v].second[i];
				r_tmp_idx[v].first.clear();
				r_tmp_idx[v].second.clear();
				r_tmp_idx[v].first.shrink_to_fit();
				r_tmp_idx[v].second.shrink_to_fit();
				r_tmp_idx_parent[v].first.clear();
				r_tmp_idx_parent[v].second.clear();
				r_tmp_idx_parent[v].first.shrink_to_fit();
				r_tmp_idx_parent[v].second.shrink_to_fit();
			}
		}

		num_of_trees = st_roots.size();

		double small_tree_num = 0;

		for (int i = 0; i < num_of_trees; ++i) {
			if (trees_pointers[i] != -1) {
				small_tree_num++;
				//small_tree_num += strees[trees_pointers[i]].size();
			}
		}

		total_resources = small_tree_num;
	}
	
	Betweenness_Ordering(NodeID k, double beta, Graph& graph, NodeID stopv, bool path_flag) {

		double start_time = GetCurrentTimeSec();

		if (k < CNT) {
			cout << " Should pick at least 16 seed!" << endl;
			return;
		}
		iteration_time.resize(numOfVertices);

		updating_time = 0;
		labeling_time = 0;
		adding_time = 0;
		selecting_time = 0;
		init_time = 0;
		u1_time = 0;
		u2_time = 0;
		u3_time = 0;
		num_of_trees = 0;
		total_resources = 0;
		get_small_parent_time = 0;

		long long bound_times = 10;
		long long trans_times = 8;

		inv.resize(numOfVertices, -1);
		rank.resize(numOfVertices, -1);

		if (DIRECTED_FLAG == false) {
			labels.index_.resize(numOfVertices);
		}
		else { //if (DIRECTED_FLAG == true) {
			dlabels.index_.resize(numOfVertices);
			dlabels.bindex_.resize(numOfVertices);
		}

		// Resource Bounding.
		NodeID new_tree_arcs = 0;
		NodeID labeled_arcs = 0;
		alive_resources = 0;
		bounded_resources = (long long)10 * (long long)k * (long long)numOfVertices;

		bool stop_flag = false;

		vector<NodeID> st_roots; // The root vertices of the sampled trees.
		vector<NodeID> st_sizes; // The size of the sampled trees.

		NodeID size_thr = numOfVertices < 8 ? numOfVertices : numOfVertices / 8; // When a tree shrink to a size smaller than the size_thr, it should be converted into a small_tree representation.
																			 //	vector<vector<int> > cover;


		init_time = GetCurrentTimeSec();
		// Building the initial k trees.
		//	srand((unsigned)time(NULL));
		srand(0x728f9d0);


		vector<NodeID> descendants_parents(numOfVertices, -1);

		vector<NodeID> parents(numOfVertices, numOfVertices);

		ltrees.reserve(numOfVertices);
		strees.reserve(numOfVertices);
		trees_pointers.reserve(numOfVertices);
		st_roots.reserve(numOfVertices);
		st_sizes.reserve(numOfVertices);
		st_roots.resize(k);
		st_sizes.resize(k);

		ltrees.resize(k);
		strees.clear();
		trees_pointers.resize(k, -1); // Point the sampling trees to its small tree representations.
									  //	cover.resize(k, vector<int>(V, 0));

		vector<NodeID> ascendants;
		ascendants.reserve(numOfVertices);
		vector<NodeID> descendants;
		descendants.reserve(numOfVertices);
		vector<NodeID> que(numOfVertices);
		vector<bool> vis(numOfVertices);
		vector<bool> has_picked(numOfVertices);
		vector<vector<cover_value_type> > cnt16(numOfVertices, vector<cover_value_type>(CNT));
		//		vector<bool> affected(V, true);
		vector<bool> picked_as_root(numOfVertices, false);
		vector<bool> trivial_root(numOfVertices, false);
		int picked_count = 0;
		int used_as_cover = 0;
		vector<cover_value_type> sum(numOfVertices);
		vector<bool> usd(numOfVertices, false);
		vector<cover_value_type> cover_tid(numOfVertices, 0);

		//unordered_set<NodeID> affected_v; //GOOGLE
		google::dense_hash_set<NodeID> affected_v;
		affected_v.set_empty_key(numOfVertices + 1);
		affected_v.set_deleted_key(numOfVertices + 2);

		vector<bool> affected_list(numOfVertices, false);
		int affected_count = 0;
		bool switch_set = false;

		//vector<vector<int> > cover_record(V, vector<int>(V, 0));

		//vector<unordered_set<int> > reverse_index(0); //GOOGLE
		//		vector<unordered_set<int> > reverse_index(0);
		vector<google::dense_hash_set<NodeID> > reverse_index;
		reverse_index.reserve(numOfVertices);

		// Preparing basic structure for pl algorithm.
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT))); // Backward labels.


		vector<pair<vector<NodeID>, vector<NodeID> > >
			tmp_idx_parent(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<NodeID>(1, numOfVertices)));

		vector<pair<vector<NodeID>, vector<NodeID> > >
			r_tmp_idx_parent(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<NodeID>(1, numOfVertices))); // Backward labels.

													 //	int best_sum = 0;
													 //	NodeID best_sum_vid = 0;
													 //vector<bool> lazy_updates(V, false);

													 // Building the initial k sampled trees. 
		for (int i = 0; i < k; ++i) {
			NodeID r = rand() % numOfVertices;
			if (picked_as_root[r] == true) { --i; continue; }

			st_roots[i] = r;
			large_tree &ltree_i = ltrees[i];
			picked_as_root[r] = true;
			picked_count++;
			used_as_cover++;

			// Array Representation
			//build_large_tree_bfs(r, ltree_i, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph.adj, 0, numOfVertices);
			build_large_tree_bfs(r, ltree_i, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph, 0, numOfVertices);

			st_sizes[i] = descendants.size();
			alive_resources += (long long)numOfVertices;

			// Turn large tree to small tree.
			// 1. Build new small tree with ltree_i and add it to strees;
			// 2. Set tree_pointers[i] to the position of the strees;
			// 3. Clear the large tree ltree_i by setting its size to 0.
			if (st_sizes[i] < size_thr)
				shrink_tree(i, st_roots[i], ltrees[i]);


			// Calculate the subtree size of each nodes			.
			//	vector<int>& cover_i = cover[i];

			bool isLarge = false;
			int stree_idx = trees_pointers[i];
			if (stree_idx == -1)
				isLarge = true;

			small_tree& stree_i = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;

			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];

				++cover_tid[q]; // Increment because it hitsthe r-q path.
				++cnt16[q][i%CNT];
				//++sum[q];
				//	if (sum[q] > best_sum) {
				//		best_sum = sum[q];
				//		best_sum_vid = q;
				//	}
				////	cover_record[i][q]++;

				int parent_node;
				if (isLarge)
					parent_node = get_large_parent(q, ltree_i);
				else
					parent_node = get_small_parent(q, stree_i);


				if (parent_node > -1) {
					cover_tid[parent_node] += cover_tid[q];
					cnt16[parent_node][i%CNT] += cover_tid[q];
					//sum[parent_node] += cover_i[q];
					//	cover_record[i][parent_node] += cover_record[i][q];
				}
			}
			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];
				cover_tid[q] = 0;
			}
		}


		//cnt16_aggreate(cnt16, has_picked, cover);

		vector<cover_value_type> picked_value;
		benchmark::heap<2, cover_value_type, NodeID> value_heap(numOfVertices);
		bool switch_small = false;

		//cnt16_2max(picked_value, sum, has_picked,affected, cnt16, switch_small, value_heap);
		//		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap, best_sum_vid, lazy_updates);
		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);

		init_time = GetCurrentTimeSec() - init_time;


		cout << "Finish initial sampling " << init_time << "secs." << endl;

		//for (NodeID v = 0; v < V; ++v)	value_heap.update(v, -picked_value[v]); // Min-heap.

		double max_value = -1;


		// Perform V iterations to rank all vertices.
		NodeID ranking = 0;
		for (; ranking < numOfVertices; ++ranking) {
			/*	if (ranking % (numOfVertices / 10) == 0)
			cout << "picking the " << ranking << "-th orders." << endl;*/

			double time_this_iteration = GetCurrentTimeSec();

			//best_sum = 0;

			NodeID pv = 0;
			cover_value_type pv_cover = -1;

			if (ranking == stopv) {
				stop_flag = true;
				cout << "topv labeling time:" << labeling_time << endl;
				break;
			}
			double selecting_acc_time = GetCurrentTimeSec();
			// In the beginning, pick the best cover vertex by a linear scan.
			if (switch_small == false) {
				double u1_acc_time = GetCurrentTimeSec();
				//while (true) {
				for (NodeID v = 0; v < numOfVertices; ++v) {
					if (has_picked[v] == true) continue;
					if (pv_cover == picked_value[v]) {
						if (sum[pv] < sum[v]) {
							pv_cover = picked_value[v];
							pv = v;
						}
					}
					else if (pv_cover < picked_value[v]) {
						pv_cover = picked_value[v];
						pv = v;
					}
				}
				
				if (sum[pv] < numOfVertices / 8 && switch_small == false) {
					cout << ranking << endl;
					cout << pv_cover << endl;
					cout << sum[pv] << endl;
					cout << picked_count << endl;
					cout << "Switch to max-heap representations." << endl;
					switch_small = true;
				}
				u1_time += (GetCurrentTimeSec() - u1_acc_time);
			}
			else {// Later, utilizing the max-heap directly
				double u2_acc_time = GetCurrentTimeSec();
				//while (true) {
				value_heap.extract_min(pv, pv_cover);
				pv_cover = -pv_cover;

				u2_time += (GetCurrentTimeSec() - u2_acc_time);

			}

			selecting_time += (GetCurrentTimeSec() - selecting_acc_time);

			has_picked[pv] = true;
			double labeling_acc_time = GetCurrentTimeSec();
			labeled_arcs = labeling_source_bfs_path(pv, que, vis, parents, dst_r, usd, ranking, tmp_idx, r_tmp_idx, tmp_idx_parent, r_tmp_idx_parent,graph);
			labeling_time += (GetCurrentTimeSec() - labeling_acc_time);
			inv[ranking] = pv;
			rank[pv] = ranking;

			//cout << ranking << "\t" << labeled_arcs << endl;


			int current_trees_num = st_roots.size();

			// 4 steps in updating:
			double updating_acc_time = GetCurrentTimeSec();
			if (switch_small == false) { // We still iterate all alive trees to calculate the cover.	
				for (int tid = 0; tid < current_trees_num; ++tid) {

					// This tree has already be taken because its root has been selected.
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					// pv has already been covered in tree tid.
					int pvpid = 0;
					if (isLarge)
						pvpid = get_large_parent(pv, ltree_tid);
					else
						pvpid = get_small_parent(pv, stree_tid);
					if (pvpid == -1 && st_roots[tid] != pv) // pv has been covered in tree tid 
						continue;


					// Get the descendants and ascendants of pv in tree tid.

					// Array Representation
					//// 1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;				 
					//get_descendants(tid, pv, descendants, descendants_parents, que, vis, graph.adj, isLarge, stree_tid, ltree_tid);

					//get_ascendants(tid, pv, ascendants, que, vis, graph.adj, graph.r_adj, isLarge, stree_tid, ltree_tid);

					get_descendants(tid, pv, descendants, descendants_parents, que, vis, graph, isLarge, stree_tid, ltree_tid);

					get_ascendants(tid, pv, ascendants, que, vis, graph, isLarge, stree_tid, ltree_tid);

					NodeID subtract_size = descendants.size();
					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
												//cover[tid][av] -= subtract_size;
						cnt16[av][tid%CNT] -= subtract_size;
						//sum[av] -= subtract_size;
						//	cover_record[tid][av] -= subtract_size;
						//	if (cover_record[tid][av] < 0)
						//	cout << "1: " << cover_record[tid][av] << endl;
						if (affected_list[av] == false) {
							affected_list[av] = true;
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
						}
					}


					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					//	vector<int> cover_tid(numOfVertices, 0);
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];

						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//sum[dv] -= cover_tid[dv];
						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];

						/*	if (cnt16[dv][tid%CNT] < 0)
						cout << "a6 " << cnt16[dv][tid%CNT] << endl;*/

						//	cover_record[tid][dv] -= cover_tid[dv];
						//	if (cover_record[tid][dv] < 0) {
						//	cout << "2: " << cover_record[tid][dv]  << " " << cover_tid[dv]  << " " << tid << " " << dv <<  endl;
						//	}

						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else
							remove_small_node(dv, stree_tid);


						//	remove_node(tid, dv);
						//affected[dv] = true;
						if (affected_list[dv] == false) {
							affected_count++;
							affected_list[dv] = true;
							if (switch_set)
								affected_v.insert(dv);
						}
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;
					}


					// Drop the size of subtree.
					st_sizes[tid] -= subtract_size;
					if (trees_pointers[tid] != -1) alive_resources -= (cover_value_type)subtract_size;
					//	if (st_sizes[tid] < 0) {
					//		cout << "a5" << st_sizes[tid] << endl;
					//	}
					if (st_sizes[tid] > 0) {
						// 3. Convert shrink trees to small trees representataions.				
						if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
							shrink_tree(tid, st_roots[tid], ltrees[tid]);
						}
					}

					if (st_roots[tid] == pv)
						if (isLarge) {
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
							//ltrees.shrink_to_fit();
						}
						else {
							stree_tid.clear();
							stree_tid.resize(0);
							//strees.shrink_to_fit();
						}

				}
			}
			else { // We use the reverse index to update trees.

				   // Build the reverse_index first (one time).
				if (reverse_index.size() == 0) {
					reverse_index.resize(numOfVertices);
					for (NodeID v = 0; v < numOfVertices; ++v) {
						reverse_index[v].set_empty_key(numOfVertices + 1);
						reverse_index[v].set_deleted_key(numOfVertices + 2);
					}
					for (int tid = 0; tid < current_trees_num; ++tid) {

						if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue; //This whole tree has already been covered. 

						bool isLarge = false;
						int stree_idx = trees_pointers[tid];
						if (stree_idx == -1)
							isLarge = true;

						small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
						large_tree& ltree_tid = ltrees[tid];

						for (int v = 0; v < numOfVertices; ++v) {

							// Get the descendants and ascendants of pv in tree tid.

							//	if (cover[tid].find(pv) == cover[tid].end()) continue; // This pv has been covered in this tree before.
							//	if (cover[tid][pv] <= 0) continue; // modified-5-14
							int vpid = 0;
							if (isLarge)
								vpid = get_large_parent(v, ltree_tid);
							else
								vpid = get_small_parent(v, stree_tid);

							// To testify whether v is still uncovered in tree tid
							// v has parent vpid != -1, but if v is the root of tree tid, it still can be uncovered.
							if ((vpid != -1 || st_roots[tid] == v))
								reverse_index[v].insert(tid);
							//if (cover[tid].find(v) != cover[tid].end())
							//	reverse_index[v].insert(tid);
						}
					}
				}

				double u4_acc_time = GetCurrentTimeSec();

				//	double u4_acc_time = GetCurrentTimeSec();
				// Access the alive trees directly by the reverse index.
				//				cout << reverse_index[pv].size() << " " << sum[pv] << endl;


				//				for (unordered_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) { // GOOGLE
				for (google::dense_hash_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) {


					int tid = *it;
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					// Get the descendants and ascendants of pv in tree tid.

					//	if (cover[tid].find(pv) == cover[tid].end()) continue; // This pv has been covered in this tree before.
					//	if (cover[tid][pv] <= 0) continue; // modified-5-14

					int pvpid = 0;
					if (isLarge)
						pvpid = get_large_parent(pv, ltree_tid);
					else {
						//	double get_small_parent_incre_time = GetCurrentTimeSec();
						pvpid = get_small_parent(pv, stree_tid);
						//	get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}

					if (pvpid == -1 && st_roots[tid] != pv)
						continue;

					// Array Representation
					//// 1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;
					//get_descendants(tid, pv, descendants, descendants_parents, que, vis, graph.adj, isLarge, stree_tid, ltree_tid);

					//get_ascendants(tid, pv, ascendants, que, vis, graph.adj, graph.r_adj, isLarge, stree_tid, ltree_tid);

					get_descendants(tid, pv, descendants, descendants_parents, que, vis, graph, isLarge, stree_tid, ltree_tid);

					get_ascendants(tid, pv, ascendants, que, vis, graph, isLarge, stree_tid, ltree_tid);

					NodeID subtract_size = descendants.size();

					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
												//						cover[tid][av] -= subtract_size;

						cnt16[av][tid%CNT] -= subtract_size;
						//sum[av] -= subtract_size;

						//	cover_record[tid][av] -= subtract_size;
						//	if (cover_record[tid][av] < 0) {
						//	cout << "4: " << cover_record[tid][av] << " " << subtract_size << " " << tid << " " << endl;
						//	}

						//affected[av] = true;
						/*if (cnt16[av][tid%CNT] < 0)
						cout << " a2 " << cnt16[av][tid%CNT]  << " " << av << " " << st_roots[tid] << " " << tid << endl;*/
						if (affected_list[av] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
							affected_list[av] = true;
						}
					}

					//double u2_acc_time = GetCurrentTimeSec();
					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					//vector<int> cover_tid(numOfVertices, 0);
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//	sum[dv] -= cover_tid[dv];
						//cover[tid].erase(dv);
						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];

						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else
							remove_small_node(dv, stree_tid);

						//	remove_node(tid, dv);
						//affected[dv] = true;

						//	cover_record[tid][dv] -= cover_tid[dv];
						//if (cover_record[tid][dv] < 0)
						//	cout << "3: " << cover_record[tid][dv] << " " << cover_tid[dv] << " " << di << " " << descendants.size() << " " << picked_as_root[st_roots[tid]] << endl;

						/*if (cnt16[dv][tid%CNT]  < 0)
						cout << " a3 " << cnt16[dv][tid%CNT] << " " << dv << " " << st_roots[tid] << " " << tid << endl;
						*/
						if (dv != pv)
							reverse_index[dv].erase(tid);
						//in_tree_count[dv]--;


						if (affected_list[dv] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(dv);
							affected_list[dv] = true;
						}
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;
					}

					if (st_roots[tid] == pv)
						if (isLarge) {
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
						}
						else {
							stree_tid.clear();
							stree_tid.resize(0);
						}

						st_sizes[tid] -= subtract_size;
						//	if (st_sizes[tid]  < 0)
						//		cout << " a4 " << st_sizes[tid] << endl;

						if (trees_pointers[tid] != -1) alive_resources -= subtract_size;


						// 3. Convert shrink trees to small trees representataions.	
						if (st_sizes[tid] > 0) {
							if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
								shrink_tree(tid, st_roots[tid], ltrees[tid]);
							}
						}

				}
				reverse_index[pv].clear();
				reverse_index[pv].resize(0);
				u4_time += (GetCurrentTimeSec() - u4_acc_time);
			}

			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			// 4. Sample new trees based on.
			// If we sample a new tree, we need to dynamically increment the following data structures:
			// 1) ltrees
			// 2) trees_pointers (-1 or pointer to strees)
			// 3) cover
			// 4) st_roots
			// 5) st_sizes 
			// 6) strees (maybe)
			// two constraints: cl > ct and 10kn total vertices
			double adding_acc_time = GetCurrentTimeSec();
			new_tree_arcs = 0;

			//cout << "ranking:" << ranking << endl;
			//cout << "pick cover:" << pv_cover << endl;
			//cout << "pick count:" << picked_count << endl;
			//cout << "labeled_arcs:" << labeled_arcs << endl;
			//cout << "new_tree_arcs:" << new_tree_arcs << endl;
			//cout << "alive_resources:" << alive_resources << endl;
			//cout << "bounded_resources:" << bounded_resources << endl;
			//cout << "#tree:" << cover.size() << endl;
			int how_many = 0;
			//cout << "beta*labeld_arcs: " << beta * labeled_arcs << " new tree _arcs: " << new_tree_arcs << " alive " << alive_resources << " bounded " << bounded_resources << endl;
			NodeID sam_this_round = 0;
			while (beta * labeled_arcs > new_tree_arcs && alive_resources < bounded_resources) {

				if (picked_count == numOfVertices) break;



				NodeID r = rand() % numOfVertices;
				if (picked_as_root[r] == true) continue;
				if (trivial_root[r] == true) continue;
				int tmp_tree_arcs = 0;

				large_tree new_ltree;
				small_tree new_stree;
				//GOOGLE DENSE MAP
				new_stree.set_empty_key(numOfVertices + 1);
				new_stree.set_deleted_key(numOfVertices + 2);

				//	BUILD_SMALL_TREE_FLAG = false;
				bool isLarge = false;
				if (!BUILD_SMALL_TREE_FLAG) {

					double u5_acc_time = GetCurrentTimeSec();
					// Array Representation
					//tmp_tree_arcs = build_large_tree_bfs(r, new_ltree, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph.adj, new_tree_arcs, beta * labeled_arcs);
					tmp_tree_arcs = build_large_tree_bfs(r, new_ltree, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph, new_tree_arcs, beta * labeled_arcs);
					u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;

					if (descendants.size() <size_thr)
						BUILD_SMALL_TREE_FLAG = true;

					if (descendants.size() < 2) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					isLarge = true;
					how_many++;
					used_as_cover++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					st_roots.push_back(r);
					trees_pointers.push_back(-1);
					st_sizes.push_back(descendants.size());

					//new_tree_arcs += descendants.size();
					alive_resources += numOfVertices;

					// Turn large tree to small tree.
					// 1. Build new small tree with ltree_i and add it to strees;
					// 2. Set tree_pointers[i] to the position of the strees;
					// 3. Clear the large tree ltree_i by setting its size to 0.
					int i = st_roots.size() - 1;
					if (st_sizes[i] < size_thr) {
						shrink_tree(i, st_roots[i], new_ltree);
						isLarge = false;
					}
					ltrees.push_back(new_ltree);
				}
				else {//Build Small Tree Directly
					double u5_acc_time = GetCurrentTimeSec();
					// Array Representation
					//	tmp_tree_arcs = build_small_tree_bfs(r, new_stree, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph.adj, new_tree_arcs, beta * labeled_arcs);
					tmp_tree_arcs = build_small_tree_bfs(r, new_stree, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph, new_tree_arcs, beta * labeled_arcs);
					u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;

					if (descendants.size() >= size_thr)
						BUILD_SMALL_TREE_FLAG = false;

					if (descendants.size() < 2) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					isLarge = false;

					how_many++;
					used_as_cover++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					st_roots.push_back(r);
					trees_pointers.push_back(strees.size()); //we assume it is inserted.
					st_sizes.push_back(descendants.size());

					alive_resources += new_stree.size();
					//strees.push_back(new_stree);
					//// new_ltree is empty.
					//ltrees.push_back(new_ltree);

					int i = st_roots.size() - 1;
					if (st_sizes[i] >= size_thr) {
						enlarge_tree(i, st_roots[i], new_stree, new_ltree);
						new_stree.clear();
						new_stree.resize(0);
						isLarge = true;
					}
					else {
						strees.push_back(new_stree);
						ltrees.push_back(new_ltree);
					}
				}

				int i = st_roots.size() - 1;

				int numOfStrees = strees.size();
				small_tree &new_stree_ref = numOfStrees != 0 ? strees[strees.size() - 1] : empty_small_tree;

				// Calculate the subtree size of each nodes.
				//				unordered_map<int, int> cover_new;
				//	vector<cover_value_type> cover_new(numOfVertices, 0); // modified-5-14

				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];

					++cover_tid[q]; // Increment because it hitsthe r-q path. 
					++cnt16[q][i%CNT];
					//++sum[q];
					/*	if (sum[q] > best_sum) {
					best_sum = sum[q];
					best_sum_vid = q;
					}*/
					//affected[q] = true;
					//	cover_record[i][q]++;
					if (affected_list[q] == false) {
						affected_count++;
						if (switch_set)
							affected_v.insert(q);
						affected_list[q] = true;
					}

					if (switch_small == true) {
						reverse_index[q].insert(i);
					}

					int parent_node = -1;

					if (isLarge)
						parent_node = get_large_parent(q, new_ltree);
					else
						parent_node = get_small_parent(q, new_stree_ref);

					if (parent_node > -1) {
						cover_tid[parent_node] += cover_tid[q];
						cnt16[parent_node][i%CNT] += cover_tid[q];
						//	sum[parent_node] += cover_new[q];
						//	cover_record[i][parent_node] += cover_record[i][q];
					}
				}
				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];
					cover_tid[q] = 0;
				}

				sam_this_round++;
			}

			//cout << ranking << "\t" << sam_this_round << endl;


			adding_time += (GetCurrentTimeSec() - adding_acc_time);
			//cout << ranking << ":  beta:" << beta << " beta * labeled_arcs " << beta * labeled_arcs << " new_tree_arcs " << new_tree_arcs << " how many " << how_many << endl;

			//double u3_acc_time = GetCurrentTimeSec();

			updating_acc_time = GetCurrentTimeSec();

			// Build the cover value heap for one time.
			if (switch_small == true && value_heap.empty())
				for (int nid = 0; nid < numOfVertices; ++nid)
					if (has_picked[nid] == false)
						value_heap.update(nid, -picked_value[nid]);

			//	cnt16_2max(picked_value, sum, has_picked, affected, cnt16, switch_small, value_heap);		
			//double u3_acc_time = GetCurrentTimeSec();

			//			cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap, best_sum_vid, lazy_updates);
			cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);
			//u3_time += (GetCurrentTimeSec() - u3_acc_time);

			if (affected_count < numOfVertices / 8)
				switch_set = true;
			else
				switch_set = false;
			affected_count = 0;

			//iteration_time[ranking] = (GetCurrentTimeSec() - updating_acc_time);
			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			time_this_iteration = GetCurrentTimeSec() - time_this_iteration;
			iteration_time[ranking] = time_this_iteration;
		}


		cout << "picked count:" << picked_count << endl;
		cout << "picked as tree covers:" << used_as_cover << endl;
		/*for (int v = 0; v < V; ++v)
		if (has_picked[v] == false) {
		cout << "wrong 2" << endl;
		cout << picked_value[v] << endl;
		}*/

		//// FIX
		//for (int v = 0; v < numOfVertices; ++v){
		//	for (int i = 0; i < CNT; ++i) {
		//		if (cnt16[v][i] != 0)
		//			cout << "mi?" << endl;
		//	}
		//}

		if (stop_flag == true) return;

		double labeling_acc_time = GetCurrentTimeSec();
		// Order the remaining vertices with degree, if any.
		if (ranking != numOfVertices) {
			vector<NodeID> degree(numOfVertices);
			vector<NodeID> r_degree(numOfVertices);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				/*vector<NodeID>& adj_v = graph.adj[v];
				for (NodeID i = 0; i < adj_v.size(); ++i) {
				if (has_picked[adj_v[i]] != true)
				degree[v]++;
				}*/

				for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid) {
					NodeID w = graph.edges[eid];
					if (has_picked[w] != true)
						degree[v]++;
				}

				
			}

			vector<pair<NodeID, NodeID> > degree_v;
			degree_v.reserve(numOfVertices - ranking);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				if (DIRECTED_FLAG == true)
					degree_v.push_back(make_pair(degree[v] * r_degree[v], v));
				else
					degree_v.push_back(make_pair(degree[v], v));
			}
			sort(degree_v.rbegin(), degree_v.rend());

			int count = 0;
			for (; ranking < numOfVertices; ++ranking) {
				NodeID tpv = degree_v[count++].second;
				inv[ranking] = tpv;
				rank[tpv] = ranking;
				labeling_source_bfs(tpv, que, vis, dst_r, usd, ranking, tmp_idx, r_tmp_idx, graph);
			}
		}
		labeling_time += (GetCurrentTimeSec() - labeling_acc_time);

		// Finish the constructions.
		for (NodeID v = 0; v < numOfVertices; ++v) {
			if (DIRECTED_FLAG == false) {
				NodeID k = tmp_idx[v].first.size();
				plabels.index_[v].spt_v.resize(k);
				plabels.index_[v].spt_d.resize(k);
				plabels.index_[v].spt_p.resize(k);
				for (NodeID i = 0; i < k; ++i) plabels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) plabels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				for (NodeID i = 0; i < k; ++i) plabels.index_[v].spt_p[i] = tmp_idx_parent[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();
				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();
				tmp_idx_parent[v].first.clear();
				tmp_idx_parent[v].second.clear();
				tmp_idx_parent[v].first.shrink_to_fit();
				tmp_idx_parent[v].second.shrink_to_fit();
			}
			else {
				NodeID k = tmp_idx[v].first.size();
				dplabels.index_[v].spt_v.resize(k);
				dplabels.index_[v].spt_d.resize(k);
				dplabels.index_[v].spt_p.resize(k);
				for (NodeID i = 0; i < k; ++i) dplabels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dplabels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				for (NodeID i = 0; i < k; ++i) dplabels.index_[v].spt_p[i] = tmp_idx_parent[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();
				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();
				tmp_idx_parent[v].first.clear();
				tmp_idx_parent[v].second.clear();
				tmp_idx_parent[v].first.shrink_to_fit();
				tmp_idx_parent[v].second.shrink_to_fit();

				k = r_tmp_idx[v].first.size();
				dplabels.bindex_[v].spt_v.resize(k);
				dplabels.bindex_[v].spt_p.resize(k);
				dplabels.bindex_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) dplabels.bindex_[v].spt_v[i] = r_tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dplabels.bindex_[v].spt_d[i] = r_tmp_idx[v].second[i];
				for (NodeID i = 0; i < k; ++i) dplabels.bindex_[v].spt_p[i] = r_tmp_idx_parent[v].second[i];
				r_tmp_idx[v].first.clear();
				r_tmp_idx[v].second.clear();
				r_tmp_idx[v].first.shrink_to_fit();
				r_tmp_idx[v].second.shrink_to_fit();
				r_tmp_idx_parent[v].first.clear();
				r_tmp_idx_parent[v].second.clear();
				r_tmp_idx_parent[v].first.shrink_to_fit();
				r_tmp_idx_parent[v].second.shrink_to_fit();
			}
		}
        long long totalNum = 0;
    	for (NodeID v = 0; v < numOfVertices; ++v) { 
    	    totalNum += plabels.index_[v].spt_v.size();
        }
        cout <<"avg label: "<< ((double)totalNum)/numOfVertices -1.0 << endl;
	/*	for (NodeID v = 0; v < numOfVertices; ++v) {
			if (DIRECTED_FLAG == false) {
				NodeID k = tmp_idx[v].first.size();
				labels.index_[v].spt_v.resize(k);
				labels.index_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();

				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();

			}
			else {
				NodeID k = tmp_idx[v].first.size();
				dlabels.index_[v].spt_v.resize(k);
				dlabels.index_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();

				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();

				k = r_tmp_idx[v].first.size();
				dlabels.bindex_[v].spt_v.resize(k);
				dlabels.bindex_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_v[i] = r_tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_d[i] = r_tmp_idx[v].second[i];
				r_tmp_idx[v].first.clear();
				r_tmp_idx[v].second.clear();

				r_tmp_idx[v].first.shrink_to_fit();
				r_tmp_idx[v].second.shrink_to_fit();
			}
		}*/

		num_of_trees = st_roots.size();

		double small_tree_num = 0;

		for (int i = 0; i < num_of_trees; ++i) {
			if (trees_pointers[i] != -1) {
				small_tree_num++;
				//small_tree_num += strees[trees_pointers[i]].size();
			}
		}

		total_resources = small_tree_num;
		/*for (int i = 0; i < cnt16.size(); ++i){
		for (int j = 0; j < cnt16[i].size(); ++j)
		if (cnt16[i][j] != 0)
		cout << "aaaa " << cnt16[i][j] << endl;
		}*/

	}

	~Betweenness_Ordering() {
		clear();
	}

	void clear() {
		ltrees.clear();
		strees.clear();
		trees_pointers.clear();
		labels.Free();
		dlabels.Free();
	}

};

template<int kNumBitParallelRoots = 50>
class BP_Betweenness_Ordering :public Ordering {

	typedef	vector<NodeID> large_tree; // A V-sized parent-pointer array representing sampled trees. -1 for those vertices which do not appear in the tree.
									   //typedef	unordered_map<NodeID, NodeID> small_tree; // A smaller hashmap, mapping a tree node vertex to its parent vertex.
	typedef	google::dense_hash_map<NodeID, NodeID> small_tree; // A smaller hashmap, mapping a tree node vertex to its parent vertex.

public:
	vector<double> iteration_time;
	Label labels;
	DLabel dlabels;
	PLabel plabels;
	DPLabel dplabels;
	BPLabel<kNumBitParallelRoots> bplabels;
	DBPLabel<kNumBitParallelRoots> dbplabels;

	vector<large_tree> ltrees;
	vector<small_tree> strees;
	vector<NodeID> trees_pointers;
	double init_time;
	double updating_time;
	double labeling_time;
	double adding_time;
	double selecting_time;
	long long bounded_resources;
	int num_of_trees;
	long total_resources;
	long long alive_resources;
	double u1_time;
	double u2_time;
	double u3_time = 0;
	double u4_time = 0;
	double u5_time = 0;
	double u6_time = 0;
	double u7_time = 0;
	double u8_time = 0;
	double u9_time = 0;
	double u10_time = 0;
	double u11_time = 0;
	double u12_time = 0;
	double u81_time = 0;
	double u82_time = 0;
	double u91_time = 0;
	double u92_time = 0;

	int skip_count = 0;

	small_tree empty_small_tree = small_tree();
	bool BUILD_SMALL_TREE_FLAG = false;
	double get_small_parent_time;
	bool switch_small;

	int build_large_tree_bfs(NodeID source, vector<int>& ltree, vector<NodeID>& descendants, vector<NodeID>& descendants_parents, vector<NodeID>& que, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<bool>& usd, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& r_tmp_idx, /*vector<vector<NodeID> >& adj*/ Graph& graph, int accum_new_arcs, int labeled_arcs_bound, NodeID ranking) {
		NodeID new_tree_arcs = 0;

		descendants.clear();
		ltree.resize(numOfVertices, -1);
		//ltree[source] = -1; // As the root, it has no parent.
		descendants_parents[source] = -1;

		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
		const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];
		const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[source];
		index_t_bp<kNumBitParallelRoots> &idx_r = (DIRECTED_FLAG == true) ? dbplabels.index_bp[source] : bplabels.index_bp[source];



		que[que_h++] = source;
		vis[source] = true;
		que_t1 = que_h;
		++new_tree_arcs;

		int currentBP = ranking - skip_count;
		if (currentBP > kNumBitParallelRoots + skip_count)
			currentBP = kNumBitParallelRoots;

		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
		}

		for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID v = que[que_i];

				if (usd[v] == true) continue;

				//	int adjv_size = adj[v].size();

				if (DIRECTED_FLAG == false) {

					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
					index_t_bp<kNumBitParallelRoots> &idx_v =  bplabels.index_bp[v];


					// Prefetch
					_mm_prefetch(&idx_v.bpspt_d[0], _MM_HINT_T0);
					_mm_prefetch(&idx_v.bpspt_s[0][0], _MM_HINT_T0);
					_mm_prefetch(&tmp_idx_v.first[0], _MM_HINT_T0);
					_mm_prefetch(&tmp_idx_v.second[0], _MM_HINT_T0);


					for (int i = 0; i < currentBP; ++i) {
						EdgeWeight td = idx_r.bpspt_d[i] + idx_v.bpspt_d[i];
						if (td - 2 <= d) {
							td +=
								(idx_r.bpspt_s[i][0] & idx_v.bpspt_s[i][0]) ? -2 :
								((idx_r.bpspt_s[i][0] & idx_v.bpspt_s[i][1]) |
									(idx_r.bpspt_s[i][1] & idx_v.bpspt_s[i][0]))
								? -1 : 0;
							if (td <= d) goto pruned;
						}
					}

					int tmp_idx_v_size = tmp_idx_v.first.size();
					for (int i = 0; i < tmp_idx_v_size; ++i) {
						NodeID w = tmp_idx_v.first[i];
						EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
						if (td <= d) {
							// num of hubs touched so far when building a new tree.
							goto pruned;
						}
					}
				} 
				else {
					pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];
					index_t_bp<kNumBitParallelRoots> &r_idx_v = dbplabels.bindex_bp[v];

					for (int i = 0; i < currentBP; ++i) {
						EdgeWeight td = idx_r.bpspt_d[i] + r_idx_v.bpspt_d[i];
						if (td - 2 <= d) {
							td +=
								(idx_r.bpspt_s[i][0] & r_idx_v.bpspt_s[i][0]) ? -2 :
								((idx_r.bpspt_s[i][0] & r_idx_v.bpspt_s[i][1]) |
									(idx_r.bpspt_s[i][1] & r_idx_v.bpspt_s[i][0]))
								? -1 : 0;
							if (td <= d) goto pruned;
						}
					}

					// Pruned by the forward labels of r and backward labels of v in the forward search from r when reaching v.
					int r_tmp_idx_v_size = r_tmp_idx_v.first.size();
					for (int i = 0; i < r_tmp_idx_v_size; ++i) {
						NodeID w = r_tmp_idx_v.first[i];
						EdgeWeight td = r_tmp_idx_v.second[i] + dst_r[w];
						if (td <= d) {
							// num of hubs touched so far when building a new tree.
							goto pruned;
						}
					}
				}

				descendants.push_back(v);
				++new_tree_arcs;
				if (accum_new_arcs + new_tree_arcs > labeled_arcs_bound) {
					goto stop_growing;
				}

				// Array Representation
				/*		for (int i = 0; i < adjv_size; ++i) {
				NodeID w = adj[v][i];*/
				for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid) {
					NodeID w = graph.edges[eid];

					if (!vis[w]) {
						que[que_h++] = w;
						descendants_parents[w] = v;
						//ltree[w] = v; // Set parent nodes. 
						vis[w] = true;		// num of hubs touched so far when building a new tree.					
					}
				}
			pruned:
				{}
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}

	stop_growing:
		{}

		for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
		}
		for (int i = 0; i < descendants.size(); ++i) {
			ltree[descendants[i]] = descendants_parents[descendants[i]];
		}

		return new_tree_arcs;
	}

	int build_small_tree_bfs(NodeID source, small_tree& stree, vector<NodeID>& descendants, vector<NodeID>& descendants_parents, vector<NodeID>& que, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<bool>& usd, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& r_tmp_idx, /*vector<vector<NodeID> >& adj*/ Graph& graph, int accum_new_arcs, int labeled_arcs_bound, NodeID ranking) {
		NodeID new_tree_arcs = 0;

		descendants.clear();
		stree.clear();
		//	stree.reserve(V / 8);
		descendants_parents[source] = -1;
		//stree.insert(make_pair(source, -1));// As the root, it has no parent.

		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
		const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];
		const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[source];
		index_t_bp<kNumBitParallelRoots> &idx_r = (DIRECTED_FLAG == true)? dbplabels.index_bp[source] : bplabels.index_bp[source];


		que[que_h++] = source;
		vis[source] = true;
		que_t1 = que_h;
		++new_tree_arcs;

		int currentBP = ranking - skip_count;
		if (currentBP > kNumBitParallelRoots + skip_count)
			currentBP = kNumBitParallelRoots;

		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
		}

		for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID v = que[que_i];

				if (usd[v] == true) continue;

				//int adjv_size = adj[v].size();
				if (DIRECTED_FLAG == false) {




					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
					index_t_bp<kNumBitParallelRoots> &idx_v = bplabels.index_bp[v];


					// Prefetch
					_mm_prefetch(&idx_v.bpspt_d[0], _MM_HINT_T0);
					_mm_prefetch(&idx_v.bpspt_s[0][0], _MM_HINT_T0);
					_mm_prefetch(&tmp_idx_v.first[0], _MM_HINT_T0);
					_mm_prefetch(&tmp_idx_v.second[0], _MM_HINT_T0);


					for (int i = 0; i < currentBP; ++i) {
						EdgeWeight td = idx_r.bpspt_d[i] + idx_v.bpspt_d[i];
						if (td - 2 <= d) {
							td +=
								(idx_r.bpspt_s[i][0] & idx_v.bpspt_s[i][0]) ? -2 :
								((idx_r.bpspt_s[i][0] & idx_v.bpspt_s[i][1]) |
									(idx_r.bpspt_s[i][1] & idx_v.bpspt_s[i][0]))
								? -1 : 0;
							if (td <= d) goto pruned;
						}
					}

					int tmp_idx_v_size = tmp_idx_v.first.size();
					for (int i = 0; i < tmp_idx_v_size; ++i) {
						NodeID w = tmp_idx_v.first[i];
						EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
						if (td <= d) {
							// num of hubs touched so far when building a new tree.
							goto pruned;
						}
					}
				}
				else {
					pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];
					index_t_bp<kNumBitParallelRoots> &r_idx_v = dbplabels.bindex_bp[v];

					for (int i = 0; i < currentBP; ++i) {
						EdgeWeight td = idx_r.bpspt_d[i] + r_idx_v.bpspt_d[i];
						if (td - 2 <= d) {
							td +=
								(idx_r.bpspt_s[i][0] & r_idx_v.bpspt_s[i][0]) ? -2 :
								((idx_r.bpspt_s[i][0] & r_idx_v.bpspt_s[i][1]) |
									(idx_r.bpspt_s[i][1] & r_idx_v.bpspt_s[i][0]))
								? -1 : 0;
							if (td <= d) goto pruned;
						}
					}

					// Pruned by the forward labels of r and backward labels of v in the forward search from r when reaching v.
					int r_tmp_idx_v_size = r_tmp_idx_v.first.size();
					for (int i = 0; i < r_tmp_idx_v_size; ++i) {
						NodeID w = r_tmp_idx_v.first[i];
						EdgeWeight td = r_tmp_idx_v.second[i] + dst_r[w];
						if (td <= d) {
							// num of hubs touched so far when building a new tree.
							goto pruned;
						}
					}
				}

				descendants.push_back(v);
				// num of hubs touched so far when building a new tree.
				++new_tree_arcs;
				if (accum_new_arcs + new_tree_arcs > labeled_arcs_bound) {
					goto stop_growing;
				}

				// Array Representation
				/*for (int i = 0; i < adjv_size; ++i) {
				NodeID w = adj[v][i];*/
				for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid) {
					NodeID w = graph.edges[eid];

					if (!vis[w]) {
						que[que_h++] = w;
						descendants_parents[w] = v;
						//stree.insert(make_pair(w, v));// Set parent nodes. 
						vis[w] = true;
					}
				}
			pruned:
				{}
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}

	stop_growing:
		{}

		for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
		}

		for (int i = 0; i < descendants.size(); ++i) {
			stree.insert(make_pair(descendants[i], descendants_parents[descendants[i]]));
		}

		return new_tree_arcs;
	}

	void get_descendants(int tree_id, NodeID picked_v, vector<NodeID>& descendants, vector<NodeID>& que, vector<bool>& vis, vector<vector<NodeID> >& adj) {
		descendants.clear();

		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;

		que[que_h++] = picked_v;
		vis[picked_v] = true;
		que_t1 = que_h;


		bool isLarge = false;
		int stree_idx = trees_pointers[tree_id];
		if (stree_idx == -1)
			isLarge = true;

		small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
		large_tree& ltree_tid = ltrees[tree_id];


		for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID v = que[que_i];
				descendants.push_back(v);
				vector<EdgeWeight>& adj_v = adj[v];
				NodeID adj_v_size = adj_v.size();
				for (size_t i = 0; i <adj_v_size; ++i) {
					NodeID w = adj_v[i];

					NodeID parentID = 0;
					if (isLarge)
						parentID = get_large_parent(w, ltree_tid);
					else
						parentID = get_small_parent(w, stree_tid);

					if (!vis[w] && parentID == v) {
						que[que_h++] = w;
						vis[w] = true;
					}
				}
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}

		for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
	}

	void get_descendants(int tree_id, NodeID picked_v, vector<NodeID>& descendants, vector<NodeID>& descendants_parents, vector<NodeID>& que, vector<bool>& vis, vector<vector<NodeID> >& adj, const bool isLarge, small_tree& stree_tid, large_tree& ltree_tid) {

		descendants.clear();
		//	descendants.reserve(V / 8);

		descendants.push_back(picked_v);
		descendants_parents[picked_v] = -1;

		vis[picked_v] = true;

		/*	bool isLarge = false;
		int stree_idx = trees_pointers[tree_id];
		if (stree_idx == -1)
		isLarge = true;

		small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
		large_tree& ltree_tid = ltrees[tree_id];*/


		for (int i = 0; i < descendants.size(); ++i) {
			NodeID v = descendants[i];
			vector<EdgeWeight>& adj_v = adj[v];
			NodeID adj_v_size = adj_v.size();
			for (int j = 0; j < adj_v_size; ++j) {
				NodeID w = adj_v[j];

				NodeID parentID = 0;


				if (isLarge)
					parentID = get_large_parent(w, ltree_tid);
				else {
					//	double get_small_parent_incre_time = GetCurrentTimeSec();
					parentID = get_small_parent(w, stree_tid);
					//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
				}


				if (!vis[w] && parentID == v) {
					descendants.push_back(w);
					descendants_parents[w] = v;
					vis[w] = true;
				}
			}
		}
		int dsize = descendants.size();
		for (int i = 0; i <dsize; ++i) vis[descendants[i]] = false;
	}

	void get_descendants(int tree_id, NodeID picked_v, vector<NodeID>& descendants, vector<NodeID>& descendants_parents, vector<NodeID>& que, vector<bool>& vis, Graph& graph, const bool isLarge, small_tree& stree_tid, large_tree& ltree_tid) {

		descendants.clear();
		//	descendants.reserve(V / 8);

		descendants.push_back(picked_v);
		descendants_parents[picked_v] = -1;

		vis[picked_v] = true;

		/*	bool isLarge = false;
		int stree_idx = trees_pointers[tree_id];
		if (stree_idx == -1)
		isLarge = true;

		small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
		large_tree& ltree_tid = ltrees[tree_id];*/


		for (int i = 0; i < descendants.size(); ++i) {
			NodeID v = descendants[i];

			//vector<EdgeWeight>& adj_v = adj[v];
			//NodeID adj_v_size = adj_v.size();

			//for (int j = 0; j < adj_v_size; ++j) {
			//NodeID w = adj_v[j];

			//Array Representation

			for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid) {

				NodeID w = graph.edges[eid];

				NodeID parentID = 0;


				if (isLarge)
					parentID = get_large_parent(w, ltree_tid);
				else {
					//	double get_small_parent_incre_time = GetCurrentTimeSec();
					parentID = get_small_parent(w, stree_tid);
					//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
				}


				if (!vis[w] && parentID == v) {
					descendants.push_back(w);
					descendants_parents[w] = v;
					vis[w] = true;
				}
			}
		}
		int dsize = descendants.size();
		for (int i = 0; i <dsize; ++i) vis[descendants[i]] = false;
	}


	void get_descendants(int tree_id, NodeID picked_v, vector<NodeID>& descendants, vector<NodeID>& descendants_parents, vector<NodeID>& que, vector<bool>& vis, WGraph& wgraph, const bool isLarge, small_tree& stree_tid, large_tree& ltree_tid) {

		descendants.clear();
		//	descendants.reserve(V / 8);

		descendants.push_back(picked_v);
		descendants_parents[picked_v] = -1;

		vis[picked_v] = true;

		/*	bool isLarge = false;
		int stree_idx = trees_pointers[tree_id];
		if (stree_idx == -1)
		isLarge = true;

		small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
		large_tree& ltree_tid = ltrees[tree_id];*/


		for (int i = 0; i < descendants.size(); ++i) {
			NodeID v = descendants[i];

			/*	vector<EdgeWeight>& adj_v = wgraph.adj[v];
			NodeID adj_v_size = adj_v.size();

			for (int j = 0; j < adj_v_size; ++j) {
			NodeID w = adj_v[j];*/

			//Array Representation

			for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid) {

				NodeID w = wgraph.edges[eid].first;

				//if(wgraph.adj[v][eid - wgraph.vertices[v]] != w) cout << "we got problems. get descent" << endl;

				NodeID parentID = 0;


				if (isLarge)
					parentID = get_large_parent(w, ltree_tid);
				else {
					//	double get_small_parent_incre_time = GetCurrentTimeSec();
					parentID = get_small_parent(w, stree_tid);
					//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
				}


				if (!vis[w] && parentID == v) {
					descendants.push_back(w);
					descendants_parents[w] = v;
					vis[w] = true;
				}
			}
		}
		int dsize = descendants.size();
		for (int i = 0; i <dsize; ++i) vis[descendants[i]] = false;
	}

	void get_ascendants(int tree_id, NodeID picked_v, vector<NodeID>& ascendants, vector<NodeID>& que, vector<bool>& vis, vector<vector<NodeID> >& adj, vector<vector<NodeID> >& r_adj) {
		ascendants.clear();

		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
		//	vector<vector<NodeID> >& adj = graph.adj;
		//vector<vector<NodeID> >& r_adj = graph.r_adj;

		que[que_h++] = picked_v;
		vis[picked_v] = true;
		que_t1 = que_h;


		bool isLarge = false;
		int stree_idx = trees_pointers[tree_id];
		if (stree_idx == -1)
			isLarge = true;

		small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
		large_tree& ltree_tid = ltrees[tree_id];

		for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID v = que[que_i];
				ascendants.push_back(v);

				NodeID parentID = 0;
				if (isLarge)
					parentID = get_large_parent(v, ltree_tid);
				else
					parentID = get_small_parent(v, stree_tid);

				if (DIRECTED_FLAG == false) {
					vector<EdgeWeight>& adj_v = adj[v];
					int adj_v_size = adj_v.size();
					for (int i = 0; i < adj_v_size; ++i) {
						NodeID w = adj_v[i];


						//if (!vis[w] && get_parent(tree_id, v) == w) {
						if (!vis[w] && parentID == w) {
							que[que_h++] = w;
							vis[w] = true;
						}
					}
				}
				else {
					vector<EdgeWeight>& r_adj_v = r_adj[v];
					int r_adj_v_size = r_adj_v.size();


					for (int i = 0; i < r_adj_v_size; ++i) {
						NodeID w = r_adj_v[i];

						//	if (!vis[w] && get_parent(tree_id, v) == w) {
						if (!vis[w] && parentID == w) {
							que[que_h++] = w;
							vis[w] = true;
						}
					}
				}
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}

		for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
	}

	void get_ascendants(int tree_id, NodeID picked_v, vector<NodeID>& ascendants, vector<NodeID>& que, vector<bool>& vis, vector<vector<NodeID> >& adj, vector<vector<NodeID> >& r_adj, const bool isLarge, small_tree& stree_tid, large_tree& ltree_tid) {
		ascendants.clear();

		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
		//	vector<vector<NodeID> >& adj = graph.adj;
		//vector<vector<NodeID> >& r_adj = graph.r_adj;

		que[que_h++] = picked_v;
		vis[picked_v] = true;
		que_t1 = que_h;


		//bool isLarge = false;
		//int stree_idx = trees_pointers[tree_id];
		//if (stree_idx == -1) 
		//	isLarge = true;
		//
		//small_tree& stree_tid = (stree_idx != -1 )? strees[stree_idx] : empty_small_tree;
		//large_tree& ltree_tid = ltrees[tree_id];

		for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID v = que[que_i];
				ascendants.push_back(v);

				NodeID parentID = 0;
				if (isLarge)
					parentID = get_large_parent(v, ltree_tid);
				else {
					//double get_small_parent_incre_time = GetCurrentTimeSec();
					parentID = get_small_parent(v, stree_tid);
					//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
				}
				if (DIRECTED_FLAG == false) {
					vector<EdgeWeight>& adj_v = adj[v];
					int adj_v_size = adj_v.size();
					for (int i = 0; i < adj_v_size; ++i) {
						NodeID w = adj_v[i];

						//if (!vis[w] && get_parent(tree_id, v) == w) {
						if (!vis[w] && parentID == w) {
							que[que_h++] = w;
							vis[w] = true;
						}
					}
				}
				else {
					vector<EdgeWeight>& r_adj_v = r_adj[v];
					int r_adj_v_size = r_adj_v.size();

					for (int i = 0; i < r_adj_v_size; ++i) {
						NodeID w = r_adj_v[i];

						if (!vis[w] && parentID == w) {
							que[que_h++] = w;
							vis[w] = true;
						}
					}
				}
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}

		for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
	}

	void get_ascendants(int tree_id, NodeID picked_v, vector<NodeID>& ascendants, vector<NodeID>& que, vector<bool>& vis, Graph& graph, const bool isLarge, small_tree& stree_tid, large_tree& ltree_tid) {
		ascendants.clear();

		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
		//	vector<vector<NodeID> >& adj = graph.adj;
		//vector<vector<NodeID> >& r_adj = graph.r_adj;

		que[que_h++] = picked_v;
		vis[picked_v] = true;
		que_t1 = que_h;


		//bool isLarge = false;
		//int stree_idx = trees_pointers[tree_id];
		//if (stree_idx == -1) 
		//	isLarge = true;
		//
		//small_tree& stree_tid = (stree_idx != -1 )? strees[stree_idx] : empty_small_tree;
		//large_tree& ltree_tid = ltrees[tree_id];

		for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID v = que[que_i];
				ascendants.push_back(v);

				NodeID parentID = 0;
				if (isLarge)
					parentID = get_large_parent(v, ltree_tid);
				else {
					//double get_small_parent_incre_time = GetCurrentTimeSec();
					parentID = get_small_parent(v, stree_tid);
					//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
				}
				if (DIRECTED_FLAG == false) {
					//vector<EdgeWeight>& adj_v = adj[v];
					//	int adj_v_size = adj_v.size();
					//	for (int i = 0; i < adj_v_size; ++i) {
					//		NodeID w = adj_v[i];

					//Array Representation
					for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid) {

						NodeID w = graph.edges[eid];

						//if (!vis[w] && get_parent(tree_id, v) == w) {
						if (!vis[w] && parentID == w) {
							que[que_h++] = w;
							vis[w] = true;
						}
					}
				}
				else {
					//vector<EdgeWeight>& r_adj_v = r_adj[v];
					//	int r_adj_v_size = r_adj_v.size();
					//for (int i = 0; i < r_adj_v_size; ++i) {
					//	NodeID w = r_adj_v[i];

					//Array Representation
					for (EdgeID eid = graph.r_vertices[v]; eid < graph.r_vertices[v + 1]; ++eid) {

						NodeID w = graph.r_edges[eid];

						if (!vis[w] && parentID == w) {
							que[que_h++] = w;
							vis[w] = true;
						}
					}
				}
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}

		for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
	}

	int labeling_source_bfs(NodeID source, vector<NodeID>& que, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<bool>& usd, vector<bool>& r_usd, NodeID ranking, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& r_tmp_idx, Graph& graph, vector<EdgeWeight>& tmp_d, vector<std::pair<uint64_t, uint64_t> >& tmp_s, vector<EdgeWeight>& r_tmp_d, vector<std::pair<uint64_t, uint64_t> >& r_tmp_s, vector<std::pair<NodeID, NodeID> >& sibling_es, vector<std::pair<NodeID, NodeID> >& child_es, vector<cover_value_type>& sum, unordered_set<NodeID>& picked_neighbors, const vector<vector<cover_value_type> >& cnt16) {
		 
		if (usd[source]) {
			if (ranking - skip_count < kNumBitParallelRoots) {
				skip_count++;
			}
			return 0;
		}

		if (ranking - skip_count < kNumBitParallelRoots) {

			index_t_bp<kNumBitParallelRoots>*& index_ = (DIRECTED_FLAG == true) ? dbplabels.index_bp : bplabels.index_bp;
			index_t_bp<kNumBitParallelRoots>*& bindex_ = dbplabels.bindex_bp;
			
			if (DIRECTED_FLAG == true) {
				fill(r_tmp_d.begin(), r_tmp_d.end(), INF_WEIGHT);
				fill(r_tmp_s.begin(), r_tmp_s.end(), std::make_pair(0, 0));
				fill(tmp_d.begin(), tmp_d.end(), INF_WEIGHT);
				fill(tmp_s.begin(), tmp_s.end(), std::make_pair(0, 0));
			}else {
				fill(tmp_d.begin(), tmp_d.end(), INF_WEIGHT);
				fill(tmp_s.begin(), tmp_s.end(), std::make_pair(0, 0));
			}

			if (DIRECTED_FLAG == true) {
				int que_t0 = 0, que_t1 = 0, que_h = 0;
				que[que_h++] = source;
				r_tmp_d[source] = 0;
				que_t1 = que_h;
				r_usd[source] = true;
				picked_neighbors.insert(source);
				int ns = 0;
				vector<int> vs;

				vector<NodeID> adj_r(graph.vertices[source + 1] - graph.vertices[source]);
				for (EdgeID eid = graph.vertices[source]; eid < graph.vertices[source + 1]; eid++) {
					adj_r[eid - graph.vertices[source]] = graph.edges[eid];
				}

				sort(adj_r.begin(), adj_r.end());

				vector<NodeID> r_adj_r(graph.r_vertices[source + 1] - graph.r_vertices[source]);
				for (EdgeID eid = graph.r_vertices[source]; eid < graph.r_vertices[source + 1]; eid++) {
					r_adj_r[eid - graph.r_vertices[source]] = graph.r_edges[eid];
				}

				sort(r_adj_r.begin(), r_adj_r.end());
				vector<NodeID> common_adj;
				set_intersection(adj_r.begin(), adj_r.end(), r_adj_r.begin(), r_adj_r.end(), back_inserter(common_adj));

				// Adj is sorted by descending sum[v] value.
				vector<pair<cover_value_type, NodeID> > sort_by_on_sum(common_adj.size());

			

				for (size_t i = 0; i < common_adj.size(); ++i) {
					NodeID v = common_adj[i];

					const vector<cover_value_type>& cnt16_v = cnt16[v];
					int sum_v = 0;
					int max1 = 0;
					int max2 = 0;
					for (int i = 0; i < CNT; ++i) {
						sum_v += cnt16_v[i];
						if (cnt16_v[i] > max2)
							max2 = cnt16_v[i];
						if (max2 > max1) {
							int tmp = max1;
							max1 = max2;
							max2 = tmp;
						}
					}

					sort_by_on_sum[i] = make_pair(sum_v - max1 - max2, v);
					//sort_by_on_sum[i] = make_pair((graph.vertices[v+1] - graph.vertices[v]) * (graph.r_vertices[v + 1] - graph.r_vertices[v]), v);
				}

				sort(sort_by_on_sum.rbegin(), sort_by_on_sum.rend());
				for (size_t i = 0; i < common_adj.size(); ++i) {
					common_adj[i] = sort_by_on_sum[i].second;
				}

				//sort(common_adj.begin(), common_adj.end());

				//forward search
				for (size_t i = 0; i < common_adj.size(); ++i) {
					NodeID v = common_adj[i];
					if (!r_usd[v]) {
						r_usd[v] = true;
						que[que_h++] = v;
						r_tmp_d[v] = 1;
						r_tmp_s[v].first = 1ULL << ns;
						vs.push_back(v);
						picked_neighbors.insert(v);
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
								sibling_es[num_sibling_es].first = v;
								sibling_es[num_sibling_es].second = tv;
								++num_sibling_es;
								//}
							}
							else {
								if (r_tmp_d[tv] == INF_WEIGHT) {
									que[que_h++] = tv;
									r_tmp_d[tv] = td;
								}
								child_es[num_child_es].first = v;
								child_es[num_child_es].second = tv;
								++num_child_es;
							}
						}
					}

					for (int i = 0; i < num_sibling_es; ++i) {
						int v = sibling_es[i].first, w = sibling_es[i].second;
						//r_tmp_s[v].second |= r_tmp_s[w].first;
						r_tmp_s[w].second |= r_tmp_s[v].first;
					}
					for (int i = 0; i < num_child_es; ++i) {
						int v = child_es[i].first, c = child_es[i].second;
						r_tmp_s[c].first |= r_tmp_s[v].first;
						r_tmp_s[c].second |= r_tmp_s[v].second;
					}

					que_t0 = que_t1;
					que_t1 = que_h;
				}

				for (NodeID v = 0; v < numOfVertices; ++v) {
					bindex_[v].bpspt_d[ranking - skip_count] = r_tmp_d[v];
					bindex_[v].bpspt_s[ranking - skip_count][0] = r_tmp_s[v].first;
					bindex_[v].bpspt_s[ranking - skip_count][1] = r_tmp_s[v].second & ~r_tmp_s[v].first;
				}

				int forward_search_space = que_h;


				//backward
				usd[source] = true;
				picked_neighbors.insert(source);
				fill(tmp_d.begin(), tmp_d.end(), INF_WEIGHT);
				fill(tmp_s.begin(), tmp_s.end(), std::make_pair(0, 0));

				que_t0 = 0, que_t1 = 0, que_h = 0;
				que[que_h++] = source;
				tmp_d[source] = 0;
				que_t1 = que_h;

				ns = 0;
				vs.clear();

				for (size_t i = 0; i < common_adj.size(); ++i) {
					NodeID v = common_adj[i];
					if (!usd[v]) {
						usd[v] = true;
						que[que_h++] = v;
						tmp_d[v] = 1;
						tmp_s[v].first = 1ULL << ns;
						vs.push_back(v);
						picked_neighbors.insert(v);
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
					index_[v].bpspt_d[ranking - skip_count] = tmp_d[v];
					index_[v].bpspt_s[ranking - skip_count][0] = tmp_s[v].first;
					index_[v].bpspt_s[ranking - skip_count][1] = tmp_s[v].second & ~tmp_s[v].first;
				}

				int backward_search_space = que_h;

				usd[source] = true;

				//cout << picked_neighbors.size() << endl;
				return (forward_search_space + backward_search_space); // * picked_neighbors.size();

			//	else
					//return forward_search_space;
				//	return forward_search_space * picked_neighbors.size();
			}
			else {

				usd[source] = true;

				fill(tmp_d.begin(), tmp_d.end(), INF_WEIGHT);
				fill(tmp_s.begin(), tmp_s.end(), std::make_pair(0, 0));

				int que_t0 = 0, que_t1 = 0, que_h = 0;
				que[que_h++] = source;
				tmp_d[source] = 0;
				que_t1 = que_h;
				picked_neighbors.insert(source);

				int ns = 0;
				vector<int> vs;

				vector<NodeID> adj_r(graph.vertices[source + 1] - graph.vertices[source]);
				for (EdgeID eid = graph.vertices[source]; eid < graph.vertices[source + 1]; eid++) {
					adj_r[eid - graph.vertices[source]] = graph.edges[eid];
				}

				vector<pair<cover_value_type, NodeID> > sort_by_on_sum(adj_r.size());
				for (size_t i = 0; i < adj_r.size(); ++i) {
					NodeID v = adj_r[i];
					
					const vector<cover_value_type>& cnt16_v = cnt16[v];
					int sum_v = 0;
					int max1 = 0;
					int max2 = 0;
					for (int i = 0; i < CNT; ++i) {
						sum_v += cnt16_v[i];
						if (cnt16_v[i] > max2)
							max2 = cnt16_v[i];
						if (max2 > max1) {
							int tmp = max1;
							max1 = max2;
							max2 = tmp;
						}
					}
					 
					sort_by_on_sum[i] = make_pair(sum_v - max1 - max2, v);

					//sort_by_on_sum[i] = make_pair(graph.vertices[v+1] - graph.vertices[v], v);
				}

				sort(sort_by_on_sum.rbegin(), sort_by_on_sum.rend());
				for (size_t i = 0; i < adj_r.size(); ++i) {
					adj_r[i] = sort_by_on_sum[i].second;
				}

				// Sort the adj_r based on current cover
				//sort(adj_r.begin(), adj_r.end());


					for (size_t i = 0; i < adj_r.size(); ++i) {
						NodeID v = adj_r[i];
						if (!usd[v]) {
							usd[v] = true;
							que[que_h++] = v;
							tmp_d[v] = 1;
							tmp_s[v].first = 1ULL << ns;
							vs.push_back(v);
							picked_neighbors.insert(v);
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
						index_[v].bpspt_d[ranking - skip_count] = tmp_d[v];
						index_[v].bpspt_s[ranking - skip_count][0] = tmp_s[v].first;
						index_[v].bpspt_s[ranking - skip_count][1] = tmp_s[v].second & ~tmp_s[v].first;
					}
					//return que_h;//as search space
					return que_h;// *picked_neighbors.size();
				}	
	}else {

			//Building normal labels.

			NodeID labeled_arcs = 0;

			//vector<vector<NodeID> >& adj = graph.adj;


			index_t_bp<kNumBitParallelRoots>*& index_ = (DIRECTED_FLAG == true) ? dbplabels.index_bp : bplabels.index_bp;
			index_t_bp<kNumBitParallelRoots>*& bindex_ = dbplabels.bindex_bp;

			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];
			const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[source];
			index_t_bp<kNumBitParallelRoots> &idx_r = index_[source];
			index_t_bp<kNumBitParallelRoots> &r_idx_r = bindex_[source];

			for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}

			NodeID que_t0 = 0, que_t1 = 0, que_h = 0;

			que[que_h++] = source;
			vis[source] = true;
			que_t1 = que_h;

			for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID v = que[que_i];
					// hubs touched so far when building this labels
					if (usd[v]) continue;


					if (DIRECTED_FLAG == false) {
						pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
						index_t_bp<kNumBitParallelRoots> &idx_v = index_[v];


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

						labeled_arcs++;
						tmp_idx_v.first.back() = ranking;
						tmp_idx_v.second.back() = d;
						tmp_idx_v.first.push_back(numOfVertices);
						tmp_idx_v.second.push_back(INF_WEIGHT);
					}
					else {
						pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];

						index_t_bp<kNumBitParallelRoots> &r_idx_v = bindex_[v];

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
								//hubs touched so far when building this labels
								goto pruned_forward;
							}
						}


						// arcs touched so far when building this labels
						labeled_arcs++;
						r_tmp_idx_v.first.back() = ranking;
						r_tmp_idx_v.second.back() = d;
						r_tmp_idx_v.first.push_back(numOfVertices);
						r_tmp_idx_v.second.push_back(INF_WEIGHT);
					}

					// Array Representation
					/*for (size_t i = 0; i < adj[v].size(); ++i) {
					NodeID w = adj[v][i];*/
					for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid) {

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
			for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
			}


			// Backward search.
			if (DIRECTED_FLAG == true) {

				for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i) {
					dst_r[r_tmp_idx_r.first[i]] = r_tmp_idx_r.second[i];
				}

				//vector<vector<NodeID> >& r_adj = graph.r_adj;

				que_t0 = 0, que_t1 = 0, que_h = 0;

				que[que_h++] = source;
				vis[source] = true;
				que_t1 = que_h;

				for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
					for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
						NodeID v = que[que_i];

						if (usd[v]) continue;

						pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
						index_t_bp<kNumBitParallelRoots> &idx_v = index_[v];

						for (int i = 0; i < kNumBitParallelRoots; ++i) {
							EdgeWeight td = r_idx_r.bpspt_d[i] + idx_v.bpspt_d[i];
							if (td - 2 <= d) {
								/*td +=
								(r_idx_r.bpspt_s[i][0] & idx_v.bpspt_s[i][0]) ? -2 :
								((r_idx_r.bpspt_s[i][0] & idx_v.bpspt_s[i][1]) |
								(r_idx_r.bpspt_s[i][1] & idx_v.bpspt_s[i][0]))
								? -1 : 0;*/
								td +=
									(idx_v.bpspt_s[i][0]) & r_idx_r.bpspt_s[i][0] ? -2 :
									((idx_v.bpspt_s[i][1] & r_idx_r.bpspt_s[i][0]) |
										(idx_v.bpspt_s[i][0] & r_idx_r.bpspt_s[i][1]))
									? -1 : 0;
								if (td <= d) goto pruned_backward;
							}
						}

						for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
							NodeID w = tmp_idx_v.first[i];
							EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
							if (td <= d) {
								goto pruned_backward;
							}
						}

						labeled_arcs++;
						tmp_idx_v.first.back() = ranking;
						tmp_idx_v.second.back() = d;
						tmp_idx_v.first.push_back(numOfVertices);
						tmp_idx_v.second.push_back(INF_WEIGHT);

						//	for (size_t i = 0; i < r_adj[v].size(); ++i) {
						//	NodeID w = r_adj[v][i];

						for (EdgeID eid = graph.r_vertices[v]; eid < graph.r_vertices[v + 1]; ++eid) {
							NodeID w = graph.r_edges[eid];
							if (!vis[w]) {
								// hubs touched so far when building this labels							
								que[que_h++] = w;
								vis[w] = true;
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
					dst_r[r_tmp_idx_r.first[i]] = INF_WEIGHT;
			}

			usd[source] = true;
			return labeled_arcs;
		}
	}

	void shrink_tree(int tid, NodeID root, large_tree& ltree_tid) {
		small_tree new_stree;
		//GOOGLE DENSE MAP
		new_stree.set_empty_key(numOfVertices + 1);
		new_stree.set_deleted_key(numOfVertices + 2);
		for (NodeID v = 0; v < ltree_tid.size(); ++v) {
			int pt = ltree_tid[v];
			if (pt != -1)
				new_stree.insert(make_pair(v, pt));
		}
		new_stree.insert(make_pair(root, -1));
		strees.push_back(new_stree);
		trees_pointers[tid] = strees.size() - 1;
		alive_resources += (long long)new_stree.size();

		//ltrees[tid].clear();
		// FIX
		ltree_tid.clear();
		ltree_tid.shrink_to_fit();
		alive_resources -= (long long)numOfVertices;
	}

	void enlarge_tree(int tid, NodeID root, small_tree& stree_tid, large_tree& ltree_tid) {
		ltree_tid.resize(numOfVertices, -1);

		for (small_tree::iterator it = stree_tid.begin(); it != stree_tid.end(); ++it) {
			NodeID v = (*it).first;
			NodeID pt = (*it).second;
			ltree_tid[v] = pt;
		}

		ltree_tid[root] = -1;

		ltrees.push_back(ltree_tid);

		trees_pointers[tid] = -1;
		alive_resources += (long long)numOfVertices;

		alive_resources -= (long long)stree_tid.size();
		stree_tid.clear();
		stree_tid.resize(0);
	}

	int get_large_parent(NodeID v, large_tree& ltree) {
		return ltree[v];
	}

	int get_small_parent_withouttest(NodeID v, small_tree& stree) {

		return stree[v];

	}

	int get_small_parent(NodeID v, small_tree& stree) {
		small_tree::iterator it = stree.find(v);

		if (it != stree.end())
			return stree[v];
		else
			return -1;
	}


	void remove_large_node(NodeID v, large_tree& ltree) {
		ltree[v] = -1;
	}
	void remove_small_node(NodeID v, small_tree& stree) {
		stree.erase(v);
		//	stree.resize(0);
		//stree[v] = -1;
	}

	//GOOGLE
	void cnt16_2max_queue(vector<cover_value_type>& picked_value, vector<cover_value_type>& sum, vector<bool>& has_picked, /*unordered_set<NodeID>& affected_v*/ google::dense_hash_set<NodeID>& affected_v, vector<bool>& affected_list, const vector<vector<cover_value_type> >& cnt16, bool& switch_small, bool& switch_set, benchmark::heap<2, cover_value_type, NodeID>& value_heap) {

		// First Round, build the entire cnt16 one time.
		if (picked_value.size() == 0) {
			picked_value.resize(numOfVertices, 0);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				NodeID max1 = 0;
				NodeID max2 = 0;
				const vector<cover_value_type>& cnt16_v = cnt16[v];
				cover_value_type &picked_value_v = picked_value[v];
				cover_value_type& sum_v = sum[v];
				sum_v = 0;
				for (int i = 0; i < CNT; ++i) {
					sum_v += cnt16_v[i];
					if (cnt16_v[i] > max2)
						max2 = cnt16_v[i];
					if (max2 > max1) {
						NodeID tmp = max1;
						max1 = max2;
						max2 = tmp;
					}
				}
				picked_value_v = sum_v - max1 - max2;
				if (switch_small == true) {
					value_heap.update(v, -picked_value_v);
				}
			}
		}

		//for (unordered_set<NodeID>::iterator it = affected_v.begin(); it != affected_v.end(); ++it) {
		if (!switch_set) {
			for (int v = 0; v < numOfVertices; ++v) {
				//NodeID v = *it;
				if (affected_list[v] == false) continue;
				affected_list[v] = false;
				if (has_picked[v] == true) continue;
				//if (v == best_vid) continue;
				cover_value_type& sum_v = sum[v];

				/*	if (sum_v < best_so_far) {
				picked_value[v] = sum_v;
				if (switch_small == true)
				value_heap.update(v, sum_v);
				lazy_updates[v] = true;
				continue;
				}*/
				NodeID max1 = 0;
				NodeID max2 = 0;
				const vector<cover_value_type>& cnt16_v = cnt16[v];
				cover_value_type& picked_value_v = picked_value[v];


				sum_v = 0;
				for (int i = 0; i < CNT; ++i) {
					sum_v += cnt16_v[i];
					if (cnt16_v[i] > max2)
						max2 = cnt16_v[i];
					if (max2 > max1) {
						int tmp = max1;
						max1 = max2;
						max2 = tmp;
					}
				}
				picked_value_v = sum_v - max1 - max2;
				//if (picked_value_v > best_so_far) best_so_far = picked_value_v;
				if (switch_small == true) {
					value_heap.update(v, -picked_value_v);
				}
			}
		}
		else {
			//GOOGLE

			//for (unordered_set<NodeID>::iterator it = affected_v.begin(); it != affected_v.end(); ++it) {
			for (google::dense_hash_set<NodeID>::iterator it = affected_v.begin(); it != affected_v.end(); ++it) {

				NodeID v = *it;
				affected_list[v] = false;
				if (has_picked[v] == true) continue;
				//if (v == best_vid) continue;
				cover_value_type& sum_v = sum[v];

				//if (sum_v < best_so_far) {
				//	picked_value[v] = sum_v;
				//	if (switch_small == true)
				//		value_heap.update(v, sum_v);
				//	lazy_updates[v] = true;
				//	continue;
				//}
				NodeID max1 = 0;
				NodeID max2 = 0;
				const vector<cover_value_type>& cnt16_v = cnt16[v];
				cover_value_type& picked_value_v = picked_value[v];


				sum_v = 0;
				for (int i = 0; i < CNT; ++i) {
					sum_v += cnt16_v[i];
					if (cnt16_v[i] > max2)
						max2 = cnt16_v[i];
					if (max2 > max1) {
						int tmp = max1;
						max1 = max2;
						max2 = tmp;
					}
				}
				picked_value_v = sum_v - max1 - max2;
				//if (picked_value_v > best_so_far) best_so_far = picked_value_v;
				if (switch_small == true) {
					value_heap.update(v, -picked_value_v);
				}
			}
			affected_v.clear();
		}
	}

	BP_Betweenness_Ordering(NodeID k, double beta, Graph& graph, NodeID stopv, long long bound_times = 10, long long trans_times = 8) {

		double start_time = GetCurrentTimeSec();



		if (k < CNT) {
			cout << " Should pick at least 16 seed!" << endl;
			return;
		}
		iteration_time.resize(numOfVertices);

		updating_time = 0;
		labeling_time = 0;
		adding_time = 0;
		selecting_time = 0;
		init_time = 0;
		u1_time = 0;
		u2_time = 0;
		u3_time = 0;
		num_of_trees = 0;
		total_resources = 0;
		get_small_parent_time = 0;

		inv.resize(numOfVertices, -1);
		rank.resize(numOfVertices, -1);


		//if (DIRECTED_FLAG == false) {
		//	labels.index_.resize(numOfVertices);
		//}
		//else { //if (DIRECTED_FLAG == true) {
		//	dlabels.index_.resize(numOfVertices);
		//	dlabels.bindex_.resize(numOfVertices);
		//}

		// Resource Bounding.
		NodeID new_tree_arcs = 0;
		NodeID labeled_arcs = 0;
		alive_resources = 0;
		bounded_resources = (long long)10 * (long long)k * (long long)numOfVertices;

		bool stop_flag = false;

		vector<NodeID> st_roots; // The root vertices of the sampled trees.
		vector<NodeID> st_sizes; // The size of the sampled trees.

		NodeID size_thr = numOfVertices < 8 ? numOfVertices : numOfVertices / 8; // When a tree shrink to a size smaller than the size_thr, it should be converted into a small_tree representation.
																				 //	vector<vector<int> > cover;


		init_time = GetCurrentTimeSec();
		// Building the initial k trees.
		//	srand((unsigned)time(NULL));
		srand(0x728f9d0);


		vector<NodeID> descendants_parents(numOfVertices, -1);

		ltrees.reserve(numOfVertices);
		strees.reserve(numOfVertices);
		trees_pointers.reserve(numOfVertices);
		st_roots.reserve(numOfVertices);
		st_sizes.reserve(numOfVertices);
		st_roots.resize(k);
		st_sizes.resize(k);

		ltrees.resize(k);
		strees.clear();
		trees_pointers.resize(k, -1); // Point the sampling trees to its small tree representations.
									  //	cover.resize(k, vector<int>(V, 0));

		vector<NodeID> ascendants;
		ascendants.reserve(numOfVertices);
		vector<NodeID> descendants;
		descendants.reserve(numOfVertices);
		vector<NodeID> que(numOfVertices);
		vector<bool> vis(numOfVertices);
		vector<bool> has_picked(numOfVertices);
		vector<vector<cover_value_type> > cnt16(numOfVertices, vector<cover_value_type>(CNT));
		//		vector<bool> affected(V, true);
		vector<bool> picked_as_root(numOfVertices, false);
		vector<bool> trivial_root(numOfVertices, false);
		int picked_count = 0;
		int used_as_cover = 0;
		vector<cover_value_type> sum(numOfVertices);
		vector<bool> usd(numOfVertices, false);
		vector<bool> r_usd(numOfVertices, false);
		vector<cover_value_type> cover_tid(numOfVertices, 0);

		//unordered_set<NodeID> affected_v; //GOOGLE
		google::dense_hash_set<NodeID> affected_v;
		affected_v.set_empty_key(numOfVertices + 1);
		affected_v.set_deleted_key(numOfVertices + 2);

		vector<bool> affected_list(numOfVertices, false);
		int affected_count = 0;
		bool switch_set = false;

		//vector<vector<int> > cover_record(V, vector<int>(V, 0));

		//vector<unordered_set<int> > reverse_index(0); //GOOGLE
		//		vector<unordered_set<int> > reverse_index(0);
		vector<google::dense_hash_set<NodeID> > reverse_index;
		reverse_index.reserve(numOfVertices);

		// Preparing basic structure for pl algorithm.
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT))); // Backward labels.


		vector<EdgeWeight> tmp_d(numOfVertices);
		vector<std::pair<uint64_t, uint64_t> > tmp_s(numOfVertices);
		vector<EdgeWeight> r_tmp_d(numOfVertices);
		vector<std::pair<uint64_t, uint64_t> > r_tmp_s(numOfVertices);

		//vector<NodeID> que(numOfVertices);
		vector<std::pair<NodeID, NodeID> > sibling_es(numOfEdges);
		vector<std::pair<NodeID, NodeID> > child_es(numOfEdges);



		index_t_bp<kNumBitParallelRoots>*& index_ = (DIRECTED_FLAG == true )?  dbplabels.index_bp : bplabels.index_bp;
		index_t_bp<kNumBitParallelRoots>*& bindex_ = dbplabels.bindex_bp;

		index_ = (index_t_bp<kNumBitParallelRoots>*)memalign(64, numOfVertices * sizeof(index_t_bp<kNumBitParallelRoots>));
		bindex_ = (index_t_bp<kNumBitParallelRoots>*)memalign(64, numOfVertices * sizeof(index_t_bp<kNumBitParallelRoots>));
													 //	int best_sum = 0;
													 //	NodeID best_sum_vid = 0;
													 //vector<bool> lazy_updates(V, false);

													 // Building the initial k sampled trees. 
		for (int i = 0; i < k; ++i) {
			NodeID r = rand() % numOfVertices;
			if (picked_as_root[r] == true) { --i; continue; }

			st_roots[i] = r;
			large_tree &ltree_i = ltrees[i];
			picked_as_root[r] = true;
			picked_count++;
			used_as_cover++;

			// Array Representation
			//build_large_tree_bfs(r, ltree_i, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph.adj, 0, numOfVertices);
			build_large_tree_bfs(r, ltree_i, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph, 0, numOfVertices, 0);

			st_sizes[i] = descendants.size();
			alive_resources += (long long)numOfVertices;

			// Turn large tree to small tree.
			// 1. Build new small tree with ltree_i and add it to strees;
			// 2. Set tree_pointers[i] to the position of the strees;
			// 3. Clear the large tree ltree_i by setting its size to 0.
			if (st_sizes[i] < size_thr)
				shrink_tree(i, st_roots[i], ltrees[i]);


			// Calculate the subtree size of each nodes			.
			//	vector<int>& cover_i = cover[i];

			bool isLarge = false;
			int stree_idx = trees_pointers[i];
			if (stree_idx == -1)
				isLarge = true;

			small_tree& stree_i = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;

			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];

				++cover_tid[q]; // Increment because it hitsthe r-q path.
				++cnt16[q][i%CNT];
				//++sum[q];
				//	if (sum[q] > best_sum) {
				//		best_sum = sum[q];
				//		best_sum_vid = q;
				//	}
				////	cover_record[i][q]++;

				int parent_node;
				if (isLarge)
					parent_node = get_large_parent(q, ltree_i);
				else
					parent_node = get_small_parent(q, stree_i);


				if (parent_node > -1) {
					cover_tid[parent_node] += cover_tid[q];
					cnt16[parent_node][i%CNT] += cover_tid[q];
					//sum[parent_node] += cover_i[q];
					//	cover_record[i][parent_node] += cover_record[i][q];
				}
			}
			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];
				cover_tid[q] = 0;
			}
		}


		//cnt16_aggreate(cnt16, has_picked, cover);

		vector<cover_value_type> picked_value;
		benchmark::heap<2, cover_value_type, NodeID> value_heap(numOfVertices);
		bool switch_small = false;

		//cnt16_2max(picked_value, sum, has_picked,affected, cnt16, switch_small, value_heap);
		//		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap, best_sum_vid, lazy_updates);
		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);

		init_time = GetCurrentTimeSec() - init_time;


		cout << "Finish initial sampling " << init_time << "secs." << endl;

		//for (NodeID v = 0; v < V; ++v)	value_heap.update(v, -picked_value[v]); // Min-heap.

		double max_value = -1;


		// Perform V iterations to rank all vertices.
		NodeID ranking = 0;
		for (; ranking < numOfVertices; ++ranking) {

			//cout << ranking << endl;
			/*	if (ranking % (numOfVertices / 10) == 0)
			cout << "picking the " << ranking << "-th orders." << endl;*/

			double time_this_iteration = GetCurrentTimeSec();

			//best_sum = 0;

			NodeID pv = 0;
			cover_value_type pv_cover = -1;

			if (ranking == stopv) {
				stop_flag = true;
				cout << "topv labeling time:" << labeling_time << endl;
				break;
			}
			double selecting_acc_time = GetCurrentTimeSec();
			// In the beginning, pick the best cover vertex by a linear scan.
			if (switch_small == false) {
				double u1_acc_time = GetCurrentTimeSec();
				//while (true) {
				for (NodeID v = 0; v < numOfVertices; ++v) {
					if (has_picked[v] == true) continue;
					if (pv_cover == picked_value[v]) {
						if (sum[pv] < sum[v]) {
							pv_cover = picked_value[v];
							pv = v;
						}
					}
					else if (pv_cover < picked_value[v]) {
						pv_cover = picked_value[v];
						pv = v;
					}
				}
				/*	if (lazy_updates[pv] == true) {
				NodeID max1 = 0;
				NodeID max2 = 0;
				const vector<int>& cnt16_pv = cnt16[pv];
				int& picked_value_v = picked_value[pv];

				int sum_v = sum[pv];
				for (int i = 0; i < CNT; ++i) {
				if (cnt16_pv[i] > max2)
				max2 = cnt16_pv[i];
				if (max2 > max1) {
				int tmp = max1;
				max1 = max2;
				max2 = tmp;
				}
				}
				picked_value_v = sum_v - max1 - max2;

				lazy_updates[pv] = false;
				}
				else
				break;*/
				//}

				if (sum[pv] < numOfVertices / 8 && switch_small == false) {
					cout << ranking << endl;
					cout << pv_cover << endl;
					cout << sum[pv] << endl;
					cout << picked_count << endl;
					cout << "Switch to max-heap representations." << endl;
					switch_small = true;
				}
				u1_time += (GetCurrentTimeSec() - u1_acc_time);
			}
			else {// Later, utilizing the max-heap directly
				double u2_acc_time = GetCurrentTimeSec();
				//while (true) {
				value_heap.extract_min(pv, pv_cover);
				pv_cover = -pv_cover;
				/*	if (lazy_updates[pv] == true) {
				NodeID max1 = 0;
				NodeID max2 = 0;
				const vector<int>& cnt16_pv = cnt16[pv];
				int& picked_value_v = picked_value[pv];

				int sum_v = sum[pv];
				for (int i = 0; i < CNT; ++i) {
				if (cnt16_pv[i] > max2)
				max2 = cnt16_pv[i];
				if (max2 > max1) {
				int tmp = max1;
				max1 = max2;
				max2 = tmp;
				}
				}
				picked_value_v = sum_v - max1 - max2;
				value_heap.update(pv, -picked_value_v);
				lazy_updates[pv] = false;
				}
				else
				break;
				}*/
				u2_time += (GetCurrentTimeSec() - u2_acc_time);

			}

			selecting_time += (GetCurrentTimeSec() - selecting_acc_time);

			//// FIX
			//if (switch_small == true && reverse_index.size() > 0 && sum[pv] <= reverse_index[pv].size()) {
			//	cout << "I think we can stop here " << GetCurrentTimeSec() - start_time << " ranking:" << ranking << endl;
			//	break;
			//}

			has_picked[pv] = true;
			unordered_set<NodeID> picked_neighbors;
			double labeling_acc_time = GetCurrentTimeSec();
		//	labeled_arcs = labeling_source_bfs(pv, que, vis, dst_r, usd, r_usd, ranking, tmp_idx, r_tmp_idx, graph);
			labeled_arcs = labeling_source_bfs(pv, que, vis, dst_r, usd, r_usd, ranking, tmp_idx, r_tmp_idx, graph, tmp_d, tmp_s, r_tmp_d, r_tmp_s, sibling_es, child_es, sum, picked_neighbors, cnt16);
			labeling_time += (GetCurrentTimeSec() - labeling_acc_time);
			inv[ranking] = pv;
			rank[pv] = ranking;

			//cout << ranking << "\t" << labeled_arcs << endl;


			int current_trees_num = st_roots.size();

			// 4 steps in updating:
			double updating_acc_time = GetCurrentTimeSec();
			if (switch_small == false) { // We still iterate all alive trees to calculate the cover.	
				for (int tid = 0; tid < current_trees_num; ++tid) {

					// This tree has already be taken because its root has been selected.
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					// pv has already been covered in tree tid.
					int pvpid = 0;
					if (isLarge)
						pvpid = get_large_parent(pv, ltree_tid);
					else
						pvpid = get_small_parent(pv, stree_tid);
					if (pvpid == -1 && st_roots[tid] != pv) // pv has been covered in tree tid 
						continue;


					// Get the descendants and ascendants of pv in tree tid.

					// Array Representation
					//// 1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;				 
					//get_descendants(tid, pv, descendants, descendants_parents, que, vis, graph.adj, isLarge, stree_tid, ltree_tid);

					//get_ascendants(tid, pv, ascendants, que, vis, graph.adj, graph.r_adj, isLarge, stree_tid, ltree_tid);

					get_descendants(tid, pv, descendants, descendants_parents, que, vis, graph, isLarge, stree_tid, ltree_tid);

					get_ascendants(tid, pv, ascendants, que, vis, graph, isLarge, stree_tid, ltree_tid);

					NodeID subtract_size = descendants.size();
					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
												//cover[tid][av] -= subtract_size;
						cnt16[av][tid%CNT] -= subtract_size;
						//sum[av] -= subtract_size;
						//	cover_record[tid][av] -= subtract_size;
						//	if (cover_record[tid][av] < 0)
						//	cout << "1: " << cover_record[tid][av] << endl;
						if (affected_list[av] == false) {
							affected_list[av] = true;
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
						}
					}


					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					//	vector<int> cover_tid(numOfVertices, 0);
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];

						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//sum[dv] -= cover_tid[dv];
						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];

						/*	if (cnt16[dv][tid%CNT] < 0)
						cout << "a6 " << cnt16[dv][tid%CNT] << endl;*/

						//	cover_record[tid][dv] -= cover_tid[dv];
						//	if (cover_record[tid][dv] < 0) {
						//	cout << "2: " << cover_record[tid][dv]  << " " << cover_tid[dv]  << " " << tid << " " << dv <<  endl;
						//	}

						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else
							remove_small_node(dv, stree_tid);


						//	remove_node(tid, dv);
						//affected[dv] = true;
						if (affected_list[dv] == false) {
							affected_count++;
							affected_list[dv] = true;
							if (switch_set)
								affected_v.insert(dv);
						}
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;
					}


					// Drop the size of subtree.
					st_sizes[tid] -= subtract_size;
					if (trees_pointers[tid] != -1) alive_resources -= (cover_value_type)subtract_size;
					//	if (st_sizes[tid] < 0) {
					//		cout << "a5" << st_sizes[tid] << endl;
					//	}
					if (st_sizes[tid] > 0) {
						// 3. Convert shrink trees to small trees representataions.				
						if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
							shrink_tree(tid, st_roots[tid], ltrees[tid]);
						}
					}

					if (st_roots[tid] == pv)
						if (isLarge) {
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
							//ltrees.shrink_to_fit();
						}
						else {
							stree_tid.clear();
							stree_tid.resize(0);
							//strees.shrink_to_fit();
						}

				}
			}
			else { // We use the reverse index to update trees.

				   // Build the reverse_index first (one time).
				if (reverse_index.size() == 0) {
					reverse_index.resize(numOfVertices);
					for (NodeID v = 0; v < numOfVertices; ++v) {
						reverse_index[v].set_empty_key(numOfVertices + 1);
						reverse_index[v].set_deleted_key(numOfVertices + 2);
					}
					for (int tid = 0; tid < current_trees_num; ++tid) {

						if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue; //This whole tree has already been covered. 

						bool isLarge = false;
						int stree_idx = trees_pointers[tid];
						if (stree_idx == -1)
							isLarge = true;

						small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
						large_tree& ltree_tid = ltrees[tid];

						for (int v = 0; v < numOfVertices; ++v) {

							// Get the descendants and ascendants of pv in tree tid.

							//	if (cover[tid].find(pv) == cover[tid].end()) continue; // This pv has been covered in this tree before.
							//	if (cover[tid][pv] <= 0) continue; // modified-5-14
							int vpid = 0;
							if (isLarge)
								vpid = get_large_parent(v, ltree_tid);
							else
								vpid = get_small_parent(v, stree_tid);

							// To testify whether v is still uncovered in tree tid
							// v has parent vpid != -1, but if v is the root of tree tid, it still can be uncovered.
							if ((vpid != -1 || st_roots[tid] == v))
								reverse_index[v].insert(tid);
							//if (cover[tid].find(v) != cover[tid].end())
							//	reverse_index[v].insert(tid);
						}
					}
				}

				double u4_acc_time = GetCurrentTimeSec();

				//	double u4_acc_time = GetCurrentTimeSec();
				// Access the alive trees directly by the reverse index.
				//				cout << reverse_index[pv].size() << " " << sum[pv] << endl;


				//				for (unordered_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) { // GOOGLE
				for (google::dense_hash_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) {


					int tid = *it;
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					// Get the descendants and ascendants of pv in tree tid.

					//	if (cover[tid].find(pv) == cover[tid].end()) continue; // This pv has been covered in this tree before.
					//	if (cover[tid][pv] <= 0) continue; // modified-5-14

					int pvpid = 0;
					if (isLarge)
						pvpid = get_large_parent(pv, ltree_tid);
					else {
						//	double get_small_parent_incre_time = GetCurrentTimeSec();
						pvpid = get_small_parent(pv, stree_tid);
						//	get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}

					if (pvpid == -1 && st_roots[tid] != pv)
						continue;

					// Array Representation
					//// 1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;
					//get_descendants(tid, pv, descendants, descendants_parents, que, vis, graph.adj, isLarge, stree_tid, ltree_tid);

					//get_ascendants(tid, pv, ascendants, que, vis, graph.adj, graph.r_adj, isLarge, stree_tid, ltree_tid);

					get_descendants(tid, pv, descendants, descendants_parents, que, vis, graph, isLarge, stree_tid, ltree_tid);

					get_ascendants(tid, pv, ascendants, que, vis, graph, isLarge, stree_tid, ltree_tid);

					NodeID subtract_size = descendants.size();

					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
												//						cover[tid][av] -= subtract_size;

						cnt16[av][tid%CNT] -= subtract_size;
						//sum[av] -= subtract_size;

						//	cover_record[tid][av] -= subtract_size;
						//	if (cover_record[tid][av] < 0) {
						//	cout << "4: " << cover_record[tid][av] << " " << subtract_size << " " << tid << " " << endl;
						//	}

						//affected[av] = true;
						/*if (cnt16[av][tid%CNT] < 0)
						cout << " a2 " << cnt16[av][tid%CNT]  << " " << av << " " << st_roots[tid] << " " << tid << endl;*/
						if (affected_list[av] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
							affected_list[av] = true;
						}
					}

					//double u2_acc_time = GetCurrentTimeSec();
					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					//vector<int> cover_tid(numOfVertices, 0);
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//	sum[dv] -= cover_tid[dv];
						//cover[tid].erase(dv);
						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];

						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else
							remove_small_node(dv, stree_tid);

						//	remove_node(tid, dv);
						//affected[dv] = true;

						//	cover_record[tid][dv] -= cover_tid[dv];
						//if (cover_record[tid][dv] < 0)
						//	cout << "3: " << cover_record[tid][dv] << " " << cover_tid[dv] << " " << di << " " << descendants.size() << " " << picked_as_root[st_roots[tid]] << endl;

						/*if (cnt16[dv][tid%CNT]  < 0)
						cout << " a3 " << cnt16[dv][tid%CNT] << " " << dv << " " << st_roots[tid] << " " << tid << endl;
						*/
						if (dv != pv)
							reverse_index[dv].erase(tid);
						//in_tree_count[dv]--;


						if (affected_list[dv] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(dv);
							affected_list[dv] = true;
						}
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;
					}

					if (st_roots[tid] == pv)
						if (isLarge) {
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
						}
						else {
							stree_tid.clear();
							stree_tid.resize(0);
						}

						st_sizes[tid] -= subtract_size;
						//	if (st_sizes[tid]  < 0)
						//		cout << " a4 " << st_sizes[tid] << endl;

						if (trees_pointers[tid] != -1) alive_resources -= subtract_size;


						// 3. Convert shrink trees to small trees representataions.	
						if (st_sizes[tid] > 0) {
							if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
								shrink_tree(tid, st_roots[tid], ltrees[tid]);
							}
						}

				}
				reverse_index[pv].clear();
				reverse_index[pv].resize(0);
				u4_time += (GetCurrentTimeSec() - u4_acc_time);
			}

			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			// 4. Sample new trees based on.
			// If we sample a new tree, we need to dynamically increment the following data structures:
			// 1) ltrees
			// 2) trees_pointers (-1 or pointer to strees)
			// 3) cover
			// 4) st_roots
			// 5) st_sizes 
			// 6) strees (maybe)
			// two constraints: cl > ct and 10kn total vertices
			double adding_acc_time = GetCurrentTimeSec();
			new_tree_arcs = 0;

			//cout << "ranking:" << ranking << endl;
			//cout << "pick cover:" << pv_cover << endl;
			//cout << "pick count:" << picked_count << endl;
			//cout << "labeled_arcs:" << labeled_arcs << endl;
			//cout << "new_tree_arcs:" << new_tree_arcs << endl;
			//cout << "alive_resources:" << alive_resources << endl;
			//cout << "bounded_resources:" << bounded_resources << endl;
			//cout << "#tree:" << cover.size() << endl;
			int how_many = 0;
			//cout << "beta*labeld_arcs: " << beta * labeled_arcs << " new tree _arcs: " << new_tree_arcs << " alive " << alive_resources << " bounded " << bounded_resources << endl;
			NodeID sam_this_round = 0;
			while (beta * labeled_arcs > new_tree_arcs && alive_resources < bounded_resources) {
			//	cout << beta * labeled_arcs << "," << new_tree_arcs << endl;
				if (picked_count == numOfVertices) break;



				NodeID r = rand() % numOfVertices;
				if (picked_as_root[r] == true) continue;
				if (trivial_root[r] == true) continue;
				int tmp_tree_arcs = 0;

				large_tree new_ltree;
				small_tree new_stree;
				//GOOGLE DENSE MAP
				new_stree.set_empty_key(numOfVertices + 1);
				new_stree.set_deleted_key(numOfVertices + 2);

				//	BUILD_SMALL_TREE_FLAG = false;
				bool isLarge = false;
				if (!BUILD_SMALL_TREE_FLAG) {

					double u5_acc_time = GetCurrentTimeSec();
					
					tmp_tree_arcs = build_large_tree_bfs(r, new_ltree, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph, new_tree_arcs, beta * labeled_arcs, ranking);

					u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;

					if (descendants.size() <size_thr)
						BUILD_SMALL_TREE_FLAG = true;

					if (descendants.size() < 2) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					isLarge = true;
					how_many++;
					used_as_cover++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					st_roots.push_back(r);
					trees_pointers.push_back(-1);
					st_sizes.push_back(descendants.size());

					//new_tree_arcs += descendants.size();
					alive_resources += numOfVertices;

					// Turn large tree to small tree.
					// 1. Build new small tree with ltree_i and add it to strees;
					// 2. Set tree_pointers[i] to the position of the strees;
					// 3. Clear the large tree ltree_i by setting its size to 0.
					int i = st_roots.size() - 1;
					if (st_sizes[i] < size_thr) {
						shrink_tree(i, st_roots[i], new_ltree);
						isLarge = false;
					}
					ltrees.push_back(new_ltree);
				}
				else {//Build Small Tree Directly
					double u5_acc_time = GetCurrentTimeSec();
					// Array Representation
					//	tmp_tree_arcs = build_small_tree_bfs(r, new_stree, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph.adj, new_tree_arcs, beta * labeled_arcs);
					tmp_tree_arcs = build_small_tree_bfs(r, new_stree, descendants, descendants_parents, que, vis, dst_r, usd, tmp_idx, r_tmp_idx, graph, new_tree_arcs, beta * labeled_arcs, ranking);
					u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;

					if (descendants.size() >= size_thr)
						BUILD_SMALL_TREE_FLAG = false;

					if (descendants.size() < 2) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					isLarge = false;

					how_many++;
					used_as_cover++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					st_roots.push_back(r);
					trees_pointers.push_back(strees.size()); //we assume it is inserted.
					st_sizes.push_back(descendants.size());

					alive_resources += new_stree.size();
					//strees.push_back(new_stree);
					//// new_ltree is empty.
					//ltrees.push_back(new_ltree);

					int i = st_roots.size() - 1;
					if (st_sizes[i] >= size_thr) {
						enlarge_tree(i, st_roots[i], new_stree, new_ltree);
						new_stree.clear();
						new_stree.resize(0);
						isLarge = true;
					}
					else {
						strees.push_back(new_stree);
						ltrees.push_back(new_ltree);
					}
				}

				int i = st_roots.size() - 1;

				int numOfStrees = strees.size();
				small_tree &new_stree_ref = numOfStrees != 0 ? strees[strees.size() - 1] : empty_small_tree;

				// Calculate the subtree size of each nodes.
				//				unordered_map<int, int> cover_new;
				//	vector<cover_value_type> cover_new(numOfVertices, 0); // modified-5-14

				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];

					++cover_tid[q]; // Increment because it hitsthe r-q path. 
					++cnt16[q][i%CNT];
					//++sum[q];
					/*	if (sum[q] > best_sum) {
					best_sum = sum[q];
					best_sum_vid = q;
					}*/
					//affected[q] = true;
					//	cover_record[i][q]++;
					if (affected_list[q] == false) {
						affected_count++;
						if (switch_set)
							affected_v.insert(q);
						affected_list[q] = true;
					}

					if (switch_small == true) {
						reverse_index[q].insert(i);
					}

					int parent_node = -1;

					if (isLarge)
						parent_node = get_large_parent(q, new_ltree);
					else
						parent_node = get_small_parent(q, new_stree_ref);

					if (parent_node > -1) {
						cover_tid[parent_node] += cover_tid[q];
						cnt16[parent_node][i%CNT] += cover_tid[q];
						//	sum[parent_node] += cover_new[q];
						//	cover_record[i][parent_node] += cover_record[i][q];
					}
				}
				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];
					cover_tid[q] = 0;
				}

				sam_this_round++;
			}

			//cout << ranking << "\t" << sam_this_round << endl;


			adding_time += (GetCurrentTimeSec() - adding_acc_time);
			//cout << ranking << ":  beta:" << beta << " beta * labeled_arcs " << beta * labeled_arcs << " new_tree_arcs " << new_tree_arcs << " how many " << how_many << endl;

			//double u3_acc_time = GetCurrentTimeSec();

			updating_acc_time = GetCurrentTimeSec();

			// Build the cover value heap for one time.
			if (switch_small == true && value_heap.empty())
				for (int nid = 0; nid < numOfVertices; ++nid)
					if (has_picked[nid] == false)
						value_heap.update(nid, -picked_value[nid]);

			//	cnt16_2max(picked_value, sum, has_picked, affected, cnt16, switch_small, value_heap);		
			//double u3_acc_time = GetCurrentTimeSec();

			//			cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap, best_sum_vid, lazy_updates);
			cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);
			//u3_time += (GetCurrentTimeSec() - u3_acc_time);

			if (affected_count < numOfVertices / 8)
				switch_set = true;
			else
				switch_set = false;
			affected_count = 0;

			//iteration_time[ranking] = (GetCurrentTimeSec() - updating_acc_time);
			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			time_this_iteration = GetCurrentTimeSec() - time_this_iteration;
			iteration_time[ranking] = time_this_iteration;
		}


		cout << "picked count:" << picked_count << endl;
		cout << "picked as tree covers:" << used_as_cover << endl;
		

		if (stop_flag == true) return;

		double labeling_acc_time = GetCurrentTimeSec();
		
		for (NodeID v = 0; v < numOfVertices; ++v) {
			if (DIRECTED_FLAG == false) {
				NodeID k = tmp_idx[v].first.size();
				index_[v].spt_v = (NodeID*)memalign(64, k * sizeof(NodeID));
				index_[v].spt_d = (EdgeWeight*)memalign(64, k * sizeof(EdgeWeight));
				for (NodeID i = 0; i < k; ++i) index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) index_[v].spt_d[i] = tmp_idx[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();

				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();

			}
			else {
				NodeID k = tmp_idx[v].first.size();
				//index_[inv[v]].spt_v.resize(k);
				//index_[inv[v]].spt_d.resize(k);
				index_[v].spt_v = (NodeID*)memalign(64, k * sizeof(NodeID));
				index_[v].spt_d = (EdgeWeight*)memalign(64, k * sizeof(EdgeWeight));

				for (NodeID i = 0; i < k; ++i) index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) index_[v].spt_d[i] = tmp_idx[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();
				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();


				k = r_tmp_idx[v].first.size();
				//index_[inv[v]].spt_v.resize(k);
				//index_[inv[v]].spt_d.resize(k);
				bindex_[v].spt_v = (NodeID*)memalign(64, k * sizeof(NodeID));
				bindex_[v].spt_d = (EdgeWeight*)memalign(64, k * sizeof(EdgeWeight));

				for (NodeID i = 0; i < k; ++i) bindex_[v].spt_v[i] = r_tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) bindex_[v].spt_d[i] = r_tmp_idx[v].second[i];
				r_tmp_idx[v].first.clear();
				r_tmp_idx[v].second.clear();
				r_tmp_idx[v].first.shrink_to_fit();
				r_tmp_idx[v].second.shrink_to_fit();
			}
		}

		num_of_trees = st_roots.size();

		double small_tree_num = 0;

		for (int i = 0; i < num_of_trees; ++i) {
			if (trees_pointers[i] != -1) {
				small_tree_num++;
				//small_tree_num += strees[trees_pointers[i]].size();
			}
		}

		total_resources = small_tree_num;
	}

	~BP_Betweenness_Ordering() {
		clear();
	}

	void clear() {
	}

};

class Given_Ordering : public Ordering{
public:
	Given_Ordering(char* order_file, Graph& graph) {
		
		inv.resize(numOfVertices);
		rank.resize(numOfVertices);
		
		ifstream ifs(order_file);

		for (NodeID v = 0; v < numOfVertices; ++v) {
			ifs >> inv[v];
		}


		ifs.close();

		Relabel(graph);
	}

	Given_Ordering(char* order_file, WGraph& wgraph) {

		inv.resize(numOfVertices);
		rank.resize(numOfVertices);

		ifstream ifs(order_file);

		for (NodeID v = 0; v < numOfVertices; ++v) {
			ifs >> inv[v];
		}


		ifs.close();
		rank.resize(inv.size());
		Relabel(wgraph);
	}

	Given_Ordering(char* order_file) {

		inv.resize(0);

		ifstream ifs(order_file);

		NodeID inv_id;
		while (ifs >> inv_id) {
			inv.push_back(inv_id);
		}
		ifs.close();

		rank.resize(inv.size());
		for (NodeID v = 0; v < inv.size(); ++v) rank[inv[v]] = v;
	}


};

#endif
