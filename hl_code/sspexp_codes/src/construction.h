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
#ifndef CONSTRUCTION_H
#define CONSTRUCTION_H

#include <queue>
#include <set>
#include "graph.h"
#include "paras.h"
#include "labels.h"
#include "ordering.h"
#include "heap.h"
#include "graph_search.h"
#include <omp.h>
#include <unordered_map>

#define numOfVertices SP_Constants::numOfVertices
#define numOfEdges SP_Constants::numOfEdges
#define INF_WEIGHT SP_Constants::INF_WEIGHT


class construction {
public:
	Label labels;
	DLabel dlabels;
	PLabel plabels;
	DPLabel dplabels;
	CLabel clabels;
	Ordering orders;
};


class PL : public construction {
public:
	
	vector<double> iteration_generated;
	vector<double> pruning_power;

	PL(Graph &graph, Ordering &orders) {
		
		iteration_generated.resize(numOfVertices);
		pruning_power.resize(numOfVertices);

		vector<index_t>& index_ = labels.index_;
		vector<NodeID> &inv = orders.inv;
		vector<NodeID> &rank = orders.rank;
		//vector<vector<NodeID> > &adj = graph.adj;

		vector<bool> usd(numOfVertices, false);

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<bool> vis(numOfVertices);
		vector<NodeID> que(numOfVertices);
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);

		for (size_t r = 0; r < numOfVertices; ++r) {
			if (usd[r]) continue;

			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[r];

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}

			NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
			que[que_h++] = r;
			vis[r] = true;
			que_t1 = que_h;

			for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID v = que[que_i];
					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
					index_t &idx_v = index_[inv[v]];

					if (usd[v]) continue;
					for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
						NodeID w = tmp_idx_v.first[i];
						EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
						if (td <= d) {
							pruning_power[w]++;
							goto pruned;
						}
					}

					// Traverse
					tmp_idx_v.first.back() = r;
					tmp_idx_v.second.back() = d;
					tmp_idx_v.first.push_back(numOfVertices);
					tmp_idx_v.second.push_back(INF_WEIGHT);
					iteration_generated[r]++;

					/*for (size_t i = 0; i < adj[v].size(); ++i) {
						NodeID w = adj[v][i];*/
					for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid){
						NodeID w = graph.edges[eid];
						if (!vis[w]) {
							que[que_h++] = w;
							vis[w] = true;
						}
					}
					pruned:
						{}
				}
				que_t0 = que_t1;
				que_t1 = que_h;
			}
			
						
			for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
			usd[r] = true;
		}

		for (size_t v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			index_[inv[v]].spt_v.resize(k);
			index_[inv[v]].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_d[i] = tmp_idx[v].second[i];
			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();
			tmp_idx[v].first.shrink_to_fit();
			tmp_idx[v].second.shrink_to_fit();
		}
	}

	PL(Graph &graph, Ordering &orders, bool D_FLAGS) {

		iteration_generated.resize(numOfVertices);
		pruning_power.resize(numOfVertices);

		vector<index_t>& index_ = dlabels.index_;
		vector<index_t>& bindex_ = dlabels.bindex_;
				
		vector<NodeID> &inv = orders.inv;
		vector<NodeID> &rank = orders.rank;
	/*	vector<vector<NodeID> > &adj = graph.adj;
		vector<vector<NodeID> > &r_adj = graph.r_adj;*/

		// Array Representation
		vector<EdgeID>& vertices = graph.vertices;
		vector<EdgeID>& r_vertices = graph.r_vertices;
		vector<NodeID>& edges = graph.edges;
		vector<NodeID>& r_edges = graph.r_edges;

		vector<bool> usd(numOfVertices, false);

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT))); // Backward labels.
		
		vector<bool> vis(numOfVertices);
		vector<NodeID> que(numOfVertices);
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT); // Forward labels of root.
		vector<EdgeWeight> r_dst_r(numOfVertices + 1, INF_WEIGHT); // Backward labels of root.

		for (size_t r = 0; r < numOfVertices; ++r) {
			if (usd[r]) continue;

			// Forward search.
			// Initialize forward labels of r.
			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[r];

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}

			NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
			que[que_h++] = r;
			vis[r] = true;
			que_t1 = que_h;

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
							pruning_power[w]++;
							goto pruned_forward;
						}
					}

					// Traverse
					r_tmp_idx_v.first.back() = r;
					r_tmp_idx_v.second.back() = d;
					r_tmp_idx_v.first.push_back(numOfVertices);
					r_tmp_idx_v.second.push_back(INF_WEIGHT);
					iteration_generated[r]++;
					/*for (size_t i = 0; i < adj[v].size(); ++i) {
						NodeID w = adj[v][i];*/
					// Array Representation
					for (EdgeID eid = vertices[v]; eid < vertices[v + 1]; ++eid){
						NodeID w = edges[eid];
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
			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;


			// Backward search.
			// Initialize backward labels of r.
			const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[r];

			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i) {
				r_dst_r[r_tmp_idx_r.first[i]] = r_tmp_idx_r.second[i];
			}


			que_t0 = 0, que_t1 = 0, que_h = 0;
			que[que_h++] = r;
			vis[r] = true;
			que_t1 = que_h;

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
							pruning_power[w]++;
							goto pruned_backward;
						}
					}

					// Traverse
					tmp_idx_v.first.back() = r;
					tmp_idx_v.second.back() = d;
					tmp_idx_v.first.push_back(numOfVertices);
					tmp_idx_v.second.push_back(INF_WEIGHT);
					iteration_generated[r]++;

					/*for (size_t i = 0; i < r_adj[v].size(); ++i) {
						NodeID w = r_adj[v][i];*/

					// Array Representation
					for (EdgeID eid = r_vertices[v]; eid < r_vertices[v + 1]; ++eid) {
						NodeID w = r_edges[eid];

						if (!vis[w]) {
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

			usd[r] = true;
		}

		for (size_t v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			index_[inv[v]].spt_v.resize(k);
			index_[inv[v]].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_d[i] = tmp_idx[v].second[i];
			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();

			tmp_idx[v].first.shrink_to_fit();
			tmp_idx[v].second.shrink_to_fit();

			k = r_tmp_idx[v].first.size();
			bindex_[inv[v]].spt_v.resize(k);
			bindex_[inv[v]].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) bindex_[inv[v]].spt_v[i] = r_tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) bindex_[inv[v]].spt_d[i] = r_tmp_idx[v].second[i];
			r_tmp_idx[v].first.clear();
			r_tmp_idx[v].second.clear();

			r_tmp_idx[v].first.shrink_to_fit();
			r_tmp_idx[v].second.shrink_to_fit();
		}


	}
	
	void TD_BP_UNDIRECTED(Graph& graph, Ordering &orderes, int kNumBitParallelRoots, bool directed = false, bool bp = true) {
		iteration_generated.resize(numOfVertices);
		pruning_power.resize(numOfVertices);

		vector<index_t>& index_ = labels.index_;
		vector<NodeID> &inv = orders.inv;
		vector<NodeID> &rank = orders.rank;
		//vector<vector<NodeID> > &adj = graph.adj;

		vector<bool> usd(numOfVertices, false);

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<bool> vis(numOfVertices);
		vector<NodeID> que(numOfVertices);
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);

		for (size_t r = 0; r < numOfVertices; ++r) {
			if (usd[r]) continue;

			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[r];

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}

			NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
			que[que_h++] = r;
			vis[r] = true;
			que_t1 = que_h;

			for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID v = que[que_i];
					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
					index_t &idx_v = index_[inv[v]];

					if (usd[v]) continue;
					for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
						NodeID w = tmp_idx_v.first[i];
						EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
						if (td <= d) {
							pruning_power[w]++;
							goto pruned;
						}
					}

					// Traverse
					tmp_idx_v.first.back() = r;
					tmp_idx_v.second.back() = d;
					tmp_idx_v.first.push_back(numOfVertices);
					tmp_idx_v.second.push_back(INF_WEIGHT);
					iteration_generated[r]++;

					/*for (size_t i = 0; i < adj[v].size(); ++i) {
					NodeID w = adj[v][i];*/
					for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid) {
						NodeID w = graph.edges[eid];
						if (!vis[w]) {
							que[que_h++] = w;
							vis[w] = true;
						}
					}
				pruned:
					{}
				}
				que_t0 = que_t1;
				que_t1 = que_h;
			}


			for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
			usd[r] = true;
		}

		for (size_t v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			index_[inv[v]].spt_v.resize(k);
			index_[inv[v]].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_d[i] = tmp_idx[v].second[i];
			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();
			tmp_idx[v].first.shrink_to_fit();
			tmp_idx[v].second.shrink_to_fit();
		}

	}

	void TP_path(Graph &graph, Ordering &orders) {

		iteration_generated.resize(numOfVertices);
		pruning_power.resize(numOfVertices);

		vector<index_t_path>& index_ = plabels.index_;
		vector<NodeID> &inv = orders.inv;
		vector<NodeID> &rank = orders.rank;
		//vector<vector<NodeID> > &adj = graph.adj;

		vector<bool> usd(numOfVertices, false);
		vector<NodeID> parents(numOfVertices, numOfVertices);

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<pair<vector<NodeID>, vector<NodeID> > >
			tmp_idx_parents(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<NodeID>(1, numOfVertices)));


		vector<bool> vis(numOfVertices);
		vector<NodeID> que(numOfVertices);
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);

		for (size_t r = 0; r < numOfVertices; ++r) {
			if (usd[r]) continue;

			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[r];

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}

			NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
			que[que_h++] = r;
			vis[r] = true;
			que_t1 = que_h;
			parents[r] = inv[r];

			for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID v = que[que_i];
					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
					pair<vector<NodeID>, vector<NodeID> > &tmp_idx_parent_v = tmp_idx_parents[v];
					index_t_path &idx_v = index_[inv[v]];

					if (usd[v]) continue;
					for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
						NodeID w = tmp_idx_v.first[i];
						EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
						if (td <= d) {
							pruning_power[w]++;
							goto pruned;
						}
					}

					// Traverse
					tmp_idx_v.first.back() = r;
					tmp_idx_v.second.back() = d;
					tmp_idx_v.first.push_back(numOfVertices);
					tmp_idx_v.second.push_back(INF_WEIGHT);
					
					tmp_idx_parent_v.first.back() = r;
					tmp_idx_parent_v.second.back() = parents[v];
					tmp_idx_parent_v.first.push_back(numOfVertices);
					tmp_idx_parent_v.second.push_back(numOfVertices);
					
					iteration_generated[r]++;



					/*for (size_t i = 0; i < adj[v].size(); ++i) {
					NodeID w = adj[v][i];*/
					for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid) {
						NodeID w = graph.edges[eid];
						if (!vis[w]) {
							parents[w] = inv[v];
							que[que_h++] = w;
							vis[w] = true;
						}
					}
				pruned:
					{}
				}
				que_t0 = que_t1;
				que_t1 = que_h;
			}


			for (size_t i = 0; i < que_h; ++i) {
				vis[que[i]] = false;
				parents[que[i]] = numOfVertices;
			}
			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
			usd[r] = true;
		}

		for (size_t v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			index_[inv[v]].spt_v.resize(k);
			index_[inv[v]].spt_d.resize(k);
			index_[inv[v]].spt_p.resize(k);
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_p[i] = tmp_idx_parents[v].second[i];
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_d[i] = tmp_idx[v].second[i];
			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();
			tmp_idx[v].first.shrink_to_fit();
			tmp_idx[v].second.shrink_to_fit();

			tmp_idx_parents[v].first.clear();
			tmp_idx_parents[v].second.clear();
			tmp_idx_parents[v].first.shrink_to_fit();
			tmp_idx_parents[v].second.shrink_to_fit();

			
		}
	}

	void TP_path_d(Graph &graph, Ordering &orders) {

		iteration_generated.resize(numOfVertices);
		pruning_power.resize(numOfVertices);

		vector<index_t_path>& index_ = dplabels.index_;
		vector<index_t_path>& bindex_ = dplabels.bindex_;

		vector<NodeID> &inv = orders.inv;
		vector<NodeID> &rank = orders.rank;
		/*	vector<vector<NodeID> > &adj = graph.adj;
		vector<vector<NodeID> > &r_adj = graph.r_adj;*/

		// Array Representation
		vector<EdgeID>& vertices = graph.vertices;
		vector<EdgeID>& r_vertices = graph.r_vertices;
		vector<NodeID>& edges = graph.edges;
		vector<NodeID>& r_edges = graph.r_edges;

		vector<bool> usd(numOfVertices, false);
		vector<NodeID> parents(numOfVertices, numOfVertices);		
		vector<NodeID> r_parents(numOfVertices, numOfVertices);


		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx_parent(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<pair<vector<NodeID>, vector<NodeID> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<NodeID>(1, numOfVertices))); // Backward labels.

		vector<pair<vector<NodeID>, vector<NodeID> > >
			r_tmp_idx_parent(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<NodeID>(1, numOfVertices))); // Backward labels.

		vector<bool> vis(numOfVertices);
		vector<NodeID> que(numOfVertices);
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT); // Forward labels of root.
		vector<EdgeWeight> r_dst_r(numOfVertices + 1, INF_WEIGHT); // Backward labels of root.

		for (size_t r = 0; r < numOfVertices; ++r) {
			if (usd[r]) continue;

			// Forward search.
			// Initialize forward labels of r.
			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[r];

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}

			NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
			que[que_h++] = r;
			vis[r] = true;
			que_t1 = que_h;
			parents[r] = inv[r];

			for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID v = que[que_i];
					pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];
					pair<vector<NodeID>, vector<NodeID> > &r_tmp_idx_parent_v = r_tmp_idx_parent[v];

					//index_t &idx_v = index_[inv[v]];

					if (usd[v]) continue;

					// Pruned by the forward labels of r and backward labels of v in the forward search from r when reaching v.
					for (size_t i = 0; i < r_tmp_idx_v.first.size(); ++i) {
						NodeID w = r_tmp_idx_v.first[i];
						EdgeWeight td = r_tmp_idx_v.second[i] + dst_r[w];
						if (td <= d) {
							pruning_power[w]++;
							goto pruned_forward;
						}
					}

					// Traverse
					r_tmp_idx_v.first.back() = r;
					r_tmp_idx_v.second.back() = d;
					r_tmp_idx_v.first.push_back(numOfVertices);
					r_tmp_idx_v.second.push_back(INF_WEIGHT);

					r_tmp_idx_parent_v.first.back() = r;
					r_tmp_idx_parent_v.second.back() = parents[v];
					r_tmp_idx_parent_v.first.push_back(numOfVertices);
					r_tmp_idx_parent_v.second.push_back(numOfVertices);

					iteration_generated[r]++;
					/*for (size_t i = 0; i < adj[v].size(); ++i) {
					NodeID w = adj[v][i];*/
					// Array Representation
					for (EdgeID eid = vertices[v]; eid < vertices[v + 1]; ++eid) {
						NodeID w = edges[eid];
						if (!vis[w]) {
							parents[w] = inv[v];
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
			for (size_t i = 0; i < que_h; ++i) {
				vis[que[i]] = false;
				parents[que[i]] = numOfVertices;
			}
			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;

			parents[r] = inv[r];

			// Backward search.
			// Initialize backward labels of r.
			const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[r];

			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i) {
				r_dst_r[r_tmp_idx_r.first[i]] = r_tmp_idx_r.second[i];
			}


			que_t0 = 0, que_t1 = 0, que_h = 0;
			que[que_h++] = r;
			vis[r] = true;
			que_t1 = que_h;

			for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID v = que[que_i];
					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
					pair<vector<NodeID>, vector<NodeID> > &tmp_idx_parent_v = tmp_idx_parent[v];
					//index_t &idx_v = index_[inv[v]];

					if (usd[v]) continue;

					// Pruned by the backward labels of r and forward labels of v in the backward search from r when reaching v (v->r path).
					for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
						NodeID w = tmp_idx_v.first[i];
						EdgeWeight td = tmp_idx_v.second[i] + r_dst_r[w];
						if (td <= d) {
							pruning_power[w]++;
							goto pruned_backward;
						}
					}

					// Traverse
					tmp_idx_v.first.back() = r;
					tmp_idx_v.second.back() = d;
					tmp_idx_v.first.push_back(numOfVertices);
					tmp_idx_v.second.push_back(INF_WEIGHT);
				
					iteration_generated[r]++;

					/*for (size_t i = 0; i < r_adj[v].size(); ++i) {
					NodeID w = r_adj[v][i];*/

					// Array Representation
					for (EdgeID eid = r_vertices[v]; eid < r_vertices[v + 1]; ++eid) {
						NodeID w = r_edges[eid];

						if (!vis[w]) {
							parents[w] = inv[v];
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
			for (size_t i = 0; i < que_h; ++i) {
				vis[que[i]] = false;
				parents[que[i]] = numOfVertices;
			}
			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i)
				r_dst_r[r_tmp_idx_r.first[i]] = INF_WEIGHT;

			usd[r] = true;
		}

		for (size_t v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			index_[inv[v]].spt_v.resize(k);
			index_[inv[v]].spt_d.resize(k);
			index_[inv[v]].spt_p.resize(k);
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_p[i] = tmp_idx_parent[v].second[i];
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_d[i] = tmp_idx[v].second[i];

			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();
			tmp_idx_parent[v].first.clear();
			tmp_idx_parent[v].second.clear();

			tmp_idx[v].first.shrink_to_fit();
			tmp_idx[v].second.shrink_to_fit();
			tmp_idx_parent[v].first.shrink_to_fit();
			tmp_idx_parent[v].second.shrink_to_fit();

			k = r_tmp_idx[v].first.size();
			bindex_[inv[v]].spt_v.resize(k);
			bindex_[inv[v]].spt_d.resize(k);
			bindex_[inv[v]].spt_p.resize(k);
			for (NodeID i = 0; i < k; ++i) bindex_[inv[v]].spt_v[i] = r_tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) bindex_[inv[v]].spt_p[i] = r_tmp_idx_parent[v].second[i];
			for (NodeID i = 0; i < k; ++i) bindex_[inv[v]].spt_d[i] = r_tmp_idx[v].second[i];

			r_tmp_idx[v].first.clear();
			r_tmp_idx[v].second.clear();

			r_tmp_idx_parent[v].first.clear();
			r_tmp_idx_parent[v].second.clear();

			r_tmp_idx[v].first.shrink_to_fit();
			r_tmp_idx[v].second.shrink_to_fit();
			r_tmp_idx_parent[v].first.shrink_to_fit();
			r_tmp_idx_parent[v].second.shrink_to_fit();
		}


	}

	PL(Graph &graph, Ordering &orders, bool path_flags, bool directed_flags) {
		if (path_flags == true) {

			if (directed_flags == true) {
				TP_path_d(graph, orders);
			}
			else {
				TP_path(graph, orders);
			}

		}
	}

/*	PL(Graph &graph, Ordering &orders, bool path_flags, bool directed_flags, bool bp_flags, int kNumBitParallelRoots) {
		if (path_flags == true) {

			if (directed_flags == true) {
				TP_path_d(graph, orders);
			}
			else {
				TP_path(graph, orders);
			}

		}
		else {
			if (directed_flags == false) {
				TP_BP(graph, orders, kNumBitParallelRoots);
			}

		}
	}
	*/


};

template<int kNumBitParallelRoots = 50>
class BPL {

	//typedef  BPLabel<kNumBitParallelRoots>::index_t_bp index_t_bp;

public:
	BPLabel<kNumBitParallelRoots> bplabels;
	DBPLabel<kNumBitParallelRoots> dbplabels;
	
	BPL(Graph &graph, Ordering &orders) {
		//bplabels = BPLabel(kNumBitParallelRoots);
		//kNumBitParallelRoots = 64;
		//bplabels.setParas(kNumBitParallelRoots);

	//	iteration_generated.resize(numOfVertices);
	//	pruning_power.resize(numOfVertices);


		index_t_bp<kNumBitParallelRoots>*& index_ = bplabels.index_bp;


		vector<NodeID> &inv = orders.inv;
		vector<NodeID> &rank = orders.rank;
		//vector<vector<NodeID> > &adj = graph.adj;

		vector<bool> usd(numOfVertices, false);

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<bool> vis(numOfVertices);
		vector<NodeID> que(numOfVertices);
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);


		{
			vector<EdgeWeight> tmp_d(numOfVertices);
			vector<std::pair<uint64_t, uint64_t> > tmp_s(numOfVertices);
			vector<NodeID> que(numOfVertices);
			vector<std::pair<NodeID, NodeID> > sibling_es(numOfEdges);
			vector<std::pair<NodeID, NodeID> > child_es(numOfEdges);

			cout << "Building BP labels" << endl;
			index_ = (index_t_bp<kNumBitParallelRoots>*)memalign(64, numOfVertices * sizeof(index_t_bp<kNumBitParallelRoots>));
			//for (NodeID v = 0; v < numOfVertices; ++v) {
			//	index_t_bp &idx = index_[v];
			//	idx.bpspt_d = (EdgeWeight*)memalign(64, kNumBitParallelRoots * sizeof(EdgeWeight));
			//	idx.bpspt_s = (uint64_t*)memalign(64, kNumBitParallelRoots * 2 * sizeof(uint64_t));
			//	/*for (int i = 0; i < kNumBitParallelRoots; ++i) {
			//	idx.bpspt_s[i] = (uint64_t*)memalign(64, 2 * sizeof(uint64_t));
			//	}*/
			//}

			int r = 0;
			for (int i_bpspt = 0; i_bpspt < kNumBitParallelRoots; ++i_bpspt) {
				while (r < numOfVertices && usd[r]) ++r;
				if (r == numOfVertices) {
					for (NodeID v = 0; v < numOfVertices; ++v) index_[v].bpspt_d[i_bpspt] = INF_WEIGHT;
					continue;
				}
				usd[r] = true;

				fill(tmp_d.begin(), tmp_d.end(), INF_WEIGHT);
				fill(tmp_s.begin(), tmp_s.end(), std::make_pair(0, 0));

				int que_t0 = 0, que_t1 = 0, que_h = 0;
				que[que_h++] = r;
				tmp_d[r] = 0;
				que_t1 = que_h;

				int ns = 0;
				vector<int> vs;

				vector<NodeID> adj_r(graph.vertices[r + 1] - graph.vertices[r]);
				for (EdgeID eid = graph.vertices[r]; eid < graph.vertices[r + 1]; eid++) {
					adj_r[eid - graph.vertices[r]] = graph.edges[eid];
				}

				sort(adj_r.begin(), adj_r.end());


				for (size_t i = 0; i < adj_r.size(); ++i) {
					NodeID v = adj_r[i];
					if (!usd[v]) {
						usd[v] = true;
						que[que_h++] = v;
						tmp_d[v] = 1;
						tmp_s[v].first = 1ULL << ns;
						vs.push_back(v);
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
					index_[inv[v]].bpspt_d[i_bpspt] = tmp_d[v];
					index_[inv[v]].bpspt_s[i_bpspt][0] = tmp_s[v].first;
					index_[inv[v]].bpspt_s[i_bpspt][1] = tmp_s[v].second & ~tmp_s[v].first;
				}
			}
		}


		cout << "Building normal labels" << endl;
		for (size_t r = 0; r < numOfVertices; ++r) {
			if (usd[r]) continue;

			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[r];
			index_t_bp<kNumBitParallelRoots> &idx_r = index_[inv[r]];

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}

			NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
			que[que_h++] = r;
			vis[r] = true;
			que_t1 = que_h;

			for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID v = que[que_i];
					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
					index_t_bp<kNumBitParallelRoots> &idx_v = index_[inv[v]];


					// Prefetch
					_mm_prefetch(&idx_v.bpspt_d[0], _MM_HINT_T0);
					_mm_prefetch(&idx_v.bpspt_s[0][0], _MM_HINT_T0);
					_mm_prefetch(&tmp_idx_v.first[0], _MM_HINT_T0);
					_mm_prefetch(&tmp_idx_v.second[0], _MM_HINT_T0);

					if (usd[v]) continue;

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
							if (td <= d) goto pruned;
						}
					}

					for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
						NodeID w = tmp_idx_v.first[i];
						EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
						if (td <= d) {
							//pruning_power[w]++;
							goto pruned;
						}
					}

					// Traverse
					tmp_idx_v.first.back() = r;
					tmp_idx_v.second.back() = d;
					tmp_idx_v.first.push_back(numOfVertices);
					tmp_idx_v.second.push_back(INF_WEIGHT);
					//iteration_generated[r]++;

					/*for (size_t i = 0; i < adj[v].size(); ++i) {
					NodeID w = adj[v][i];*/
					for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid) {
						NodeID w = graph.edges[eid];
						if (!vis[w]) {
							que[que_h++] = w;
							vis[w] = true;
						}
					}
				pruned:
					{}
				}
				que_t0 = que_t1;
				que_t1 = que_h;
			}


			for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
			usd[r] = true;
		}

		for (size_t v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			//index_[inv[v]].spt_v.resize(k);
			//index_[inv[v]].spt_d.resize(k);
			index_[inv[v]].spt_v = (NodeID*)memalign(64, k * sizeof(NodeID));
			index_[inv[v]].spt_d = (EdgeWeight*)memalign(64, k * sizeof(EdgeWeight));

			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_d[i] = tmp_idx[v].second[i];
			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();
			tmp_idx[v].first.shrink_to_fit();
			tmp_idx[v].second.shrink_to_fit();
		}
	}

	BPL(Graph &graph, Ordering &orders, bool directed_flags) {
		//bplabels = BPLabel(kNumBitParallelRoots);
		//kNumBitParallelRoots = 64;
		//bplabels.setParas(kNumBitParallelRoots);

		//	iteration_generated.resize(numOfVertices);
		//	pruning_power.resize(numOfVertices);

		if (directed_flags == false)
			return;

		index_t_bp<kNumBitParallelRoots>*& index_ = dbplabels.index_bp;
		index_t_bp<kNumBitParallelRoots>*& bindex_ = dbplabels.bindex_bp;


		vector<NodeID> &inv = orders.inv;
		vector<NodeID> &rank = orders.rank;
		//vector<vector<NodeID> > &adj = graph.adj;

		vector<bool> usd(numOfVertices, false);
		vector<bool> r_usd(numOfVertices, false); 
	//	vector<bool> bp_usd(numOfVertices, false);
	//	vector<bool> r_bp_usd(numOfVertices, false);

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<bool> vis(numOfVertices);
		vector<NodeID> que(numOfVertices);
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<EdgeWeight> r_dst_r(numOfVertices + 1, INF_WEIGHT); // Backward labels of root.


		{
			vector<EdgeWeight> tmp_d(numOfVertices);
			vector<std::pair<uint64_t, uint64_t> > tmp_s(numOfVertices);
			vector<EdgeWeight> r_tmp_d(numOfVertices);
			vector<std::pair<uint64_t, uint64_t> > r_tmp_s(numOfVertices);

			vector<NodeID> que(numOfVertices);
			vector<std::pair<NodeID, NodeID> > sibling_es(numOfEdges);
			vector<std::pair<NodeID, NodeID> > child_es(numOfEdges);
			vector<std::pair<NodeID, NodeID> > r_sibling_es(numOfEdges);
			vector<std::pair<NodeID, NodeID> > r_child_es(numOfEdges);

			cout << "Building BP labels" << endl;
			index_ = (index_t_bp<kNumBitParallelRoots>*)memalign(64, numOfVertices * sizeof(index_t_bp<kNumBitParallelRoots>));
			bindex_ = (index_t_bp<kNumBitParallelRoots>*)memalign(64, numOfVertices * sizeof(index_t_bp<kNumBitParallelRoots>));
			//for (NodeID v = 0; v < numOfVertices; ++v) {
			//	index_t_bp &idx = index_[v];
			//	idx.bpspt_d = (EdgeWeight*)memalign(64, kNumBitParallelRoots * sizeof(EdgeWeight));
			//	idx.bpspt_s = (uint64_t*)memalign(64, kNumBitParallelRoots * 2 * sizeof(uint64_t));
			//	/*for (int i = 0; i < kNumBitParallelRoots; ++i) {
			//	idx.bpspt_s[i] = (uint64_t*)memalign(64, 2 * sizeof(uint64_t));
			//	}*/
			//}

			int r = 0;
			for (int i_bpspt = 0; i_bpspt < kNumBitParallelRoots; ++i_bpspt) {
				while (r < numOfVertices && usd[r] ) ++r;
				if (r == numOfVertices) {
					for (NodeID v = 0; v < numOfVertices; ++v) {
						index_[v].bpspt_d[i_bpspt] = INF_WEIGHT;
						bindex_[v].bpspt_d[i_bpspt] = INF_WEIGHT;
					}
					continue;
				}

				r_usd[r] = true;
				//r_bp_usd[r] = true;

				//forward search
				fill(r_tmp_d.begin(), r_tmp_d.end(), INF_WEIGHT);
				fill(r_tmp_s.begin(), r_tmp_s.end(), std::make_pair(0, 0));

				int que_t0 = 0, que_t1 = 0, que_h = 0;
				que[que_h++] = r;
				r_tmp_d[r] = 0;
				que_t1 = que_h;

				int ns = 0;
				vector<int> vs;

				vector<NodeID> adj_r(graph.vertices[r + 1] - graph.vertices[r]);
				for (EdgeID eid = graph.vertices[r]; eid < graph.vertices[r + 1]; eid++) {
					adj_r[eid - graph.vertices[r]] = graph.edges[eid];
				}

				sort(adj_r.begin(), adj_r.end());

				vector<NodeID> r_adj_r(graph.r_vertices[r + 1] - graph.r_vertices[r]);
				for (EdgeID eid = graph.r_vertices[r]; eid < graph.r_vertices[r + 1]; eid++) {
					r_adj_r[eid - graph.r_vertices[r]] = graph.r_edges[eid];
				}

				sort(r_adj_r.begin(), r_adj_r.end());
				 
				vector<NodeID> common_adj;
				set_intersection(adj_r.begin(), adj_r.end(), r_adj_r.begin(), r_adj_r.end(), back_inserter(common_adj));
				sort(common_adj.begin(), common_adj.end());

				for (size_t i = 0; i < common_adj.size(); ++i) {
					NodeID v = common_adj[i];
					if (!r_usd[v]) {
						r_usd[v] = true;
						que[que_h++] = v;
						r_tmp_d[v] = 1;
						r_tmp_s[v].first = 1ULL << ns;
						vs.push_back(v);
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
								//}
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
						//r_tmp_s[v].second |= r_tmp_s[w].first;
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
					bindex_[inv[v]].bpspt_d[i_bpspt] = r_tmp_d[v];
					bindex_[inv[v]].bpspt_s[i_bpspt][0] = r_tmp_s[v].first;
					bindex_[inv[v]].bpspt_s[i_bpspt][1] = r_tmp_s[v].second & ~r_tmp_s[v].first;
				}
				//forward search end

				//backward

				usd[r] = true;
				fill(tmp_d.begin(), tmp_d.end(), INF_WEIGHT);
				fill(tmp_s.begin(), tmp_s.end(), std::make_pair(0, 0));

				que_t0 = 0, que_t1 = 0, que_h = 0;
				que[que_h++] = r;
				tmp_d[r] = 0;
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
					index_[inv[v]].bpspt_d[i_bpspt] = tmp_d[v];
					index_[inv[v]].bpspt_s[i_bpspt][0] = tmp_s[v].first;
					index_[inv[v]].bpspt_s[i_bpspt][1] = tmp_s[v].second & ~tmp_s[v].first;
				}
			}
		}


		cout << "Building normal labels" << endl;
		//for (size_t r = 0; r < numOfVertices; ++r) {
		//	if (usd[r]) continue;

		//	const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[r];
		//	index_t_bp<kNumBitParallelRoots> &idx_r = index_[inv[r]];

		//	for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
		//		dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
		//	}

		//	NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
		//	que[que_h++] = r;
		//	vis[r] = true;
		//	que_t1 = que_h;

		//	for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
		//		for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
		//			NodeID v = que[que_i];
		//			pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
		//			index_t_bp<kNumBitParallelRoots> &idx_v = index_[inv[v]];


		//			// Prefetch
		//			_mm_prefetch(&idx_v.bpspt_d[0], _MM_HINT_T0);
		//			_mm_prefetch(&idx_v.bpspt_s[0][0], _MM_HINT_T0);
		//			_mm_prefetch(&tmp_idx_v.first[0], _MM_HINT_T0);
		//			_mm_prefetch(&tmp_idx_v.second[0], _MM_HINT_T0);

		//			if (usd[v]) continue;

		//			for (int i = 0; i < kNumBitParallelRoots; ++i) {
		//				EdgeWeight td = idx_r.bpspt_d[i] + idx_v.bpspt_d[i];
		//				if (td - 2 <= d) {
		//					td +=
		//						(idx_r.bpspt_s[i][0] & idx_v.bpspt_s[i][0]) ? -2 :
		//						((idx_r.bpspt_s[i][0] & idx_v.bpspt_s[i][1]) |
		//							(idx_r.bpspt_s[i][1] & idx_v.bpspt_s[i][0]))
		//						? -1 : 0;
		//					if (td <= d) goto pruned;
		//				}
		//			}

		//			for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
		//				NodeID w = tmp_idx_v.first[i];
		//				EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
		//				if (td <= d) {
		//					//pruning_power[w]++;
		//					goto pruned;
		//				}
		//			}

		//			// Traverse
		//			tmp_idx_v.first.back() = r;
		//			tmp_idx_v.second.back() = d;
		//			tmp_idx_v.first.push_back(numOfVertices);
		//			tmp_idx_v.second.push_back(INF_WEIGHT);
		//			//iteration_generated[r]++;

		//			/*for (size_t i = 0; i < adj[v].size(); ++i) {
		//			NodeID w = adj[v][i];*/
		//			for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid) {
		//				NodeID w = graph.edges[eid];
		//				if (!vis[w]) {
		//					que[que_h++] = w;
		//					vis[w] = true;
		//				}
		//			}
		//		pruned:
		//			{}
		//		}
		//		que_t0 = que_t1;
		//		que_t1 = que_h;
		//	}


		//	for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
		//	for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
		//		dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
		//	usd[r] = true;
		//}

		for (size_t r = 0; r < numOfVertices; ++r) {
			if (usd[r]) continue;

			// Forward search.
			// Initialize forward labels of r.
			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[r];
			index_t_bp<kNumBitParallelRoots> &idx_r = index_[inv[r]];


			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}

			NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
			que[que_h++] = r;
			vis[r] = true;
			que_t1 = que_h;

			for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID v = que[que_i];
					pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];
					index_t_bp<kNumBitParallelRoots> &r_idx_v = bindex_[inv[v]];

					//index_t &idx_v = index_[inv[v]];
					// Prefetch
					_mm_prefetch(&r_idx_v.bpspt_d[0], _MM_HINT_T0);
					_mm_prefetch(&r_idx_v.bpspt_s[0][0], _MM_HINT_T0);
					_mm_prefetch(&r_tmp_idx_v.first[0], _MM_HINT_T0);
					_mm_prefetch(&r_tmp_idx_v.second[0], _MM_HINT_T0);

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
					r_tmp_idx_v.first.back() = r;
					r_tmp_idx_v.second.back() = d;
					r_tmp_idx_v.first.push_back(numOfVertices);
					r_tmp_idx_v.second.push_back(INF_WEIGHT);
					/*for (size_t i = 0; i < adj[v].size(); ++i) {
					NodeID w = adj[v][i];*/
					// Array Representation
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
			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;


			// Backward search.
			// Initialize backward labels of r.
			const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[r];
			index_t_bp<kNumBitParallelRoots> &r_idx_r = bindex_[inv[r]];

			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i) {
				r_dst_r[r_tmp_idx_r.first[i]] = r_tmp_idx_r.second[i];
			}


			que_t0 = 0, que_t1 = 0, que_h = 0;
			que[que_h++] = r;
			vis[r] = true;
			que_t1 = que_h;

			for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID v = que[que_i];
					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
					index_t_bp<kNumBitParallelRoots> &idx_v = index_[inv[v]];
					//index_t &idx_v = index_[inv[v]];
					_mm_prefetch(&idx_v.bpspt_d[0], _MM_HINT_T0);
					_mm_prefetch(&idx_v.bpspt_s[0][0], _MM_HINT_T0);
					_mm_prefetch(&tmp_idx_v.first[0], _MM_HINT_T0);
					_mm_prefetch(&tmp_idx_v.second[0], _MM_HINT_T0);

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
					tmp_idx_v.first.back() = r;
					tmp_idx_v.second.back() = d;
					tmp_idx_v.first.push_back(numOfVertices);
					tmp_idx_v.second.push_back(INF_WEIGHT);

					/*for (size_t i = 0; i < r_adj[v].size(); ++i) {
					NodeID w = r_adj[v][i];*/

					// Array Representation
					for (EdgeID eid = graph.r_vertices[v]; eid < graph.r_vertices[v + 1]; ++eid) {
						NodeID w = graph.r_edges[eid];

						if (!vis[w]) {
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

			usd[r] = true;
		}


		for (size_t v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			//index_[inv[v]].spt_v.resize(k);
			//index_[inv[v]].spt_d.resize(k);
			index_[inv[v]].spt_v = (NodeID*)memalign(64, k * sizeof(NodeID));
			index_[inv[v]].spt_d = (EdgeWeight*)memalign(64, k * sizeof(EdgeWeight));

			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_d[i] = tmp_idx[v].second[i];
			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();
			tmp_idx[v].first.shrink_to_fit();
			tmp_idx[v].second.shrink_to_fit();


			k = r_tmp_idx[v].first.size();
			//index_[inv[v]].spt_v.resize(k);
			//index_[inv[v]].spt_d.resize(k);
			bindex_[inv[v]].spt_v = (NodeID*)memalign(64, k * sizeof(NodeID));
			bindex_[inv[v]].spt_d = (EdgeWeight*)memalign(64, k * sizeof(EdgeWeight));

			for (NodeID i = 0; i < k; ++i) bindex_[inv[v]].spt_v[i] = r_tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) bindex_[inv[v]].spt_d[i] = r_tmp_idx[v].second[i];
			r_tmp_idx[v].first.clear();
			r_tmp_idx[v].second.clear();
			r_tmp_idx[v].first.shrink_to_fit();
			r_tmp_idx[v].second.shrink_to_fit();

		}
	
}
};


class PL_W : public construction {
public:

	vector<double> iteration_generated;
	vector<double> pruning_power;

	PL_W(WGraph &wgraph, Ordering &orders) {
 
		iteration_generated.resize(numOfVertices);
		pruning_power.resize(numOfVertices);

		vector<index_t>& index_ = labels.index_;
		vector<NodeID> &inv = orders.inv;
		vector<NodeID> &rank = orders.rank;
	//	vector<vector<NodeID> > &adj = wgraph.adj;
	//	vector<vector<EdgeWeight> > &adj_weight = wgraph.adj_weight;

		vector<EdgeID>& vertices = wgraph.vertices;
		vector<NodeEdgeWeightPair>& edges = wgraph.edges;

		vector<bool> usd(numOfVertices, false);

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<bool> vis(numOfVertices);
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		queue<NodeID> visited_que;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);


		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);
		long pop = 0;
		double hsize = 0;
		for (size_t r = 0; r < numOfVertices; ++r) {
			if (usd[r]) continue;

			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[r];

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}

			pqueue.update(r, 0);
			//vis[r] = true;

			long max_heap_size = 0;
			long heap_size = 1;
			while (!pqueue.empty()) {


				pop++;
				heap_size--;

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
						pruning_power[w]++;
						goto pruned;
					}
				}

				// Traverse
				tmp_idx_v.first.back() = r;
				tmp_idx_v.second.back() = v_d;
				tmp_idx_v.first.push_back(numOfVertices);
				tmp_idx_v.second.push_back(INF_WEIGHT);

				iteration_generated[r]++;

			//	for (size_t i = 0; i < adj[v].size(); ++i) {
			//		NodeID w = adj[v][i];
			//		EdgeWeight w_d = adj_weight[v][i] + v_d;

			for (EdgeID eid = vertices[v]; eid < vertices[v + 1]; ++eid) {
				NodeID w = edges[eid].first;
				EdgeWeight w_d = edges[eid].second + v_d;
					if (!vis[w]) {
						if (distances[w] == INF_WEIGHT) {
							heap_size++;
							if (max_heap_size < heap_size)
								max_heap_size = heap_size;
						}

						if( distances[w] > w_d ){
							pqueue.update(w, w_d);
							distances[w] = w_d;
						}
					}
				}
				pruned: 
					{}
			}
			hsize = hsize + max_heap_size;
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
			usd[r] = true;
		}

		//cout << "total pop:" << pop << endl;
		//cout << "heap size:" << (double)hsize / (double)numOfVertices << endl;
		
		double count = 0;
		for (size_t v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			count = count + k - 1;
			index_[inv[v]].spt_v.resize(k);
			index_[inv[v]].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_d[i] = tmp_idx[v].second[i];
			tmp_idx[v].first.clear(); 
			tmp_idx[v].second.clear();
			tmp_idx[v].first.shrink_to_fit();
			tmp_idx[v].second.shrink_to_fit();
		}
		cout << "Average Label Size:" << count / numOfVertices << endl;
		
	}

	PL_W(WGraph &wgraph, Ordering &orders, bool DIRECTED, bool PATH_QUERY) {
		
		// Generating Path Labels

		iteration_generated.resize(numOfVertices);
		pruning_power.resize(numOfVertices);

		vector<index_t_path>& index_ = plabels.index_;
		vector<NodeID> &inv = orders.inv;
		vector<NodeID> &rank = orders.rank;
		//	vector<vector<NodeID> > &adj = wgraph.adj;
		//	vector<vector<EdgeWeight> > &adj_weight = wgraph.adj_weight;

		vector<EdgeID>& vertices = wgraph.vertices;
		vector<NodeEdgeWeightPair>& edges = wgraph.edges;

		vector<bool> usd(numOfVertices, false);

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));
		vector<pair<vector<NodeID>, vector<NodeID> > > tmp_idx_parent(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices), vector<NodeID>(1, numOfVertices)));

		vector<bool> vis(numOfVertices);
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		queue<NodeID> visited_que;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
		vector<NodeID> parents(numOfVertices, numOfVertices);		

		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);

		for (size_t r = 0; r < numOfVertices; ++r) {
			if (usd[r]) continue;

			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[r];

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}

			parents[r] = inv[r]; 
			pqueue.update(r, 0);
			//vis[r] = true;

			while (!pqueue.empty()) {
				NodeID v;
				EdgeWeight v_d;
				pqueue.extract_min(v, v_d);
				pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
				pair<vector<NodeID>, vector<NodeID> > &tmp_idx_parent_v = tmp_idx_parent[v];

				vis[v] = true;
				visited_que.push(v);

				if (usd[v]) continue;
				for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
					NodeID w = tmp_idx_v.first[i];
					EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
					if (td <= v_d) {
						pruning_power[w]++;
						goto pruned;
					}
				}

				// Traverse
				tmp_idx_v.first.back() = r;
				tmp_idx_v.second.back() = v_d;
				tmp_idx_v.first.push_back(numOfVertices);
				tmp_idx_v.second.push_back(INF_WEIGHT);


				tmp_idx_parent_v.first.back() = r;
				tmp_idx_parent_v.second.back() = parents[v];
				tmp_idx_parent_v.first.push_back(numOfVertices);
				tmp_idx_parent_v.second.push_back(numOfVertices);

				iteration_generated[r]++;

				//	for (size_t i = 0; i < adj[v].size(); ++i) {
				//		NodeID w = adj[v][i];
				//		EdgeWeight w_d = adj_weight[v][i] + v_d;

				for (EdgeID eid = vertices[v]; eid < vertices[v + 1]; ++eid) {
					NodeID w = edges[eid].first;
					EdgeWeight w_d = edges[eid].second + v_d;
					if (!vis[w]) {
						if (distances[w] > w_d) {
							pqueue.update(w, w_d);
							distances[w] = w_d;
							parents[w] = inv[v];
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
				parents[vis_v] = numOfVertices;
				pqueue.clear(vis_v);
			}

			pqueue.clear_n();

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
			usd[r] = true;
		}

		for (size_t v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			index_[inv[v]].spt_v.resize(k);
			index_[inv[v]].spt_d.resize(k);
			index_[inv[v]].spt_p.resize(k);
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_d[i] = tmp_idx[v].second[i];
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_p[i] = tmp_idx_parent[v].second[i];
			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();
			tmp_idx_parent[v].first.clear();
			tmp_idx_parent[v].second.clear();
		}
	}


	PL_W(WGraph &wgraph, Ordering &orders, bool D_FLAGS) {

		iteration_generated.resize(numOfVertices);
		pruning_power.resize(numOfVertices);

		vector<index_t>& index_ = dlabels.index_;

		vector<NodeID> &inv = orders.inv;
		vector<NodeID> &rank = orders.rank;
		/*vector<vector<NodeID> > &adj = wgraph.adj;
		vector<vector<EdgeWeight> > &adj_weight = wgraph.adj_weight;*/

		vector<index_t>& bindex_ = dlabels.bindex_;/*
		vector<vector<NodeID> > &r_adj = wgraph.r_adj;
		vector<vector<EdgeWeight> > &r_adj_weight = wgraph.r_adj_weight;*/

		//Array Representation
		vector<EdgeID>& vertices = wgraph.vertices;
		vector<EdgeID>& r_vertices = wgraph.r_vertices;
		vector<NodeEdgeWeightPair>& edges = wgraph.edges;
		vector<NodeEdgeWeightPair>& r_edges = wgraph.r_edges;

		vector<bool> usd(numOfVertices, false);

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<bool> vis(numOfVertices);
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<EdgeWeight> r_dst_r(numOfVertices + 1, INF_WEIGHT);
		queue<NodeID> visited_que;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);


		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);

		// Forward search from r.
		for (size_t r = 0; r < numOfVertices; ++r) {
			if (usd[r]) continue;

			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[r];

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}

			pqueue.update(r, 0);
			//vis[r] = true;

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
						pruning_power[w]++;
						goto pruned_forward;
					}
				}

				// Traverse
				r_tmp_idx_v.first.back() = r;
				r_tmp_idx_v.second.back() = v_d;
				r_tmp_idx_v.first.push_back(numOfVertices);
				r_tmp_idx_v.second.push_back(INF_WEIGHT);
				iteration_generated[r]++;
				
			
			/*	for (size_t i = 0; i < adj[v].size(); ++i) {
					NodeID w = adj[v][i];
					EdgeWeight w_d = adj_weight[v][i] + v_d;*/

					//Array Representation

			for (EdgeID eid = vertices[v]; eid < vertices[v + 1]; ++eid) {
				NodeID w = edges[eid].first;
				EdgeWeight w_d = edges[eid].second + v_d;
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

			while (!visited_que.empty()) {
				NodeID vis_v = visited_que.front();
				visited_que.pop();
				vis[vis_v] = false;
				distances[vis_v] = INF_WEIGHT;
				//pqueue.clear(vis_v);
			}

			//pqueue.clear_n();

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;

			
			// Backward search from r.
			const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[r];

			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i) {
				r_dst_r[r_tmp_idx_r.first[i]] = r_tmp_idx_r.second[i];
			}

			pqueue.update(r, 0);

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
					EdgeWeight td = tmp_idx_v.second[i] + r_dst_r[w];
					if (td <= v_d) {
						pruning_power[w]++;
						goto pruned_backward;
					}
				}

				// Traverse
				tmp_idx_v.first.back() = r;
				tmp_idx_v.second.back() = v_d;
				tmp_idx_v.first.push_back(numOfVertices);
				tmp_idx_v.second.push_back(INF_WEIGHT);

				iteration_generated[r]++;

				/*for (size_t i = 0; i < r_adj[v].size(); ++i) {
					NodeID w = r_adj[v][i];
					EdgeWeight w_d = r_adj_weight[v][i] + v_d;*/

				//Array Representation
				for (EdgeID eid = r_vertices[v]; eid < r_vertices[v + 1]; ++eid) {
					NodeID w = r_edges[eid].first;
					EdgeWeight w_d = r_edges[eid].second + v_d;

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
				//pqueue.clear(vis_v);
			}

		//	pqueue.clear_n();

			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i)
				r_dst_r[r_tmp_idx_r.first[i]] = INF_WEIGHT;

			usd[r] = true;			
		}

		for (size_t v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			index_[inv[v]].spt_v.resize(k);
			index_[inv[v]].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_d[i] = tmp_idx[v].second[i];
			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();
			tmp_idx[v].first.shrink_to_fit();
			tmp_idx[v].second.shrink_to_fit();

			k = r_tmp_idx[v].first.size();
			bindex_[inv[v]].spt_v.resize(k);
			bindex_[inv[v]].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) bindex_[inv[v]].spt_v[i] = r_tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) bindex_[inv[v]].spt_d[i] = r_tmp_idx[v].second[i];
			r_tmp_idx[v].first.clear();
			r_tmp_idx[v].second.clear();		
			r_tmp_idx[v].first.shrink_to_fit();
			r_tmp_idx[v].second.shrink_to_fit();	
		}

	}
	
	PL_W(WGraph &wgraph, Ordering &orders, bool D_FLAGS, bool PATH_QUERY, bool dwpath) {

		iteration_generated.resize(numOfVertices);
		pruning_power.resize(numOfVertices);

		vector<index_t_path>& index_ = dplabels.index_;

		vector<NodeID> &inv = orders.inv;
		vector<NodeID> &rank = orders.rank;
		/*vector<vector<NodeID> > &adj = wgraph.adj;
		vector<vector<EdgeWeight> > &adj_weight = wgraph.adj_weight;*/

		vector<index_t_path>& bindex_ = dplabels.bindex_;/*
		vector<vector<NodeID> > &r_adj = wgraph.r_adj;
		vector<vector<EdgeWeight> > &r_adj_weight = wgraph.r_adj_weight;*/

		//Array Representation
		vector<EdgeID>& vertices = wgraph.vertices;
		vector<EdgeID>& r_vertices = wgraph.r_vertices;
		vector<NodeEdgeWeightPair>& edges = wgraph.edges;
		vector<NodeEdgeWeightPair>& r_edges = wgraph.r_edges;

		vector<bool> usd(numOfVertices, false);

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		
		vector<NodeID> parents(numOfVertices, numOfVertices);
		vector<pair<vector<NodeID>, vector<NodeID> > > tmp_idx_parent(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices), vector<NodeID>(1, numOfVertices)));
		
		
		vector<NodeID> r_parents(numOfVertices, numOfVertices);
		vector<pair<vector<NodeID>, vector<NodeID> > > r_tmp_idx_parent(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices), vector<NodeID>(1, numOfVertices)));
		
				
		vector<bool> vis(numOfVertices);
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<EdgeWeight> r_dst_r(numOfVertices + 1, INF_WEIGHT);
		queue<NodeID> visited_que;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);


		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);

		// Forward search from r.
		for (size_t r = 0; r < numOfVertices; ++r) {
			if (usd[r]) continue;

			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[r];

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}
			pqueue.update(r, 0);
			distances[r] = 0;
			parents[r] = inv[r];
			//vis[r] = true;

			while (!pqueue.empty()) {
				NodeID v;
				EdgeWeight v_d;
				pqueue.extract_min(v, v_d);
				pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];
				pair<vector<NodeID>, vector<NodeID> > &r_tmp_idx_parent_v = r_tmp_idx_parent[v];
				vis[v] = true;
				visited_que.push(v);

				if (usd[v]) continue;
				for (size_t i = 0; i < r_tmp_idx_v.first.size(); ++i) {
					NodeID w = r_tmp_idx_v.first[i];
					EdgeWeight td = r_tmp_idx_v.second[i] + dst_r[w];
					if (td <= v_d) {
						pruning_power[w]++;
						goto pruned_forward;
					}
				}

				// Traverse
				r_tmp_idx_v.first.back() = r;
				r_tmp_idx_v.second.back() = v_d;
				r_tmp_idx_v.first.push_back(numOfVertices);
				r_tmp_idx_v.second.push_back(INF_WEIGHT);
				iteration_generated[r]++;
				
				r_tmp_idx_parent_v.first.back() = r;
				r_tmp_idx_parent_v.second.back() = parents[v];
				r_tmp_idx_parent_v.first.push_back(numOfVertices);
				r_tmp_idx_parent_v.second.push_back(numOfVertices);

			
			/*	for (size_t i = 0; i < adj[v].size(); ++i) {
					NodeID w = adj[v][i];
					EdgeWeight w_d = adj_weight[v][i] + v_d;*/

					//Array Representation

			for (EdgeID eid = vertices[v]; eid < vertices[v + 1]; ++eid) {
				NodeID w = edges[eid].first;
				EdgeWeight w_d = edges[eid].second + v_d;
					if (!vis[w]) {
						if (distances[w] > w_d) {
							parents[w] = inv[v];
							pqueue.update(w, w_d);
							distances[w] = w_d;
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
				parents[vis_v] = numOfVertices;
				//pqueue.clear(vis_v);
			}

			//pqueue.clear_n();

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;

			
			// Backward search from r.
			const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[r];
			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i) {
				r_dst_r[r_tmp_idx_r.first[i]] = r_tmp_idx_r.second[i];
			}

			pqueue.update(r, 0);
			r_parents[r] = inv[r];


			while (!pqueue.empty()) {
				NodeID v;
				EdgeWeight v_d;
				pqueue.extract_min(v, v_d);
				pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
				pair<vector<NodeID>, vector<NodeID> > &tmp_idx_parent_v = tmp_idx_parent[v];
				vis[v] = true;
				visited_que.push(v);

				if (usd[v]) continue;
				for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
					NodeID w = tmp_idx_v.first[i];
					EdgeWeight td = tmp_idx_v.second[i] + r_dst_r[w];
					if (td <= v_d) {
						pruning_power[w]++;
						goto pruned_backward;
					}
				}

				// Traverse
				tmp_idx_v.first.back() = r;
				tmp_idx_v.second.back() = v_d;
				tmp_idx_v.first.push_back(numOfVertices);
				tmp_idx_v.second.push_back(INF_WEIGHT);

				tmp_idx_parent_v.first.back() = r;
				tmp_idx_parent_v.second.back() = r_parents[v];
				tmp_idx_parent_v.first.push_back(numOfVertices);
				tmp_idx_parent_v.second.push_back(numOfVertices);

			
				iteration_generated[r]++;

				/*for (size_t i = 0; i < r_adj[v].size(); ++i) {
					NodeID w = r_adj[v][i];
					EdgeWeight w_d = r_adj_weight[v][i] + v_d;*/

				//Array Representation
				for (EdgeID eid = r_vertices[v]; eid < r_vertices[v + 1]; ++eid) {
					NodeID w = r_edges[eid].first;
					EdgeWeight w_d = r_edges[eid].second + v_d;

					if (!vis[w]) {
						if (distances[w] > w_d) {
							pqueue.update(w, w_d);
							distances[w] = w_d;
							r_parents[w] = inv[v];
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
				r_parents[vis_v] = numOfVertices;
				//pqueue.clear(vis_v);
			}

		//	pqueue.clear_n();

			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i)
				r_dst_r[r_tmp_idx_r.first[i]] = INF_WEIGHT;

			usd[r] = true;			
		}

		for (size_t v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			index_[inv[v]].spt_v.resize(k);
			index_[inv[v]].spt_d.resize(k);
			index_[inv[v]].spt_p.resize(k);
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_d[i] = tmp_idx[v].second[i];
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_p[i] = tmp_idx_parent[v].second[i];
			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();
			tmp_idx[v].first.shrink_to_fit();
			tmp_idx[v].second.shrink_to_fit();

			tmp_idx_parent[v].first.clear();
			tmp_idx_parent[v].second.clear();
			tmp_idx_parent[v].first.shrink_to_fit();
			tmp_idx_parent[v].second.shrink_to_fit();

			k = r_tmp_idx[v].first.size();
			bindex_[inv[v]].spt_v.resize(k);
			bindex_[inv[v]].spt_d.resize(k);
			bindex_[inv[v]].spt_p.resize(k);
			for (NodeID i = 0; i < k; ++i) bindex_[inv[v]].spt_v[i] = r_tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) bindex_[inv[v]].spt_d[i] = r_tmp_idx[v].second[i];
			for (NodeID i = 0; i < k; ++i) bindex_[inv[v]].spt_p[i] = r_tmp_idx_parent[v].second[i];
		
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
		
};
	
	

class CPL : public construction {
public:
	
	vector<double> iteration_generated;
	vector<double> pruning_power;
	bool SECOND_LEVEL = false;
	long children_size;
	long r_children_size;


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

	CPL(Graph &graph, Ordering &orders, bool slevel) {
		SECOND_LEVEL = slevel;
		iteration_generated.resize(numOfVertices);
		pruning_power.resize(numOfVertices);

		vector<index_t>& index_ = labels.index_;
		vector<NodeID> &inv = orders.inv;
		vector<NodeID> &rank = orders.rank;
		//vector<vector<NodeID> > &adj = graph.adj;

		vector<bool> usd(numOfVertices, false);

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

	
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
		tmp_idx_token_parents(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
			vector<EdgeWeight>(1, numOfVertices)));

		vector<bool> vis(numOfVertices);
		vector<NodeID> que(numOfVertices);
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);

		for (size_t r = 0; r < numOfVertices; ++r) {
			if (usd[r]) continue;

			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[r];

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}

			NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
			que[que_h++] = r;
			vis[r] = true;
			que_t1 = que_h;

			for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID v = que[que_i];
					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_token_parents_v = tmp_idx_token_parents[v];
					index_t &idx_v = index_[inv[v]];

					if (usd[v]) continue;
					for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
						NodeID w = tmp_idx_v.first[i];
						EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
						if(tmp_idx_v.second[i] == d + dst_r[w] && tmp_idx_token_parents_v.first[i] == numOfVertices){
							tmp_idx_token_parents_v.first[i] = r;//tmp_idx_token_parents_v.first.size();
							tmp_idx_token_parents_v.second[i] = dst_r[w];
						}	
						if (td <= d) {
							pruning_power[w]++;
							goto pruned;
						}
					}

					// Traverse
					tmp_idx_v.first.back() = r;
					tmp_idx_v.second.back() = d;
					tmp_idx_v.first.push_back(numOfVertices);
					tmp_idx_v.second.push_back(INF_WEIGHT);

					tmp_idx_token_parents_v.first.back() = numOfVertices;
					tmp_idx_token_parents_v.second.back() = INF_WEIGHT;
					tmp_idx_token_parents_v.first.push_back(numOfVertices);
					tmp_idx_token_parents_v.second.push_back(INF_WEIGHT);

					iteration_generated[r]++;

					/*for (size_t i = 0; i < adj[v].size(); ++i) {
						NodeID w = adj[v][i];*/
					for (EdgeID eid = graph.vertices[v]; eid < graph.vertices[v + 1]; ++eid){
						NodeID w = graph.edges[eid];
						if (!vis[w]) {
							que[que_h++] = w;
							vis[w] = true;
						}
					}
					pruned:
						{}
				}
				que_t0 = que_t1;
				que_t1 = que_h;
			}
			
						
			for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
			usd[r] = true;
		}
		

		converttokens(tmp_idx,  tmp_idx_token_parents, false);

		vector<NodeID> change_anchor(numOfVertices);
		for (size_t v = 0; v < numOfVertices; ++v) {
			change_anchor[inv[v]] = clabels.anchor_p[v];
		}
		for (size_t v = 0; v < numOfVertices; ++v) {
			clabels.anchor_p[v] = change_anchor[v];
		}
		
/*
		for (size_t v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			index_[inv[v]].spt_v.resize(k);
			index_[inv[v]].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_d[i] = tmp_idx[v].second[i];
			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();
			tmp_idx[v].first.shrink_to_fit();
			tmp_idx[v].second.shrink_to_fit();
		}
		*/
	}

	CPL(Graph &graph, Ordering &orders, bool slevel, bool D_FLAGS) {
		SECOND_LEVEL = slevel;
		iteration_generated.resize(numOfVertices);
		pruning_power.resize(numOfVertices);

		vector<index_t>& index_ = dlabels.index_;
		vector<index_t>& bindex_ = dlabels.bindex_;
				
		vector<NodeID> &inv = orders.inv;
		vector<NodeID> &rank = orders.rank;
	/*	vector<vector<NodeID> > &adj = graph.adj;
		vector<vector<NodeID> > &r_adj = graph.r_adj;*/

		// Array Representation
		vector<EdgeID>& vertices = graph.vertices;
		vector<EdgeID>& r_vertices = graph.r_vertices;
		vector<NodeID>& edges = graph.edges;
		vector<NodeID>& r_edges = graph.r_edges;

		vector<bool> usd(numOfVertices, false);

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

		vector<bool> vis(numOfVertices);
		vector<NodeID> que(numOfVertices);
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT); // Forward labels of root.
		vector<EdgeWeight> r_dst_r(numOfVertices + 1, INF_WEIGHT); // Backward labels of root.

		for (size_t r = 0; r < numOfVertices; ++r) {
			if (usd[r]) continue;

			// Forward search.
			// Initialize forward labels of r.
			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[r];

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}
			
			const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[r];
			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i) {
				r_dst_r[r_tmp_idx_r.first[i]] = r_tmp_idx_r.second[i];
			}
			

			NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
			que[que_h++] = r;
			vis[r] = true;
			que_t1 = que_h;

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
							r_tmp_idx_token_parents_v.first[i] = r;//tmp_idx_token_parents_v.first.size();
							r_tmp_idx_token_parents_v.second[i] = r_dst_r[w];
						}
						if (td <= d) {
							pruning_power[w]++;
							goto pruned_forward;
						}
					}

					// Traverse
					r_tmp_idx_v.first.back() = r;
					r_tmp_idx_v.second.back() = d;
					r_tmp_idx_v.first.push_back(numOfVertices);
					r_tmp_idx_v.second.push_back(INF_WEIGHT);
					
					
					r_tmp_idx_token_parents_v.first.back() = numOfVertices;
					r_tmp_idx_token_parents_v.second.back() = INF_WEIGHT;
					r_tmp_idx_token_parents_v.first.push_back(numOfVertices);
					r_tmp_idx_token_parents_v.second.push_back(INF_WEIGHT);
					
					iteration_generated[r]++;
					/*for (size_t i = 0; i < adj[v].size(); ++i) {
						NodeID w = adj[v][i];*/
					// Array Representation
					for (EdgeID eid = vertices[v]; eid < vertices[v + 1]; ++eid){
						NodeID w = edges[eid];
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
			// Initialize backward labels of r.
			
 

			que_t0 = 0, que_t1 = 0, que_h = 0;
			que[que_h++] = r;
			vis[r] = true;
			que_t1 = que_h;

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
							tmp_idx_token_parents_v.first[i] = r;//tmp_idx_token_parents_v.first.size();
							tmp_idx_token_parents_v.second[i] = dst_r[w];
						}		
						if (td <= d) {
							pruning_power[w]++;
							goto pruned_backward;
						}
					}

					// Traverse
					tmp_idx_v.first.back() = r;
					tmp_idx_v.second.back() = d;
					tmp_idx_v.first.push_back(numOfVertices);
					tmp_idx_v.second.push_back(INF_WEIGHT);
					
					
					tmp_idx_token_parents_v.first.back() = numOfVertices;
					tmp_idx_token_parents_v.second.back() = INF_WEIGHT;
					tmp_idx_token_parents_v.first.push_back(numOfVertices);
					tmp_idx_token_parents_v.second.push_back(INF_WEIGHT);
					
					iteration_generated[r]++;

					/*for (size_t i = 0; i < r_adj[v].size(); ++i) {
						NodeID w = r_adj[v][i];*/

					// Array Representation
					for (EdgeID eid = r_vertices[v]; eid < r_vertices[v + 1]; ++eid) {
						NodeID w = r_edges[eid];

						if (!vis[w]) {
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
			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i)
				r_dst_r[r_tmp_idx_r.first[i]] = INF_WEIGHT;

			
			usd[r] = true;
		}
		
		
		converttokens(tmp_idx,  tmp_idx_token_parents, false);
		//converttokens(tmp_idx,  tmp_idx_token_parents, r_tmp_idx, r_tmp_idx_token_parents);
		cout << clabels.numOfTokens << " Tokens in total" << endl;
		cout << (double)children_size / (double) clabels.numOfTokens << " average children number" << endl;
		
		converttokens(r_tmp_idx,  r_tmp_idx_token_parents, true);
		cout << clabels.r_numOfTokens << " Tokens in total" << endl;
		cout << (double)r_children_size / (double) clabels.r_numOfTokens << " average children number" << endl;
		
		
		
		vector<NodeID> change_anchor(numOfVertices);
		for (size_t v = 0; v < numOfVertices; ++v) {
			change_anchor[inv[v]] = clabels.anchor_p[v];
		}
		for (size_t v = 0; v < numOfVertices; ++v) {
			clabels.anchor_p[v] = change_anchor[v];
		}
		
		for (size_t v = 0; v < numOfVertices; ++v) {
			change_anchor[inv[v]] = clabels.r_anchor_p[v];
		}
		for (size_t v = 0; v < numOfVertices; ++v) {
			clabels.r_anchor_p[v] = change_anchor[v];
		}
		
		
	}

};

class CPL_W : public construction {
public:

	vector<double> iteration_generated;
	vector<double> pruning_power;
	bool SECOND_LEVEL = false;
	long children_size;
	long r_children_size;
 
 

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

	

	CPL_W(WGraph &wgraph, Ordering &orders, bool slevel) {
 
		SECOND_LEVEL = slevel;
		iteration_generated.resize(numOfVertices);
		pruning_power.resize(numOfVertices);

		vector<index_t>& index_ = labels.index_;
		vector<NodeID> &inv = orders.inv;
		vector<NodeID> &rank = orders.rank;
	//	vector<vector<NodeID> > &adj = wgraph.adj;
	//	vector<vector<EdgeWeight> > &adj_weight = wgraph.adj_weight;

		vector<EdgeID>& vertices = wgraph.vertices;
		vector<NodeEdgeWeightPair>& edges = wgraph.edges;

		vector<bool> usd(numOfVertices, false);

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
		tmp_idx_token_parents(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
			vector<EdgeWeight>(1, numOfVertices)));

				
		vector<bool> vis(numOfVertices);
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		queue<NodeID> visited_que;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);


		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);
		long pop = 0;
		double hsize = 0;
		for (size_t r = 0; r < numOfVertices; ++r) {
			if (usd[r]) continue;

			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[r];

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}

			pqueue.update(r, 0);
			//vis[r] = true;

			long max_heap_size = 0;
			long heap_size = 1;
			while (!pqueue.empty()) {


				pop++;
				heap_size--;

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
					EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
					if(tmp_idx_v.second[i] == v_d + dst_r[w] && tmp_idx_token_parents_v.first[i] == numOfVertices){
							tmp_idx_token_parents_v.first[i] = r;//tmp_idx_token_parents_v.first.size();
							tmp_idx_token_parents_v.second[i] = dst_r[w];
					}
					if (td <= v_d) {
						pruning_power[w]++;
						goto pruned;
					}
				}

				// Traverse
				tmp_idx_v.first.back() = r;
				tmp_idx_v.second.back() = v_d;
				tmp_idx_v.first.push_back(numOfVertices);
				tmp_idx_v.second.push_back(INF_WEIGHT);

				tmp_idx_token_parents_v.first.back() = numOfVertices;
				tmp_idx_token_parents_v.second.back() = INF_WEIGHT;
				tmp_idx_token_parents_v.first.push_back(numOfVertices);
				tmp_idx_token_parents_v.second.push_back(INF_WEIGHT);

				
				iteration_generated[r]++;

			//	for (size_t i = 0; i < adj[v].size(); ++i) {
			//		NodeID w = adj[v][i];
			//		EdgeWeight w_d = adj_weight[v][i] + v_d;

			for (EdgeID eid = vertices[v]; eid < vertices[v + 1]; ++eid) {
				NodeID w = edges[eid].first;
				EdgeWeight w_d = edges[eid].second + v_d;
					if (!vis[w]) {
						if (distances[w] == INF_WEIGHT) {
							heap_size++;
							if (max_heap_size < heap_size)
								max_heap_size = heap_size;
						}

						if( distances[w] > w_d ){
							pqueue.update(w, w_d);
							distances[w] = w_d;
						}
					}
				}
				pruned: 
					{}
			}
			hsize = hsize + max_heap_size;
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
			usd[r] = true;
		}

		//cout << "total pop:" << pop << endl;
		//cout << "heap size:" << (double)hsize / (double)numOfVertices << endl;
		
		converttokens(tmp_idx,  tmp_idx_token_parents, false);

		vector<NodeID> change_anchor(numOfVertices);
		for (size_t v = 0; v < numOfVertices; ++v) {
			change_anchor[inv[v]] = clabels.anchor_p[v];
		}
		for (size_t v = 0; v < numOfVertices; ++v) {
			clabels.anchor_p[v] = change_anchor[v];
		}
		
		/*
		double count = 0;
		for (size_t v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			count = count + k - 1;
			index_[inv[v]].spt_v.resize(k);
			index_[inv[v]].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_d[i] = tmp_idx[v].second[i];
			tmp_idx[v].first.clear(); 
			tmp_idx[v].second.clear();
			tmp_idx[v].first.shrink_to_fit();
			tmp_idx[v].second.shrink_to_fit();
		}
		cout << "Average Label Size:" << count / numOfVertices << endl;
		*/
		
	}
	
	CPL_W(WGraph &wgraph, Ordering &orders, bool slevel, bool D_FLAGS) {
		SECOND_LEVEL = slevel;
		iteration_generated.resize(numOfVertices);
		pruning_power.resize(numOfVertices);

		vector<index_t>& index_ = dlabels.index_;

		vector<NodeID> &inv = orders.inv;
		vector<NodeID> &rank = orders.rank;
		/*vector<vector<NodeID> > &adj = wgraph.adj;
		vector<vector<EdgeWeight> > &adj_weight = wgraph.adj_weight;*/

		vector<index_t>& bindex_ = dlabels.bindex_;/*
		vector<vector<NodeID> > &r_adj = wgraph.r_adj;
		vector<vector<EdgeWeight> > &r_adj_weight = wgraph.r_adj_weight;*/

		//Array Representation
		vector<EdgeID>& vertices = wgraph.vertices;
		vector<EdgeID>& r_vertices = wgraph.r_vertices;
		vector<NodeEdgeWeightPair>& edges = wgraph.edges;
		vector<NodeEdgeWeightPair>& r_edges = wgraph.r_edges;

		vector<bool> usd(numOfVertices, false);

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

				
		vector<bool> vis(numOfVertices);
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<EdgeWeight> r_dst_r(numOfVertices + 1, INF_WEIGHT);
		queue<NodeID> visited_que;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);


		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);

		// Forward search from r.
		for (size_t r = 0; r < numOfVertices; ++r) {
			if (usd[r]) continue;

			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[r];

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}

			const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[r];

			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i) {
				r_dst_r[r_tmp_idx_r.first[i]] = r_tmp_idx_r.second[i];
			}

			
			pqueue.update(r, 0);
			//vis[r] = true;

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
						r_tmp_idx_token_parents_v.first[i] = r;//tmp_idx_token_parents_v.first.size();
						r_tmp_idx_token_parents_v.second[i] = r_dst_r[w];
					}
					if (td <= v_d) {
						pruning_power[w]++;
						goto pruned_forward;
					}
				}

				// Traverse
				r_tmp_idx_v.first.back() = r;
				r_tmp_idx_v.second.back() = v_d;
				r_tmp_idx_v.first.push_back(numOfVertices);
				r_tmp_idx_v.second.push_back(INF_WEIGHT);
				iteration_generated[r]++;
				
				r_tmp_idx_token_parents_v.first.back() = numOfVertices;
				r_tmp_idx_token_parents_v.second.back() = INF_WEIGHT;
				r_tmp_idx_token_parents_v.first.push_back(numOfVertices);
				r_tmp_idx_token_parents_v.second.push_back(INF_WEIGHT);
				
			/*	for (size_t i = 0; i < adj[v].size(); ++i) {
					NodeID w = adj[v][i];
					EdgeWeight w_d = adj_weight[v][i] + v_d;*/

					//Array Representation

			for (EdgeID eid = vertices[v]; eid < vertices[v + 1]; ++eid) {
				NodeID w = edges[eid].first;
				EdgeWeight w_d = edges[eid].second + v_d;
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

			while (!visited_que.empty()) {
				NodeID vis_v = visited_que.front();
				visited_que.pop();
				vis[vis_v] = false;
				distances[vis_v] = INF_WEIGHT;
				//pqueue.clear(vis_v);
			}

			//pqueue.clear_n();
			// Backward search from r.


			
			pqueue.update(r, 0);

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
					EdgeWeight td = tmp_idx_v.second[i] + r_dst_r[w];
					if(tmp_idx_v.second[i] == v_d + dst_r[w] && tmp_idx_token_parents_v.first[i] == numOfVertices){
						tmp_idx_token_parents_v.first[i] = r;//tmp_idx_token_parents_v.first.size();
						tmp_idx_token_parents_v.second[i] = dst_r[w];
					}	
					if (td <= v_d) {
						pruning_power[w]++;
						goto pruned_backward;
					}
				}

				// Traverse
				tmp_idx_v.first.back() = r;
				tmp_idx_v.second.back() = v_d;
				tmp_idx_v.first.push_back(numOfVertices);
				tmp_idx_v.second.push_back(INF_WEIGHT);

				tmp_idx_token_parents_v.first.back() = numOfVertices;
				tmp_idx_token_parents_v.second.back() = INF_WEIGHT;
				tmp_idx_token_parents_v.first.push_back(numOfVertices);
				tmp_idx_token_parents_v.second.push_back(INF_WEIGHT);
					
				iteration_generated[r]++;

				/*for (size_t i = 0; i < r_adj[v].size(); ++i) {
					NodeID w = r_adj[v][i];
					EdgeWeight w_d = r_adj_weight[v][i] + v_d;*/

				//Array Representation
				for (EdgeID eid = r_vertices[v]; eid < r_vertices[v + 1]; ++eid) {
					NodeID w = r_edges[eid].first;
					EdgeWeight w_d = r_edges[eid].second + v_d;

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
				//pqueue.clear(vis_v);
			}

		//	pqueue.clear_n();

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i)
				r_dst_r[r_tmp_idx_r.first[i]] = INF_WEIGHT;

			usd[r] = true;			
		}

		converttokens(tmp_idx,  tmp_idx_token_parents, false);
		//converttokens(tmp_idx,  tmp_idx_token_parents, r_tmp_idx, r_tmp_idx_token_parents);
		cout << clabels.numOfTokens << " Tokens in total" << endl;
		cout << (double)children_size / (double) clabels.numOfTokens << " average children number" << endl;
		
		converttokens(r_tmp_idx,  r_tmp_idx_token_parents, true);
		cout << clabels.r_numOfTokens << " Tokens in total" << endl;
		cout << (double)r_children_size / (double) clabels.r_numOfTokens << " average children number" << endl;
				
		
		vector<NodeID> change_anchor(numOfVertices);
		for (size_t v = 0; v < numOfVertices; ++v) {
			change_anchor[inv[v]] = clabels.anchor_p[v];
		}
		for (size_t v = 0; v < numOfVertices; ++v) {
			clabels.anchor_p[v] = change_anchor[v];
		}
		
		for (size_t v = 0; v < numOfVertices; ++v) {
			change_anchor[inv[v]] = clabels.r_anchor_p[v];
		}
		for (size_t v = 0; v < numOfVertices; ++v) {
			clabels.r_anchor_p[v] = change_anchor[v];
		}
		
		/*
		for (size_t v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			index_[inv[v]].spt_v.resize(k);
			index_[inv[v]].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_d[i] = tmp_idx[v].second[i];
			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();
			tmp_idx[v].first.shrink_to_fit();
			tmp_idx[v].second.shrink_to_fit();

			k = r_tmp_idx[v].first.size();
			bindex_[inv[v]].spt_v.resize(k);
			bindex_[inv[v]].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) bindex_[inv[v]].spt_v[i] = r_tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) bindex_[inv[v]].spt_d[i] = r_tmp_idx[v].second[i];
			r_tmp_idx[v].first.clear();
			r_tmp_idx[v].second.clear();		
			r_tmp_idx[v].first.shrink_to_fit();
			r_tmp_idx[v].second.shrink_to_fit();	
		}*/

	}
		
};
	
class Bottomup : public construction {


public:
	vector<double> iteration_generated;
	vector<double> pruning_power;

	typedef vector<pair<NodeID, EdgeWeight> >  Bucket;
	vector<bool> contracted;
	vector<Bucket> mtmBucket;
	vector<EdgeWeight> possibleWitness;

	int relabelByOrder(CHGraph &chgraph, Ordering &orders) {

		int totalEdge = 0;

		vector<NodeID>& inv = orders.inv;
		vector<NodeID>& rank = orders.rank;
		for (NodeID v = 0; v < numOfVertices; ++v) rank[inv[v]] = v;

		vector<vector<CHGraph::CH_Edge> > new_adj(numOfVertices);
		vector<vector<CHGraph::CH_Edge> > new_r_adj(numOfVertices);

		//for (int i = 0; i < numOfVertices; ++i) {
		//	for (int j = 0; j < chgraph.adj[i].size(); ++j) {
		//		if (j != chgraph.adj[i].size() - 1)
		//			if (chgraph.adj[i][j].target == chgraph.adj[i][j + 1].target) {
		//				cout << i << " hahaah:" << chgraph.adj[i][j].weight << "," << chgraph.adj[i][j + 1].weight << endl;
		//			}
		//	}
		//}

		for (NodeID v = 0; v < numOfVertices; ++v) {
			for (NodeID i = 0; i < chgraph.adj[v].size(); ++i)
				new_adj[rank[v]].push_back(CHGraph::CH_Edge(rank[chgraph.adj[v][i].target], 0, chgraph.adj[v][i].weight));

			totalEdge += chgraph.adj[v].size();

			if (DIRECTED_FLAG == true) {
				totalEdge += chgraph.r_adj[v].size();

				for (NodeID i = 0; i < chgraph.r_adj[v].size(); ++i)
					new_r_adj[rank[v]].push_back(CHGraph::CH_Edge(rank[chgraph.r_adj[v][i].target], 0, chgraph.r_adj[v][i].weight));
			}
		}
		chgraph.adj.swap(new_adj);
		if (DIRECTED_FLAG == true) {
			chgraph.r_adj.swap(new_r_adj);
		}

		for (int i = 0; i < numOfVertices; ++i) {
			if (DIRECTED_FLAG == true) {
				sort(chgraph.r_adj[i].begin(), chgraph.r_adj[i].end());
			}
		}

		new_adj.clear(); 
		if (DIRECTED_FLAG == true) {
			new_r_adj.clear();
		}

		return totalEdge;
	}

	int witness_search(NodeID v, CHGraph& chgraph, const int hopLimitsParameter, Ordering& orders, vector<bool>& vis, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, unordered_set<NodeID>& visited_set, vector<EdgeWeight>& distances, vector<bool>& hitTarget) {


		
		int addShortcuts = 0;

		vector<CHGraph::CH_Edge>& inEdges = chgraph.adj[v];
		vector<CHGraph::CH_Edge>& outEdges = chgraph.adj[v];

		vector<pair<NodeID, CHGraph::CH_Edge> > possibleShortcuts;

		// 1-hop witness search
		if	(hopLimitsParameter == 1) {			
			for (int i = 0; i < inEdges.size(); ++i) {
				//if (inEdges[i].level == V) continue; 
				NodeID inNode = inEdges[i].target;
				if (inNode >= v) continue; // Skip the nodes that have been already contracted.
				EdgeWeight inWeight = inEdges[i].weight;
				vector<CHGraph::CH_Edge>& outEdgesOfInNode = chgraph.adj[inNode];

				for (int j = 0; j < outEdges.size(); ++j) {
					//if (outEdges[j].level == V) continue;
					NodeID outNode = outEdges[j].target;
					if (outNode > v ) continue; // Skip the nodes that have been already contracted.
					if (outNode >= inNode) continue; // For undirected case, only test each pair once and we also skip the case inNode == outNode 
					EdgeWeight outWeight = outEdges[j].weight;

					EdgeWeight walkThroughWeight = inWeight + outWeight; // Distance for the path (inNode - v - outNode).

					bool foundWitness = false;
					for (int k = 0; k < outEdgesOfInNode.size(); ++k) {
						NodeID outNeighborOfInNode = outEdgesOfInNode[k].target;
					//	if (outNeighborLevelOfInNode == V) continue;
						if (outNeighborOfInNode >= v) continue; // Skip the ondes that have been already contracted.

						if (outNeighborOfInNode == outNode) {
							EdgeWeight outNeighborWeightOfInNode = outEdgesOfInNode[k].weight; // Distance for the direct path (inNode - outNode).



							if (outNeighborWeightOfInNode <= walkThroughWeight) {
								foundWitness = true;
								//walkThroughWeight = outNeighborWeightOfInNode;
							}
							break;
						}
					}
					if (foundWitness == false) {
						possibleShortcuts.push_back(make_pair(inNode, CHGraph::CH_Edge(outNode, v, walkThroughWeight)));
					}
				}
			}

		}

		// 2-hop witness search.
		if (hopLimitsParameter == 2) {

			// Init the many-to-many bucket.
			for (int j = 0; j < outEdges.size(); ++j) {
				NodeID outNode = outEdges[j].target;
				if (outNode >= v) continue; // Skip the nodes that have been already contracted.
				EdgeWeight outWeight = outEdges[j].weight;
				
				vector<CHGraph::CH_Edge>& inEdgesOfOutNode = chgraph.adj[outNode];
				for (int k = 0; k < inEdgesOfOutNode.size(); ++k) {
					NodeID inNeighborOfOutNode = inEdgesOfOutNode[k].target;
					NodeID inNeighborLevelOfOutNode = inEdgesOfOutNode[k].level;
					EdgeWeight inNeighborWeightOfOutNode = inEdgesOfOutNode[k].weight;

					if (inNeighborOfOutNode >= v) continue; // Skip the nodes that have been already contracted. And skip the current contracted node (this is important).
					mtmBucket[inNeighborOfOutNode].push_back(make_pair(outNode, inNeighborWeightOfOutNode));
				}
			}
			// Finish the init.
			
			for (int i = 0; i < inEdges.size(); ++i) {
				NodeID inNode = inEdges[i].target;
				if (inNode >= v) continue; // Skip the nodes that have bee	n already contracted.
				EdgeWeight inWeight = inEdges[i].weight;
				vector<CHGraph::CH_Edge>& outEdgesOfInNode = chgraph.adj[inNode];

				// One-hop witness search from the bucket of inNode.
				Bucket& bucketInNode = mtmBucket[inNode];
				for (int k = 0; k < bucketInNode.size(); ++k) {
					NodeID target = bucketInNode[k].first;
					EdgeWeight targetWeight = bucketInNode[k].second;
					if (target >= inNode) continue;

					if (possibleWitness[target] > targetWeight)
						possibleWitness[target] = targetWeight;
				}
				
				// Two-hop witness search from inNode.
				// 1-hop forward search from inNode and scan the buckets of the reached nodes to find the length of 2-hop witness.
				for (int k = 0; k < outEdgesOfInNode.size(); ++k) {
					NodeID reachNode = outEdgesOfInNode[k].target;
					if (reachNode >= v) continue;// Skip the nodes that have been already contracted.
					NodeID reachNodeLevel = outEdgesOfInNode[k].level;
					EdgeWeight reachNodeWeight = outEdgesOfInNode[k].weight;

					Bucket& bucketReachNode = mtmBucket[reachNode];
					for (int q = 0; q < bucketReachNode.size(); ++q) {
						NodeID target = bucketReachNode[q].first;
						if (target >= inNode) continue;
						EdgeWeight newTargetWeight = bucketReachNode[q].second + reachNodeWeight;
						if (possibleWitness[target] > newTargetWeight)
							possibleWitness[target] = newTargetWeight;
					}
				}

				// Scan the outNode of v, to check whether shortcuts are needed for (inNode - v - outNode).
				for (int j = 0; j < outEdges.size(); ++j) {
					NodeID outNode = outEdges[j].target;
					if (outNode > v) {
						possibleWitness[outNode] = INF_WEIGHT;
						continue;
					}
					if (outNode >= inNode) {
						possibleWitness[outNode] = INF_WEIGHT;
						continue;
					}// undirected, scan each pair only once.
					EdgeWeight outWeight = outEdges[j].weight;

					EdgeWeight witnessWeight = possibleWitness[outNode];
					possibleWitness[outNode] = INF_WEIGHT;
					
					EdgeWeight walkThroughWeight = inWeight + outWeight;

					if(witnessWeight > walkThroughWeight)
						possibleShortcuts.push_back(make_pair(inNode, CHGraph::CH_Edge(outNode, v, walkThroughWeight)));
				}
			}
			
			// Cleanup the many-to-many bucket.
			for (int j = 0; j < outEdges.size(); ++j) {
				NodeID outNode = outEdges[j].target;
				if (outNode > v) continue;
				EdgeWeight outWeight = outEdges[j].weight;

				vector<CHGraph::CH_Edge>& inEdgesOfOutNode = chgraph.adj[outNode];
				for (int k = 0; k < inEdgesOfOutNode.size(); ++k) {
					NodeID inNeighborOfOutNode = inEdgesOfOutNode[k].target;
					NodeID inNeighborLevelOfOutNode = inEdgesOfOutNode[k].level;
					EdgeWeight inNeighborWeightOfOutNode = inEdgesOfOutNode[k].weight;

					if (inNeighborOfOutNode >= v) continue;
					if(!mtmBucket[inNeighborOfOutNode].empty())
						mtmBucket[inNeighborOfOutNode].clear();
				}
			}
	
		}

		// Dijkstra Local Search.
		if (hopLimitsParameter > 2) {

			// Init the target table.
			int noOfTarget = 0;
			for (int j = 0; j < outEdges.size(); ++j) {
				NodeID outNode = outEdges[j].target;
				if (outNode >= v) continue; // Skip the nodes that have been already contracted.
								
				hitTarget[outNode] = true;
				noOfTarget++;
			}

			// Loop the incoming node for witness search.
			for (int i = 0; i < inEdges.size(); ++i) {
				//if (inEdges[i].level == V) continue; 
				NodeID inNode = inEdges[i].target;
				int visitedTarget = 0;
				if (inNode >= v) continue; // Skip the nodes that have been already contracted.
				EdgeWeight inWeight = inEdges[i].weight;
				vector<CHGraph::CH_Edge>& outEdgesOfInNode = chgraph.adj[inNode];

				pqueue.update(inNode, 0);
				visited_set.insert(inNode);
				distances[inNode] = 0;

				while (!pqueue.empty()) {
					NodeID u;
					EdgeWeight u_d;
					pqueue.extract_min(u, u_d);
					vis[u] = true;

					if (hitTarget[u] == true)
						visitedTarget++;
					if (visitedTarget == noOfTarget)
						break;

					for (int j = 0; j < chgraph.adj[u].size(); ++j) {
						NodeID w = chgraph.adj[u][j].target;
						EdgeWeight w_d = chgraph.adj[u][j].weight + u_d;
						if (w >= v) continue; // Can not visit contracted nodes.
						if (!vis[w]) {
							if (distances[w] > w_d) {
								pqueue.update(w, w_d);
								distances[w] = w_d;
								visited_set.insert(w);
							}
						}
					}
				}

				// Test witness for all outNode.
				for (int j = 0; j < outEdges.size(); ++j) {
					NodeID outNode = outEdges[j].target;
					if (outNode >= v) continue;
					if (outNode >= inNode) continue; // undirected, scan each pair only once.
					EdgeWeight outWeight = outEdges[j].weight;
					EdgeWeight walkThroughWeight = inWeight + outWeight;
					if (distances[outNode] > walkThroughWeight) {
						possibleShortcuts.push_back(make_pair(inNode, CHGraph::CH_Edge(outNode, v, walkThroughWeight)));
					}
				}

				// Clean up the dijkstra structures otherwise the trash will manipulate the next inNode's dijkstra search.
				for (unordered_set<NodeID>::iterator it = visited_set.begin(); it != visited_set.end(); ++it) {
					NodeID cv = *it;
					vis[cv] = false;
					distances[cv] = INF_WEIGHT;		
				}
				while (!pqueue.empty()) {
					NodeID tmpv;
					EdgeWeight tmpweight;
					pqueue.extract_min(tmpv, tmpweight);
				}
				visited_set.clear();

			}

			// Clean the target table for the next contracted node.
			for (int j = 0; j < outEdges.size(); ++j) {
				NodeID outNode = outEdges[j].target;
				if (outNode >= v) continue; // Skip the nodes that have been already contracted.
				hitTarget[outNode] = false;
			}
		}

		// Append the shortcuts.
		for (int i = 0; i < possibleShortcuts.size(); ++i) {
			NodeID fromNode = possibleShortcuts[i].first;
			NodeID toNode = possibleShortcuts[i].second.target;
			NodeID level = possibleShortcuts[i].second.level;
			EdgeWeight weight = possibleShortcuts[i].second.weight;
			
			addShortcuts++;
			addShortcuts++;

			int fromAdjSize = chgraph.adj[fromNode].size();
			bool skipfrom = false;
			for (int j = fromAdjSize - 1; j + 1 > 0; --j) {
				if (chgraph.adj[fromNode][j].target == toNode) {
					if (weight > chgraph.adj[fromNode][j].weight) break;
					chgraph.adj[fromNode][j].weight = weight;
					chgraph.adj[fromNode][j].level = level;
					skipfrom = true;
					addShortcuts--;
					break;
				}
			}
			if (!skipfrom) {
				chgraph.adj[fromNode].push_back(CHGraph::CH_Edge(toNode, level, weight));
			}
			
			int toAdjSize = chgraph.adj[toNode].size();
			bool skipto = false;
			for (int j = toAdjSize - 1; j + 1 > 0; --j) {
				if (chgraph.adj[toNode][j].target == fromNode) {
					if (weight > chgraph.adj[toNode][j].weight) break;
					chgraph.adj[toNode][j].weight = weight;
					chgraph.adj[toNode][j].level = level;
					skipto = true;
					addShortcuts--;
					break;
				}
			}

			if (!skipto)
				chgraph.adj[toNode].push_back(CHGraph::CH_Edge(fromNode, level, weight));

		}
		addShortcuts -= chgraph.adj[v].size(); // Two times because it will appear in others' adj.

		possibleShortcuts.clear();

		return addShortcuts;
	}
	
	int witness_search_directed(NodeID v, CHGraph& chgraph, const int hopLimitsParameter, Ordering& orders, vector<bool>& vis, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, unordered_set<NodeID>& visited_set, vector<EdgeWeight>& distances, vector<bool>& hitTarget) {
				
		int addShortcuts = 0;

		vector<CHGraph::CH_Edge>& inEdges = chgraph.r_adj[v];
		vector<CHGraph::CH_Edge>& outEdges = chgraph.adj[v];

		vector<pair<NodeID, CHGraph::CH_Edge> > possibleShortcuts;

		// 1-hop witness search
		if (hopLimitsParameter == 1) {
			for (int i = 0; i < inEdges.size(); ++i) {
				//if (inEdges[i].level == V) continue; 
				NodeID inNode = inEdges[i].target;
				if (inNode >= v) continue; // Skip the nodes that have been already contracted.
				EdgeWeight inWeight = inEdges[i].weight;
				vector<CHGraph::CH_Edge>& outEdgesOfInNode = chgraph.adj[inNode];

				for (int j = 0; j < outEdges.size(); ++j) {
					//if (outEdges[j].level == V) continue;
					NodeID outNode = outEdges[j].target;
					if (outNode >= v) continue; // Skip the nodes that have been already contracted.
					if (outNode == inNode) continue;
					EdgeWeight outWeight = outEdges[j].weight;

					EdgeWeight walkThroughWeight = inWeight + outWeight; // Distance for the path (inNode - v - outNode).

					bool foundWitness = false;
					for (int k = 0; k < outEdgesOfInNode.size(); ++k) {
						NodeID outNeighborOfInNode = outEdgesOfInNode[k].target;
						//	if (outNeighborLevelOfInNode == V) continue;
						if (outNeighborOfInNode >= v) continue; // Skip the ondes that have been already contracted.

						if (outNeighborOfInNode == outNode) {
							EdgeWeight outNeighborWeightOfInNode = outEdgesOfInNode[k].weight; // Distance for the direct path (inNode - outNode).
							
							if (outNeighborWeightOfInNode <= walkThroughWeight) {
								foundWitness = true;
								//walkThroughWeight = outNeighborWeightOfInNode;
							}
							break;
						}
					}
					if (foundWitness == false) {
						possibleShortcuts.push_back(make_pair(inNode, CHGraph::CH_Edge(outNode, v, walkThroughWeight)));
					}
				}
			}

		}

		// 2-hop witness search.
		if (hopLimitsParameter == 2) {
			// Init the many-to-many bucket.
			for (int j = 0; j < outEdges.size(); ++j) {
				NodeID outNode = outEdges[j].target;
				if (outNode >= v) continue; // Skip the nodes that have been already contracted.
				EdgeWeight outWeight = outEdges[j].weight;

				vector<CHGraph::CH_Edge>& inEdgesOfOutNode = chgraph.r_adj[outNode];
				for (int k = 0; k < inEdgesOfOutNode.size(); ++k) {
					NodeID inNeighborOfOutNode = inEdgesOfOutNode[k].target;
					NodeID inNeighborLevelOfOutNode = inEdgesOfOutNode[k].level;
					EdgeWeight inNeighborWeightOfOutNode = inEdgesOfOutNode[k].weight;

					if (inNeighborOfOutNode >= v) continue; // Skip the nodes that have been already contracted. And skip the current contracted node (this is important).
					mtmBucket[inNeighborOfOutNode].push_back(make_pair(outNode, inNeighborWeightOfOutNode));
				}
			}
			// Finish the init.

			for (int i = 0; i < inEdges.size(); ++i) {
				NodeID inNode = inEdges[i].target;
				if (inNode >= v) continue; // Skip the nodes that have bee	n already contracted.
				EdgeWeight inWeight = inEdges[i].weight;
				vector<CHGraph::CH_Edge>& outEdgesOfInNode = chgraph.adj[inNode];

				// One-hop witness search from the bucket of inNode.
				Bucket& bucketInNode = mtmBucket[inNode];
				for (int k = 0; k < bucketInNode.size(); ++k) {
					NodeID target = bucketInNode[k].first;
					EdgeWeight targetWeight = bucketInNode[k].second;
					//if (target >= inNode) continue;

					if (possibleWitness[target] > targetWeight)
						possibleWitness[target] = targetWeight;
				}

				// Two-hop witness search from inNode.
				// 1-hop forward search from inNode and scan the buckets of the reached nodes to find the length of 2-hop witness.
				for (int k = 0; k < outEdgesOfInNode.size(); ++k) {
					NodeID reachNode = outEdgesOfInNode[k].target;
					if (reachNode >= v) continue;// Skip the nodes that have been already contracted.
					NodeID reachNodeLevel = outEdgesOfInNode[k].level;
					EdgeWeight reachNodeWeight = outEdgesOfInNode[k].weight;

					Bucket& bucketReachNode = mtmBucket[reachNode];
					for (int q = 0; q < bucketReachNode.size(); ++q) {
						NodeID target = bucketReachNode[q].first;
					//	if (target >= inNode) continue;
						EdgeWeight newTargetWeight = bucketReachNode[q].second + reachNodeWeight;
						if (possibleWitness[target] > newTargetWeight)
							possibleWitness[target] = newTargetWeight;
					}
				}

				// Scan the outNode of v, to check whether shortcuts are needed for (inNode - v - outNode).
				for (int j = 0; j < outEdges.size(); ++j) {
					NodeID outNode = outEdges[j].target;
					if (outNode > v) {
						possibleWitness[outNode] = INF_WEIGHT;
						continue;
					}
					if (outNode == inNode) {
						possibleWitness[outNode] = INF_WEIGHT;
						continue;
					}
					EdgeWeight outWeight = outEdges[j].weight;

					EdgeWeight witnessWeight = possibleWitness[outNode];
					possibleWitness[outNode] = INF_WEIGHT;

					EdgeWeight walkThroughWeight = inWeight + outWeight;

					if (witnessWeight > walkThroughWeight)
						possibleShortcuts.push_back(make_pair(inNode, CHGraph::CH_Edge(outNode, v, walkThroughWeight)));
				}


			}

			// Cleanup the many-to-many bucket.
			for (int j = 0; j < outEdges.size(); ++j) {
				NodeID outNode = outEdges[j].target;
				if (outNode > v) continue;
				EdgeWeight outWeight = outEdges[j].weight;

				vector<CHGraph::CH_Edge>& inEdgesOfOutNode = chgraph.r_adj[outNode];
				for (int k = 0; k < inEdgesOfOutNode.size(); ++k) {
					NodeID inNeighborOfOutNode = inEdgesOfOutNode[k].target;
					NodeID inNeighborLevelOfOutNode = inEdgesOfOutNode[k].level;
					EdgeWeight inNeighborWeightOfOutNode = inEdgesOfOutNode[k].weight;

					if (inNeighborOfOutNode >= v) continue;
					if (!mtmBucket[inNeighborOfOutNode].empty())
						mtmBucket[inNeighborOfOutNode].clear();
				}
			}			
		}

		// Dijkstra Local Search.
		if (hopLimitsParameter > 2) {

			// Init the target table.
			int noOfTarget = 0;
			for (int j = 0; j < outEdges.size(); ++j) {
				NodeID outNode = outEdges[j].target;
				if (outNode >= v) continue; // Skip the nodes that have been already contracted.

				hitTarget[outNode] = true;
				noOfTarget++;
			}

			// Loop the incoming node for witness search.
			for (int i = 0; i < inEdges.size(); ++i) {
				//if (inEdges[i].level == V) continue; 
				NodeID inNode = inEdges[i].target;
				int visitedTarget = 0;
				if (inNode >= v) continue; // Skip the nodes that have been already contracted.
				EdgeWeight inWeight = inEdges[i].weight;
				vector<CHGraph::CH_Edge>& outEdgesOfInNode = chgraph.adj[inNode];

				pqueue.update(inNode, 0);
				visited_set.insert(inNode);
				distances[inNode] = 0;

				while (!pqueue.empty()) {
					NodeID u;
					EdgeWeight u_d;
					pqueue.extract_min(u, u_d);
					vis[u] = true;

					if (hitTarget[u] == true)
						visitedTarget++;
					if (visitedTarget == noOfTarget)
						break;

					for (int j = 0; j < chgraph.adj[u].size(); ++j) {
						NodeID w = chgraph.adj[u][j].target;
						EdgeWeight w_d = chgraph.adj[u][j].weight + u_d;
						if (w >= v) continue; // Can not visit contracted nodes.
						if (!vis[w]) {
							if (distances[w] > w_d) {
								pqueue.update(w, w_d);
								distances[w] = w_d;
								visited_set.insert(w);
							}
						}
					}
				}

				// Test witness for all outNode.
				for (int j = 0; j < outEdges.size(); ++j) {
					NodeID outNode = outEdges[j].target;
					if (outNode >= v) continue;
					if (outNode == inNode) continue;
					EdgeWeight outWeight = outEdges[j].weight;
					EdgeWeight walkThroughWeight = inWeight + outWeight;
					if (distances[outNode] > walkThroughWeight) {
						possibleShortcuts.push_back(make_pair(inNode, CHGraph::CH_Edge(outNode, v, walkThroughWeight)));
					}
				}

				// Clean up the dijkstra structures otherwise the trash will manipulate the next inNode's dijkstra search.
				for (unordered_set<NodeID>::iterator it = visited_set.begin(); it != visited_set.end(); ++it) {
					NodeID cv = *it;
					vis[cv] = false;
					distances[cv] = INF_WEIGHT;
				}
				while (!pqueue.empty()) {
					NodeID tmpv;
					EdgeWeight tmpweight;
					pqueue.extract_min(tmpv, tmpweight);
				}
				visited_set.clear();

			}

			// Clean the target table for the next contracted node.
			for (int j = 0; j < outEdges.size(); ++j) {
				NodeID outNode = outEdges[j].target;
				if (outNode >= v) continue; // Skip the nodes that have been already contracted.
				hitTarget[outNode] = false;
			}
		}

		// Append the shortcuts.
		for (int i = 0; i < possibleShortcuts.size(); ++i) {

			addShortcuts += 2;

			NodeID fromNode = possibleShortcuts[i].first;
			NodeID toNode = possibleShortcuts[i].second.target;
			NodeID level = possibleShortcuts[i].second.level;
			EdgeWeight weight = possibleShortcuts[i].second.weight;
			int fromAdjSize = chgraph.adj[fromNode].size();
			bool skipfrom = false;
			for (int j = fromAdjSize - 1; j + 1 > 0; --j) {
				if (chgraph.adj[fromNode][j].target == toNode) {
					if (weight > chgraph.adj[fromNode][j].weight) break;
					chgraph.adj[fromNode][j].weight = weight;
					chgraph.adj[fromNode][j].level = level;
					skipfrom = true;
					addShortcuts--;
					break;
				}
			}
			if (!skipfrom) {
				chgraph.adj[fromNode].push_back(CHGraph::CH_Edge(toNode, level, weight));
			}

			int toAdjSize = chgraph.r_adj[toNode].size();
			bool skipto = false;
			for (int j = toAdjSize - 1; j + 1 > 0; --j) {
				if (chgraph.r_adj[toNode][j].target == fromNode) {
					if (weight > chgraph.r_adj[toNode][j].weight) break;
					chgraph.r_adj[toNode][j].weight = weight;
					chgraph.r_adj[toNode][j].level = level;
					skipto = true;
					addShortcuts--;
					break;
				}
			}

			if (!skipto)
				chgraph.r_adj[toNode].push_back(CHGraph::CH_Edge(fromNode, level, weight));

		}
	
		addShortcuts -= 2 * chgraph.adj[v].size(); // 2 times because these out (in) edges will appear in others' in (out) edges.
		addShortcuts -= 2 * chgraph.r_adj[v].size();

		possibleShortcuts.clear();

		return addShortcuts;
	}

	int build_shortcuts_dij(NodeID source, vector<vector<CHGraph::CH_Edge> >& adj, vector<bool>& vis, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, unordered_set<NodeID>& visited_set, vector<EdgeWeight>& distances, vector<NodeID>& max_parents, vector<bool>& isNeighbor, vector<bool>& waitForPop, vector<vector<CHGraph::CH_Edge> >& shortcuts) {
		
		// All dsitances[] is INF_WEIGHT, all vis[] is false,  all max_parents[] is numOfVertices
		if (source == 0) return 0; // No shortcut for the most important vertex

		int uncover_count = source;
		pqueue.update(source, 0);
		distances[source] = 0;
		//isNeighbor[source] = true;
		visited_set.insert(source);

		int in_que_cover_wait_for_pop = -1;
		int in_que_cover = 0;
		int in_que = 1;
		/*if (source == 3)
			cout << "one search" << endl;*/
		//int sspace = 0;
		while (!pqueue.empty()) {
			NodeID v;
			EdgeWeight v_d;
			pqueue.extract_min(v, v_d);
			in_que--;
			vis[v] = true;
		//	sspace++;
		//	cout << "source: " << source << " visiting:" << v  << " inque:" << in_que << " inque_cover:" << in_que_cover << " inque_cover_wait_for_pop:" << in_que_cover_wait_for_pop << endl;
			if (max_parents[v] != numOfVertices) { 
				// If v is the first vertex having higher rank than source along its shortest path, add a shortcut.
				if (max_parents[v] == v && isNeighbor[v] == false) {
					shortcuts[source].push_back(CHGraph::CH_Edge(v, source, v_d));
					in_que_cover_wait_for_pop--;
				}
				else if (max_parents[v] == v && isNeighbor[v] == true) {
					in_que_cover_wait_for_pop--;
				}
				if (v < source)
					uncover_count--;
				if (uncover_count == 0) { // All vertices higher rank than source have been visited.
				//	cout << "good 1" << endl;
					break;
				}
				in_que_cover--; 
			}


		//	if (in_que == in_que_cover && v != source && in_que != 0)// All elments in queue have been covered.
		//		continue;
			 
			//if (source == 3)
			//	cout << "v: " << v << endl

			//if (in_que == in_que_cover && v != source && in_que_cover_wait_for_pop != -1 && in_que_cover_wait_for_pop != 0) {
			//		cout << "source: " << v << " inque:" << in_que << " inque_cover:" << in_que_cover << " inque_cover_wait_for_pop:" << in_que_cover_wait_for_pop << endl;
			//	//	cout << "good 2" << endl;  
			//	continue;
			//};

			for (int i = 0; i < adj[v].size(); ++i) {	
				NodeID w = adj[v][i].target;
				EdgeWeight w_d = adj[v][i].weight + v_d;				
				if (!vis[w]) {
					// When it is first inserted into the queue. 
					if (distances[w] == INF_WEIGHT && w_d < distances[w]) {
						in_que++;
						visited_set.insert(w);
					}
					if (distances[w] > w_d) {
						distances[w] = w_d;
						pqueue.update(w, w_d);

						//if (source == 3)
						//	cout << v << "," << max_parents[v] << "," << w << "," << max_parents[w] << endl;

						if (v == source)
							isNeighbor[w] = true;

						int uncover_flag = false;
						if (max_parents[w] == numOfVertices)
							uncover_flag = true;

					//	// w is the first vertex has higher rank than source.
					//	if (w < source && max_parents[v] == numOfVertices) {
					//		// Newly covered vertex in queue.
					//		/*if (max_parents[w] == numOfVertices) {
					//			in_que_cover++;  */
					//		if(waitForPop[w] == false){
					//			if (in_que_cover_wait_for_pop == -1)
					//				in_que_cover_wait_for_pop = 0;
					//			in_que_cover_wait_for_pop++;
					//			waitForPop[w] = true;
					//		}
					//		max_parents[w] = w;
					//	}
					//	
					//	// w-source has already been hit by previous max_parents[v].
					//	if (max_parents[v] != numOfVertices) {
					//		// Newly covered vertex in queue.
					///*		if(max_parents[w] == numOfVertices)
					//			in_que_cover++;*/
					//		if (waitForPop[w] == true) {
					//			if (in_que_cover_wait_for_pop == -1)
					//				in_que_cover_wait_for_pop = 0;
					//			else
					//				in_que_cover_wait_for_pop--;
					//			waitForPop[w] = false;
					//		}
					//		max_parents[w] = max_parents[v];
					//	}

					//	// w has already been covered previously but it is not the shortest path. it has not been covered yet.
					//	if (w > source && max_parents[v] == numOfVertices && max_parents[w] != numOfVertices) {
					//		max_parents[w] = numOfVertices;
					//		in_que_cover--; 
					//	}

					//	// w has not yet been covered. it is the first time it has been ocvered.
					//	if (uncover_flag == true && max_parents[w] != numOfVertices)
					//		in_que_cover++;
					//	



						// The shorter path comes from u-x-v-w rather than u(v)-w.
						if (isNeighbor[w] == true && v != source) {
							isNeighbor[w] = false;
						}

						// w has not been covered yet.
						if (max_parents[w] == numOfVertices) {
							// w is the first vertex has higher rank than source.
							if (w < source && max_parents[v] == numOfVertices) {
								// Newly covered vertex in queue.
								/*if (max_parents[w] == numOfVertices) {
									in_que_cover++;  */
								if(waitForPop[w] == false){
									if (in_que_cover_wait_for_pop == -1)
										in_que_cover_wait_for_pop = 0;
									in_que_cover_wait_for_pop++;
									waitForPop[w] = true;									
								}
								max_parents[w] = w;
							}

							// w-source has already been hit by previous max_parents[v].
							if (max_parents[v] != numOfVertices) {
								// Newly covered vertex in queue.
						/*		if(max_parents[w] == numOfVertices)
									in_que_cover++;*/
								if (waitForPop[w] == true) {
									if (in_que_cover_wait_for_pop == -1)
										in_que_cover_wait_for_pop = 0;
									else
										in_que_cover_wait_for_pop--;
									waitForPop[w] = false;
								}
								max_parents[w] = max_parents[v];
							}

							// w has not yet been covered. it is the first time it has been ocvered.
							if (uncover_flag == true && max_parents[w] != numOfVertices)
								in_que_cover++;
						}
						else { // max_parents[w] != numOfVertices // w has been covered previously

							// v has not been covered
							if (max_parents[v] == numOfVertices) {
								if (w < source) {//w is  the first higher rank vertex in this shortest path
									if (waitForPop[w] == false) {
										if (in_que_cover_wait_for_pop == -1)
											in_que_cover_wait_for_pop = 0;
										in_que_cover_wait_for_pop++;
										waitForPop[w] = true;
									}
									max_parents[w] = w;
								}
								else {//w has not been covered yet, remove its covered information
									max_parents[w] = numOfVertices;
									in_que_cover--;
								}
							}

							// w has already been covered by v. if w is the first higher rank vertex previously, remove this information
							if (max_parents[v] != numOfVertices) {
								if (max_parents[w] == w) {
									if (waitForPop[w] == true) {
										if (in_que_cover_wait_for_pop == -1)
											in_que_cover_wait_for_pop = 0;
										else
											in_que_cover_wait_for_pop--;
										waitForPop[w] = false;
									}
									max_parents[w] = max_parents[v];
								}
							}
						}
					}
				}
			}
		

			if (in_que == in_que_cover && v != source && in_que_cover_wait_for_pop != -1 && in_que_cover_wait_for_pop == 0) {
				//	cout << "source: " << v << " inque:" << in_que << " inque_cover:" << in_que_cover << " inque_cover_wait_for_pop:" << in_que_cover_wait_for_pop << endl;
				//	cout << "good 2" << endl;  
				break;
			}

		

		}	
		
	//	cout << "search space:" << sspace << endl;

		// Clean up the dijkstra structures otherwise the trash will manipulate the next inNode's dijkstra search.
		for (unordered_set<NodeID>::iterator it = visited_set.begin(); it != visited_set.end(); ++it) {
			NodeID cv = *it;
			vis[cv] = false;
			distances[cv] = INF_WEIGHT;
			max_parents[cv] = numOfVertices;
			isNeighbor[cv] = false;
			pqueue.clear(cv);
			waitForPop[cv] = false;
		}
	/*	while (!pqueue.empty()) {
			NodeID tmpv;
			EdgeWeight tmpweight;
			pqueue.extract_min(tmpv, tmpweight);
		}*/
		pqueue.clear_n();
		visited_set.clear();
		
		return 0;
	}

	//int build_shortcuts_dij(NodeID source, vector<vector<CHGraph::CH_Edge> >& adj, vector<bool>& vis, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, unordered_set<NodeID>& visited_set, vector<EdgeWeight>& distances, vector<NodeID>& max_parents, vector<bool>& isNeighbor, vector<vector<CHGraph::CH_Edge> >& shortcuts) {

	//	int uncover_count = source;
	//	pqueue.update(source, 0);
	//	distances[source] = 0;
	//	//isNeighbor[source] = true;
	//	visited_set.insert(source);

	//	while (!pqueue.empty()) {
	//		NodeID v;
	//		EdgeWeight v_d;
	//		pqueue.extract_min(v, v_d);
	//		vis[v] = true;

	//		if (max_parents[v] != numOfVertices) { // Either one of the ancestors of v or v has higher rank than source, we can stop the search from here. 
	//			if (max_parents[v] == v && isNeighbor[v] == false)
	//				shortcuts[source].push_back(CHGraph::CH_Edge(v, source, v_d));
	//			if (v < source)
	//				uncover_count--;
	//			if (uncover_count == 0) // All vertices higher rank than source have been covered.
	//				break;
	//			continue;
	//		}

	//		for (int i = 0; i < adj[v].size(); ++i) {
	//			NodeID w = adj[v][i].target;
	//			EdgeWeight w_d = adj[v][i].weight + v_d;
	//			if (!vis[w]) {
	//				if (distances[w] > w_d) {
	//					distances[w] = w_d;
	//					pqueue.update(w, w_d);
	//					visited_set.insert(w);

	//					if (v == source)
	//						isNeighbor[w] = true;


	//					// w is the first vertex has higher rank than source.
	//					if (w < source && max_parents[v] == numOfVertices)
	//						max_parents[w] = w;

	//					// w-source has already been hit by previous max_parents[v].
	//					if (max_parents[v] != numOfVertices)
	//						max_parents[w] = max_parents[v];

	//					// The shorter path comes from u-x-v-w rather than u(v)-w.
	//					if (isNeighbor[w] == true && v != source) {
	//						isNeighbor[w] = false;
	//					}


	//					/*			if (max_parents[w] > w)
	//					max_parents[w] = w;

	//					if(max_parents[w] > max_parents[v])
	//					max_parents[w] = max_parents[v];*/

	//				}
	//			}
	//		}
	//	}

	//	// Clean up the dijkstra structures otherwise the trash will manipulate the next inNode's dijkstra search.
	//	for (unordered_set<NodeID>::iterator it = visited_set.begin(); it != visited_set.end(); ++it) {
	//		NodeID cv = *it;
	//		vis[cv] = false;
	//		distances[cv] = INF_WEIGHT;
	//		max_parents[cv] = numOfVertices;
	//		isNeighbor[cv] = false;
	//		pqueue.clear(cv);
	//	}
	//	/*	while (!pqueue.empty()) {
	//	NodeID tmpv;
	//	EdgeWeight tmpweight;
	//	pqueue.extract_min(tmpv, tmpweight);
	//	}*/
	//	pqueue.clear_n();
	//	visited_set.clear();

	//	return 0;
	//}


	//int build_shortcuts_bfs(NodeID source, vector<vector<CHGraph::CH_Edge> >& adj, vector<bool>& vis, vector<NodeID>& max_parents, vector<bool>& isNeighbor, vector<vector<CHGraph::CH_Edge> >& shortcuts, vector<NodeID>& que) {

	//	int uncover_count = source;

	//	//vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);

	//		
	//		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
	//		que[que_h++] = source;
	//		vis[source] = true;
	//		que_t1 = que_h;

	//		for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
	//			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
	//				NodeID v = que[que_i];

	//				if (max_parents[v] != numOfVertices) { // Either one of the ancestors of v or v has higher rank than source, we can stop the search from here. 
	//					
	//					if (max_parents[v] == v && isNeighbor[v] == false) 
	//						shortcuts[source].push_back(CHGraph::CH_Edge(v, source, d));
	//						
	//					if (v < source)
	//						uncover_count--;
	//					if (uncover_count == 0) { // All vertices higher rank than source have been covered.
	//						break; 
	//					}
	//					continue;
	//				}

	//				for (size_t i = 0; i < adj[v].size(); ++i) {
	//					NodeID w = adj[v][i].target;
	//					if (!vis[w]) {
	//						vis[w] = true;

	//						if (v == source)
	//							isNeighbor[w] = true;

	//						// w is the first vertex has higher rank than source.
	//						if (w < source )//&& max_parents[v] == numOfVertices)
	//							max_parents[w] = w;

	//						//// w-source has already been hit by previous max_parents[v].
	//						//if (max_parents[v] != numOfVertices) {
	//						//	max_parents[w] = max_parents[v];
	//						//	continue;
	//						//}
	//												

	//						//// The shorter path comes from u-x-v-w rather than u(v)-w. *We wont have this case in unweigted graphs.
	//						//if (isNeighbor[w] == true && v != source) {
	//						//	isNeighbor[w] = false;
	//						//}


	//						que[que_h++] = w;

	//					}
	//				}
	//			pruned:
	//				{}
	//			}
	//			que_t0 = que_t1;
	//			que_t1 = que_h;
	//		}

	//		for (size_t i = 0; i < que_h; ++i) {
	//			vis[que[i]] = false;
	//			max_parents[que[i]] = numOfVertices;
	//			isNeighbor[que[i]] = false;
	//		}		

	//}
	
	int build_shortcuts_bfs(NodeID source, vector<vector<CHGraph::CH_Edge> >& adj, vector<bool>& vis, vector<NodeID>& max_parents, vector<bool>& isNeighbor, vector<vector<CHGraph::CH_Edge> >& shortcuts, vector<NodeID>& que) {

		int uncover_count = source;

		//vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);


		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
		que[que_h++] = source;
		vis[source] = true;
		que_t1 = que_h;

		int in_que = 1;
		int in_que_cover = 0;
		int in_que_wait_for_pop = -1;

		for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID v = que[que_i];
				in_que--;

				if (max_parents[v] != numOfVertices) { // Either one of the ancestors of v or v has higher rank than source, we can stop the search from here. 

					if (max_parents[v] == v && isNeighbor[v] == false) {
						shortcuts[source].push_back(CHGraph::CH_Edge(v, source, d));
						in_que_wait_for_pop--;
					}

					if(max_parents[v] == v && isNeighbor[v] ==  true)
						in_que_wait_for_pop--;

					if (v < source)
						uncover_count--;
					if (uncover_count == 0) { // All vertices higher rank than source have been covered.
						goto jumpout;
					} 

					in_que_cover--;
				}


				for (size_t i = 0; i < adj[v].size(); ++i) {
					NodeID w = adj[v][i].target;
					if (!vis[w]) {
						vis[w] = true;
						in_que++;

						if (v == source)
							isNeighbor[w] = true;

						// w is the first vertex has higher rank than source.
						if (w < source && max_parents[v] == numOfVertices) {//&& max_parents[v] == numOfVertices)
							max_parents[w] = w;
							if (in_que_wait_for_pop == -1)
								in_que_wait_for_pop = 0;
							in_que_wait_for_pop++;
							in_que_cover++;
						}

						if (max_parents[v] != numOfVertices) {
							in_que_cover++;
							max_parents[w] = max_parents[v];
						}

						//// w-source has already been hit by previous max_parents[v].
						//if (max_parents[v] != numOfVertices) {
						//	max_parents[w] = max_parents[v];
						//	continue;
						//}


						//// The shorter path comes from u-x-v-w rather than u(v)-w. *We wont have this case in unweigted graphs.
						//if (isNeighbor[w] == true && v != source) {
						//	isNeighbor[w] = false;
						//}


						que[que_h++] = w;

					}
				}


				if (in_que == in_que_cover && in_que_wait_for_pop == 0) {
					goto jumpout;
				}
			pruned:
				{}
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}

	jumpout:
		{}

		for (size_t i = 0; i < que_h; ++i) {
			vis[que[i]] = false;
			max_parents[que[i]] = numOfVertices;
			isNeighbor[que[i]] = false;
		}

	}


	void labeling(CHGraph &chgraph, Ordering &orders) {

		for (int i = 0; i < numOfVertices; ++i) {
			if(chgraph.adj[i].size() != 0)
				sort(chgraph.adj[i].begin(), chgraph.adj[i].end());
			if (DIRECTED_FLAG == true) {
				if (chgraph.r_adj[i].size() != 0)
					sort(chgraph.r_adj[i].begin(), chgraph.r_adj[i].end());
			}
		}

		vector<index_t>& index_ = labels.index_;
		vector<NodeID> &inv = orders.inv;
		vector<NodeID> &rank = orders.rank;

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(),
				vector<EdgeWeight>()));

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(),
				vector<EdgeWeight>()));

		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<EdgeWeight> r_dst_r(numOfVertices + 1, INF_WEIGHT);

		set<NodeID> appendCandidate;

	//	vector<CHGraph::CH_Edge>& inEdges = chgraph.r_adj[v];
		//vector<CHGraph::CH_Edge>& outEdges = chgraph.adj[v];


		// Start to process every nodes.
		for (NodeID v = 0; v < numOfVertices; ++v) {
			vector<CHGraph::CH_Edge>&  adj_v = chgraph.adj[v];
			pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
			
			if (DIRECTED_FLAG == false) {
				//for (int i = 0; i < tmp_idx_v.first.size() - 1; ++i)
				//	dst_r[tmp_idx_v.first[i]] = tmp_idx_v.second[i];

				// Search all the neighrbor u of v.
				for (int i = 0; i < adj_v.size(); ++i) {
					NodeID u = adj_v[i].target;
					NodeID level = adj_v[i].level;
					if (u >= v ) continue;
					EdgeWeight weight = adj_v[i].weight;

					// Check all the existing labels w in the labels of u.
					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_u = tmp_idx[u];
					for (int j = 0; j < tmp_idx_u.first.size(); ++j) {
						NodeID w = tmp_idx_u.first[j];
						EdgeWeight vuw_distance = tmp_idx_u.second[j] + weight; // Distance from v to u to w.
						// Check whether the path (v->u->w) should be added into the labels of v. We prune it if the existing labels of v and w can answer the distance directly.
						//pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_w = tmp_idx[w];
						//for (int k = 0; k < tmp_idx_w.first.size(); ++k) {
						//	NodeID labelOfW = tmp_idx_w.first[k];
						//	EdgeWeight distanceToIt = tmp_idx_w.second[k] + dst_r[labelOfW];

						//	// Prune this label, we do not add (w, vuw_distance) to the label set of v.
						//	if (distanceToIt <= vuw_distance) goto pruning;
						//}

						if (dst_r[w] > vuw_distance) {
							if (dst_r[w] == INF_WEIGHT) appendCandidate.insert(w);				
							dst_r[w] = vuw_distance;
						}
						//pruning: {}
					}
				}

				// Post prune the candidate by running query test, start from the max rank candidate. It takes O(|M|^2) time, |M| is the average label size.
				for (set<NodeID>::iterator it = appendCandidate.begin(); it != appendCandidate.end(); ++it) {
					NodeID w = *it;
					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_w = tmp_idx[w];
					for (int k = 0; k < tmp_idx_w.first.size(); ++k) {
						NodeID labelOfW = tmp_idx_w.first[k];
						EdgeWeight distanceToIt = tmp_idx_w.second[k] + dst_r[labelOfW];
						if (dst_r[w] >= distanceToIt) {
							if (labelOfW != w) 
								dst_r[w] = INF_WEIGHT; // We prune it because the this path has been hit by another higher rank vertex labelOfW.
							break;
						}
					}
				}

				for (set<NodeID>::iterator it = appendCandidate.begin(); it != appendCandidate.end(); ++it) {
					NodeID w = *it;
					if (dst_r[w] == INF_WEIGHT) continue;					
					tmp_idx_v.first.push_back(w);
					tmp_idx_v.second.push_back(dst_r[w]);
					dst_r[w] = INF_WEIGHT;
				}
				appendCandidate.clear();
				tmp_idx_v.first.push_back(v);
				tmp_idx_v.second.push_back(0);

			}
		}

		if (DIRECTED_FLAG == false) {
			for (size_t v = 0; v < numOfVertices; ++v) {
				tmp_idx[v].first.push_back(numOfVertices);
				tmp_idx[v].second.push_back(INF_WEIGHT);
				NodeID k = tmp_idx[v].first.size();
				index_[inv[v]].spt_v.resize(k);
				index_[inv[v]].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_d[i] = tmp_idx[v].second[i];
//				tmp_idx[v].first.clear();
				vector<NodeID>().swap(tmp_idx[v].first);
				vector<EdgeWeight>().swap(tmp_idx[v].second);
//				tmp_idx[v].second.clear(); 
			}
		}
	}
	
	void td_labeling(CHGraph &chgraph, Ordering &orders) {

		for (int i = 0; i < numOfVertices; ++i) {
			if (chgraph.adj[i].size() != 0)
				sort(chgraph.adj[i].rbegin(), chgraph.adj[i].rend());
			if (chgraph.r_adj[i].size() != 0)
				sort(chgraph.r_adj[i].rbegin(), chgraph.r_adj[i].rend());
		}

		vector<index_t>& index_ = labels.index_;
		vector<NodeID> &inv = orders.inv;
		vector<NodeID> &rank = orders.rank;

	/*	vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(),
				vector<EdgeWeight>()));*/

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(),
				vector<EdgeWeight>()));

		//vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<EdgeWeight> r_dst_r(numOfVertices + 1, INF_WEIGHT);

		set<NodeID> appendCandidate;

		//	vector<CHGraph::CH_Edge>& inEdges = chgraph.r_adj[v];
		//vector<CHGraph::CH_Edge>& outEdges = chgraph.adj[v];

		//vector<index_t>& index_ = labels.index_;
		//vector<NodeID> &inv = orders.inv;
		//vector<NodeID> &rank = orders.rank;
		//	vector<vector<NodeID> > &adj = wgraph.adj;
		//	vector<vector<EdgeWeight> > &adj_weight = wgraph.adj_weight;

		//vector<EdgeID>& vertices = wgraph.vertices;
		//vector<NodeEdgeWeightPair>& edges = wgraph.edges;

		vector<bool> usd(numOfVertices, false);

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<bool> vis(numOfVertices);
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		queue<NodeID> visited_que;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);


		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);
		long pop = 0;
		double hsize = 0;
		for (size_t r = 0; r < numOfVertices; ++r) {
			if (usd[r]) continue;

			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[r];

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}

			pqueue.update(r, 0);
			//vis[r] = true;
			long max_heap_size = 0;
			long heap_size = 1;
			while (!pqueue.empty()) {

				pop++;
				heap_size--;
				NodeID v;
				EdgeWeight v_d;
				pqueue.extract_min(v, v_d);
				pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
				vis[v] = true;
				visited_que.push(v);
				vector<CHGraph::CH_Edge>&  adj_v = chgraph.adj[v];

				if (usd[v]) continue;
				for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
					NodeID w = tmp_idx_v.first[i];
					EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
					if (td <= v_d) {
					//	pruning_power[w]++;
						goto pruned;
					}
				}

				// Traverse
				tmp_idx_v.first.back() = r;
				tmp_idx_v.second.back() = v_d;
				tmp_idx_v.first.push_back(numOfVertices);
				tmp_idx_v.second.push_back(INF_WEIGHT);

				iteration_generated[r]++;

				//	for (size_t i = 0; i < adj[v].size(); ++i) {
				//		NodeID w = adj[v][i];
				//		EdgeWeight w_d = adj_weight[v][i] + v_d;


				//for (EdgeID eid = vertices[v]; eid < vertices[v + 1]; ++eid) {
				//	NodeID w = edges[eid].first;
				//	EdgeWeight w_d = edges[eid].second + v_d;

			for (int i = 0; i < adj_v.size(); ++i) {
				NodeID w = adj_v[i].target;
				EdgeWeight w_d = adj_v[i].weight + v_d;
				if (w < v) break; //Only the neighbors with lower ranks will be visited.
				if (!vis[w]) {
		
					if (distances[w] == INF_WEIGHT) {
						heap_size++;
						if (max_heap_size < heap_size)
							max_heap_size = heap_size;
					}
					if (distances[w] > w_d) {
						pqueue.update(w, w_d);
						distances[w] = w_d;
					}
					
				}
			}
			pruned:
				{}
			}
			hsize = hsize + max_heap_size;
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
			usd[r] = true;
		}

		cout << "total pop:" << pop << endl;
		cout << "average max heap size " << (double)hsize / (double)numOfVertices << endl;

		for (size_t v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			index_[inv[v]].spt_v.resize(k);
			index_[inv[v]].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_d[i] = tmp_idx[v].second[i];
			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();
		}

	}


	void labeling_directed(CHGraph &chgraph, Ordering &orders) {
		for (int i = 0; i < numOfVertices; ++i) {
			/* cout << i << endl;
			 if(i == 5073)
				 cout << chgraph.adj[i].size() << endl;
			for (int j = 0; j < chgraph.adj[i].size(); ++j) {
				if (i == 5073)
					cout << "," << j;
				NodeID hh = chgraph.adj[i][j].target;
			}
			 cout << "test ok " << chgraph.adj[i].size();
			for (int j = 0; j < chgraph.r_adj[i].size(); ++j) {
				NodeID hh = chgraph.r_adj[i][j].target;
			}*/
			sort(chgraph.adj[i].begin(), chgraph.adj[i].end());
			sort(chgraph.r_adj[i].begin(), chgraph.r_adj[i].end());
			//cout << "," << chgraph.r_adj[i].size() << endl;
		}

		vector<index_t>& index_ = dlabels.index_;
		vector<index_t>& bindex_ = dlabels.bindex_;
		vector<NodeID> &inv = orders.inv;
		vector<NodeID> &rank = orders.rank;

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(),
				vector<EdgeWeight>()));

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(),
				vector<EdgeWeight>()));

		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<EdgeWeight> r_dst_r(numOfVertices + 1, INF_WEIGHT);

		set<NodeID> appendCandidate;
		
		// The main process:
		// Forward search: fetch necessary outLabel of the forward adj. (all outgoing labels must pass out neighbors)
		// Backward search: fetch necessary inLabel of the backward adj.(all incoming labels must pass in neighbors)
		// Start to process every nodes.
		for (NodeID v = 0; v < numOfVertices; ++v) {

			// Forward search.
			{
				vector<CHGraph::CH_Edge>&  adj_v = chgraph.adj[v];
				pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
				
				// Search all the neighrbor u of v.
				for (int i = 0; i < adj_v.size(); ++i) {
					NodeID u = adj_v[i].target;
					NodeID level = adj_v[i].level;
					if (u >= v) continue;
					EdgeWeight weight = adj_v[i].weight;
					// Check all the existing labels w in the labels of u.
					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_u = tmp_idx[u];
					for (int j = 0; j < tmp_idx_u.first.size(); ++j) {
						NodeID w = tmp_idx_u.first[j];
						EdgeWeight vuw_distance = tmp_idx_u.second[j] + weight; // Distance from v to u to w.
						if (dst_r[w] > vuw_distance) {
							if (dst_r[w] == INF_WEIGHT) appendCandidate.insert(w);
							dst_r[w] = vuw_distance;
						}
					}
				}

				// Post prune the candidate by running query test, start from the max rank candidate. It takes O(|M|^2) time, |M| is the average label size.
				for (set<NodeID>::iterator it = appendCandidate.begin(); it != appendCandidate.end(); ++it) {
					NodeID w = *it;
					pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_w = r_tmp_idx[w];
					for (int k = 0; k < r_tmp_idx_w.first.size(); ++k) {
						NodeID labelOfW = r_tmp_idx_w.first[k];
						EdgeWeight distanceToIt = r_tmp_idx_w.second[k] + dst_r[labelOfW];
						if (dst_r[w] >= distanceToIt) {
							if (labelOfW != w)
								dst_r[w] = INF_WEIGHT; // We prune it because the this path has been hit by another higher rank vertex labelOfW.
							break;
						}
					}
				}

				for (set<NodeID>::iterator it = appendCandidate.begin(); it != appendCandidate.end(); ++it) {
					NodeID w = *it;
					if (dst_r[w] == INF_WEIGHT) continue;
					tmp_idx_v.first.push_back(w);
					tmp_idx_v.second.push_back(dst_r[w]);
					dst_r[w] = INF_WEIGHT;
				}
				appendCandidate.clear();
				tmp_idx_v.first.push_back(v);
				tmp_idx_v.second.push_back(0);
			}

			// Backward search.
			{
				vector<CHGraph::CH_Edge>&  r_adj_v = chgraph.r_adj[v];
				pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];
				
				// Search all the neighrbor u of v.
				for (int i = 0; i < r_adj_v.size(); ++i) {
					NodeID u = r_adj_v[i].target;
					NodeID level = r_adj_v[i].level;
					if (u >= v) continue;
					EdgeWeight weight = r_adj_v[i].weight;
					// Check all the existing labels w in the labels of u.
					pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_u = r_tmp_idx[u];
					for (int j = 0; j < r_tmp_idx_u.first.size(); ++j) {
						NodeID w = r_tmp_idx_u.first[j];
						EdgeWeight wuv_distance = r_tmp_idx_u.second[j] + weight; // Distance from w to u to v.																
						if (r_dst_r[w] > wuv_distance) {
							if (r_dst_r[w] == INF_WEIGHT) appendCandidate.insert(w);
							r_dst_r[w] = wuv_distance;
						}
						//pruning: {}
					}
				}

				// Post prune the candidate by running query test, start from the max rank candidate. It takes O(|M|^2) time, |M| is the average label size.
				for (set<NodeID>::iterator it = appendCandidate.begin(); it != appendCandidate.end(); ++it) {
					NodeID w = *it;
					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_w = tmp_idx[w];
					for (int k = 0; k < tmp_idx_w.first.size(); ++k) {
						NodeID labelOfW = tmp_idx_w.first[k];
						EdgeWeight distanceToIt = tmp_idx_w.second[k] + r_dst_r[labelOfW];
						if (r_dst_r[w] >= distanceToIt) {
							if (labelOfW != w)
								r_dst_r[w] = INF_WEIGHT; // We prune it because the this path has been hit by another higher rank vertex labelOfW.
							break;
						}
					} 
				}
				for (set<NodeID>::iterator it = appendCandidate.begin(); it != appendCandidate.end(); ++it) {
					NodeID w = *it;
					if (r_dst_r[w] == INF_WEIGHT) continue;
					r_tmp_idx_v.first.push_back(w);
					r_tmp_idx_v.second.push_back(r_dst_r[w]);
					r_dst_r[w] = INF_WEIGHT;
				}

				appendCandidate.clear();
				r_tmp_idx_v.first.push_back(v);
				r_tmp_idx_v.second.push_back(0);
			}
		}

		int added = 0;
		for (size_t v = 0; v < numOfVertices; ++v ){
			tmp_idx[v].first.push_back(numOfVertices);
			tmp_idx[v].second.push_back(INF_WEIGHT);
			NodeID k = tmp_idx[v].first.size();
			index_[inv[v]].spt_v.resize(k);
			index_[inv[v]].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_d[i] = tmp_idx[v].second[i];
			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();
			added += k;

			r_tmp_idx[v].first.push_back(numOfVertices);
			r_tmp_idx[v].second.push_back(INF_WEIGHT);
			k = r_tmp_idx[v].first.size();
			bindex_[inv[v]].spt_v.resize(k);
			bindex_[inv[v]].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) bindex_[inv[v]].spt_v[i] = r_tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) bindex_[inv[v]].spt_d[i] = r_tmp_idx[v].second[i];
			r_tmp_idx[v].first.clear();
			r_tmp_idx[v].second.clear();
			added += k;
		}
		cout << added << endl;


	}
	
	Bottomup(CHGraph &chgraph, Ordering &orders, double& time_contracting, const double SWITCH_DEGREE_PARA, const int HOP_LIMIT_PARA) {
		iteration_generated.resize(numOfVertices);

		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);
		unordered_set<NodeID> visited_set;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
		vector<bool> hitTarget(numOfVertices);// false
		vector<bool> vis(numOfVertices);

		
		contracted.resize(numOfVertices, false);
		possibleWitness.resize(numOfVertices, INF_WEIGHT);
		mtmBucket.resize(numOfVertices);
				
		// Temporary structures for Dijkstra searches.
		NodeID remainNode = numOfVertices;
		int currentEdges = relabelByOrder(chgraph, orders);

		int hopLimitsParameter = 1;

		double currentDegree = (double)currentEdges / (double)numOfVertices;
		time_contracting = GetCurrentTimeSec();
		// A main loop to pick all vertices from the least to the most important ones. 
		for (NodeID v = numOfVertices - 1; v > 0; --v) {
			if(v % (numOfVertices /10) == 0)
				cout << "Building shortcuts for " << v << "th vertices("<< orders.inv[v] <<"). Total shortcuts ratio so far " << currentEdges / numOfEdges << " time:" << GetCurrentTimeSec() - time_contracting << "s" << endl;
			contracted[v] = true;

			double 	time_thisround = GetCurrentTimeSec();
			int addShortcuts = witness_search(v, chgraph, hopLimitsParameter, orders, vis, pqueue, visited_set, distances, hitTarget);
			iteration_generated[v] = GetCurrentTimeSec() - time_thisround;

			currentEdges += addShortcuts;
			currentDegree = (double)currentEdges / (double)v;
			if (HOP_LIMIT_PARA == 2 || HOP_LIMIT_PARA == 1 || HOP_LIMIT_PARA == 5) {
				hopLimitsParameter = HOP_LIMIT_PARA;
				continue;
			}
			if (currentDegree > 3.3 && currentDegree <= SWITCH_DEGREE_PARA && hopLimitsParameter != 2) {
				hopLimitsParameter = 2;
				cout << "At Vertex " << v << ". Switch to 2-hop limit with average degree:" << currentDegree << " #edge:" << currentEdges << endl;
			}else if (currentDegree <= 3.3 && hopLimitsParameter !=1) {
				hopLimitsParameter = 1;
				cout << "At Vertex " << v << ". Switch to 1-hop limit with average degree:" << currentDegree << " #edge:" << currentEdges << endl;
			}else if(hopLimitsParameter != 5 && currentDegree > SWITCH_DEGREE_PARA){
				hopLimitsParameter = 5;
				cout << "At Vertex " << v << ". Switch to multi-hop limit with average degree:" << currentDegree << " #edge:" << currentEdges << endl;
			}
		}
		time_contracting = GetCurrentTimeSec() - time_contracting;
		cout << "Spent " << time_contracting << " s on creating shortcuts..." << endl;
		cout << "Start to label vertices..." << endl;
		labeling(chgraph, orders);
	}

	Bottomup(CHGraph &chgraph, Ordering &orders, bool directed_flag, double& time_contracting, const double SWITCH_DEGREE_PARA, const int HOP_LIMIT_PARA) {

		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);
		unordered_set<NodeID> visited_set;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
		vector<bool> hitTarget(numOfVertices);// false
		vector<bool> vis(numOfVertices);
		iteration_generated.resize(numOfVertices);

		contracted.resize(numOfVertices, false);
		possibleWitness.resize(numOfVertices, INF_WEIGHT);
		mtmBucket.resize(numOfVertices);

		// Temporary structures for Dijkstra searches.
		NodeID remainNode = numOfVertices;
		int currentEdges = relabelByOrder(chgraph, orders);

		int hopLimitsParameter = 1;

		double currentDegree = (double)currentEdges / (double)numOfVertices;
		time_contracting = GetCurrentTimeSec();
		// A main loop to pick all vertices from the least to the most important ones. 
		for (NodeID v = numOfVertices - 1; v > 0; --v) {
			
			if (v % (numOfVertices / 10) == 0)
				cout << "Building shortcuts for " << v << "th vertices(" << orders.inv[v] << "). Total shortcuts ratio so far " << currentEdges / numOfEdges << " time:" << GetCurrentTimeSec()- time_contracting << "s" <<  endl;
			contracted[v] = true;			

			double 	time_thisround = GetCurrentTimeSec();
			int addShortcuts = witness_search_directed(v, chgraph, hopLimitsParameter, orders, vis, pqueue, visited_set, distances, hitTarget);
			iteration_generated[v] = GetCurrentTimeSec() - time_thisround;

			currentEdges += addShortcuts;
			currentDegree = (double)currentEdges / (double)v;			
			
			if (HOP_LIMIT_PARA == 2 || HOP_LIMIT_PARA == 1 || HOP_LIMIT_PARA == 5) {
				hopLimitsParameter = HOP_LIMIT_PARA;
				continue;
			}

			if (currentDegree > 3.3 && currentDegree <= SWITCH_DEGREE_PARA && hopLimitsParameter != 2) {
				hopLimitsParameter = 2;
				cout << "At Vertex " << v << ". Switch to 2-hop limit with average degree:" << currentDegree << " #edge:" << currentEdges << endl;
			}
			else if (currentDegree <= 3.3 && hopLimitsParameter != 1) {
				hopLimitsParameter = 1;
				cout << "At Vertex " << v << ". Switch to 1-hop limit with average degree:" << currentDegree << " #edge:" << currentEdges << endl;
			}
			else if (hopLimitsParameter != 5 && currentDegree > SWITCH_DEGREE_PARA) {
				hopLimitsParameter = 5;
				cout << "At Vertex " << v << ". Switch to multi-hop limit with average degree:" << currentDegree << " #edge:" << currentEdges << endl;
			}
		}
		time_contracting = GetCurrentTimeSec() - time_contracting;
		cout << "Spent " << time_contracting << " s on creating shortcuts..." << endl;
		cout << "Start to label vertices..." << endl;
		labeling_directed(chgraph, orders);
	}

	Bottomup(CHGraph &chgraph, Ordering &orders, double& time_contracting) {
		vector<vector<CHGraph::CH_Edge> > shortcuts(numOfVertices);
		vector<vector<CHGraph::CH_Edge> > r_shortcuts(numOfVertices);
		//benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);
		//unordered_set<NodeID> visited_set;
		//vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
		//vector<bool> hitTarget(numOfVertices);// false
		//vector<bool> vis(numOfVertices);
		//iteration_generated.resize(numOfVertices);
		//vector<NodeID> max_parents(numOfVertices, numOfVertices);
		//vector<bool> isNeighbor(numOfVertices, false);
		//vector<bool> waitForPop(numOfVertices, false);
		//
		//vector<NodeID> que(numOfVertices);

		relabelByOrder(chgraph, orders);
		long total_shortcut = 0;
		time_contracting = GetCurrentTimeSec();
		double time_shortcuting = GetCurrentTimeSec();

		int num_threads = 10;
		omp_set_num_threads(num_threads);
		vector<vector<vector<CHGraph::CH_Edge> > > shortcuts_omp(num_threads, vector<vector<CHGraph::CH_Edge> >(numOfVertices));
		vector<vector<vector<CHGraph::CH_Edge> > > r_shortcuts_omp(num_threads, vector<vector<CHGraph::CH_Edge> >(numOfVertices));


		vector<benchmark::heap<2, EdgeWeight, NodeID> > pqueue(num_threads, benchmark::heap<2, EdgeWeight, NodeID>(numOfVertices) );
		vector<unordered_set<NodeID> > visited_set(num_threads);
		vector<vector<EdgeWeight> > distances(num_threads, vector<EdgeWeight>(numOfVertices, INF_WEIGHT));
		vector<vector<bool> > hitTarget(num_threads, vector<bool>(numOfVertices) );// false
		vector<vector<bool> > vis(num_threads, vector<bool>(numOfVertices));
		iteration_generated.resize(numOfVertices);
		vector<vector<NodeID> > max_parents(num_threads, vector<NodeID>(numOfVertices, numOfVertices));
		vector<vector<bool> > isNeighbor(num_threads, vector<bool>(numOfVertices, false));
		vector<vector<bool> > waitForPop(num_threads, vector<bool>(numOfVertices, false));
		vector<vector<NodeID> > que(num_threads, vector<NodeID>(numOfVertices));


		#pragma omp parallel for schedule(dynamic)
		for (NodeID v = numOfVertices - 1; v > 0; --v) {
			if (WEIGHTED_FLAG == true) {
				//build_shortcuts_dij(v, chgraph.adj, vis, pqueue, visited_set, distances, max_parents, isNeighbor, shortcuts);
				build_shortcuts_dij(v, chgraph.adj, vis[omp_get_thread_num()], pqueue[omp_get_thread_num()], visited_set[omp_get_thread_num()], distances[omp_get_thread_num()], max_parents[omp_get_thread_num()], isNeighbor[omp_get_thread_num()], waitForPop[omp_get_thread_num()], shortcuts_omp[omp_get_thread_num()]);
				if (DIRECTED_FLAG == true)
					//	build_shortcuts_dij(v, chgraph.r_adj, vis, pqueue, visited_set, distances, max_parents, isNeighbor, r_shortcuts);
					build_shortcuts_dij(v, chgraph.r_adj, vis[omp_get_thread_num()], pqueue[omp_get_thread_num()], visited_set[omp_get_thread_num()], distances[omp_get_thread_num()], max_parents[omp_get_thread_num()], isNeighbor[omp_get_thread_num()], waitForPop[omp_get_thread_num()], r_shortcuts_omp[omp_get_thread_num()]);
			}
			else {
				build_shortcuts_bfs(v, chgraph.adj, vis[omp_get_thread_num()], max_parents[omp_get_thread_num()], isNeighbor[omp_get_thread_num()], shortcuts_omp[omp_get_thread_num()], que[omp_get_thread_num()]);
		//		build_shortcuts_bfs(v, chgraph.adj, vis[omp_get_thread_num()], max_parents[omp_get_thread_num()], isNeighbor[omp_get_thread_num()], shortcuts_omp[omp_get_thread_num()], que[omp_get_thread_num()]);
				if (DIRECTED_FLAG == true)
					build_shortcuts_bfs(v, chgraph.r_adj, vis[omp_get_thread_num()], max_parents[omp_get_thread_num()], isNeighbor[omp_get_thread_num()], r_shortcuts_omp[omp_get_thread_num()], que[omp_get_thread_num()]);
					//build_shortcuts_bfs(v, chgraph.r_adj, vis[omp_get_thread_num()], max_parents[omp_get_thread_num()], isNeighbor[omp_get_thread_num()], r_shortcuts_omp[omp_get_thread_num()], que[omp_get_thread_num()]);
			}
		}

		for (int t = 0; t < num_threads; ++t) {
			for (int v = 0; v < numOfVertices; ++v) {
				for (int i = 0; i < shortcuts_omp[t][v].size(); ++i) {
					shortcuts[v].push_back(shortcuts_omp[t][v][i]);
				}
				for (int i = 0; i < r_shortcuts_omp[t][v].size(); ++i) {
					r_shortcuts[v].push_back(r_shortcuts_omp[t][v][i]);
				}
			}
		}



		//for (NodeID v = numOfVertices - 1; v > 0; --v) {
		//	if (WEIGHTED_FLAG == true) {
		//		//build_shortcuts_dij(v, chgraph.adj, vis, pqueue, visited_set, distances, max_parents, isNeighbor, shortcuts);
		//		build_shortcuts_dij(v, chgraph.adj, vis, pqueue, visited_set, distances, max_parents, isNeighbor, waitForPop, shortcuts);
		//		if (DIRECTED_FLAG == true)
		//		//	build_shortcuts_dij(v, chgraph.r_adj, vis, pqueue, visited_set, distances, max_parents, isNeighbor, r_shortcuts);
		//			build_shortcuts_dij(v, chgraph.r_adj, vis, pqueue, visited_set, distances, max_parents, isNeighbor, waitForPop, r_shortcuts);
		//	}	else {
		//		build_shortcuts_bfs(v, chgraph.adj, vis,  max_parents, isNeighbor, shortcuts, que);
		//		if (DIRECTED_FLAG == true)
		//			build_shortcuts_bfs(v, chgraph.r_adj, vis, max_parents, isNeighbor, r_shortcuts, que);
		//	}
		//}

		double time_searching = GetCurrentTimeSec() - time_contracting;

		cout << "Spent " << time_searching << " s on searching shortcuts..." << endl;


		if (DIRECTED_FLAG == false)
			total_shortcut += process_shortcuts(chgraph, shortcuts);
		else
			total_shortcut += process_shortcuts_directed(chgraph, shortcuts, r_shortcuts);

		time_contracting = GetCurrentTimeSec() - time_contracting - time_searching;
		cout << "Spent " << time_contracting << " s on processing shortcuts..." << endl;
		cout << "Start to label vertices..." << endl;
		cout << total_shortcut << " shortcuts in total." << endl;

		time_shortcuting = GetCurrentTimeSec() - time_shortcuting;
		time_contracting = time_shortcuting;
		if (DIRECTED_FLAG == false)
			td_labeling(chgraph, orders);
			//labeling(chgraph, orders);
		else
			labeling_directed(chgraph, orders);

	}
	
	Bottomup(CHGraph &chgraph, Ordering &orders, double& time_contracting, bool CH_ORDER_FLAGS) {


		relabelByOrder(chgraph, orders);
		labeling(chgraph, orders);
	}

	long process_shortcuts(CHGraph& chgraph, vector<vector<CHGraph::CH_Edge> > shortcuts) {
		vector<vector<bool> > add_out(numOfVertices);
		vector<vector<bool> > add_in(numOfVertices);
		long added_shortcuts = 0;
		// Test whether shortcut can replace the original edges.
		for (NodeID v = 0; v < numOfVertices; ++v) {
			vector<CHGraph::CH_Edge>& shortcuts_v = shortcuts[v];
			vector<bool>& add_out_v = add_out[v];
			vector<bool>& add_in_v = add_in[v];
			add_out_v.resize(shortcuts_v.size(), true);
			add_in_v.resize(shortcuts_v.size(), true);
			for (int i = 0; i < shortcuts_v.size(); ++i) {
				NodeID w = shortcuts_v[i].target;
				vector<CHGraph::CH_Edge>::iterator pos = lower_bound(chgraph.adj[v].begin(), chgraph.adj[v].end(), shortcuts_v[i]);
				if (pos != chgraph.adj[v].end() && (*pos).target == shortcuts_v[i].target) {
					if ((*pos).weight > shortcuts_v[i].weight)
						(*pos).weight = shortcuts_v[i].weight;
					add_out_v[i] = false;
				}

				CHGraph::CH_Edge dummy(v, 0, 0);
				pos = lower_bound(chgraph.adj[w].begin(), chgraph.adj[w].end(), dummy);
				if (pos != chgraph.adj[w].end() && (*pos).target == v) {
					if ((*pos).weight > shortcuts_v[i].weight)
						(*pos).weight = shortcuts_v[i].weight;
					add_in_v[i] = false;
				}
			}
		}

		// Append the actual shortcuts.
		for (NodeID v = 0; v < numOfVertices; ++v) {
			vector<CHGraph::CH_Edge>& shortcuts_v = shortcuts[v];
			vector<bool>& add_out_v = add_out[v];
			vector<bool>& add_in_v = add_in[v];
			for (int i = 0; i < shortcuts_v.size(); ++i) {
				NodeID w = shortcuts_v[i].target;
				NodeID level = shortcuts_v[i].level;
				EdgeWeight weight = shortcuts_v[i].weight;
				if (add_out_v[i] == true)
					chgraph.adj[v].push_back(CHGraph::CH_Edge(w, level, weight));
				if (add_in_v[i] == true)
					chgraph.adj[w].push_back(CHGraph::CH_Edge(v, level, weight));
				if (add_out_v[i]==false && add_in_v[i]==false)
					added_shortcuts--;
				added_shortcuts++;
			}
		}

		return added_shortcuts;
	}

	long process_shortcuts_directed(CHGraph& chgraph, vector<vector<CHGraph::CH_Edge> > shortcuts, vector<vector<CHGraph::CH_Edge> > r_shortcuts) {
		vector<vector<bool> > add_out(numOfVertices);
		vector<vector<bool> > add_in(numOfVertices);
		vector<vector<bool> > r_add_out(numOfVertices);
		vector<vector<bool> > r_add_in(numOfVertices);

		long added_shortcuts = 0;
		// Test whether shortcut can replace the original edges.
		for (NodeID v = 0; v < numOfVertices; ++v) {
			{// Forward shortcut from (v->w), so add this shortcut into adj[v] and r_adj[w].
				vector<CHGraph::CH_Edge>& shortcuts_v = shortcuts[v];
				vector<bool>& add_out_v = add_out[v];
				vector<bool>& add_in_v = add_in[v];
				add_out_v.resize(shortcuts_v.size(), true);
				add_in_v.resize(shortcuts_v.size(), true);
				for (int i = 0; i < shortcuts_v.size(); ++i) {
					NodeID w = shortcuts_v[i].target;
					vector<CHGraph::CH_Edge>::iterator pos = lower_bound(chgraph.adj[v].begin(), chgraph.adj[v].end(), shortcuts_v[i]);
					if (pos != chgraph.adj[v].end() && (*pos).target == shortcuts_v[i].target) {
						if ((*pos).weight > shortcuts_v[i].weight)
							(*pos).weight = shortcuts_v[i].weight;
						add_out_v[i] = false;
					}

					CHGraph::CH_Edge dummy(v, 0, 0);
					pos = lower_bound(chgraph.r_adj[w].begin(), chgraph.r_adj[w].end(), dummy);
					if (pos != chgraph.r_adj[w].end() && (*pos).target == v) {
						if ((*pos).weight > shortcuts_v[i].weight)
							(*pos).weight = shortcuts_v[i].weight;
						add_in_v[i] = false;
					}
				}
			}


			{// Backward shortcuts, add (w->v) into adj[w] and r_adj[v].
				vector<CHGraph::CH_Edge>& r_shortcuts_v = r_shortcuts[v];
				vector<bool>& r_add_out_v = r_add_out[v];
				vector<bool>& r_add_in_v = r_add_in[v];
				r_add_out_v.resize(r_shortcuts_v.size(), true);
				r_add_in_v.resize(r_shortcuts_v.size(), true);
				for (int i = 0; i < r_shortcuts_v.size(); ++i) {
					NodeID w = r_shortcuts_v[i].target;

					vector<CHGraph::CH_Edge>::iterator pos = lower_bound(chgraph.r_adj[v].begin(), chgraph.r_adj[v].end(), r_shortcuts_v[i]);
					if (pos != chgraph.r_adj[v].end() && (*pos).target == r_shortcuts_v[i].target) {
						if ((*pos).weight > r_shortcuts_v[i].weight)
							(*pos).weight = r_shortcuts_v[i].weight;
						r_add_in_v[i] = false;
					}

					CHGraph::CH_Edge dummy(v, 0, 0);
					pos = lower_bound(chgraph.adj[w].begin(), chgraph.adj[w].end(), dummy);
					if (pos != chgraph.adj[w].end() && (*pos).target == v) {
						if ((*pos).weight > r_shortcuts_v[i].weight)
							(*pos).weight = r_shortcuts_v[i].weight;
						r_add_out_v[i] = false;
					}
				}
			}
		}

		// Append the actual shortcuts.
		for (NodeID v = 0; v < numOfVertices; ++v) {

			{//Forward shortcuts, add (v->w) to adj[v] and r_adj[w].
				vector<CHGraph::CH_Edge>& shortcuts_v = shortcuts[v];
				vector<bool>& add_out_v = add_out[v];
				vector<bool>& add_in_v = add_in[v];
				for (int i = 0; i < shortcuts_v.size(); ++i) {
					NodeID w = shortcuts_v[i].target;
					NodeID level = shortcuts_v[i].level;
					EdgeWeight weight = shortcuts_v[i].weight;
					if (add_out_v[i] == true) {
						chgraph.adj[v].push_back(CHGraph::CH_Edge(w, level, weight));
					//	cout << v << "," << w << "," << weight << endl;
					}
					if (add_in_v[i] == true) {
						chgraph.r_adj[w].push_back(CHGraph::CH_Edge(v, level, weight));
					//	cout << w << "," << v << "," << weight << endl;
					}
					if (add_out_v[i] == false && add_in_v[i] == false)
						added_shortcuts--;
					added_shortcuts++;
				}
			}

			{// Backward shortcuts, add (w->v) to adj[w] and r_adj[v].
				vector<CHGraph::CH_Edge>& r_shortcuts_v = r_shortcuts[v];
				vector<bool>& r_add_out_v = r_add_out[v];
				vector<bool>& r_add_in_v = r_add_in[v];
				for (int i = 0; i < r_shortcuts_v.size(); ++i) {
					NodeID w = r_shortcuts_v[i].target;
					NodeID level = r_shortcuts_v[i].level;
					EdgeWeight weight = r_shortcuts_v[i].weight;
					if (r_add_out_v[i] == true) {
						chgraph.adj[w].push_back(CHGraph::CH_Edge(v, level, weight));
					//	cout << w << "," << v << "," << weight << endl;
					}
					if (r_add_in_v[i] == true) {
						chgraph.r_adj[v].push_back(CHGraph::CH_Edge(w, level, weight));
					//	cout << v << "," << w << "," << weight << endl;
					}
					if (r_add_out_v[i] == false && r_add_in_v[i] == false)
						added_shortcuts--;
					added_shortcuts++;
				}
			}

		}
		return added_shortcuts;
	}

};

#endif
