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
#ifndef GRAPH_H
#define GRAPH_H

#include<vector>
#include<fstream>
#include <algorithm>
#include<cstdint> 
#include "paras.h"
#include<iostream>
#include<assert.h>
#include<map>
#include<cmath>


using namespace std;

//typedef unsigned int NodeID;
typedef int NodeID;
typedef long long EdgeID;
//typedef unsigned int EdgeID;
typedef pair<NodeID, NodeID> Edge;
typedef double EdgeWeight;
//typedef  int EdgeWeight;
typedef pair<NodeID, EdgeWeight> NodeEdgeWeightPair;



#define numOfVertices SP_Constants::numOfVertices
#define numOfEdges SP_Constants::numOfEdges
#define DIRECTED_FLAG SP_Constants::DIRECTED_FLAG

// Unweighted Graph.
class Graph {
public:
	

	vector<vector<NodeID> > adj; // Adjacent Lists.
	vector<vector<NodeID> > r_adj; // Reverse adjacent lists for directed search.
	vector<EdgeID> vertices;
	vector<EdgeID> r_vertices;
	vector<NodeID> edges;
	vector<NodeID> r_edges;


	Graph() {
		numOfVertices = 0;
		numOfEdges = 0;
		adj.clear();
		r_adj.clear();
	}

	~Graph() {
		adj.clear();
		r_adj.clear();
	}

	// Load Graph file with format:
	// u	v
	// node id in the graph file should start with 0.
	// The graph file should imply whether it is directed. If it is directed, we will duplicaate (u ,v) with (v, u).
	bool load_graph(const char* graph_file) {
		
		ifstream ifs(graph_file);

		vector<Edge> es;


		for (NodeID u, v; ifs >> u >> v;) {
			numOfVertices = max(numOfVertices, max(u, v) + 1);
			if (u == v) continue; // Delete self loop.
			es.push_back(make_pair(u, v));
		}
		numOfEdges = es.size();
		if (DIRECTED_FLAG == true)
			numOfEdges += numOfEdges;

		if (ifs.bad()) return false;
		ifs.close();

		adj.resize(numOfVertices);
		if (DIRECTED_FLAG == true) 	r_adj.resize(numOfVertices);
		

		for (size_t i = 0; i < es.size(); ++i) {
			NodeID u = es[i].first, v = es[i].second;
			adj[u].push_back(v);
			if (DIRECTED_FLAG == false)
				adj[v].push_back(u);
			else {
				r_adj[v].push_back(u);
			}
		}

		for (NodeID v = 0; v < numOfVertices; ++v) {
			sort(adj[v].begin(), adj[v].end());
			if(DIRECTED_FLAG == true)
				sort(r_adj[v].begin(), r_adj[v].end());
		}

		long long sum_of_edges = 0;
		for (NodeID v = 0; v < numOfVertices; ++v) {
			sum_of_edges += adj[v].size();
		}

		vertices.resize(numOfVertices + 1);
		edges.reserve(sum_of_edges);
		for (NodeID v = 0; v < numOfVertices; ++v) {
			vertices[v] = edges.size();

			for (NodeID i = 0; i < adj[v].size(); ++i)
				edges.push_back(adj[v][i]);
		}

		vertices[numOfVertices] = edges.size();

		if (DIRECTED_FLAG == true) {
			sum_of_edges = 0;
			for (NodeID v = 0; v < numOfVertices; ++v) {
				sum_of_edges += r_adj[v].size();
			}
			r_vertices.resize(numOfVertices + 1);
			r_edges.reserve(sum_of_edges);

			for (NodeID v = 0; v < numOfVertices; ++v) {
				r_vertices[v] = r_edges.size();

				for (NodeID i = 0; i < r_adj[v].size(); ++i)
					r_edges.push_back(r_adj[v][i]);
			}
			r_vertices[numOfVertices] = r_edges.size();
		}

		adj.clear();
		r_adj.clear();
		adj.shrink_to_fit();
		r_adj.shrink_to_fit();

		return true;
	}

};

// Weighted Graph.
class WGraph : public Graph {

public:
	vector< vector<EdgeWeight> > adj_weight; // Weights of adjacent lists;
	vector< vector<EdgeWeight> > r_adj_weight; // Weights of reverse adjacent lists;
	vector<NodeEdgeWeightPair> edges;
	vector<NodeEdgeWeightPair> r_edges;


	WGraph() {
		numOfVertices = 0;
		numOfEdges = 0;
		adj.clear();
		adj_weight.clear();
		r_adj.clear();
		r_adj_weight.clear();
	}

	~WGraph() {
		adj.clear();
		adj_weight.clear();
		r_adj.clear();
		r_adj_weight.clear();
	}

	// Load Graph file with format:
	// u	v	w
	// node id in the graph file should start with 0.
	bool load_graph(const char* graph_file) {

		ifstream ifs(graph_file);

		vector<Edge> es;
		vector<EdgeWeight> es_weight;

		NodeID u, v;
		EdgeWeight w;
		for (; ifs >> u >> v >> w;) {
			numOfVertices = max(numOfVertices, max(u, v) + 1);
			if (u == v) continue; // Delete self loop.
			es.push_back(make_pair(u, v));
			es_weight.push_back(w);
		}

		numOfEdges = es.size();
		
		if (ifs.bad()) return false;
		ifs.close();
		
		adj.resize(numOfVertices);
		adj_weight.resize(numOfVertices);
		assert(adj.size() == numOfVertices && adj_weight.size() == numOfVertices && "No enough memory for adj lists!");
		if (DIRECTED_FLAG == true) {
			r_adj.resize(numOfVertices);
			r_adj_weight.resize(numOfVertices);
			assert(r_adj.size() == numOfVertices && r_adj_weight.size() == numOfVertices  && "No enough memory for adj lists!");
		}

		for (size_t i = 0; i < es.size(); ++i) {
			NodeID u = es[i].first, v = es[i].second;
			EdgeWeight w = es_weight[i];
			adj[u].push_back(v);
			adj_weight[u].push_back(w);
			if (DIRECTED_FLAG == false) {
				adj[v].push_back(u);
				adj_weight[v].push_back(w);
			}
			else {
				r_adj[v].push_back(u);
				r_adj_weight[v].push_back(w);
			}
		}

		for (NodeID v = 0; v < numOfVertices; ++v) {
			vector<NodeID>& adj_v = adj[v];
			vector<EdgeWeight>& adj_weight_v = adj_weight[v];
			vector<pair<NodeID, EdgeWeight> > adj_vw_v(adj_v.size());
			for (size_t i = 0; i < adj_v.size(); ++i) {
				NodeID u = adj_v[i];
				EdgeWeight w = adj_weight_v[i];
				adj_vw_v[i] = make_pair(u, w);
			}
			sort(adj_vw_v.begin(), adj_vw_v.end());
			for (size_t i = 0; i < adj_v.size(); ++i) {
				adj_v[i] = adj_vw_v[i].first;
				adj_weight_v[i] = adj_vw_v[i].second;
			}
			adj_vw_v.clear();

		}

		long long sum_of_edges = 0;
		for (NodeID v = 0; v < numOfVertices; ++v) {
			sum_of_edges += adj[v].size();
		}


		vertices.resize(numOfVertices + 1);
		edges.reserve(sum_of_edges);
		EdgeID duplicated = 0;
		for (NodeID v = 0; v < numOfVertices; ++v) {
			vertices[v] = edges.size();

			NodeID last_v = numOfVertices;
			EdgeWeight last_weight = -1;
			for (NodeID i = 0; i < adj[v].size(); ++i) {
				if (adj[v][i] == last_v) {
					if (last_weight <= adj_weight[v][i]) {
						last_v = adj[v][i];
						last_weight = adj_weight[v][i];
						duplicated++;
						continue;
					}
				}
				edges.push_back(NodeEdgeWeightPair(adj[v][i], adj_weight[v][i]));
				last_v = adj[v][i];
				last_weight = adj_weight[v][i];
			}
		}

		vertices[numOfVertices] = edges.size();
		cout << sum_of_edges << " edges in forward adj." << endl;
		cout << edges.size() << " edges inserted." << endl;
		cout << duplicated << " edges deleted." << endl;

		adj.clear();
		r_adj.clear();
		adj.shrink_to_fit();
		r_adj.shrink_to_fit();

		adj_weight.clear();
		r_adj_weight.clear();
		adj_weight.shrink_to_fit();
		r_adj_weight.shrink_to_fit();

		/*for (NodeID v = 0; v < numOfVertices; ++v) {
			vector<NodeID>& adj_v = adj[v];
			vector<EdgeWeight>& adj_weight_v = adj_weight[v];
			for (size_t i = 0; i < adj_v.size(); ++i) {
				NodeID u = adj_v[i];
				EdgeWeight w = adj_weight_v[i];
				int idx = vertices[v] + i;
				NodeID au = edges[idx].first;
				EdgeWeight aw = edges[idx].second;
				if (u != au || w != aw)
					cout << "we got problems." << endl;
			}
			
			if (DIRECTED_FLAG == true) {
				vector<NodeID>& r_adj_v = r_adj[v];
				vector<EdgeWeight>& r_adj_weight_v = r_adj_weight[v];
				for (size_t i = 0; i < r_adj_v.size(); ++i) {
					NodeID u = r_adj_v[i];
					EdgeWeight w = r_adj_weight_v[i];
					int idx = r_vertices[v] + i;
					NodeID au = r_edges[idx].first;
					EdgeWeight aw = r_edges[idx].second;
					if (u != au || w != aw)
						cout << "we got problems." << endl;
				}
			}
		}*/

		return true;
	}
};

// CH Graph.
class CHGraph : public Graph {



public:
	struct CH_Edge {
		NodeID target;
		NodeID level; // 0: the edges in original graphs otherwise they are added because the level-th node is contracted.
		EdgeWeight weight;

		CH_Edge(NodeID v, NodeID level, EdgeWeight weight) :target(v), level(level), weight(weight) {}
		
		bool operator<(const CH_Edge& rhs) const {			
	//		if (target < rhs.target)
				return target < rhs.target;
			/*	return true;
			else if (target == rhs.target)
				if (weight < rhs.weight)
					return true;
				else
					return false;	*/	
		}

		//bool operator=(const CH_Edge& rhs) const {
		//	//return (target == rhs.target);
		//	return (weight < rhs.weight);
		//}

	};

	vector<vector<CH_Edge> > adj;
	vector<vector<CH_Edge> > r_adj;


	CHGraph() {
		numOfVertices = 0;
		numOfEdges = 0;
		adj.clear();
		r_adj.clear();
	}

	~CHGraph() {
		adj.clear();
		r_adj.clear();
	}

	bool load_ch_overlay_graph(const char* ch_graph_file) {
		ifstream ifs(ch_graph_file);
		vector<EdgeID> veid;
		vector<EdgeID> vsid;
		vector<NodeEdgeWeightPair> edge_tmp;
		int tmp_edges_num;
		ifs >> numOfVertices;
		veid.resize(numOfVertices);
		vsid.resize(numOfVertices);
		for (int i = 0; i < numOfVertices; ++i) {
			ifs >> veid[i];
			ifs >> vsid[i];
		}

		ifs >> tmp_edges_num;
		edge_tmp.resize(tmp_edges_num);
		cout << numOfVertices << "\t" << tmp_edges_num << endl;
		numOfEdges = tmp_edges_num;
		for (int i = 0; i < tmp_edges_num; ++i) {
			NodeID target;
			EdgeWeight weight;
			ifs >> target >> weight;
			edge_tmp[i] = make_pair(target, weight);
		}

		adj.resize(numOfVertices);
		for (NodeID i = 0; i < numOfVertices; ++i) {
			EdgeID st;
			EdgeID et;
			st = veid[i];
			et = veid[i] + vsid[i];
			int cc = 0;
			map<NodeID, NodeID> in;
			for (EdgeID j = st; j < et; ++j) {
				NodeID target = edge_tmp[j].first;
				EdgeWeight weight = edge_tmp[j].second;
				if (target < numOfVertices) {
					if (in.find(target) == in.end()) {
						adj[i].push_back(CH_Edge(target, cc++, weight));
						in[target] = adj[i].size() - 1;
					}
					else {
						if (weight < adj[i][in[target]].weight)
							adj[i][in[target]].weight = weight;
					}
				}
			}
			in.clear();
		
			if(adj[i].size() != 0 )
				sort(adj[i].begin(), adj[i].end());
		}

	}

	bool load_graph(const char* graph_file) {
		ifstream ifs(graph_file);

		vector<Edge> es;
		
		for (NodeID u, v; ifs >> u >> v;) {
			numOfVertices = max(numOfVertices, max(u, v) + 1);
			if (u == v) continue;
			es.push_back(make_pair(u, v));
		}
		numOfEdges = es.size();

		if (ifs.bad()) return false;
		ifs.close();

		adj.resize(numOfVertices);
		r_adj.resize(numOfVertices);


		for (size_t i = 0; i < es.size(); ++i) {
			NodeID u = es[i].first, v = es[i].second;

			adj[u].push_back(CH_Edge(v,0,1));
			if (DIRECTED_FLAG == false) {
				adj[v].push_back(CH_Edge(u, 0, 1));
				r_adj[v].push_back(CH_Edge(u, 0, 1));
				r_adj[u].push_back(CH_Edge(v, 0, 1));
			} else {
				r_adj[v].push_back(CH_Edge(u,0,1));
			}
		}
		return true;

	}

	bool load_wgraph(const char* graph_file) {
		ifstream ifs(graph_file);

		vector<Edge> es;
		vector<EdgeWeight> es_weight;

		NodeID u, v;
		EdgeWeight w;
		for (; ifs >> u >> v >> w;) {
			numOfVertices = max(numOfVertices, max(u, v) + 1);
			if (u == v) continue;
			es.push_back(make_pair(u, v));
			es_weight.push_back(w);
		}

		numOfEdges = es.size();

		if (ifs.bad()) return false;
		ifs.close();

		adj.resize(numOfVertices);
		r_adj.resize(numOfVertices);
		

		for (size_t i = 0; i < es.size(); ++i) {
			NodeID u = es[i].first, v = es[i].second;
			EdgeWeight w = es_weight[i];
			adj[u].push_back(CH_Edge(v, 0, w));
			if (DIRECTED_FLAG == false) {
				adj[v].push_back(CH_Edge(u, 0, w));
				r_adj[v].push_back(CH_Edge(u, 0, w));
				r_adj[u].push_back(CH_Edge(v, 0, w));
			}
			else {
				r_adj[v].push_back(CH_Edge(u, 0, w));
			}
		}
		return true;
	}

};

#endif
