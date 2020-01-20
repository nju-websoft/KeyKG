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
#ifndef GRAPH_SEARCH_H
#define GRAPH_SEARCH_H

#include "graph.h"
#include "heap.h"
#include "paras.h"
#include <queue>
#include <map>
#define INF_WEIGHT SP_Constants::INF_WEIGHT

namespace graph_search {

	void BFS(vector<EdgeWeight> &distances, vector<NodeID> &que, NodeID& que_h, vector<bool> &vis, NodeID &source, Graph &graph) {
		
		NodeID que_t0 = 0, que_t1 = 0;
		que_h = 0;

		//vector<vector<NodeID> >& adj = graph.adj;
		vector<EdgeID>& vertices = graph.vertices;
		vector<NodeID>& edges = graph.edges;

		que[que_h++] = source;
		vis[source] = true;
		distances[source] = 0;
		que_t1 = que_h;
		
		for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID v = que[que_i];
				distances[v] = d;

	/*			for (size_t i = 0; i < adj[v].size(); ++i) {
					NodeID w = adj[v][i];*/
				//Array Representation
				for (EdgeID eid = vertices[v]; eid < vertices[v + 1]; ++eid){

					NodeID w = edges[eid];

					if (!vis[w]) {
						que[que_h++] = w;
						vis[w] = true;
					}
				}
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}
	}	

	void BFS_with_target(vector<EdgeWeight> &distances, vector<NodeID> &que, NodeID& que_h, vector<bool> &vis, NodeID &source, Graph &graph, NodeID target, EdgeWeight& dis) {

		NodeID que_t0 = 0, que_t1 = 0;
		que_h = 0;

		//vector<vector<NodeID> >& adj = graph.adj;
		vector<EdgeID>& vertices = graph.vertices;
		vector<NodeID>& edges = graph.edges;

		que[que_h++] = source;
		vis[source] = true;
		distances[source] = 0;
		que_t1 = que_h;

		for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID v = que[que_i];
				distances[v] = d;

				if (v == target) {
					dis = d;
					return;
				}

				/*			for (size_t i = 0; i < adj[v].size(); ++i) {
				NodeID w = adj[v][i];*/
				//Array Representation
				for (EdgeID eid = vertices[v]; eid < vertices[v + 1]; ++eid) {

					NodeID w = edges[eid];

					if (!vis[w]) {
						que[que_h++] = w;
						vis[w] = true;
					}
				}
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}
	}
	
	void backward_BFS(vector<EdgeWeight> &distances, vector<NodeID> &que, NodeID& que_h, vector<bool> &vis, NodeID &source, Graph &graph) {

		NodeID que_t0 = 0, que_t1 = 0;
		que_h = 0;

		//vector<vector<NodeID> >& adj = graph.adj;
		vector<EdgeID>& r_vertices = graph.r_vertices;
		vector<NodeID>& r_edges = graph.r_edges;

		que[que_h++] = source;
		vis[source] = true;
		distances[source] = 0;
		que_t1 = que_h;

		for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID v = que[que_i];
				distances[v] = d;

				/*			for (size_t i = 0; i < adj[v].size(); ++i) {
				NodeID w = adj[v][i];*/
				//Array Representation
				for (EdgeID eid = r_vertices[v]; eid < r_vertices[v + 1]; ++eid) {

					NodeID w = r_edges[eid];

					if (!vis[w]) {
						que[que_h++] = w;
						vis[w] = true;
					}
				}
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}
	}


	void BFS_init(vector<EdgeWeight> &distances, vector<NodeID> &que, NodeID& que_h, vector<bool> &vis) {
		distances.resize(numOfVertices, INF_WEIGHT);
		que.resize(numOfVertices);
		que_h = 0;
		vis.resize(numOfVertices);
	}
	
	// Clear after each single-source shortest path. 
	void BFS_clear(vector<EdgeWeight> &distances, vector<NodeID> &que, NodeID& que_h, vector<bool> &vis) {
		for (int i = 0; i < que_h; ++i) {
			vis[que[i]] = false;
			distances[que[i]] = INF_WEIGHT;
		}
		que_h = 0;
	}
	
	EdgeWeight BFS(NodeID& source, NodeID& target, Graph& graph) {
		EdgeWeight d = INF_WEIGHT;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
		vector<NodeID> que;
		vector<bool> vis;
		NodeID que_h = 0;

		BFS_init(distances, que, que_h, vis);
		BFS_with_target(distances, que, que_h, vis, source, graph, target, d);

		BFS_clear(distances, que, que_h, vis);
		return d;
	}


	EdgeWeight BiBFS(NodeID& source, NodeID& target, Graph& graph) {
		EdgeWeight d = INF_WEIGHT;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
		vector<NodeID> que;
		vector<bool> vis;
		NodeID que_h = 0;

		vector<EdgeWeight> r_distances(numOfVertices, INF_WEIGHT);
		vector<NodeID> r_que;
		vector<bool> r_vis;
		NodeID r_que_h = 0;


		BFS_init(distances, que, que_h, vis);
		BFS_init(r_distances, r_que, r_que_h, r_vis);

		NodeID que_t0 = 0, que_t1 = 0;
		NodeID r_que_t0 = 0, r_que_t1 = 0;
		que_h = 0;
		r_que_h = 0;

		//vector<vector<NodeID> >& adj = graph.adj;
		vector<EdgeID>& vertices = graph.vertices;
		vector<NodeID>& edges = graph.edges;
		vector<EdgeID>& r_vertices = graph.r_vertices;
		vector<NodeID>& r_edges = graph.r_edges;

		que[que_h++] = source;
		vis[source] = true;
		distances[source] = 0;
		que_t1 = que_h;

		r_que[r_que_h++] = target;
		r_vis[target] = true;
		r_distances[target] = 0;
		r_que_t1 = r_que_h;

		EdgeWeight meet_distance = INF_WEIGHT;
		bool fin_flag = false;

		for (EdgeWeight d = 0; que_t0 < que_h && r_que_t0 < r_que_h; d = d + 1) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID v = que[que_i];
				distances[v] = d;

				if (v == target) {
					meet_distance = distances[v];
					fin_flag = true;
					break;
				}

				if (r_distances[v] < INF_WEIGHT) {
					if (meet_distance > distances[v] + r_distances[v])
						meet_distance = distances[v] + r_distances[v];
				}

				/*			for (size_t i = 0; i < adj[v].size(); ++i) {
				NodeID w = adj[v][i];*/
				//Array Representation
				for (EdgeID eid = vertices[v]; eid < vertices[v + 1]; ++eid) {

					NodeID w = edges[eid];

					if (!vis[w]) {
						que[que_h++] = w;
						vis[w] = true;
					}
				}
			}
			que_t0 = que_t1;
			que_t1 = que_h;

			if (fin_flag == true)
				break;

			for (NodeID r_que_i = r_que_t0; r_que_i < r_que_t1; ++r_que_i) {
				NodeID v = r_que[r_que_i];
				r_distances[v] = d;

				if (v == source) {
					meet_distance = r_distances[v];
					fin_flag = true;
					break;
				}

				if (distances[v] < INF_WEIGHT) {
					if (meet_distance > distances[v] + r_distances[v])
						meet_distance = distances[v] + r_distances[v];
				}

				/*			for (size_t i = 0; i < adj[v].size(); ++i) {
				NodeID w = adj[v][i];*/
				//Array Representation
				for (EdgeID eid = vertices[v]; eid < vertices[v + 1]; ++eid) {

					NodeID w = edges[eid];

					if (!r_vis[w]) {
						r_que[que_h++] = w;
						r_vis[w] = true;
					}
				}
			}
			r_que_t0 = r_que_t1;
			r_que_t1 = r_que_h;


			if (fin_flag == true)
				break;

		}




		BFS_clear(distances, que, que_h, vis);
		BFS_clear(r_distances, r_que, r_que_h, r_vis);
		return meet_distance;
	}

	void BFS_rank(const char* query_file, int num_sources, Graph& graph) {

	
		srand(time(NULL));

		map<EdgeWeight, vector<pair<NodeID, NodeID> > > distance_pairs;
		map<EdgeWeight, vector<pair<NodeID, NodeID> > > bfs_rank;

		for (NodeID i = 0; i < num_sources; ++i) {
			NodeID source = rand() % numOfVertices;

			EdgeWeight d = INF_WEIGHT;
			vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
			vector<NodeID> que;
			vector<bool> vis;
			NodeID que_h = 0;
			BFS_init(distances, que, que_h, vis);
			BFS(distances, que, que_h, vis, source, graph);


			for (NodeID j = 0; j < numOfVertices; ++j) {
				if (distances[j] > 0 && distances[j] < 999999) {
					distance_pairs[distances[j]].push_back(make_pair(source, j));					
				}
			}
			BFS_clear(distances, que, que_h, vis);
		}
		cout << "Collecting Info..." << endl;
		srand(time(NULL));
		for (map<EdgeWeight, vector<pair<NodeID, NodeID> > >::iterator it = distance_pairs.begin(); it != distance_pairs.end(); ++it) {
			int limit = (*it).second.size() > num_sources ? num_sources : (*it).second.size();
			if (limit < num_sources) {
				for (NodeID i = 0; i < (*it).second.size(); ++i)
					bfs_rank[(*it).first].push_back((*it).second[i]);
			}
			else {
				while (limit > 0) {
					NodeID pid = rand() % (*it).second.size();
					bfs_rank[(*it).first].push_back((*it).second[pid]);
					limit--;
				}
			}
		}

		ofstream ofs(query_file);
		for (map<EdgeWeight, vector<pair<NodeID, NodeID> > >::iterator it = bfs_rank.begin(); it != bfs_rank.end(); ++it) {
			ofs << (*it).first << "\t" << (*it).second.size();
			for (NodeID i = 0; i < (*it).second.size(); ++i) {
				ofs << "\t" << (*it).second[i].first << "\t" << (*it).second[i].second;
			}
			ofs << endl;
		}
	/*	for (NodeID i = 0; i < num_sources; ++i) {
			ofs << sources[i] << endl;
			ofs << bfs_rank[i].size() << endl;
			for (map<EdgeWeight, unordered_set<NodeID> >::iterator it = bfs_rank[i].begin(); it != bfs_rank[i].end(); ++it) {
				ofs << (*it).first << "\t" << (*it).second.size();
				for (unordered_set<NodeID>::iterator sit = (*it).second.begin(); sit != (*it).second.end(); ++sit) {
					ofs << "\t" << (*sit);
				}
				ofs << endl;
			}
		}*/
		ofs.close();

	}

	

	void Dijkstra(vector<EdgeWeight> &distances, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, queue<NodeID> &visited_que, vector<bool>& vis, NodeID& source, WGraph& wgraph) {

	/*	vector<vector<NodeID> >& adj = wgraph.adj;
		vector<vector<EdgeWeight> >& adj_weight = wgraph.adj_weight;*/
		vector<EdgeID>& vertices = wgraph.vertices;
		vector<NodeEdgeWeightPair>& edges = wgraph.edges;

		pqueue.update(source, 0);
		
		while (!pqueue.empty()) {
			NodeID v;
			EdgeWeight v_d;
			pqueue.extract_min(v, v_d);
			distances[v] = v_d;
			vis[v] = true;
			visited_que.push(v);



		/*	for (size_t i = 0; i < adj[v].size(); ++i) {
				NodeID w = adj[v][i];

				EdgeWeight w_d = adj_weight[v][i] + v_d;*/
			
			// Array Representation
			for (EdgeID eid = vertices[v]; eid < vertices[v + 1]; ++eid) {
				NodeID w = edges[eid].first;
				EdgeWeight w_d = edges[eid].second + v_d;
				if (!vis[w]) {
					if (distances[w] > w_d) {
						distances[w] = w_d;
						pqueue.update(w, w_d);
					}
				}
			}
		}
	}

	void backward_Dijkstra(vector<EdgeWeight> &distances, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, queue<NodeID> &visited_que, vector<bool>& vis, NodeID& source, WGraph& wgraph) {

		//vector<vector<NodeID> >& adj = wgraph.r_adj;
		//vector<vector<EdgeWeight> >& adj_weight = wgraph.r_adj_weight;
		
		vector<EdgeID>& r_vertices = wgraph.r_vertices;
		vector<NodeEdgeWeightPair>& r_edges = wgraph.r_edges;

		pqueue.update(source, 0);

		while (!pqueue.empty()) {
			NodeID v;
			EdgeWeight v_d;
			pqueue.extract_min(v, v_d);
			distances[v] = v_d;
			vis[v] = true;
			visited_que.push(v);

			//for (size_t i = 0; i < adj[v].size(); ++i) {
			//	NodeID w = adj[v][i];

			//	EdgeWeight w_d = adj_weight[v][i] + v_d;

			// Array Representation

			for (EdgeID eid = r_vertices[v]; eid < r_vertices[v + 1]; ++eid) {

				NodeID w = r_edges[eid].first;
				EdgeWeight w_d = r_edges[eid].second + v_d;

				if (!vis[w]) {
					if (distances[w] > w_d) {
						distances[w] = w_d;
						pqueue.update(w, w_d);
					}
				}
			}
		}
	}

	void Dijkstra_init(vector<EdgeWeight> &distances, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, queue<NodeID> &visited_que, vector<bool>& vis) {
		distances.resize(numOfVertices, INF_WEIGHT);
		pqueue.clear();
		while (!visited_que.empty()) visited_que.pop();
		vis.resize(numOfVertices);
	}

	void Dijkstra_clear(vector<EdgeWeight> &distances, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, queue<NodeID> &visited_que, vector<bool>& vis) {
		while (!visited_que.empty()) {
			NodeID vis_v = visited_que.front();
			distances[vis_v] = INF_WEIGHT;
			vis[vis_v] = false;
			pqueue.clear(vis_v);
			visited_que.pop();
		}
	}

	EdgeWeight Dijkstra(NodeID& source, NodeID& target, WGraph& wgraph) {
		EdgeWeight d = INF_WEIGHT;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
		vector<NodeID> que;
		vector<bool> vis;
		NodeID que_h = 0;
		queue<NodeID> visited_que;
		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);


		Dijkstra_init(distances, pqueue, visited_que, vis);
		
		
		
		Dijkstra(distances, pqueue, visited_que, vis, source, wgraph);




		d = distances[target];
		Dijkstra_clear(distances, pqueue, visited_que, vis);
		return d;
	}

	EdgeWeight BiDijkstra(NodeID& source, NodeID& target, WGraph& wgraph) {
		EdgeWeight d = INF_WEIGHT;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
		vector<bool> vis;
		queue<NodeID> visited_que;
		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);

		vector<EdgeWeight> r_distances(numOfVertices, INF_WEIGHT);
		vector<bool> r_vis;
		queue<NodeID> r_visited_que;
		benchmark::heap<2, EdgeWeight, NodeID> r_pqueue(numOfVertices);

		Dijkstra_init(distances, pqueue, visited_que, vis);
		Dijkstra_init(r_distances, r_pqueue, r_visited_que, r_vis);

		vector<EdgeID>& vertices = wgraph.vertices;
		vector<NodeEdgeWeightPair>& edges = wgraph.edges;


		vector<EdgeID>& r_vertices = wgraph.r_vertices;
		vector<NodeEdgeWeightPair>& r_edges = wgraph.r_edges;

		pqueue.update(source, 0);
		r_pqueue.update(target, 0);

		while (pqueue.top_value() + r_pqueue.top_value() <= d) {
			NodeID v;
			EdgeWeight v_d;
			pqueue.extract_min(v, v_d);
			distances[v] = v_d;

			if (r_distances[v] < INF_WEIGHT)
				if (d > distances[v] + r_distances[v])
					d = distances[v] + r_distances[v];

			if (v == target) break;

			vis[v] = true;
			visited_que.push(v);
			// Array Representation
			for (EdgeID eid = vertices[v]; eid < vertices[v + 1]; ++eid) {
				NodeID w = edges[eid].first;
				EdgeWeight w_d = edges[eid].second + v_d;
				if (!vis[w]) {
					if (distances[w] > w_d) {
						distances[w] = w_d;
						pqueue.update(w, w_d);
					}
				}
			}

			//Backward Search from target.
			r_pqueue.extract_min(v, v_d);
			r_distances[v] = v_d;

			if (distances[v] < INF_WEIGHT)
				if (d > distances[v] + r_distances[v])
					d = distances[v] + r_distances[v];

			if (v == source) break;

			r_vis[v] = true;
			r_visited_que.push(v);
			if (DIRECTED_FLAG == false) {
				for (EdgeID eid = vertices[v]; eid < vertices[v + 1]; ++eid) {
					NodeID w = edges[eid].first;
					EdgeWeight w_d = edges[eid].second + v_d;
					if (!r_vis[w]) {
						if (r_distances[w] > w_d) {
							r_distances[w] = w_d;
							r_pqueue.update(w, w_d);
						}
					}
				}
			}
			else {
				for (EdgeID eid = r_vertices[v]; eid < r_vertices[v + 1]; ++eid) {
					NodeID w = r_edges[eid].first;
					EdgeWeight w_d = r_edges[eid].second + v_d;
					if (!r_vis[w]) {
						if (r_distances[w] > w_d) {
							r_distances[w] = w_d;
							r_pqueue.update(w, w_d);
						}
					}
				}
			}

		}
				
	//	Dijkstra(distances, pqueue, visited_que, vis, source, wgraph);

		Dijkstra_clear(distances, pqueue, visited_que, vis);
		Dijkstra_clear(r_distances, r_pqueue, r_visited_que, r_vis);
		return d;
	}


	void Dijkstra_rank(const char* query_file, int num_sources, WGraph& wgraph) {


		srand(time(NULL));

		map<EdgeWeight, vector<pair<NodeID, NodeID> > > dijrank_pairs;
		map<EdgeWeight, vector<pair<NodeID, NodeID> > > dij_rank;

		for (NodeID i = 0; i < num_sources; ++i) {
			NodeID source = rand() % numOfVertices;
			EdgeWeight d = INF_WEIGHT;

			vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
			vector<NodeID> que;
			vector<bool> vis;
			NodeID que_h = 0;
			queue<NodeID> visited_que;
			benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);


			Dijkstra_init(distances, pqueue, visited_que, vis);
			Dijkstra(distances, pqueue, visited_que, vis, source, wgraph);

			map < EdgeWeight, vector<pair<NodeID, NodeID> > > this_dij_rank;
			for (NodeID j = 0; j < numOfVertices; ++j) {
				if (distances[j] > 0 && distances[j] < 99999999) {
					this_dij_rank[distances[j]].push_back(make_pair(source, j));
				}
			}

			int rank = 0;
			for (map < EdgeWeight, vector<pair<NodeID, NodeID> > >::iterator it = this_dij_rank.begin(); it != this_dij_rank.end(); ++it) {
				rank++;
				for (NodeID j = 0; j < (*it).second.size(); ++j)
					dijrank_pairs[rank].push_back((*it).second[j]);

			}
			Dijkstra_clear(distances, pqueue, visited_que, vis);
		}
		cout << "Collecting Info..." << endl;
		srand(time(NULL));
		for (map<EdgeWeight, vector<pair<NodeID, NodeID> > >::iterator it = dijrank_pairs.begin(); it != dijrank_pairs.end(); ++it) {
			int limit = (*it).second.size() > num_sources ? num_sources : (*it).second.size();
			if (limit < num_sources) {
				for (NodeID i = 0; i < (*it).second.size(); ++i)
					dij_rank[(*it).first].push_back((*it).second[i]);
			}
			else {
				while (limit > 0) {
					NodeID pid = rand() % (*it).second.size();
					dij_rank[(*it).first].push_back((*it).second[pid]);
					limit--;
				}
			}
		}

		ofstream ofs(query_file);
		for (map<EdgeWeight, vector<pair<NodeID, NodeID> > >::iterator it = dij_rank.begin(); it != dij_rank.end(); ++it) {
			ofs << (*it).first << "\t" << (*it).second.size();
			for (NodeID i = 0; i < (*it).second.size(); ++i) {
				ofs << "\t" << (*it).second[i].first << "\t" << (*it).second[i].second;
			}
			ofs << endl;
		}
		/*	for (NodeID i = 0; i < num_sources; ++i) {
		ofs << sources[i] << endl;
		ofs << bfs_rank[i].size() << endl;
		for (map<EdgeWeight, unordered_set<NodeID> >::iterator it = bfs_rank[i].begin(); it != bfs_rank[i].end(); ++it) {
		ofs << (*it).first << "\t" << (*it).second.size();
		for (unordered_set<NodeID>::iterator sit = (*it).second.begin(); sit != (*it).second.end(); ++sit) {
		ofs << "\t" << (*sit);
		}
		ofs << endl;
		}
		}*/
		ofs.close();

	}

	void CH_Dijkstra(vector<EdgeWeight> &distances, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, queue<NodeID> &visited_que, vector<bool>& vis, NodeID& source, const vector<vector<CHGraph::CH_Edge> > adj) {

		pqueue.update(source, 0);

		while (!pqueue.empty()) {
			NodeID v;
			EdgeWeight v_d;
			pqueue.extract_min(v, v_d);
			distances[v] = v_d;
			vis[v] = true;
			visited_que.push(v);

			for (size_t i = 0; i < adj[v].size(); ++i) {
				NodeID w = adj[v][i].target;
				if (w > v) continue; // CH search, only search the upper graph.
				
				EdgeWeight w_d = adj[v][i].weight + v_d;
				if (!vis[w]) {
					if (distances[w] > w_d) {
						distances[w] = w_d;
						pqueue.update(w, w_d);
					}
				}
			}
		}
	}

	EdgeWeight Dijkstra(NodeID& source, NodeID& target, CHGraph& chgraph) {
		EdgeWeight d = INF_WEIGHT;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
		vector<bool> vis;
		queue<NodeID> visited_que;
		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);

		vector<EdgeWeight> r_distances(numOfVertices, INF_WEIGHT);
		vector<bool> r_vis;
		queue<NodeID> r_visited_que;
		benchmark::heap<2, EdgeWeight, NodeID> r_pqueue(numOfVertices);

		Dijkstra_init(distances, pqueue, visited_que, vis);

		Dijkstra_init(r_distances, r_pqueue, r_visited_que, r_vis);

		CH_Dijkstra(distances, pqueue, visited_que, vis, source, chgraph.adj);

		CH_Dijkstra(r_distances, r_pqueue, r_visited_que, r_vis, target, chgraph.r_adj);

		for (int i = 0; i < numOfVertices; ++i) {
			if (distances[i] < INF_WEIGHT && r_distances[i] < INF_WEIGHT)
				if (d > distances[i] + r_distances[i])
					d = distances[i] + r_distances[i];
		}

		Dijkstra_clear(distances, pqueue, visited_que, vis);
		Dijkstra_clear(r_distances, r_pqueue, r_visited_que, r_vis);
		return d;
	}

	void read_query_rank(const char* rank_query_file, map<NodeID, vector<pair<NodeID, NodeID> > >& rank_query_pair) {
		ifstream ifs(rank_query_file);
		NodeID rank = 0;
		NodeID num = 0;
		while (ifs >> rank >> num) {
			NodeID s, t;
			for (NodeID i = 0; i < num; ++i) {
				ifs >> s >> t;
				rank_query_pair[rank].push_back(make_pair(s, t));
			}
		}
		ifs.close();
	}

};

#endif