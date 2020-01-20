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

#ifndef LABELS_H
#define LABELS_H


#include <limits>
#include <climits>
#include <stdlib.h>
#include <iostream>
#include <sys/time.h>
#include "graph.h"
#include "paras.h" 
#include <malloc.h>
#include <xmmintrin.h>
//typedef unsigned __int64 BPSeed;
#include <omp.h>
#include<bitset>

#define numOfVertices SP_Constants::numOfVertices
#define numOfEdges SP_Constants::numOfEdges
#define INF_WEIGHT SP_Constants::INF_WEIGHT

struct index_t {
	vector<NodeID> spt_v;
	vector<EdgeWeight> spt_d;

	NodeID size() {
		return spt_v.size();
	}

};

struct index_t_p {
	NodeID* spt_v;
	EdgeWeight* spt_d;
}__attribute__((aligned(64)));  // Aligned for cache lines;


struct two_index_t_p {
	NodeID* spt_v;
	EdgeWeight* spt_d;
	uint8_t* spt_lv;
	EdgeWeight* spt_ld;
}__attribute__((aligned(64)));  // Aligned for cache lines;

struct index_t_path {
	vector<NodeID> spt_v;
	vector<NodeID> spt_p;//parent nodes
	vector<EdgeWeight> spt_d;

	NodeID size() {
		return spt_v.size();
	}

};

struct index_t_path_p {
	NodeID* spt_v;
	NodeID* spt_p;
	EdgeWeight* spt_d;
};

struct query_info {
	NodeID meet_node;
	NodeID search_len;
	double time_cost;
	EdgeWeight distance;
};

template<int kNumBitParallelRoots = 50>
struct index_t_bp {
	NodeID* spt_v;
	EdgeWeight* spt_d;
	EdgeWeight bpspt_d[kNumBitParallelRoots];
	uint64_t bpspt_s[kNumBitParallelRoots][2];
}__attribute__((aligned(64)));  // Aligned for cache lines;


struct token_t {
	NodeID* sptc_v; // sptc_v[0] is the root
	EdgeWeight* sptc_d;	 // |*| = k + 1, sptc_d[0] is the number of children - k
	unsigned char* sptc_fbv; // first-level bit vector
	unsigned char* sptc_sbv; // second-level bit vector
	NodeID* sptc_pathv; // intermediate point for a path
}__attribute__((aligned(64)));

class CLabel {

public:
	token_t* supertokenindex_p;
	token_t* tokenindex_p;
	NodeID* anchor_p;
	NodeID numOfTokens;
	long total_children;
	
	
	token_t* r_supertokenindex_p;
	token_t* r_tokenindex_p;
	NodeID* r_anchor_p;
	NodeID r_numOfTokens;
	long r_total_children;
	
	void save_labels(const char* save_filename) {
		ofstream ofs(save_filename, ios::binary | ios::out);

		ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
		for (NodeID v = 0; v < numOfVertices; ++v) {
			ofs.write((const char*)&anchor_p[v], sizeof(anchor_p[v]));
		//	ofs.write((const char*)&index_[v].spt_d[i], sizeof(index_[v].spt_d[i]));
		}
		
		ofs.write((const char*)&numOfTokens, sizeof(numOfTokens));
		for (NodeID t = 0; t < numOfTokens; ++t) {
			token_t& tt = tokenindex_p[t];
			EdgeWeight tsize = tt.sptc_d[0];
			ofs.write((const char*)&tt.sptc_v[0], sizeof(tt.sptc_v[0]));
			ofs.write((const char*)&tsize, sizeof(tsize));
			for(NodeID c = 0; c < tsize; ++c){				
				ofs.write((const char*)&tt.sptc_v[1 + c], sizeof(tt.sptc_v[1 + c]));
				ofs.write((const char*)&tt.sptc_d[1 + c], sizeof(tt.sptc_d[1 + c]));
			} 
		} 
		
		ofs.close();
	}
	
	
	void save_labels_path(const char* save_filename) {
		ofstream ofs(save_filename, ios::binary | ios::out);

		ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
		for (NodeID v = 0; v < numOfVertices; ++v) {
			ofs.write((const char*)&anchor_p[v], sizeof(anchor_p[v]));
		//	ofs.write((const char*)&index_[v].spt_d[i], sizeof(index_[v].spt_d[i]));
		}
		
		ofs.write((const char*)&numOfTokens, sizeof(numOfTokens));
		for (NodeID t = 0; t < numOfTokens; ++t) {
			token_t& tt = tokenindex_p[t];
			EdgeWeight tsize = tt.sptc_d[0];
			ofs.write((const char*)&tt.sptc_v[0], sizeof(tt.sptc_v[0]));
			ofs.write((const char*)&tsize, sizeof(tsize));
			for(NodeID c = 0; c < tsize; ++c){				
				ofs.write((const char*)&tt.sptc_v[1 + c], sizeof(tt.sptc_v[1 + c]));
				ofs.write((const char*)&tt.sptc_d[1 + c], sizeof(tt.sptc_d[1 + c]));
				ofs.write((const char*)&tt.sptc_pathv[1 + c], sizeof(tt.sptc_pathv[1 + c]));
			} 
		} 
		
		ofs.close();
	}
	
	void save_labels_d(const char* save_filename) {
		ofstream ofs(save_filename, ios::binary | ios::out);

		ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
		for (NodeID v = 0; v < numOfVertices; ++v) {
			ofs.write((const char*)&anchor_p[v], sizeof(anchor_p[v]));
		//	ofs.write((const char*)&index_[v].spt_d[i], sizeof(index_[v].spt_d[i]));
		}
		for (NodeID v = 0; v < numOfVertices; ++v) {
			ofs.write((const char*)&r_anchor_p[v], sizeof(r_anchor_p[v]));
		//	ofs.write((const char*)&index_[v].spt_d[i], sizeof(index_[v].spt_d[i]));
		}
		
		ofs.write((const char*)&numOfTokens, sizeof(numOfTokens));
		for (NodeID t = 0; t < numOfTokens; ++t) {
			token_t& tt = tokenindex_p[t];
			EdgeWeight tsize = tt.sptc_d[0];
			ofs.write((const char*)&tt.sptc_v[0], sizeof(tt.sptc_v[0]));
			ofs.write((const char*)&tsize, sizeof(tsize));
			for(NodeID c = 0; c < tsize; ++c){				
				ofs.write((const char*)&tt.sptc_v[1 + c], sizeof(tt.sptc_v[1 + c]));
				ofs.write((const char*)&tt.sptc_d[1 + c], sizeof(tt.sptc_d[1 + c]));
			} 
		}
		
		ofs.write((const char*)&r_numOfTokens, sizeof(r_numOfTokens));
		for (NodeID t = 0; t < r_numOfTokens; ++t) {
			token_t& tt = r_tokenindex_p[t];
			EdgeWeight tsize = tt.sptc_d[0];
			ofs.write((const char*)&tt.sptc_v[0], sizeof(tt.sptc_v[0]));
			ofs.write((const char*)&tsize, sizeof(tsize));
			for(NodeID c = 0; c < tsize; ++c){				
				ofs.write((const char*)&tt.sptc_v[1 + c], sizeof(tt.sptc_v[1 + c]));
				ofs.write((const char*)&tt.sptc_d[1 + c], sizeof(tt.sptc_d[1 + c]));
			} 
		}
		
		ofs.close();
	}
	
	
	void load_labels_path(const char* load_filename) {
		
		total_children = 0;
		
		tokenindex_p = NULL;
		anchor_p = NULL;

		ifstream ifs(load_filename);
		NodeID isize = 0;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;

		anchor_p = (NodeID*)memalign(64, numOfVertices * sizeof(NodeID));
		NodeID anchor_id;
		for (NodeID v = 0; v < numOfVertices; ++v) {
				ifs.read((char*)&anchor_id, sizeof(anchor_id));
				anchor_p[v] = anchor_id;
		}
		
		ifs.read((char*)&isize, sizeof(isize));
		numOfTokens = isize;		
		tokenindex_p = (token_t*)memalign(64, numOfTokens * sizeof(token_t));
		
		EdgeWeight csize;
		NodeID cid;
		EdgeWeight cd;
		for (NodeID v = 0; v < numOfTokens; ++v) {
			token_t& tt = tokenindex_p[v];
		
			ifs.read((char*)&cid, sizeof(cid));
			ifs.read((char*)&csize, sizeof(csize));
			
			tt.sptc_v = (NodeID*)memalign(64, (csize + 1) * sizeof(NodeID));
			tt.sptc_d = (EdgeWeight*)memalign(64, (csize + 1 ) * sizeof(EdgeWeight));
			
			total_children += (csize + 1);
			
			tt.sptc_v[0] = cid;
			tt.sptc_d[0] = csize;
			for (NodeID i = 0; i < csize; ++i) {
				ifs.read((char*)&cid, sizeof(cid));
				ifs.read((char*)&cd, sizeof(cd));
				tt.sptc_v[i + 1] = cid;
				tt.sptc_d[i + 1] = cd;			
			}
		}
		ifs.close();
	}
	
	void load_labels(const char* load_filename) {
		
		total_children = 0;
		
		tokenindex_p = NULL;
		anchor_p = NULL;

		ifstream ifs(load_filename);
		NodeID isize = 0;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;

		anchor_p = (NodeID*)memalign(64, numOfVertices * sizeof(NodeID));
		NodeID anchor_id;
		for (NodeID v = 0; v < numOfVertices; ++v) {
				ifs.read((char*)&anchor_id, sizeof(anchor_id));
				anchor_p[v] = anchor_id;
		}
		
		ifs.read((char*)&isize, sizeof(isize));
		numOfTokens = isize;		
		tokenindex_p = (token_t*)memalign(64, numOfTokens * sizeof(token_t));
		
		EdgeWeight csize;
		NodeID cid;
		EdgeWeight cd;
		for (NodeID v = 0; v < numOfTokens; ++v) {
			token_t& tt = tokenindex_p[v];
		
			ifs.read((char*)&cid, sizeof(cid));
			ifs.read((char*)&csize, sizeof(csize));
			
			tt.sptc_v = (NodeID*)memalign(64, (csize + 1) * sizeof(NodeID));
			tt.sptc_d = (EdgeWeight*)memalign(64, (csize + 1 ) * sizeof(EdgeWeight));
			
			total_children += (csize + 1);
			
			tt.sptc_v[0] = cid;
			tt.sptc_d[0] = csize;
			
			for (NodeID i = 0; i < csize; ++i) {
				ifs.read((char*)&cid, sizeof(cid));
				ifs.read((char*)&cd, sizeof(cd));
				tt.sptc_v[i + 1] = cid;
				tt.sptc_d[i + 1] = cd;		
			}
		}
		ifs.close();
	}
	
	
	void load_labels_d(const char* load_filename) {
		
		total_children = 0;
		r_total_children = 0;
		
		tokenindex_p = NULL;
		anchor_p = NULL;

		r_tokenindex_p = NULL;
		r_anchor_p = NULL;
		
		ifstream ifs(load_filename);
		NodeID isize = 0;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;

		anchor_p = (NodeID*)memalign(64, numOfVertices * sizeof(NodeID));
		r_anchor_p = (NodeID*)memalign(64, numOfVertices * sizeof(NodeID));
		NodeID anchor_id;
		for (NodeID v = 0; v < numOfVertices; ++v) {
				ifs.read((char*)&anchor_id, sizeof(anchor_id));
				anchor_p[v] = anchor_id;
		}
		for (NodeID v = 0; v < numOfVertices; ++v) {
				ifs.read((char*)&anchor_id, sizeof(anchor_id));
				r_anchor_p[v] = anchor_id;
		}
		 
		ifs.read((char*)&isize, sizeof(isize));
		numOfTokens = isize;		
		tokenindex_p = (token_t*)memalign(64, numOfTokens * sizeof(token_t));
		
		EdgeWeight csize;
		NodeID cid;
		EdgeWeight cd;
		for (NodeID v = 0; v < numOfTokens; ++v) {
			token_t& tt = tokenindex_p[v];
		
			ifs.read((char*)&cid, sizeof(cid));
			ifs.read((char*)&csize, sizeof(csize));
			
			tt.sptc_v = (NodeID*)memalign(64, (csize + 1) * sizeof(NodeID));
			tt.sptc_d = (EdgeWeight*)memalign(64, (csize + 1 ) * sizeof(EdgeWeight));
			
			total_children += (csize + 1);
			
			tt.sptc_v[0] = cid;
			tt.sptc_d[0] = csize;
			for (NodeID i = 0; i < csize; ++i) {
				ifs.read((char*)&cid, sizeof(cid));
				ifs.read((char*)&cd, sizeof(cd));
				tt.sptc_v[i + 1] = cid;
				tt.sptc_d[i + 1] = cd;			
			}
		}
		
		ifs.read((char*)&isize, sizeof(isize));
		r_numOfTokens = isize;		
		r_tokenindex_p = (token_t*)memalign(64, r_numOfTokens * sizeof(token_t));
		
		for (NodeID v = 0; v < r_numOfTokens; ++v) {
			token_t& tt = r_tokenindex_p[v];
		 
			ifs.read((char*)&cid, sizeof(cid));
			ifs.read((char*)&csize, sizeof(csize));
			
			tt.sptc_v = (NodeID*)memalign(64, (csize + 1) * sizeof(NodeID));
			tt.sptc_d = (EdgeWeight*)memalign(64, (csize + 1 ) * sizeof(EdgeWeight));
			
			r_total_children += (csize + 1);
			
			tt.sptc_v[0] = cid;
			tt.sptc_d[0] = csize;
			for (NodeID i = 0; i < csize; ++i) {
				ifs.read((char*)&cid, sizeof(cid));
				ifs.read((char*)&cd, sizeof(cd));
				tt.sptc_v[i + 1] = cid;
				tt.sptc_d[i + 1] = cd;			
			}
		}
		cout << "finish loading" << endl;
		ifs.close();
	}
	
	 
	void print_stat() {
		cout << "Total Token #: " << numOfTokens << endl;
		cout << "Average Children (Super) Token #: " << (double)total_children/(double)numOfTokens << endl;
		//cout << "Maximum Label Size: " << max_size() << endl;
	}
	
	void print_stat_d() {
		cout << "Total Token #: " << numOfTokens << endl;
		cout << "Total r_Token #: " << r_numOfTokens << endl;
		cout << "Average Children (Super) Token #: " << (double)total_children/(double)numOfTokens << endl;
		cout << "Average Children (Super) Token #: " << (double)r_total_children/(double)r_numOfTokens << endl;
	//	cout << "Maximum Label Size: " << max_size() << endl;
	} 
	
	EdgeWeight query_p(NodeID s, NodeID t, long ts, vector<NodeID>& dis_vec, vector<long>& ts_vec, vector<NodeID>& que, vector<EdgeWeight>& que_d) {
		if(s==t) return 0;
		
		EdgeWeight distance = INF_WEIGHT;
		
		NodeID anchor_s = anchor_p[s];
		NodeID anchor_t = anchor_p[t];

		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
		
		que_d[que_h] = 0;
		que[que_h++] = anchor_s;
		que_t1 = que_h;
		
		if(anchor_s < numOfVertices){
			if(ts_vec[anchor_s] != ts){
				ts_vec[anchor_s] = ts;
				dis_vec[anchor_s] = 0;
			}			
		}
		else{
			for (; que_t0 < que_h;) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID tid = que[que_i];
					EdgeWeight tdis = que_d[que_i];
					
					const token_t& token_v = tokenindex_p[tid - numOfVertices];
					
					_mm_prefetch(&token_v.sptc_v[0], _MM_HINT_T0);
					_mm_prefetch(&token_v.sptc_d[0], _MM_HINT_T0);
					
					NodeID r = token_v.sptc_v[0];
					EdgeWeight csize = token_v.sptc_d[0];
					
					// hashing, can be replaced by 1024 linear probing for efficiency.
					if(ts_vec[r] != ts){
						ts_vec[r] = ts;
						dis_vec[r] = tdis;
					}
					
					for (EdgeWeight i = 0; i < csize; ++i){
						NodeID w = token_v.sptc_v[i+1];
						EdgeWeight w_d = token_v.sptc_d[i+1] + tdis;
						if( w < numOfVertices){// hashing, can be replaced by 1024 linear probing for efficiency.
							if(ts_vec[w] != ts){
								ts_vec[w] = ts;
								dis_vec[w] = w_d;
							}
						}else{
							que_d[que_h] = w_d;
							que[que_h++] = w;
						}
					}
				}
				que_t0 = que_t1;
				que_t1 = que_h;
			}
		}
		
		
		que_t0 = 0, que_t1 = 0, que_h = 0;
		que_d[que_h] = 0;
		que[que_h++] = anchor_t;
		
		if(anchor_t < numOfVertices){
			if(ts_vec[anchor_t] == ts){
				EdgeWeight current_dis = dis_vec[anchor_t] + 0;
				if(current_dis < distance)
					distance = current_dis;
			}
		}else{
			que_t1 = que_h;
			for (; que_t0 < que_h;) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID tid = que[que_i];
					EdgeWeight tdis = que_d[que_i];
					
					const token_t& token_v = tokenindex_p[tid - numOfVertices];
					
					_mm_prefetch(&token_v.sptc_v[0], _MM_HINT_T0);
					_mm_prefetch(&token_v.sptc_d[0], _MM_HINT_T0);
					
					NodeID r = token_v.sptc_v[0];
					EdgeWeight csize = token_v.sptc_d[0];
					
					// hashing, can be replaced by 1024 linear probing for efficiency.
					if(ts_vec[r] == ts){
						EdgeWeight current_dis = dis_vec[r] + tdis;
						if(current_dis < distance)
							distance = current_dis;
					}
					
					for (EdgeWeight i = 0; i < csize; ++i){
						NodeID w = token_v.sptc_v[i+1];
						EdgeWeight w_d = token_v.sptc_d[i+1] + tdis;
						if( w < numOfVertices){
						// hashing, can be replaced by 1024 linear probing for efficiency.
							if(ts_vec[w] == ts){
								EdgeWeight current_dis = dis_vec[w] + w_d;
								if(current_dis < distance)
									distance = current_dis;
							}
						}else{							
							que_d[que_h] = w_d;
							que[que_h++] = w;
						}
					}
				}
				que_t0 = que_t1;
				que_t1 = que_h;
			}
		}		
		return distance;
	}
	
	EdgeWeight query_p_d(NodeID s, NodeID t, long ts, vector<NodeID>& dis_vec, vector<long>& ts_vec, vector<NodeID>& que, vector<EdgeWeight>& que_d) {
		if(s==t) return 0;
		
		EdgeWeight distance = INF_WEIGHT;
		
		NodeID anchor_s = anchor_p[s];
		NodeID anchor_t = r_anchor_p[t];

		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
		
		que_d[que_h] = 0;
		que[que_h++] = anchor_s;
		que_t1 = que_h;
		
		if(anchor_s < numOfVertices){
			if(ts_vec[anchor_s] != ts){
				ts_vec[anchor_s] = ts;
				dis_vec[anchor_s] = 0;
			}			
		}
		else{
			for (; que_t0 < que_h;) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID tid = que[que_i];
					EdgeWeight tdis = que_d[que_i];
					
					const token_t& token_v = tokenindex_p[tid - numOfVertices];
					
					_mm_prefetch(&token_v.sptc_v[0], _MM_HINT_T0);
					_mm_prefetch(&token_v.sptc_d[0], _MM_HINT_T0);
					
					NodeID r = token_v.sptc_v[0];
					EdgeWeight csize = token_v.sptc_d[0];
					
					// hashing, can be replaced by 1024 linear probing for efficiency.
					if(ts_vec[r] != ts){
						ts_vec[r] = ts;
						dis_vec[r] = tdis;
					}
					
					for (EdgeWeight i = 0; i < csize; ++i){
						NodeID w = token_v.sptc_v[i+1];
						EdgeWeight w_d = token_v.sptc_d[i+1] + tdis;
						if( w < numOfVertices){// hashing, can be replaced by 1024 linear probing for efficiency.
							if(ts_vec[w] != ts){
								ts_vec[w] = ts;
								dis_vec[w] = w_d;
							}
						}else{
							que_d[que_h] = w_d;
							que[que_h++] = w;
						}
					}
				}
				que_t0 = que_t1;
				que_t1 = que_h;
			} 
		}
		
		
		que_t0 = 0, que_t1 = 0, que_h = 0;
		que_d[que_h] = 0;
		que[que_h++] = anchor_t;
		
		if(anchor_t < numOfVertices){
			if(ts_vec[anchor_t] == ts){
				EdgeWeight current_dis = dis_vec[anchor_t] + 0;
				if(current_dis < distance)
					distance = current_dis;
			}
		}else{
			que_t1 = que_h;
			for (; que_t0 < que_h;) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID tid = que[que_i];
					EdgeWeight tdis = que_d[que_i];
					
					const token_t& token_v = r_tokenindex_p[tid - numOfVertices];
					
					_mm_prefetch(&token_v.sptc_v[0], _MM_HINT_T0);
					_mm_prefetch(&token_v.sptc_d[0], _MM_HINT_T0);
					
					NodeID r = token_v.sptc_v[0];
					EdgeWeight csize = token_v.sptc_d[0];
					
					// hashing, can be replaced by 1024 linear probing for efficiency.
					if(ts_vec[r] == ts){
						EdgeWeight current_dis = dis_vec[r] + tdis;
						if(current_dis < distance)
							distance = current_dis;
					}
					
					for (EdgeWeight i = 0; i < csize; ++i){
						NodeID w = token_v.sptc_v[i+1];
						EdgeWeight w_d = token_v.sptc_d[i+1] + tdis;
						if( w < numOfVertices){
						// hashing, can be replaced by 1024 linear probing for efficiency.
							if(ts_vec[w] == ts){
								EdgeWeight current_dis = dis_vec[w] + w_d;
								if(current_dis < distance)
									distance = current_dis;
							}
						}else{							
							que_d[que_h] = w_d;
							que[que_h++] = w;
						}
					}
				}
				que_t0 = que_t1;
				que_t1 = que_h;
			}
		}
		
		
		
		return distance;
	}
	
	
	void save_two_level_labels(const char* save_filename) {
		ofstream ofs(save_filename, ios::binary | ios::out);

		ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
		for (NodeID v = 0; v < numOfVertices; ++v) {
			ofs.write((const char*)&anchor_p[v], sizeof(anchor_p[v]));
		//	ofs.write((const char*)&index_[v].spt_d[i], sizeof(index_[v].spt_d[i]));
		}
		
		// Store supertokens 
		for (NodeID v = 0; v < numOfVertices; ++v) {
			token_t& supertoken_v = supertokenindex_p[v];
			
			NodeID isize = supertoken_v.sptc_v[0];			
			ofs.write((const char*)&isize, sizeof(isize));
			
			for(NodeID i = 0; i < isize; ++i){
				NodeID tid = supertoken_v.sptc_v[i + 1];
				EdgeWeight ew = supertoken_v.sptc_d[i + 1];	
				
				ofs.write((const char*)&tid, sizeof(tid));
				ofs.write((const char*)&ew, sizeof(ew));
				
			}			
		//	ofs.write((const char*)&index_[v].spt_d[i], sizeof(index_[v].spt_d[i]));
		}
		
		// Store normal tokens
		ofs.write((const char*)&numOfTokens, sizeof(numOfTokens));
		for (NodeID t = 0; t < numOfTokens; ++t) {
			token_t& tt = tokenindex_p[t];
			NodeID sid = tt.sptc_v[0];
			EdgeWeight ssize = tt.sptc_d[0];
			EdgeWeight fsize = supertokenindex_p[sid].sptc_d[0];
			ofs.write((const char*)&sid, sizeof(sid));			
			ofs.write((const char*)&ssize, sizeof(ssize));	
			if(ssize == 0) continue;
			//ofs.write((const char*)&fsize, sizeof(fsize));
			//if(t < 10) 
			//	cout << sid << "vs" << fsize << "vs" << ssize << endl;
			
			for(NodeID c = 0; c < fsize; ++c){		
				//char a = tt.sptc_fbv[c];
				//ofs.write((const char*)&a, sizeof(a));
				ofs.write((const char*)&tt.sptc_fbv[c], sizeof(tt.sptc_fbv[c]));
			//	if(t < 10){
		//			bitset<8> s(tt.sptc_fbv[c]);
			//		cout << s;
			//	}
			} 
			//if(t < 10)
			//	cout << endl;
					
			for(NodeID c = 0; c < ssize; ++c){			
				//char a = tt.sptc_sbv[c];
				//ofs.write((const char*)&a, sizeof(a));	
				ofs.write((const char*)&tt.sptc_sbv[c], sizeof(tt.sptc_sbv[c]));
			//	if(t < 10){
			//		bitset<8> s(tt.sptc_sbv[c]);
			//		cout << s;
			//	}
			} 
			//if(t < 10)
			//	cout << endl;
		}
		
		ofs.close();
	}
	 
	void load_two_level_labels(const char* load_filename) {
		total_children = 0;
		
		tokenindex_p = NULL;
		anchor_p = NULL;
		supertokenindex_p = NULL;

		ifstream ifs(load_filename);
		NodeID isize = 0;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;

		anchor_p = (NodeID*)memalign(64, numOfVertices * sizeof(NodeID));
		NodeID anchor_id;
		for (NodeID v = 0; v < numOfVertices; ++v) {
				ifs.read((char*)&anchor_id, sizeof(anchor_id));
				anchor_p[v] = anchor_id;
		}
		
		//load supertokens
		NodeID cid;
		EdgeWeight cd;
		supertokenindex_p = (token_t*)memalign(64, numOfVertices * sizeof(token_t));
		for (NodeID v = 0; v < numOfVertices; ++v) {
				
				token_t& supertoken_v = supertokenindex_p[v];					
			
				NodeID csize;
				ifs.read((char*)&csize, sizeof(csize));
				
				supertoken_v.sptc_v = (NodeID*)memalign(64, (csize + 1) * sizeof(NodeID));
				supertoken_v.sptc_d = (EdgeWeight*)memalign(64, (csize + 1 ) * sizeof(EdgeWeight));
				
				supertoken_v.sptc_v[0] = csize;
				
				NodeID intsize = ceil((double)ceil((double)csize / (double)8) / (double)8);				
				supertoken_v.sptc_d[0] = intsize;
				
				total_children += csize;
				for(EdgeWeight i = 0; i < csize; ++i){
					ifs.read((char*)&cid, sizeof(cid));
					ifs.read((char*)&cd, sizeof(cd));
					supertoken_v.sptc_v[i + 1] = cid;
					supertoken_v.sptc_d[i + 1] = cd;	
				}				
		}
		cout << "loaded supertokens" << endl;
		cout << "Average Children Super Token #: " << (double)total_children/(double)numOfVertices << endl;
		
		ifs.read((char*)&isize, sizeof(isize));
		numOfTokens = isize;			
		NodeID sid;
		EdgeWeight ssize;
		EdgeWeight fsize;
		cout<< numOfTokens << " tokens in total." << endl;
		
		tokenindex_p = (token_t*)memalign(64, numOfTokens * sizeof(token_t));	
		for (NodeID v = 0; v < numOfTokens; ++v) {
			token_t& tt = tokenindex_p[v];
		
			ifs.read((char*)&sid, sizeof(sid));
			ifs.read((char*)&ssize, sizeof(ssize));
			
			tt.sptc_v = (NodeID*)memalign(64, 1 * sizeof(NodeID));
			tt.sptc_d = (EdgeWeight*)memalign(64, 1 * sizeof(EdgeWeight));
			
			tt.sptc_v[0] = sid;
			tt.sptc_d[0] = ssize;			
			fsize = supertokenindex_p[sid].sptc_d[0];			
			
			if(ssize == 0) continue;
			//if(v < 10) 
			//	cout << sid << "vs" << fsize << "vs" << ssize << endl;
			
			tt.sptc_fbv = (unsigned char*)memalign(64, fsize * sizeof(unsigned char));		
		//	unsigned char fb;
			char fb;
			for (NodeID i = 0; i < fsize; ++i) {
				ifs.read((char*)&(tt.sptc_fbv[i]), sizeof(tt.sptc_fbv[i]));
				//ifs.read((char*)&fb, sizeof(fb));
			//	if(v < 10){
			//		bitset<8> s(tt.sptc_fbv[i]);
			//		cout << s;
			//	}
			}
			//if(v < 10)
			//	cout << endl;
			
			  
			tt.sptc_sbv = (unsigned char*)memalign(64, ssize * sizeof(unsigned char));		
			//unsigned char sb;
			char sb;
			for (NodeID i = 0; i < ssize; ++i) {
				ifs.read((char*)&(tt.sptc_sbv[i]), sizeof(tt.sptc_sbv[i]));
				//ifs.read((char*)&sb, sizeof(sb));
			//	if(v < 10){
			//		bitset<8> s(tt.sptc_sbv[i]);
			//		cout << s;
				//}
			}
			//if(v < 10)
			//	cout << endl;
			//
		}
		cout << "loaded standard tokens" << endl;
		ifs.close();
	} 
	
	void save_two_level_labels_path(const char* save_filename) {
		ofstream ofs(save_filename, ios::binary | ios::out);

		ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
		for (NodeID v = 0; v < numOfVertices; ++v) {
			ofs.write((const char*)&anchor_p[v], sizeof(anchor_p[v]));
		//	ofs.write((const char*)&index_[v].spt_d[i], sizeof(index_[v].spt_d[i]));
		}
		
		// Store supertokens 
		for (NodeID v = 0; v < numOfVertices; ++v) {
			token_t& supertoken_v = supertokenindex_p[v];
			
			NodeID isize = supertoken_v.sptc_v[0];			
			ofs.write((const char*)&isize, sizeof(isize));
			
			for(NodeID i = 0; i < isize; ++i){
				NodeID tid = supertoken_v.sptc_v[i + 1];
				EdgeWeight ew = supertoken_v.sptc_d[i + 1];	
				NodeID pid = supertoken_v.sptc_pathv[i + 1];
				
				ofs.write((const char*)&tid, sizeof(tid));
				ofs.write((const char*)&ew, sizeof(ew));
				ofs.write((const char*)&pid, sizeof(pid));
				
			}			
		//	ofs.write((const char*)&index_[v].spt_d[i], sizeof(index_[v].spt_d[i]));
		}
		
		// Store normal tokens
		ofs.write((const char*)&numOfTokens, sizeof(numOfTokens));
		for (NodeID t = 0; t < numOfTokens; ++t) {
			token_t& tt = tokenindex_p[t];
			NodeID sid = tt.sptc_v[0];
			EdgeWeight ssize = tt.sptc_d[0];
			EdgeWeight fsize = supertokenindex_p[sid].sptc_d[0];
			ofs.write((const char*)&sid, sizeof(sid));			
			ofs.write((const char*)&ssize, sizeof(ssize));	
			if(ssize == 0) continue;
			//ofs.write((const char*)&fsize, sizeof(fsize));
			//if(t < 10) 
			//	cout << sid << "vs" << fsize << "vs" << ssize << endl;
			
			for(NodeID c = 0; c < fsize; ++c){		
				//char a = tt.sptc_fbv[c];
				//ofs.write((const char*)&a, sizeof(a));
				ofs.write((const char*)&tt.sptc_fbv[c], sizeof(tt.sptc_fbv[c]));
			//	if(t < 10){
		//			bitset<8> s(tt.sptc_fbv[c]);
			//		cout << s;
			//	}
			} 
			//if(t < 10)
			//	cout << endl;
					
			for(NodeID c = 0; c < ssize; ++c){			
				//char a = tt.sptc_sbv[c];
				//ofs.write((const char*)&a, sizeof(a));	
				ofs.write((const char*)&tt.sptc_sbv[c], sizeof(tt.sptc_sbv[c]));
			//	if(t < 10){
			//		bitset<8> s(tt.sptc_sbv[c]);
			//		cout << s;
			//	}
			} 
			//if(t < 10)
			//	cout << endl;
		}
		
		ofs.close();
	}
	
	
	void load_two_level_labels_path(const char* load_filename) {
		total_children = 0;
		
		tokenindex_p = NULL;
		anchor_p = NULL;
		supertokenindex_p = NULL;

		ifstream ifs(load_filename);
		NodeID isize = 0;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;

		anchor_p = (NodeID*)memalign(64, numOfVertices * sizeof(NodeID));
		NodeID anchor_id;
		for (NodeID v = 0; v < numOfVertices; ++v) {
				ifs.read((char*)&anchor_id, sizeof(anchor_id));
				anchor_p[v] = anchor_id;
		}
		
		//load supertokens
		NodeID cid;
		EdgeWeight cd;
		supertokenindex_p = (token_t*)memalign(64, numOfVertices * sizeof(token_t));
		for (NodeID v = 0; v < numOfVertices; ++v) {
				
				token_t& supertoken_v = supertokenindex_p[v];					
			
				NodeID csize;
				ifs.read((char*)&csize, sizeof(csize));
				
				supertoken_v.sptc_v = (NodeID*)memalign(64, (csize + 1) * sizeof(NodeID));
				supertoken_v.sptc_d = (EdgeWeight*)memalign(64, (csize + 1 ) * sizeof(EdgeWeight));
				supertoken_v.sptc_pathv = (EdgeWeight*)memalign(64, (csize + 1 ) * sizeof(EdgeWeight));
				
				supertoken_v.sptc_v[0] = csize;
				
				NodeID intsize = ceil((double)ceil((double)csize / (double)8) / (double)8);				
				supertoken_v.sptc_d[0] = intsize;
				
				supertoken_v.sptc_pathv[0] = numOfVertices;
				
				total_children += csize;
				for(EdgeWeight i = 0; i < csize; ++i){
					ifs.read((char*)&cid, sizeof(cid));
					ifs.read((char*)&cd, sizeof(cd));
					supertoken_v.sptc_v[i + 1] = cid;
					supertoken_v.sptc_d[i + 1] = cd;	
					ifs.read((char*)&cid, sizeof(cid));
					supertoken_v.sptc_pathv[i + 1] = cid;
				}				
		}
		cout << "loaded supertokens" << endl;
		cout << "Average Children Super Token #: " << (double)total_children/(double)numOfVertices << endl;
		
		ifs.read((char*)&isize, sizeof(isize));
		numOfTokens = isize;			
		NodeID sid;
		EdgeWeight ssize;
		EdgeWeight fsize;
		cout<< numOfTokens << " tokens in total." << endl;
		
		tokenindex_p = (token_t*)memalign(64, numOfTokens * sizeof(token_t));	
		for (NodeID v = 0; v < numOfTokens; ++v) {
			token_t& tt = tokenindex_p[v];
		
			ifs.read((char*)&sid, sizeof(sid));
			ifs.read((char*)&ssize, sizeof(ssize));
			
			tt.sptc_v = (NodeID*)memalign(64, 1 * sizeof(NodeID));
			tt.sptc_d = (EdgeWeight*)memalign(64, 1 * sizeof(EdgeWeight));
			
			tt.sptc_v[0] = sid;
			tt.sptc_d[0] = ssize;			
			fsize = supertokenindex_p[sid].sptc_d[0];			
			
			if(ssize == 0) continue;
			//if(v < 10) 
			//	cout << sid << "vs" << fsize << "vs" << ssize << endl;
			
			tt.sptc_fbv = (unsigned char*)memalign(64, fsize * sizeof(unsigned char));		
		//	unsigned char fb;
			char fb;
			for (NodeID i = 0; i < fsize; ++i) {
				ifs.read((char*)&(tt.sptc_fbv[i]), sizeof(tt.sptc_fbv[i]));
				//ifs.read((char*)&fb, sizeof(fb));
			//	if(v < 10){
			//		bitset<8> s(tt.sptc_fbv[i]);
			//		cout << s;
			//	}
			}
			//if(v < 10)
			//	cout << endl;
			
			  
			tt.sptc_sbv = (unsigned char*)memalign(64, ssize * sizeof(unsigned char));		
			//unsigned char sb;
			char sb;
			for (NodeID i = 0; i < ssize; ++i) {
				ifs.read((char*)&(tt.sptc_sbv[i]), sizeof(tt.sptc_sbv[i]));
				//ifs.read((char*)&sb, sizeof(sb));
			//	if(v < 10){
			//		bitset<8> s(tt.sptc_sbv[i]);
			//		cout << s;
				//}
			}
			//if(v < 10)
			//	cout << endl;
			//
		}
		cout << "loaded standard tokens" << endl;
		ifs.close();
	} 
	
	void save_two_level_labels_d(const char* save_filename) {
		ofstream ofs(save_filename, ios::binary | ios::out);
		//cout << "1" << endl;
		ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
		for (NodeID v = 0; v < numOfVertices; ++v) {
			ofs.write((const char*)&anchor_p[v], sizeof(anchor_p[v]));
		//	ofs.write((const char*)&index_[v].spt_d[i], sizeof(index_[v].spt_d[i]));
		}
		for (NodeID v = 0; v < numOfVertices; ++v) {
			ofs.write((const char*)&r_anchor_p[v], sizeof(r_anchor_p[v]));
		//	ofs.write((const char*)&index_[v].spt_d[i], sizeof(index_[v].spt_d[i]));
		}
				
		
		// Store supertokens 
	//	cout << "2" << endl;
		for (NodeID v = 0; v < numOfVertices; ++v) {
			token_t& supertoken_v = supertokenindex_p[v];
			
			NodeID isize = supertoken_v.sptc_v[0];			
			ofs.write((const char*)&isize, sizeof(isize));
			
			for(NodeID i = 0; i < isize; ++i){
				NodeID tid = supertoken_v.sptc_v[i + 1];
				EdgeWeight ew = supertoken_v.sptc_d[i + 1];	
				
				ofs.write((const char*)&tid, sizeof(tid));
				ofs.write((const char*)&ew, sizeof(ew));
				
			}			
		}
		
		for (NodeID v = 0; v < numOfVertices; ++v) {
			token_t& supertoken_v = r_supertokenindex_p[v];
			
			NodeID isize = supertoken_v.sptc_v[0];			
			ofs.write((const char*)&isize, sizeof(isize));
			
			for(NodeID i = 0; i < isize; ++i){
				NodeID tid = supertoken_v.sptc_v[i + 1];
				EdgeWeight ew = supertoken_v.sptc_d[i + 1];	
				
				ofs.write((const char*)&tid, sizeof(tid));
				ofs.write((const char*)&ew, sizeof(ew));
				
			}			
		}
		
		 
		// Store normal tokens
		//cout << "3" << endl;
		ofs.write((const char*)&numOfTokens, sizeof(numOfTokens));
		for (NodeID t = 0; t < numOfTokens; ++t) {
			
		//	cout << "31:" << t << endl;
			token_t& tt = tokenindex_p[t]; 
			NodeID sid = tt.sptc_v[0];
			EdgeWeight ssize = tt.sptc_d[0];
			EdgeWeight fsize = supertokenindex_p[sid].sptc_d[0];
			ofs.write((const char*)&sid, sizeof(sid));			
			ofs.write((const char*)&ssize, sizeof(ssize));	
			
			
		//	cout << "32:" << t << endl;
			if(ssize == 0) continue;
			//ofs.write((const char*)&fsize, sizeof(fsize));
			//if(t < 10) 
			//	cout << sid << "vs" << fsize << "vs" << ssize << endl;
			
			
		//	cout << "33:" << t << endl;
			for(NodeID c = 0; c < fsize; ++c){		
				//char a = tt.sptc_fbv[c];
				//ofs.write((const char*)&a, sizeof(a));
				ofs.write((const char*)&tt.sptc_fbv[c], sizeof(tt.sptc_fbv[c]));
			//	if(t < 10){
		//			bitset<8> s(tt.sptc_fbv[c]);
			//		cout << s;
			//	}
			} 
			//if(t < 10) 
			//	cout << endl;
					
		//	cout << "34:" << t << endl;
			for(NodeID c = 0; c < ssize; ++c){			
				//char a = tt.sptc_sbv[c];
				//ofs.write((const char*)&a, sizeof(a));	
				ofs.write((const char*)&tt.sptc_sbv[c], sizeof(tt.sptc_sbv[c]));
			//	if(t < 10){
			//		bitset<8> s(tt.sptc_sbv[c]);
			//		cout << s;
			//	}
			} 
			//if(t < 10)
			//	cout << endl;
		} 
		
		//cout << "4" << endl;
		ofs.write((const char*)&r_numOfTokens, sizeof(r_numOfTokens));
		for (NodeID t = 0; t < r_numOfTokens; ++t) {
					
			//cout << "41:" << t << endl;
			token_t& tt = r_tokenindex_p[t];
			NodeID sid = tt.sptc_v[0];
			EdgeWeight ssize = tt.sptc_d[0];
			EdgeWeight fsize = r_supertokenindex_p[sid].sptc_d[0];
			ofs.write((const char*)&sid, sizeof(sid));			
			ofs.write((const char*)&ssize, sizeof(ssize));	
			if(ssize == 0) continue; 
			//ofs.write((const char*)&fsize, sizeof(fsize));
			//if(t < 10) 
			//	cout << sid << "vs" << fsize << "vs" << ssize << endl;
			
			//cout << "42:" << t << "," << fsize <<  endl;
			for(NodeID c = 0; c < fsize; ++c){		
				ofs.write((const char*)&tt.sptc_fbv[c], sizeof(tt.sptc_fbv[c]));
			} 

					
			//cout << "43:" << t << "," << ssize << endl;
			for(NodeID c = 0; c < ssize; ++c){		
				ofs.write((const char*)&tt.sptc_sbv[c], sizeof(tt.sptc_sbv[c]));	
			} 
		}	
		
		ofs.close();
	}
	 
	void load_two_level_labels_d(const char* load_filename) {
		total_children = 0;		
		tokenindex_p = NULL;
		anchor_p = NULL;
		supertokenindex_p = NULL;
		
		r_total_children = 0;		
		r_tokenindex_p = NULL;
		r_anchor_p = NULL;
		r_supertokenindex_p = NULL;

		ifstream ifs(load_filename);
		NodeID isize = 0;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;

		anchor_p = (NodeID*)memalign(64, numOfVertices * sizeof(NodeID));
		r_anchor_p = (NodeID*)memalign(64, numOfVertices * sizeof(NodeID));
		NodeID anchor_id;
		for (NodeID v = 0; v < numOfVertices; ++v) {
				ifs.read((char*)&anchor_id, sizeof(anchor_id));
				anchor_p[v] = anchor_id;
		}
		for (NodeID v = 0; v < numOfVertices; ++v) {
				ifs.read((char*)&anchor_id, sizeof(anchor_id));
				r_anchor_p[v] = anchor_id;
		}
		
		//load supertokens
		NodeID cid;
		EdgeWeight cd;
		supertokenindex_p = (token_t*)memalign(64, numOfVertices * sizeof(token_t));
		for (NodeID v = 0; v < numOfVertices; ++v) {
				
				token_t& supertoken_v = supertokenindex_p[v];					
			
				NodeID csize;
				ifs.read((char*)&csize, sizeof(csize));
				
				supertoken_v.sptc_v = (NodeID*)memalign(64, (csize + 1) * sizeof(NodeID));
				supertoken_v.sptc_d = (EdgeWeight*)memalign(64, (csize + 1 ) * sizeof(EdgeWeight));
				
				supertoken_v.sptc_v[0] = csize;
				
				NodeID intsize = ceil((double)ceil((double)csize / (double)8) / (double)8);				
				supertoken_v.sptc_d[0] = intsize;
				
				total_children += csize;
				for(EdgeWeight i = 0; i < csize; ++i){
					ifs.read((char*)&cid, sizeof(cid));
					ifs.read((char*)&cd, sizeof(cd));
					supertoken_v.sptc_v[i + 1] = cid;
					supertoken_v.sptc_d[i + 1] = cd;	
				}				
		}	
		
		r_supertokenindex_p = (token_t*)memalign(64, numOfVertices * sizeof(token_t));
		for (NodeID v = 0; v < numOfVertices; ++v) {
				
				token_t& supertoken_v = r_supertokenindex_p[v];					
			
				NodeID csize;
				ifs.read((char*)&csize, sizeof(csize));
				
				supertoken_v.sptc_v = (NodeID*)memalign(64, (csize + 1) * sizeof(NodeID));
				supertoken_v.sptc_d = (EdgeWeight*)memalign(64, (csize + 1 ) * sizeof(EdgeWeight));
				
				supertoken_v.sptc_v[0] = csize;
				
				NodeID intsize = ceil((double)ceil((double)csize / (double)8) / (double)8);				
				supertoken_v.sptc_d[0] = intsize;
				
				r_total_children += csize;
				for(EdgeWeight i = 0; i < csize; ++i){
					ifs.read((char*)&cid, sizeof(cid));
					ifs.read((char*)&cd, sizeof(cd));
					supertoken_v.sptc_v[i + 1] = cid;
					supertoken_v.sptc_d[i + 1] = cd;	
				}				
		}
		cout << "loaded supertokens" << endl;
		cout << "Average Children Super Token #: " << (double)total_children/(double)numOfVertices << endl;
		cout << "Average Children Super Token #: " << (double)r_total_children/(double)numOfVertices << endl;
	
	
		ifs.read((char*)&isize, sizeof(isize));
		numOfTokens = isize;			
		NodeID sid;
		EdgeWeight ssize;
		EdgeWeight fsize;
		cout<< numOfTokens << " tokens in total." << endl;		
		tokenindex_p = (token_t*)memalign(64, numOfTokens * sizeof(token_t));	
		for (NodeID v = 0; v < numOfTokens; ++v) {
			token_t& tt = tokenindex_p[v];
		
			ifs.read((char*)&sid, sizeof(sid));
			ifs.read((char*)&ssize, sizeof(ssize));
			
			tt.sptc_v = (NodeID*)memalign(64, 1 * sizeof(NodeID));
			tt.sptc_d = (EdgeWeight*)memalign(64, 1 * sizeof(EdgeWeight));
			
			tt.sptc_v[0] = sid;
			tt.sptc_d[0] = ssize;			
			fsize = supertokenindex_p[sid].sptc_d[0];			
			
			if(ssize == 0) continue;
			//if(v < 10) 
			//	cout << sid << "vs" << fsize << "vs" << ssize << endl;
			
			tt.sptc_fbv = (unsigned char*)memalign(64, fsize * sizeof(unsigned char));		
		//	unsigned char fb;
			char fb;
			for (NodeID i = 0; i < fsize; ++i) {
				ifs.read((char*)&(tt.sptc_fbv[i]), sizeof(tt.sptc_fbv[i]));
				//ifs.read((char*)&fb, sizeof(fb));
			//	if(v < 10){
			//		bitset<8> s(tt.sptc_fbv[i]);
			//		cout << s;
			//	}
			}
			//if(v < 10)
			//	cout << endl;
			
			  
			tt.sptc_sbv = (unsigned char*)memalign(64, ssize * sizeof(unsigned char));		
			//unsigned char sb;
			char sb;
			for (NodeID i = 0; i < ssize; ++i) {
				ifs.read((char*)&(tt.sptc_sbv[i]), sizeof(tt.sptc_sbv[i]));
				//ifs.read((char*)&sb, sizeof(sb));
			//	if(v < 10){
			//		bitset<8> s(tt.sptc_sbv[i]);
			//		cout << s;
				//}
			}
			//if(v < 10)
			//	cout << endl;
			//
		}
		
		ifs.read((char*)&isize, sizeof(isize));
		r_numOfTokens = isize;			
		cout<< r_numOfTokens << " tokens in total." << endl;		
		r_tokenindex_p = (token_t*)memalign(64, r_numOfTokens * sizeof(token_t));	
		for (NodeID v = 0; v < r_numOfTokens; ++v) {
			token_t& tt = r_tokenindex_p[v];
		
			ifs.read((char*)&sid, sizeof(sid));
			ifs.read((char*)&ssize, sizeof(ssize));
			
			tt.sptc_v = (NodeID*)memalign(64, 1 * sizeof(NodeID));
			tt.sptc_d = (EdgeWeight*)memalign(64, 1 * sizeof(EdgeWeight));
			
			tt.sptc_v[0] = sid;
			tt.sptc_d[0] = ssize;			
			fsize = r_supertokenindex_p[sid].sptc_d[0];			
			
			if(ssize == 0) continue;
			//if(v < 10) 
			//	cout << sid << "vs" << fsize << "vs" << ssize << endl;
			
			tt.sptc_fbv = (unsigned char*)memalign(64, fsize * sizeof(unsigned char));		
		//	unsigned char fb;
			char fb;
			for (NodeID i = 0; i < fsize; ++i) {
				ifs.read((char*)&(tt.sptc_fbv[i]), sizeof(tt.sptc_fbv[i]));
				//ifs.read((char*)&fb, sizeof(fb));
			//	if(v < 10){
			//		bitset<8> s(tt.sptc_fbv[i]);
			//		cout << s;
			//	}
			}
			//if(v < 10)
			//	cout << endl;
			
			  
			tt.sptc_sbv = (unsigned char*)memalign(64, ssize * sizeof(unsigned char));		
			//unsigned char sb;
			char sb;
			for (NodeID i = 0; i < ssize; ++i) {
				ifs.read((char*)&(tt.sptc_sbv[i]), sizeof(tt.sptc_sbv[i]));
				//ifs.read((char*)&sb, sizeof(sb));
			//	if(v < 10){
			//		bitset<8> s(tt.sptc_sbv[i]);
			//		cout << s;
				//}
			}
			//if(v < 10)
			//	cout << endl;
			//
		}
		 
		cout << "loaded standard tokens" << endl;
		ifs.close();
	}
	 
	 

	EdgeWeight query_p_two_level(NodeID s, NodeID t, long ts, vector<NodeID>& dis_vec, vector<long>& ts_vec, vector<NodeID>& que, vector<EdgeWeight>& que_d) {
		if(s==t) return 0;
		
		EdgeWeight distance = INF_WEIGHT;
		
		NodeID anchor_s = anchor_p[s];
		NodeID anchor_t = anchor_p[t];

		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
		
		que_d[que_h] = 0;
		que[que_h++] = anchor_s;
		que_t1 = que_h;
		
		if(anchor_s < numOfVertices){
			if(ts_vec[anchor_s] != ts){
				ts_vec[anchor_s] = ts;
				dis_vec[anchor_s] = 0;
			}			
		}
		else{
			for (; que_t0 < que_h;) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID tid = que[que_i];
					EdgeWeight tdis = que_d[que_i];
					
					const token_t& token_v = tokenindex_p[tid - numOfVertices];
					
					_mm_prefetch(&token_v.sptc_v[0], _MM_HINT_T0);
					_mm_prefetch(&token_v.sptc_d[0], _MM_HINT_T0);
					
					NodeID r = token_v.sptc_v[0];					
					EdgeWeight ssize = token_v.sptc_d[0];
					
					token_t& supertoken_r = supertokenindex_p[r];
					EdgeWeight fsize = supertoken_r.sptc_d[0];
					
					// hashing, can be replaced by 1024 linear probing for efficiency.
					if(ts_vec[r] != ts){
						ts_vec[r] = ts;
						dis_vec[r] = tdis;
					}
					
					EdgeWeight spos = 0;
					
					for(EdgeWeight i = 0; i < fsize; ++i){
						unsigned char fmask = token_v.sptc_fbv[i];						
						bitset<8> fbs(fmask);
						for(NodeID j = 0; j < 8; ++j){							
							if(fbs[ 7 - j]){
								unsigned char smask = token_v.sptc_sbv[spos++];
								bitset<8> sbs(smask);
								for(NodeID k = 0; k < 8; ++k){
									if(sbs[7 - k]){
										NodeID w = supertoken_r.sptc_v[ (i * 8 + j) * 8 + k  +  1];
										EdgeWeight w_d = supertoken_r.sptc_d[(i * 8 + j) * 8 + k  +  1] + tdis;											
										if( w < numOfVertices){// hashing, can be replaced by 1024 linear probing for efficiency.
											if(ts_vec[w] != ts){
												ts_vec[w] = ts;
												dis_vec[w] = w_d;
											}
										}else{
											que_d[que_h] = w_d;
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
		}
		
		
		que_t0 = 0, que_t1 = 0, que_h = 0;
		que_d[que_h] = 0;
		que[que_h++] = anchor_t;
		
		if(anchor_t < numOfVertices){
			if(ts_vec[anchor_t] == ts){
				EdgeWeight current_dis = dis_vec[anchor_t] + 0;
				if(current_dis < distance)
					distance = current_dis;
			}
		}else{
			que_t1 = que_h;
			for (; que_t0 < que_h;) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID tid = que[que_i];
					EdgeWeight tdis = que_d[que_i];
					
					const token_t& token_v = tokenindex_p[tid - numOfVertices];
					
					_mm_prefetch(&token_v.sptc_v[0], _MM_HINT_T0);
					_mm_prefetch(&token_v.sptc_d[0], _MM_HINT_T0);
					
					
					
					NodeID r = token_v.sptc_v[0];					
					EdgeWeight ssize = token_v.sptc_d[0];					
					token_t& supertoken_r = supertokenindex_p[r];
					EdgeWeight fsize = supertoken_r.sptc_d[0];
					
					
					// hashing, can be replaced by 1024 linear probing for efficiency.
					if(ts_vec[r] == ts){
						EdgeWeight current_dis = dis_vec[r] + tdis;
						if(current_dis < distance)
							distance = current_dis;
					}
					
					EdgeWeight spos = 0;
					
					for(EdgeWeight i = 0; i < fsize; ++i){
						unsigned char fmask = token_v.sptc_fbv[i];						
						bitset<8> fbs(fmask);
						for(NodeID j = 0; j < 8; ++j){							
							if(fbs[7 - j]){
								unsigned char smask = token_v.sptc_sbv[spos++];
								bitset<8> sbs(smask);
								for(NodeID k = 0; k < 8; ++k){
									if(sbs[7 - k]){
										NodeID w = supertoken_r.sptc_v[ (i * 8 + j) * 8 + k  +  1];
										EdgeWeight w_d = supertoken_r.sptc_d[(i * 8 + j) * 8 + k  +  1] + tdis;											
										if( w < numOfVertices){// hashing, can be replaced by 1024 linear probing for efficiency.
											if(ts_vec[w] == ts){
												EdgeWeight current_dis = dis_vec[w] + w_d;
												if(current_dis < distance)
													distance = current_dis;
											}
										}else{
											que_d[que_h] = w_d;
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
		}
		
		
		
		return distance;
	}
	 
	EdgeWeight query_p_two_level_d(NodeID s, NodeID t, long ts, vector<NodeID>& dis_vec, vector<long>& ts_vec, vector<NodeID>& que, vector<EdgeWeight>& que_d) {
		if(s==t) return 0;
		
		EdgeWeight distance = INF_WEIGHT;
		
		NodeID anchor_s = anchor_p[s];
		NodeID anchor_t = r_anchor_p[t];

		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
		
		que_d[que_h] = 0;
		que[que_h++] = anchor_s;
		que_t1 = que_h;
		
		if(anchor_s < numOfVertices){
			if(ts_vec[anchor_s] != ts){
				ts_vec[anchor_s] = ts;
				dis_vec[anchor_s] = 0;
			}			
		}
		else{
			for (; que_t0 < que_h;) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID tid = que[que_i];
					EdgeWeight tdis = que_d[que_i];
					
					const token_t& token_v = tokenindex_p[tid - numOfVertices];
					
					_mm_prefetch(&token_v.sptc_v[0], _MM_HINT_T0);
					_mm_prefetch(&token_v.sptc_d[0], _MM_HINT_T0);
					
					NodeID r = token_v.sptc_v[0];					
					EdgeWeight ssize = token_v.sptc_d[0];
					
					token_t& supertoken_r = supertokenindex_p[r];
					EdgeWeight fsize = supertoken_r.sptc_d[0];
					
					// hashing, can be replaced by 1024 linear probing for efficiency.
					if(ts_vec[r] != ts){
						ts_vec[r] = ts;
						dis_vec[r] = tdis;
					}
					
					EdgeWeight spos = 0;
					
					for(EdgeWeight i = 0; i < fsize; ++i){
						unsigned char fmask = token_v.sptc_fbv[i];						
						bitset<8> fbs(fmask);
						for(NodeID j = 0; j < 8; ++j){							
							if(fbs[ 7 - j]){
								unsigned char smask = token_v.sptc_sbv[spos++];
								bitset<8> sbs(smask);
								for(NodeID k = 0; k < 8; ++k){
									if(sbs[7 - k]){
										NodeID w = supertoken_r.sptc_v[ (i * 8 + j) * 8 + k  +  1];
										EdgeWeight w_d = supertoken_r.sptc_d[(i * 8 + j) * 8 + k  +  1] + tdis;											
										if( w < numOfVertices){// hashing, can be replaced by 1024 linear probing for efficiency.
											if(ts_vec[w] != ts){
												ts_vec[w] = ts;
												dis_vec[w] = w_d;
											}
										}else{
											que_d[que_h] = w_d;
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
		}
		
		
		que_t0 = 0, que_t1 = 0, que_h = 0;
		que_d[que_h] = 0;
		que[que_h++] = anchor_t;
		
		if(anchor_t < numOfVertices){
			if(ts_vec[anchor_t] == ts){
				EdgeWeight current_dis = dis_vec[anchor_t] + 0;
				if(current_dis < distance)
					distance = current_dis;
			}
		}else{
			que_t1 = que_h;
			for (; que_t0 < que_h;) {
				for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
					NodeID tid = que[que_i];
					EdgeWeight tdis = que_d[que_i];
					
					const token_t& token_v = r_tokenindex_p[tid - numOfVertices];
					
					_mm_prefetch(&token_v.sptc_v[0], _MM_HINT_T0);
					_mm_prefetch(&token_v.sptc_d[0], _MM_HINT_T0);
					
					
					
					NodeID r = token_v.sptc_v[0];					
					EdgeWeight ssize = token_v.sptc_d[0];					
					token_t& supertoken_r = r_supertokenindex_p[r];
					EdgeWeight fsize = supertoken_r.sptc_d[0];
					
					
					// hashing, can be replaced by 1024 linear probing for efficiency.
					if(ts_vec[r] == ts){
						EdgeWeight current_dis = dis_vec[r] + tdis;
						if(current_dis < distance)
							distance = current_dis;
					}
					
					EdgeWeight spos = 0;
					
					for(EdgeWeight i = 0; i < fsize; ++i){
						unsigned char fmask = token_v.sptc_fbv[i];						
						bitset<8> fbs(fmask);
						for(NodeID j = 0; j < 8; ++j){							
							if(fbs[7 - j]){
								unsigned char smask = token_v.sptc_sbv[spos++];
								bitset<8> sbs(smask);
								for(NodeID k = 0; k < 8; ++k){
									if(sbs[7 - k]){
										NodeID w = supertoken_r.sptc_v[ (i * 8 + j) * 8 + k  +  1];
										EdgeWeight w_d = supertoken_r.sptc_d[(i * 8 + j) * 8 + k  +  1] + tdis;											
										if( w < numOfVertices){// hashing, can be replaced by 1024 linear probing for efficiency.
											if(ts_vec[w] == ts){
												EdgeWeight current_dis = dis_vec[w] + w_d;
												if(current_dis < distance)
													distance = current_dis;
											}
										}else{
											que_d[que_h] = w_d;
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
		}
		
		
		
		return distance;
	}
	
	
};


class Label {

public:
	vector<index_t> index_;	
	index_t_p* index_p;
	two_index_t_p* two_index_p;


	double GetCurrentTimeSec() {
		struct timeval tv;
		gettimeofday(&tv, NULL);
		return tv.tv_sec + tv.tv_usec * 1e-6;
	}


	Label() {
		index_.resize(numOfVertices);
	}

	~Label() {
		Free();
	}

	EdgeWeight query_p(NodeID s, NodeID t) {
		//
		//EdgeWeight distance = INF_WEIGHT;

		//NodeID *vs = index_p[s].spt_v;
		//NodeID *vt = index_p[t].spt_v;
		//EdgeWeight* ws = index_p[s].spt_d;
		//EdgeWeight* wt = index_p[t].spt_d;

		//_mm_prefetch(vs, _MM_HINT_T0);
		//_mm_prefetch(vt, _MM_HINT_T0);
		//_mm_prefetch(ws, _MM_HINT_T0);
		//_mm_prefetch(wt, _MM_HINT_T0);

		//for (unsigned i = 0, j = 0; ; ) {
		//	if (*(vs + i) == *(vt + j)) {
		//		if (*(vs + i) == numOfVertices) break;  // Sentinel
		//		EdgeWeight td = *(ws + i) + *(wt + j);
		//		if (td < distance) distance = td;
		//		++i;
		//		++j;
		//	}
		//	else {
		//		i += *(vs + i) < *(vt + j) ? 1 : 0;
		//		j += *(vs + i) > *(vt + j) ? 1 : 0;
		//	}
		//}
		//return distance;

		EdgeWeight distance = INF_WEIGHT;

		const index_t_p &idx_s = index_p[s];
		const index_t_p &idx_t = index_p[t];

		_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);

		for (int i = 0, j = 0; ; ) {
			NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

			if (v1 == numOfVertices) break;  // Sentinel

			if (v1 == v2) {
				EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
				if (td < distance) distance = td;
				++i;
				++j;
			} 
			else {
				i += v1 < v2 ? 1 : 0;
				j += v1 > v2 ? 1 : 0;
			}
		}
		return distance;
	}

	EdgeWeight two_query_p_sequential(NodeID s, NodeID t) {
		
		EdgeWeight distance = INF_WEIGHT;
		EdgeWeight ldistance = INF_WEIGHT;

		const two_index_t_p &idx_s = two_index_p[s];
		const two_index_t_p &idx_t = two_index_p[t];

		_mm_prefetch(&idx_s.spt_lv[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_lv[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_ld[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_ld[0], _MM_HINT_T0);

		for (uint8_t i = 0, j = 0; ; ) {
			uint8_t uv8_1 = idx_s.spt_lv[i], uv8_2 = idx_t.spt_lv[j];

			if (uv8_1 == UCHAR_MAX) break;  // Sentinel

			if (uv8_1 == uv8_2) {
				EdgeWeight td = idx_s.spt_ld[i] + idx_t.spt_ld[j];
				if (td < ldistance) ldistance = td;
				++i;
				++j;
			}
			else {
				i += uv8_1 < uv8_2 ? 1 : 0;
				j += uv8_1 > uv8_2 ? 1 : 0;
			}
		}
	
		_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);

		for (int i = 0, j = 0; ; ) {
			NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

			if (v1 == numOfVertices) break;  // Sentinel

			if (v1 == v2) {
				EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
				if (td < distance) distance = td;
				++i;
				++j;
			}
			else {
				i += v1 < v2 ? 1 : 0;
				j += v1 > v2 ? 1 : 0;
			}
		}
		
		if(distance < ldistance) 
			return distance;
		else
			return ldistance;		
	}
	
	EdgeWeight two_query_p_parallel(NodeID s, NodeID t) {
		
		EdgeWeight distance = INF_WEIGHT;
		EdgeWeight ldistance = INF_WEIGHT;

		const two_index_t_p &idx_s = two_index_p[s];
		const two_index_t_p &idx_t = two_index_p[t];

		
		#pragma omp parallel sections
		{
			#pragma omp section
			{
				_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
				_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
				_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
				_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);

				for (int i = 0, j = 0; ; ) {
					NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

					if (v1 == numOfVertices) break;  // Sentinel

					if (v1 == v2) {
						EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
						if (td < distance) distance = td;
						++i;
						++j;
					}
					else {
						i += v1 < v2 ? 1 : 0;
						j += v1 > v2 ? 1 : 0;
					}
				}
			}
			
			#pragma omp section
			{
				_mm_prefetch(&idx_s.spt_lv[0], _MM_HINT_T0);
				_mm_prefetch(&idx_t.spt_lv[0], _MM_HINT_T0);
				_mm_prefetch(&idx_s.spt_ld[0], _MM_HINT_T0);
				_mm_prefetch(&idx_t.spt_ld[0], _MM_HINT_T0);

				for (uint8_t i = 0, j = 0; ; ) {
					uint8_t uv8_1 = idx_s.spt_lv[i], uv8_2 = idx_t.spt_lv[j];

					if (uv8_1 == UCHAR_MAX) break;  // Sentinel

					if (uv8_1 == uv8_2) {
						EdgeWeight td = idx_s.spt_ld[i] + idx_t.spt_ld[j];
						if (td < ldistance) ldistance = td;
						++i;
						++j;
					}
					else {
						i += uv8_1 < uv8_2 ? 1 : 0;
						j += uv8_1 > uv8_2 ? 1 : 0;
					}
				}
			}
		}
		if(distance < ldistance) 
			return distance;
		else
			return ldistance;		
	}
	
	EdgeWeight query_p_with_nums(NodeID s, NodeID t, int k) {
		//
		//EdgeWeight distance = INF_WEIGHT;

		//NodeID *vs = index_p[s].spt_v;
		//NodeID *vt = index_p[t].spt_v;
		//EdgeWeight* ws = index_p[s].spt_d;
		//EdgeWeight* wt = index_p[t].spt_d;

		//_mm_prefetch(vs, _MM_HINT_T0);
		//_mm_prefetch(vt, _MM_HINT_T0);
		//_mm_prefetch(ws, _MM_HINT_T0);
		//_mm_prefetch(wt, _MM_HINT_T0);

		//for (unsigned i = 0, j = 0; ; ) {
		//	if (*(vs + i) == *(vt + j)) {
		//		if (*(vs + i) == numOfVertices) break;  // Sentinel
		//		EdgeWeight td = *(ws + i) + *(wt + j);
		//		if (td < distance) distance = td;
		//		++i;
		//		++j;
		//	}
		//	else {
		//		i += *(vs + i) < *(vt + j) ? 1 : 0;
		//		j += *(vs + i) > *(vt + j) ? 1 : 0;
		//	}
		//}
		//return distance;


		EdgeWeight distance = INF_WEIGHT;

		const index_t_p &idx_s = index_p[s];
		const index_t_p &idx_t = index_p[t];

		_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);
		int k1 = k, k2 = k;
		for (int i = 0, j = 0; ; ) {
			NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

			if (v1 == numOfVertices) break;  // Sentinel

			if (v1 == v2) {
				EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
				if (td < distance) distance = td;
				++i;
				++j;
			}
			else {
				i += v1 < v2 ? 1 : 0;
				j += v1 > v2 ? 1 : 0;
			}

			if (i > k1 || j > k2) break;
		}
		return distance;
	}


	EdgeWeight query(NodeID s, NodeID t) {
		EdgeWeight distance = INF_WEIGHT;
		vector<NodeID>& index_s = index_[s].spt_v;
		vector<EdgeWeight>& index_s_d = index_[s].spt_d;

		vector<NodeID>& index_t = index_[t].spt_v;
		vector<EdgeWeight>& index_t_d = index_[t].spt_d;

		for (int i = 0, j = 0; i < index_s.size(), j < index_t.size(); ) {
			if (index_s[i] == index_t[j]) 
				distance = min(distance, (EdgeWeight)(index_s_d[i++] + index_t_d[j++]));
			else {
				if (index_s[i] < index_t[j])
					++i;
				else
					++j;
			}
		}
		return distance;
	}
	
	

	EdgeWeight query(NodeID s, NodeID t, NodeID& meet, EdgeWeight& dis1, EdgeWeight& dis2) {
		EdgeWeight distance = INF_WEIGHT;
		vector<NodeID>& index_s = index_[s].spt_v;
		vector<EdgeWeight>& index_s_d = index_[s].spt_d;

		vector<NodeID>& index_t = index_[t].spt_v;
		vector<EdgeWeight>& index_t_d = index_[t].spt_d;
		meet = numeric_limits<NodeID>::max();
		dis1 = numeric_limits<EdgeWeight>::max();
		dis2 = numeric_limits<EdgeWeight>::max();
		for (int i = 0, j = 0; i < index_s.size(), j < index_t.size(); ) {
			if (index_s[i] == index_t[j]) {
				if (distance > (EdgeWeight)(index_s_d[i] + index_t_d[j])) {
					distance = (EdgeWeight)(index_s_d[i] + index_t_d[j]);
					meet = index_s[i];
					dis1 = index_s_d[i];
					dis2 = index_t_d[j];
				}
				++i; ++j;
			}
			else {
				if (index_s[i] < index_t[j])
					++i;
				else
					++j;
			}
		}
		return distance;
	}

	/*EdgeWeight query_new(NodeID s, NodeID t, Ordering& ordering) {
		EdgeWeight distance = INF_WEIGHT;
		vector<NodeID>& index_s = index_[s].spt_v;
		vector<EdgeWeight>& index_s_d = index_[s].spt_d;

		vector<NodeID>& index_t = index_[t].spt_v;
		vector<EdgeWeight>& index_t_d = index_[t].spt_d;



		for (int i = 0, j = 0; i < index_s.size(), j < index_t.size(); ) {
			if (index_s[i] == index_t[j])
				distance = min(distance, (EdgeWeight)(index_s_d[i++] + index_t_d[j++]));
			else {
				if (index_s[i] < index_t[j])
					++i;
				else
					++j;
			}
		}
		return distance;
	}
	*/
	double avg_size() {
		
			double total = 0;
		if(index_.size()!=0){
			for (int i = 0; i < numOfVertices; ++i) total += index_[i].spt_v.size();

			double avg = total / numOfVertices - 1; // We do not count the trivial label (V, INF_WEIGHT).
			return avg;
		}
		
		total = 0;
		for (int i = 0; i < numOfVertices; ++i) {
			int unit_count = 0;
			const index_t_p &idx_s = index_p[i];
			for(int j = 0; ;){
				NodeID v = idx_s.spt_v[j++];
				++unit_count;
				if( v == numOfVertices) break;
			}
			total += unit_count;
		}

		double avg = total / numOfVertices - 1; // We do not count the trivial label (V, INF_WEIGHT).

		return avg;
	}
	/*
	NodeID max_size() {
		NodeID maxsize = numeric_limits<NodeID>::min();
		for (int i = 0; i < V; ++i) maxsize = max(maxsize, index_[i].spt_v.size());
		return maxsize;
	}*/

	void append(NodeID v, NodeID root, EdgeWeight distance) {
		index_[v].spt_v.push_back(root);
		index_[v].spt_d.push_back(distance);
	}

	void print_stat() {
		cout << "Average Label Size: " << avg_size() << endl;
		//cout << "Maximum Label Size: " << max_size() << endl;
	}

	void Free() {
		if (index_.size() == 0) return;
		for (int v = 0; v < numOfVertices; ++v) {
			index_[v].spt_v.clear();
			index_[v].spt_d.clear();
		}
		index_.clear();
	}

	void save_labels(const char* save_filename) {
		ofstream ofs(save_filename, ios::binary | ios::out);

		ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
		for (NodeID v = 0; v < numOfVertices; ++v) {
			NodeID isize = index_[v].size();
			ofs.write((const char*)&isize, sizeof(isize));
			for (NodeID i = 0; i < index_[v].size(); ++i) {
				ofs.write((const char*)&index_[v].spt_v[i], sizeof(index_[v].spt_v[i]));
				ofs.write((const char*)&index_[v].spt_d[i], sizeof(index_[v].spt_d[i]));
			}
		}
		ofs.close();
	}

	
	void load_labels(const char* load_filename) {
	/*	for (NodeID v = 0; v < numOfVertices; ++v) {
			free(index_p[v].spt_v);
			free(index_p[v].spt_d);
		}
		*/
		//free(index_p);
		index_p = NULL;

		ifstream ifs(load_filename);
		NodeID isize = 0;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;

		index_p = (index_t_p*)memalign(64, numOfVertices * sizeof(index_t_p));
		


		for (NodeID v = 0; v < numOfVertices; ++v) {
			index_t_p &idx = index_p[v];
			ifs.read((char*)&isize, sizeof(isize));

			idx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
			idx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));

		//	index_[v].spt_v.resize(isize);
		//	index_[v].spt_d.resize(isize);

			for (NodeID i = 0; i < isize; ++i) {
				NodeID hub;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;

				idx.spt_v[i] = hub;
				idx.spt_d[i] = hub_weight;

			}
		}
		ifs.close();

		/*
		index_.clear();
		ifstream ifs(load_filename);
		NodeID isize = 0;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;
		index_.resize(numOfVertices);
		for (NodeID v = 0; v < numOfVertices; ++v) {

			ifs.read((char*)&isize, sizeof(isize));
			index_[v].spt_v.resize(isize);
			index_[v].spt_d.resize(isize);

			for (NodeID i = 0; i < index_[v].size(); ++i) {
				NodeID hub;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));
				index_[v].spt_v[i] = hub;
				index_[v].spt_d[i] = hub_weight;
			}
		}
		ifs.close();
		*/
	}

	void convert_to_fewerbit(){
		
		two_index_p = NULL;
		two_index_p = (two_index_t_p*)memalign(64, numOfVertices * sizeof(two_index_t_p));

		double compressed_size = 0;
		double total_size = 0;
		for (NodeID v = 0; v < numOfVertices; ++v) {
			two_index_t_p &idx = two_index_p[v];
			
			index_t_p &idx_original = index_p[v];
			
			NodeID isize = 0;
			for(NodeID i = 0; idx_original.spt_v[i] < UCHAR_MAX; ++i){
				++isize;
			}
			

			idx.spt_lv = (uint8_t*)memalign(64, (isize + 1) * sizeof(uint8_t));
			idx.spt_ld = (EdgeWeight*)memalign(64, (isize + 1) * sizeof(EdgeWeight));

		//	index_[v].spt_v.resize(isize);
		//	index_[v].spt_d.resize(isize);

			for (NodeID i = 0; i < isize; ++i) {
				uint8_t hub;
				EdgeWeight hub_weight;
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;

				idx.spt_lv[i] = idx_original.spt_v[i];
				idx.spt_ld[i] = idx_original.spt_d[i];

			}
			
			compressed_size += 4 * (isize - 1)  - isize;
			
			idx.spt_lv[isize] = UCHAR_MAX;
			idx.spt_ld[isize] = INF_WEIGHT;
			
			NodeID larger_size = 0;
			for(NodeID i = isize; idx_original.spt_v[i] != numOfVertices; ++i){
				++larger_size;
			}
			
			larger_size++;
			
			idx.spt_v = (NodeID*)memalign(64, larger_size * sizeof(NodeID));
			idx.spt_d = (EdgeWeight*)memalign(64, larger_size * sizeof(EdgeWeight));
			
			for (NodeID i = 0; i < larger_size; ++i) {
				uint8_t hub;
				EdgeWeight hub_weight;
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;

				idx.spt_v[i] = idx_original.spt_v[i + isize];
				idx.spt_d[i] = idx_original.spt_d[i + isize];
			}			
			
			total_size += 4 * (isize - 1 + larger_size) * 2;
			
		}		
		cout << "reduce size :" << compressed_size << " out of " << total_size << " saving " << int(compressed_size * 100 / total_size) << "%" << endl;

	}
	
	void load_labels_with_k(const char* load_filename, int k) {
		/*	for (NodeID v = 0; v < numOfVertices; ++v) {
		free(index_p[v].spt_v);
		free(index_p[v].spt_d);
		}
		*/
		//free(index_p);

		long total_amount = 0;
		long actual_amount = 0;

		index_p = NULL;

		ifstream ifs(load_filename);
		NodeID isize = 0;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;

		index_p = (index_t_p*)memalign(64, numOfVertices * sizeof(index_t_p));



		for (NodeID v = 0; v < numOfVertices; ++v) {
			index_t_p &idx = index_p[v];
			ifs.read((char*)&isize, sizeof(isize));
			int actual_isize = k;
			if (isize > k) actual_isize = k;
			else actual_isize = isize;

			total_amount += isize;
			actual_amount += actual_isize;

			idx.spt_v = (NodeID*)memalign(64, actual_isize * sizeof(NodeID));
			idx.spt_d = (EdgeWeight*)memalign(64, actual_isize * sizeof(EdgeWeight));

			//	index_[v].spt_v.resize(isize);
			//	index_[v].spt_d.resize(isize);

			for (NodeID i = 0; i < isize; ++i) {
				
				NodeID hub;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;

				if (i > actual_isize) continue;
				if (i == actual_isize - 1) {
					idx.spt_v[i] = numOfVertices;
					idx.spt_d[i] = INF_WEIGHT;
				}else {
					idx.spt_v[i] = hub;
					idx.spt_d[i] = hub_weight;
				}
			}
		}
		ifs.close();

		cout << "Total Labels:" << total_amount << endl;
		cout << "Actual Labels:" << actual_amount << endl;

		/*
		index_.clear();
		ifstream ifs(load_filename);
		NodeID isize = 0;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;
		index_.resize(numOfVertices);
		for (NodeID v = 0; v < numOfVertices; ++v) {

		ifs.read((char*)&isize, sizeof(isize));
		index_[v].spt_v.resize(isize);
		index_[v].spt_d.resize(isize);

		for (NodeID i = 0; i < index_[v].size(); ++i) {
		NodeID hub;
		EdgeWeight hub_weight;
		ifs.read((char*)&hub, sizeof(hub));
		ifs.read((char*)&hub_weight, sizeof(hub_weight));
		index_[v].spt_v[i] = hub;
		index_[v].spt_d[i] = hub_weight;
		}
		}
		ifs.close();
		*/
	}


	void save_labels_iteration_stats(const char* save_filename) {

		vector<NodeID> stat(numOfVertices);
		for (NodeID v = 0; v < numOfVertices; ++v) {
			for (NodeID i = 0; i < index_[v].size(); ++i)
				stat[index_[v].spt_v[i]]++;
		}

		ofstream ofs(save_filename);

		for (NodeID v = 0; v < numOfVertices; ++v) {
			ofs << stat[v] << endl;
		}
		ofs.close();
	}


	EdgeWeight query_with_info(NodeID s, NodeID t, query_info& q_info) {

		double stime = GetCurrentTimeSec();

		EdgeWeight distance = INF_WEIGHT;
		vector<NodeID>& index_s = index_[s].spt_v;
		vector<EdgeWeight>& index_s_d = index_[s].spt_d;

		vector<NodeID>& index_t = index_[t].spt_v;
		vector<EdgeWeight>& index_t_d = index_[t].spt_d;		

		q_info.meet_node = numOfVertices;
		double meet_distance;

		for (int i = 0, j = 0; i < index_s.size(), j < index_t.size(); ) {
			if (index_s[i] == index_t[j]) {
				meet_distance = (EdgeWeight)(index_s_d[i++] + index_t_d[j++]);
				if ( distance >  meet_distance) {
					distance = meet_distance;
					q_info.meet_node = index_s[i];
				}
			}
			else {
				if (index_s[i] < index_t[j])
					++i;
				else
					++j;
			}
		};

		stime = GetCurrentTimeSec() - stime;

		q_info.time_cost = stime;

		if (index_s.size() < index_t.size())
			q_info.search_len = index_s.size();
		else
			q_info.search_len = index_t.size();

		return distance;
	}


};

class PLabel {

public:
	vector<index_t_path> index_;
	index_t_path_p* index_p;


	double GetCurrentTimeSec() {
		struct timeval tv;
		gettimeofday(&tv, NULL);
		return tv.tv_sec + tv.tv_usec * 1e-6;
	}


	PLabel() {
		index_.resize(numOfVertices);
	}

	~PLabel() {
		Free();
	}

	EdgeWeight query_p(NodeID s, NodeID t) {

		//EdgeWeight distance = INF_WEIGHT;

		//NodeID *vs = index_p[s].spt_v;
		//NodeID *vt = index_p[t].spt_v;
		//EdgeWeight* ws = index_p[s].spt_d;
		//EdgeWeight* wt = index_p[t].spt_d;

		//_mm_prefetch(vs, _MM_HINT_T0);
		//_mm_prefetch(vt, _MM_HINT_T0);
		//_mm_prefetch(ws, _MM_HINT_T0);
		//_mm_prefetch(wt, _MM_HINT_T0);

		//for (unsigned i = 0, j = 0; ; ) {
		//	if (*(vs + i) == *(vt + j)) {
		//		if (*(vs + i) == numOfVertices) break;  // Sentinel
		//		EdgeWeight td = *(ws + i) + *(wt + j);
		//		if (td < distance) distance = td;
		//		++i;
		//		++j;
		//	}
		//	else {
		//		i += *(vs + i) < *(vt + j) ? 1 : 0;
		//		j += *(vs + i) > *(vt + j) ? 1 : 0;
		//	}
		//}
		//return distance;


		EdgeWeight distance = INF_WEIGHT;
		NodeID meet;

		const index_t_path_p &idx_s = index_p[s];
		const index_t_path_p &idx_t = index_p[t];

		_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);

		for (int i = 0, j = 0; ; ) {
			NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

			if (v1 == numOfVertices) break;  // Sentinel

			if (v1 == v2) {
				EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
				if (td < distance) {
					distance = td;
				}
				++i;
				++j;
			}
			else {
				i += v1 < v2 ? 1 : 0;
				j += v1 > v2 ? 1 : 0;
			}
		}

		return distance;
	}

	EdgeWeight query(NodeID s, NodeID t) {
		EdgeWeight distance = INF_WEIGHT;
		vector<NodeID>& index_s = index_[s].spt_v;
		vector<EdgeWeight>& index_s_d = index_[s].spt_d;

		vector<NodeID>& index_t = index_[t].spt_v;
		vector<EdgeWeight>& index_t_d = index_[t].spt_d;

		for (int i = 0, j = 0; i < index_s.size(), j < index_t.size(); ) {
			if (index_s[i] == index_t[j])
				distance = min(distance, (EdgeWeight)(index_s_d[i++] + index_t_d[j++]));
			else {
				if (index_s[i] < index_t[j])
					++i;
				else
					++j;
			}
		}
		return distance;
	}

	EdgeWeight query(NodeID s, NodeID t, NodeID& meet, EdgeWeight& dis1, EdgeWeight& dis2) {
		EdgeWeight distance = INF_WEIGHT;
		vector<NodeID>& index_s = index_[s].spt_v;
		vector<EdgeWeight>& index_s_d = index_[s].spt_d;

		vector<NodeID>& index_t = index_[t].spt_v;
		vector<EdgeWeight>& index_t_d = index_[t].spt_d;
		meet = numeric_limits<NodeID>::max();
		dis1 = numeric_limits<EdgeWeight>::max();
		dis2 = numeric_limits<EdgeWeight>::max();
		for (int i = 0, j = 0; i < index_s.size(), j < index_t.size(); ) {
			if (index_s[i] == index_t[j]) {
				if (distance >(EdgeWeight)(index_s_d[i] + index_t_d[j])) {
					distance = (EdgeWeight)(index_s_d[i] + index_t_d[j]);
					meet = index_s[i];
					dis1 = index_s_d[i];
					dis2 = index_t_d[j];
				}
				++i; ++j;
			}
			else {
				if (index_s[i] < index_t[j])
					++i;
				else
					++j;
			}
		}
		return distance;
	}


	EdgeWeight query_path(NodeID s, NodeID t, vector<NodeID>& rank, vector<NodeID>& inv) {


		EdgeWeight distance = INF_WEIGHT;
		NodeID meetnode = numOfVertices;
		NodeID s_parent;
		NodeID t_parent;

		const index_t_path_p &idx_s = index_p[s];
		const index_t_path_p &idx_t = index_p[t];

		_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_p[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_p[0], _MM_HINT_T0);

		for (int i = 0, j = 0; ; ) {
			NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

			if (v1 == numOfVertices) break;  // Sentinel

			if (v1 == v2) {
				EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
				if (td < distance) {
					distance = td;
				//	if (v1 < meetnode) {
						meetnode = v1;
						s_parent = idx_s.spt_p[i];
						t_parent = idx_t.spt_p[j];
					//}
				}
				++i;
				++j;
			}
			else {
				i += v1 < v2 ? 1 : 0;
				j += v1 > v2 ? 1 : 0;
			}
		}


		//Next, retrieve path from s - meetnode and meetnode - t.
		vector<NodeID> path_from_s;
		vector<NodeID> path_to_t;
		path_from_s.push_back(s_parent);
		path_to_t.push_back(t_parent);
		
		int operation = 0;

		/*	if (s == 194569 && t == 20072)
		cout << "debug." << " meet: " << meetnode << " sparent:" << s_parent << " tparent:" << t_parent <<  endl;*/

		NodeID inv_meetnode = inv[meetnode];

		while (path_from_s.back() != inv_meetnode) {
			/*if (s == 194569 && t == 20072)
			cout << "s meet:" << path_from_s.back() << endl;*/
			const index_t_path_p &idx_from_s = index_p[path_from_s.back()];

			_mm_prefetch(&idx_from_s.spt_v[0], _MM_HINT_T0);
			_mm_prefetch(&idx_from_s.spt_p[0], _MM_HINT_T0);
			 
		//	vector<NodeID>& index_from_s = index_[path_from_s.back()].spt_v;
			for (int i = 0; ; ++i) {
				operation++;
				if (idx_from_s.spt_v[i] == numOfVertices) break;
				if (idx_from_s.spt_v[i] == meetnode) {
					path_from_s.push_back(idx_from_s.spt_p[i]);
					break;
				}
			}
		}

		while (path_to_t.back() != inv_meetnode) {
			/*if (s == 194569 && t == 20072)
			cout << "t meet:" << path_to_t.back() << endl;*/
		//	vector<NodeID>& index_to_t = index_[path_to_t.back()].spt_v;
			const index_t_path_p &idx_to_t = index_p[path_to_t.back()];
			
			_mm_prefetch(&idx_to_t.spt_v[0], _MM_HINT_T0);
			_mm_prefetch(&idx_to_t.spt_p[0], _MM_HINT_T0);
			for (int i = 0; ; ++i) {
				operation++;
				if (idx_to_t.spt_v[i] == numOfVertices) break;
				if (idx_to_t.spt_v[i] == meetnode) {
					path_to_t.push_back(idx_to_t.spt_p[i]);
					break;
				}
			}
		}

		distance = 0;
		distance += path_from_s.size() + path_to_t.size();

//		return distance;
		return distance;
		//EdgeWeight distance = INF_WEIGHT;
		//vector<NodeID>& index_s = index_[s].spt_v;
		//vector<EdgeWeight>& index_s_d = index_[s].spt_d;

		//vector<NodeID>& bindex_t = index_[t].spt_v;
		//vector<EdgeWeight>& bindex_t_d = index_[t].spt_d;


		//NodeID meetnode = numOfVertices;
		//int s_parent;
		//int t_parent;
		//for (int i = 0, j = 0; i < index_s.size(), j < bindex_t.size(); ) {
		//	if (index_s[i] == bindex_t[j]) {
		//		if (distance >(EdgeWeight)(index_s_d[i] + bindex_t_d[j])) {
		//			distance = (EdgeWeight)(index_s_d[i] + bindex_t_d[j]);
		//			if (index_s[i] < meetnode) {
		//				meetnode = index_s[i];
		//				s_parent = index_[s].spt_p[i];
		//				t_parent = index_[t].spt_p[j];
		//			}
		//		}
		//		//distance = min(distance, (EdgeWeight)(index_s_d[i] + bindex_t_d[j]));
		//		++i;
		//		++j;
		//	}
		//	else {
		//		if (index_s[i] < bindex_t[j])
		//			++i;
		//		else
		//			++j;
		//	}
		//}

		////Next, retrieve path from s - meetnode and meetnode - t.
		//vector<NodeID> path_from_s;
		//vector<NodeID> path_to_t;
		//path_from_s.push_back(s_parent);
		//path_to_t.push_back(t_parent);

		///*	if (s == 194569 && t == 20072)
		//cout << "debug." << " meet: " << meetnode << " sparent:" << s_parent << " tparent:" << t_parent <<  endl;*/

		//while (path_from_s.back() != inv[meetnode]) {
		//	/*if (s == 194569 && t == 20072)
		//	cout << "s meet:" << path_from_s.back() << endl;*/
		//	vector<NodeID>& index_from_s = index_[path_from_s.back()].spt_v;
		//	for (int i = 0; i < index_from_s.size(); ++i) {
		//		if (index_from_s[i] == meetnode) {
		//			path_from_s.push_back(index_[path_from_s.back()].spt_p[i]);
		//			break;
		//		}
		//	}
		//}

		//while (path_to_t.back() != inv[meetnode]) {
		//	/*if (s == 194569 && t == 20072)
		//	cout << "t meet:" << path_to_t.back() << endl;*/
		//	vector<NodeID>& index_to_t = index_[path_to_t.back()].spt_v;
		//	for (int i = 0; i < index_to_t.size(); ++i) {
		//		if (index_to_t[i] == meetnode) {
		//			path_to_t.push_back(index_[path_to_t.back()].spt_p[i]);
		//			break;
		//		}
		//	}
		//}

		////for (int i = 0; i < path_from_s.size(); ++i)
		////	path_from_s[i] = inv[path_from_s[i]];
		////for (int i = 0; i < path_to_t.size(); ++i)
		////	path_to_t[i] = inv[path_to_t[i]];

		//return path_from_s.size() + path_to_t.size();
	}

	EdgeWeight query_path_check(NodeID s, NodeID t, vector<NodeID>& rank, vector<NodeID>& inv) {


		EdgeWeight distance = INF_WEIGHT;
		NodeID meetnode = numOfVertices;
		NodeID s_parent;
		NodeID t_parent;

		const index_t_path_p &idx_s = index_p[s];
		const index_t_path_p &idx_t = index_p[t];

		_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_p[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_p[0], _MM_HINT_T0);

		for (int i = 0, j = 0; ; ) {
			NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

			if (v1 == numOfVertices) break;  // Sentinel

			if (v1 == v2) {
				EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
				if (td < distance) {
					distance = td;
				//	if (v1 < meetnode) {
						meetnode = v1;
						s_parent = idx_s.spt_p[i];
						t_parent = idx_t.spt_p[j];
					//}
				}
				++i;
				++j;
			}
			else {
				i += v1 < v2 ? 1 : 0;
				j += v1 > v2 ? 1 : 0;
			}
		}
		
		NodeID inv_meetnode = inv[meetnode];
		//Next, retrieve path from s - meetnode and meetnode - t.
		vector<NodeID> path_from_s;
		vector<NodeID> path_to_t;
		if(s !=inv_meetnode)
			path_from_s.push_back(s);
		path_from_s.push_back(s_parent);
		if (t != inv_meetnode)	
			path_to_t.push_back(t);

		path_to_t.push_back(t_parent);

		/*	if (s == 194569 && t == 20072)
		cout << "debug." << " meet: " << meetnode << " sparent:" << s_parent << " tparent:" << t_parent <<  endl;*/


		while (path_from_s.back() != inv_meetnode) {
			/*if (s == 194569 && t == 20072)
			cout << "s meet:" << path_from_s.back() << endl;*/
			const index_t_path_p &idx_from_s = index_p[path_from_s.back()];

			_mm_prefetch(&idx_from_s.spt_v[0], _MM_HINT_T0);
			_mm_prefetch(&idx_from_s.spt_p[0], _MM_HINT_T0);

			//	vector<NodeID>& index_from_s = index_[path_from_s.back()].spt_v;
			for (int i = 0; ; ++i) {
				if (idx_from_s.spt_v[i] == numOfVertices) break;
				if (idx_from_s.spt_v[i] == meetnode) {
					path_from_s.push_back(idx_from_s.spt_p[i]);
					break;
				}
			}
		}

		while (path_to_t.back() != inv_meetnode) {
			/*if (s == 194569 && t == 20072)
			cout << "t meet:" << path_to_t.back() << endl;*/
			//	vector<NodeID>& index_to_t = index_[path_to_t.back()].spt_v;
			const index_t_path_p &idx_to_t = index_p[path_to_t.back()];
			_mm_prefetch(&idx_to_t.spt_v[0], _MM_HINT_T0);
			_mm_prefetch(&idx_to_t.spt_p[0], _MM_HINT_T0);
			for (int i = 0; ; ++i) {
				if (idx_to_t.spt_v[i] == numOfVertices) break;
				if (idx_to_t.spt_v[i] == meetnode) {
					path_to_t.push_back(idx_to_t.spt_p[i]);
					break;
				}
			}
		}

		//return distance;

		EdgeWeight alldis = 0;

		if (path_from_s.size() == 1)
			if (s != inv_meetnode)
				alldis += query_p(s, inv_meetnode);

		if (path_to_t.size() == 1)
			if (t != inv_meetnode)
				alldis += query_p(t, inv_meetnode);

		for (int i = 0; i < path_from_s.size() - 1; ++i) {
			alldis += query_p(path_from_s[i], path_from_s[i + 1]);
			//cout << "s " << path_from_s[i] << "," << path_from_s[i + 1] << endl;
		}
		for (int i = 0; i < path_to_t.size() - 1; ++i) {
			alldis += query_p(path_to_t[i], path_to_t[i + 1]);

			//cout <<"t " <<  path_to_t[i] << "," << path_to_t[i + 1] << endl;
		}
		/*if (distance != alldis)
			cout << "a?" << endl;*/
		//cout << distance << "," << alldis << "," << path_from_s.size() + path_to_t.size() << endl;
//		cout << s << "," << t << "," << inv_meetnode << "   " << distance << "vs." << alldis << endl;

		return distance;
	}


	//EdgeWeight query_path_check(NodeID s, NodeID t, vector<NodeID>& rank, vector<NodeID>& inv) {


	//	EdgeWeight distance = INF_WEIGHT;
	//	NodeID meetnode = numOfVertices;
	//	NodeID s_parent;
	//	NodeID t_parent;

	//	const index_t_path_p &idx_s = index_p[s];
	//	const index_t_path_p &idx_t = index_p[t];

	//	_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
	//	_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
	//	_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
	//	_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);
	//	_mm_prefetch(&idx_s.spt_p[0], _MM_HINT_T0);
	//	_mm_prefetch(&idx_t.spt_p[0], _MM_HINT_T0);

	//	for (int i = 0, j = 0; ; ) {
	//		NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

	//		if (v1 == numOfVertices) break;  // Sentinel

	//		if (v1 == v2) {
	//			EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
	//			if (td < distance) {
	//				distance = td;
	//				if (v1 < meetnode) {
	//					meetnode = v1;
	//					s_parent = idx_s.spt_p[i];
	//					t_parent = idx_t.spt_p[j];
	//				}
	//			}
	//			++i;
	//			++j;
	//		}
	//		else {
	//			i += v1 < v2 ? 1 : 0;
	//			j += v1 > v2 ? 1 : 0;
	//		}
	//	}

	//	//Next, retrieve path from s - meetnode and meetnode - t.
	//	vector<NodeID> path_from_s;
	//	vector<NodeID> path_to_t;
	//	path_from_s.push_back(s_parent);
	//	path_to_t.push_back(t_parent);

	//	/*	if (s == 194569 && t == 20072)
	//	cout << "debug." << " meet: " << meetnode << " sparent:" << s_parent << " tparent:" << t_parent <<  endl;*/

	//	NodeID inv_meetnode = inv[meetnode];

	//	while (path_from_s.back() != inv_meetnode) {
	//		/*if (s == 194569 && t == 20072)
	//		cout << "s meet:" << path_from_s.back() << endl;*/
	//		const index_t_path_p &idx_from_s = index_p[path_from_s.back()];

	//		_mm_prefetch(&idx_from_s.spt_v[0], _MM_HINT_T0);
	//		_mm_prefetch(&idx_from_s.spt_p[0], _MM_HINT_T0);

	//		//	vector<NodeID>& index_from_s = index_[path_from_s.back()].spt_v;
	//		for (int i = 0; ; ++i) {
	//			if (idx_from_s.spt_v[i] == numOfVertices) break;
	//			if (idx_from_s.spt_v[i] == meetnode) {
	//				path_from_s.push_back(idx_from_s.spt_p[i]);
	//				break;
	//			}
	//		}
	//	}

	//	while (path_to_t.back() != inv_meetnode) {
	//		/*if (s == 194569 && t == 20072)
	//		cout << "t meet:" << path_to_t.back() << endl;*/
	//		//	vector<NodeID>& index_to_t = index_[path_to_t.back()].spt_v;
	//		const index_t_path_p &idx_to_t = index_p[path_to_t.back()];
	//		_mm_prefetch(&idx_to_t.spt_v[0], _MM_HINT_T0);
	//		_mm_prefetch(&idx_to_t.spt_p[0], _MM_HINT_T0);
	//		for (int i = 0; ; ++i) {
	//			if (idx_to_t.spt_v[i] == numOfVertices) break;
	//			if (idx_to_t.spt_v[i] == meetnode) {
	//				path_to_t.push_back(idx_to_t.spt_p[i]);
	//				break;
	//			}
	//		}
	//	}

	//	EdgeWeight path_from_s = 0;
	//	for (int i = 0; i < path_from_s.size(); ++i) {

	//	}

	//	

	//	return distance;

	//	
	//}


	/*EdgeWeight query_new(NodeID s, NodeID t, Ordering& ordering) {
	EdgeWeight distance = INF_WEIGHT;
	vector<NodeID>& index_s = index_[s].spt_v;
	vector<EdgeWeight>& index_s_d = index_[s].spt_d;

	vector<NodeID>& index_t = index_[t].spt_v;
	vector<EdgeWeight>& index_t_d = index_[t].spt_d;



	for (int i = 0, j = 0; i < index_s.size(), j < index_t.size(); ) {
	if (index_s[i] == index_t[j])
	distance = min(distance, (EdgeWeight)(index_s_d[i++] + index_t_d[j++]));
	else {
	if (index_s[i] < index_t[j])
	++i;
	else
	++j;
	}
	}
	return distance;
	}
	*/
	double avg_size() {
		double total = 0;
		for (int i = 0; i < numOfVertices; ++i) total += index_[i].spt_v.size();

		double avg = total / numOfVertices - 1; // We do not count the trivial label (V, INF_WEIGHT).

		return avg;
	}
	/*
	NodeID max_size() {
	NodeID maxsize = numeric_limits<NodeID>::min();
	for (int i = 0; i < V; ++i) maxsize = max(maxsize, index_[i].spt_v.size());
	return maxsize;
	}*/

	void append(NodeID v, NodeID root, EdgeWeight distance) {
		index_[v].spt_v.push_back(root);
		index_[v].spt_d.push_back(distance);
	}

	void print_stat() {
		cout << "Average Label Size: " << avg_size() << endl;
		//cout << "Maximum Label Size: " << max_size() << endl;
	}

	void Free() {
		if (index_.size() == 0) return;
		for (int v = 0; v < numOfVertices; ++v) {
			index_[v].spt_v.clear();
			index_[v].spt_d.clear();
		}
		index_.clear();
	}

	void save_labels(const char* save_filename) {
		ofstream ofs(save_filename, ios::binary | ios::out);

		ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
		for (NodeID v = 0; v < numOfVertices; ++v) {
			NodeID isize = index_[v].size();
			ofs.write((const char*)&isize, sizeof(isize));
			for (NodeID i = 0; i < index_[v].size(); ++i) {
				ofs.write((const char*)&index_[v].spt_v[i], sizeof(index_[v].spt_v[i]));
				ofs.write((const char*)&index_[v].spt_p[i], sizeof(index_[v].spt_p[i]));
				ofs.write((const char*)&index_[v].spt_d[i], sizeof(index_[v].spt_d[i]));
			}
		}
		ofs.close();
	}

	void load_labels(const char* load_filename) {
		/*	for (NodeID v = 0; v < numOfVertices; ++v) {
		free(index_p[v].spt_v);
		free(index_p[v].spt_d);
		}
		*/
		//free(index_p);
		index_p = NULL;

		ifstream ifs(load_filename);
		NodeID isize = 0;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;

		index_p = (index_t_path_p*)memalign(64, numOfVertices * sizeof(index_t_path_p));

		for (NodeID v = 0; v < numOfVertices; ++v) {
			index_t_path_p &idx = index_p[v];
			ifs.read((char*)&isize, sizeof(isize));

			idx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
			idx.spt_p = (NodeID*)memalign(64, isize * sizeof(NodeID));
			idx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));

			//	index_[v].spt_v.resize(isize);
			//	index_[v].spt_d.resize(isize);

			for (NodeID i = 0; i < isize; ++i) {
				NodeID hub;
				NodeID hub_parent;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_parent, sizeof(hub_parent));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;
				idx.spt_v[i] = hub;
				idx.spt_p[i] = hub_parent;
				idx.spt_d[i] = hub_weight;
			}
		}
		ifs.close();

		/*
		index_.clear();
		ifstream ifs(load_filename);
		NodeID isize = 0;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;
		index_.resize(numOfVertices);
		for (NodeID v = 0; v < numOfVertices; ++v) {

		ifs.read((char*)&isize, sizeof(isize));
		index_[v].spt_v.resize(isize);
		index_[v].spt_d.resize(isize);

		for (NodeID i = 0; i < index_[v].size(); ++i) {
		NodeID hub;
		EdgeWeight hub_weight;
		ifs.read((char*)&hub, sizeof(hub));
		ifs.read((char*)&hub_weight, sizeof(hub_weight));
		index_[v].spt_v[i] = hub;
		index_[v].spt_d[i] = hub_weight;
		}
		}
		ifs.close();
		*/
	}

	void save_labels_iteration_stats(const char* save_filename) {

		vector<NodeID> stat(numOfVertices);
		for (NodeID v = 0; v < numOfVertices; ++v) {
			for (NodeID i = 0; i < index_[v].size(); ++i)
				stat[index_[v].spt_v[i]]++;
		}

		ofstream ofs(save_filename);

		for (NodeID v = 0; v < numOfVertices; ++v) {
			ofs << stat[v] << endl;
		}
		ofs.close();
	}

	EdgeWeight query_with_info(NodeID s, NodeID t, query_info& q_info) {

		double stime = GetCurrentTimeSec();

		EdgeWeight distance = INF_WEIGHT;
		vector<NodeID>& index_s = index_[s].spt_v;
		vector<EdgeWeight>& index_s_d = index_[s].spt_d;

		vector<NodeID>& index_t = index_[t].spt_v;
		vector<EdgeWeight>& index_t_d = index_[t].spt_d;

		q_info.meet_node = numOfVertices;
		double meet_distance;

		for (int i = 0, j = 0; i < index_s.size(), j < index_t.size(); ) {
			if (index_s[i] == index_t[j]) {
				meet_distance = (EdgeWeight)(index_s_d[i++] + index_t_d[j++]);
				if (distance >  meet_distance) {
					distance = meet_distance;
					q_info.meet_node = index_s[i];
				}
			}
			else {
				if (index_s[i] < index_t[j])
					++i;
				else
					++j;
			}
		};

		stime = GetCurrentTimeSec() - stime;

		q_info.time_cost = stime;

		if (index_s.size() < index_t.size())
			q_info.search_len = index_s.size();
		else
			q_info.search_len = index_t.size();

		return distance;
	}


};


class DLabel : public Label {
public:
	vector<index_t> bindex_; // Backward labels.

	index_t_p* bindex_p;
	
	
	two_index_t_p* b_two_index_p;


	DLabel() {
		index_.resize(numOfVertices);
		bindex_.resize(numOfVertices);
	}
	
	~DLabel() {  
		Free();
	}

	EdgeWeight query(NodeID s, NodeID t) {

		EdgeWeight distance = INF_WEIGHT;
		vector<NodeID>& index_s = index_[s].spt_v;
		vector<EdgeWeight>& index_s_d = index_[s].spt_d;

		vector<NodeID>& bindex_t = bindex_[t].spt_v;
		vector<EdgeWeight>& bindex_t_d = bindex_[t].spt_d;

		for (int i = 0, j = 0; i < index_s.size(), j < bindex_t.size(); ) {
			if (index_s[i] == bindex_t[j]) {
				distance = min(distance, (EdgeWeight)(index_s_d[i] + bindex_t_d[j]));
				++i;
				++j;
			}
			else {
				if (index_s[i] < bindex_t[j])
					++i;
				else
					++j;
			}
		}

		return distance;
	}

	EdgeWeight query(NodeID s, NodeID t, NodeID& meet, EdgeWeight& dis1, EdgeWeight& dis2) {
		EdgeWeight distance = INF_WEIGHT;
		vector<NodeID>& index_s = index_[s].spt_v;
		vector<EdgeWeight>& index_s_d = index_[s].spt_d;

		vector<NodeID>& bindex_t = bindex_[t].spt_v;
		vector<EdgeWeight>& bindex_t_d = bindex_[t].spt_d;

		meet = numeric_limits<NodeID>::max();
		dis1 = numeric_limits<EdgeWeight>::max();
		dis2 = numeric_limits<EdgeWeight>::max();
		for (int i = 0, j = 0; i < index_s.size(), j < bindex_t.size(); ) {
			if (index_s[i] == bindex_t[j]) {
				if (distance > (EdgeWeight)(index_s_d[i] + bindex_t_d[j])) {
					distance = (EdgeWeight)(index_s_d[i] + bindex_t_d[j]);
					meet = index_s[i];
					dis1 = index_s_d[i];
					dis2 = bindex_t_d[j];
				}
				++i;
				++j;
			}
			else {
				if (index_s[i] < bindex_t[j])
					++i;
				else
					++j;
			}
		}

		return distance;
	}
	
	inline EdgeWeight query_p(NodeID s, NodeID t) {

		//EdgeWeight distance = INF_WEIGHT;
		//
		////const index_t_p &idx_s = index_p[s];
		////const index_t_p &idx_t = bindex_p[t];

		//NodeID *vs = index_p[s].spt_v;
		//NodeID *vt = bindex_p[t].spt_v;
		//EdgeWeight* ws = index_p[s].spt_d;
		//EdgeWeight* wt = bindex_p[t].spt_d;


		//_mm_prefetch(vs, _MM_HINT_T0);
		//_mm_prefetch(vt, _MM_HINT_T0);
		//_mm_prefetch(ws, _MM_HINT_T0); 
		//_mm_prefetch(wt, _MM_HINT_T0);

		//for (unsigned i = 0, j = 0; ; ) {
		//	if (*(vs + i) == *(vt + j)) {
		//		if (*(vs + i) == numOfVertices) break;  // Sentinel
		//		EdgeWeight td = *(ws + i) + *(wt + j);
		//		if (td < distance) distance = td;
		//		++i;
		//		++j;
		//	}
		//	else {
		//		i += *(vs + i) < *(vt + j) ? 1 : 0;
		//		j += *(vs + i) > *(vt + j) ? 1 : 0;
		//	}
		//}
		//return distance;

		
		EdgeWeight distance = INF_WEIGHT;

		const index_t_p &idx_s = index_p[s];
		const index_t_p &idx_t = bindex_p[t];


		_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);

		for (int i = 0, j = 0; ; ) {
			NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

			if (v1 == v2) {
				if (v1 == numOfVertices) break;  // Sentinel
				EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
				if (td < distance) distance = td;
				++i;
				++j;
			}
			else {
				i += v1 < v2 ? 1 : 0;
				j += v1 > v2 ? 1 : 0;
			}
		}
		return distance; 
	}

	EdgeWeight query_with_info(NodeID s, NodeID t, query_info& q_info) {

		double stime = GetCurrentTimeSec();

		EdgeWeight distance = INF_WEIGHT;
		vector<NodeID>& index_s = index_[s].spt_v;
		vector<EdgeWeight>& index_s_d = index_[s].spt_d;

	//	vector<NodeID>& index_t = index_[t].spt_v;
	//	vector<EdgeWeight>& index_t_d = index_[t].spt_d;	
		vector<NodeID>& bindex_t = bindex_[t].spt_v;
		vector<EdgeWeight>& bindex_t_d = bindex_[t].spt_d;

		q_info.meet_node = numOfVertices;
		double meet_distance;

		for (int i = 0, j = 0; i < index_s.size(), j < bindex_t.size(); ) {
			if (index_s[i] == bindex_t[j]) {
				meet_distance = (EdgeWeight)(index_s_d[i++] + bindex_t[j++]);
				if (distance >  meet_distance) {
					distance = meet_distance;
					q_info.meet_node = index_s[i];
				}
			}
			else {
				if (index_s[i] < bindex_t[j])
					++i;
				else
					++j;
			}
		};

		stime = GetCurrentTimeSec() - stime;

		q_info.time_cost = stime;

		if (index_s.size() < bindex_t.size())
			q_info.search_len = index_s.size();
		else
			q_info.search_len = bindex_t.size();

		return distance;
	}


	void append(NodeID v, NodeID root, EdgeWeight distance, bool forward) { // forward(backward) search from root to vertex v.
		if (forward) { // forward search from root to vertex v, hence append (root, distance) to backward index of vertex v.
			bindex_[v].spt_v.push_back(root);
			bindex_[v].spt_d.push_back(distance);
		}
		else { // backward search from root to vertex v, hence append (root, distance) to forward index of vertex v.
			index_[v].spt_v.push_back(root);
			index_[v].spt_d.push_back(distance);
		}
	}

	void Free() {
		if (index_.size() == 0 || bindex_.size() == 0) return;
		for (int v = 0; v < numOfVertices; ++v) {
			index_[v].spt_v.clear();
			index_[v].spt_d.clear();
			if (DIRECTED_FLAG == true) {
				bindex_[v].spt_v.clear();
				bindex_[v].spt_d.clear();
			}
		}
		index_.clear();
		bindex_.clear();
	}

	double avg_size() {
		double total = 0;
		for (int i = 0; i < numOfVertices; ++i) {
			total += index_[i].spt_v.size() ;
			total += bindex_[i].spt_v.size();
		}

		double avg = total / numOfVertices / 2 - 1; // We do not count the trivial labels (V, INF_WEIGHT).

		return avg;
	}

	void print_stat() {
		cout << "Average Label Size: " << avg_size() << endl;
		//cout << "Maximum Label Size: " << max_size() << endl;
	}

	void save_labels(const char* save_filename) {
		ofstream ofs(save_filename, ios::binary | ios::out);

		ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
		for (NodeID v = 0; v < numOfVertices; ++v) {
			int isize = index_[v].size();
			ofs.write((const char*)&isize, sizeof(isize));
			for (NodeID i = 0; i < index_[v].size(); ++i) {
				ofs.write((const char*)&index_[v].spt_v[i], sizeof(index_[v].spt_v[i]));
				ofs.write((const char*)&index_[v].spt_d[i], sizeof(index_[v].spt_d[i]));
			}
			int bisize = bindex_[v].size();
			ofs.write((const char*)&bisize, sizeof(bisize));
			for (NodeID i = 0; i < bindex_[v].size(); ++i) {
				ofs.write((const char*)&bindex_[v].spt_v[i], sizeof(bindex_[v].spt_v[i]));
				ofs.write((const char*)&bindex_[v].spt_d[i], sizeof(bindex_[v].spt_d[i]));
			}
		}
		ofs.close();
	}

	void load_labels(const char* load_filename) {
		cout << "Loading Labels" << endl;
/*
		for (NodeID v = 0; v < numOfVertices; ++v) {
			free(index_p[v].spt_v);
			free(index_p[v].spt_d);
		}*/

		//free(index_p);
		index_p = NULL;
		bindex_p = NULL;

		ifstream ifs(load_filename);
		NodeID isize = 0;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;

		index_p = (index_t_p*)memalign(64, numOfVertices * sizeof(index_t_p));
		bindex_p = (index_t_p*)memalign(64, numOfVertices * sizeof(index_t_p));

		cout << numOfVertices << " vertices." << endl;

		for (NodeID v = 0; v < numOfVertices; ++v) {
			index_t_p &idx = index_p[v];
			ifs.read((char*)&isize, sizeof(isize));

			idx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
			idx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));

			for (NodeID i = 0; i < isize; ++i) {
				NodeID hub;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;

				idx.spt_v[i] = hub;
				idx.spt_d[i] = hub_weight;

			}

			//	index_[v].spt_v.resize(isize);
			//	index_[v].spt_d.resize(isize);


			index_t_p &bidx = bindex_p[v];
			ifs.read((char*)&isize, sizeof(isize));
			bidx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
			bidx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));

			for (NodeID i = 0; i < isize; ++i) {
				NodeID hub;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;
				bidx.spt_v[i] = hub;
				bidx.spt_d[i] = hub_weight;

			}

		}
		ifs.close();

		/*
		index_.clear();
		bindex_.clear();
		ifs.open(load_filename, ios::binary | ios::in);
		
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;
		index_.resize(numOfVertices);
		bindex_.resize(numOfVertices);
		for (NodeID v = 0; v < numOfVertices; ++v) {

			ifs.read((char*)&isize, sizeof(isize));
			index_[v].spt_v.resize(isize);
			index_[v].spt_d.resize(isize);
			for (NodeID i = 0; i < index_[v].size(); ++i) {
				NodeID hub;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));
				index_[v].spt_v[i] = hub;
				index_[v].spt_d[i] = hub_weight;
			}

			ifs.read((char*)&isize, sizeof(isize));
			bindex_[v].spt_v.resize(isize);
			bindex_[v].spt_d.resize(isize); 			
			for (NodeID i = 0; i < bindex_[v].size(); ++i) {
				NodeID hub;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));
				bindex_[v].spt_v[i] = hub;
				bindex_[v].spt_d[i] = hub_weight;
			}
		}
		ifs.close();
		*/
	/*	for (int i = 0; i < numOfVertices; ++i) {
			for (int j = 0; j < index_[i].size(); ++j)
				if (index_[i].spt_v[j] != index_p[i].spt_v[j])
					cout << "warning." << endl;
		}*/
		
	}


	
	void convert_to_fewerbit(){
		
		two_index_p = NULL;
		b_two_index_p = NULL;
		two_index_p = (two_index_t_p*)memalign(64, numOfVertices * sizeof(two_index_t_p));
		b_two_index_p = (two_index_t_p*)memalign(64, numOfVertices * sizeof(two_index_t_p));

		for (NodeID v = 0; v < numOfVertices; ++v) {
			two_index_t_p &idx = two_index_p[v];
			
			index_t_p &idx_original = index_p[v];
			
			NodeID isize = 0;
			for(NodeID i = 0; idx_original.spt_v[i] < UCHAR_MAX; ++i){
				++isize;
			}
			

			idx.spt_lv = (uint8_t*)memalign(64, (isize + 1) * sizeof(uint8_t));
			idx.spt_ld = (EdgeWeight*)memalign(64, (isize + 1) * sizeof(EdgeWeight));

		//	index_[v].spt_v.resize(isize);
		//	index_[v].spt_d.resize(isize);

			for (NodeID i = 0; i < isize; ++i) {
				uint8_t hub;
				EdgeWeight hub_weight;
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;

				idx.spt_lv[i] = idx_original.spt_v[i];
				idx.spt_ld[i] = idx_original.spt_d[i];

			}
			 
			idx.spt_lv[isize] = UCHAR_MAX;
			idx.spt_ld[isize] = INF_WEIGHT;
			
			NodeID larger_size = 0;
			for(NodeID i = isize; idx_original.spt_v[i] != numOfVertices; ++i){
				++larger_size;
			}
			
			
			idx.spt_v = (NodeID*)memalign(64, larger_size * sizeof(NodeID));
			idx.spt_d = (EdgeWeight*)memalign(64, larger_size * sizeof(EdgeWeight));
			
			for (NodeID i = 0; i < larger_size; ++i) {
				uint8_t hub;
				EdgeWeight hub_weight;
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;

				idx.spt_v[i] = idx_original.spt_v[i + isize];
				idx.spt_d[i] = idx_original.spt_d[i + isize];

			}			
			
			two_index_t_p &b_idx = b_two_index_p[v];
			
			index_t_p &b_idx_original = bindex_p[v];
			
			isize = 0;
			for(NodeID i = 0; b_idx_original.spt_v[i] < UCHAR_MAX; ++i){
				++isize;
			}
			

			b_idx.spt_lv = (uint8_t*)memalign(64, (isize + 1) * sizeof(uint8_t));
			b_idx.spt_ld = (EdgeWeight*)memalign(64, (isize + 1) * sizeof(EdgeWeight));

		//	index_[v].spt_v.resize(isize);
		//	index_[v].spt_d.resize(isize);

			for (NodeID i = 0; i < isize; ++i) {
				uint8_t hub;
				EdgeWeight hub_weight;
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;

				b_idx.spt_lv[i] = b_idx_original.spt_v[i];
				b_idx.spt_ld[i] = b_idx_original.spt_d[i];

			}
			 
			b_idx.spt_lv[isize] = UCHAR_MAX;
			b_idx.spt_ld[isize] = INF_WEIGHT;
			
			larger_size = 0;
			for(NodeID i = isize; b_idx_original.spt_v[i] != numOfVertices; ++i){
				++larger_size;
			}
			
			
			b_idx.spt_v = (NodeID*)memalign(64, larger_size * sizeof(NodeID));
			b_idx.spt_d = (EdgeWeight*)memalign(64, larger_size * sizeof(EdgeWeight));
			
			for (NodeID i = 0; i < larger_size; ++i) {
				uint8_t hub;
				EdgeWeight hub_weight;
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;
 
				b_idx.spt_v[i] = b_idx_original.spt_v[i + isize];
				b_idx.spt_d[i] = b_idx_original.spt_d[i + isize];

			}			
		}		
	}
	
	
	void save_labels_iteration_stats(const char* save_filename) {

		vector<NodeID> stat(numOfVertices);
		for (NodeID v = 0; v < numOfVertices; ++v) {
			for (NodeID i = 0; i < index_[v].size(); ++i)
				stat[index_[v].spt_v[i]]++;
			for (NodeID i = 0; i < bindex_[v].size(); ++i)
				stat[bindex_[v].spt_v[i]]++;
		}

		ofstream ofs(save_filename);

		for (NodeID v = 0; v < numOfVertices; ++v) {
			ofs << stat[v] << endl;
		}
		ofs.close();
	}

};

class DPLabel{
public:
	vector<index_t_path> index_;
	vector<index_t_path> bindex_; // Backward labels.

	index_t_path_p* index_p;
	index_t_path_p* bindex_p;


	DPLabel() {
		index_.resize(numOfVertices);
		bindex_.resize(numOfVertices);
	}

	~DPLabel() {
		Free();
	}

	inline EdgeWeight query_path(NodeID s, NodeID t, vector<NodeID>& rank, vector<NodeID>& inv) {



		EdgeWeight distance = INF_WEIGHT;
		NodeID meetnode = numOfVertices;
		NodeID s_parent;
		NodeID t_parent;

		const index_t_path_p &idx_s = index_p[s];
		const index_t_path_p &idx_t = bindex_p[t];

		_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_p[0], _MM_HINT_T0); 
		_mm_prefetch(&idx_t.spt_p[0], _MM_HINT_T0);

		for (int i = 0, j = 0; ; ) {
			NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

			if (v1 == numOfVertices) break;  // Sentinel

			if (v1 == v2) {
				EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
				if (td < distance) {
					distance = td;
					//if (v1 < meetnode) {
						meetnode = v1;
						s_parent = idx_s.spt_p[i];
						t_parent = idx_t.spt_p[j];
				//	}
				}
				++i;
				++j;
			}
			else {
				i += v1 < v2 ? 1 : 0;
				j += v1 > v2 ? 1 : 0;
			}
		}

		//Next, retrieve path from s - meetnode and meetnode - t.
		vector<NodeID> path_from_s;
		vector<NodeID> path_to_t;
		path_from_s.push_back(s_parent);
		path_to_t.push_back(t_parent);

		/*	if (s == 194569 && t == 20072)
		cout << "debug." << " meet: " << meetnode << " sparent:" << s_parent << " tparent:" << t_parent <<  endl;*/

		NodeID inv_meetnode = inv[meetnode];

		while (path_from_s.back() != inv_meetnode) {
			/*if (s == 194569 && t == 20072)
			cout << "s meet:" << path_from_s.back() << endl;*/
			const index_t_path_p &idx_from_s = index_p[path_from_s.back()];

			_mm_prefetch(&idx_from_s.spt_v[0], _MM_HINT_T0);
			_mm_prefetch(&idx_from_s.spt_p[0], _MM_HINT_T0);

			//	vector<NodeID>& index_from_s = index_[path_from_s.back()].spt_v;
			for (int i = 0; ; ++i) {
				if (idx_from_s.spt_v[i] == numOfVertices) break;
				if (idx_from_s.spt_v[i] == meetnode) {
					path_from_s.push_back(idx_from_s.spt_p[i]);
					break;
				}
			}
		}

		while (path_to_t.back() != inv_meetnode) {
			/*if (s == 194569 && t == 20072)
			cout << "t meet:" << path_to_t.back() << endl;*/
			//	vector<NodeID>& index_to_t = index_[path_to_t.back()].spt_v;
			const index_t_path_p &idx_to_t = bindex_p[path_to_t.back()];
			_mm_prefetch(&idx_to_t.spt_v[0], _MM_HINT_T0);
			_mm_prefetch(&idx_to_t.spt_p[0], _MM_HINT_T0);
			for (int i = 0; ; ++i) {
				if (idx_to_t.spt_v[i] == numOfVertices) break;
				if (idx_to_t.spt_v[i] == meetnode) {
					path_to_t.push_back(idx_to_t.spt_p[i]);
					break;
				}
			}
		}

		return distance;
		
	}

	EdgeWeight query_path_p(NodeID s, NodeID t, vector<NodeID>& rank, vector<NodeID>& inv) {

		EdgeWeight distance = INF_WEIGHT;
		vector<NodeID>& index_s = index_[s].spt_v;
		vector<EdgeWeight>& index_s_d = index_[s].spt_d;

		vector<NodeID>& bindex_t = bindex_[t].spt_v;
		vector<EdgeWeight>& bindex_t_d = bindex_[t].spt_d;

		
		NodeID meetnode = numOfVertices;
		int s_parent;
		int t_parent;
		for (int i = 0, j = 0; i < index_s.size(), j < bindex_t.size(); ) {
			if (index_s[i] == bindex_t[j]) {
				if (distance >(EdgeWeight)(index_s_d[i] + bindex_t_d[j])) {
					distance = (EdgeWeight)(index_s_d[i] + bindex_t_d[j]);
				//	if (index_s[i] < meetnode) {
						meetnode = index_s[i];
						s_parent = index_[s].spt_p[i];
						t_parent = index_[t].spt_p[j];
				//	}
				}
				//distance = min(distance, (EdgeWeight)(index_s_d[i] + bindex_t_d[j]));
				++i;
				++j;
			}
			else {
				if (index_s[i] < bindex_t[j])
					++i;
				else
					++j;
			}
		}

		//Next, retrieve path from s - meetnode and meetnode - t.
		vector<NodeID> path_from_s;
		vector<NodeID> path_to_t;
		path_from_s.push_back(s_parent);
		path_to_t.push_back(t_parent);

		/*	if (s == 194569 && t == 20072)
		cout << "debug." << " meet: " << meetnode << " sparent:" << s_parent << " tparent:" << t_parent <<  endl;*/

		while (path_from_s.back() != inv[meetnode]) {
			/*if (s == 194569 && t == 20072)
			cout << "s meet:" << path_from_s.back() << endl;*/
			vector<NodeID>& index_from_s = index_[path_from_s.back()].spt_v;
			for (int i = 0; i < index_from_s.size(); ++i) {
				if (index_from_s[i] == meetnode) {
					path_from_s.push_back(index_[path_from_s.back()].spt_p[i]);
					break;
				}
			}
		}

		while (path_to_t.back() != inv[meetnode]) {
			/*if (s == 194569 && t == 20072)
			cout << "t meet:" << path_to_t.back() << endl;*/
			vector<NodeID>& index_to_t = bindex_[path_to_t.back()].spt_v;
			for (int i = 0; i < index_to_t.size(); ++i) {
				if (index_to_t[i] == meetnode) {
					path_to_t.push_back(bindex_[path_to_t.back()].spt_p[i]);
					break;
				}
			}
		}

		//for (int i = 0; i < path_from_s.size(); ++i)
		//	path_from_s[i] = inv[path_from_s[i]];
		//for (int i = 0; i < path_to_t.size(); ++i)
		//	path_to_t[i] = inv[path_to_t[i]];

		return path_from_s.size() + path_to_t.size();
	}


	void Free() {
		if (index_.size() == 0 || bindex_.size() == 0) return;
		for (int v = 0; v < numOfVertices; ++v) {
			index_[v].spt_v.clear();
			index_[v].spt_d.clear();
			if (DIRECTED_FLAG == true) {
				bindex_[v].spt_v.clear();
				bindex_[v].spt_d.clear();
			}
		}
		index_.clear();
		bindex_.clear();
	}

	double avg_size() {
		double total = 0;
		for (int i = 0; i < numOfVertices; ++i) {
			total += index_[i].spt_v.size();
			total += bindex_[i].spt_v.size();
		}
		double avg = total / numOfVertices / 2 - 1; // We do not count the trivial labels (V, INF_WEIGHT).

		return avg;
	}

	void print_stat() {
		cout << "Average Label Size: " << avg_size() << endl;
		//cout << "Maximum Label Size: " << max_size() << endl;
	}

	void save_labels(const char* save_filename) {
		ofstream ofs(save_filename, ios::binary | ios::out);

		ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
		for (NodeID v = 0; v < numOfVertices; ++v) {
			int isize = index_[v].size();
			ofs.write((const char*)&isize, sizeof(isize));
			for (NodeID i = 0; i < index_[v].size(); ++i) {
				ofs.write((const char*)&index_[v].spt_v[i], sizeof(index_[v].spt_v[i]));
				ofs.write((const char*)&index_[v].spt_p[i], sizeof(index_[v].spt_p[i]));
				ofs.write((const char*)&index_[v].spt_d[i], sizeof(index_[v].spt_d[i]));
			}
			int bisize = bindex_[v].size();
			ofs.write((const char*)&bisize, sizeof(bisize));
			for (NodeID i = 0; i < bindex_[v].size(); ++i) {
				ofs.write((const char*)&bindex_[v].spt_v[i], sizeof(bindex_[v].spt_v[i]));
				ofs.write((const char*)&bindex_[v].spt_p[i], sizeof(bindex_[v].spt_p[i]));
				ofs.write((const char*)&bindex_[v].spt_d[i], sizeof(bindex_[v].spt_d[i]));
			}
		}
		ofs.close();
	}

	void load_labels(const char* load_filename) {

		index_p = NULL;
		bindex_p = NULL;

		ifstream ifs(load_filename, ios::binary | ios::in);

		NodeID isize;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;

		index_p = (index_t_path_p*)memalign(64, numOfVertices * sizeof(index_t_path_p));
		bindex_p = (index_t_path_p*)memalign(64, numOfVertices * sizeof(index_t_path_p));

		cout << numOfVertices << " vertices." << endl;

		for (NodeID v = 0; v < numOfVertices; ++v) {
			index_t_path_p &idx = index_p[v];
			ifs.read((char*)&isize, sizeof(isize));

			idx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
			idx.spt_p = (NodeID*)memalign(64, isize * sizeof(NodeID));
			idx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));

			for (NodeID i = 0; i < isize; ++i) {
				NodeID hub;
				EdgeWeight hub_weight;
				NodeID hub_parent;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_parent, sizeof(hub_parent));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;

				idx.spt_v[i] = hub;
				idx.spt_d[i] = hub_weight;
				idx.spt_p[i] = hub_parent;
			}

			//	index_[v].spt_v.resize(isize);
			//	index_[v].spt_d.resize(isize);


			index_t_path_p &bidx = bindex_p[v];
			ifs.read((char*)&isize, sizeof(isize));
			bidx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
			bidx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));
			bidx.spt_p = (NodeID*)memalign(64, isize * sizeof(NodeID));

			for (NodeID i = 0; i < isize; ++i) {
				NodeID hub;
				EdgeWeight hub_weight;
				NodeID hub_parent;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_parent, sizeof(hub_parent));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;
				bidx.spt_v[i] = hub;
				bidx.spt_d[i] = hub_weight;
				bidx.spt_p[i] = hub_parent;
			}

		}
		ifs.close();

		/*index_.clear();
		bindex_.clear();
		ifstream ifs(load_filename, ios::binary | ios::in);
		NodeID isize;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;
		index_.resize(numOfVertices);
		bindex_.resize(numOfVertices);
		for (NodeID v = 0; v < numOfVertices; ++v) {
			ifs.read((char*)&isize, sizeof(isize));
			index_[v].spt_v.resize(isize);
			index_[v].spt_p.resize(isize);
			index_[v].spt_d.resize(isize);
			for (NodeID i = 0; i < index_[v].size(); ++i) {
				NodeID hub;
				NodeID parent;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&parent, sizeof(parent));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));
				index_[v].spt_v[i] = hub;
				index_[v].spt_p[i] = parent;
				index_[v].spt_d[i] = hub_weight;
			}

			ifs.read((char*)&isize, sizeof(isize));
			bindex_[v].spt_v.resize(isize);
			bindex_[v].spt_d.resize(isize);
			for (NodeID i = 0; i < bindex_[v].size(); ++i) {
				NodeID hub;
				NodeID parent;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&parent, sizeof(parent));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));
				bindex_[v].spt_v[i] = hub;
				bindex_[v].spt_p[i] = parent;
				bindex_[v].spt_d[i] = hub_weight;
			}
		}
		ifs.close();*/
	}
	
	inline EdgeWeight query_p(NodeID s, NodeID t) {

		//EdgeWeight distance = INF_WEIGHT;
		//
		////const index_t_p &idx_s = index_p[s];
		////const index_t_p &idx_t = bindex_p[t];

		//NodeID *vs = index_p[s].spt_v;
		//NodeID *vt = bindex_p[t].spt_v;
		//EdgeWeight* ws = index_p[s].spt_d;
		//EdgeWeight* wt = bindex_p[t].spt_d;


		//_mm_prefetch(vs, _MM_HINT_T0);
		//_mm_prefetch(vt, _MM_HINT_T0);
		//_mm_prefetch(ws, _MM_HINT_T0); 
		//_mm_prefetch(wt, _MM_HINT_T0);

		//for (unsigned i = 0, j = 0; ; ) {
		//	if (*(vs + i) == *(vt + j)) {
		//		if (*(vs + i) == numOfVertices) break;  // Sentinel
		//		EdgeWeight td = *(ws + i) + *(wt + j);
		//		if (td < distance) distance = td;
		//		++i;
		//		++j;
		//	}
		//	else {
		//		i += *(vs + i) < *(vt + j) ? 1 : 0;
		//		j += *(vs + i) > *(vt + j) ? 1 : 0;
		//	}
		//}
		//return distance;


		EdgeWeight distance = INF_WEIGHT;

		const index_t_path_p &idx_s = index_p[s];
		const index_t_path_p &idx_t = bindex_p[t];


		_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);

		for (int i = 0, j = 0; ; ) {
			NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

			if (v1 == v2) {
				if (v1 == numOfVertices) break;  // Sentinel
				EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
				if (td < distance) distance = td;
				++i;
				++j;
			}
			else {
				i += v1 < v2 ? 1 : 0;
				j += v1 > v2 ? 1 : 0;
			}
		}
		return distance;
	}


};

template<int kNumBitParallelRoots = 50>
class BPLabel {

public:
	index_t_bp<kNumBitParallelRoots>* index_bp;
	
	BPLabel() {
	}

	~BPLabel() {
		//Free();
	}

	EdgeWeight query_p(NodeID s, NodeID t) {

		EdgeWeight distance = INF_WEIGHT;

		NodeID *vs = index_bp[s].spt_v;
		NodeID *vt = index_bp[t].spt_v;
		EdgeWeight* ws = index_bp[s].spt_d;
		EdgeWeight* wt = index_bp[t].spt_d;

		_mm_prefetch(vs, _MM_HINT_T0);
		_mm_prefetch(vt, _MM_HINT_T0);
		_mm_prefetch(ws, _MM_HINT_T0);
		_mm_prefetch(wt, _MM_HINT_T0);

		for (int i = 0; i < kNumBitParallelRoots; ++i) {
			EdgeWeight td = index_bp[s].bpspt_d[i] + index_bp[t].bpspt_d[i];
			if (td - 2 <= distance) {
				td +=
					(index_bp[s].bpspt_s[i][0] & index_bp[t].bpspt_s[i][0]) ? -2 :
					((index_bp[s].bpspt_s[i][0] & index_bp[t].bpspt_s[i][1]) | (index_bp[s].bpspt_s[i][1] & index_bp[t].bpspt_s[i][0]))
					? -1 : 0;

				if (td < distance) distance = td;
			}
		}

		for (unsigned i = 0, j = 0; ; ) {
			if (*(vs + i) == *(vt + j)) {
				if (*(vs + i) == numOfVertices) break;  // Sentinel
				EdgeWeight td = *(ws + i) + *(wt + j);
				if (td < distance) distance = td;
				++i;
				++j;
			}
			else {
				i += *(vs + i) < *(vt + j) ? 1 : 0;
				j += *(vs + i) > *(vt + j) ? 1 : 0;
			}
		}
		return distance;
	}
	
	EdgeWeight query_p(NodeID s, NodeID t, bool& isBP) {

		EdgeWeight distance = INF_WEIGHT;

		const index_t_bp<kNumBitParallelRoots> &idx_s = index_bp[s];
		const index_t_bp<kNumBitParallelRoots> &idx_t = index_bp[t];

		_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);

		isBP = false;

		for (int i = 0; i < kNumBitParallelRoots; ++i) {
			EdgeWeight td = index_bp[s].bpspt_d[i] + index_bp[t].bpspt_d[i];
			if (td - 2 <= distance) {
				td +=
					(index_bp[s].bpspt_s[i][0] & index_bp[t].bpspt_s[i][0]) ? -2 :
					((index_bp[s].bpspt_s[i][0] & index_bp[t].bpspt_s[i][1]) | (index_bp[s].bpspt_s[i][1] & index_bp[t].bpspt_s[i][0]))
					? -1 : 0;

				if (td < distance) {
					distance = td;
					isBP = true;
				}
			}
		}

		for (int i = 0, j = 0; ; ) {
			NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

			if (v1 == numOfVertices) break;  // Sentinel

			if (v1 == v2) {
				EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
				if (td < distance) {
					distance = td;
					isBP = false;
				}
				++i;
				++j;
			}
			else {
				i += v1 < v2 ? 1 : 0;
				j += v1 > v2 ? 1 : 0;
			}
		}
		return distance;
	}


	/*
	NodeID max_size() {
	NodeID maxsize = numeric_limits<NodeID>::min();
	for (int i = 0; i < V; ++i) maxsize = max(maxsize, index_[i].spt_v.size());
	return maxsize;
	}*/

	void print_stat() {
		cout << "Average Label Size: " << avg_size() << endl;
		//cout << "Maximum Label Size: " << max_size() << endl;
	}

	double avg_size() {
		double lab_count = 0;
		for (NodeID v = 0; v < numOfVertices; ++v) {
			NodeID isize;
			for (isize = 1; index_bp[v].spt_v[isize - 1] != numOfVertices; ++isize) continue;
			lab_count += isize;
		}
		lab_count = (double)lab_count / (double)numOfVertices - 1;
		return lab_count;
	}

	void Free() {
		for (int v = 0; v < numOfVertices; ++v) {
			free(index_bp[v].spt_v);
			free(index_bp[v].spt_d);
	
		}

		free(index_bp);
		index_bp = NULL;
	}
	
	void save_labels(const char* save_filename) {
		ofstream ofs(save_filename, ios::binary | ios::out);

		ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
		int knumbit = kNumBitParallelRoots;
		ofs.write((const char*)&knumbit, sizeof(knumbit));

		for (NodeID v = 0; v < numOfVertices; ++v) {
			index_t_bp<kNumBitParallelRoots> &idx = index_bp[v];

			for (int i = 0; i < kNumBitParallelRoots; ++i) {
				EdgeWeight d = idx.bpspt_d[i];
				uint64_t a = idx.bpspt_s[i][0];
				uint64_t b = idx.bpspt_s[i][1];
				ofs.write((const char*)&d, sizeof(d));
				ofs.write((const char*)&a, sizeof(a));
				ofs.write((const char*)&b, sizeof(b));
			}

			NodeID isize;
			for (isize = 1; idx.spt_v[isize - 1] != numOfVertices; ++isize) continue;  // Find the sentinel
			ofs.write((const char*)&isize, sizeof(isize));
			for (NodeID i = 0; i < isize; ++i) {
				ofs.write((const char*)&idx.spt_v[i], sizeof(idx.spt_v[i]));
				ofs.write((const char*)&idx.spt_d[i], sizeof(idx.spt_d[i]));
			}
		}

		ofs.close();
	}
	
	 
	void load_labels(const char* load_filename){
		
		index_bp = NULL;

		int knumbit;

		ifstream ifs(load_filename);
		NodeID isize = 0;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;		
		ifs.read((char*)&knumbit, sizeof(isize));

		if (knumbit != kNumBitParallelRoots) {
			cout << knumbit << "!=" << kNumBitParallelRoots << endl;
			return;
		}
		

		index_bp = (index_t_bp<kNumBitParallelRoots>*)memalign(64, numOfVertices * sizeof(index_t_bp<kNumBitParallelRoots>));
		

		for (NodeID v = 0; v < numOfVertices; ++v) {
			index_t_bp<kNumBitParallelRoots> &idx = index_bp[v];

			for (int i = 0; i < kNumBitParallelRoots; ++i) {

				//idx.bpspt_s[i] = (uint64_t*)memalign(64, 2 * sizeof(uint64_t));

				EdgeWeight d;
				uint64_t a, b;

				ifs.read((char*)&d, sizeof(EdgeWeight));
				ifs.read((char*)&a, sizeof(uint64_t));
				ifs.read((char*)&b, sizeof(uint64_t));
				idx.bpspt_d[i] = d;
				idx.bpspt_s[i][0] = a;
				idx.bpspt_s[i][1] = b;
			}

			ifs.read((char*)&isize, sizeof(isize));

			idx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
			idx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));


			for (NodeID i = 0; i < isize; ++i) {
				NodeID hub;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));

				idx.spt_v[i] = hub;
				idx.spt_d[i] = hub_weight;

			}
		}
		ifs.close();
	}


};

template<int kNumBitParallelRoots = 50>
class DBPLabel {

public:
	index_t_bp<kNumBitParallelRoots>* index_bp;
	index_t_bp<kNumBitParallelRoots>* bindex_bp;
	DBPLabel() {
	}

	~DBPLabel() {
	}

	/*EdgeWeight query_p(NodeID s, NodeID t) {

	EdgeWeight distance = INF_WEIGHT;

	NodeID *vs = index_p[s].spt_v;
	NodeID *vt = index_p[t].spt_v;
	EdgeWeight* ws = index_p[s].spt_d;
	EdgeWeight* wt = index_p[t].spt_d;

	_mm_prefetch(vs, _MM_HINT_T0);
	_mm_prefetch(vt, _MM_HINT_T0);
	_mm_prefetch(ws, _MM_HINT_T0);
	_mm_prefetch(wt, _MM_HINT_T0);

	for (unsigned i = 0, j = 0; ; ) {
	if (*(vs + i) == *(vt + j)) {
	if (*(vs + i) == numOfVertices) break;  // Sentinel
	EdgeWeight td = *(ws + i) + *(wt + j);
	if (td < distance) distance = td;
	++i;
	++j;
	}
	else {
	i += *(vs + i) < *(vt + j) ? 1 : 0;
	j += *(vs + i) > *(vt + j) ? 1 : 0;
	}
	}
	return distance;


	//EdgeWeight distance = INF_WEIGHT;

	//const index_t_p &idx_s = index_p[s];
	//const index_t_p &idx_t = index_p[t];

	//_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
	//_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
	//_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
	//_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);

	//for (int i = 0, j = 0; ; ) {
	//	NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

	//	if (v1 == numOfVertices) break;  // Sentinel

	//	if (v1 == v2) {
	//		EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
	//		if (td < distance) distance = td;
	//		++i;
	//		++j;
	//	}
	//	else {
	//		i += v1 < v2 ? 1 : 0;
	//		j += v1 > v2 ? 1 : 0;
	//	}
	//}
	//return distance;
	}
	*/

	/*
	NodeID max_size() {
	NodeID maxsize = numeric_limits<NodeID>::min();
	for (int i = 0; i < V; ++i) maxsize = max(maxsize, index_[i].spt_v.size());
	return maxsize;
	}*/

	EdgeWeight query_p(NodeID s, NodeID t) {

		EdgeWeight distance = INF_WEIGHT;

		NodeID *vs = index_bp[s].spt_v;
		NodeID *vt = bindex_bp[t].spt_v;
		EdgeWeight* ws = index_bp[s].spt_d;
		EdgeWeight* wt = bindex_bp[t].spt_d;

		_mm_prefetch(vs, _MM_HINT_T0);
		_mm_prefetch(vt, _MM_HINT_T0);
		_mm_prefetch(ws, _MM_HINT_T0);
		_mm_prefetch(wt, _MM_HINT_T0);

		for (int i = 0; i < kNumBitParallelRoots; ++i) {
			EdgeWeight td = index_bp[s].bpspt_d[i] + bindex_bp[t].bpspt_d[i];
			if (td - 2 <= distance) {
				td +=
					(index_bp[s].bpspt_s[i][0] & bindex_bp[t].bpspt_s[i][0]) ? -2 :
					((index_bp[s].bpspt_s[i][0] & bindex_bp[t].bpspt_s[i][1]) | (index_bp[s].bpspt_s[i][1] & bindex_bp[t].bpspt_s[i][0]))
					? -1 : 0;

				if (td < distance) distance = td;
			}
		}

		for (unsigned i = 0, j = 0; ; ) {
			if (*(vs + i) == *(vt + j)) {
				if (*(vs + i) == numOfVertices) break;  // Sentinel
				EdgeWeight td = *(ws + i) + *(wt + j);
				if (td < distance) distance = td;
				++i;
				++j;
			}
			else {
				i += *(vs + i) < *(vt + j) ? 1 : 0;
				j += *(vs + i) > *(vt + j) ? 1 : 0;
			}
		}
		return distance;
	}

	EdgeWeight query_p(NodeID s, NodeID t, bool& isBP) {
		isBP = false;
		EdgeWeight distance = INF_WEIGHT;

		NodeID *vs = index_bp[s].spt_v;
		NodeID *vt = bindex_bp[t].spt_v;
		EdgeWeight* ws = index_bp[s].spt_d;
		EdgeWeight* wt = bindex_bp[t].spt_d;

		_mm_prefetch(vs, _MM_HINT_T0);
		_mm_prefetch(vt, _MM_HINT_T0);
		_mm_prefetch(ws, _MM_HINT_T0);
		_mm_prefetch(wt, _MM_HINT_T0);

		for (int i = 0; i < kNumBitParallelRoots; ++i) {
			EdgeWeight td = index_bp[s].bpspt_d[i] + bindex_bp[t].bpspt_d[i];
			if (td - 2 <= distance) {
				td +=
					(index_bp[s].bpspt_s[i][0] & bindex_bp[t].bpspt_s[i][0]) ? -2 :
					((index_bp[s].bpspt_s[i][0] & bindex_bp[t].bpspt_s[i][1]) | (index_bp[s].bpspt_s[i][1] & bindex_bp[t].bpspt_s[i][0]))
					? -1 : 0;

				if (td < distance) {
					distance = td;
					isBP = true;
				}
			}
		}

		for (unsigned i = 0, j = 0; ; ) {
			if (*(vs + i) == *(vt + j)) {
				if (*(vs + i) == numOfVertices) break;  // Sentinel
				EdgeWeight td = *(ws + i) + *(wt + j);
				if (td < distance) {
					distance = td;
					isBP = false;
				}
				++i;
				++j;
			}
			else {
				i += *(vs + i) < *(vt + j) ? 1 : 0;
				j += *(vs + i) > *(vt + j) ? 1 : 0;
			}
		}
		return distance;
	}



	void print_stat() {
		cout << "Average Label Size: " << avg_size() << endl;
		//cout << "Maximum Label Size: " << max_size() << endl;
	}

	double avg_size() {
		double lab_count = 0;
		for (NodeID v = 0; v < numOfVertices; ++v) {
			NodeID isize;
			for (isize = 1; index_bp[v].spt_v[isize - 1] != numOfVertices; ++isize) continue;
			lab_count += isize;
			for (isize = 1; bindex_bp[v].spt_v[isize - 1] != numOfVertices; ++isize) continue;
		}
		lab_count = (double)lab_count / (double)numOfVertices - 1 / (double)2;
		return lab_count;
	}

	void Free() {
		for (int v = 0; v < numOfVertices; ++v) {
			free(index_bp[v].spt_v);
			free(index_bp[v].spt_d);
			free(index_bp[v].bpspt_d);
			free(index_bp[v].bpspt_s);
			free(bindex_bp[v].spt_v);
			free(bindex_bp[v].spt_d);
			free(bindex_bp[v].bpspt_d);
			free(bindex_bp[v].bpspt_s);
		}

		free(index_bp);
		free(bindex_bp);
		index_bp = NULL;
		bindex_bp = NULL;
	}

	void save_labels(const char* save_filename) {
		ofstream ofs(save_filename, ios::binary | ios::out);

		ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
		int knumbit = kNumBitParallelRoots;
		ofs.write((const char*)&knumbit, sizeof(knumbit));

		for (NodeID v = 0; v < numOfVertices; ++v) {
			index_t_bp<kNumBitParallelRoots> &idx = index_bp[v];
			index_t_bp<kNumBitParallelRoots> &r_idx = bindex_bp[v];

			for (int i = 0; i < kNumBitParallelRoots; ++i) {
				EdgeWeight d = idx.bpspt_d[i];
				uint64_t a = idx.bpspt_s[i][0];
				uint64_t b = idx.bpspt_s[i][1];
				ofs.write((const char*)&d, sizeof(d));
				ofs.write((const char*)&a, sizeof(a));
				ofs.write((const char*)&b, sizeof(b));
			}

			for (int i = 0; i < kNumBitParallelRoots; ++i) {
				EdgeWeight d = r_idx.bpspt_d[i];
				uint64_t a = r_idx.bpspt_s[i][0];
				uint64_t b = r_idx.bpspt_s[i][1];
				ofs.write((const char*)&d, sizeof(d));
				ofs.write((const char*)&a, sizeof(a));
				ofs.write((const char*)&b, sizeof(b));
			}

			NodeID isize;
			for (isize = 1; idx.spt_v[isize - 1] != numOfVertices; ++isize) continue;  // Find the sentinel
			ofs.write((const char*)&isize, sizeof(isize));
			for (NodeID i = 0; i < isize; ++i) {
				ofs.write((const char*)&idx.spt_v[i], sizeof(idx.spt_v[i]));
				ofs.write((const char*)&idx.spt_d[i], sizeof(idx.spt_d[i]));
			}
			for (isize = 1; r_idx.spt_v[isize - 1] != numOfVertices; ++isize) continue;  // Find the sentinel
			ofs.write((const char*)&isize, sizeof(isize));
			for (NodeID i = 0; i < isize; ++i) {
				ofs.write((const char*)&r_idx.spt_v[i], sizeof(r_idx.spt_v[i]));
				ofs.write((const char*)&r_idx.spt_d[i], sizeof(r_idx.spt_d[i]));
			}
		}

		ofs.close();
	}


	void load_labels(const char* load_filename) {

		index_bp = NULL;

		int knumbit;

		ifstream ifs(load_filename);
		NodeID isize = 0;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;
		ifs.read((char*)&knumbit, sizeof(isize));

		if (knumbit != kNumBitParallelRoots) {
			cout << knumbit << "!=" << kNumBitParallelRoots << endl;
			return;
		}


		index_bp = (index_t_bp<kNumBitParallelRoots>*)memalign(64, numOfVertices * sizeof(index_t_bp<kNumBitParallelRoots>));
		bindex_bp = (index_t_bp<kNumBitParallelRoots>*)memalign(64, numOfVertices * sizeof(index_t_bp<kNumBitParallelRoots>));


		for (NodeID v = 0; v < numOfVertices; ++v) {
			index_t_bp<kNumBitParallelRoots> &idx = index_bp[v];
			index_t_bp<kNumBitParallelRoots> &r_idx = bindex_bp[v];

			for (int i = 0; i < kNumBitParallelRoots; ++i) {

				//idx.bpspt_s[i] = (uint64_t*)memalign(64, 2 * sizeof(uint64_t));

				EdgeWeight d;
				uint64_t a, b;

				ifs.read((char*)&d, sizeof(EdgeWeight));
				ifs.read((char*)&a, sizeof(uint64_t));
				ifs.read((char*)&b, sizeof(uint64_t));

				idx.bpspt_d[i] = d;
				idx.bpspt_s[i][0] = a;
				idx.bpspt_s[i][1] = b;
			}
			for (int i = 0; i < kNumBitParallelRoots; ++i) {

				//idx.bpspt_s[i] = (uint64_t*)memalign(64, 2 * sizeof(uint64_t));

				EdgeWeight d;
				uint64_t a, b;

				ifs.read((char*)&d, sizeof(EdgeWeight));
				ifs.read((char*)&a, sizeof(uint64_t));
				ifs.read((char*)&b, sizeof(uint64_t));

				r_idx.bpspt_d[i] = d;
				r_idx.bpspt_s[i][0] = a;
				r_idx.bpspt_s[i][1] = b;
			}

			ifs.read((char*)&isize, sizeof(isize));

			idx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
			idx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));


			for (NodeID i = 0; i < isize; ++i) {
				NodeID hub;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));

				idx.spt_v[i] = hub;
				idx.spt_d[i] = hub_weight;

			}


			ifs.read((char*)&isize, sizeof(isize));

			r_idx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
			r_idx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));
			for (NodeID i = 0; i < isize; ++i) {
				NodeID hub;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));

				r_idx.spt_v[i] = hub;
				r_idx.spt_d[i] = hub_weight;
			}
		}
		ifs.close();
	}


};

#endif