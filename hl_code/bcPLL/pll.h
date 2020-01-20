#ifndef __PLL__
#define __PLL__
#include <algorithm>
#include <limits> 
#include <map>
#include <set>
#include "UsedSt.h"
#define ISEQUAL(a,b) fabs(a - b)<1e-8
using namespace std;

typedef struct valPll {
	int name;
	double f;
	valPll(int name, double f) {
		this->name = name;
		this->f = f;
	}
} valPll;

static bool lessSortPll(valPll v1, valPll v2) {
	if (v1.f != v2.f)
		return v1.f > v2.f;
	else
		return v1.name > v2.name;
}

typedef struct HeapNode {
	double d;
	int u;
	HeapNode(double d, int u){
		this->d = d;
		this->u = u;
	}
	bool operator  < (const HeapNode & rhs) const {
		return d > rhs.d;
	}
} HeapNode;

class pll {

public:
	pll(){
		maxDouble = numeric_limits<double>::max();
	};
	pll(string st){
		file = st;
		maxDouble = numeric_limits<double>::max();
	};

	void tryAddNode(int);
	void fileAssign(string);
	void graphInput();
	void graphInit();
	void labelOutput(string);
	void pruneBfs(int);
	double query(int, int);
	void prune();	//prune using the order accroding to it's degree
	void disOutput(string);	//print distance matrix
	long long getLabelSize(); //get the total size of labels
	void restart();	//clear everyting except graph
	void reinput();	//read a new graph
	void orderSort();
	//void pruneBfs(int , vector<int>&, vector<int>&, vector<int>&);	//used to prune a bfs
	void pruneBC(vector<int>&, vector<int>&);	//prune use bc
	void betweenness(int);
	void graphRDGen(string);
	void graphKGGen(string);
	int tryAddHashNodeSet(int);
	int tryAddHashNode(int);
	vector<set<int> > graphSet;
	vector<vector<hop> > graph;
	vector<vector<hop> > label;
	vector<int> newName;
	vector<int> oldName;
	int nodeNum;
	long long labelSize;
	double *centrality;
private:
	double maxDouble;
	double *dist;	
	int *fai;
	bool *visited;
	double *delta;

	string file;
	map<int, int> hashnode;	//hash node to name
	vector<double> distPrune;	//speed up vk-u distance	
	vector<int> changedList;
	vector<valPll> valList;
	vector<int> stacklist;
};

#endif
