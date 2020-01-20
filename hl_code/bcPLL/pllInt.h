#ifndef __PLLINT__
#define __PLLINT__
#include <algorithm>
#include <limits> 
#include <map>
#include "UsedSt.h"
using namespace std;

typedef struct valpllInt {
	int name;
	int f;
	valpllInt(int name, int f) {
		this->name = name;
		this->f = f;
	}
} valpllInt;

static bool lessSortpllInt(valpllInt v1, valpllInt v2) {
	if (v1.f != v2.f)
		return v1.f > v2.f;
	else
		return v1.name > v2.name;
}

typedef struct HeapNodeInt {
	int d;
	int u;
	HeapNodeInt(int d, int u){
		this->d = d;
		this->u = u;
	}
	bool operator  < (const HeapNodeInt & rhs) const {
		return d > rhs.d;
	}
} HeapNodeInt;

class pllInt {

public:
	pllInt(){
		maxDouble = numeric_limits<double>::max();
	};
	pllInt(string st){
		file = st;
		maxDouble = numeric_limits<double>::max();
	};

	void tryAddNode(int);
	void fileAssign(string);
	void fileAssign(string, string);
	void graphInput();
	void graphInput2();
	void graphInit();
	void labelOutput(string);
	void pruneBfs(int);
	int query(int, int);
	void prune();	//prune using the order accroding to it's degree
	void disOutput(string);	//print distance matrix
	long long getLabelSize(); //get the total size of labels
	void restart();	//clear everyting except graph
	void reinput();	//read a new graph
	void orderSort();
	//void pruneBfs(int , vector<int>&, vector<int>&, vector<int>&);	//used to prune a bfs
	void pruneBC(vector<int>&, vector<int>&);	//prune use bc
	void betweenness(int);
	void graphGen();
	int tryAddHashNode(int);
	vector<vector<hopInt> > graph;
	vector<vector<hopInt> > label;
	vector<int> newName;
	vector<int> oldName;
	int nodeNum;
	long long labelSize;
	double *centrality;
private:
	double maxDouble;
	int *dist;
	int *fai;
	bool *visited;
	double *delta;

	string file;
	string fileW;
	map<int, int> hashnode;	//hash node to name
	vector<int> distPrune;	//speed up vk-u distance	
	vector<int> changedList;
	vector<valpllInt> valList;
	vector<int> stacklist;
};

#endif
