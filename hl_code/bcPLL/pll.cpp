#include "pll.h"

void pll::tryAddNode(int x){
	x++;
	if (nodeNum < x){
		for (int i = nodeNum; i < x; i++)
			graph.push_back(vector<hop>());
		nodeNum = x;
	}
}

int pll::tryAddHashNode(int x){
	if (hashnode.find(x) == hashnode.end()){		
        int k = hashnode.size();
		hashnode[x] = k;
        if (k + 1 > nodeNum){
		    nodeNum++;
		    graph.push_back(vector<hop>());
        }
	}
	return hashnode.find(x)->second;
}

int pll::tryAddHashNodeSet(int x){
	if (hashnode.find(x) == hashnode.end()){		
        int k = hashnode.size();
		hashnode[x] = k;
        if (k + 1 > nodeNum){
		    nodeNum++;
		    graphSet.push_back(set<int>());
        }
	}
	return hashnode.find(x)->second;
}
void pll::graphInput(){
	int x, y;
	double w;
	ifstream in1;
	in1.open(file.c_str());
	nodeNum = 0;
	for (int i = 0; i < nodeNum; i++)
		graph.push_back(vector<hop>());	
	while (in1 >> x >> y >> w) {
		tryAddNode(x);
		tryAddNode(y);
		graph[x].push_back(hop(y, w));
		graph[y].push_back(hop(x, w));
	}
	in1.close();
}

void pll::graphRDGen(string outfile){
	hashnode.clear();
	int x, y;
	srand((unsigned)time(0));
	ifstream in1;
	in1.open(file.c_str());
    
	nodeNum = 0;
	in1 >> nodeNum;
	for (int i = 0; i < nodeNum; i++)
		graphSet.push_back(set<int>());
	while (in1 >> x >> y) {		
		x = tryAddHashNodeSet(x);
		y = tryAddHashNodeSet(y);
		graphSet[x].insert(y);
		graphSet[y].insert(x);
	}
	cout << "input end" << endl;
	in1.close();
	visited = new bool[nodeNum];
	for (int i = 0; i < nodeNum; i++)visited[i] = false;
	int loc = 0;
	for (int i = 1; i < (int)graphSet.size(); i++)
	if (graphSet[i].size() > graphSet[loc].size())
		loc = i;

	//bfs from node with biggest degree
	visited[loc] = true;
	queue<int> ll;
	ll.push(loc);
	while (!ll.empty()){
		loc = ll.front();
        for (int it : graphSet[loc]){
    		if (!visited[it]){
	    		ll.push(it);
		    	visited[it] = true;
    		}
        }
		ll.pop();
	}
	int nodecnt = 0;
	int *linkName = new int[nodeNum];
	for (int i = 0; i < nodeNum; i++)
	if (visited[i]){
		linkName[i] = nodecnt;
		nodecnt++;
	}
	cout << "bfs end!" << endl;

	ofstream out1;
	out1.open((outfile).c_str());
	for (int i = 0; i < (int)graphSet.size(); i++){
		if (!visited[i]) continue;
        for (int it : graphSet[i])
		    if (visited[it] && linkName[i] < linkName[it])
			out1 << linkName[i]<<' '<<linkName[it]<<' '<<(rand()%100 + 1)<< endl;
			//out1 << linkName[i]<<' '<<linkName[it]<<' '<< 1<< endl;
	}
	out1.close();
	delete(linkName);
	delete(visited);
	hashnode.clear();
}

void pll::graphKGGen(string outfile){
	hashnode.clear();
	int x, y;
	srand((unsigned)time(0));
	ifstream in1;
	in1.open(file.c_str());

   	nodeNum = 0;
	in1 >> nodeNum;
	for (int i = 0; i < nodeNum; i++)
		graphSet.push_back(set<int>());
	while (in1 >> x >> y) {		
		x = tryAddHashNodeSet(x);
		y = tryAddHashNodeSet(y);
		graphSet[x].insert(y);
		graphSet[y].insert(x);
	}
	cout << "input end" << endl;
	in1.close();
	visited = new bool[nodeNum];
	for (int i = 0; i < nodeNum; i++)visited[i] = false;
	int loc = 0;
	for (int i = 1; i < (int)graphSet.size(); i++)
	if (graphSet[i].size() > graphSet[loc].size())
		loc = i;

	//bfs from node with biggest degree
	visited[loc] = true;
	queue<int> ll;
	ll.push(loc);
	while (!ll.empty()){
		loc = ll.front();
        for (int it : graphSet[loc]){
    		if (!visited[it]){
	    		ll.push(it);
		    	visited[it] = true;
    		}
        }
		ll.pop();
	}
	int nodecnt = 0;
	int *linkName = new int[nodeNum];
	for (int i = 0; i < nodeNum; i++)
	if (visited[i]){
		linkName[i] = nodecnt;
		nodecnt++;
	}
	cout << "bfs end!" << endl;

    int *edgeNum = new int[nodeNum];
    for (int i = 0; i < nodeNum; i++) edgeNum[i] = 0;
	for (int i = 0; i < (int)graphSet.size(); i++){
		if (!visited[i]) continue;
        for (int it : graphSet[i])
		    if (visited[it]) edgeNum[i]++;
	}

	ofstream out1;
	out1.open((outfile).c_str());
	for (int i = 0; i < (int)graphSet.size(); i++){
		if (!visited[i]) continue;
        for (int it : graphSet[i])
		    if (visited[it] && linkName[i] < linkName[it]){
                double nn = edgeNum[i] + edgeNum[it];
                //if (nn < edgeNum[it.goal]) nn = edgeNum[it.goal];
                //nn = log(1+it.goal)/log(2);
                //int nn1 = (int)(log(nn)/log(2)*100);
                int nn1 = round(log(edgeNum[i])/log(2)+log(edgeNum[it])/log(2));
			    out1 << linkName[i]<<' '<<linkName[it]<<' '<< nn1<< endl;
            }
	}
	out1.close();
	delete(linkName);
	delete(visited);
	hashnode.clear();
}

void pll::graphInit() {
	graphInput();
	stacklist.clear();
	oldName.assign(nodeNum, 0);
	newName.assign(nodeNum, -1);
	label.assign(nodeNum, vector<hop>());
	distPrune.assign(nodeNum, -1);	
	dist = new double[nodeNum];
	centrality = new double[nodeNum];
	fai = new int[nodeNum];
	visited = new bool[nodeNum];
	delta = new double[nodeNum];
	labelSize = 0;
	for (int i = 0; i < nodeNum; i++) centrality[i] = 0;
}


void pll::fileAssign(string st) {
	file = st;
}

void pll::restart() {
	oldName.clear();
	newName.clear();
	stacklist.clear();
	for (int i = 0; i < nodeNum; i++)
		label[i].clear();
	label.clear();
	distPrune.clear();	
	delete(dist);	
	delete(fai);
	delete(centrality);
	delete(visited);
	delete(delta);

	oldName.assign(nodeNum, 0);
	newName.assign(nodeNum, -1);
	label.assign(nodeNum, vector<hop>());
	distPrune.assign(nodeNum, -1);	
	dist = new double[nodeNum];	
	fai = new int[nodeNum];
	centrality = new double[nodeNum];
	visited = new bool[nodeNum];
	delta = new double[nodeNum];
	labelSize = 0;
	for (int i = 0; i < nodeNum; i++) centrality[i] = 0;
}

void pll::reinput() {
	int x, y;
	double w;
	ifstream in1;
	in1.open(file.c_str());

	in1 >> nodeNum;
	while (in1 >> x >> y >> w) {
		tryAddNode(x);
		tryAddNode(y);
		graph[x].push_back(hop(y, w));
		graph[y].push_back(hop(x, w));
	}
	in1.close();
}

void pll::orderSort() {
	valList.assign(nodeNum, valPll(0, 0));
	for (int i = 0; i < nodeNum; i++) {
		valList[i].name = i;
		valList[i].f = graph[i].size();
	}
	sort(valList.begin(), valList.end(), lessSortPll);
	for (int i = 0; i < nodeNum; i++)
		newName[valList[i].name] = i;
	for (int i = 0; i < nodeNum; i++)
		oldName[i] = valList[i].name;
}

void pll::prune() {
	//ofstream f1;
	//f1.open((file + ".size").c_str());
	for (int i = 0; i < nodeNum; i++) {
		pruneBfs(oldName[i]);
		if (i % 100 == 0)
			printf("%d has been pruned: %d\n", i, labelSize);
		//f1 << i << " " << labelSize << endl;
	}
	//f1.close();
}

long long pll::getLabelSize() {
	return labelSize;
}

void pll::pruneBC(vector<int>& newName, vector<int>& oldName){
	this->oldName = oldName;
	this->newName = newName;
	/*
    delete(fai);	
	delete(delta);
	delete(centrality);
    */
    int cnt = 0;
	for (int i = 0; i < nodeNum; i++) dist[i] = -1;
	for (int i = 0; i < nodeNum; i++) visited[i] = false;	
	for (int i = 0; i < nodeNum; i++) {
		pruneBfs(oldName[i]);
		if (labelSize / nodeNum > cnt){
            cnt = labelSize / nodeNum;
            if (cnt % 100 == 0)
    			printf("%d has been pruned: %d\n", i, cnt);	
        }
	}
}
void pll::pruneBfs(int u) {
	int newU = newName[u];
	for (hop it : label[u]) distPrune[it.goal] = it.len;

	priority_queue<HeapNode> q; /// -æ‡¿Î,µ„
	q.push(HeapNode{ 0, u });
	dist[u] = 0;
	changedList.push_back(u);
	while (!q.empty()) {
		int vk = q.top().u;
		q.pop();
		if (visited[vk]) continue;
		visited[vk] = true;		
		double query = maxDouble;		
		//get distance u-vk using 2 hop cover
		for (hop it : label[vk]) {
			if (distPrune[it.goal] > 0 && distPrune[it.goal] + it.len < query)
				query = distPrune[it.goal] + it.len;
		}
		//if 2 hop cover u-vk is not shortest u-vk
		if (query <= dist[vk]) continue;

		label[vk].push_back(hop(newU, dist[vk]));
		labelSize++;
		for (hop it : graph[vk]){
			if (visited[it.goal]) continue;
			if (dist[it.goal] < 0 || (dist[it.goal] > dist[vk] + it.len)) {
				changedList.push_back(it.goal);
				dist[it.goal] = dist[vk] + it.len;
				q.push(HeapNode{ dist[it.goal], it.goal });
			}
		}
	}
	//withdraw distance from u to other node
	for (hop it : label[u]) distPrune[it.goal] = -1;
	for (int i : changedList){
		dist[i] = -1;
		visited[i] = false;
	}
	changedList.clear();
}

void pll::betweenness(int s){
	for (int i = 0; i < nodeNum; i++) fai[i] = 0;
	for (int i = 0; i < nodeNum; i++) delta[i] = 0;	
	for (int i = 0; i < nodeNum; i++) dist[i] = -1;
	for (int i = 0; i < nodeNum; i++) visited[i] = false;
	fai[s] = 1;
	dist[s] = 0;
	
	priority_queue<HeapNode> q; /// -æ‡¿Î,µ„
	q.push(HeapNode{ 0, s });

	while (!q.empty()){
		HeapNode tp = q.top(); q.pop();
		int u = tp.u;
		if (visited[u]) continue;		
		visited[u] = true;
		//cout << u <<endl;
		stacklist.push_back(u);
		for (hop edge : graph[u]){
			int v = edge.goal;
			double len = edge.len;
			if (visited[v]){
				if (dist[u] >= dist[v] + len)
					fai[u] += fai[v];
				continue;
			}
			if (dist[v] < 0 || dist[v] > dist[u] + len){
				dist[v] = dist[u] + len;
				q.push(HeapNode{ dist[v], v });
			}			
		}
	}	
		
	int order = stacklist.size();
	while (order > 1) {
		order--;
		int w = stacklist[order];
		//cout << w << " ";
		for (hop edge : graph[w]) {
			int v = edge.goal;
			if (ISEQUAL(dist[v] + edge.len, dist[w]))
				delta[v] += ((double)fai[v]) / fai[w] * (1 + delta[w]);
		}
		centrality[w] = centrality[w] + delta[w];
	}	
	stacklist.clear();
}

double pll::query(int s, int t) {
	double ans = maxDouble;
	int is = 0;
	int it = 0;
	while (is < (int)label[s].size() && it < (int)label[t].size()) {
		if (label[s][is].goal > label[t][it].goal) {
			it++;
			continue;
		}
		if (label[s][is].goal < label[t][it].goal) {
			is++;
			continue;
		}
		if (label[s][is].goal == label[t][it].goal) {
			if (ans > label[s][is].len + label[t][it].len)
				ans = label[s][is].len + label[t][it].len;
			is++;
			it++;
		}
	}
	return ans;
}

void pll::labelOutput(string st) {
	ofstream f1;
	f1.open((file + "." + st).c_str());
	for (int i = 0; i < nodeNum; i++)
		for (int j = 0; j < (int)label[i].size(); j++)
			f1 << newName[i] << label[i][j].goal << ' ' << label[i][j].len << endl;
}

void pll::disOutput(string st) {
	ofstream f1;
	f1.open((file + "." + st).c_str());
	for (int i = 0; i < nodeNum; i++) {
		for (int j = 0; j < nodeNum; j++) {
			double t = query(i, j);
			if (t >= maxDouble - 1)
				t = -1;
			f1 << t << " ";
		}
		f1 << endl;
	}
}
