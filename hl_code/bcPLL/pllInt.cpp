#include "pllInt.h"

void pllInt::tryAddNode(int x){
	x++;
	if (nodeNum < x){
		for (int i = nodeNum; i < x; i++)
			graph.push_back(vector<hopInt>());
		nodeNum = x;
	}
}

int pllInt::tryAddHashNode(int x){
	if (hashnode.find(x) == hashnode.end()){
		int k = (int)hashnode.size();
        hashnode[x] = k;
		if (k + 1 > nodeNum){
			nodeNum++;
			graph.push_back(vector<hopInt>());
			cout << x << " " << nodeNum << endl;
		}
	}
	return hashnode.find(x)->second;
}
void pllInt::graphInput(){
	int x, y;
	int w = 1;
	ifstream in1;
	in1.open(file.c_str());
	nodeNum = 0;
	//in1 >> nodeNum;
	for (int i = 0; i < nodeNum; i++)
		graph.push_back(vector<hopInt>());
	//while (in1 >> x >> y) {
	while (in1 >> x >> y >> w) {
		tryAddNode(x);
		tryAddNode(y);
		graph[x].push_back(hopInt(y, w));
		graph[y].push_back(hopInt(x, w));
	}
	in1.close();
}

void pllInt::graphGen(){
	hashnode.clear();
	int x, y, w;
	srand((unsigned)time(0));
	ifstream in1;
	in1.open(file.c_str());

	nodeNum = 0;
	in1 >> nodeNum;
	for (int i = 0; i < nodeNum; i++)
		graph.push_back(vector<hopInt>());
	while (in1 >> x >> y) {
		x = tryAddHashNode(x);
		y = tryAddHashNode(y);
		w = rand() % 100 + 1;
		graph[x].push_back(hopInt(y, w));
		graph[y].push_back(hopInt(x, w));
	}
	cout << "input end" << endl;
	cout << nodeNum << endl;
	in1.close();
	visited = new bool[nodeNum];
	for (int i = 0; i < nodeNum; i++)visited[i] = false;
	int loc = 0;
	for (int i = 1; i < (int)graph.size(); i++)
	if (graph[i].size() > graph[loc].size())
		loc = i;

	//bfs from node with biggest degree
	visited[loc] = true;
	queue<int> ll;
	ll.push(loc);
	while (!ll.empty()){
		loc = ll.front();
		for (hopInt it : graph[loc])
		if (!visited[it.goal]){
			ll.push(it.goal);
			visited[it.goal] = true;
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
	cout << "bfs end" << endl;
	ofstream out1;
	out1.open((file + ".weight").c_str());
	//out1 << nodecnt << endl;
	for (int i = 0; i < (int)graph.size(); i++){
		if (!visited[i]) continue;
		for (hopInt it : graph[i])
		if (visited[it.goal] && linkName[i] < linkName[it.goal])
			out1 << linkName[i] << ' ' << linkName[it.goal] << ' ' << it.len << endl;
	}
	out1.close();
	delete(linkName);
	delete(visited);
	hashnode.clear();
}

void pllInt::graphInit() {
	graphInput();
	//graphInput2();
	stacklist.clear();
	oldName.assign(nodeNum, 0);
	newName.assign(nodeNum, -1);
	label.assign(nodeNum, vector<hopInt>());
	distPrune.assign(nodeNum, -1);
	dist = new int[nodeNum];
	centrality = new double[nodeNum];
	fai = new int[nodeNum];
	visited = new bool[nodeNum];
	delta = new double[nodeNum];
	labelSize = 0;
	for (int i = 0; i < nodeNum; i++) centrality[i] = 0;
}


void pllInt::fileAssign(string st) {
	file = st;
}

void pllInt::fileAssign(string st, string st2) {
	file = st;
	fileW = st2;
}
void pllInt::restart() {
	oldName.clear();
	newName.clear();
	stacklist.clear();
	for (int i = 0; i < nodeNum; i++)
		vector<hopInt>().swap(label[i]);
	//label[i].clear();
	label.clear();
	distPrune.clear();
	delete(dist);
	delete(fai);
	delete(centrality);
	delete(delta);
	delete(visited);

	oldName.assign(nodeNum, 0);
	newName.assign(nodeNum, -1);
	label.assign(nodeNum, vector<hopInt>());
	distPrune.assign(nodeNum, -1);
	dist = new int[nodeNum];
	fai = new int[nodeNum];
	centrality = new double[nodeNum];
	visited = new bool[nodeNum];
	delta = new double[nodeNum];
	labelSize = 0;
	for (int i = 0; i < nodeNum; i++) centrality[i] = 0;
}

void pllInt::reinput() {
	int x, y, w;
	ifstream in1;
	in1.open(file.c_str());

	in1 >> nodeNum;
	while (in1 >> x >> y >> w) {
		tryAddNode(x);
		tryAddNode(y);
		graph[x].push_back(hopInt(y, w));
		graph[y].push_back(hopInt(x, w));
	}
	in1.close();
}

void pllInt::orderSort() {
	valList.assign(nodeNum, valpllInt(0, 0));
	for (int i = 0; i < nodeNum; i++) {
		valList[i].name = i;
		valList[i].f = (int)graph[i].size();
	}
	sort(valList.begin(), valList.end(), lessSortpllInt);
	for (int i = 0; i < nodeNum; i++)
		newName[valList[i].name] = i;
	for (int i = 0; i < nodeNum; i++)
		oldName[i] = valList[i].name;
}

void pllInt::prune() {
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

long long pllInt::getLabelSize() {
	return labelSize;
}

void pllInt::pruneBC(vector<int>& newName, vector<int>& oldName){
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
			cnt = (int)labelSize / nodeNum;
			if (cnt % 100 == 0)
				printf("%d has been pruned: %d\n", i, cnt);
		}
	}
}
void pllInt::pruneBfs(int u) {
	int newU = newName[u];
	for (hopInt it : label[u]) distPrune[it.goal] = it.len;

	priority_queue<HeapNodeInt> q; /// -æ‡¿Î,µ„
	q.push(HeapNodeInt{ 0, u });
	dist[u] = 0;
	changedList.push_back(u);
	while (!q.empty()) {
		int vk = q.top().u;
		q.pop();
		if (visited[vk]) continue;
		visited[vk] = true;
		double query = maxDouble;
		//get distance u-vk using 2 hop cover
		for (hopInt it : label[vk]) {
			if (distPrune[it.goal] > 0 && distPrune[it.goal] + it.len < query)
				query = distPrune[it.goal] + it.len;
		}
		//if 2 hop cover u-vk is not shortest u-vk
		if (query <= dist[vk]) continue;

		label[vk].push_back(hopInt(newU, dist[vk]));
		labelSize++;
		for (hopInt it : graph[vk]){
			if (visited[it.goal]) continue;
			if (dist[it.goal] < 0 || (dist[it.goal] > dist[vk] + it.len)) {
				changedList.push_back(it.goal);
				dist[it.goal] = dist[vk] + it.len;
				q.push(HeapNodeInt{ dist[it.goal], it.goal });
			}
		}
	}
	//withdraw distance from u to other node
	for (hopInt it : label[u]) distPrune[it.goal] = -1;
	for (int i : changedList){
		dist[i] = -1;
		visited[i] = false;
	}
	changedList.clear();
}

void pllInt::betweenness(int s){
	for (int i = 0; i < nodeNum; i++) fai[i] = 0;
	for (int i = 0; i < nodeNum; i++) delta[i] = 0;
	for (int i = 0; i < nodeNum; i++) dist[i] = -1;
	for (int i = 0; i < nodeNum; i++) visited[i] = false;
	fai[s] = 1;
	dist[s] = 0;

	priority_queue<HeapNodeInt> q; /// -æ‡¿Î,µ„
	q.push(HeapNodeInt{ 0, s });

	while (!q.empty()){
		HeapNodeInt tp = q.top(); q.pop();
		int u = tp.u;
		if (visited[u]) continue;
		visited[u] = true;
		//cout << u <<endl;
		stacklist.push_back(u);
		for (hopInt edge : graph[u]){
			int v = edge.goal;
			int len = edge.len;
			if (visited[v]){
				if (dist[u] >= dist[v] + len)
					fai[u] += fai[v];
				continue;
			}
			if (dist[v] < 0 || dist[v] > dist[u] + len){
				dist[v] = dist[u] + len;
				q.push(HeapNodeInt{ dist[v], v });
			}
		}
	}

	int order = (int)stacklist.size();
	while (order > 1) {
		order--;
		int w = stacklist[order];
		//cout << w << " ";
		for (hopInt edge : graph[w]) {
			int v = edge.goal;
			if (dist[v] + edge.len == dist[w])
				delta[v] += ((double)fai[v]) / fai[w] * (1 + delta[w]);
		}
		centrality[w] = centrality[w] + delta[w];
	}
	stacklist.clear();
}

int pllInt::query(int s, int t) {
	int ans = maxDouble;
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

void pllInt::labelOutput(string st) {
	ofstream f1;
	f1.open((file + "." + st).c_str());
	for (int i = 0; i < nodeNum; i++)
	for (int j = 0; j < (int)label[i].size(); j++)
		f1 << newName[i] << label[i][j].goal << ' ' << label[i][j].len << endl;
}

void pllInt::disOutput(string st) {
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
