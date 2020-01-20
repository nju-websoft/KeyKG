#ifndef __TEST__
#define __TEST__
#include "UsedSt.h"
using namespace std;
#define node 1000
class test{
public:
	test(string);
	string file;
	void graphName(string);
	void graphGen(){
		srand((unsigned)time(0));
		int **graph = new int *[node];
		for (int i = 0; i < node; i++){
			graph[i] = new int[node];
			for (int j = 0; j < node; j++)
				graph[i][j] = 0;
		}
		int times = node * 2;
		for (int i = 0; i < times; i++){
			int x = rand() % node;
			int y = rand() % node;
			if (x > y){
				int t = x;
				x = y;
				y = t;
			}
			graph[x][y] = graph[y][x] = rand() % 100 + 1;
		}

		ofstream f1;
		f1.open(file.c_str());
		f1 << node << endl;
		for (int i = 0; i < node; i++)
		for (int j = i + 1; j < node; j++)
		if (graph[i][j] != 0)
			f1 << i << ' ' << j << ' ' << graph[i][j] << endl;
		for (int i = 0; i < node; i++)
			delete[]graph[i];
		delete[]graph;
		f1.close();
	};
	void testBfs(){
		vector<vector<int> > graph;
		vector <int> queue;
		vector <int> dis;
		ifstream in1;
		in1.open(file.c_str());
		int nodeNum;
		in1 >> nodeNum;
		int x, y;
		for (int i = 0; i < nodeNum; i++)
			graph.push_back(vector<int>());
		while (in1 >> x >> y) {
			graph[x].push_back(y);
			graph[y].push_back(x);
		}
		in1.close();
		ofstream f1;
		f1.open((file + ".bfs").c_str());
		for (int t = 0; t < nodeNum; t++){
			queue.clear();
			dis.assign(nodeNum, -1);
			dis[t] = 0;
			int cnt = 0;
			queue.push_back(t);
			while (cnt < (int)queue.size()){
				int vk = queue[cnt];
				for (int i = 0; i < (int)graph[vk].size(); i++)
				if (dis[graph[vk][i]] == -1){
					queue.push_back(graph[vk][i]);
					dis[graph[vk][i]] = dis[vk] + 1;
				}
				cnt++;
			}
			for (int i = 0; i < (int)dis.size(); i++)
				f1 << dis[i] << ' ';
			f1 << endl;
		}
		f1.close();
	};
	void testFloyd() {
		int nodeNum;
		ifstream in1;
		in1.open(file.c_str());

		in1 >> nodeNum;

		int **graphFloyd;
		graphFloyd = new int*[nodeNum];
		for (int i = 0; i < nodeNum; i++) {
			graphFloyd[i] = new int[nodeNum];
			for (int j = 0; j < nodeNum; j++)
				graphFloyd[i][j] = 10000000;
		}
		int x, y, w;
		while (in1 >> x >> y >> w) {
			graphFloyd[x][y] = graphFloyd[y][x] = w;
		}
		for (int k = 0; k < nodeNum; k++)
		for (int i = 0; i < nodeNum; i++)
		for (int j = 0; j < nodeNum; j++)
		if (i != j)
		if (graphFloyd[i][j] > graphFloyd[i][k] + graphFloyd[k][j])
			graphFloyd[i][j] = graphFloyd[i][k] + graphFloyd[k][j];
		ofstream f2;
		f2.open((file + ".floyd").c_str());
		for (int i = 0; i < nodeNum; i++) graphFloyd[i][i] = 0;
		for (int i = 0; i < nodeNum; i++) {
			for (int j = 0; j < nodeNum; j++) {
				if (graphFloyd[i][j] == 10000000) graphFloyd[i][j] = -1;
				f2 << graphFloyd[i][j] << " ";
			}
			f2 << endl;
		}
		f2.close();

	};
	void testFile(string file1, string file2){
		ifstream f1, f2;
		f1.open(file1.c_str());
		f2.open(file2.c_str());
		double x, y, flag = 0;
		while (f1 >> x &&f2 >> y){
			if (x != y){
				flag = 1; break;
			}
		}
		if (f1 >> x) flag = 1;
		if (f2 >> y) flag = 1;
		if (flag == 1)
			cout << file1 + " " + file2 << " have difference\n";
		else
			printf("no difference\n");
		f1.close();
		f2.close();
	};
private:
};
#endif
