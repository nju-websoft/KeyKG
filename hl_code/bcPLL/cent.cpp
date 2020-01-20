#include "cent.h"

cent::cent() {
}

void cent::chlV1(int runcase) {	
	int topNum = runcase;
	if (topNum > nodeNum) topNum = nodeNum;
	//get the order according to it's degree
	valList.assign(nodeNum, valCent(0, 0));
	for (int i = 0; i < nodeNum; i++) {
		valList[i].name = i;
		valList[i].f = (double)p1.graph[i].size();
	}
	sort(valList.begin(), valList.end(), lessSortCent);
	for (int i = 0; i < topNum; i++){
		p1.betweenness(valList[i].name);
		if (i % 50 == 49)
			printf("betweenness %d has done\n", i);
	}
	for (int i = 0; i < nodeNum; i++) {
		valList[i].name = i;
		valList[i].f = p1.centrality[i];
	}	
	sort(valList.begin(), valList.end(), lessSortCent);	
	
	//for (int i = 0; i < 100; i++)
		//cout << valList[i].name << " " << valList[i].f << endl;

	vector<int> newName(nodeNum, 0);
	vector<int> oldName(nodeNum, 0);
	for (int i = 0; i < nodeNum; i++) {
		newName[valList[i].name] = i;
		oldName[i] = valList[i].name;
	}
	p1.pruneBC(newName, oldName);
}

void cent::chlV2(int runcase) {
    int topNum = runcase;
    if (topNum > nodeNum) topNum = nodeNum;

    vector<int> v(nodeNum);
    for (int i = 0; i<v.size(); ++i) v[i] = i;
    random_device rd;
    mt19937 g(rd());
    shuffle(v.begin(), v.end(), g);

    srand(time(NULL));
    for (int i = 0; i < topNum; i++){
        p1.betweenness(v[i]);
        if (i % 50 == 49)
            printf("betweenness %d has done\n", i);
    }
    valList.assign(nodeNum, valCent(0, 0));
    for (int i = 0; i < nodeNum; i++) {
        valList[i].name = i;
        valList[i].f = p1.centrality[i];
    }
    sort(valList.begin(), valList.end(), lessSortCent);

    vector<int> newName(nodeNum, 0);
    vector<int> oldName(nodeNum, 0);
    for (int i = 0; i < nodeNum; i++) {
        newName[valList[i].name] = i;
        oldName[i] = valList[i].name;
    }
    p1.pruneBC(newName, oldName);
}
void cent::deal(string st, bool degree, int runcase = 200) {
	file = st;
	p1.fileAssign(file);
	p1.graphInit();
	nodeNum = p1.nodeNum;

	clock_t b3;
	b3 = clock();
	if (degree) chlV1(runcase); else chlV2(runcase);
	long long size = p1.getLabelSize();
	cout <<"running time:"<< (double)(clock() - b3)/CLOCKS_PER_SEC << "s" << endl;
	printf("chl label size %lld\n", size);
	printf("chl ave label size %f\n", (double)size/nodeNum);
	//p1.disOutput("pl");
}

void cent::deal(string st, string st2, int runcase = 200) {
	file = st;
	p1.graphInit();
	nodeNum = p1.nodeNum;

	clock_t b3;
	b3 = clock();
	chlV1(runcase);
	long long size = p1.getLabelSize();
	cout << clock() - b3 << "ms" << endl;
	printf("chl label size %lld\n", size);

}
