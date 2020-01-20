package DPBF;

import java.util.*;

public class MST {
	int parent[];
	int tag[];
	int nodeNum;
	int tagNum;
	MST(int nodeNum){
		this.nodeNum = nodeNum;
		
		parent = new int[nodeNum];
		for (int i = 0; i < nodeNum; i++) parent[i] = i;
		tagNum =-1;
		tag = new int[nodeNum];
		for (int i = 0; i < nodeNum; i++) tag[i] = 0;
	}
	
	int find(int x) {
		if (tag[x] != tagNum) {
			tag[x] = tagNum;
			parent[x] = x;
		}
		if (parent[x] != x)
			parent[x] = find(parent[x]);
		return parent[x];
	}
	
	void union(int x, int y) {
		parent[find(x)] = find(y);
	}
	
	DTreeInfo prim(ArrayList<EdgeDPBF> edges1){
		PriorityQueue<EdgeDPBF> edges = new PriorityQueue<EdgeDPBF>();
		for (EdgeDPBF it : edges1) edges.add(it);
	
		tagNum++;
		
		double ans = 0;
		ArrayList<EdgeDPBF> TreeEdges = new ArrayList<EdgeDPBF>();
		while (!edges.isEmpty()) {
			EdgeDPBF it = edges.poll();
			if (find(it.u) != find(it.v)) {
				ans += it.weight;
				TreeEdges.add(it);
				union(it.u, it.v);
			}
		}
		return new DTreeInfo(ans, TreeEdges);
	}
}
