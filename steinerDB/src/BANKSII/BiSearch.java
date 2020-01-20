package BANKSII;

import java.util.*;

import DOGST.Hop;


public class BiSearch {
	int dmax = 8;
	double dist[][];
	int sp[][];
	boolean visited[];	//used for remove duplicate visit
	int realName[];
	int nodeNum, keyNum;
	public PriorityQueue<AnsTree> ansQueue;	//record ans tree
	public ArrayList<AnsTree> ansList;
	int topk = 1;
	double mu = 0.5;
	public double totalweight = 0;
	double [][]at;
	double []atv;
	int []depth;
	MaxHeap Qin, Qout;
	GraphInfo g1;
	void para_init() {
		nodeNum = g1.nodeNum;
		keyNum = g1.keyNum;
		visited = new boolean[nodeNum];
		ansQueue = new PriorityQueue<AnsTree>();
		ansList = new ArrayList<AnsTree>();
		for (int i = 0; i < nodeNum; i++)
			visited[i] = false;
		
		Qin = new MaxHeap();
		Qin.heapInit(g1.nodeNum, g1.keyNum);
		
		Qout = new MaxHeap();
		Qout.heapInit(g1.nodeNum, g1.keyNum);
		
		dist = new double[nodeNum][keyNum];
		for (int i = 0; i < nodeNum; i++)
			for (int j = 0; j < keyNum; j++)
				dist[i][j] = -1;
		sp = new int[nodeNum][keyNum];
		for (int i = 0; i < nodeNum; i++)
			for (int j = 0; j < keyNum; j++)
				sp[i][j] = -1;
		
		at = new double[nodeNum][keyNum];
		for (int i = 0; i < nodeNum; i++)
			for (int j = 0; j < keyNum; j++)
				at[i][j] = 0;
		
		atv = new double[nodeNum];
		for (int i = 0; i < nodeNum; i++)
			atv[i] = 0;
		
		depth = new int[nodeNum];
		for (int i = 0; i < nodeNum; i++)
			depth[i] = -1;
		
		realName = new int[keyNum];
		for (int i = 0; i < keyNum; i++)
			realName[i] = g1.HashKeywordTable.get(g1.keywordList.get(i));			
		
		for (int i = 0; i < keyNum; i++) {
			ArrayList<Integer> Si=  g1.invTable.get(realName[i]);
			for (int it : Si) {
				at[it][i] = 1.0 / Si.size() * (1 - mu);
				atv[it] += at[it][i];
				dist[it][i] = 0;
				sp[it][i] = it;
				depth[it] = 0;
				Qin.tryHigher(it, atv[it]);
			}
		}
	}
	double bidir_exp_search(){		
		para_init();		
		while (true) {
			//get node with highest activation
			double maxatv = -1;
			int goal = -1;
			if (!Qin.isempty()) {				
				if (maxatv < Qin.topatv()) {
					maxatv = Qin.topatv();
					goal = 1;
				}
			}
			if (!Qout.isempty()) {				
				if (maxatv < Qout.topatv()) {
					maxatv = Qout.topatv();
					goal = 2;
				}
			}
			if (maxatv < 0)
				break;
			/*
			cnt++;
			if (cnt % 100 == 0)
				System.out.println(cnt + " has been poped");
			*/
			switch (goal) {
			case 1:
				extend(Qin, true, Qout);
				break;
			case 2:
				extend(Qout, false, null);
				break;
			default:
				System.out.println("Nope, find no node to extend!");
				break;
			}
							
		}
		
		while (!ansQueue.isEmpty())
			ansList.add(ansQueue.poll());
		for (AnsTree it : ansList) {
			totalweight += it.weight;
			//System.out.println(it.weight);
		}		
		double keynode = 0;
		for (int i = 1; i < keyNum; i++)	//record all 
			keynode = keynode + g1.invTable.get(realName[i]).size();
		return keynode;		
	}
	
	//merge line 7-14 and 16-23 in one function
	void extend(MaxHeap q1, boolean isIn, MaxHeap q2) {
		//Pop best v from Q and insert in X
		HeapNode topnode = q1.pop();
		//System.out.println("extend node : "+topnode.v + " atv: " + topnode.activation);
		//if is-Complete(u) then EMIT(u)
		if (is_complete(topnode.v))	emit(topnode.v);
		//if depth u < d max then
		if (depth[topnode.v] < dmax) {
			//∀sh.v ∈ outgoing[topnode.v]
			for (Hop sh : g1.graph.get(topnode.v)) {
				//ExploreEdge(u,v)
				exploreEdge(sh.v, topnode.v, sh.dis);
				//if v /∈  X
				if (depth[sh.v] < 0)
					depth[sh.v] = depth[topnode.v] + 1;				
				q1.tryHigher(sh.v, atv[sh.v]);
			}
			if (!isIn) return;
			if (!q2.inX(topnode.v))
				q2.tryHigher(topnode.v, atv[topnode.v]);
		}
	}
	
	void exploreEdge(int u, int v, double w) {
		//for each keyword i
		for (int i = 0; i < keyNum; i++) {
			//if u has a better path to t i via v then
			if (dist[v][i] >= 0 && (dist[u][i] < 0 || dist[u][i] > dist[v][i] + w)) {
				//sp u,i ← v
				sp[u][i] = v;
				//update dist u,i with this new dist
				dist[u][i] = dist[v][i] + w;
				attach(u,i);
				//if is-Complete(u) then EMIT(u)
				if (is_complete(u)) emit(u);	
			}
			//if v spreads more activation to u from t i then
			if (at[u][i]/ (1-mu) < at[v][i]/(1-mu) * mu *(1/w) / g1.nodeInWeight.get(v)) {
				//update a u,i with this new activation
				atv[u] = atv[u] - at[u][i];
				at[u][i] = at[v][i]/(1-mu) * mu *(1/w) / g1.nodeInWeight.get(v) * (1-mu);
				atv[u] = atv[u] + at[u][i];
				
				activate(u, i);	//Activate(u,i)
			}							
		}
	}
		
	void attach(int v, int k) {
		//update priority of v if it is present in Qin
		Qin.tryHigher(v, atv[k]);
		
		//dijkstra to propagate distance
		ArrayList<Integer> list = new ArrayList<Integer>();		
		PriorityQueue<Hop> nodes = new PriorityQueue<Hop>();		
		nodes.add(new Hop(v, dist[v][k]));		
		while (!nodes.isEmpty()) {
			Hop sh = nodes.poll();
			if (visited[sh.v])
				continue;
			list.add(sh.v);
			visited[sh.v] = true;
			for (Hop edge : g1.graph.get(sh.v)) {
				if (depth[edge.v] < 0)	//withdraw none reached ancestor
					continue;
				if (dist[sh.v][k] >= 0 && (dist[edge.v][k] < 0 || edge.dis + dist[sh.v][k] < dist[edge.v][k])) {
					dist[edge.v][k] = edge.dis + dist[sh.v][k];
					sp[edge.v][k] = sh.v; 
					nodes.add(new Hop(edge.v, dist[edge.v][k]));
				}
			}
		}
		for (int i : list) visited[i] = false;
		list.clear();
	}
	
	void activate(int v, int k) {
		//update priority of v if it is present in Qin
		Qin.tryHigher(v, atv[k]);
		
		ArrayList<Integer> list = new ArrayList<Integer>();		
		PriorityQueue<Hop> nodes = new PriorityQueue<Hop>();		
		nodes.add(new Hop(v, -at[v][k]/(1-mu)*mu));		
		while (!nodes.isEmpty()) {
			Hop sh = nodes.poll();	//node with biggest atv to spread
			if (visited[sh.v])
				continue;
			list.add(sh.v);
			atv[sh.v] = atv[sh.v] - at[sh.v][k];  
			visited[sh.v] = true;
			
			for (Hop edge : g1.graph.get(sh.v)) {				
				if (depth[edge.v] < 0 || dist[sh.v][k] <= dist[edge.v][k]) continue;
				double newat = -sh.dis* (1/edge.dis) / g1.nodeInWeight.get(sh.v);
				if (newat*(1-mu) > at[edge.v][k]) {
					at[edge.v][k] = newat*(1-mu);
					nodes.add(new Hop(edge.v, newat*mu));
				}					
			}
		}
		for (int i : list) {
			visited[i] = false;
			atv[i] = atv[i] + at[i][k];
			Qin.tryHigher(i, atv[i]);
		};
		list.clear();
	}
	
	void emit(int v) {
		//System.out.println("emit: "+ v);
		AnsTree anst = new AnsTree();
		ArrayList<Integer> path = new ArrayList<Integer>();
		visited[v] = true;
		//generate tree
		for (int i = 0; i < keyNum; i++) {
			path.add(v);
			int tv = v;
			while (sp[tv][i] != tv) {
				tv = sp[tv][i];
				path.add(tv);
			}
			
			//merge path to tree
			for (int j = path.size() - 1; j >= 0; j--) {
				if (visited[path.get(j)]) {
					anst.weight += dist[path.get(j)][i];	//add path length to tree
					break;
				}
				visited[path.get(j)] = true;
				anst.addedge(path.get(j), path.get(j - 1));
			}
			path.clear();
		}		
		
		//clear all visted path
		for (int i = 0; i <anst.edge.size(); i++) {
			visited[anst.edge.get(i).u] = false;
			visited[anst.edge.get(i).v] = false;
		}
		visited[v] = false;
		//store candidate tree in queue 
		if (ansQueue.size() < topk) {	//topk tree
			ansQueue.add(anst);
			return;
		}
		//a smaller tree
		if (ansQueue.peek().weight > anst.weight) {
			ansQueue.poll();
			ansQueue.add(anst);
		}
	}
	
	boolean is_complete(int v) {
		for (double it : dist[v])
			if (it < 0)
				return false;
		return true;
	}
	
	
	public double kdealRandom(GraphInfo gg){
		g1 = gg;
		totalweight = 0;
		return bidir_exp_search();		
	}
	
	public double kdealRandom(GraphInfo gg, List<String> keys){
		g1 = gg;
		g1.givenQueryWord(keys);
		totalweight = 0;
		return bidir_exp_search();		
	}
}
