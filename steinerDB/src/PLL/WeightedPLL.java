package PLL;

import java.io.*;
import java.util.*;

import database.dbReader;
import database.dbWriter;
import graph.WeightedGraph;
import graph.Hop;
import graph.OrderGen;

public class WeightedPLL {
	static String testFile = "E:\\data\\exampleU.txt";
	//static String testFile = "E:\\data\\undirectedGraph(linkedmdb).txt";;
	//static String testFile = "E:\\data\\undirectedGraph(dbpedia201510).txt";	
	//static String testFile = "E:\\data\\3\\com-youtube.ungraph.txt";
	
	int nodeNum;
	long labelSize = 0;
	ArrayList<ArrayList<Hop>> graph; // raw graph edge
	ArrayList<Integer> list; //	
	ArrayList<Double> queryarray; // record nodes visited during prunedBfs	
	ArrayList<ArrayList<EHop>> label; // 2 hop label
	double errorNum = Math.pow(10 , -6);
	int[] oldName, newName;
	EHopMinHeap heap;
	OrderGen og1;
	
	void PrunedBfs(int vk) {
		
		//start form vk
		int newvk = newName[vk];
		heap.push(vk, 0, vk);
			
		// to speed up Query(vk,u)
		Iterator<EHop> it = label.get(vk).iterator();
		while (it.hasNext()) {
			EHop path = (EHop) it.next();
			queryarray.set(path.v, path.dis);
		}

		int u;
		while (!heap.empty()) {
			EHop hu = heap.pop();
			u = hu.v;
			list.add(u);	//u has been changed
			
			// find the nearest 2 hop cover u-vk path and record the length in query
			double query = Double.MAX_VALUE;
			it = label.get(u).iterator();
			while (it.hasNext()) {
				EHop path = (EHop) it.next();
				if (queryarray.get(path.v) > -1 && query > queryarray.get(path.v) + path.dis)
					query = queryarray.get(path.v) + path.dis;
			}

			// if u-vk is not a shorter path than query
			if (query < heap.dist[u] + 1e-8)			
				continue;

			// u need to record a hop to vk
			label.get(u).add(new EHop(newvk, heap.dist[u], hu.par));
			labelSize++;

			// add potential neighbour of u			
			for (int i = 0; i < graph.get(u).size(); i++) {
				Hop edge = graph.get(u).get(i);	
				heap.tryLower(edge.v, heap.dist[u] + edge.dis, u);
			}
		}

		// clear for speed up Query(vk+1,u)
		it = label.get(vk).iterator();
		while (it.hasNext()) {
			EHop path = (EHop) it.next();
			queryarray.set(path.v, -2.0);
		}

		// repair modified p value
		for (int i = 0; i < list.size(); i++) {
			heap.heapLoc[list.get(i)] = -1;
			heap.dist[list.get(i)] = -1;
		}
		list.clear();

		//System.out.println(vk + " is pruned");
	}

	void Prune(WeightedGraph w1) throws IOException {						
		graph = w1.graph;
		nodeNum = w1.nodeNum;
		heap = new EHopMinHeap();
		heap.heapInit(nodeNum);
		queryarray = new ArrayList<Double>();
		for (int i = 0; i < nodeNum; i++)		
			queryarray.add(-2.0);
		
		label = new ArrayList<ArrayList<EHop>>();
		for (int i = 0; i < nodeNum; i++)
			label.add(new ArrayList<EHop>());
		
		list = new ArrayList<Integer>();
		
		//label generate
		int cnt = 0;
		for (int i = 0; i < nodeNum; i++) {
			PrunedBfs(oldName[i]);
			if (labelSize / nodeNum > cnt) {
				cnt = (int) (labelSize / nodeNum);
				if (cnt % 10 == 0)
					System.out.println(i + ":" + cnt);
			}
		}		
				
	}
	
	//query distance between s and t
	double query(int s, int t) {
		double ans = Double.MAX_VALUE;
		int move1 = 1;	//move indicate whether choose next hop
		int move2 = 1;
		EHop path1 = null, path2 = null;
		Iterator<EHop> it1 = label.get(s).iterator();
		Iterator<EHop> it2 = label.get(t).iterator();
		while ((it1.hasNext()||move1==0) && (it2.hasNext()||move2==0)) {
			if (move1 == 1)
				path1 = it1.next();
			if (move2 == 1)
				path2 = it2.next();
			
			//if s-t is in a hop then the result must be right
			if (path1.v == newName[t])
				return path1.dis;
			if (path2.v == newName[s])
				return path2.dis;
			move1 = move2 = 0;
			
			//the program produces hop while targer in increasing order 
			if (path1.v < path2.v)
				move1 = 1;
			if (path1.v > path2.v)
				move2 = 1;
			
			//a potential 2 hop cover
			if (path1.v == path2.v) {
				if (ans > path1.dis + path2.dis)
					ans = path1.dis + path2.dis;
				move1 = move2 = 1;
			}
		}
		return ans;
	}
	
	double getPath(int s, int t) {
		double ans = Double.MAX_VALUE;
		int move1 = 1;	//move indicate whether choose next hop
		int move2 = 1;
		int goal = 0;
		EHop path1 = null, path2 = null;
		Iterator<EHop> it1 = label.get(s).iterator();
		Iterator<EHop> it2 = label.get(t).iterator();
		while ((it1.hasNext()||move1==0) && (it2.hasNext()||move2==0)) {
			if (move1 == 1)
				path1 = it1.next();
			if (move2 == 1)
				path2 = it2.next();
			
			//if s-t is in a hop then the result must be right
			if (path1.v == newName[t]) {
				ans = path1.dis;
				goal = t;
				break;
			}
			if (path2.v == newName[s]) {
				ans = path2.dis;
				goal = s;
				break;
			}
			move1 = move2 = 0;
			
			//the program produces hop while targer in increasing order 
			if (path1.v < path2.v)
				move1 = 1;
			if (path1.v > path2.v)
				move2 = 1;
			
			//a potential 2 hop cover
			if (path1.v == path2.v) {
				if (ans > path1.dis + path2.dis) {
					ans = path1.dis + path2.dis;
					goal = oldName[path1.v];
				}
				move1 = move2 = 1;
			}
		}
		
		sumPath(s,goal);
		sumPath(t,goal);
		
		return ans;
	}
	
	//used for debug CHL
	//get path len from t to s
	void sumPath(int s, int t) {
		int it = s;
		int newt = newName[t];
		int left, right, mid;
		while (it != t) {
			
			//Binary Search for hop (it, newt)
			left = 0;
			right = label.get(it).size() - 1;
			while (left < right) {
				mid = (left + right)/2;
				if (label.get(it).get(mid).v < newt)
					left = mid + 1;
				else
					right = mid;
			}
			if (label.get(it).get(left).v != newt)
				System.out.println("binary wrong");
			
			//now left.par is it's parent
			it = label.get(it).get(left).par;
		}		
	}
		
	void pllDBDeal(String DBName, String mode) throws IOException{
		WeightedGraph w1 = new WeightedGraph();
		og1 = new OrderGen();	//new order generate
		w1.graphDBRead(DBName);
		newName = new int[w1.nodeNum];
		oldName = new int[w1.nodeNum];

		//start chl
		labelSize  = 0;
		if (!mode.equals("PLL")) {
			og1.orderGen2(w1, newName, oldName);
			Prune(w1);
			System.out.println("avg CHL: "+(((double)labelSize)/w1.nodeNum -  1));
			orderOutput(w1, DBName);
			labelOutput(w1, DBName);
		}			
		else {
			og1.orderGen1(w1, newName, oldName);
			Prune(w1);
			System.out.println("avg pll: "+(((double)labelSize)/w1.nodeNum -  1));		
			orderOutputPLL(w1, DBName);
			labelOutputPLL(w1, DBName);
		}
		
	}
	
	//have read the graph
	public void pllDBDeal(String DBName, WeightedGraph ww, String mode) throws IOException{
		WeightedGraph w1 = ww;
		og1 = new OrderGen();	//new order generate		
		newName = new int[w1.nodeNum];
		oldName = new int[w1.nodeNum];
		
		//start chl
		labelSize  = 0;
		if (mode.equals("PLL")) {
			og1.orderGen1(w1, newName, oldName);
			Prune(w1);
			System.out.println("avg PLL: "+(((double)labelSize)/w1.nodeNum + 1));		
			orderOutputPLL(w1, DBName);
			labelOutputPLL(w1, DBName);
		}			
		else {
			og1.orderGen2(w1, newName, oldName);
			Prune(w1);
			System.out.println("avg CHL: "+(((double)labelSize)/w1.nodeNum));
			orderOutput(w1, DBName);
			labelOutput(w1, DBName);			
		}
	}
	
	void orderOutput(WeightedGraph w1, String dbName) {
		dbWriter d1= new dbWriter();
		d1.dbInit(dbName, "subName");
		for (int i = 0; i < nodeNum; i++)
			d1.insertSubName(w1.reHashTable.get(oldName[i]), i);
		d1.dbWriterEnd();
		System.out.println("Sub Name Generate end!");
		
		ArrayList<String> aa = new ArrayList<String>();
		ArrayList<Integer> bb = new ArrayList<Integer>();
		dbReader dh = new dbReader(); 
		dh.dbInit(dbName, "invertedTable");		
		String []a = new String[1];
		int []b = new int [1];
		while(dh.readInvertedTable(a,b)) {
			aa.add(a[0]); 
			bb.add(b[0]);
		}
		dbWriter di = new dbWriter();
		di.dbInit(dbName, "invertedTable");
		for (int i = 0; i < aa.size(); i++)
			di.insertInvertedTable(aa.get(i), newName[bb.get(i)]);
		di.dbWriterEnd();
		
		System.out.println("Inverted Table Generate end!");
	}
	
	void orderOutputPLL(WeightedGraph w1, String dbName) {
		dbWriter d1= new dbWriter();
		d1.dbInit(dbName, "subNamePLL");
		for (int i = 0; i < nodeNum; i++)
			d1.insertSubName(w1.reHashTable.get(oldName[i]), i);
		d1.dbWriterEnd();
		System.out.println("PLL Sub Name Generate end!");

		ArrayList<String> aa = new ArrayList<String>();
		ArrayList<Integer> bb = new ArrayList<Integer>();
		dbReader dh = new dbReader(); 
		dh.dbInit(dbName, "invertedTable");		
		String []a = new String[1];
		int []b = new int [1];
		while(dh.readInvertedTable(a,b)) {
			aa.add(a[0]); 
			bb.add(b[0]);
		}		
		dbWriter di = new dbWriter();
		di.dbInit(dbName, "invertedTablePLL");
		for (int i = 0; i < aa.size(); i++)
			di.insertInvertedTable(aa.get(i), newName[bb.get(i)]);
		di.dbWriterEnd();
		
		System.out.println("PLL Inverted Table Generate end!");
	}
	
	void labelOutput(WeightedGraph w1, String dbName){
		dbWriter d1= new dbWriter();
		d1.dbInit(dbName, "hubLabel");
		// label output
		
		for (int i = 0; i < nodeNum; i++) {
			int u = oldName[i];
			for (int j = 0; j < label.get(u).size(); j++) {
				EHop ehv = label.get(u).get(j);
				d1.insertHubLabel(i, ehv.v, ehv.dis, newName[ehv.par]);
			}			
		}
		d1.dbWriterEnd();
		
		System.out.println("hub Label Generate end!");
	}
	
	void labelOutputPLL(WeightedGraph w1, String dbName){
		dbWriter d1= new dbWriter();
		d1.dbInit(dbName, "hubLabelPLL");
		// label output
		
		for (int i = 0; i < nodeNum; i++) {
			int u = oldName[i];
			for (int j = 0; j < label.get(u).size(); j++) {
				EHop ehv = label.get(u).get(j);
				d1.insertHubLabel(i, ehv.v, ehv.dis, newName[ehv.par]);
			}			
		}
		d1.dbWriterEnd();		
		System.out.println("PLL hub Label Generate end!");
	}
	
	public static void main(String[] args) throws IOException {
	}
}
