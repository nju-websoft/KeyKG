package prunedDP;

import java.io.*;
import java.util.*;



public class MST {
	int parent[];
	int tag[];
	int nodeNum;
	int tagNum;
	static String testFile = "E:\\data\\test.txt";
	PrintWriter pw;
	boolean fileoutput = false;
	boolean primoutput = false;
	boolean treetest = false;
	public MST(int nodeNum){
		
		this.nodeNum = nodeNum;
		
		parent = new int[nodeNum];
		for (int i = 0; i < nodeNum; i++) parent[i] = i;
		tagNum =-1;
		tag = new int[nodeNum];
		for (int i = 0; i < nodeNum; i++) tag[i] = 0;
		
		try {
			if (fileoutput)
				pw = new PrintWriter(new File(testFile));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
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
	
	void testTree(ArrayList<EdgeDP> edges1, ArrayList<EdgeDP> edges2) {
		ArrayList<EdgeDP> edges = new ArrayList<EdgeDP>();
		if (edges1 != null)
		for (EdgeDP it : edges1) edges.add(it);
		if (edges2 != null)
			for (EdgeDP it : edges2) edges.add(it);
		int fa = find(edges.get(0).u);
		boolean flag = true;
		for (EdgeDP it : edges) {
			if (fa != find(it.u)) flag = false;
			if (fa != find(it.v)) flag = false;
		}
		if (!flag)
			System.out.println("Not a tree");
	}
	
	double prim(ArrayList<EdgeDP> edges1, ArrayList<EdgeDP> edges2) {
		PriorityQueue<EdgeDP> edges = new PriorityQueue<EdgeDP>();
		if (edges1 != null)
			for (EdgeDP it : edges1)
				edges.add(it);
		if (edges2 != null)
			for (EdgeDP it : edges2)
				edges.add(it);
		tagNum++;
		double ans = 0;
		if (primoutput) {
			if (fileoutput)
				pw.println("Prim start");
			else
				System.out.println("Prim start");
		}
		
		while (!edges.isEmpty()) {
			EdgeDP it = edges.poll();
			if (find(it.u) != find(it.v)) {
				ans += it.weight;
				union(it.u, it.v);				
				if (primoutput) {
					if (fileoutput)
						pw.println(it.u + " " + it.v + " " + it.weight);
					else
						System.out.println(it.u + " " + it.v + " " + it.weight);
				}
			}
		}
		if (primoutput) {
			if (fileoutput)
				pw.println("Prim : " + ans);
			else
				System.out.println("Prim : " + ans);			
		}
		if (treetest) testTree(edges1, edges2);
		return ans;
	}

}
