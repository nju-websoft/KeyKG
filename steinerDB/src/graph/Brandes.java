package graph;

import java.io.*;
import java.util.*;


public class Brandes {
	static String testFile = "E:\\data\\exampleU.txt";
	//static String testFile = "E:\\data\\test.txt";
	int nodeNum;	
	double []centrality;
	double []distance;
	int []fai;	
	
	double []delta;
	ArrayList<ArrayList<Hop>> graph;
	
	public Brandes() {
	}
	public Brandes(WeightedGraph w1) {
		graph = w1.graph;
		nodeNum = w1.nodeNum;

		centrality = new double[nodeNum];
		distance = new double[nodeNum];
		fai = new int[nodeNum];
		delta = new double[nodeNum];		

		for (int i = 0; i < nodeNum; i++)
			centrality[i] = 0;
	}
	
	//use dfs to
	void dfs(int s){		
		for (int i = 0; i < nodeNum; i++)
			fai[i] = 0;
		for (int i = 0; i < nodeNum; i++)
			delta[i] = 0;		
		
		HopMinHeap heap = new HopMinHeap();
		heap.heapInit(nodeNum);
		heap.push(s, 0);		
		fai[s] = 1;
		ArrayList<Integer> list = new ArrayList<Integer>();				
		
		//form v to u
		while (!heap.empty()) {
			Hop sv = heap.pop();
			int v = sv.v;
			list.add(v);
			for (int i = 0; i < graph.get(v).size(); i++) {
				Hop edge = graph.get(v).get(i);
				int u = edge.v;			
				heap.tryLower(u, sv.dis + edge.dis);	
				if (heap.dist[u] + edge.dis <= sv.dis + 1e-8)	//shortest path to v view u
					fai[v] = fai[u] + fai[v];
			}
		}
		/*
		for (int i = 0; i < nodeNum; i++)
			System.out.println(i + " "+ heap.dist[i]);
		*/
		int order = list.size();
		while (order > 1) {
			order--;
			int w = list.get(order);
			for (int i = 0; i < graph.get(w).size(); i++) {
				Hop edge = graph.get(w).get(i);
				int v = edge.v;				
				if (heap.dist[v] + edge.dis <= heap.dist[w])
					delta[v] = delta[v] + ((double)fai[v])/fai[w]*(1 + delta[w]);
			}
			centrality[w] = centrality[w] + delta[w];
		}
		list.clear();		
	}
	
	void betweenness(WeightedGraph w1) {
		graph = w1.graph;
		nodeNum = w1.nodeNum;
		
		centrality = new double[nodeNum];
		distance = new double[nodeNum];
		fai = new int[nodeNum];
		delta = new double[nodeNum];
		
		for (int i = 0; i < nodeNum; i++)
			centrality[i] = 0;		
		for (int i = 0; i < nodeNum; i++) {			
			dfs(i);
			System.out.println(i + " end");
		}
		for (int i = 0; i < nodeNum; i++)
			System.out.println(centrality[i]);
	}
	
	public static void main(String[] args) throws IOException {
	}
}