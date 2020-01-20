package graph;

import java.io.*;
import java.util.*;



public class GraphChange {

	class Edge implements Comparable<Edge>{
		String e;
		int v;
		
		Edge(String y, int z){
			e = y;
			v =z;
		}
		public int compareTo(Edge e2) {		
			int result = this.v > e2.v ? 1 : (this.v < e2.v ? -1 : 0);			
			return result;
		}
	}
	
	Map<String, Integer> hashNodeTable;
	ArrayList<String> reHashTable;
	ArrayList<Set<Edge>> graph;
	boolean[] visited;
	
	int tryAdd(String x) {
		if (!hashNodeTable.containsKey(x)) {
			hashNodeTable.put(x, hashNodeTable.size());
			reHashTable.add(x);
			graph.add(new TreeSet<Edge>());
			//System.out.println(hashNodeTable.size());
		}
		return hashNodeTable.get(x);
	}
		
	public void writeToOut(String fileName, String fileout) throws IOException {
		hashNodeTable = new HashMap<String, Integer>();		
		reHashTable = new ArrayList<String>();
		graph = new ArrayList<Set<Edge>>();				
        		
		Scanner input1 = new Scanner(new File(fileName));
		//input1.nextInt();
		while (input1.hasNext()) {
			String x = ((Integer) input1.nextInt()).toString();
			String y = ((Integer) input1.nextInt()).toString();
			//input1.nextInt();
			int ux = tryAdd(x);
			int uz = tryAdd(y);
			graph.get(ux).add(new Edge(null, uz));
			graph.get(uz).add(new Edge(null, ux));
		}
		input1.close();
		
		int loc = 0;
		for (int i = 1; i < graph.size(); i++)
			if (graph.get(i).size() > graph.get(loc).size())
				loc = i;
				
		visited = new boolean[graph.size()];
		for (int i = 1; i < graph.size(); i++) visited[i] = false;
		int nodeNum = 0;
		visited[loc] = true;
		ArrayList<Integer> vlist = new ArrayList<Integer>();
		vlist.add(loc);
		while (nodeNum < vlist.size()) {
			int u = vlist.get(nodeNum);
			for (Edge it : graph.get(u))
				if(!visited[it.v]) {
					visited[it.v] = true;
					vlist.add(it.v);
				}
			nodeNum++;
		}
		System.out.println(nodeNum);			
		int name[] = new int[graph.size()];
		int cnt = 0;
		for (int i = 0; i < graph.size(); i++)
			if (visited[i]) {
				name[i] = cnt;
				cnt++;
		}
		
		int[] edgenum = new int[graph.size()];
		for (int i = 0; i < graph.size(); i++) edgenum[i] = 0;
		for (int i = 0; i < graph.size(); i++)
			if (visited[i]) {
				for (Edge it : graph.get(i))
					if (i < it.v && visited[it.v])			
						edgenum[i]++;
			}
		
		PrintWriter fopinfo;
		//File file = new File(fileName+".regular");
		System.out.println(fileout);
		File file = new File(fileout);
		if (!file.exists()) file.createNewFile();
		fopinfo = new PrintWriter(file);
		//fopinfo.println(cnt);
		for (int i = 0; i < graph.size(); i++)
			if (visited[i]) {
				for (Edge it : graph.get(i))
					if (i < it.v && visited[it.v]) {
						double nn = edgenum[i];
						if (edgenum[it.v] > nn) nn = edgenum[it.v];
						nn = Math.log(1 + nn)/Math.log(2.0);
						fopinfo.println(name[i] + " " + name[it.v] + " " + nn);
						//fopinfo.println(name[i] + " " + name[it.v]);
					}
			}
		fopinfo.close();
		System.out.println("Edge weights write end!");
	}
	
	public void writelinkedmdb(String fileName, String fileout) throws IOException {
		hashNodeTable = new HashMap<String, Integer>();		
		reHashTable = new ArrayList<String>();
		graph = new ArrayList<Set<Edge>>();				
        		
		Scanner input1 = new Scanner(new File(fileName));
		//input1.nextInt();
		while (input1.hasNext()) {
			String x =  input1.next();
			String y =  input1.next();
			String z =  input1.next();		
			int ux = tryAdd(x);
			int uz = tryAdd(z);
			graph.get(ux).add(new Edge(y, uz));
			graph.get(uz).add(new Edge(y, ux));
		}
		input1.close();
		
		int loc = 0;
		for (int i = 1; i < graph.size(); i++)
			if (graph.get(i).size() > graph.get(loc).size())
				loc = i;
				
		visited = new boolean[graph.size()];
		for (int i = 1; i < graph.size(); i++) visited[i] = false;
		int nodeNum = 0;
		visited[loc] = true;
		ArrayList<Integer> vlist = new ArrayList<Integer>();
		vlist.add(loc);
		while (nodeNum < vlist.size()) {
			int u = vlist.get(nodeNum);
			for (Edge it : graph.get(u))
				if(!visited[it.v]) {
					visited[it.v] = true;
					vlist.add(it.v);
				}
			nodeNum++;
		}
		System.out.println(nodeNum);			
		int name[] = new int[graph.size()];
		int cnt = 0;
		for (int i = 0; i < graph.size(); i++)
			if (visited[i]) {
				name[i] = cnt;
				cnt++;
		}
		
		int[] edgenum = new int[graph.size()];
		for (int i = 0; i < graph.size(); i++) edgenum[i] = 0;
		for (int i = 0; i < graph.size(); i++)
			if (visited[i]) {
				for (Edge it : graph.get(i))
					if (i < it.v && visited[it.v])			
						edgenum[i]++;
			}
		
		Map<String, Integer> mapIt = new HashMap<String, Integer>();
		
		for (int i = 0; i < graph.size(); i++)
			if (visited[i]) {
				for (Edge it : graph.get(i))
					if (i < it.v && visited[it.v]) {
						if (!mapIt.containsKey(it.e)) mapIt.put(it.e, 0);
						mapIt.put(it.e, mapIt.get(it.e) + 1);
					}
						
			}
		
		PrintWriter fopinfo;
		File file = new File(fileout);
		if (!file.exists()) file.createNewFile();
		fopinfo = new PrintWriter(file);
		//fopinfo.println(cnt);
		for (int i = 0; i < graph.size(); i++)
			if (visited[i]) {
				for (Edge it : graph.get(i))
					if (i < it.v && visited[it.v]) {
						double nn = edgenum[i];
						if (edgenum[it.v] > nn) nn = edgenum[it.v];
						nn = Math.log(1 + nn)/Math.log(2.0);
						//fopinfo.println(name[i] + " " + name[it.v]);
						//fopinfo.println(name[i] + " " + name[it.v] + " " + nn);
						//fopinfo.println(name[i] + " " + name[it.v] + " " + 1);
						fopinfo.println(name[i] + " " + name[it.v] + " " + (Math.log(mapIt.get(it.e)+1)/Math.log(2.0)));
					}
			}
		fopinfo.close();
		System.out.println("Edge weights write end!");
	}
	
	
	
	void dbdeal() throws IOException  {
		String []fileName = new String[5]; 			
		fileName[0] = "E:\\data\\undirectedGraph(dbpedia201510).txt";
		fileName[1] = "E:\\data\\undirectedGraph(linkedmdb).txt";
		fileName[2] = "E:\\data\\pedia16.txt";
		fileName[3] = "E:\\data\\linkedmdb.txt";
		fileName[4] = "E:\\data\\mondial.txt";
		
		String []outfile = new String[5];
		outfile[0] = "E:\\data\\pedia16node.txt";
		outfile[2] = "E:\\data\\linkedmdbnode.txt";
		outfile[3] = "E:\\data\\linkedmdbedge.txt";
		outfile[4] = "E:\\data\\mondialgg.txt";
		
		
		String st = fileName[1];
		String ss = outfile[2];
		
		writeToOut(st, ss);
		//writelinkedmdb(st, ss);

		
	}
	
	public static void main(String[] args) throws IOException {
		GraphChange link = new GraphChange();
		link.dbdeal();		
	}
	
}
