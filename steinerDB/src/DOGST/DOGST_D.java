package DOGST;

import java.io.*;
import java.util.*;

import PLL.EHop;

public class DOGST_D {
	public ArrayList<String> keywordList;
	ArrayList<ArrayList<Hop>> disList;
	int keyNum;	
	DiskCommonStruct com1;
	int[] realName;	//node name of keyword	
	ArrayList<Hop> hopTerminal;
	public ArrayList<AnsTree> ansList;	//result tree
	public double totalweight = 0.0;
	
	public DOGST_D(DiskCommonStruct c1){
		com1 = c1;		
		hopTerminal = new ArrayList<Hop>();
		for (int i = 0; i < c1.nodeNum; i++) hopTerminal.add(null);
	}
	public void randQuery(DiskCommonStruct c1,int kNum){
		keywordList = new ArrayList<String>();
		
		Random r1 = new Random();	
		int kM;
		if (kNum == 0)
			kM = r1.nextInt(6)+2;
		else kM = kNum;
		int gg;
		int []gu = new int[kM];
		for (int i = 0; i < kM; i++) {
			while (true) {
				gg = r1.nextInt(c1.keywordNum);
				int flag = 0;
				for (int j  = 0; j < i -1 ; j++)	//not generated keyword
					if (gu[j] == gg)
						flag = 1;
				if (flag == 0)
					break;
			}
			keywordList.add(c1.reHashKeywordTable.get(gg));
			gu[i] = gg;
		}
		keyNum = keywordList.size();
		
	}
	
	public void biggestQuery(int kNum){
		keywordList = new ArrayList<String>();
		PriorityQueue<Hop> nodes = new PriorityQueue<Hop>();
		for (int i = 0; i < com1.invTable.size(); i++)
			nodes.add(new Hop(i, -com1.invTable.get(i).size()));
		
		for (int i = 0; i < kNum; i++) {
			Hop it = nodes.poll();
			System.out.println(-it.dis);
			keywordList.add(com1.reHashKeywordTable.get(it.v));
		}
		keyNum = keywordList.size();
	}
	
	//get keyword from file
	void fileQuery(String testfile) throws IOException{
		File kFile = new File(testfile);
		if (!kFile.exists())
			return;
		Scanner input = new Scanner(kFile);
		keywordList = new ArrayList<String>();
		while (input.hasNext()) {
			keywordList.add(input.next());
		}			
		input.close();
		keyNum = keywordList.size();
		for (int i = 0; i < keyNum; i++)
			System.out.print(keywordList.get(i) + " ");
		System.out.println();
	}
	
	double getdisList(int chins) {
		realName = new int[keyNum];
		for (int i = 0; i < keyNum; i++)
			realName[i] = com1.HashKeywordTable.get(keywordList.get(i));
		/*
		{
		//make set 0 has least nodes
				int minloc = 0;
				for (int i = 1; i < keyNum; i++)			
					if (com1.invTable.get(realName[i]).size() < com1.invTable.get(realName[minloc]).size())	//get node set with smallest set size
						minloc = i;
				if (minloc != 0) {	
					int tmp = realName[minloc];
					realName[minloc] = realName[0];
					realName[0] = tmp;
				}
				if (com1.invTable.get(realName[0]).size() != 0) {
					int gg = com1.invTable.get(realName[0]).size();
					if (gg > 200) gg = 200;
					System.out.println(gg);
					return 0;
				}
		}
		*/
		for (int i = 0; i < keyNum; i++)	//get hub label for speed up
			com1.fetchHubLabel(com1.invTable.get(realName[i]));
		
		double keynode = 0;
		for (int i = 1; i < keyNum; i++)	//record all 
			keynode = keynode + com1.invTable.get(realName[i]).size();		
		
		//make set 0 has least nodes
		int minloc = 0;
		for (int i = 1; i < keyNum; i++)			
			if (com1.invTable.get(realName[i]).size() < com1.invTable.get(realName[minloc]).size())	//get node set with smallest set size
				minloc = i;
		if (minloc != 0) {	
			int tmp = realName[minloc];
			realName[minloc] = realName[0];
			realName[0] = tmp;
		}
		if (com1.invTable.get(realName[0]).size() == 0) {
			System.out.println("no such keyword");
			com1.clearHubLabel();
			return keynode;
		}
				
		disList = new ArrayList<ArrayList<Hop>>();
		for (int i = 0; i < com1.invTable.get(realName[0]).size(); i++)
			disList.add(new ArrayList<Hop>());					
		ArrayList<Integer> changedTer = new ArrayList<Integer>();
		for (int i = 1; i < keyNum; i++) {	//for each set i
			for (int v : com1.invTable.get(realName[i])) 
				for (EHop it : com1.label.get(v)) {
					//System.out.println(it.v);
					if (hopTerminal.get(it.v) == null) {
						hopTerminal.set(it.v, new Hop(v, it.dis));
						changedTer.add(it.v);
					} 
					else if (hopTerminal.get(it.v).dis > it.dis)
						hopTerminal.set(it.v, new Hop(v, it.dis));
				}			
			
			for (int j = 0; j < com1.invTable.get(realName[0]).size(); j++) {
				int u = com1.invTable.get(realName[0]).get(j);
				double minhopdis = -1;
				int minhoploc = -1;
				for (EHop it : com1.label.get(u)) {
					if (hopTerminal.get(it.v) == null) continue;
					if (minhopdis < 0 || minhopdis > it.dis + hopTerminal.get(it.v).dis) {
						minhopdis = it.dis + hopTerminal.get(it.v).dis;
						minhoploc = hopTerminal.get(it.v).v;
					}						
				}
				if (minhopdis < 0)
					disList.get(j).add(new Hop(-1, Double.MAX_VALUE));			
				else {					
					disList.get(j).add(new Hop(minhoploc, minhopdis));
				}
			}
			for (int it : changedTer)
				hopTerminal.set(it, null);
			changedTer.clear();
		}
		
		//heap
		Queue<TreeStr> unusedHeap = new PriorityQueue<TreeStr>(new TreeComparator());
		
		//insert all possible top-k trees into heap
		for (int it = 0; it < com1.invTable.get(realName[0]).size(); it++) {
			TreeStr t1 = new TreeStr(keyNum);
			t1.root = it;
			double sum = 0;
			for (int i = 0; i < keyNum - 1; i++)
				sum = sum + disList.get(it).get(i).dis;
			t1.weight = sum;
			unusedHeap.add(t1);
		}
			
		TreeStr t1 = unusedHeap.peek();		
		if (t1.weight >= Double.MAX_VALUE - 2) {
			System.out.println("keyword can't be linked!");
			com1.clearHubLabel();
			return keynode;
		}
		//System.out.println("chings: "+t1.weight);
		//try to add the tree to ans list 
		AnsTree canT = new AnsTree();
		canT.root = com1.invTable.get(realName[0]).get(t1.root);
		//canT.terminal.add(canT.root);
		for (int i = 0; i < keyNum - 1; i++)	//generate candidate ans
			canT.terminal.add(disList.get(t1.root).get(i).v);
		ansList.add(canT);
		//now we generate result
		for (int k = 0; k < ansList.size(); k++) {					
			if (chins == 1) {
				com1.fastchins(ansList.get(k));					
			}
		
			if (chins == 2) {
			//get all terminals
				ArrayList<Integer> terList = new ArrayList<Integer>();
				for(int it : ansList.get(k).terminal)
					terList.add(it);
				
				AnsTree atbest = null;				
				for (int t = 0; t < terList.size(); t++) {	//choose t as root
					AnsTree at1 = new AnsTree();
					for (int it = 0; it < terList.size(); it++)
						if (it != t)
							at1.terminal.add(terList.get(it));
					at1.root = terList.get(t);
					at1.terminal.add(ansList.get(k).root);
					
					com1.fastchins(at1);
					if (atbest == null) {
						atbest= at1;
						continue;
					}
					//System.out.println("top " + t +" "+at1.weight);
					if (atbest.weight > at1.weight)
						atbest = at1;						
				}
				
				com1.fastchins(ansList.get(k));
				if (atbest != null && ansList.get(k).weight > atbest.weight)
					ansList.set(k, atbest);
			}
			
			if (chins == 0){
				AnsTree anst = ansList.get(k);
				Iterator<Integer> it = anst.terminal.iterator();
				while (it.hasNext()) {	//link root to all terminals
					int v = it.next();
					if (v != anst.root)
						com1.merge(anst.root, v, anst);
				}
				for (int i = 0; i <anst.edge.size(); i++) {
					com1.visited[anst.edge.get(i).u] = false;
					com1.visited[anst.edge.get(i).v] = false;
				}							
			}			
			//System.out.println(ansList.get(k).weight);
			totalweight += ansList.get(k).weight;
		}
		
		totalweight = Double.MAX_VALUE;
		for (AnsTree anst : ansList)
			if (totalweight > anst.weight)
				totalweight = anst.weight;
		/*
		for (TreeEdge it1 : ansList.get(0).edge)
			System.out.println(it1.u + " " + it1.v + " "+ 0);
		System.out.println();
		*/
		com1.clearHubLabel();
		return keynode;
	}
	
	double kdealRandom(DiskCommonStruct c1, int chins){
		totalweight = 0.0;
		com1 = c1;		
		randQuery(c1, 10);
		return getdisList(chins)/keyNum;		
	}
	
	public double kdealRandom(DiskCommonStruct c1, ArrayList<String> keys, int chins) {
		totalweight = 0.0;
		com1 = c1;
		keywordList = keys;
		keyNum = keywordList.size();
		ansList = new ArrayList<AnsTree>();		
		return getdisList(chins)/keyNum;
	}
	/*
	void kdealFile(CommonStruct c1, String testfile) throws IOException{
		com1 = c1;
		fileQuery(testfile);
		getdisList();
	}
	*/
}
