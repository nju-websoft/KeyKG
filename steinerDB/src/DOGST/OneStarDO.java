package DOGST;

import java.util.*;

public class OneStarDO {
	public ArrayList<String> keywordList;
	ArrayList<ArrayList<Hop>> disList;
	int keyNum;	
	CommonStruct com1;
	int[] realName;	//node name of keyword
	int topk = 100;
	public ArrayList<AnsTree> ansList;	//result tree
	public double totalweight = 0.0;
	
	public OneStarDO(CommonStruct c1){
		com1 = c1;
	}
	
	double getdisList(int chins) {
		realName = new int[keyNum];
		for (int i = 0; i < keyNum; i++)
			realName[i] = com1.HashKeywordTable.get(keywordList.get(i));
		
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
			return keynode;
		}
				
		/*
		long x = 0, y = 0;
		for (int v : com1.invTable.get(realName[0])) {
			x = x + com1.label.get(v).size();
		}
		for (int i = 1; i < keyNum; i++)
			for (int v : com1.invTable.get(realName[i]))
				y = y + com1.label.get(v).size();
		System.out.println((x+y)+ " " + x + " " + y + " " + x*y + " " + (x*(keyNum-1)+y));
		if (true) return keynode;
		*/
		//get dis from u to all v,where u is in set 0, v is node of any set 
		disList = new ArrayList<ArrayList<Hop>>();
		for (int i = 0; i < com1.invTable.get(realName[0]).size(); i++)
			disList.add(new ArrayList<Hop>());
				
		
		for (int it = 0; it < com1.invTable.get(realName[0]).size(); it++) {	//for each node in set 0
			int u = com1.invTable.get(realName[0]).get(it);			
			
			for (int i = 1; i < keyNum; i++) {	//for each set i
				double minhopdis = -1;
				int minhoploc = -1;
				for (int j = 0; j < com1.invTable.get(realName[i]).size(); j++) {	//for each node in set i
					int v = com1.invTable.get(realName[i]).get(j);
					double dis = com1.query(u, v);
					if (dis < Double.MAX_VALUE - 1) {
						if (minhopdis < 0 || minhopdis > dis) {
							minhopdis = dis;
							minhoploc = v;
						}							
					}							
				}				
				if (minhopdis < 0)
					disList.get(it).add(new Hop(-1, Double.MAX_VALUE));			
				else					
					disList.get(it).add(new Hop(minhoploc, minhopdis));				
			}
		}
			
		TreeStr tsmall = null;				
		for (int it = 0; it < com1.invTable.get(realName[0]).size(); it++) {
			TreeStr t1 = new TreeStr(keyNum);
			t1.root = it;
			double sum = 0;
			for (int i = 0; i < keyNum - 1; i++)
				sum = sum + disList.get(it).get(i).dis;
			t1.weight = sum;
			if (tsmall == null|| t1.weight < tsmall.weight)
				tsmall = t1;			
		}
		
		if (tsmall.weight >= Double.MAX_VALUE - 2) {
			System.out.println("keyword can't be linked!");
			return keynode;
		}
				
		AnsTree canT = new AnsTree();
		canT.root = com1.invTable.get(realName[0]).get(tsmall.root);
		//canT.terminal.add(canT.root);
		for (int i = 0; i < keyNum - 1; i++)	//generate candidate ans
			canT.terminal.add(disList.get(tsmall.root).get(i).v);
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
		return keynode;
	}
	
	
	public double kdealRandom(CommonStruct c1, ArrayList<String> keys, int chins) {
		totalweight = 0.0;
		com1 = c1;
		keywordList = keys;
		keyNum = keywordList.size();
		ansList = new ArrayList<AnsTree>();
		return getdisList(chins)/keyNum;
	}
}
