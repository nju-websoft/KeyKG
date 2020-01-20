package DOGST;

import java.util.*;

import PLL.EHop;

public class DOGST {
	public ArrayList<String> keywordList;
	ArrayList<ArrayList<Hop>> disList;
	int keyNum;	
	CommonStruct com1;
	int[] realName;	//node name of keyword	
	ArrayList<hopage> hopTerminal;
	public ArrayList<AnsTree> ansList;	//result tree
	public double totalweight = 0.0;
	int epoch;
	class hopage{		
		int v;
		double dis;
		int age;
		hopage(int age){
			this.age = age;			
		}
		
		void trychange(int v, double dis){
			if (age < epoch) {
				age = epoch;
				this.v = v;
				this.dis =dis;
				return;
			}
			if (dis < this.dis) {
				this.v = v;
				this.dis =dis;
			}			
			return;
		}
	}
	public DOGST(CommonStruct c1){
		com1 = c1;		
		hopTerminal = new ArrayList<hopage>();
		epoch = 0;
		for (int i = 0; i < c1.nodeNum; i++) hopTerminal.add(new hopage(epoch));
	}
	public void randQuery(CommonStruct c1,int kNum){
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
		
		/*
		for (int i = 0; i < keyNum; i++)
			System.out.print(keywordList.get(i) + " ");
		System.out.println();
		*/
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
		
		disList = new ArrayList<ArrayList<Hop>>();
		for (int i = 0; i < com1.invTable.get(realName[0]).size(); i++)					
			disList.add(new ArrayList<Hop>());
		
		//ArrayList<Integer> changedTer = new ArrayList<Integer>();
		for (int i = 1; i < keyNum; i++) {	//for each set i
			epoch++;
			for (int v : com1.invTable.get(realName[i])) 
				for (EHop it : com1.label.get(v))
					hopTerminal.get(it.v).trychange(v, it.dis);			
			
			for (int j = 0; j < com1.invTable.get(realName[0]).size(); j++) {
				int u = com1.invTable.get(realName[0]).get(j);
				double minhopdis = -1;
				int minhoploc = -1;
				for (EHop it : com1.label.get(u)) {
					if (hopTerminal.get(it.v).age < epoch) continue;
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
			return keynode;
		}
		//try to add the tree to ans list 
		AnsTree canT = new AnsTree();
		canT.root = com1.invTable.get(realName[0]).get(t1.root);
		//canT.terminal.add(canT.root);
		for (int i = 0; i < keyNum - 1; i++)	//generate candidate ans
			canT.terminal.add(disList.get(t1.root).get(i).v);
		ansList.add(canT);
		/*
		System.out.print(canT.root);
		for (int v : canT.terminal) System.out.print(" " + v);*/
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
				if (atbest!= null && ansList.get(k).weight > atbest.weight)
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
		return keynode;
	}
	
	double kdealRandom(CommonStruct c1, int chins){
		totalweight = 0.0;
		com1 = c1;		
		randQuery(c1, 10);
		return getdisList(chins)/keyNum;		
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
