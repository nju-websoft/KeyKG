package DOGST;

import java.util.*;

import PLL.EHop;
import database.dbReader;

public class DiskCommonStruct {
	int nodeNum = 0;
	int keywordNum = 0;
	public ArrayList<String> reHashTable;	//rehash Integer to Entity Name
	ArrayList<ArrayList<Integer>> invTable;	//keyword inverted table
	public HashMap<String, Integer> HashKeywordTable;	//keyword hash
	public ArrayList<String> reHashKeywordTable;	//rehash Integer to keyword
	ArrayList<ArrayList<EHop>> label; // 2 hop label
	boolean []visited;
	String database;
	dbReader dbDisk;
	ArrayList<Integer> changeList;
	ArrayList<Hop> treeLabel;
	ArrayList<Integer> treeLabelChanged;
	ArrayList<PriorityQueue<Hop>> terminalLabel;
	
	//read entity hash table
	void readhashTable(String DBName) {
		//nodeNum = 0;
		reHashTable = new ArrayList<String>();
		dbReader dh = new dbReader(); 
		dh.dbInit(DBName, "subName");		
		String []a = new String[1];
		int []b = new int [1];
		while(dh.readSubName(a,b)){			
			reHashTable.add(a[0]);						
		}
		nodeNum = reHashTable.size();
	}
	
	//read keywrod hash table
	void readKeyMap(String DBName) {
		HashKeywordTable = new HashMap<String, Integer>();
		reHashKeywordTable = new ArrayList<String>();
		dbReader dh = new dbReader(); 
		dh.dbInit(DBName, "keyMap");				
		int []a = new int [1];
		String []b = new String[1];
		while(dh.readKeyMap(a,b)){
			reHashKeywordTable.add(b[0]);
			if (!HashKeywordTable.containsKey(b[0])) {
				HashKeywordTable.put(b[0], a[0]);				
			}
		}
		keywordNum = reHashKeywordTable.size();
	}
	
	//read inverted table
	void readInvTable(String DBName) {
		invTable = new ArrayList<ArrayList<Integer>>();
		for (int i = 0; i < HashKeywordTable.size(); i++)
			invTable.add(new ArrayList<Integer>());
		
		dbReader dh = new dbReader(); 
		dh.dbInit(DBName, "invertedTable");		
		String []a = new String[1];
		int []b = new int [1];
		while(dh.readInvertedTable(a,b)){
			invTable.get(HashKeywordTable.get(a[0])).add(b[0]);
		}
	}
	
	//read hub label
	void initHubLabel(){
		label = new ArrayList<ArrayList<EHop>>();
		for (int i = 0; i < nodeNum; i++) label.add(null);
		changeList = new ArrayList<Integer>();
	}
	
	void clearHubLabel() {
		for (int it : changeList) label.set(it, null);
		changeList.clear();
			//label.get(need.get(i).intValue()).clear();
	}
	
	void fetchHubLabel(ArrayList<Integer> need) {		
		for (int i : need) 
			fetchHubLabel(i);
	}
	
	void fetchHubLabel(int i) {
		if (label.get(i) == null){
			ArrayList<EHop> newla = new ArrayList<EHop>();
			dbDisk.dbGet(i, newla);
			label.set(i, newla);
			changeList.add(i);
		}
	}
	
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
			if (path1.v == t)
				return ans = path1.dis;				
			if (path2.v == s)
				ans = path2.dis;
			move1 = move2 = 0;
			
			//the program produces hop while target in increasing order 
			if (path1.v < path2.v)
				move1 = 1;
			if (path1.v > path2.v)
				move2 = 1;
			
			//a potential 2 hop cover
			if (path1.v == path2.v) {
				if (ans > path1.dis + path2.dis) {
					ans = path1.dis + path2.dis;					
				}
				move1 = move2 = 1;
			}
		}
		return ans;
	}
	
	//merge shortest path between s and t in the result tree
	void merge(int root, int t, AnsTree tree) {	//s is the root
		int s = root;
		double ans = Double.MAX_VALUE;
		int move1 = 1;	//move indicate whether choose next hop
		int move2 = 1;
		int goal = t;
		EHop path1 = null, path2 = null;
		Iterator<EHop> it1 = label.get(s).iterator();
		Iterator<EHop> it2 = label.get(t).iterator();
		while ((it1.hasNext()||move1==0) && (it2.hasNext()||move2==0)) {
			if (move1 == 1)
				path1 = it1.next();
			if (move2 == 1)
				path2 = it2.next();
			
			//if s-t is in a hop then the result must be right
			if (path1.v == t) {
				ans = path1.dis;
				goal = t;
				break;
			}
			if (path2.v == s) {
				ans = path2.dis;
				goal = s;
				break;
			}
			move1 = move2 = 0;
			
			//the program produces hop while target in increasing order 
			if (path1.v < path2.v)
				move1 = 1;
			if (path1.v > path2.v)
				move2 = 1;
			
			//a potential 2 hop cover
			if (path1.v == path2.v) {
				if (ans > path1.dis + path2.dis) {
					ans = path1.dis + path2.dis;
					goal = path1.v;
				}
				move1 = move2 = 1;
			}
		}
		if (t == goal)	//root to t
			mergeRevPath(root, goal, tree);
		if (root == goal)	//t to root
			mergePath(t, goal, tree);
		if (t != goal&&root != goal) {	//t to goal, root to goal
			if (!mergePath(t, goal, tree))	//don't reach visited node
				mergeRevPath(root, goal, tree);	//root to goal should have reverse edge
		}
				
	}
	
	
	void mergeRevPath(int root, int t, AnsTree tree) {	//path from root to t
		ArrayList<Hop> tmpPath = new ArrayList<Hop>();		
		int it = root;
		int[] a = new int[1];
		double[] b = new double[1];
		while (it != t) {
			dbDisk.dbGet(it, t, a, b);
			tmpPath.add(new Hop(it, b[0])); // the dis between it and t
			it = a[0];
		}
		tmpPath.add(new Hop(t, 0));		
		for (int i = tmpPath.size() -  1; i > 0; i--) 
			if (!visited[tmpPath.get(i).v]){	//not visited node
				visited[tmpPath.get(i).v] = true;
				tree.edge.add(new TreeEdge(tmpPath.get(i).v, tmpPath.get(i - 1).v));				
			}
			else {	//visited node
				tree.weight = tree.weight + tmpPath.get(i).dis;
				return;
			}
		tree.weight = tree.weight + tmpPath.get(0).dis;
	}
	
	//merge shortest path from s to t
	//return true if reaches visited node
	boolean mergePath(int s, int t, AnsTree tree) {	//there s is the end point
		ArrayList<Hop> tmpPath = new ArrayList<Hop>();
		int it = s;	
		int []a = new int[1];
		double []b = new double[1];
		while (it != t) {			
			dbDisk.dbGet(it, t, a, b);
			tmpPath.add(new Hop(it, b[0]));	//the dis between it and t
			it = a[0];
		}
		tmpPath.add(new Hop(t, 0));
		tree.weight = tree.weight + tmpPath.get(0).dis;
		for (int i = 0; i< tmpPath.size() - 1; i++) 
			if (!visited[tmpPath.get(i).v]){	//not visited node
				visited[tmpPath.get(i).v] = true;
				tree.edge.add(new TreeEdge(tmpPath.get(i).v, tmpPath.get(i + 1).v));				
			}
			else {	//visited node
				tree.weight = tree.weight - tmpPath.get(i).dis;
				return true;
			}
		return false;
	}
	
	//add new to tree
	void treeLabelAdd(int u) {		
		if (label.get(u)== null) fetchHubLabel(u);
		for (EHop it : label.get(u)) {
			if (treeLabel.get(it.v) == null) {	//a new node labeled
				treeLabel.set(it.v, new Hop(u, it.dis));
				treeLabelChanged.add(it.v);
			}
			else
				if (it.dis < treeLabel.get(it.v).dis)
					treeLabel.set(it.v, new Hop(u, it.dis));
		}
	}
	
	
	void chinsmerge(int s, int t, AnsTree anst) {
		int it = s;		
		int left, right, mid;
		while (it != t) {
			treeLabelAdd(it);
			//Binary Search for hop (it, newt)
			left = 0;
			right = label.get(it).size() - 1;
			while (left < right) {
				mid = (left + right)/2;
				if (label.get(it).get(mid).v < t)
					left = mid + 1;
				else
					right = mid;
			}			
			if (label.get(it).get(left).v != t)
				System.out.println("binary wrong");			
			
			anst.edge.add(new TreeEdge(it, label.get(it).get(left).par));
			//now left.par is it's parent			
			it = label.get(it).get(left).par;
		}
	}
	
	//geterate answer tree 
	void fastchins(AnsTree anst) {
		//anst.weight = 0;
		treeLabelAdd(anst.root);
		
		Map<Integer, Integer> terTable = new HashMap<Integer, Integer>();
		for (int u : anst.terminal)
			terTable.put(u, terTable.size());	
		
		//get label set of terminals waitting to be linked to answer tree
		Map<Integer, Integer> nodeTable = new HashMap<Integer, Integer>();
		ArrayList<Integer> nameNode = new ArrayList<Integer>();	//the real name of i 
		for (int u : anst.terminal) {
			for (EHop it : label.get(u)) {
				if (!nodeTable.containsKey(it.v)) {
					nodeTable.put(it.v, nodeTable.size());
					nameNode.add(it.v);
					terminalLabel.add(new PriorityQueue<Hop>());
				}
				terminalLabel.get(nodeTable.get(it.v)).offer(new Hop(u, it.dis));	
			}
		}		
		nodeTable.clear();
		
		while (!terminalLabel.isEmpty()) {
			
			//choose next terminal to be linked to tree			
			double minLen = -1;
			int t = 0, u = 0;
			int cnt = 0;
			while (cnt < terminalLabel.size()) {
				Queue<Hop> relay = terminalLabel.get(cnt);
				while (!relay.isEmpty() && terTable.get(relay.peek().v) == -1) relay.poll();
				if (relay.isEmpty()) {	//terminalLabel[cnt] is dropped
					int gg = terminalLabel.size() - 1;
					terminalLabel.set(cnt, terminalLabel.get(gg));
					nameNode.set(cnt, nameNode.get(gg));	//change the name of cnt
					terminalLabel.remove(gg);
					continue;
				}
				Hop candi = relay.peek();
				if (treeLabel.get(nameNode.get(cnt)) != null) 
					if (minLen < 0 || treeLabel.get(nameNode.get(cnt)).dis + candi.dis < minLen) {					
						minLen = treeLabel.get(nameNode.get(cnt)).dis + candi.dis;					
						t = candi.v;						
						u = nameNode.get(cnt);
					}
				cnt++;
			}
			
			if (minLen >= 0) 
			{
				int s = treeLabel.get(u).v;
				//now we konw the shorest path is s-u-t
				treeLabelAdd(u);	//add u to tree
				anst.weight += minLen;	//update tree weight			
				chinsmerge(s, u, anst);	//merge s-u to tree			
				chinsmerge(t, u, anst);	//merge t-u to tree	
				//System.out.println("tree find:"+s + " "+ t +" "+minLen);
				terTable.put(t, -1);	//t is in the tree, not in waiting list
			}
		}
		
		for (int i : treeLabelChanged){
			treeLabel.set(i, null);
		}
		treeLabelChanged.clear();
		
	}
	
	public void Init(String DBName){
		database = DBName;	
		readhashTable(DBName);
		initHubLabel();
		readKeyMap(DBName);
		readInvTable(DBName);		
		
		visited = new boolean[nodeNum];
		for (int i = 0; i < nodeNum; i++)
			visited[i] = false;
		
		treeLabel = new ArrayList<Hop>();
		for (int i = 0; i < nodeNum; i++)
			treeLabel.add(null);
		treeLabelChanged = new ArrayList<Integer>();		
		terminalLabel = new ArrayList<PriorityQueue<Hop>>();
		
		System.out.println("graph init end!");
		dbDisk = new dbReader();
		dbDisk.dbInitDisk(database, "hubLabel");
	}
	
	void output() {
		System.out.println(nodeNum);
	}
}
