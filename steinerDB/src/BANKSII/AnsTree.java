package BANKSII;

import java.util.*;

public class AnsTree implements Comparable<AnsTree>{
	public class TreeEdge {
		public int u;
		public int v;
		TreeEdge(int u, int v){
			this.u = u;
			this.v = v;
		}
	}
	 	
	public ArrayList<TreeEdge> edge;	//record tree edge
	public double weight;
	
	AnsTree(){		
		edge = new ArrayList<TreeEdge>(); 
	}
	
	void addedge(int u, int v) {
		edge.add(new TreeEdge(u,v));
	}
	
	void outputedge() {
		for (TreeEdge it : edge)
			System.out.println(it.u + " " + it.v);
	}
	public int compareTo(AnsTree a2) {	
		if (this.weight < a2.weight)
			return 1;
		return -1;
	}
}
