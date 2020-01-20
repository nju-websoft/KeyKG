package DPBF;

import java.util.*;

public class CandiTree implements Comparable<CandiTree>{
	int v;
	BitSet X;
	double weight;	
	ArrayList<EdgeDPBF> edges;
	
	CandiTree(int v, int xs){		
		this.v = v;
		X = new BitSet(xs);
		weight = 0;		
		edges = new ArrayList<EdgeDPBF>();
	}
	
	CandiTree(int v, int xs, int loc){		
		this.v = v;
		X = new BitSet(xs);
		X.set(loc);
		weight = 0;
		edges = new ArrayList<EdgeDPBF>();
	}
	
	CandiTree(int v, BitSet X, double weight, ArrayList<EdgeDPBF> edges){
		this.v = v;
		this.X = (BitSet) X.clone();
		this.weight = weight;		
		this.edges = new ArrayList<EdgeDPBF>();
		for (EdgeDPBF it : edges)
			this.edges.add(it);
	}
	
	public int compareTo(CandiTree c2) {		
		int result = this.weight > c2.weight ? 1 : (this.weight < c2.weight ? -1 : 0);
		if (result == 0)
			result = (this.v < c2.v) ? 1 : (this.v == c2.v ? 0 : -1);
		return result;
	}
}
