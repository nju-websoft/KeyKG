package prunedDP;

import java.util.*;

public class CandiDTree implements Comparable<CandiDTree>{
	int v;
	BitSet X;
	double weight;
	double lb;
	ArrayList<EdgeDP> edges;
	DTree dt;
	
	CandiDTree(int v, int xs){		
		this.v = v;
		X = new BitSet(xs);
		weight = 0;
		lb = 0;
		edges = new ArrayList<EdgeDP>();
		dt = new DTree(v, X);
	}
	
	CandiDTree(int v, int xs, int loc){		
		this.v = v;
		X = new BitSet(xs);
		X.set(loc);
		weight = 0;
		lb = 0;
		edges = new ArrayList<EdgeDP>();
		dt = new DTree(v, X);
	}
	
	CandiDTree(int v, BitSet X, double weight, double lb, ArrayList<EdgeDP> edges){
		this.v = v;
		this.X = (BitSet) X.clone();
		this.weight = weight;
		this.lb = lb;
		this.edges = new ArrayList<EdgeDP>();
		for (EdgeDP it : edges)
			this.edges.add(it);
		dt = new DTree(v, this.X);
	}
	
	public int compareTo(CandiDTree c2) {		
		int result  = 0;
		result = this.lb > c2.lb ? 1 : (this.lb < c2.lb ? -1 : 0);		
		//if (result == 0) result = (this.weight > c2.weight) ? 1 : (this.weight == c2.weight ? 0 : -1);
		if (result == 0)
			result = (this.v < c2.v) ? 1 : (this.v == c2.v ? 0 : -1);
		return result;
	}
}
