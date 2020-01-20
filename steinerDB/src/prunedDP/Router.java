package prunedDP;

import java.util.*;


public class Router implements Comparable<Router>{
	int u, v;
	BitSet X;
	double weight;
	
	Router(int u, int v, int xs){
		this.u = u;
		this.v = v;
		X = new BitSet(xs);
		weight = 0;
	}
	
	Router(int u, int v, int xs, int loc){
		this.u = u;
		this.v = v;
		X = new BitSet(xs);
		X.set(loc);
		weight = 0;
	}
	
	public int compareTo(Router r2) {		
		int result = this.weight > r2.weight ? 1 : (this.weight < r2.weight ? -1 : 0);
		if (result == 0)
			result = (this.u < r2.u) ? 1 : (this.u == r2.u ? 0 : -1);
		if (result == 0)
			result = (this.v < r2.v) ? 1 : (this.v == r2.v ? 0 : -1);
		return result;
	}
}
