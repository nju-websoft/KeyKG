package prunedDP;

import java.util.BitSet;

public class DTree implements Comparable<DTree>{
	int v;
	BitSet X;
	
	DTree(int v, int xs){		
		this.v = v;
		X = new BitSet(xs);		
	}
	
	DTree(int v, int xs, int loc){		
		this.v = v;
		X = new BitSet(xs);
		X.set(loc);
	}
	
	DTree(int v, BitSet X){		
		this.v = v;
		this.X = (BitSet) X.clone();		
	}
	
	int getloc() {
		int ans = 0;
		for (int i = X.nextSetBit(0); i >= 0; i = X.nextSetBit(i + 1))
		      ans += X.get(i) ? (1 << i) : 0;
		return ans;
	}
	
	public int compareTo(DTree r2) {		
		int result = this.v > r2.v ? 1 : (this.v < r2.v ? -1 : 0);		
		if (result == 0)
			result = (this.getloc() < r2.getloc()) ? 1 : (this.getloc() == r2.getloc() ? 0 : -1);
		return result;
	}
}
