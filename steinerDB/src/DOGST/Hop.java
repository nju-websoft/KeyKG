package DOGST;

public class Hop  implements Comparable<Hop>{
	public int v;
	public double dis;
	public Hop(int v, double dis) {
		this.v = v;
		this.dis = dis;
	}
	
	public int compareTo(Hop s2) {		
		int result = this.dis > s2.dis ? 1 : (this.dis < s2.dis ? -1 : 0);
		if (result == 0)
			result = (this.v < s2.v) ? 1 : (this.v == s2.v ? 0 : -1);
		return result;
	}
}
