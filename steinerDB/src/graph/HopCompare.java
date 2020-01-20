package graph;

import java.util.*;

public class HopCompare  implements Comparator<Hop>{
	public int compare(Hop s1, Hop s2) {		
		int result = s1.dis < s2.dis ? 1 : (s1.dis > s2.dis ? -1 : 0);
		if (result == 0)
			result = (s1.v > s2.v) ? 1 : (s1.v == s2.v ? 0 : -1);
		return result;
	}
}
