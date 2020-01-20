package DOGST;

import java.util.Comparator;

public class TreeComparator  implements Comparator<TreeStr>{
	public int compare(TreeStr s1, TreeStr s2) {
		double diff = s1.weight - s2.weight;
		int result = diff > Math.pow(10, -8) ? 1 : (diff < -Math.pow(10, -8) ? -1 : 0);
		
		if (result == 0)
			result = s1.root > s2.root ? 1 : (s1.root < s2.root ? -1 : 0);
		if (result == 0)
			for (int i = 0; i < s1.child.length; i++) {
				result = s1.child[i] > s2.child[i] ? 1 : (s1.child[i] < s2.child[i] ? -1 : 0);
				if (result != 0)
					break;
			}
		return result;
	}
}