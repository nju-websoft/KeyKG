package DOGST;

import java.util.*;


public class AnsTree  implements Comparable<AnsTree>{
	public int root;
	public Set<Integer> terminal;
	public ArrayList<TreeEdge> edge;	//record tree edge
	public double weight;
	
	AnsTree() {
		this.terminal = new HashSet<Integer>();
		edge = new ArrayList<TreeEdge>(); 
	}
	
	public int compareTo(AnsTree s2) {		
		int result = this.weight > s2.weight ? 1 : (this.weight < s2.weight ? -1 : 0);
		if (result == 0)
			result = (this.root < s2.root) ? 1 : (this.root == s2.root ? 0 : -1);
		return result;
	}
	
	//test whether this ansTree has the same terminal with te2
	public boolean isSetEqual(AnsTree te2) {
		Set<Integer> set1 = this.terminal;
		Set<Integer> set2 = te2.terminal;

		if (set1 == null && set2 == null) {
			return true; // Both are null
		}

		if (set1 == null || set2 == null || set1.size() != set2.size() || set1.size() == 0 || set2.size() == 0) {
			return false;
		}

		Iterator<Integer> it2 = set2.iterator();

		boolean isFullEqual = true;

		while (it2.hasNext()) {
			if (!set1.contains(it2.next())) {
				isFullEqual = false;
				break;
			}
		}

		return isFullEqual;
	}
	
}
