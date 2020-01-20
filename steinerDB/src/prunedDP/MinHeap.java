package prunedDP;

import java.util.ArrayList;
import java.util.Map;
import java.util.TreeMap;

public class MinHeap {
	
	ArrayList<CandiDTree> heap;	//start from 1
	Map<DTree, Integer>heapLoc;	
	int heapSize = 0;
	MinHeap(){		
		heap = new ArrayList<CandiDTree>();
		heap.add(null);		
		heapLoc = new TreeMap<DTree, Integer>();
	}
	
	//swap two heap element
	void swap(int i, int j) {
		CandiDTree tmp = heap.get(i);
		heap.set(i, heap.get(j));
		heap.set(j, tmp);
	}
	
	void push(CandiDTree ct) {		
		if (!heapLoc.containsKey(ct.dt)) {
			heap.add(ct);
			heapSize++;
			DTree dt = new DTree(ct.v, ct.X);
			heapLoc.put(dt, heapSize);
			upModify(heapSize);
			return;
		}		
		int loc = heapLoc.get(ct.dt);
		if (loc < 1) return;
		if (heap.get(loc).lb <= ct.lb) return;
		heap.set(loc, ct);
		upModify(loc);		
		
	}
		
	CandiDTree poll() {
		if (heapSize == 0) return null;
		CandiDTree ans = heap.get(1);		
		heapLoc.put(ans.dt, 0);
		heapLoc.put(heap.get(heapSize).dt, 1);		 
		swap(1, heapSize);
		heap.remove(heapSize);
		heapSize--;
		downModify(1);		
		return ans;
	}
	
	
	boolean isEmpty() {
		if (heapSize == 0) return true;
		return false;
	}	
		
	void downModify(int loc) {
		int target = loc;
		int left = loc * 2;
		int right = loc * 2 + 1;	
		if (left <= heapSize && heap.get(left).lb < heap.get(target).lb)			
			target = left;
		if (right <= heapSize && heap.get(right).lb < heap.get(target).lb)			
			target = right;
		if (target != loc) {					
			heapLoc.put(heap.get(target).dt, loc);
			heapLoc.put(heap.get(loc).dt, target);		
			swap(target, loc);
			downModify(target);
		}
	}
	
	void upModify(int loc) {
		if (loc > 1) {
			int target = loc / 2;
			if (heap.get(target).lb > heap.get(loc).lb) {
				heapLoc.put(heap.get(target).dt, loc);
				heapLoc.put(heap.get(loc).dt, target);				
				swap(target, loc);
				upModify(target);
			}
		}
	}	
}
