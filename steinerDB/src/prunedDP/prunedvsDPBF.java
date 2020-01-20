package prunedDP;

import java.io.*;
import java.util.*;

public class prunedvsDPBF {
	
	class inst implements Comparable<inst>{
		int v;
		int X;
		
		inst(int v, int X){
			this.v = v;
			this.X = X;
		}
		public int compareTo(inst r2) {		
			int result = this.v > r2.v ? 1 : (this.v < r2.v ? -1 : 0);		
			if (result == 0)
				result = (this.X < r2.X) ? 1 : (this.X == r2.X ? 0 : -1);
			return result;
		}
	}
	
	void testDeal() throws FileNotFoundException{
		Map<inst, Double> hashm = new TreeMap<inst, Double>();
		Scanner input1 = new Scanner(new File("E:\\data\\test.txt"));
		while (input1.hasNext()) {
			inst i = new inst(input1.nextInt(), input1.nextInt());			
			double wei = input1.nextDouble();
			hashm.put(i, wei);			
		}
		input1.close();
		System.out.println("start");
		Scanner input2 = new Scanner(new File("E:\\data\\testpr.txt"));
		while (input2.hasNext()) {
			inst i = new inst(input2.nextInt(), input2.nextInt());			
			double wei = input2.nextDouble();
			if (!hashm.containsKey(i)) {
				//System.out.println(i.v + " " + i.X + " " + wei);
				continue;
			}
			if (hashm.get(i) != wei)
				System.out.println(i.v + " " + i.X + " " + wei + " " + hashm.get(i));			
		}
		input2.close();
	}
	
	public static void main(String[] args) throws IOException {
		prunedvsDPBF t1 = new prunedvsDPBF();
		t1.testDeal();	
	}
}
