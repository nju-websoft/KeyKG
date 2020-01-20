package PLL;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

public class TestPLL {
	int nodeNum = 400;
	double[][] distance;
	WeightedPLL ww;
	
	void graphGenerate(String st) throws IOException {
		double [][]Matrix = new double[nodeNum][nodeNum];
		for(int i=0; i<nodeNum; i++)
			for(int j=0; j<nodeNum; j++)
				Matrix[i][j] = Double.MAX_VALUE;
		Random r1 = new Random();
		for(int i=0; i<nodeNum; i++)
			for (int j = 0; j < nodeNum/10; j++) {
				double x1 = r1.nextDouble()*100;
				int s = (r1.nextInt()%nodeNum+nodeNum)%nodeNum;
				Matrix[i][s] = x1;
				Matrix[s][i] = x1;
			}
		File fileG = new File(st+".regular");
		if (!fileG.exists()) { // if file doesn't exists, then create it
			fileG.createNewFile();
		}
		File fileW = new File(st+".weight");
		if (!fileW.exists()) { // if file doesn't exists, then create it
			fileW.createNewFile();
		}
		PrintWriter fopG = new PrintWriter(fileG);
		PrintWriter fopW = new PrintWriter(fileW);
		fopG.println(nodeNum);
		for(int i=0; i<nodeNum; i++)
			for(int j=i+1; j<nodeNum; j++)
				if (Matrix[i][j]<Double.MAX_VALUE-1) {
					fopG.println(i + " " + j);
					fopW.println(Matrix[i][j]);
				}
		fopG.close();
		fopW.close();
	}
	
	public void floyd(String st) throws FileNotFoundException{
		distance = new double[nodeNum][nodeNum];
	
		for(int i=0; i<nodeNum; i++)
			for(int j=0; j<nodeNum; j++)
				distance[i][j] = Double.MAX_VALUE;
		for(int i=0; i<nodeNum; i++)
			distance[i][i]=0;
		File regularFile = new File(st+".regular");
		Scanner inputR = new Scanner(regularFile);
		
		File weightFile = new File(st+".weight");
		Scanner inputW = new Scanner(weightFile);
		int x, y;
		double w;
		inputR.nextInt();
		while (inputR.hasNext()) {
			x = inputR.nextInt();
			y = inputR.nextInt();
			w = inputW.nextDouble();
			distance[x][y] = w;
			distance[y][x] = w;
		}
		inputR.close();
		inputW.close();
		
		//循环更新矩阵的值
		for(int k=0; k<nodeNum; k++){
			for(int i=0; i<nodeNum; i++){
				for(int j=0; j<nodeNum; j++){
					double temp = (distance[i][k] >= Double.MAX_VALUE-1 || distance[k][j] >=Double.MAX_VALUE-1) 
							? Double.MAX_VALUE : distance[i][k] + distance[k][j];
					if(distance[i][j] > temp){
						distance[i][j] = temp;
					}
				}
			}
		}
					
	}
	
	void testDis() {		
		for (int i = 0; i < nodeNum; i++)
			for (int j = 0; j < nodeNum; j++)
				if(Math.abs(distance[i][j]-ww.getPath(i,j))>Math.pow(10, -3))
					System.out.println(i+" "+j+" "+distance[i][j]+" "+ww.query(i,j));
	}
	void testDeal(String st) throws IOException {
		graphGenerate(st);
		floyd(st);		
	}
	static String testFile = "E:\\data\\test.txt";
	public static void main(String[] args) throws IOException {
		TestPLL t1 = new TestPLL();
		t1.testDeal(testFile);
	}
}
