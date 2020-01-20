//#include "centThread.h"
#include "centInt.h"
#include "cent.h"
#include "pll.h"
#include "test.h"
#include "UsedSt.h"
using namespace std;
Mode CodeMode;

int main(int argc, char*argv[]) {
	string testFile,goalFile;
	CodeMode = RUN_PLLD;
    int runcase = 200;
    int threadnum = 1;
	for (int i = 1; i < argc; i++)
		if (argv[i][0] == '-'){
			if (!strcmp("-PD",argv[i])){          
				CodeMode = RUN_PLLD;
                if (i + 1 < argc)
    			    testFile = argv[i+1];
                else{
                    cout<<"wrong parameter"<<endl;
                    return 0;
                }
            }
			if (!strcmp("-PR",argv[i])){
				CodeMode = RUN_PLLR;
                if (i + 1 < argc)
    			    testFile = argv[i+1];
                else{
                    cout<<"wrong parameter"<<endl;
                    return 0;
                }
            }
            if (!strcmp("-GK",argv[i])){
				CodeMode = GENKG;
                if (i + 2 < argc){
    			    testFile = argv[i+1];
                    goalFile = argv[i+2];
                }
                else{
                    cout<<"wrong parameter!"<<endl;
                    return 0;
                }
            }
            if (!strcmp("-GR",argv[i])){
				CodeMode = GENRD;
                if (i + 2 < argc){
    			    testFile = argv[i+1];
                    goalFile = argv[i+2];
                }
                else{
                    cout<<"wrong parameter!"<<endl;
                    return 0;
                }
            }
            if (!strcmp("-T", argv[i]))
                if (i + 1 < argc)
                    runcase = atoi(argv[i + 1]);
                else{
                    cout << "wrong parameter!" << endl;
                    return 0;
                }
            if (!strcmp("-PP", argv[i])){
                CodeMode = RUN_PLLPAL;
                if (i + 1 < argc)
                    threadnum = atoi(argv[i + 1]);
                else{
                    cout << "wrong parameter!" << endl;
                    return 0;
                }
            }
         }
	if (CodeMode == GENKG) {
		pll p1 = pll(testFile);
		p1.graphKGGen(goalFile);
	}

	if (CodeMode == GENRD) {
		pll p1 = pll(testFile);
		p1.graphRDGen(goalFile);
	}
    if (CodeMode == RUN_PLLD) {
        centInt c1 = centInt();
		c1.deal(testFile, true, runcase);
	}
    if (CodeMode == RUN_PLLR) {
        centInt c1 = centInt();
		c1.deal(testFile, false, runcase);
	}
	if (CodeMode == RUN_PLLPAL) {
        //centThread c1;
		//c1.deal(testFile, true, runcase, threadnum);
	}/*
	testFile = "test.txt";
	test t1 = test(testFile);
	t1.graphGen();
	t1.testFloyd();
	cent c2 = cent();
	c2.deal(testFile);
	t1.testFile(testFile+".floyd" ,testFile+".pl");
	return 0;
	*/
}
