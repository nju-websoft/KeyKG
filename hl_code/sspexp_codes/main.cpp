// umac_ssp_exp.cpp : Defines the entry point for the console application.
//

//#include "command/ConstructionParadigms.h"
//#include "command/IndexProcessing.h"
//#include "command/QueryProcessing.h"

#include <iostream>
#include <fstream>

#include "command/IndexProcessing.h"
#include "command/ConstructionParadigms.h"
#include "command/QueryProcessing.h"

using namespace std;

void exit_with_help(){
    printf("Command Usage: (the first parameter specifies the procedure be executed)\n");
    printf("-------------------------------------------------------------------\n");
    printf("(1) -x: Indexing by different ordering schemes:\n");
    printf("\tsspexp_run -x -d [directedGraphFlag] -w [weightedGraphFlag] -s [specialFlag] -m [orderingSchemes] [-a [s_beta]] -g [graphFileName] -e [exportLabelFileName] \n");
    printf("-------------------------------------------------------------------\n");
    printf("(2) -q: Query testing for a label file:\n");
    printf("\tsspexp_run -q -d [directedGraphFlag] -s [specialFlag] [-o [orderFileName]] -l [inputLabelFileName] -n [queryNum] \n");
    printf("-------------------------------------------------------------------\n");
    printf("(3) -y: Indexing by different construction paradigms (given ordering file):\n");
    printf("\tsspexp_run -y -d [directedGraphFlag] -w [weightedGraphFlag] -m [consturctionParadigms] -g [graphFileName] -o [orderFileName] -e [exportLabelFileName] \n");
    printf("-------------------------------------------------------------------\n");
    printf("Parameter explanation:\n");
    printf("\t[directedGraphFlag]: 0 or 1, for undirected and directed graphs, default is 0\n");
    printf("\t[weightedGraphFlag] = 0 or 1, for unweighted and weighted graphs, default is 0\n");
    printf("\t[specialFlag] = 0: default label\n \t\t\t1: path label\n \t\t\t2: bp label\n \t\t\t3: HLC label\n \t\t\t4: HLCM label\n");
    printf("\t[orderingSchemes] = 0: DHP \n \t\t\t1: BHP\n \t\t\t2: SHP\n");
    printf("\t[s_beta]: [Optional]trade-off parameter for BHP, default is 1\n");
    printf("\t[queryNum]: the number of randomly generated non-duplicated reachable queries. Extra half of queryNum will be generated for warm-up\n");
    printf("\t[orderFileName]: only required for path query processing\n");
    printf("\t[consturctionParadigms] = 0: Hub Pushing Algorithm \n \t\t\t\t1: Hub Pulling Algorithm\n");

    
    printf("-------------------------------------------------------------------\n");
    printf("Examples:\n");
    printf("Indexing SHP with BP optimization for directed unweighted graph a.txt, outputing alabel.label and alabel.order\n");
    printf("\tsspexp_run -x -d 1 -w 0 -s 2 -m 2 -g a.txt -e alabel \n");
    printf("Indexing Hub Pushing Algorithm given odering file alabel.order for directed unweighted graph a.txt, outputing alabel.label\n");
    printf("\tsspexp_run -y -d 1 -w 0 -m 0 -g a.txt -o alabel.order -e alabel \n");
    printf("Query processing for 1000000 random queries for directed path label alabel_path.label\n");
    printf("\tsspexp_run -q -d 1 -s 1 -l alabel_path.label -n 1000000 \n");
    printf("Query processing for 1000000 random queries for directed HLC label alabel_hlc.label with alabel_hlc.order\n");
    printf("\tsspexp_run -q -d 1 -s 1 -o alabel_hlc.order -l alabel_hlc.label -n 1000000 \n");
    printf("-------------------------------------------------------------------\n");
    exit(1);
}
                                
    /** The main program. */
                                
    int main(int argc, char *argv[]){
        // The program is controlled by command-line arguments. The order of those
        // arguments is important. The first argument specifies the Command-
        // class that is used.
        
        int opt = 'm';
        if (argc > 1){
            switch (argv[1][1]){
                case 'x':
                    opt = 'x';
                    break;
                case 'q':
                    opt = 'q';
                    break;
                case 'y':
                    opt = 'y';
                    break;
                    
                default:
                    exit_with_help();
                    break;
            }
        }
        
        if (argc == 1)
            exit_with_help();
        
        Command* m = NULL;
        int result = 0;
        switch (opt) {
                
            case 'x':
                m = new command::IndexProcessing();
                break;
            case 'y':
                m = new command::ConstructionParadigms();
                break;
            case 'q':
                m = new command::QueryProcessing();
                break;
            default:
                break;
                
        }
        if (m != NULL) {
            result = m->main(argc, argv);
            delete m;
        } else {
            exit_with_help();
        }
        return result;
        
        return 0;
    }

