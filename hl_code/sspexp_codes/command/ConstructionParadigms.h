#ifndef _COMMAND_CONSTRUCTION_PARADIGMS
#define _COMMAND_CONSTRUCTION_PARADIGMS

#include "../src/graph.h"
#include "../src/ordering.h"
#include "../src/construction.h"
#include "../src/coverage_ordering.h"
#include "../src/coverage_ordering_path.h"
#include "../src/coverage_ordering_bp.h"
#include "../src/coverage_ordering_compress.h"
#include "../src/labels.h"
#include "../src/time_util.h"
#include "../command.h"
#include<cstring>

using namespace time_util;

namespace command{
    class ConstructionParadigms: public Command{
    public:
        void exit_with_help(){
            printf("Usage:\n");
            printf("\tsspexp_run -y -d [directedGraphFlag] -w [weightedGraphFlag] -m [consturctionParadigms] -g [graphFileName] -o [orderFileName] -e [exportLabelFileName] \n");
            printf("-------------------------------------------------------------------\n");
            printf("Parameter explanation:\n");
            printf("\t[directedGraphFlag]: 0 or 1, for undirected and directed graphs, default is 0\n");
            printf("\t[weightedGraphFlag] = 0 or 1, for unweighted and weighted graphs, default is 0\n");
            printf("\t[specialFlag] = 0: default label\n \t\t\t1: path label\n \t\t\t2: bp label\n \t\t\t3: HLC label\n \t\t\t4: HLCM label\n");
            printf("\t[consturctionParadigms] = 0: Hub Pushing Algorithm \n \t\t\t\t1: Hub Pulling Algorithm\n");
            printf("-------------------------------------------------------------------\n");
            printf("Examples:\n");
            printf("Indexing Hub Pushing Algorithm given odering file alabel.order for directed unweighted graph a.txt, outputing alabel.label\n");
            printf("\tsspexp_run -y -d 1 -w 0 -m 0 -g a.txt -o alabel.order -e alabel \n");
            printf("-------------------------------------------------------------------\n");
            exit(1);
        }
        int main(int argc, char *argv[])
        {
            char graphFileName[255] = "";
            char labelFileName[255] = "";
            char orderFileName[255] = "";
            int t_directed_flag = 0;
            int t_weighted_flag = 0;
            int t_construction_flag = 0;
            double beta = 1;
            DIRECTED_FLAG = false;
            WEIGHTED_FLAG = false;
            
            if(argc < 14)
                exit_with_help();
            
            for(int i = 2; i < argc; i++){
                if(argv[i][0] != '-') break;
                if(++i >= argc)
                    exit_with_help();
                switch (argv[i-1][1]){
                    case 'g':
                        strcpy(graphFileName,argv[i]);
                        break;
                    case 'o':
                        strcpy(orderFileName,argv[i]);
                        break;
                    case 'd':
                        t_directed_flag = atoi(argv[i]);
                        break;
                    case 'w':
                        t_weighted_flag = atoi(argv[i]);
                        break;
                    case 'm':
                        t_construction_flag = atoi(argv[i]);
                        break;
                    case 'e':
                        strcpy(labelFileName, argv[i]);
                        break;
                    default:
                        exit_with_help();
                }
            }
            
            if (t_directed_flag == 1)
                DIRECTED_FLAG = true;
            if (t_weighted_flag == 1)
                WEIGHTED_FLAG = true;
            
            if(t_directed_flag != 1 && t_directed_flag != 0)
                exit_with_help();
            if(t_weighted_flag != 1 && t_weighted_flag != 0)
                exit_with_help();
            if(t_construction_flag <0 || t_construction_flag > 2)
                exit_with_help();
            
            Graph graph;
            WGraph wgraph;
            CHGraph chgraph;
            if(t_construction_flag == 0){
                if(WEIGHTED_FLAG == true){
                    wgraph.load_graph(graphFileName);
                }else{
                    graph.load_graph(graphFileName);
                }
            }else if (t_construction_flag == 1){
                if(WEIGHTED_FLAG == true){
                    chgraph.load_wgraph(argv[1]);  
                }else{
                    chgraph.load_graph(argv[1]);  
                }
            }
            cout << numOfVertices << " nodes and " << numOfEdges << " arcs " << endl;
            
            if (numOfVertices == 0 || numOfEdges == 0){
                cout << "Corruptted graph file" << endl;
                return 0;
            }
            
            // Indexing
            if (DIRECTED_FLAG == true){
                if( WEIGHTED_FLAG == true){
                    if(t_construction_flag == 0){
                        double _labeling_time = GetCurrentTimeSec();
                        Given_Ordering given_order(orderFileName, wgraph);
                    	PL_W pl_w(wgraph, given_order, DIRECTED_FLAG);
                        _labeling_time = GetCurrentTimeSec() - _labeling_time;
                        cout << "Indexing time:" << _labeling_time << endl;
                        
                        string labelFile(labelFileName);
                        labelFile.append(".label");
                        pl_w.dlabels.save_labels(labelFile.c_str());
                        return 0;
                    }else if (t_construction_flag == 1){
                        double _labeling_time = GetCurrentTimeSec();
                        Given_Ordering given_order(orderFileName);
                        double _contracting_time = 0; 
                        Bottomup bu(chgraph, given_order, _contracting_time);
                        _labeling_time = GetCurrentTimeSec() - _labeling_time;
                        cout << "Indexing time:" << _labeling_time << endl;
                        cout << "Contracting time:" << _contracting_time << endl;
                        
                        string labelFile(labelFileName);
                        labelFile.append(".label");
                        bu.dlabels.save_labels(labelFile.c_str());
                        return 0;
                    }
                }else{
                    if(t_construction_flag == 0){
                        double _labeling_time = GetCurrentTimeSec();
                        Given_Ordering given_order(orderFileName, graph);
                    	PL pl(graph, given_order, DIRECTED_FLAG);
                        _labeling_time = GetCurrentTimeSec() - _labeling_time;
                        cout << "Indexing time:" << _labeling_time << endl;
                        
                        string labelFile(labelFileName);
                        labelFile.append(".label");
                        pl.dlabels.save_labels(labelFile.c_str());
                        return 0;
                    }else if (t_construction_flag == 1){
                        double _labeling_time = GetCurrentTimeSec();
                        Given_Ordering given_order(orderFileName);
                        double _contracting_time = 0; 
                        Bottomup bu(chgraph, given_order, _contracting_time);
                        _labeling_time = GetCurrentTimeSec() - _labeling_time;
                        cout << "Indexing time:" << _labeling_time << endl;
                        cout << "Contracting time:" << _contracting_time << endl;
                        
                        string labelFile(labelFileName);
                        labelFile.append(".label");
                        bu.dlabels.save_labels(labelFile.c_str());
                        return 0;
                    }
                }
            }else{
                if( WEIGHTED_FLAG == true){
                    if(t_construction_flag == 0){
                        double _labeling_time = GetCurrentTimeSec();
                        Given_Ordering given_order(orderFileName, wgraph);
                    	PL_W pl_w(wgraph, given_order);
                        _labeling_time = GetCurrentTimeSec() - _labeling_time;
                        cout << "Indexing time:" << _labeling_time << endl;
                        
                        string labelFile(labelFileName);
                        labelFile.append(".label");
                        pl_w.labels.save_labels(labelFile.c_str());
                        return 0;
                    }else if (t_construction_flag == 1){
                        double _labeling_time = GetCurrentTimeSec();
                        Given_Ordering given_order(orderFileName);
                        double _contracting_time = 0; 
                        Bottomup bu(chgraph, given_order, _contracting_time);
                        _labeling_time = GetCurrentTimeSec() - _labeling_time;
                        cout << "Indexing time:" << _labeling_time << endl;
                        cout << "Contracting time:" << _contracting_time << endl;
                        
                        string labelFile(labelFileName);
                        labelFile.append(".label");
                        bu.labels.save_labels(labelFile.c_str());
                        return 0;
                    }
                }else{
                    if(t_construction_flag == 0){
                        double _labeling_time = GetCurrentTimeSec();
                        Given_Ordering given_order(orderFileName, graph);
                    	PL pl(graph, given_order);
                        _labeling_time = GetCurrentTimeSec() - _labeling_time;
                        cout << "Indexing time:" << _labeling_time << endl;
                        
                        string labelFile(labelFileName);
                        labelFile.append(".label");
                        pl.labels.save_labels(labelFile.c_str());
                        return 0;
                    }else if (t_construction_flag == 1){
                        double _labeling_time = GetCurrentTimeSec();
                        Given_Ordering given_order(orderFileName);
                        double _contracting_time = 0; 
                        Bottomup bu(chgraph, given_order, _contracting_time);
                        _labeling_time = GetCurrentTimeSec() - _labeling_time;
                        cout << "Indexing time:" << _labeling_time << endl;
                        cout << "Contracting time:" << _contracting_time << endl;
                        
                        string labelFile(labelFileName);
                        labelFile.append(".label");
                        bu.labels.save_labels(labelFile.c_str());
                        return 0;
                    }
                }
            }
             
        }
    };
}
#endif 

