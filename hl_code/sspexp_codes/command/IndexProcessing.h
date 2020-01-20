#ifndef _COMMAND_INDEXING_PROCESSING
#define _COMMAND_INDEXING_PROCESSING


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
    class IndexProcessing: public Command{
    public:
        void exit_with_help(){
            printf("Usage:\n");
            printf("\tsspexp_run -x -d [directedGraphFlag] -w [weightedGraphFlag] -s [specialFlag] \n -m [orderingSchemes] [-a [s_beta]] -g [graphFileName] -e [exportLabelFileName] \n");
            printf("-------------------------------------------------------------------\n");
            printf("Parameter explanation:\n");
            printf("\t[directedGraphFlag]: 0 or 1, for undirected and directed graphs, default is 0\n");
            printf("\t[weightedGraphFlag] = 0 or 1, for unweighted and weighted graphs, default is 0\n");
            printf("\t[specialFlag] = 0: default label\n \t\t\t1: path label\n \t\t\t2: bp label\n \t\t\t3: HLC label\n \t\t\t4: HLCM label\n");
            printf("\t[orderingSchemes] = 0: DHP \n \t\t\t1: BHP\n \t\t\t2: SHP\n");
            printf("\t[s_beta]: [Optional]trade-off parameter for BHP, default is 1\n");
            
            printf("-------------------------------------------------------------------\n");
            printf("Examples:\n");
            printf("Indexing SHP with BP optimization for directed unweighted graph a.txt, outputing alabel.label and alabel.order\n");
            printf("\tsspexp_run -x -d 1 -w 0 -s 2 -m 2 -g a.txt -e alabel \n\n");
        
            printf("Indexing BHP for undirected weighted graph r.txt with s_\beta=1.5, outputing rlabel.label and rlabel.order\n");
            printf("\tsspexp_run -x -d 0 -w 1 -s 0 -m 1 -a 1.5 -g r.txt -e rlabel \n");
            printf("-------------------------------------------------------------------\n");
            exit(1);
        }
        int main(int argc, char *argv[])
        {
            char graphFileName[255] = "";
            char labelFileName[255] = "";
            int t_directed_flag = 0;
            int t_weighted_flag = 0;
            int t_special_flag = 0;
            int t_ordering_flag = 0;
            double beta = 1;
            
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
                    case 'd':
                        t_directed_flag = atoi(argv[i]);
                        break;
                    case 'w':
                        t_weighted_flag = atoi(argv[i]);
                        break;
                    case 's':
                        t_special_flag = atoi(argv[i]);
                        break;
                    case 'm':
                        t_ordering_flag = atoi(argv[i]);
                        break;
                    case 'e':
                        strcpy(labelFileName, argv[i]);
                        break;
                    case 'a':
                        beta = atof(argv[i]);
                        if(argc < 16)
                            exit_with_help();
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
            if(t_special_flag <0 || t_special_flag > 4)
                exit_with_help();
            if(t_ordering_flag <0 || t_ordering_flag > 2)
                exit_with_help();
            if(beta<=0)
                exit_with_help();
            
            Graph graph;
            WGraph wgraph;
            if(WEIGHTED_FLAG == true){
                wgraph.load_graph(graphFileName);
            }else{
                graph.load_graph(graphFileName);
            }
            cout << numOfVertices << " nodes and " << numOfEdges << " arcs " << endl;
            
            if (numOfVertices == 0 || numOfEdges == 0){
                cout << "Corruptted graph file" << endl;
                return 0;
            }
            
            // Indexing
            if(t_special_flag == 0){ // default labels
                if (DIRECTED_FLAG == true){
                    if( WEIGHTED_FLAG == true){
                        if(t_ordering_flag == 0){
                            double _labeling_time = GetCurrentTimeSec();
                            Degree_Ordering degree_order(wgraph);
                            PL_W pl_w(wgraph, degree_order, DIRECTED_FLAG);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            degree_order.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            pl_w.dlabels.save_labels(labelFile.c_str());
                            return 0;
                        }else if (t_ordering_flag == 1){
                            double _labeling_time = GetCurrentTimeSec();
                            Betweenness_Ordering betweenness_ordering(16, beta, wgraph, numOfVertices);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            betweenness_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            betweenness_ordering.dlabels.save_labels(labelFile.c_str());
                            return 0;
                        }else{ // t_ordering_flag == 2
                            double _labeling_time = GetCurrentTimeSec();
                            Coverage_Ordering coverage_ordering(wgraph, DIRECTED_FLAG, true);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            coverage_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            coverage_ordering.dlabels.save_labels(labelFile.c_str());
                            return 0;
                        }
                    }else{
                        if(t_ordering_flag == 0){
                            double _labeling_time = GetCurrentTimeSec();
                            Degree_Ordering degree_order(graph);
                            PL pl(graph, degree_order, DIRECTED_FLAG);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            degree_order.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            pl.dlabels.save_labels(labelFile.c_str());
                            return 0;
                        }else if (t_ordering_flag == 1){
                            double _labeling_time = GetCurrentTimeSec();
                            Betweenness_Ordering betweenness_ordering(16, beta, graph, numOfVertices);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            betweenness_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            betweenness_ordering.dlabels.save_labels(labelFile.c_str());
                            return 0;
                        }else{ // t_ordering_flag == 2
                            double _labeling_time = GetCurrentTimeSec();
                            Coverage_Ordering coverage_ordering(graph, DIRECTED_FLAG);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            coverage_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            coverage_ordering.dlabels.save_labels(labelFile.c_str());
                            return 0;
                        }
                    }
                }else{
                    if( WEIGHTED_FLAG == true){
                        if(t_ordering_flag == 0){
                            double _labeling_time = GetCurrentTimeSec();
                            Degree_Ordering degree_order(wgraph);
                            PL_W pl_w(wgraph, degree_order);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                                                        
                            /*
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            degree_order.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            pl_w.labels.save_labels(labelFile.c_str());
                            */
                            return 0;
                        }else if (t_ordering_flag == 1){
                            double _labeling_time = GetCurrentTimeSec();
                            Betweenness_Ordering betweenness_ordering(16, beta, wgraph, numOfVertices);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;

                            /*
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            betweenness_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            betweenness_ordering.labels.save_labels(labelFile.c_str());
                            */
                            cout << "Avg label: "<<betweenness_ordering.labels.avg_size()<< endl;
                            return 0;
                        }else{ // t_ordering_flag == 2
                            double _labeling_time = GetCurrentTimeSec();
                            Coverage_Ordering coverage_ordering(wgraph);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            /*
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            coverage_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            coverage_ordering.labels.save_labels(labelFile.c_str());
                            */
                            return 0;
                        }
                    }else{
                        if(t_ordering_flag == 0){
                            double _labeling_time = GetCurrentTimeSec();
                            Degree_Ordering degree_order(graph);
                            PL pl(graph, degree_order);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            degree_order.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            pl.labels.save_labels(labelFile.c_str());
                            return 0;
                        }else if (t_ordering_flag == 1){
                            double _labeling_time = GetCurrentTimeSec();
                            Betweenness_Ordering betweenness_ordering(16, beta, graph, numOfVertices);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            betweenness_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            betweenness_ordering.labels.save_labels(labelFile.c_str());
                            return 0;
                        }else{ // t_ordering_flag == 2
                            double _labeling_time = GetCurrentTimeSec();
                            vector<NodeID> border(numOfVertices);
                            Coverage_Ordering coverage_ordering(graph, border);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            coverage_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            coverage_ordering.labels.save_labels(labelFile.c_str());
                            return 0;
                        }
                    }
                }
            } else if(t_special_flag == 1){ // path labels
                if (DIRECTED_FLAG == true){
                    if( WEIGHTED_FLAG == true){
                        if(t_ordering_flag == 0){
                            double _labeling_time = GetCurrentTimeSec();
	                        Degree_Ordering degree_order(wgraph);
                            PL_W pl_w(wgraph, degree_order, true, true,true);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            degree_order.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            pl_w.dplabels.save_labels(labelFile.c_str());
                            return 0;
                        }else if (t_ordering_flag == 1){
                            double _labeling_time = GetCurrentTimeSec();
                            Betweenness_Ordering betweenness_ordering(16, beta, wgraph, numOfVertices, true);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            betweenness_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            betweenness_ordering.dplabels.save_labels(labelFile.c_str());
                            return 0;
                        }else{ // t_ordering_flag == 2
                            double _labeling_time = GetCurrentTimeSec();
                        	Coverage_Ordering_Path coverage_ordering(wgraph, DIRECTED_FLAG);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            coverage_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            coverage_ordering.dplabels.save_labels(labelFile.c_str());
                            return 0;
                        }
                    }else{
                        if(t_ordering_flag == 0){
                            double _labeling_time = GetCurrentTimeSec();
                        	Degree_Ordering degree_order(graph);
                            PL pl(graph, degree_order, true, true);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            degree_order.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            pl.dplabels.save_labels(labelFile.c_str());
                            return 0;
                        }else if (t_ordering_flag == 1){
                            double _labeling_time = GetCurrentTimeSec();
                            Betweenness_Ordering betweenness_ordering(16, beta, graph, numOfVertices, true);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            betweenness_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            betweenness_ordering.dplabels.save_labels(labelFile.c_str());
                            return 0;
                        }else{ // t_ordering_flag == 2
                            double _labeling_time = GetCurrentTimeSec();
                            Coverage_Ordering_Path coverage_ordering(graph, DIRECTED_FLAG);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            coverage_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            coverage_ordering.dplabels.save_labels(labelFile.c_str());
                            return 0;
                        }
                    }
                }else{
                    if( WEIGHTED_FLAG == true){
                        if(t_ordering_flag == 0){
                            double _labeling_time = GetCurrentTimeSec();
                            Degree_Ordering degree_order(wgraph);
                            PL_W pl_w(wgraph, degree_order, true, true);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            degree_order.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            pl_w.plabels.save_labels(labelFile.c_str());
                            return 0;
                        }else if (t_ordering_flag == 1){
                            double _labeling_time = GetCurrentTimeSec();
                            Betweenness_Ordering betweenness_ordering(16, beta, wgraph, numOfVertices, true);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            betweenness_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            betweenness_ordering.plabels.save_labels(labelFile.c_str());
                            return 0;
                        }else{ // t_ordering_flag == 2
                            double _labeling_time = GetCurrentTimeSec();
                            Coverage_Ordering_Path coverage_ordering(wgraph);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            coverage_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            coverage_ordering.plabels.save_labels(labelFile.c_str());
                            return 0;
                        }
                    }else{
                        if(t_ordering_flag == 0){
                            double _labeling_time = GetCurrentTimeSec();
                            Degree_Ordering degree_order(graph);
                            PL pl(graph, degree_order, true, false);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            degree_order.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            pl.plabels.save_labels(labelFile.c_str());
                            return 0;
                        }else if (t_ordering_flag == 1){
                            double _labeling_time = GetCurrentTimeSec();
                            Betweenness_Ordering betweenness_ordering(16, beta, graph, numOfVertices, true);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            betweenness_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            betweenness_ordering.plabels.save_labels(labelFile.c_str());
                            return 0;
                        }else{ // t_ordering_flag == 2
                            double _labeling_time = GetCurrentTimeSec();
                            vector<NodeID> border(numOfVertices);
                            Coverage_Ordering_Path coverage_ordering(graph, border);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            coverage_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            coverage_ordering.plabels.save_labels(labelFile.c_str());
                            return 0;
                        }
                    }
                }
            } else if(t_special_flag == 2){ // BP labels
                if (DIRECTED_FLAG == true){
                    if( WEIGHTED_FLAG == true){
                        cout << "BP does not support weighted graphs" << endl;
                        return 0;
                    }else{
                        if(t_ordering_flag == 0){
                            double _labeling_time = GetCurrentTimeSec();
                            Degree_Ordering degree_order(graph);
                            BPL<50> bpl(graph, degree_order, DIRECTED_FLAG);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            bpl.dbplabels.save_labels(labelFile.c_str());
                            return 0;
                        }else if (t_ordering_flag == 1){
                            double _labeling_time = GetCurrentTimeSec();
                            BP_Betweenness_Ordering<50> betweenness_ordering(16, beta, graph, numOfVertices);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            betweenness_ordering.dbplabels.save_labels(labelFile.c_str());
                            return 0;
                        }else{ // t_ordering_flag == 2
                            double _labeling_time = GetCurrentTimeSec();
                            Coverage_Ordering_BP<50> coverage_ordering(graph, DIRECTED_FLAG);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            coverage_ordering.dbplabels.save_labels(labelFile.c_str());
                            return 0;
                        }
                    }
                }else{
                    if( WEIGHTED_FLAG == true){
                        cout << "BP does not support weighted graphs" << endl;
                        return 0;
                    }else{
                        if(t_ordering_flag == 0){
                            double _labeling_time = GetCurrentTimeSec();
                            Degree_Ordering degree_order(graph);
                            BPL<50> bpl(graph, degree_order);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            bpl.bplabels.save_labels(labelFile.c_str());
                            return 0;
                        }else if (t_ordering_flag == 1){
                            double _labeling_time = GetCurrentTimeSec();
	                        BP_Betweenness_Ordering<50> betweenness_ordering(16, beta, graph, numOfVertices);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            betweenness_ordering.bplabels.save_labels(labelFile.c_str());
                            return 0;
                        }else{ // t_ordering_flag == 2
                            double _labeling_time = GetCurrentTimeSec();
                            vector<NodeID> border(numOfVertices);
                            Coverage_Ordering_BP<50> coverage_ordering(graph, border);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            coverage_ordering.bplabels.save_labels(labelFile.c_str());
                            return 0;
                        }
                    }
                }
            } else if(t_special_flag == 3){ // HLC labels
                if (DIRECTED_FLAG == true){
                    if( WEIGHTED_FLAG == true){
                       cout << "This part of the codes has been removed for code simplicity. It involves no experiments conducted in the paper. We will add this function back upon acceptance." << endl;
                        return 0;
                    }else{
                        if(t_ordering_flag == 0){
                            double _labeling_time = GetCurrentTimeSec();
                            Degree_Ordering degree_order(graph);
	                        CPL pl(graph, degree_order, false, DIRECTED_FLAG);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            degree_order.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            pl.clabels.save_labels_d(labelFile.c_str());
                            return 0;
                        }else if (t_ordering_flag == 1){
                            double _labeling_time = GetCurrentTimeSec();
	                        Betweenness_Ordering betweenness_ordering(16, beta, graph, numOfVertices, 10, 8 , false, true);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            betweenness_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            betweenness_ordering.clabels.save_labels_d(labelFile.c_str());
                            return 0;
                        }else{ // t_ordering_flag == 2
                            double _labeling_time = GetCurrentTimeSec();
                        	Coverage_Ordering_Compress coverage_ordering(graph, true, false); 
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            coverage_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            coverage_ordering.clabels.save_labels_d(labelFile.c_str());
                            return 0;
                        }
                    }
                }else{
                    if( WEIGHTED_FLAG == true){
                        if(t_ordering_flag == 0){
                            double _labeling_time = GetCurrentTimeSec();
                            Degree_Ordering degree_order(wgraph);
                           	CPL_W pl(wgraph, degree_order, false); 
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            degree_order.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            pl.clabels.save_labels(labelFile.c_str());
                            return 0;
                        }else if (t_ordering_flag == 1){
                            double _labeling_time = GetCurrentTimeSec();
                           	Betweenness_Ordering betweenness_ordering(16, beta, wgraph, numOfVertices, 10, 8 , false, false);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            betweenness_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            betweenness_ordering.clabels.save_labels(labelFile.c_str());
                            return 0;
                        }else{ // t_ordering_flag == 2
                            double _labeling_time = GetCurrentTimeSec();
                        	Coverage_Ordering_Compress coverage_ordering(wgraph, false, false);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            coverage_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            coverage_ordering.clabels.save_labels(labelFile.c_str());
                            return 0;
                        }
                    }else{
                        if(t_ordering_flag == 0){
                            double _labeling_time = GetCurrentTimeSec();
                            Degree_Ordering degree_order(graph);
                            CPL pl(graph, degree_order, false); 
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            degree_order.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            pl.clabels.save_labels(labelFile.c_str());
                            return 0;
                        }else if (t_ordering_flag == 1){
                            double _labeling_time = GetCurrentTimeSec();
	                        Betweenness_Ordering betweenness_ordering(16, beta, graph, numOfVertices, 10, 8 , false, false);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            betweenness_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            betweenness_ordering.clabels.save_labels(labelFile.c_str());
                            return 0;
                        }else{ // t_ordering_flag == 2
                            double _labeling_time = GetCurrentTimeSec();
                            vector<NodeID> border(numOfVertices);
                            Coverage_Ordering_Compress coverage_ordering(graph);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            coverage_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            coverage_ordering.clabels.save_labels(labelFile.c_str());
                            return 0;
                        }
                    } 
                }
            } else if(t_special_flag == 4){ // HLCM labels
                if (DIRECTED_FLAG == true){
                    if( WEIGHTED_FLAG == true){
                       cout << "This part of the codes has been removed for code simplicity. It involves no experiments conducted in the paper. We will add this function back upon acceptance." << endl;
                        return 0;
                    }else{
                        if(t_ordering_flag == 0){
                            double _labeling_time = GetCurrentTimeSec();
                            Degree_Ordering degree_order(graph);
	                        CPL pl(graph, degree_order, true, DIRECTED_FLAG);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            degree_order.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            pl.clabels.save_two_level_labels_d(labelFile.c_str());
                            return 0;
                        }else if (t_ordering_flag == 1){
                            double _labeling_time = GetCurrentTimeSec();
	                        Betweenness_Ordering betweenness_ordering(16, beta, graph, numOfVertices, 10, 8 , true, true);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            betweenness_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            betweenness_ordering.clabels.save_two_level_labels_d(labelFile.c_str());
                            return 0;
                        }else{ // t_ordering_flag == 2
                            double _labeling_time = GetCurrentTimeSec();
                        	Coverage_Ordering_Compress coverage_ordering(graph, true, true); 
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            coverage_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            coverage_ordering.clabels.save_two_level_labels_d(labelFile.c_str());
                            return 0;
                        }
                    }
                }else{
                    if( WEIGHTED_FLAG == true){
                        if(t_ordering_flag == 0){
                            double _labeling_time = GetCurrentTimeSec();
                            Degree_Ordering degree_order(wgraph);
                           	CPL_W pl(wgraph, degree_order, true); 
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            degree_order.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            pl.clabels.save_two_level_labels(labelFile.c_str());
                            return 0;
                        }else if (t_ordering_flag == 1){
                            double _labeling_time = GetCurrentTimeSec();
                           	Betweenness_Ordering betweenness_ordering(16, beta, wgraph, numOfVertices, 10, 8 , true, false);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            betweenness_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            betweenness_ordering.clabels.save_two_level_labels(labelFile.c_str());
                            return 0;
                        }else{ // t_ordering_flag == 2
                            double _labeling_time = GetCurrentTimeSec();
                        	Coverage_Ordering_Compress coverage_ordering(wgraph, false, true);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            coverage_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            coverage_ordering.clabels.save_two_level_labels(labelFile.c_str());
                            return 0;
                        }
                    }else{
                        if(t_ordering_flag == 0){
                            double _labeling_time = GetCurrentTimeSec();
                            Degree_Ordering degree_order(graph);
                            CPL pl(graph, degree_order, true); 
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            degree_order.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            pl.clabels.save_two_level_labels(labelFile.c_str());
                            return 0;
                        }else if (t_ordering_flag == 1){
                            double _labeling_time = GetCurrentTimeSec();
	                        Betweenness_Ordering betweenness_ordering(16, beta, graph, numOfVertices, 10, 8 , true, false);
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            betweenness_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            betweenness_ordering.clabels.save_two_level_labels(labelFile.c_str());
                            return 0;
                        }else{ // t_ordering_flag == 2
                            double _labeling_time = GetCurrentTimeSec();
                            vector<NodeID> border(numOfVertices);
                            Coverage_Ordering_Compress coverage_ordering(graph, false, true); 
                            _labeling_time = GetCurrentTimeSec() - _labeling_time;
                            cout << "Indexing time:" << _labeling_time << endl;
                            
                            string orderFileName(labelFileName);
                            orderFileName.append(".order");
                            coverage_ordering.save_rank(orderFileName.c_str());
                            
                            string labelFile(labelFileName);
                            labelFile.append(".label");
                            coverage_ordering.clabels.save_two_level_labels(labelFile.c_str());
                            return 0;
                        }
                    }
                }
            }
        }
    };
}
#endif 
