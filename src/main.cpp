#include "fpga_gr.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string.h>

using namespace std;

int main(int argc, char **argv)
{
     FPGA_Gr fgr;
     auto t1 = clock();
     char *k = argv[1];
     char num[20];
     char output[100] = "../output/result_";
     strcpy(num, k);

     char f11[100] = "../../benchmark/testcase_new/sy";
     char f12[20] = "r.csv";
     strcat(f11, num);
     strcat(f11, f12);

     char f21[100] = "../../benchmark/testcase_new/syn";
     char f22[20] = ".csv";
     strcat(f21, num);
     strcat(f21, f22);

     strcat(output, num);
     strcat(output, ".out");

     //test random
     int multi_round = 3; //跑幾次init route(含第一次)

     fgr.getfile(f11, f21);
     cout << "Loading files : " << f11 << endl;
     cout << "Loading files : " << f21 << endl;

     fgr.breakdown();
     cout << "construct path table...";

     /*---------construct path table---------*/
     fgr.construct_table_ver2();
     cout << "initial routing...";

     /*---------global routing---------*/
     auto init_time = clock();
     fgr.global_routing_ver3();
     //fgr.global_routing_ver2();
     cout << "OK" << endl;
     double init_cost = fgr.compute_TDM_cost();
     double initt, rrt;
     initt = (double)(clock() - init_time) / (double)CLOCKS_PER_SEC;
     cout << "initial cost = " << fixed << setprecision(0) << init_cost << ", avg TDM = " << fgr.avg_tdm_ratio << ", "
          << ", time = " << fixed << setprecision(2)
          << initt << " seconds" << endl;

     //fgr.show_congestion_map();

     /*---------Rip up and reroute---------*/
     double old_cost = init_cost;
     for (int i = 0; i < 5; i++)
     {
          cout << "iter " << i + 1 << " : ";
          //fgr.update_history_cost();
          auto rrtime = clock();
          fgr.congestion_RR();
          //fgr.subtree_sink_RR();
          //fgr.max_subpath_RR();
          double rr_cost = fgr.compute_TDM_cost();
          fgr.set_after_conj_cost();
          rrt = (double)(clock() - rrtime) / (double)CLOCKS_PER_SEC;

          cout << "\treroute cost = " << fixed << setprecision(0) << rr_cost
               << "\n\ttime = " << fixed << setprecision(2) << rrt << " seconds" << "\n\tavg TDM = " << fgr.avg_tdm_ratio
               << "\n\timprove = " << fixed << setprecision(2)
               << (old_cost - rr_cost) / old_cost * 100 << "%" << endl;

          old_cost = rr_cost;          
     }

     fgr.check_result();
     cout << "runtime = " << fixed << setprecision(2) << initt + rrt << " seconds\n";
     fgr.output_file(output, t1);
}