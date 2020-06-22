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
    fgr.global_routing_ver3();
    //fgr.global_routing_ver2();
    cout << "OK" << endl;
    cout << "initial cost = " << fixed << setprecision(0) << fgr.compute_TDM_cost() << endl;

    for (int i = 0; i < 3; i++)
    {
        cout << "iter " << i << " : ";
        fgr.max_subpath_RR();
        cout << "reroute cost = " << fgr.compute_TDM_cost() << endl;
    }

    getchar();

    cout << "round 1 cost = ";
    fgr.total_cost = fgr.compute_TDM_cost();
    cout << fixed << setprecision(0) << fgr.total_cost;
    cout << ", runtime = " << fixed << setprecision(2) << (double)(clock() - t1) / (double)CLOCKS_PER_SEC << " seconds\n";

    //getchar();

    double start = fgr.total_cost;
    double best = start;
    int best_round = 1;

    for (int i = 1; i < multi_round; i++)
    {
        cout << "round " << i + 1 << " cost = ";
        int count = 0;
        fgr.initial_route_result();
        fgr.global_routing_ver3();
        fgr.total_cost = fgr.compute_TDM_cost();
        cout << fixed << setprecision(0) << fgr.total_cost;
        cout << ", runtime = " << fixed << setprecision(2) << (double)(clock() - t1) / (double)CLOCKS_PER_SEC << " seconds\n";

        if (fgr.total_cost < best)
        {
            best = fgr.total_cost;
            best_round = i + 1;
        }

        fgr.round++;
    }

    double end = fgr.total_cost;
    fgr.check_result();
    cout << "runtime = " << fixed << setprecision(2) << (double)(clock() - t1) / (double)CLOCKS_PER_SEC << " seconds\n";
    cout << "best round = " << best_round << endl;
    cout << "improve = " << fixed << setprecision(2) << (start - best) / start * 100 << "%" << endl;
    fgr.output_file(output, t1);
}