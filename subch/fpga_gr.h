#ifndef FPGA_GR_HPP
#define FPGA_GR_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <deque>
#include <queue>
#include <map>
#include <cmath>
#include <climits>
#include <iomanip>
#include <time.h>
#include "node.h"

#define LIMIT_HOP 1

using namespace std;

class Channel;

class FPGA
{
public:
    int id;
    vector<pair<int, int>> nbr_pair; //nbr_fpga pairs
    map<int, int> nbr_index;         //input nbr fpga id --> output fpga index in vector nbr_pair
};

class Sink
{
public:
    int id;
    int weight;
    Sink() {}
    ~Sink() {}
};

class SubNet
{
public:
    int parent_net; //record parent net id
    int source, sink;
    int weight;
    vector<int> path;
    //map<vector<int>, >
};

class Ch_nets
{
public:
    Channel *ch; //point to channel
    //double tdm;
    //double edge_weight;
};

class Net
{
public:
    int id;
    string name;
    int source;
    double cost;
    bool sorted;
    vector<Sink> sink;
    vector<SubNet> sbnet;
    vector<pair<int, int>> channels;
    Tree_Node *rtree_root; //routing tree root
    map<pair<int, int>, int> edge_crit;
    int total_tree_edge;
    double signal_weight;

    Net()
    {
        sorted = false;
        rtree_root = NULL;
        total_tree_edge = 0;
    }
};

class Channel
{
public:
    pair<int, int> name;
    list<pair<Net *, double>> net_ch_weight; //net list and edge weight
    double history_used[2];                  // min-->max : index=0
    int capacity;
    vector<Net *> net_sigweight[2];

    Channel()
    {
        capacity = 0;
        history_used[0] = history_used[1] = 0.0;
    }

    ~Channel();
};

class Table_content
{
public:
    int hops;
    vector<int> parent;
    Table_content()
    {
        hops = 0;
    }
};

class Path_table_ver2
{
public:
    vector<Table_content> cand;
};

class FPGA_Gr
{
public:
    int round;
    int fpga_num;
    int sink_num;
    int total_demand;
    double total_cost, avg_sk_weight;
    int avg_tdm_ratio;
    int maxtdm, mintdm;
    double maxsgw, minsgw; // max and min signal weight
    vector<FPGA> fpga;
    vector<Net> net;
    vector<SubNet> subnet;
    vector<vector<Path_table_ver2>> path_table_ver2;
    map<pair<int, int>, int> channel_demand; //2 fpga --> demand signals
    map<pair<int, int>, Channel *> map_to_channel;
    map<pair<int, int>, int> channel_capacity; //2 fpga --> channel capacity

    FPGA_Gr()
    {
        round = 1;
        sink_num = total_cost = total_demand = 0;
        maxsgw = maxtdm = 0;
        minsgw = mintdm = INT_MAX;
    }
    ~FPGA_Gr() {}
    void getfile(char *, char *);
    void breakdown(); //break down all net into 2 pin subnet
    void construct_table();
    void show_path_table();
    void add_channel_demand(const int &, const int &);
    void sub_channel_demand(const int &, const int &);
    void global_routing();
    void routing_tree(Net &, const vector<vector<int>> &);
    void compute_edge_weight(Net &, Tree_Node *);
    bool find_rip_up_edge(Net &, int &, int &); //input net root and get edge with parent and child
    Tree_Node *rip_up_edge(Net &, const int &, const int &);
    void reroute_edge(Net &, const int &, Tree_Node *, double);
    void rip_up_reroute(time_t);
    int channel_used(int, int);
    double channel_TDM(int, int);
    double compute_TDM_cost();
    double comptue_tree_TDM_cost(Tree_Node *);
    void record_net_channel_used();
    void show_net_channel_table();
    void output_file(char *output, time_t);
    void check_result();

    //another global routing
    void construct_table_ver2(); //考慮hops數多1~2的可能
    void global_routing_ver2();  //考慮tdm(orcd  congestion)
    void show_path_table_ver2();
    double compute_cost_for_gr2(Net &, const vector<int> &, const SubNet &, int routed_net);

    //history cost 2020/03/29
    void initial_route_result();

    //channel direct 2020/04/08
    //void distribute_channel_capacity(); //依比例分配channel的capacity
};

#endif