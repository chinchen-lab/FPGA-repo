#ifndef NODE_HPP
#define NODE_HPP

#include <iostream>
#include <queue>
#include <list>
using namespace std;

class Tree_Node
{
public:
    int fpga_id;
    int max_value; //max value
    int sink_weight;
    int edge_weight;
    bool flag;
    Tree_Node *parent;
    list<Tree_Node *> children;

    Tree_Node() 
    {
        sink_weight = 0;
        flag = false;
    }
    ~Tree_Node() {}
};

#endif