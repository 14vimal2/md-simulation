// kd_tree.h

#include "para.h"
#include <algorithm> 
#include <unordered_set>
#include <vector>

#ifndef _KD_TREE_
#define _KD_TREE_



class kd_tree
{
private:
    struct Node
    {
        int ind;
        double vals[DIMENSION];
        Node *left, *right;
    };
    Node *newNode(int ind, double *val, Node *l, Node *r);    

    Node **data;
    Node *root;
public:
    kd_tree();
    ~kd_tree();
    void build(double *r);
    bool verifyBuild(Node *node, int axis);
    Node * recBuild(int i, int j, int axis);
    void searchNeighbour(int ind, double *r, std::vector<int> &indices);
    void recSearchNeighbour(Node *key_node, std::unordered_set<int> &visited,  Node *node, int axis);
    void recCount(Node *node, int &number_of_points);
    void displayTree(Node *node);
};

#endif

