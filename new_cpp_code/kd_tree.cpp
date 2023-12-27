// kd_tree.cpp

#include "kd_tree.h"
#include "para.h"
#include <unordered_set>
#include <iostream>
#include <queue>

kd_tree::Node * kd_tree::newNode(int i, double *r = nullptr, Node *left=nullptr, Node *right=nullptr) {
    Node * temp = new Node();
    temp->ind = i;
    for (size_t j = 0; j < DIMENSION; j++)
    {
        temp->vals[ j] = r[ i *DIMENSION +  j];
    }
    temp->left = left;
    temp->right = right;
    return temp;
}

void kd_tree::displayTree(kd_tree::Node * node) {
    if (node)
    {
        std::cout << "(" << node->ind << ", ";
        displayTree(node->left);
        std::cout << ", ";
        displayTree(node->right);
        std::cout << ")";
    }
    
}


kd_tree::kd_tree()
{
    root = nullptr;
    data = new Node*[NUMBER_OF_PARTICLES];
    for (size_t i = 0; i < NUMBER_OF_PARTICLES; i++)
    {
        data[i] = new Node();
    }
    
}


void kd_tree::build(double *r) {

    for (size_t i = 0; i < NUMBER_OF_PARTICLES; i++)
    {
        data[i]->ind = i;
        for (size_t j = 0; j < DIMENSION; j++)
        {
            data[i]->vals[j] = r[i*DIMENSION + j];
        }

        // set both left and right ptr to nullptr 
        data[i]->left = nullptr;
        data[i]->right = nullptr;
    }

    root = recBuild(0, NUMBER_OF_PARTICLES-1, 0);    

    if (!verifyBuild(root, 0))
        std::cout << "build failed" << std::endl;    


    int number_of_points = 0;
    recCount(root, number_of_points);

    if (number_of_points != NUMBER_OF_PARTICLES)
    {
        std::cout << "Error: " << number_of_points << " != " << NUMBER_OF_PARTICLES << std::endl;
    }         
}

bool kd_tree::verifyBuild(Node *node, int axis) {    
    if (node == nullptr)
        return true;

    std::queue<Node*> q;    
    if (node->left) {
        q.push(node->left);
        Node *t;
        while (q.size()) {
            t = q.front();
            if (t->vals[axis] >= node->vals[axis]) return false; 
            q.pop();
            if (t->left) q.push(t->left);
            if (t->right) q.push(t->right);
        }
    }

    while (q.size()) q.pop();
    
    if (node->right) {
        q.push(node->right);
        Node *t;
        while (q.size()) {
            t = q.front();
            if (t->vals[axis] < node->vals[axis]) return false;
            q.pop();
            if (t->left) q.push(t->left);
            if (t->right) q.push(t->right);
        }
    }
    return verifyBuild(node->left, (axis + 1) % DIMENSION) && verifyBuild(node->right, (axis+1) % DIMENSION);


}


void kd_tree::recCount(Node *node, int &number_of_points) {
    if (node == nullptr) return;
    number_of_points++;
    recCount(node->left, number_of_points);
    recCount(node->right, number_of_points);
}

kd_tree::Node * kd_tree::recBuild(int i, int j, int axis) {

    if (i > j) return nullptr;

    if (i == j) return data[i];
    
    // find the middle index
    int mid = (i + j) / 2;

    // sort the array of Node pointers based on the current axis
    // both i and j index are included
    std::sort(data + i, data + j + 1, [&axis](auto *a, auto *b) {
        return  a->vals[axis] < b->vals[axis];
    });

    // handle the case if it has duplicate values
    while (mid > 0 && mid > i && data[mid]->vals[axis] == data[mid-1]->vals[axis] ) mid--;
        

    // create a new Node with the middle value
    // recursively build the left and right subtrees
    data[mid]->left = recBuild( i, mid-1, (axis + 1) % DIMENSION);
    data[mid]->right = recBuild(mid + 1, j, (axis + 1) % DIMENSION);

    return data[mid];
}


void kd_tree::searchNeighbour(int ind, double *r, std::vector<int> &indices){
    Node * key_node = newNode(ind, r);
    std::unordered_set<int> visited = {ind};
    recSearchNeighbour(key_node , visited,  root, 0);
    auto it = visited.begin();

    while (it != visited.end())
    {
        if (*it != ind)
        {
            if (ind < *it) {
                indices.push_back(ind);
                indices.push_back(*it);
            } else {
                indices.push_back(*it);
                indices.push_back(ind);    
            }
        }        
        it++;
    }
    delete key_node;
}

void kd_tree::recSearchNeighbour(Node *key_node, std::unordered_set<int> &visited,  Node *node, int axis){
    if (node == nullptr) return;

    if ( visited.count(node->ind) > 0 ) return;

    
    double r2 = 0.;
    for (size_t j = 0; j < DIMENSION; j++)
        r2 += (key_node->vals[j] - node->vals[j]) * (key_node->vals[j] - node->vals[j]);

    if (r2 < CUTOFF_DISTANCE * CUTOFF_DISTANCE)
        visited.insert(node->ind);

    if (key_node->vals[axis] - CUTOFF_DISTANCE <= node->vals[axis])
    {
        recSearchNeighbour(key_node, visited, node->left, (axis + 1) % DIMENSION);
    }
    
    if (key_node->vals[axis] + CUTOFF_DISTANCE > node->vals[axis])
    {
        recSearchNeighbour(key_node, visited, node->right, (axis + 1) % DIMENSION);
    }


}



kd_tree::~kd_tree()
{
    // delete pos
    for (size_t i = 0; i < NUMBER_OF_PARTICLES; i++)
    {
        delete data[i];
    }
    
    delete[] data;
    delete root;
}

