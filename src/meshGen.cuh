#ifndef MESHGEN_H
#define MESHGEN_H

#include <stdlib.h>
#include "defConstants.cuh"

// Create a node class
class Node
{
    protected:
        int n_Id;
        float n_X;
        float n_Val;
    public:
        Node(int n, float f, float v);
        ~Node();
        int get_nId();
        float get_nX();
        float get_nVal();
        void set_nId(int n);
        void set_nX(float f);
        void set_nVal(float f);
};

// Create a Gauss point class that inherits from the Node class
class GaussPoint : public Node
{
    private:
        float gp_Weight;
    public:
        float get_gpWeight();
        void set_gpWeight(float w);
};

// Create an Element class
class Element
{
    private:
        int e_Id;
        int e_Order;
        int e_NumNodes;
        int e_NumGPs;
        Node *e_Nodes;
        GaussPoint *e_GPs;
    public:
        Element(int id, int p);
        ~Element();
        int get_eId();
        int get_eOrder();
        int get_eNumNodes();
        int get_eNumGP();
        Node *get_eNodes();
        GaussPoint *get_eGPs();
        void set_eId(int id);
        void set_eOrder(int p);
};

// Create a Mesh class
class Mesh
{
    private:
        int m_NumElements;
        int m_NumNodes;
    public:
};

#endif // MESHGEN_H