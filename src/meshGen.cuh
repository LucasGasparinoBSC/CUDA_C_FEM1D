#ifndef MESHGEN_H
#define MESHGEN_H

#include <stdlib.h>
#include <stdio.h>
#include "defConstants.cuh"

// Create a node class
class Node
{
    protected:
        int n_Id;
        float n_X;
        float n_Val;
    public:
        Node();
        Node(Node &in);
        Node(int &n, float &f, float &v);
        ~Node();
        inline const int& get_nId();
        inline const float& get_nX();
        inline const float& get_nVal();
        inline int& extract_nId();
        inline float& extract_nX();
        inline float& extract_nVal();
        inline void set_nId(int &n);
        inline void set_nX(float &f);
        inline void set_nVal(float &v);
};

// Create a Gauss point class that inherits from the Node class
class GaussPoint : public Node
{
    private:
        float gp_Weight;
    public:
        inline GaussPoint();
        inline GaussPoint(GaussPoint &in);
        inline GaussPoint(int &n, float &f, float &v, float &w);
        inline ~GaussPoint();
        inline const float& get_gpWeight();
        inline float& extract_gpWeight();
        inline void set_gpWeight(float &w);
};

// Create a Master element class
class MasterElement
{
    protected:
        int e_Order;
        int e_NumNodes;
        int e_NumGPs;
        int e_hoNodes;
        int *e_NodeOrdering;
        float *e_Ngp;
        float *e_dNgp;
    public:
        MasterElement(int p){};
        ~MasterElement();
        int get_eOrder();
        int get_eNumNodes();
        int get_eNumGPs();
        int get_ehoNodes();
        int *get_eMasterNodes();
};

#endif // MESHGEN_H