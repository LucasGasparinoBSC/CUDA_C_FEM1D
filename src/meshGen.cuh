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
        GaussPoint *e_GaussPoints;
    private:
        inline void createGaussPoints();
        inline void createNodeOrdering();
        inline void createShapeFunctions();
    public:
        inline MasterElement();
        inline MasterElement(MasterElement &in);
        inline MasterElement(int &p);
        inline ~MasterElement();
        inline const int& get_eOrder();
        inline const int& get_eNumNodes();
        inline const int& get_eNumGPs();
        inline const int& get_ehoNodes();
        inline const int* get_eNodeOrdering();
        inline const float* get_eNgp();
        inline const float* get_edNgp();
        inline const GaussPoint* get_eGaussPoints();
        inline int& extract_eOrder();
        inline int& extract_eNumNodes();
        inline int& extract_eNumGPs();
        inline int& extract_ehoNodes();
        inline int* extract_eNodeOrdering();
        inline float* extract_eNgp();
        inline float* extract_edNgp();
        inline GaussPoint* extract_eGaussPoints();
};

// Create an Element class that extends the MasterElement class
class Element : public MasterElement
{
    private:
        int e_Id;
        float e_Length;
        Node *e_Nodes;
    public:
        inline Element();
        inline Element(Element &in);
        inline Element(int &p, int &id);
        inline ~Element();
};

#endif // MESHGEN_H