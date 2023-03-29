#include "meshGen.cuh"

// Node class methods

// Constructor
Node::Node(int n, float f, float v)
{
    n_Id = n;
    n_X = f;
    n_Val = v;
}

// Destructor
Node::~Node()
{
}

// Getters
int Node::get_nId()
{
    return n_Id;
}

float Node::get_nX()
{
    return n_X;
}

float Node::get_nVal()
{
    return n_Val;
}

// Setters
void Node::set_nId(int n)
{
    n_Id = n;
}

void Node::set_nX(float f)
{
    n_X = f;
}

void Node::set_nVal(float v)
{
    n_Val = v;
}

// GaussPoint class methods

// Getters
float GaussPoint::get_gpWeight()
{
    return gp_Weight;
}

// Setters
void GaussPoint::set_gpWeight(float w)
{
    gp_Weight = w;
}

// Element class methods

// Constructor
Element::Element(int id, int p)
{
    e_Id = id;
    e_Order = p;
    e_NumNodes = p + 1;
    e_NumGPs = p + 1;
    e_Nodes = (Node*)malloc(e_NumNodes * sizeof(Node));
    e_GPs = (GaussPoint*)malloc(e_NumGPs * sizeof(GaussPoint));
}

// Destructor
Element::~Element()
{
    free(e_Nodes);
    free(e_GPs);
}