#include "meshGen.cuh"

// Node class methods

// Constructors: Defualt, Copy, and Parameterized
Node::Node() : n_Id(-1), n_X(0.0f), n_Val(0.0f){}
Node::Node(Node &in) : n_Id(in.n_Id), n_X(in.n_X), n_Val(in.n_Val){}
Node::Node(int &n, float &f, float &v) : n_Id(n), n_X(f), n_Val(v){}

// Destructor
Node::~Node(){}

// Getters and Extractors
const int& Node::get_nId() {return n_Id;}
const float& Node::get_nX() {return n_X;}
const float& Node::get_nVal() {return n_Val;}
int& Node::extract_nId() {return n_Id;}
float& Node::extract_nX() {return n_X;}
float& Node::extract_nVal() {return n_Val;}

// Setters
void Node::set_nId(int &n)  {n_Id = n;}
void Node::set_nX(float &f) {n_X = f;}
void Node::set_nVal(float &v) {n_Val = v;}

// GaussPoint class methods

// Constructors: Defualt, Copy, and Parameterized
GaussPoint::GaussPoint() : Node(), gp_Weight(0.0f){}
GaussPoint::GaussPoint(GaussPoint &in) : Node(in), gp_Weight(in.gp_Weight){}
GaussPoint::GaussPoint(int &n, float &f, float &v, float &w) : Node(n, f, v), gp_Weight(w){}

// Destructor
GaussPoint::~GaussPoint(){}

// Getters and Extractors
const float& GaussPoint::get_gpWeight() {return gp_Weight;}
float& GaussPoint::extract_gpWeight() {return gp_Weight;}

// Setters
void GaussPoint::set_gpWeight(float &w) {gp_Weight = w;}

// MasterElement class methods

// Constructors: Defualt, Copy, and Parameterized
MasterElement::MasterElement() : e_Order(-1), e_NumNodes(-1), e_NumGPs(-1), e_hoNodes(-1){}
MasterElement::MasterElement(MasterElement &in) : e_Order(in.e_Order), e_NumNodes(in.e_NumNodes), e_NumGPs(in.e_NumGPs), e_hoNodes(in.e_hoNodes), e_NodeOrdering(in.e_NodeOrdering), e_Ngp(in.e_Ngp), e_dNgp(in.e_dNgp), e_GaussPoints(in.e_GaussPoints){}
MasterElement::MasterElement(int &p) : e_Order(p), e_NumNodes(p+1), e_NumGPs(p+1), e_hoNodes(e_NumNodes-2)
{
    // If order p greater than MAX_ORDER, abort
    if (e_Order > MAX_ORDER)
    {
        printf("ERROR: Order of element greater than MAX_ORDER. Aborting.\n");
        exit(1);
    }
    // Allocate arrays for node ordering, shape functions, shape function derivatives and Gauss points
    e_NodeOrdering = (int*)malloc(e_NumNodes*sizeof(int));
    e_Ngp = (float*)malloc(e_NumNodes*e_NumGPs*sizeof(float));
    e_dNgp = (float*)malloc(e_NumNodes*e_NumGPs*sizeof(float));
    e_GaussPoints = (GaussPoint*)malloc(e_NumGPs*sizeof(GaussPoint));

    // Create Gauss points, shape functions and derivatives
    createGaussPoints();
    createNodeOrdering();
    createShapeFunctions();
}

// Destructor
MasterElement::~MasterElement()
{
    // Free arrays for node ordering, shape functions, shape function derivatives and Gauss points
    free(e_NodeOrdering);
    free(e_Ngp);
    free(e_dNgp);
    free(e_GaussPoints);
}

// Define Gauss points
void MasterElement::createGaussPoints()
{
    // Create a table of possible Gauss points based on order p
    float auxXgp[MAX_NGAUS];
    float auxWgp[MAX_NGAUS];
    memset(auxXgp, 0.0f, MAX_NGAUS*sizeof(float));
    memset(auxWgp, 0.0f, MAX_NGAUS*sizeof(float));
    if (e_Order == 0) {
        auxXgp[0] = 0.0f;
        auxWgp[0] = 2.0f;
    }
    else if (e_Order == 1) {
        auxXgp[0] = -0.577350269189626f;
        auxXgp[1] = 0.577350269189626f;
        auxWgp[0] = 1.0f;
        auxWgp[1] = 1.0f;
    }
    else if (e_Order == 2) {
        auxXgp[0] = -0.774596669241483f;
        auxXgp[1] = 0.0f;
        auxXgp[2] = 0.774596669241483f;
        auxWgp[0] = 0.555555555555556f;
        auxWgp[1] = 0.888888888888889f;
        auxWgp[2] = 0.555555555555556f;
    }
    else if (e_Order == 3) {
        auxXgp[0] = -0.861136311594053f;
        auxXgp[1] = -0.339981043584856f;
        auxXgp[2] = 0.339981043584856f;
        auxXgp[3] = 0.861136311594053f;
        auxWgp[0] = 0.347854845137454f;
        auxWgp[1] = 0.652145154862546f;
        auxWgp[2] = 0.652145154862546f;
        auxWgp[3] = 0.347854845137454f;
    }
    else if (e_Order == 4) {
        auxXgp[0] = -0.906179845938664f;
        auxXgp[1] = -0.538469310105683f;
        auxXgp[2] = 0.0f;
        auxXgp[3] = 0.538469310105683f;
        auxXgp[4] = 0.906179845938664f;
        auxWgp[0] = 0.236926885056189f;
        auxWgp[1] = 0.478628670499366f;
        auxWgp[2] = 0.568888888888889f;
        auxWgp[3] = 0.478628670499366f;
        auxWgp[4] = 0.236926885056189f;
    }
    else if (e_Order == 5) {
        auxXgp[0] = -0.932469514203152f;
        auxXgp[1] = -0.661209386466265f;
        auxXgp[2] = -0.238619186083197f;
        auxXgp[3] = 0.238619186083197f;
        auxXgp[4] = 0.661209386466265f;
        auxXgp[5] = 0.932469514203152f;
        auxWgp[0] = 0.171324492379170f;
        auxWgp[1] = 0.360761573048139f;
        auxWgp[2] = 0.467913934572691f;
        auxWgp[3] = 0.467913934572691f;
        auxWgp[4] = 0.360761573048139f;
        auxWgp[5] = 0.171324492379170f;
    }
    else if (e_Order == 6) {
        auxXgp[0] = -0.949107912342759f;
        auxXgp[1] = -0.741531185599394f;
        auxXgp[2] = -0.405845151377397f;
        auxXgp[3] = 0.0f;
        auxXgp[4] = 0.405845151377397f;
        auxXgp[5] = 0.741531185599394f;
        auxXgp[6] = 0.949107912342759f;
        auxWgp[0] = 0.129484966168870f;
        auxWgp[1] = 0.279705391489277f;
        auxWgp[2] = 0.381830050505119f;
        auxWgp[3] = 0.417959183673469f;
        auxWgp[4] = 0.381830050505119f;
        auxWgp[5] = 0.279705391489277f;
        auxWgp[6] = 0.129484966168870f;
    }
    else if (e_Order == 7) {
        auxXgp[0] = -0.960289856497536f;
        auxXgp[1] = -0.796666477413627f;
        auxXgp[2] = -0.525532409916329f;
        auxXgp[3] = -0.183434642495650f;
        auxXgp[4] = 0.183434642495650f;
        auxXgp[5] = 0.525532409916329f;
        auxXgp[6] = 0.796666477413627f;
        auxXgp[7] = 0.960289856497536f;
        auxWgp[0] = 0.101228536290376f;
        auxWgp[1] = 0.222381034453374f;
        auxWgp[2] = 0.313706645877887f;
        auxWgp[3] = 0.362683783378362f;
        auxWgp[4] = 0.362683783378362f;
        auxWgp[5] = 0.313706645877887f;
        auxWgp[6] = 0.222381034453374f;
        auxWgp[7] = 0.101228536290376f;
    }
    else if (e_Order == 8) {
        printf("Gauss quadrature order 8 not implemented yet!\n");
        exit(1);
    }

    // Set Gauss points
    for (int i = 0; i < e_NumGPs; i++)
    {
        // Set Gauss point ID
        e_GaussPoints[i].set_nId(i);
        // Set Gauss point location
        e_GaussPoints[i].set_nX(auxXgp[i]);
        // Set Gauss point weight
        e_GaussPoints[i].set_gpWeight(auxWgp[i]);
    }

    // Free memory
    free(auxXgp);
    free(auxWgp);
}

// Define element orrdering
void MasterElement::createNodeOrdering()
{
    e_NodeOrdering[0] = 0;
    e_NodeOrdering[1] = 1;
    if (e_Order > 1) {
        for (int i = 2; i < e_NumNodes; i++)
        {
            e_NodeOrdering[i] = i;
        }
    }
}

// Define element shape functions and derivatives using Lagrangian polynomials
void MasterElement::createShapeFunctions()
{
    // Form the element grid (equispaced)
    float *xgrid = (float *)malloc(e_NumNodes * sizeof(float));
    xgrid[0] = -1.0f;
    xgrid[1] = 1.0f;
    for (int inode = 2; inode < e_NumNodes; inode++)
    {
        xgrid[inode] = -1.0f + 2.0f * ((float)(inode - 1) / (float)(e_NumNodes - 1));
    }

    // Compute shape functions and derivatives at every gauss point
    for (int igaus = 0; igaus < e_NumGPs; igaus++)
    {
        // Get Gauss point location
        const float xgp = e_GaussPoints[igaus].get_nX();

        // Compute Lagrangian polynomial of order p
        for (int inode = 0; inode < e_NumNodes; inode++)
        {
            for (int j = 0; j < e_NumNodes; j++)
            {
                float N = 1.0f;
                if (j != inode)
                {
                    N *= (xgp - xgrid[j]) / (xgrid[inode] - xgrid[j]);
                }
                e_Ngp[igaus* e_NumNodes + inode] = N;
            }
        }

        // Compute derivatives of Lagrangian polynomial of order p
        for (int inode = 0; inode < e_NumNodes; inode++)
        {
            for (int j = 0; j < e_NumNodes; j++)
            {
                float dN = 0.0f;
                if (j != inode)
                {
                    float l = 1.0f;
                    for (int k = 0; k < e_NumNodes; k++)
                    {
                        if (k != inode && k != j)
                        {
                            l *= (xgp - xgrid[k]) / (xgrid[j] - xgrid[k]);
                        }
                    }
                    dN += l/(xgrid[inode] - xgrid[j]);
                }
                e_dNgp[igaus* e_NumNodes + inode] = dN;
            }
        }
    }

    // Freee the memory
    free(xgrid);
}

// Getters and extractors
const int& MasterElement::get_eOrder() {return e_Order;}
const int& MasterElement::get_eNumNodes() {return e_NumNodes;}
const int& MasterElement::get_eNumGPs() {return e_NumGPs;}
const int& MasterElement::get_ehoNodes() {return e_hoNodes;}
const int* MasterElement::get_eNodeOrdering() {return e_NodeOrdering;}
const float* MasterElement::get_eNgp() {return e_Ngp;}
const float* MasterElement::get_edNgp() {return e_dNgp;}
const GaussPoint* MasterElement::get_eGaussPoints() {return e_GaussPoints;}
int& MasterElement::extract_eOrder() {return e_Order;}
int& MasterElement::extract_eNumNodes() {return e_NumNodes;}
int& MasterElement::extract_eNumGPs() {return e_NumGPs;}
int& MasterElement::extract_ehoNodes() {return e_hoNodes;}
int* MasterElement::extract_eNodeOrdering() {return e_NodeOrdering;}
float* MasterElement::extract_eNgp() {return e_Ngp;}
float* MasterElement::extract_edNgp() {return e_dNgp;}
GaussPoint* MasterElement::extract_eGaussPoints() {return e_GaussPoints;}

// Element class methods

// Constructors: default, copy and parameterized
Element::Element() : MasterElement(), e_Id(-1), e_Length(0.0f){}
Element::Element(Element &in) : MasterElement(in), e_Id(in.e_Id), e_Length(in.e_Length), e_Nodes(in.e_Nodes){}
Element::Element(int &p, int &id) : MasterElement(p), e_Id(id)
{
    // Create array of element nodes
    e_Nodes = (Node *)malloc(e_NumNodes * sizeof(Node));
}