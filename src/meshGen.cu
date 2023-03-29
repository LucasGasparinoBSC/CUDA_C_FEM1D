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