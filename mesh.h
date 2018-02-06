#ifndef MESH_H
#define MESH_H

#include "vertex.h"
#include "triangle.h"
#include "colors.h"
#include <vector>

using namespace std;

struct point {
    double x, y, z;
};

struct IndexPair {
    int a, b;
};

class Mesh
{
public:
    vector<Vertex> vertices;
    vector<Triangle> triangles;
    Colors colors;

    Mesh();
    void build_mesh();
    void add_random_vertex();
    void add_vertex(double x, double y, double z);
    int get_target_triangle(double x, double y, double z);
    void add_vertex_inside(int target, double x, double y, double z);
    bool is_inside(Triangle t, double x, double y, double z);
    void flip(int a_index, int b_index);
    int is_delaunay(int t_index);
    int get_opposite_vertex(int t_index, int v_index);
    void flip_algorithm();
    void lawson_algorithm();
    bool should_flip(int a_index, int b_index);
    void log_triangles();
    void draw();
};

#endif // MESH_H
