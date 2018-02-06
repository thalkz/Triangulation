#ifndef TRIANGLE_H
#define TRIANGLE_H

#include"vertex.h"
#include <vector>

class Triangle
{
public:
    int t[3];
    int v[3];
    bool is_box = false;

    Triangle();
    Triangle(int v0, int v1, int v2, int t0, int t1, int t2);

    bool contains(std::vector<Vertex>& vertices, double x, double y, double z);
    double get_light(std::vector<Vertex>& vertices);
};

#endif // TRIANGLE_H
