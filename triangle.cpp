#include "triangle.h"
#include "vertex.h"
#include "math.h"
#include <vector>

Triangle::Triangle()
{
    v[0] = -1;
    v[1] = -1;
    v[2] = -1;

    t[0] = -1;
    t[1] = -1;
    t[2] = -1;
}

Triangle::Triangle(int v0, int v1, int v2, int t0, int t1, int t2)
{
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;

    t[0] = t0;
    t[1] = t1;
    t[2] = t2;
}

bool is_clockwise_rotation(Vertex a, Vertex b, Vertex c)
{
    double ab [2];
    ab[0] = b.x - a.x;
    ab[1] = b.y - a.y;

    double ac [2];
    ac[0] = c.x - a.x;
    ac[1] = c.y - a.y;

    return (ab[0] * ac[1] - ab[1] * ac[0] > 0.0);
}

bool Triangle::contains(std::vector<Vertex>& vertices, double x, double y, double z)
{
    Vertex a = vertices[this->v[0]];
    Vertex b = vertices[this->v[1]];
    Vertex c = vertices[this->v[2]];

    Vertex vertex = Vertex(x, y, z, 0);

    return is_clockwise_rotation(a, b, vertex)
        && is_clockwise_rotation(b, c, vertex)
        && is_clockwise_rotation(c, a, vertex);
}

double Triangle::get_light(std::vector<Vertex>& vertices)
{
    Vertex a = vertices[this->v[0]];
    Vertex b = vertices[this->v[1]];
    Vertex c = vertices[this->v[2]];

    double u [3];
    u[0] = b.x - a.x;
    u[1] = b.y - a.y;
    u[2] = b.z - a.z;

    double v [3];
    v[0] = c.x - a.x;
    v[1] = c.y - a.y;
    v[2] = c.z - a.z;

    double normal [3];
    normal[0] = u[1] * v[2] - u[2] * v[1];
    normal[1] = u[2] * v[0] - u[0] * v[2];
    normal[2] = u[0] * v[1] - u[1] * v[0];

    // normalize
    double kn = 1 / sqrt(pow(normal[0],2) + pow(normal[1],2) + pow(normal[2],2));
    normal[0] = kn * normal[0];
    normal[1] = kn * normal[1];
    normal[2] = kn * normal[2];

    double light[] = {0.0, 0.0, 1.0};

    return normal[0] * light[0] + normal[1] * light[1] + normal[2] * light[2];
}
