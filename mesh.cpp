#include "mesh.h"
#include <gl.h>
#include <cstdlib>
#include <time.h>
#include <math.h>
#include <iostream>
#include "colors.h"

Mesh::Mesh()
{
    colors = Colors();

    build_mesh();
    log_triangles();
}

void Mesh::build_mesh()
{
    vertices.clear();
    triangles.clear();

    double mesh_size = 0.5;

    vertices.push_back(Vertex(-1.0 * mesh_size, -1.0 * mesh_size, 0.0, 0)); // s0
    vertices.push_back(Vertex(1.0 * mesh_size, -1.0 * mesh_size, 0.0, 0)); // s1
    vertices.push_back(Vertex(1.0 * mesh_size, 1.0 * mesh_size, 0.0, 0)); // s2
    vertices.push_back(Vertex(-1.0 * mesh_size, 1.0 * mesh_size, 0.0, 1)); // s3
    vertices.push_back(Vertex(0.0, 0.0, -3.0 * mesh_size, 2)); // s4

    triangles.push_back(Triangle(0, 1, 2, 3, 1, 4)); // t0
    triangles.push_back(Triangle(0, 2, 3, 2, 5, 0)); // t1
    triangles.push_back(Triangle(3, 2, 4, 3, 5, 1)); // t2
    triangles.push_back(Triangle(2, 1, 4, 4, 2, 0)); // t3
    triangles.push_back(Triangle(0, 4, 1, 3, 0, 5)); // t4
    triangles.push_back(Triangle(0, 3, 4, 2, 4, 1)); // t5

    triangles[2].is_box = true;
    triangles[3].is_box = true;
    triangles[4].is_box = true;
    triangles[5].is_box = true;
}

void Mesh::log_triangles()
{
    cout << "Log Triangles : " << endl;
    for (int i = 0; i < (int) triangles.size(); i++)
    {
        std::cout << "Triangle " << i << " (" << triangles[i].t[0] << ", "
                                              << triangles[i].t[1] << ", "
                                              << triangles[i].t[2] << ")" << endl;
        for (int j = 0; j < 3; j++)
        {
            Triangle other = triangles[triangles[i].t[j]];
            cout << "=> Ajdacent Triangle " << triangles[i].t[j] << "... ";
            if (other.t[0] == i || other.t[1] == i || other.t[2] == i)
            {
                cout << "OK !" << endl;
            }
            else
            {
                cout << "Error !" << endl;
            }
        }
    }

    cout << endl;
}

void Mesh::add_random_vertex()
{
    double x = (rand() % 1000) / 1000.0 - 0.5;
    double y = (rand() % 1000) / 1000.0 - 0.5;
    //double z = (rand() % 1000) / 10000.0 - 0.05;
    double z = 0;

    int target = get_target_triangle(x, y, z);

    if (target != -1)
    {
        Triangle old = triangles[target];

        int new_vertex = vertices.size();

        int t_index [3];
        t_index[0] = target;
        t_index[1] = triangles.size();
        t_index[2] = triangles.size() + 1;

        Triangle t0 = Triangle(
                    old.v[0],
                    old.v[1],
                    new_vertex,
                    t_index[1],
                    t_index[2],
                    old.t[2]);

        Triangle t1 = Triangle(
                    new_vertex,
                    old.v[1],
                    old.v[2],
                    old.t[0],
                    t_index[2],
                    t_index[0]);

        Triangle t2 = Triangle(
                    old.v[0],
                    new_vertex,
                    old.v[2],
                    t_index[1],
                    old.t[1],
                    t_index[0]);

        vertices.push_back(Vertex(x, y, z, target));

        triangles[target] = t0;
        triangles.push_back(t1);
        triangles.push_back(t2);

        // Raccorder le triangle en face de v0
        for (int k = 0; k < 3; k++)
        {
            if (triangles[old.t[0]].t[k] == target)
            {
                triangles[old.t[0]].t[k] = t_index[1];
            }
        }

        // Raccorder le triangle en face de v1
        for (int k = 0; k < 3; k++)
        {
            if (triangles[old.t[1]].t[k] == target)
            {
                triangles[old.t[1]].t[k] = t_index[2];
            }
        }

        // Raccorder le triangle en face de v2
        for (int k = 0; k < 3; k++)
        {
            if (triangles[old.t[2]].t[k] == target)
            {
                triangles[old.t[2]].t[k] = t_index[0];
            }
        }
    }
    else
    {
        std::cout << "Cannot add Vertex at : " << x << ", " << y << ", " << z << endl;
    }

    log_triangles();
}

void Mesh::log_queue(vector<IndexPair> &to_flip)
{
    cout << endl << "Initial Queue contains : " << endl;

    for (int i = 0; i < (int) to_flip.size(); i++)
    {
        cout << to_flip.at(i).a << ", " << to_flip.at(i).b << endl;
    }
    cout << endl;
}

void Mesh::add_delaunay_vertex()
{
    double x = (rand() % 1000) / 1000.0 - 0.5;
    double y = (rand() % 1000) / 1000.0 - 0.5;
    //double z = (rand() % 1000) / 10000.0 - 0.05;
    double z = 0;

    int target = get_target_triangle(x, y, z);

    if (target != -1)
    {
        Triangle old = triangles[target];

        int new_vertex = vertices.size();

        int t_index [3];
        t_index[0] = target;
        t_index[1] = triangles.size();
        t_index[2] = triangles.size() + 1;

        Triangle t0 = Triangle(
                    old.v[0],
                    old.v[1],
                    new_vertex,
                    t_index[1],
                    t_index[2],
                    old.t[2]);

        Triangle t1 = Triangle(
                    new_vertex,
                    old.v[1],
                    old.v[2],
                    old.t[0],
                    t_index[2],
                    t_index[0]);

        Triangle t2 = Triangle(
                    old.v[0],
                    new_vertex,
                    old.v[2],
                    t_index[1],
                    old.t[1],
                    t_index[0]);

        vertices.push_back(Vertex(x, y, z, target));

        triangles[target] = t0;
        triangles.push_back(t1);
        triangles.push_back(t2);

        // Raccorder le triangle en face de v0
        for (int k = 0; k < 3; k++)
        {
            if (triangles[old.t[0]].t[k] == target)
            {
                triangles[old.t[0]].t[k] = t_index[1];
            }
        }

        // Raccorder le triangle en face de v1
        for (int k = 0; k < 3; k++)
        {
            if (triangles[old.t[1]].t[k] == target)
            {
                triangles[old.t[1]].t[k] = t_index[2];
            }
        }

        // Raccorder le triangle en face de v2
        for (int k = 0; k < 3; k++)
        {
            if (triangles[old.t[2]].t[k] == target)
            {
                triangles[old.t[2]].t[k] = t_index[0];
            }
        }

        // Ajouter les 3 arrêtes et faire l'algo de flip
        vector<IndexPair> to_flip;

        if (!triangles[old.t[2]].is_box) {
            to_flip.push_back({t_index[0], old.t[2]});
        }

        if (!triangles[old.t[0]].is_box) {
            to_flip.push_back({t_index[1], old.t[0]});
        }

        if (!triangles[old.t[1]].is_box) {
            to_flip.push_back({t_index[2], old.t[1]});
        }

        log_queue(to_flip);
        perform_flips(to_flip);
    }
    else
    {
        std::cout << "Cannot add Vertex at : " << x << ", " << y << ", " << z << endl;
    }


    log_triangles();
}

int Mesh::get_target_triangle(double x, double y, double z)
{
    for (int i = 0; i < (int) triangles.size(); i++)
    {
        if (triangles[i].contains(vertices, x, y, z))
        {
            return i;
        }
    }

    return -1;
}

void Mesh::flip(int a_index, int b_index)
{
    Triangle a = triangles[a_index];
    Triangle b = triangles[b_index];

    int opposite_a = -1; // Index of the vector opposite to a
    int opposite_b = -1; // Index of the vector opposite to b

    for (int i = 0; i < 3; i++)
    {
        if (a.v[i] != b.v[0] && a.v[i] != b.v[1] && a.v[i] != b.v[2])
        {
            opposite_a = i;
        }
    }

    for (int i = 0; i < 3; i++)
    {
        if (b.v[i] != a.v[0] && b.v[i] != a.v[1] && b.v[i] != a.v[2])
        {
            opposite_b = i;
        }
    }

    // Update vertices in a and b
    a.v[(opposite_a + 1) % 3] = b.v[opposite_b];
    b.v[(opposite_b + 1) % 3] = a.v[opposite_a];

    // Update triangles in
    a.t[opposite_a] = triangles[b_index].t[(opposite_b + 2) % 3];
    a.t[(opposite_a + 2) % 3] = b_index;

    // Update triangles in b
    b.t[opposite_b] = triangles[a_index].t[(opposite_a + 2) % 3];
    b.t[(opposite_b + 2) % 3] = a_index;

    // Update outside triangles
    for (int i = 0; i < 3; i++)
    {
        if (triangles[a.t[opposite_a]].t[i] == b_index)
        {
            triangles[a.t[opposite_a]].t[i] = a_index;
        }
    }

    for (int i = 0; i < 3; i++)
    {
        if (triangles[b.t[opposite_b]].t[i] == a_index)
        {
            triangles[b.t[opposite_b]].t[i] = b_index;
        }
    }

    // Save triangles
    triangles[a_index] = a;
    triangles[b_index] = b;
}

double get_angle(Vertex a, Vertex b, Vertex c)
{
    double ab[2];
    ab[0] = b.x - a.x;
    ab[1] = b.y - a.y;

    // normalize
    double kab = 1 / sqrt(pow(ab[0],2) + pow(ab[1],2));
    ab[0] = kab * ab[0];
    ab[1] = kab * ab[1];

    double ac[2];
    ac[0] = c.x - a.x;
    ac[1] = c.y - a.y;

    // normalize
    double kac = 1 / sqrt(pow(ac[0],2) + pow(ac[1],2));
    ac[0] = kac * ac[0];
    ac[1] = kac * ac[1];

    return acos(ab[0] * ac[0] + ab[1] * ac[1]);
}

Vertex get_center(Vertex a, Vertex b, Vertex c)
{
    double angle_a = get_angle(a, b, c);
    double angle_b = get_angle(b, c, a);
    double angle_c = get_angle(c, a, b);

//    cout << "angles : " << 180.0 * angle_a / 3.14 << ", " << 180.0 * angle_b / 3.14 << ", " << 180.0 * angle_c / 3.14 << endl;
//    cout << "somme des angles : " << 180.0 * (angle_a + angle_b + angle_c) / 3.14 << endl;

    double i = tan(angle_b) + tan(angle_c);
    double j = tan(angle_a) + tan(angle_c);
    double k = tan(angle_a) + tan(angle_b);

    double x = (a.x * i + b.x * j + c.x * k) / (i + j + k);
    double y = (a.y * i + b.y * j + c.y * k) / (i + j + k);

    return Vertex(x, y, 0.0);
}

double get_2d_distance(Vertex va, Vertex center)
{
    return sqrt(pow(va.x - center.x, 2) + pow(va.y - center.y, 2));
}

int Mesh::get_opposite_vertex(int t_index, int v_index)
{
    int opposite_t_index = triangles[t_index].t[v_index];

    //cout << "Vertex #" << v_index << " is opposite to Triangle #" << opposite_t_index << endl;

    Triangle opposite_t = triangles[opposite_t_index];

    for (int i = 0; i < 3; i++)
    {
        if (opposite_t.v[i] != triangles[t_index].v[0]
                && opposite_t.v[i] != triangles[t_index].v[1]
                && opposite_t.v[i] != triangles[t_index].v[2])
        {
            return opposite_t.v[i];
        }
    }
    cout << "Error : Could not find opposite Vertex : t_index=" << t_index << ", v_index=" << v_index << endl;
    return -1;
}

bool is_inside_circle(Vertex center, double radius, Vertex vertex)
{
    double dist = sqrt(pow(center.x - vertex.x, 2) + pow(center.y - vertex.y, 2));
    return dist < radius;
}

int Mesh::is_delaunay(int t_index)
{
    Triangle triangle = triangles[t_index];

    Vertex va = vertices[triangle.v[0]];
    Vertex vb = vertices[triangle.v[1]];
    Vertex vc = vertices[triangle.v[2]];

    Vertex center = get_center(va, vb, vc);
    double radius = get_2d_distance(va, center);

    // cout << "center (" << center.x << ", " << center.y << ", " << center.z << ")" << endl;
    // cout << "radius " << radius << endl;

    int opposite_va = get_opposite_vertex(t_index, 0);
    int opposite_vb = get_opposite_vertex(t_index, 1);
    int opposite_vc = get_opposite_vertex(t_index, 2);

    if (is_inside_circle(center, radius, vertices[opposite_va]) && opposite_va != 4)
    {
        return triangle.t[0];
    }
    else if (is_inside_circle(center, radius, vertices[opposite_vb]) && opposite_vb != 4)
    {
        return triangle.t[1];
    }
    else if (is_inside_circle(center, radius, vertices[opposite_vc]) && opposite_vc != 4)
    {
        return triangle.t[2];
    }
    else
    {
        return -1;
    }
}

bool is_above_3d_plane(Vertex &va, Vertex &vb, Vertex &vc, Vertex &opposite_v)
{
    double ab[3];
    ab[0] = vb.x - va.x;
    ab[1] = vb.y - va.y;
    ab[2] = vb.z - va.z;

    // normalize
    double kab = 1 / sqrt(pow(ab[0],2) + pow(ab[1],2) + pow(ab[2],2) );
    ab[0] = kab * ab[0];
    ab[1] = kab * ab[1];
    ab[2] = kab * ab[2];

    double ac[3];
    ac[0] = vc.x - va.x;
    ac[1] = vc.y - va.y;
    ac[2] = vc.z - va.z;

    // normalize
    double kac = 1 / sqrt(pow(ac[0],2) + pow(ac[1],2) + pow(ac[2],2) );
    ac[0] = kac * ac[0];
    ac[1] = kac * ac[1];
    ac[2] = kac * ac[2];

    double normal[3];
    normal[0] = ab[1] * ac[2] - ab[2] * ac[1];
    normal[1] = ab[2] * ac[0] - ab[0] * ac[2];
    normal[2] = ab[0] * ac[1] - ab[1] * ac[0];

    //cout << "Normal : " << normal[0] << ", " << normal[1] << ", " << normal[2] << endl;

    double opp[3];
    opp[0] = opposite_v.x - va.x;
    opp[1] = opposite_v.y - va.y;
    opp[2] = opposite_v.z - va.z;

    return (normal[0] * opp[0] + normal[1] * opp[1] + normal[2] * opp[2] < 0);
}

int Mesh::is_delaunay_predicat(int t_index)
{
    Triangle triangle = triangles[t_index];

    Vertex va = vertices[triangle.v[0]];
    Vertex vb = vertices[triangle.v[1]];
    Vertex vc = vertices[triangle.v[2]];

    // Projeter sur la parabole
    va.z = va.x * va.x + va.y * va.y;
    vb.z = vb.x * vb.x + vb.y * vb.y;
    vc.z = vc.x * vc.x + vc.y * vc.y;

    int opposite_va = get_opposite_vertex(t_index, 0);
    int opposite_vb = get_opposite_vertex(t_index, 1);
    int opposite_vc = get_opposite_vertex(t_index, 2);

    // Projeter sur la parabole
    Vertex opp_va = vertices[opposite_va];
    opp_va.z = opp_va.x * opp_va.x + opp_va.y * opp_va.y;

    Vertex opp_vb = vertices[opposite_vb];
    opp_vb.z = opp_vb.x * opp_vb.x + opp_vb.y * opp_vb.y;

    Vertex opp_vc = vertices[opposite_vc];
    opp_vc.z = opp_vc.x * opp_vc.x + opp_vc.y * opp_vc.y;

    if (is_above_3d_plane(va, vb, vc, opp_va) && opposite_va != 4)
    {
        return triangle.t[0];
    }
    else if (is_above_3d_plane(va, vb, vc, opp_vb) && opposite_vb != 4)
    {
        return triangle.t[1];
    }
    else if (is_above_3d_plane(va, vb, vc, opp_vc) && opposite_vc != 4)
    {
        return triangle.t[2];
    }
    else
    {
        return -1;
    }
}

bool Mesh::should_flip(int a_index, int b_index)
{
    cout << "Circle Method : " << is_delaunay(b_index) << endl;
    cout << "Predicat Method : " << is_delaunay_predicat(b_index) << endl;

    //return (a_index == is_delaunay(b_index));
    return (a_index == is_delaunay_predicat(b_index));
}

void Mesh::lawson_algorithm()
{
    vector<IndexPair> to_flip;

    // Find every triangles pair to flip and put them in 'to_flip'

    for (int t_index = 0; t_index < (int) triangles.size(); t_index++)
    {
        if (triangles[t_index].is_box)
        {
            continue;
        }

        cout << endl << "Is Triangle " << t_index << " Delaunay ?" << endl;
        int other_index = is_delaunay(t_index);

        if (other_index == -1)
        {
            cout << "Yes" << endl;
        }
        else
        {
            if (t_index < other_index) // This removes duplicate pairs
            {
                cout << "No, add Triangles (" << t_index << ", " << other_index << ") to queue." << endl;
                to_flip.push_back({t_index, other_index});
            }
            else
            {
                cout << "No, but indices aren't in ascending order." << endl;
            }
        }
    }
    //to_flip = fetch_pairs_to_flip();
    log_queue(to_flip);
    perform_flips(to_flip);
}

void Mesh::perform_flips(vector<IndexPair> &to_flip)
{
    int flip_count = 0;

    while(!to_flip.empty())
    {
        IndexPair pair = to_flip.back();
        to_flip.pop_back();
        cout << "Pop : " << pair.a << ", " << pair.b << endl;

        if (should_flip(pair.a, pair.b))
        {
            cout << "=> Flip" << endl;
            flip(pair.a, pair.b);
            flip_count++;

            // Ajouter les 2 arrêtes autour de a
            for (int k = 0; k < 3; k++)
            {
                if (triangles[pair.a].t[k] != pair.b)
                {
                    cout << "Push : " << triangles[pair.a].t[k] << ", " << pair.a << endl;
                    to_flip.push_back({triangles[pair.a].t[k], pair.a});
                }
            }

            // Ajouter les 2 arrêtes autour de b
            for (int k = 0; k < 3; k++)
            {
                if (triangles[pair.b].t[k] != pair.a)
                {
                    cout << "Push : " << triangles[pair.b].t[k] << ", " << pair.b << endl;
                    to_flip.push_back({triangles[pair.b].t[k], pair.b});
                }
            }
        }
        else
        {
            cout << "=> Don't Flip" << endl;
        }
    }

    cout << endl << "Finished with " << flip_count << " flips" << endl;
}

void Mesh::draw()
{
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < (int)triangles.size(); i++)
    {
        if (triangles[i].is_box)
        {
            point grey = {0.2, 0.2, 0.2};
            glColor3dv(&grey.x);
        }
        else
        {
            double light = triangles[i].get_light(vertices);
            rgb c = colors.get(i, light);
            glColor3dv(&c.x);
        }

        if (i == 1)
        {
            point pink = {0.8, 0.0, 0.8};
            glColor3dv(&pink.x);
        }

        Triangle* triangle = &(triangles[i]);

        for (int j = 0; j < 3; j++)
        {
            Vertex* vertex = &vertices[triangle->v[j]];

            // specify vertex
            glVertex3d(vertex->x, vertex->y, vertex->z);
        }

    }
    glEnd();

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < (int)triangles.size(); i++)
    {
        point black = {0.0, 0.0, 0.0};
        glColor3dv(&black.x);
        Triangle* triangle = &(triangles[i]);

        for (int j = 0; j < 3; j++)
        {
            Vertex* vertex = &vertices[triangle->v[j]];

            // specify vertex
            glVertex3d(vertex->x, vertex->y, vertex->z);
        }

    }
    glEnd();
}

void Mesh::draw_voronoi()
{
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glBegin(GL_LINES);
    for (int i = 0; i < (int)triangles.size(); i++)
    {
        point black = {1.0, 1.0, 1.0};
        glColor3dv(&black.x);

        Triangle* triangle = &(triangles[i]);

        for (int j = 0; j < 3; j++)
        {
            Vertex centre = get_center(vertices[triangle->v[0]], vertices[triangle->v[1]], vertices[triangle->v[2]]);

            Triangle* adj_triangle = &triangles[triangle->t[j]];
            Vertex adj_centre = get_center(vertices[adj_triangle->v[0]], vertices[adj_triangle->v[1]], vertices[adj_triangle->v[2]]);

            // specify vertex
            if (!triangle->is_box)
            {
                glVertex3d(centre.x, centre.y, 0.001);
                if (adj_triangle->is_box)
                {
                    glVertex3d(adj_centre.x * 1000, adj_centre.y * 1000, 0.001);
                }
                else
                {
                    glVertex3d(adj_centre.x, adj_centre.y, 0.001);
                }
            }
        }

    }
    glEnd();
}
