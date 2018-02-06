#ifndef COLORS_H
#define COLORS_H

#include <string>
#include <vector>

struct rgb {
    double x, y, z;
};

class Colors
{
public:
    std::vector<rgb> list;

    Colors();
    void add(int x, int y, int z);
    rgb get(int index);
    rgb get(int index, double light);
};

#endif // COLORS_H
