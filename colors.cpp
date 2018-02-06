#include "colors.h"
#include <vector>
#include <string>
#include <math.h>

Colors::Colors()
{
    add(85, 239, 196);
    add(129, 236, 236);
    add(116, 185, 255);
    add(162, 155, 254);
    add(223, 230, 233);
    add(255, 234, 167);
    add(255, 234, 167);
    add(255, 118, 117);
    add(253, 203, 110);
    add(99, 110, 114);
}

void Colors::add(int x, int y, int z)
{
    list.push_back({(double) x / 255, (double) y / 255, (double) z / 255});
}

rgb Colors::get(int index)
{
    return list[index % (int) list.size()];
}

rgb Colors::get(int index, double light)
{
    rgb color = list[index % (int) list.size()];
    return {color.x * light, color.y * light, color.z * light};
}
