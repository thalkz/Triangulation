#include "colors.h"
#include <vector>
#include <string>
#include <math.h>

Colors::Colors()
{
    add(85, 239, 196);  // 0 = green
    add(129, 236, 236); // 1 = light-blue
    add(116, 185, 255); // 2 = blue
    add(162, 155, 254); // 3 = mauve
    add(223, 230, 233); // 4 = light-grey
    add(255, 234, 167); // 5 = yellow
    add(255, 118, 117); // 6 = red
    add(253, 203, 110); // 7 = orange
    add(99, 110, 114);  // 8 = grey
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
