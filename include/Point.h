/* Little struct that describes a point in a cartesian xy plane */

#ifndef POINT_H
#define POINT_H

struct Point
{
public:
    Point(double x, double y);
    ~Point();
    double x, y;
};
#endif // POINT_H
