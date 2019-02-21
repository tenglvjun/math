#include "math_tools.h"
#include "math_def.h"
#include <cmath>

MathTools::MathTools()
{
}

MathTools::~MathTools()
{
}

double MathTools::Radia2Degree(const double r)
{
    return r * 180 / PI;
}

double MathTools::Degree2Radian(const double d)
{
    return ((d * PI) / (double)180);
}

bool MathTools::IsZero(const double v)
{
    return (fabs(v) < EPSILON) ? true : false;
}
