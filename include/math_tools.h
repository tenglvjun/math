#if !defined(__MATH_TOOLS_HEAD_FILE__)
#define __MATH_TOOLS_HEAD_FILE__

class MathTools final
{
  public:
    MathTools();
    ~MathTools();

  public:
    static double Radia2Degree(const double r);
    static double Degree2Radian(const double d);
    static bool IsZero(const double v);
};

#endif // __MATH_TOOLS_HEAD_FILE__
