#include "vector.h"
#include <assert.h>
#include <cmath>
#include <iostream>
#include "math_tools.h"
#include "math_def.h"

GeoVector2D::GeoVector2D()
{
    m_coord[0] = 0.0f;
    m_coord[1] = 0.0f;
}

GeoVector2D::GeoVector2D(const double x, const double y)
{
    m_coord[0] = x;
    m_coord[1] = y;
}

GeoVector2D::GeoVector2D(const GeoVector2D &v)
{
    m_coord[0] = v[0];
    m_coord[1] = v[1];
}

GeoVector2D::~GeoVector2D()
{
}

GeoVector2D &GeoVector2D::operator=(const GeoVector2D &v)
{
    if (&v == this)
    {
        return *this;
    }

    m_coord[0] = v[0];
    m_coord[1] = v[1];

    return *this;
}

double GeoVector2D::operator[](const unsigned int idx) const
{
    assert(idx < GeoVector2D::Size());

    return m_coord[idx];
}

double &GeoVector2D::operator[](const unsigned int idx)
{
    assert(idx < GeoVector2D::Size());

    return m_coord[idx];
}

GeoVector2D &GeoVector2D::operator+=(const GeoVector2D &v)
{
    m_coord[0] += v[0];
    m_coord[1] += v[1];

    return *this;
}

GeoVector2D GeoVector2D::operator+(const GeoVector2D &v)
{
    GeoVector2D ret;

    ret[0] = m_coord[0] + v[0];
    ret[1] = m_coord[1] + v[1];

    return ret;
}

GeoVector2D GeoVector2D::operator*(const double scale) const
{
    GeoVector2D ret;

    ret[0] = m_coord[0] * scale;
    ret[1] = m_coord[1] * scale;

    return ret;
}

GeoVector2D &GeoVector2D::operator*=(const double scale)
{
    m_coord[0] *= scale;
    m_coord[1] *= scale;

    return *this;
}

GeoVector2D GeoVector2D::operator-(const GeoVector2D &v) const
{
    GeoVector2D ret;

    ret[0] = m_coord[0] - v[0];
    ret[1] = m_coord[1] - v[1];

    return ret;
}

GeoVector2D &GeoVector2D::operator-=(const GeoVector2D &v)
{
    m_coord[0] -= v[0];
    m_coord[1] -= v[1];

    return *this;
}

bool GeoVector2D::operator==(const GeoVector2D &v)
{

    if (MathTools::IsZero(m_coord[0] - v[0]) &&
        MathTools::IsZero(m_coord[1] - v[1]))
    {
        return true;
    }

    return false;
}

double GeoVector2D::operator%(const GeoVector2D &v) const
{
    return m_coord[0] * v[0] + m_coord[1] * v[1];
}

void GeoVector2D::Normalize()
{
    double magnitude = Magnitude();

    m_coord[0] /= magnitude;
    m_coord[1] /= magnitude;
}

double GeoVector2D::Magnitude() const
{
    return sqrt(m_coord[0] * m_coord[0] + m_coord[1] * m_coord[1]);
}

double GeoVector2D::Magnitude2() const
{
    return m_coord[0] * m_coord[0] + m_coord[1] * m_coord[1];
}

GeoVector2D GeoVector2D::ProjectTo(const GeoVector2D &v)
{
    GeoVector2D tmp = v;
    tmp.Normalize();

    return tmp * (*this % tmp);
}

void GeoVector2D::Flatten(std::vector<float> &data) const
{
    data.push_back((float)m_coord[0]);
    data.push_back((float)m_coord[1]);
}

void GeoVector2D::Dump() const
{
    std::cout.precision(5);
    for (unsigned int i = 0; i < GeoVector2D::Size(); i++)
    {
        std::cout << m_coord[i] << "    ";
    }
    std::cout << std::endl;
}

unsigned int GeoVector2D::Size()
{
    return 2;
}

GeoVector3D::GeoVector3D()
{
    m_coord[0] = 0.0f;
    m_coord[1] = 0.0f;
    m_coord[2] = 0.0f;
}

GeoVector3D::GeoVector3D(const double x, const double y, const double z)
{
    m_coord[0] = x;
    m_coord[1] = y;
    m_coord[2] = z;
}

GeoVector3D::GeoVector3D(const GeoVector3D &v)
{
    m_coord[0] = v[0];
    m_coord[1] = v[1];
    m_coord[2] = v[2];
}

GeoVector3D::~GeoVector3D()
{
}

GeoVector3D &GeoVector3D::operator=(const GeoVector3D &v)
{
    if (&v == this)
    {
        return *this;
    }

    m_coord[0] = v[0];
    m_coord[1] = v[1];
    m_coord[2] = v[2];

    return *this;
}

double GeoVector3D::operator[](const unsigned int idx) const
{
    assert(idx < GeoVector3D::Size());

    return m_coord[idx];
}

double &GeoVector3D::operator[](const unsigned int idx)
{
    assert(idx < GeoVector3D::Size());

    return m_coord[idx];
}

GeoVector3D &GeoVector3D::operator+=(const GeoVector3D &v)
{
    m_coord[0] += v[0];
    m_coord[1] += v[1];
    m_coord[2] += v[2];

    return *this;
}

GeoVector3D GeoVector3D::operator+(const GeoVector3D &v)
{
    GeoVector3D ret;

    ret[0] = m_coord[0] + v[0];
    ret[1] = m_coord[1] + v[1];
    ret[2] = m_coord[2] + v[2];

    return ret;
}

GeoVector3D GeoVector3D::operator*(const double scale) const
{
    GeoVector3D ret;

    ret[0] = m_coord[0] * scale;
    ret[1] = m_coord[1] * scale;
    ret[2] = m_coord[2] * scale;

    return ret;
}

GeoVector3D &GeoVector3D::operator*=(const double scale)
{
    m_coord[0] *= scale;
    m_coord[1] *= scale;
    m_coord[2] *= scale;

    return *this;
}

GeoVector3D GeoVector3D::operator-(const GeoVector3D &v) const
{
    GeoVector3D ret;

    ret[0] = m_coord[0] - v[0];
    ret[1] = m_coord[1] - v[1];
    ret[2] = m_coord[2] - v[2];

    return ret;
}

GeoVector3D &GeoVector3D::operator-=(const GeoVector3D &v)
{
    m_coord[0] -= v[0];
    m_coord[1] -= v[1];
    m_coord[2] -= v[2];

    return *this;
}

double GeoVector3D::operator%(const GeoVector3D &v) const
{
    return m_coord[0] * v[0] + m_coord[1] * v[1] + m_coord[2] * v[2];
}

GeoVector3D GeoVector3D::operator*(const GeoVector3D &v) const
{
    GeoVector3D ret;

    ret[0] = m_coord[1] * v[2] - m_coord[2] * v[1];
    ret[1] = m_coord[2] * v[0] - m_coord[0] * v[2];
    ret[2] = m_coord[0] * v[1] - m_coord[1] * v[0];

    return ret;
}

bool GeoVector3D::operator==(const GeoVector3D &v)
{

    if (MathTools::IsZero(m_coord[0] - v[0]) &&
        MathTools::IsZero(m_coord[1] - v[1]) &&
        MathTools::IsZero(m_coord[2] - v[2]))
    {
        return true;
    }

    return false;
}

void GeoVector3D::Normalize()
{
    double magnitude = Magnitude();

    m_coord[0] /= magnitude;
    m_coord[1] /= magnitude;
    m_coord[2] /= magnitude;
}

double GeoVector3D::Magnitude() const
{
    return sqrt(m_coord[0] * m_coord[0] + m_coord[1] * m_coord[1] + m_coord[2] * m_coord[2]);
}

double GeoVector3D::Magnitude2() const
{
    return m_coord[0] * m_coord[0] + m_coord[1] * m_coord[1] + m_coord[2] * m_coord[2];
}

GeoVector3D GeoVector3D::ProjectTo(const GeoVector3D &v)
{
    GeoVector3D tmp = v;
    tmp.Normalize();

    return tmp * (*this % tmp);
}

void GeoVector3D::Flatten(std::vector<float> &data) const
{
    data.push_back((float)m_coord[0]);
    data.push_back((float)m_coord[1]);
    data.push_back((float)m_coord[2]);
}

void GeoVector3D::Dump() const
{
    std::cout.precision(5);
    for (unsigned int i = 0; i < GeoVector3D::Size(); i++)
    {
        std::cout << m_coord[i] << "    ";
    }
    std::cout << std::endl;
}

unsigned int GeoVector3D::Size()
{
    return 3;
}

GeoVector4D::GeoVector4D()
{
    m_coord[0] = 0.0f;
    m_coord[1] = 0.0f;
    m_coord[2] = 0.0f;
    m_coord[3] = 0.0f;
}

GeoVector4D::GeoVector4D(const double x, const double y, const double z, const double w)
{
    m_coord[0] = x;
    m_coord[1] = y;
    m_coord[2] = z;
    m_coord[3] = w;
}

GeoVector4D::GeoVector4D(const GeoVector4D &v)
{
    m_coord[0] = v[0];
    m_coord[1] = v[1];
    m_coord[2] = v[2];
    m_coord[3] = v[3];
}

GeoVector4D::GeoVector4D(const GeoVector3D &v, const double w)
{
    m_coord[0] = v[0];
    m_coord[1] = v[1];
    m_coord[2] = v[2];
    m_coord[3] = w;
}

GeoVector4D::~GeoVector4D()
{
}

GeoVector4D &GeoVector4D::operator=(const GeoVector4D &v)
{
    if (&v == this)
    {
        return *this;
    }

    m_coord[0] = v[0];
    m_coord[1] = v[1];
    m_coord[2] = v[2];
    m_coord[3] = v[3];

    return *this;
}

double GeoVector4D::operator[](const unsigned int idx) const
{
    assert(idx < GeoVector4D::Size());

    return m_coord[idx];
}

double &GeoVector4D::operator[](const unsigned int idx)
{
    assert(idx < GeoVector4D::Size());

    return m_coord[idx];
}

GeoVector4D &GeoVector4D::operator+=(const GeoVector4D &v)
{
    m_coord[0] += v[0];
    m_coord[1] += v[1];
    m_coord[2] += v[2];
    m_coord[3] += v[3];

    return *this;
}

GeoVector4D GeoVector4D::operator+(const GeoVector4D &v)
{
    GeoVector4D ret;

    ret[0] = m_coord[0] + v[0];
    ret[1] = m_coord[1] + v[1];
    ret[2] = m_coord[2] + v[2];
    ret[3] = m_coord[3] + v[3];

    return ret;
}

GeoVector4D GeoVector4D::operator*(const double scale) const
{
    GeoVector4D ret;

    ret[0] = m_coord[0] * scale;
    ret[1] = m_coord[1] * scale;
    ret[2] = m_coord[2] * scale;
    ret[3] = m_coord[3] * scale;

    return ret;
}

GeoVector4D &GeoVector4D::operator*=(const double scale)
{
    m_coord[0] *= scale;
    m_coord[1] *= scale;
    m_coord[2] *= scale;
    m_coord[3] *= scale;

    return *this;
}

GeoVector4D GeoVector4D::operator-(const GeoVector4D &v) const
{
    GeoVector4D ret;

    ret[0] = m_coord[0] - v[0];
    ret[1] = m_coord[1] - v[1];
    ret[2] = m_coord[2] - v[2];
    ret[3] = m_coord[3] - v[3];

    return ret;
}

GeoVector4D &GeoVector4D::operator-=(const GeoVector4D &v)
{
    m_coord[0] -= v[0];
    m_coord[1] -= v[1];
    m_coord[2] -= v[2];
    m_coord[3] -= v[3];

    return *this;
}

bool GeoVector4D::operator==(const GeoVector4D &v)
{

    if (MathTools::IsZero(m_coord[0] - v[0]) &&
        MathTools::IsZero(m_coord[1] - v[1]) &&
        MathTools::IsZero(m_coord[2] - v[2]) &&
        MathTools::IsZero(m_coord[3] - v[3]))
    {
        return true;
    }

    return false;
}

void GeoVector4D::Normalize()
{
    double magnitude = Magnitude();

    m_coord[0] /= magnitude;
    m_coord[1] /= magnitude;
    m_coord[2] /= magnitude;
    m_coord[3] /= magnitude;
}

double GeoVector4D::Magnitude() const
{
    return sqrt(m_coord[0] * m_coord[0] + m_coord[1] * m_coord[1] + m_coord[2] * m_coord[2] + m_coord[3] * m_coord[3]);
}

double GeoVector4D::Magnitude2() const
{
    return m_coord[0] * m_coord[0] + m_coord[1] * m_coord[1] + m_coord[2] * m_coord[2] + m_coord[3] * m_coord[3];
}

void GeoVector4D::Flatten(std::vector<float> &data) const
{
    data.push_back((float)m_coord[0]);
    data.push_back((float)m_coord[1]);
    data.push_back((float)m_coord[2]);
    data.push_back((float)m_coord[3]);
}

void GeoVector4D::Dump() const
{
    std::cout.precision(5);
    for (unsigned int i = 0; i < GeoVector4D::Size(); i++)
    {
        std::cout << m_coord[i] << "    ";
    }
    std::cout << std::endl;
}

unsigned int GeoVector4D::Size()
{
    return 4;
}

GeoVector::GeoVector()
    : m_coord(nullptr), m_dim(0)
{
}

GeoVector::GeoVector(const unsigned int len, const double *data /* = nullptr */)
{
    assert(len > 0);

    m_dim = len;
    m_coord = new double[m_dim];
    memset(m_coord, 0, sizeof(double) * m_dim);

    if (data)
    {
        for (unsigned int i = 0; i < m_dim; i++)
        {
            m_coord[i] = data[i];
        }
    }
}

GeoVector::GeoVector(const GeoVector &v)
{
    assert((v.m_dim > 0) && v.m_coord);

    m_dim = v.m_dim;
    m_coord = new double[m_dim];

    for (unsigned int i = 0; i < m_dim; i++)
    {
        m_coord[i] = v.m_coord[i];
    }
}

GeoVector::~GeoVector()
{
    Clear();
}

GeoVector &GeoVector::operator=(const GeoVector &v)
{
    if (this == &v)
    {
        return *this;
    }

    Clear();

    m_dim = v.m_dim;
    m_coord = new double[m_dim];

    for (unsigned int i = 0; i < m_dim; i++)
    {
        m_coord[i] = v.m_coord[i];
    }

    return *this;
}

double GeoVector::operator[](const unsigned int idx) const
{
    assert(idx < m_dim);

    return m_coord[idx];
}

double &GeoVector::operator[](const unsigned int idx)
{
    assert(idx < m_dim);

    return m_coord[idx];
}

GeoVector &GeoVector::operator+=(const GeoVector &v)
{
    assert(m_dim == v.Dim());

    for (unsigned int i = 0; i < m_dim; i++)
    {
        m_coord[i] += v[i];
    }

    return *this;
}

GeoVector GeoVector::operator+(const GeoVector &v)
{
    GeoVector ret(m_dim);

    for (unsigned int i = 0; i < m_dim; i++)
    {
        ret[i] = m_coord[i] + v[i];
    }

    return ret;
}

GeoVector GeoVector::operator*(const double scale) const
{
    GeoVector ret(m_dim);

    for (unsigned int i = 0; i < m_dim; i++)
    {
        ret[i] = m_coord[i] * scale;
    }

    return ret;
}

GeoVector &GeoVector::operator*=(const double scale)
{
    for (unsigned int i = 0; i < m_dim; i++)
    {
        m_coord[i] *= scale;
    }

    return *this;
}

GeoVector GeoVector::operator-(const GeoVector &v) const
{
    GeoVector ret(m_dim);

    for (unsigned int i = 0; i < m_dim; i++)
    {
        ret[i] = m_coord[i] - v[i];
    }

    return ret;
}

GeoVector &GeoVector::operator-=(const GeoVector &v)
{
    assert(m_dim == v.Dim());

    for (unsigned int i = 0; i < m_dim; i++)
    {
        m_coord[i] -= v[i];
    }

    return *this;
}

void GeoVector::Normalize()
{
    double magnitude = Magnitude();
    for (unsigned int i = 0; i < m_dim; i++)
    {
        m_coord[i] /= magnitude;
    }
}

double GeoVector::Magnitude() const
{
    double sum = Magnitude2();

    return sqrt(sum);
}

double GeoVector::Magnitude2() const
{
    double sum = 0.0f;
    for (unsigned int i = 0; i < m_dim; i++)
    {
        sum += m_coord[i] * m_coord[i];
    }

    return sum;
}

void GeoVector::Flatten(std::vector<float> &data) const
{
    for (unsigned int i = 0; i < m_dim; i++)
    {
        data.push_back((float)m_coord[i]);
    }
}

unsigned int GeoVector::Dim() const
{
    return m_dim;
}

void GeoVector::Dump() const
{
    std::cout.precision(5);
    for (unsigned int i = 0; i < m_dim; i++)
    {
        std::cout << m_coord[i] << "    ";
    }
    std::cout << std::endl;
}

void GeoVector::Clear()
{
    m_dim = 0;
    SAFE_DELETE_ARRAY(m_coord);
}