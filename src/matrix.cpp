#include "matrix.h"
#include <memory.h>
#include <assert.h>
#include <cmath>
#include <iostream>
#include "math_tools.h"
#include "math_def.h"

GeoMatrix::GeoMatrix(const unsigned int row, const unsigned int col)
    : m_data(nullptr), m_row(0), m_col(0)
{
    Init(row, col);
}

GeoMatrix::GeoMatrix(const unsigned int row, const unsigned int col, const double *data)
    : m_data(nullptr), m_row(0), m_col(0)
{
    Init(row, col);

    for (unsigned int i = 0; i < m_row; i++)
    {
        for (unsigned int j = 0; j < m_col; j++)
        {
            m_data[i][j] = data[m_col * i + j];
        }
    }
}

GeoMatrix::GeoMatrix(const GeoMatrix &m)
{
    Init(m.m_row, m.m_col);

    for (unsigned int row = 0; row < m_row; row++)
    {
        for (unsigned int col = 0; col < m_col; col++)
        {
            m_data[row][col] = m[row][col];
        }
    }
}

GeoMatrix::~GeoMatrix()
{
    Clear();
}

GeoMatrix &GeoMatrix::operator=(const GeoMatrix &m)
{
    if (&m == this)
    {
        return *this;
    }

    Clear();
    Init(m.m_row, m.m_col);

    for (unsigned int row = 0; row < m_row; row++)
    {
        for (unsigned int col = 0; col < m_col; col++)
        {
            m_data[row][col] = m[row][col];
        }
    }

    return *this;
}

double *GeoMatrix::operator[](const unsigned int idx) const
{
    assert(m_data && (idx < m_row));

    return m_data[idx];
}

double *GeoMatrix::operator[](const unsigned int idx)
{
    assert(m_data && (idx < m_row));

    return m_data[idx];
}

GeoVector3D GeoMatrix::operator*(const GeoVector3D &v) const
{
    assert((GeoVector3D::Size() == m_col) && (m_col == m_row));

    GeoVector3D ret;
    ret[0] = m_data[0][0] * v[0] + m_data[0][1] * v[1] + m_data[0][2] * v[2];
    ret[1] = m_data[1][0] * v[0] + m_data[1][1] * v[1] + m_data[1][2] * v[2];
    ret[2] = m_data[2][0] * v[0] + m_data[2][1] * v[1] + m_data[2][2] * v[2];

    return ret;
}

GeoVector2D GeoMatrix::operator*(const GeoVector2D &v) const
{
    assert((GeoVector2D::Size() == m_col) && (m_col == m_row));

    GeoVector2D ret;
    ret[0] = m_data[0][0] * v[0] + m_data[0][1] * v[1];
    ret[1] = m_data[1][0] * v[0] + m_data[1][1] * v[1];

    return ret;
}

GeoVector4D GeoMatrix::operator*(const GeoVector4D &v) const
{
    assert((GeoVector4D::Size() == m_col) && (m_col == m_row));

    GeoVector4D ret;

    ret[0] = m_data[0][0] * v[0] + m_data[0][1] * v[1] + m_data[0][2] * v[2] + m_data[0][3] * v[3];
    ret[1] = m_data[1][0] * v[0] + m_data[1][1] * v[1] + m_data[1][2] * v[2] + m_data[1][3] * v[3];
    ret[2] = m_data[2][0] * v[0] + m_data[2][1] * v[1] + m_data[2][2] * v[2] + m_data[2][3] * v[3];
    ret[3] = m_data[3][0] * v[0] + m_data[3][1] * v[1] + m_data[3][2] * v[2] + m_data[3][3] * v[3];

    return ret;
}

GeoMatrix GeoMatrix::operator*(const GeoMatrix &m) const
{
    assert(m_col == m.Rows());

    GeoMatrix ret(m_row, m.Cols());

    for (unsigned int i = 0; i < m_row; i++)
    {
        for (unsigned int j = 0; j < m.Cols(); j++)
        {
            for (unsigned int k = 0; k < m_col; k++)
            {
                ret[i][j] += (m_data[i][k] * m[k][j]);
            }
        }
    }

    return ret;
}

void GeoMatrix::operator*=(const double s)
{
    for (unsigned int i = 0; i < m_row; i++)
    {
        for (unsigned int j = 0; j < m_col; j++)
        {
            m_data[i][j] *= s;
        }
    }
}

GeoMatrix GeoMatrix::operator*(const double s)
{
    GeoMatrix ret(m_row, m_col);

    for (unsigned int i = 0; i < m_row; i++)
    {
        for (unsigned int j = 0; j < m_col; j++)
        {
            ret[i][j] *= s;
        }
    }

    return ret;
}

void GeoMatrix::operator+=(const GeoMatrix &m)
{
    assert((m_col == m.Cols()) && (m_row == m.Rows()));

    for (unsigned int i = 0; i < m_row; i++)
    {
        for (unsigned int j = 0; j < m_col; j++)
        {
            m_data[i][j] += m[i][j];
        }
    }
}

void GeoMatrix::SetIdentity()
{
    assert((m_row == m_col) && (m_row > 0));

    for (unsigned int row = 0; row < m_row; row++)
    {
        for (unsigned int col = 0; col < m_col; col++)
        {
            if (row == col)
            {
                m_data[row][col] = 1.0f;
            }
            else
            {
                m_data[row][col] = 0.0f;
            }
        }
    }
}

void GeoMatrix::Zeros()
{
    for (unsigned int row = 0; row < m_row; row++)
    {
        m_data[row] = new double[m_col];
        for (unsigned int col = 0; col < m_col; col++)
        {
            m_data[row][col] = 0.0f;
        }
    }
}

void GeoMatrix::Resharp(const unsigned int row, const unsigned int col)
{
    if (m_row == row && m_col == col)
    {
        Zeros();
        return;
    }

    Clear();
    Init(row, col);
}

void GeoMatrix::Flatten(std::vector<float> &data) const
{
    for (unsigned int j = 0; j < m_col; j++)
    {
        for (unsigned int i = 0; i < m_row; i++)
        {
            data.push_back((float)(m_data[i][j]));
        }
    }
}

GeoMatrix GeoMatrix::SubMatrix(const unsigned int sRow, const unsigned int eRow, const unsigned int sCol, const unsigned int eCol)
{
    GeoMatrix m(eRow - sRow, eCol - sCol);

    for (unsigned int i = sRow; i < eRow; i++)
    {
        for (unsigned int j = sCol; j < eCol; j++)
        {
            m[i - sRow][j - sCol] = m_data[i][j];
        }
    }

    return m;
}

unsigned int GeoMatrix::Rows() const
{
    return m_row;
}

unsigned int GeoMatrix::Cols() const
{
    return m_col;
}

bool GeoMatrix::LUDecompose(GeoMatrix &up, GeoMatrix &low) const
{
    if (!IsSquare())
    {
        return false;
    }

    up.Resharp(m_row, m_col);
    low.Resharp(m_row, m_col);

    unsigned int n = m_row;
    double sum;

    for (unsigned int i = 0; i < n; i++)
    {
        // Upper Triangular
        for (unsigned int k = i; k < n; k++)
        {
            // Summation of L(i, j) * U(j, k)
            sum = 0.0f;
            for (unsigned int j = 0; j < i; j++)
                sum += (low[i][j] * up[j][k]);

            // Evaluating U(i, k)
            up[i][k] = m_data[i][k] - sum;

            if ((i == k) && (MathTools::IsZero(up[i][k])))
            {
                return false;
            }
        }

        // Lower Triangular
        for (unsigned int k = i; k < n; k++)
        {
            if (i == k)
            {
                low[i][i] = 1; // Diagonal as 1
            }
            else
            {

                // Summation of L(k, j) * U(j, i)
                sum = 0;
                for (unsigned int j = 0; j < i; j++)
                {
                    sum += (low[k][j] * up[j][i]);
                }

                // Evaluating L(k, i)
                low[k][i] = (m_data[k][i] - sum) / up[i][i];
            }
        }
    }

    return true;
}

double GeoMatrix::Det() const
{
    GeoMatrix up(m_row, m_col);
    GeoMatrix low(m_row, m_col);

    if (!LUDecompose(up, low))
    {
        return 0;
    }

    double det = 1.0f;

    for (unsigned int i = 0; i < m_row; i++)
    {
        det *= up[i][i];
    }

    return det;
}

bool GeoMatrix::Inverse(GeoMatrix &inverse) const
{
    GeoMatrix up(m_row, m_col);
    GeoMatrix low(m_row, m_col);

    if (!LUDecompose(up, low))
    {
        return false;
    }

    inverse.Resharp(m_row, m_col);

    for (unsigned int i = 0; i < m_row; i++)
    {
        GeoVector b(m_col);

        for (unsigned int j = 0; j < m_col; j++)
        {
            b[j] = (i == j) ? 1 : 0;
        }

        GeoVector ret = GeoMatrix::SolveLinearEquation(up, low, b);

        inverse.SetVector(i, ret, false);
    }

    return true;
}

bool GeoMatrix::IsSquare() const
{
    assert((m_col > 0) && (m_col > 0));

    return (m_col == m_row);
}

GeoMatrix &GeoMatrix::Transpose()
{
    assert((m_col > 0) && (m_col > 0));

    GeoMatrix tmp(m_col, m_row);
    Transpose(tmp);

    (*this) = (tmp);

    return *this;
}

void GeoMatrix::Transpose(GeoMatrix &transpose) const
{
    transpose.Resharp(m_col, m_row);

    for (unsigned int row = 0; row < m_row; row++)
    {
        for (unsigned int col = 0; col < m_col; col++)
        {
            transpose[col][row] = m_data[row][col];
        }
    }
}

void GeoMatrix::SetVector(const unsigned int idx, const GeoVector &v, bool isRow)
{
    if (isRow)
    {
        assert(m_col == v.Dim());
        for (unsigned int col = 0; col < m_col; col++)
        {
            m_data[idx][col] = v[col];
        }
    }
    else
    {
        assert(m_row == v.Dim());
        for (unsigned int row = 0; row < m_row; row++)
        {
            m_data[row][idx] = v[row];
        }
    }
}

bool GeoMatrix::SolveLinearEquation(const GeoVector &b, GeoVector &x) const
{
    assert(b.Dim() == x.Dim());

    GeoMatrix up(m_row, m_col);
    GeoMatrix low(m_row, m_col);

    if (!LUDecompose(up, low))
    {
        return false;
    }

    x = GeoMatrix::SolveLinearEquation(up, low, b);

    return true;
}

void GeoMatrix::Dump() const
{
    std::cout.precision(5);
    for (unsigned int i = 0; i < m_row; i++)
    {
        for (unsigned int j = 0; j < m_col; j++)
        {
            std::cout << m_data[i][j] << "  ";
        }

        std::cout << std::endl;
    }
}

void GeoMatrix::Clear()
{
    for (unsigned int row = 0; row < m_row; row++)
    {
        SAFE_DELETE_ARRAY(m_data[row]);
    }

    SAFE_DELETE(m_data);

    m_row = 0;
    m_col = 0;
}

void GeoMatrix::Init(const unsigned int row, const unsigned int col)
{
    m_row = row;
    m_col = col;

    m_data = new double *[m_row];

    Zeros();
}

GeoMatrix GeoMatrix::TranslateMatrix(const GeoVector3D &trans)
{
    GeoMatrix matrix(4, 4);

    matrix.SetIdentity();

    matrix[0][3] = trans[0];
    matrix[1][3] = trans[1];
    matrix[2][3] = trans[2];

    return matrix;
}

GeoMatrix GeoMatrix::TranslateMatrix(const GeoVector4D &trans)
{
    GeoMatrix matrix(4, 4);

    matrix.SetIdentity();

    matrix[0][3] = trans[0];
    matrix[1][3] = trans[1];
    matrix[2][3] = trans[2];

    return matrix;
}

GeoMatrix GeoMatrix::RotateMatrix(const double angle, const GeoVector3D &axis)
{
    GeoMatrix matrix(4, 4);

    double c = cos(angle);
    double s = sin(angle);

    matrix[0][0] = c + (1 - c) * axis[0] * axis[0];
    matrix[0][1] = (1 - c) * axis[0] * axis[1] - s * axis[2];
    matrix[0][2] = (1 - c) * axis[0] * axis[2] + s * axis[1];

    matrix[1][0] = (1 - c) * axis[0] * axis[1] + s * axis[2];
    matrix[1][1] = c + (1 - c) * axis[1] * axis[1];
    matrix[1][2] = (1 - c) * axis[1] * axis[2] - s * axis[0];

    matrix[2][0] = (1 - c) * axis[0] * axis[2] - s * axis[1];
    matrix[2][1] = (1 - c) * axis[1] * axis[2] + s * axis[0];
    matrix[2][2] = c + (1 - c) * axis[2] * axis[2];

    matrix[3][3] = 1.0f;

    return matrix;
}

GeoMatrix GeoMatrix::ScaleMatrix(const double s)
{
    GeoMatrix ret(4, 4);
    ret.SetIdentity();

    ret[0][0] = s;
    ret[1][1] = s;
    ret[2][2] = s;

    return ret;
}

GeoVector GeoMatrix::SolveLinearEquation(const GeoMatrix &up, const GeoMatrix &low, const GeoVector &b)
{
    double sum;

    // Solve Lower
    GeoVector y(b.Dim());
    for (int row = 0; row < (int)low.Rows(); row++)
    {
        sum = 0.0f;
        for (int col = 0; col < row; col++)
        {
            sum += low[row][col] * y[col];
        }
        y[row] = (b[row] - sum) / low[row][row];
    }

    GeoVector ret(y.Dim());
    //Solve Upper
    for (int row = (int)(up.Rows() - 1); row >= 0; row--)
    {
        sum = 0.0f;
        for (int col = (int)(up.Cols() - 1); col > row; col--)
        {
            sum += up[row][col] * ret[col];
        }
        ret[row] = (y[row] - sum) / up[row][row];
    }

    return ret;
}

GeoVector4D operator*(const GeoVector4D &v, const GeoMatrix &m)
{
    assert((GeoVector4D::Size() == m.Cols()) && (m.IsSquare()));

    GeoVector4D ret;

    ret[0] = m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2] + m[0][3] * v[3];
    ret[1] = m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2] + m[1][3] * v[3];
    ret[2] = m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2] + m[2][3] * v[3];
    ret[3] = m[3][0] * v[0] + m[3][1] * v[1] + m[3][2] * v[2] + m[3][3] * v[3];

    return ret;
}

GeoVector3D operator*(const GeoVector3D &v, const GeoMatrix &m)
{
    assert((GeoVector3D::Size() == m.Cols()) && (m.IsSquare()));

    GeoVector3D ret;

    ret[0] = m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2];
    ret[1] = m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2];
    ret[2] = m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2];

    return ret;
}

GeoVector2D operator*(const GeoVector2D &v, const GeoMatrix &m)
{
    assert((GeoVector2D::Size() == m.Cols()) && (m.IsSquare()));

    GeoVector2D ret;

    ret[0] = m[0][0] * v[0] + m[0][1] * v[1];
    ret[1] = m[1][0] * v[0] + m[1][1] * v[1];

    return ret;
}