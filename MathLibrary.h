#ifndef MATH_LIBRARY_H
#define MATH_LIBRARY_H

#include <cmath>
#include <cassert>
#include <iostream>

constexpr float tolerance = 0.0001f;

//------------------------------------------------------------------------------------------------------
//                                  Base Vect
//------------------------------------------------------------------------------------------------------
template<int vectSize> struct Vect
{
    Vect() = default;
    float &operator[](const int i) { assert(i >= 0 && i < vectSize); return data[i]; }
    float  operator[](const int i) const { assert(i >= 0 && i < vectSize); return data[i]; }
    Vect<vectSize> norm(){norm(*this);};
    float data[vectSize] = { 0.0f };
};

//------------------------------------------------------------------------------------------------------
//                                  Operator Overloading
//------------------------------------------------------------------------------------------------------
template<int vectSize> Vect<vectSize> operator *(const Vect<vectSize> &lhs, const Vect<vectSize> &rhs)
{
    Vect<vectSize> res = lhs;

    for (int i = vectSize; i--;)
    {
        res[i] *= rhs[i];
    }

    return res;
}

template<int vectSize> Vect<vectSize> operator *=(const Vect<vectSize> &lhs, const Vect<vectSize> &rhs)
{
    for (int i = vectSize; i--;)
    {
        lhs[i] *= rhs[i];
    }

    return lhs;
}

template<int vectSize> Vect<vectSize> operator +(const Vect<vectSize> &lhs, const Vect<vectSize> &rhs)
{
    Vect<vectSize> res = lhs;

    for (int i = vectSize; i--;)
    {
        res[i] += rhs[i];
    }

    return res;
}

template<int vectSize> Vect<vectSize> operator -(const Vect<vectSize> &lhs, const Vect<vectSize> &rhs)
{
    Vect<vectSize> res = lhs;

    for (int i = vectSize; i--;)
    {
        res[i] -= rhs[i];
    }

    return res;
}

template<int vectSize> Vect<vectSize> operator *(const Vect<vectSize> &lhs, const float scalar)
{
    Vect<vectSize> res = lhs;

    for (int i = vectSize; i--;)
    {
        res[i] *= scalar;
    }

    return res;
}

template<int vectSize> Vect<vectSize> operator-(Vect<vectSize> &lhs)
{
    return lhs*(-1.f);
}

template<int vectSize> Vect<vectSize> operator-(const Vect<vectSize> &lhs)
{
    return lhs * (-1.f);
}

//------------------------------------------------------------------------------------------------------
//                                  Math Functions
//------------------------------------------------------------------------------------------------------
template<int vectSize> float t_dot(const Vect<vectSize> &lhs, const Vect<vectSize> &rhs)
{
    float a = 0;
    for (int i = vectSize; i--;)
    {
        a += lhs[i] * rhs[i];
    }

    return a;
}

template<int vectSize> float t_dot(const Vect<vectSize> &lhs, const float f)
{
    float a = 0;
    for (int i = vectSize; i--;)
    {
        a += lhs[i] * f;
    }

    return a;
}

template<int vectSize> void t_norm(Vect<vectSize> &lhs)
{
    float mag = std::sqrtf(t_dot<vectSize>(lhs, lhs));
    if(fabsf(mag) <= tolerance) 
    {
        std::cout << "Division by 0 in t_norm" << std::endl;
        return;
    }
    for (int i = vectSize; i--;)
    {
        lhs[i] /= mag;
    }
}

template<int vectSize> Vect<vectSize> t_getNorm(const Vect<vectSize> &lhs)
{
    Vect<vectSize> res = lhs;
    float mag = std::sqrtf(t_dot<vectSize>(lhs, lhs));
    for (int i = vectSize; i--;)
    {
        res[i] /= mag;
    }

    return res;
}

template<int vectSize> float t_mag(const Vect<vectSize> &lhs)
{
    return std::sqrtf(t_dot<vectSize>(lhs, lhs));
}

//------------------------------------------------------------------------------------------------------
//                                  Custom Vect Types
//------------------------------------------------------------------------------------------------------

template<> struct Vect<2>
{
    Vect() : x(0.0f), y(0.0f)
    {
    }
    Vect(const float x, const float y) : x(x), y(y)
    {
    }

    float &operator[](const int i) { assert(i >= 0 && i < 2); return i ? y : x; }
    float  operator[](const int i) const { assert(i >= 0 && i < 2); return i ? y : x; }
    float dot(const Vect<2> &v) const{return t_dot<2>(*this, v);}
    float dot(const float f) const { return t_dot<2>(*this, f);}
    void norm() { return t_norm(*this);}
    Vect<2> getNorm() const { return t_getNorm(*this);}
    float mag() const { return t_mag(*this);}

    float x{};
    float y{};
};

template<> struct Vect<3>
{
    Vect() : x(0.0f), y(0.0f), z(0.0f)
    {
    }
    Vect(const float x, const float y, const float z) : x(x), y(y), z(z)
    {
    }

    float &operator[](const int i) { assert(i >= 0 && i < 3); return i ? (1 == i ? y : z) : x; }
    float  operator[](const int i) const { assert(i >= 0 && i < 3); return i ? (1 == i ? y : z) : x; }
    float dot(const Vect<3> &v) const{ return t_dot<3>(*this, v); }
    float dot(const float f) const { return t_dot<3>(*this, f); }
    void norm() { return t_norm(*this); }
    Vect<3> getNorm() const{ return t_getNorm(*this);}
    float mag() const{ return t_mag(*this);}

    float x{};
    float y{};
    float z{};
};

template<> struct Vect<4>
{
    Vect() : x(0.0f), y(0.0f), z(0.0f), w(0.0f)
    {
    }
    Vect(const float x, const float y, const float z, const float w) : x(x), y(y), z(z), w(w)
    {
    }

    float &operator[](const int i)
    {
        assert(i >= 0 && i < 4);
        switch (i)
        {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            case 3: return w;
            default: assert(false);
        }
    }
    float  operator[](const int i) const
    {
        assert(i >= 0 && i < 4);
        switch (i)
        {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            case 3: return w;
            default: assert(false);
        }
    }
    float dot(const Vect<4> &v) const { return t_dot<4>(*this, v); }
    float dot(const float f) const { return t_dot<4>(*this, f); }
    void norm() { return t_norm(*this); }
    Vect<4> getNorm() const { return t_getNorm(*this); }
    float mag() const { return t_mag(*this);}

    float x{};
    float y{};
    float z{};
    float w{};
};

typedef Vect<2> Vec2f;
typedef Vect<3> Vec3f;
typedef Vect<4> Vec4f;

#endif // !MATH_LIBRARY_H
