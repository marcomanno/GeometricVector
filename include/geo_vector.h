
#pragma once

#include "geo_iter.h"
#include "geo_pow.h"

//#include <type_traits>
//#include <limits>
//#include <assert.h>
//#include <algorithm>
//#include <array>

namespace Geo
{
typedef double real;
template <size_t N, typename T> struct Vect_base
{
protected:
  T m_el[N];
  void set_null()
  {
    T * a = m_el;
    iterate_forw<N>::eval([&a](int i)
    {
      a[i] = static_cast<T> (0); } );
  }
  bool is_null()
  {
    T * a = m_el;
    return iterate_forw_until<N, false>::eval([&a](int i)
    {
      return a[i] == static_cast<T> (0);
    });
  }
  template<typename Arg1>
  static void assign(T * p, const Arg1 & a1)
  {
    p[0] = a1;
  }
  template<typename Arg1, typename ... Args>
  static void assign(T * p, const Arg1 & a1, const Args & ... a2)
  {
    p[0] = a1;
    assign(p + 1, a2...);
  }

public:
  const T& operator [] (int i) const
  {
    return m_el[i];
  }
  T len_sq() const
  {
    T v = 0;
    const T* a = m_el;
    iterate_forw<N>::eval([a, &v](int i)
    {
      v += gk_sq<T>(a[i]); } );
    return v;
  }
  T * coord() { return m_el; }
  const T * coord() const { return m_el; }
  typedef T mygeotype;
  enum
  {
    Dimension = N
  } ;
} ;
template <size_t N=3, typename T=real> class Vect : public Vect_base<N, T>
{
public:
  Vect()
  {
    /* Not init */
  }
  Vect(const T coo[N])
  {
    init(coo);
  }
  Vect(const Vect_base<N, T> & vb)
  {
    init(vb.coord());
  }
  Vect(const T& val)
  {
    fill(val);
  }
  template<typename... Args> explicit Vect(const Args & ... a)
  {
    static_assert(sizeof...(Args) == N, "Wrong parameter length!");
    this->assign(this->m_el, a...);
  }
  void init(const T * coo)
  {
    iterate_forw<N>::eval([this, &coo](int i)
    {
      m_el[i] = coo[i];
    } );
  }
  void fill(const T& val)
  {
    iterate_forw<N>::eval([this, &val](int i)
    {
      m_el[i] = val;
    });
  }
  Vect<N, T> & operator +=(const Vect_base<N, T> & a)
  {
    struct plus_eq { void operator()(T & a, const T & b) { a += b; }; };
    oper_self<plus_eq>(a);
    return *this;
  }
  Vect<N, T> & operator -=(const Vect_base<N, T> & a)
  {
    struct minus_eq { void operator()(T & a, const T & b) { a -= b; }; };
    oper_self<minus_eq>(a);
    return *this;
  }
  T & operator [] (int i)
  {
    return m_el[i];
  }
private:
  template <typename op> void oper_self(const Vect_base<N, T> & x)
  {
    T * a = this->m_el;
    const T * b = x.coord();
    iterate_forw<N>::eval([&a, &b](int i) { op()(a[i], b[i]); } );
  }
} ;
template <size_t N=3, typename T=real> class Vers : public Vect_base<N, T>
{
public:
  Vers()
  {
    /* Not init */
  }
  Vers(const Vect<N, T> & oth)
  {
    init(oth.coord());
  }
  template<typename... Args> explicit Vers(const Args & ... a)
  {
    static_assert(sizeof...(Args) == N, "Wrong parameter length!");
    this->assign(this->m_el, a...);
    normalize();
  }
  Vers(const T coo[N])
  {
    init(coo);
  }
  inline double init(const Vect<N, T> & v)
  {
    return init(v.coord()) > 0;
  }
  double init(const T coo[N])
  {
    T * p = this->m_el;
    iterate_forw<N>::eval([&p, coo](int i)
    {
      p[i] = coo[i]; } );
    return normalize();
  }
  T normalize()
  {
    T mod = 0, * p = this->m_el;
    iterate_forw<N>::eval([&mod, &p](int i)
    {
      mod += gk_sq<T>(p[i]);
    } );
    if (mod == 0)
      Vect_base<N, T>::set_null();
    else
    {
      mod = sqrt(mod);
      iterate_forw<N>::eval([&p, mod](int i)
      {
        p[i] /= mod; } );
    }
    return mod;
  }
} ;
template <size_t N, typename T>
Vect<N, T> operator +(const Vect_base<N, T> & a, const Vect_base<N, T> & b)
{
  Vect<N, T> x(a);
  return x += b;
}
template <size_t N, typename T>
Vect<N, T> operator -(const Vect_base<N, T> & a, const Vect_base<N, T> & b)
{
  Vect<N, T> x(a);
  return x -= b;
}
template <size_t N, typename T, typename T1>
Vect<N, T> operator*(const Vect_base<N, T> & a, const T1 & b)
{
  Vect<N, T> r;
  iterate_forw<N>::eval([&r, &a, b](size_t i)
  {
    r[i] = a[i] * b;
  } );
  return r;
}
template <size_t N, typename T, typename T1>
Vect<N, T> operator*(const T1 & a, const Vect_base<N, T> & b)
{
  return b * a;
}
template <size_t N, typename T>
T v_dot(const Vect_base<N, T> & a, const Vect_base<N, T> & b)
{
  T v = 0;
  iterate_forw<N>::eval([&a, &b, &v](size_t i)
  {
    v += a[i] * b[i];
  } );
  return v;
}
template <typename T>
Vect<3, T> v_cross(const Vect_base<3, T> & a, const Vect_base<3, T> & b)
{
  Vect<3, T> c;
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
  return c;
}
template <size_t N, typename T> T v_sq_len(const Vect<N, T> & a)
{
  return v_dot(a, a);
}
template <size_t N, typename T> T v_sq_len(const Vers<N, T> & a)
{
  return static_cast<T> (1);
}
template <size_t N, typename T> T v_len(const Vect<N, T> & a)
{
  return sqrt(v_sq_len(a));
}
template <size_t N, typename T> T v_len(const Vers<N, T> & a)
{
  return static_cast<T> (1);
}
template <typename T1, typename T2>
typename T1::mygeotype v_cos_angle(const T1 & a, const T2 & b)
{
  return v_dot(a, b) / sqrt(v_sq_len(a) * v_sq_len(b));
}
template <size_t N, typename T>
T v_sq_cos_angle(const Vect<N, T> & a, const Vect<N, T> & b)
{
  return gk_sq(v_dot(a, b)) / (v_sq_len(a) * v_sq_len(b));
}
template <typename T>
T v_cross(const Vect_base<2, T> & a, const Vect_base<2, T> & b)
{
  return a[0] * b[1] - a[1] * b[0];
}


};

