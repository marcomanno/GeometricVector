
#pragma once

#include "geo_iter.h"
#include "geo_pow.h"

namespace Geo
{
typedef double real;
template <size_t N, typename T> struct Vect_base
{
protected:
  T m_el[N];
  void set_null()
  {
    iterate_forw<N>::eval([this](int i)
    {
      m_el[i] = static_cast<T> (0);
    } );
  }
  bool is_null() const
  {
    return iterate_forw_until<N, false>::eval([this](int i)
    {
      return m_el[i] == static_cast<T> (0);
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
  T len_sq() const
  {
    T v = 0;
    iterate_forw<N>::eval([this, &v](size_t i)
    {
      v += gk_sq<T>(m_el[i]);
    });
    return v;
  }
  T len() const
  {
    return sqrt(len_sq());
  }
  T normalize()
  {
    T mod = Vect_base::len_sq();
    if (mod != 0 && mod != 1.)
    {
      mod = sqrt(mod);
      iterate_forw<N>::eval([this, mod](size_t i)
      {
        m_el[i] /= mod;
      });
    }
    return mod;
  }

public:
  const T& operator [] (size_t i) const
  {
    return m_el[i];
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
  // Make some base class public
  using Vect_base::len;
  using Vect_base::len_sq;
  using Vect_base::normalize;

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
  T & operator [] (size_t i)
  {
    return m_el[i];
  }
private:
  template <typename op> void oper_self(const Vect_base<N, T> & x)
  {
    T * a = this->m_el;
    const T * b = x.coord();
    iterate_forw<N>::eval([&a, &b](size_t i) { op()(a[i], b[i]); } );
  }
};

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
    iterate_forw<N>::eval([this, coo](size_t i)
    {
      m_el[i] = coo[i]
    } );
    return normalize();
  }
  T len_sq() const
  {
    return static_cast<T>(is_null() ? 0. : 1.);
  }
  T len() const
  {
    return len_sq();
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
template <typename T>
T v_cross(const Vect_base<2, T> & a, const Vect_base<2, T> & b)
{
  return a[0] * b[1] - a[1] * b[0];
}
template <typename T1, typename T2>
typename T1::mygeotype v_cos_angle(const T1 & a, const T2 & b)
{
  return v_dot(a, b) / sqrt(a.len_sq() * b.len_sq());
}
template <size_t N, typename T>
T v_sq_cos_angle(const Vect<N, T> & a, const Vect<N, T> & b)
{
  return gk_sq(v_cos_angle(a, b));
}
};

