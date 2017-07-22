#pragma once

namespace Geo
{
// Class interval
template <typename T = real> class Interval
{
  T m_limits[2];
  static const T min_max(bool at_end)
  {
    return at_end ?
      std::numeric_limits<T>::min() :
      std::numeric_limits<T>::max();
  }
  void check()
  {
    if (m_limits[1] < m_limits[0] && !is_empty())
      * this = Empty();
  }
public:
  typedef T mygeotype;
  Interval()
  {
    *this = Infinite();
  }
  Interval(const T vals[2])
  {
    m_limits[0] = vals[0];
    m_limits[1] = vals[1];
    check();
  }
  Interval(const T & start, const T & end)
  {
    m_limits[0] = start;
    m_limits[1] = end;
    check();
  }
  static Interval Empty()
  {
    return Interval(min_max(false), min_max(true));
  }
  static Interval Infinite()
  {
    return Interval(min_max);
  }
  bool operator == (const Interval oth) const
  {
    return m_limits[0] == oth[0] && m_limits[1] == oth[1];
  }
  bool is_empty() const
  {
    return *this == Empty();
  }
  bool is_infinite() const
  {
    return *this == Infinite();
  }
  const T & operator[](int i) const
  {
    return m_limits[i & 1];
  }
  void set(bool up_lim, const T & val)
  {
    m_limits[up_lim] = val;
    check();
  }
  T length()
  {
    return m_limits[1] - m_limits[0];
  }
  bool contain(const T & val, const T tol = 0) const
  {
    return !(val < m_limits[0] - tol) &&
      !(m_limits[1] + tol < val);
  }
  bool contain(const Interval & in, const T tol) const
  {
    return !(in[0] < m_limits[0] - tol) &&
      !(m_limits[1] + tol < in[1]);
  }
  bool disjoint(const Interval & in, const T tol) const
  {
    return (in[1] < m_limits[0] - tol) ||
      (m_limits[1] + tol < in[0]);
  }
  bool add(const T & el)
  {
    if (el < m_limits[0])
    {
      m_limits[0] = el;
      return true;
    }
    if (m_limits[1] < el)
    {
      m_limits[1] = el;
      return true;
    }
    return false;
  }
  bool add(const Interval & in)
  {
    bool achange = in[0] < m_limits[0];
    if (achange)
      m_limits[0] = in[0];
    if (m_limits[1] < in[1])
    {
      m_limits[1] = in[1];
      return true;
    }
    else
      return achange;
  }
  bool subtract(const Interval & in)
  {
    if (m_limits[0] < in[0])
    {
      m_limits[0] = in[0];
      return true;
    }
    if (in[1] < m_limits[1])
    {
      m_limits[1] = in[1];
      return true;
    }
    return false;
  }
  void grow(const T & delta)
  {
    if (m_limits[0] != min_max(false))
      m_limits[0] -= delta;
    if (m_limits[1] != min_max(true))
      m_limits[1] += delta;
  }
  void shift(const T & delta)
  {
    if (m_limits[0] != min_max(false))
      m_limits[0] += delta;
    if (m_limits[1] != min_max(true))
      m_limits[1] += delta;
  }
  Interval intersect(Interval & in, const Interval & in2) const
  {
    Interval res(*this);
    res.subtract(in);
    return res;
  }
  Interval unite(const Interval & in) const
  {
    Interval res(*this);
    res.add(in);
    return res;
  }

};
template <typename T = real> struct Range : public Interval<T>
{
  enum TY
  {
    OPEN, CLOSED, PERIODIC
  };
  Range(TY ty) : Interval<T>(), m_ty(ty)
  {
  }
  Range(const T vals[2], TY ty) : Interval<T>(vals), m_ty(ty)
  {
  }
  Range(const T & start, const T & end, TY ty) : Interval<T>(start, end), m_ty(ty)
  {
  }
  bool contain(Interval<T> & in, const T & val, const T tol = 0) const
  {
    if (whole_interval(in, tol))
      return true;
    if (in.contain(val, tol))
      return true;
    if (m_ty != PERIODIC)
      return false;
    T sh = find_peridic_shift(in, val);
    val += sh;
    return in.contain(val, tol);
  }
  bool contain(const Interval<T> & in, const Interval<T> & in2, const T tol) const
  {
    if (whole_interval(in, tol))
      return true;
    if (in.contain(in2, tol))
      return true;
    if (m_ty != PERIODIC)
      return false;
    Interval<T> in3 = align_periodic_intervals(in, in2, tol);
    return in.contain(in3, tol);
  }
  bool disjoint(const Interval<T> & in, const Interval<T> & in2, const T tol) const
  {
    if (whole_interval(in, tol) || whole_interval(in2, tol))
      return false;
    if (!in.disjoint(in2, tol))
      return false;
    if (m_ty != PERIODIC)
      return true;
    Interval<T> in3 = align_periodic_intervals(in, in2, tol);
    return in.disjoint(in3, tol);
  }
  bool add(Interval<T> & in, const T & el) const
  {
    if (m_ty == OPEN)
      in.add(el);
    else if (m_ty == CLOSED)
    {
      if (!this->contain(el, 0))
        return false;
      in.add(el);
    }
    else // if ( m_ty == PERIODIC )
    {
      T sh = find_peridic_shift(in, el, 0);
      if (sh == 0) return false;
      T val[2];
      if (sh < 0)
      {
        val[0] = sh;
        val[1] = sh + this->length();
      }
      else
      {
        val[1] = sh;
        val[0] = sh - this->length();
      }
      if (val[1] - in[1] < in[0] - val[0])
        return in.add(val[1]);
      else
        return in.add(val[0]);
    }
  }
  bool add(Interval<T> & in, const Interval<T> & in2) const
  {
    if (m_ty != PERIODIC)
      return in.add(in2);
    Interval<T> in3 = align_periodic_intervals(in, in2, 0, true);
    return in.add(in3);
  }
  bool subtract(Interval<T> & in, const Interval<T> & in2) const
  {
    if (m_ty != PERIODIC)
      return in.subtract(in2);
    Interval<T> in3 = align_periodic_intervals(in, in2, 0, true);
    return in.subtract(in3);
  }
  void grow(Interval<T> & in, const T & delta) const
  {
    in.grow(delta);
    check(in);
  }
  void shift(Interval<T> & in, const T & delta) const
  {
    in.shift(delta);
    check(in);
  }
  Interval<T> intersect(const Interval<T> & in, const Interval<T> & in2) const
  {
    Interval<T> res(in);
    subtract(res, in2);
    return res;
  }
  Interval<T> unite(Interval<T> & in, const Interval<T> & in2) const
  {
    Interval<T> res(in);
    add(in, in2);
    return res;
  }

private:
  TY m_ty;
  T find_peridic_shift(Interval<T> & in, const T & val, const T & tol)
  {
    T diff;
    bool chsign;
    if ((diff = in[0] - tol - val) > 0)      chsign = false;
    else if ((diff = val - in[1] + tol) > 0) chsign = true;
    else return false;
    T ratio = diff / this->length();
    T n = floor(ratio);
    if (ratio - n > 0) n += 1;
    if (chsign) n = -n;
    return n * this->length();
  }
  Interval<T> align_periodic_intervals(
    const Interval<T> & in, const Interval<T> & in2,
    T tol, bool near = false) const
  {
    T shift;
    if (in2[near] < in[0] - tol)
      shift = find_peridic_shift(in, in2[near], tol);
    else if (in[1] + tol < in2[!near])
      shift = find_peridic_shift(in, in2[!near], tol);
    else
      return in2;
    Interval<T> out_in(in2);
    out_in.shift(shift);
    return out_in;
  }
  bool whole_interval(const Interval<T> & in, T tol) const
  {
    if (in == Interval<T>::Infinite())
      return true;
    if (m_ty == OPEN)
      return false;
    if (!(in.length() < this->length() - 2 * tol))
      return true; // Whole interval.
    return false;
  }
  void check(Interval<T> &in) const
  {
    if (m_ty == CLOSED)
    {
      if (in[0] < this->m_limits[0])
        in.set(false, this->m_limits[0]);
      if (this->m_limits[1] < in[1])
        in.set(true, this->m_limits[1]);
    }
    else if (m_ty == PERIODIC)
    {
      if (in.length() > this->length())
        in = *this;
    }
  }
};
template <typename Interval_ty, size_t N> class Par_box
{
  Interval_ty m_intervals[N];
public:
  const Interval_ty & operator[](size_t i) const
  {
    return m_intervals[i];
  }
  void set(size_t i, const Interval_ty & new_i_int)
  {
    m_intervals[i] = new_i_int;
  }
  void set(size_t i, bool up_lim, const typename Interval_ty::mygeotype & val)
  {
    m_intervals[i].set(up_lim, val);
  }
  bool is_empty() const
  {
    Interval_ty * intervs = m_intervals;
    return iterate_forw_until<N, false>::eval([&intervs](int i)
    {
      return intervs[i].is_empty();
    });
  }
  bool is_infinite() const
  {
    Interval_ty * intervs = m_intervals;
    return !is_empty() &&
      iterate_forw_until<N, false>::eval([&intervs](int i)
    {
      return intervs[i].is_infinite();
    });
  }
  bool contain(const Interval_ty & in, const typename Interval_ty::mygeotype tol) const
  {
    Interval_ty * intervs = m_intervals;
    return iterate_forw_until<N>::eval([&intervs, &in, &tol](int i)
    {
      return intervs[i].contain(in[i], tol);
    }
    );
  }
  bool disjoint(const Par_box<Interval_ty, N> & in, const typename Interval_ty::mygeotype tol)const
  {
    Interval_ty * intervs = m_intervals;
    return iterate_forw_until<N, false>::eval([&intervs, &in, &tol](int i)
    {
      return intervs[i].disjoint(in[i], tol);
    }
    );
  }
};


enum class geo_type;

// Base class of objects with reference counting.
class Ref_count
{
  friend void intrusive_ptr_add_ref(Ref_count * p);
  friend void intrusive_ptr_release(Ref_count * p);
  size_t  m_count;
protected:
  Ref_count() : m_count(0)
  {
  };
  virtual ~Ref_count()
  {
  }
public:
  template<class T> T * get_interface()
  {
    return dynamic_cast<T *> (this);
  }
  template<class T> const T * get_interface() const
  {
    return dynamic_cast<const T *> (this);
  }
};
inline void intrusive_ptr_add_ref(Ref_count * p)
{
  ++p->m_count;
}
inline void intrusive_ptr_release(Ref_count * p)
{
  if (--p->m_count == 0u)
    delete p;
}

//template <class Obj> using smart_ptr = boost::intrusive_ptr<Obj>;
template<typename In, typename Out>
class gk_func : public Ref_count
{
public:
  virtual bool ev(const In &, Out &) = 0;
  virtual bool ev_der(const In &, const char nder[In::Dimension], Out &) = 0;
};

template <class Out> using gk_curve_fun = gk_func<real, Out>;
template <typename Out> class gk_surface : public gk_func<Vect<2>, Out>
{
  // Compute the determinant:
  // coord1  coord2  coord3  .. coordM-1
  // Der_coord1
  // Der_coord2
  // ...
  // Der_coordM-1
  virtual bool ev_norm(const Vect<2> &, Out &) = 0;
};
}