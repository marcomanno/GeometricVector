#define CATCH_CONFIG_MAIN 
#include "catch.hpp"

#include <geo_vector.h>

TEST_CASE("base", "base")
{
  Geo::Vect<4> v0(0.);
  Geo::Vect<4> v1(1.);
  REQUIRE(Geo::v_dot(v0, v1) == 0);
  Geo::Vect<4> v2( 1, 2, 3, 4 );
  Geo::Vect<4> v3{ 1., 2., 3., 4. };
  REQUIRE((v2 - v3).len_sq() == 0);
  Geo::Vers<> a(1, 4, 5);
  REQUIRE(a.len_sq() == 1);
  v0[1] = 3;
  Geo::Vers<4> a1(-1, 4, 5, 3);
  REQUIRE(Geo::v_dot(v0, a1) > 1.6803361007);
  Geo::Vers<3> a2(1, 4, 5);
  REQUIRE(Geo::v_cross(a, a2).len_sq() == 0);
 }

TEST_CASE("pow", "pow")
{
  auto val = Geo::gk_pow<double>::to<6>(4.);
}