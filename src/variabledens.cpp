// -*-c++-*-
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011,2012,2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#include <masa_internal.h>

#ifdef HAVE_METAPHYSICL

#include <ad_masa.h>

typedef ShadowNumber<double, long double> RawScalar;
const unsigned int NDIM = 3;
typedef DualNumber<RawScalar, NumberVector<NDIM, RawScalar> > FirstDerivType;
typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
typedef SecondDerivType ADType;

using namespace MASA;

template <typename Scalar>
MASA::ad_cns_3d_crossterms<Scalar>::ad_cns_3d_crossterms()
{
  this->mmsname = "ad_cns_3d_crossterms";
  this->dimension = 3;

  this->register_var("R",&R);
  this->register_var("k",&k);
  this->register_var("u_0",&u_0);
  this->register_var("u_x",&u_x);
  this->register_var("u_y",&u_y);
  this->register_var("u_z",&u_z);
  this->register_var("v_0",&v_0);
  this->register_var("v_x",&v_x);
  this->register_var("v_y",&v_y);
  this->register_var("v_z",&v_z);
  this->register_var("w_0",&w_0);
  this->register_var("w_x",&w_x);
  this->register_var("w_y",&w_y);
  this->register_var("w_z",&w_z);
  this->register_var("rho_0",&rho_0);
  this->register_var("rho_x",&rho_x);
  this->register_var("rho_y",&rho_y);
  this->register_var("rho_z",&rho_z);
  this->register_var("p_0",&p_0);
  this->register_var("p_x",&p_x);
  this->register_var("p_y",&p_y);
  this->register_var("p_z",&p_z);
  this->register_var("a_px",&a_px);
  this->register_var("a_py",&a_py);
  this->register_var("a_pz",&a_pz);
  this->register_var("a_rhox",&a_rhox);
  this->register_var("a_rhoy",&a_rhoy);
  this->register_var("a_rhoz",&a_rhoz);
  this->register_var("a_ux",&a_ux);
  this->register_var("a_uy",&a_uy);
  this->register_var("a_uz",&a_uz);
  this->register_var("a_vx",&a_vx);
  this->register_var("a_vy",&a_vy);
  this->register_var("a_vz",&a_vz);
  this->register_var("a_wx",&a_wx);
  this->register_var("a_wy",&a_wy);
  this->register_var("a_wz",&a_wz);
  this->register_var("Gamma",&Gamma);
  this->register_var("mu",&mu);
  this->register_var("L",&L);

  this->init_var();

} // done with constructor

template <typename Scalar>
int MASA::ad_cns_3d_crossterms<Scalar>::init_var()
{
  int err = 0;

  err += this->set_var("R",1.01);
  err += this->set_var("k",1.38);
  err += this->set_var("u_0",1.0);
  err += this->set_var("v_0",1.0);
  err += this->set_var("w_0",1.0);
  err += this->set_var("rho_0",10.02);
  err += this->set_var("p_0",1.0);
  err += this->set_var("a_",20.0);
  err += this->set_var("ah_",10.0);
  err += this->set_var("Cp_",0.01);
  err += this->set_var("Tref_",300.0);
  err += this->set_var("Gamma",1.01);
  err += this->set_var("mu",.00125);
  err += this->set_var("L",3.02);

  return err;

} // done with init_var

// ----------------------------------------
// Source Terms
// ----------------------------------------

// public static method, that can be called from eval_q_t
template <typename Scalar>
Scalar MASA::ad_cns_3d_crossterms<Scalar>::eval_q_u(Scalar x1, Scalar y1, Scalar z1) const
{
  using std::cos;

  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef SecondDerivType ADScalar;

  const ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  const ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  const ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U;

  // Arbitrary manufactured solution
  U[0] = -u_0 * cos(a_ * PI * x) * sin(a_ * PI * y) * sin(a_ * PI * z);
  U[1] =  v_0 * sin(a_ * PI * x) * cos(a_ * PI * y) * sin(a_ * PI * z);
  U[2] = -w_0 * sin(a_ * PI * x) * sin(a_ * PI * y) * cos(a_ * PI * z);
  ADScalar RHO = 1.0; 
  ADScalar P = -p_0/4.0 * ((cos(2.0 * a_ * PI * x) + cos(2.0 * a_ * PI * y) + cos(2.0 * a_ * PI * z));

  // Temperature
  ADScalar T = (cos(ah_ * PI * x) * cos(ah_ * PI * z) * cos(ah_ * PI * z))/Cp_ + Tref_;

  // Perfect gas energies
  ADScalar E = 1./(Gamma-1.)*P/RHO;
  ADScalar ET = E + .5 * U.dot(U);

  // The shear strain tensor
  NumberVector<NDIM, typename ADScalar::derivatives_type> GradU = gradient(U);

  // The identity tensor I
  NumberVector<NDIM, NumberVector<NDIM, Scalar> > Identity = 
    NumberVector<NDIM, Scalar>::identity();

  // The shear stress tensor
  NumberVector<NDIM, NumberVector<NDIM, ADScalar> > Tau = mu * (GradU + transpose(GradU) - 2./3.*divergence(U)*Identity);

  // Temperature flux
  NumberVector<NDIM, ADScalar> q = -k * T.derivatives();

  // Euler equation residuals
  // Scalar Q_rho = raw_value(divergence(RHO*U));
  NumberVector<NDIM, Scalar> Q_rho_u = 
    raw_value(divergence(RHO*U.outerproduct(U) - Tau) + P.derivatives());

  return Q_rho_u[0];
}

// public, static method
template <typename Scalar>
Scalar MASA::ad_cns_3d_crossterms<Scalar>::eval_q_v(Scalar x1, Scalar y1, Scalar z1) const
{
  using std::cos;

  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef SecondDerivType ADScalar;

  const ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  const ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  const ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U;

  // Arbitrary manufactured solution
  U[0] = u_0 + u_x * cos(a_ux * PI * x / L) * u_y * cos(a_uy * PI * y / L) * cos(a_uy * PI * z / L);
  U[1] = v_0 + v_x * cos(a_vx * PI * x / L) * v_y * cos(a_vy * PI * y / L) * cos(a_vy * PI * z / L);
  U[2] = w_0 + w_x * cos(a_wx * PI * x / L) * w_y * cos(a_wy * PI * y / L) * cos(a_wy * PI * z / L);
  ADScalar RHO = rho_0 + rho_x * cos(a_rhox * PI * x / L) * rho_y * cos(a_rhoy * PI * y / L) * cos(a_rhoz * PI * z / L);
  ADScalar P = p_0 + p_x * cos(a_px * PI * x / L) * p_y * cos(a_py * PI * y / L) * cos(a_pz * PI * z / L);

  // Temperature
  ADScalar T = P / RHO / R;

 // Perfect gas energies
  ADScalar E = 1./(Gamma-1.)*P/RHO;
  ADScalar ET = E + .5 * U.dot(U);

  // The shear strain tensor
  NumberVector<NDIM, typename ADScalar::derivatives_type> GradU = gradient(U);

  // The identity tensor I
  NumberVector<NDIM, NumberVector<NDIM, Scalar> > Identity = 
    NumberVector<NDIM, Scalar>::identity();

  // The shear stress tensor
  NumberVector<NDIM, NumberVector<NDIM, ADScalar> > Tau = mu * (GradU + transpose(GradU) - 2./3.*divergence(U)*Identity);

  // Temperature flux
  NumberVector<NDIM, ADScalar> q = -k * T.derivatives();

  // Euler equation residuals
  // Scalar Q_rho = raw_value(divergence(RHO*U));
  NumberVector<NDIM, Scalar> Q_rho_u = 
    raw_value(divergence(RHO*U.outerproduct(U) - Tau) + P.derivatives());

  return Q_rho_u[1];

}

// public, static method
template <typename Scalar>
Scalar MASA::ad_cns_3d_crossterms<Scalar>::eval_q_w(Scalar x1, Scalar y1, Scalar z1) const
{
  using std::cos;

  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef SecondDerivType ADScalar;

  const ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  const ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  const ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U;

  // Arbitrary manufactured solution
  U[0] = u_0 + u_x * cos(a_ux * PI * x / L) * u_y * cos(a_uy * PI * y / L) * cos(a_uy * PI * z / L);
  U[1] = v_0 + v_x * cos(a_vx * PI * x / L) * v_y * cos(a_vy * PI * y / L) * cos(a_vy * PI * z / L);
  U[2] = w_0 + w_x * cos(a_wx * PI * x / L) * w_y * cos(a_wy * PI * y / L) * cos(a_wy * PI * z / L);
  ADScalar RHO = rho_0 + rho_x * cos(a_rhox * PI * x / L) * rho_y * cos(a_rhoy * PI * y / L) * cos(a_rhoz * PI * z / L);
  ADScalar P = p_0 + p_x * cos(a_px * PI * x / L) * p_y * cos(a_py * PI * y / L) * cos(a_pz * PI * z / L);

  // Temperature
  ADScalar T = P / RHO / R;

 // Perfect gas energies
  ADScalar E = 1./(Gamma-1.)*P/RHO;
  ADScalar ET = E + .5 * U.dot(U);

  // The shear strain tensor
  NumberVector<NDIM, typename ADScalar::derivatives_type> GradU = gradient(U);

  // The identity tensor I
  NumberVector<NDIM, NumberVector<NDIM, Scalar> > Identity = 
    NumberVector<NDIM, Scalar>::identity();

  // The shear stress tensor
  NumberVector<NDIM, NumberVector<NDIM, ADScalar> > Tau = mu * (GradU + transpose(GradU) - 2./3.*divergence(U)*Identity);

  // Temperature flux
  NumberVector<NDIM, ADScalar> q = -k * T.derivatives();

  // Euler equation residuals
  // Scalar Q_rho = raw_value(divergence(RHO*U));
  NumberVector<NDIM, Scalar> Q_rho_u = 
    raw_value(divergence(RHO*U.outerproduct(U) - Tau) + P.derivatives());

  return Q_rho_u[2];

}

// public, static method
template <typename Scalar>
Scalar MASA::ad_cns_3d_crossterms<Scalar>::eval_q_e(Scalar x1, Scalar y1, Scalar z1) const
{
  using std::cos;

  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef SecondDerivType ADScalar;

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U;

  const ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  const ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  const ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  // Arbitrary manufactured solution
  U[0] = u_0 + u_x * cos(a_ux * PI * x / L) * u_y * cos(a_uy * PI * y / L) * cos(a_uy * PI * z / L);
  U[1] = v_0 + v_x * cos(a_vx * PI * x / L) * v_y * cos(a_vy * PI * y / L) * cos(a_vy * PI * z / L);
  U[2] = w_0 + w_x * cos(a_wx * PI * x / L) * w_y * cos(a_wy * PI * y / L) * cos(a_wy * PI * z / L);
  ADScalar RHO = rho_0 + rho_x * cos(a_rhox * PI * x / L) * rho_y * cos(a_rhoy * PI * y / L) * cos(a_rhoz * PI * z / L);
  ADScalar P = p_0 + p_x * cos(a_px * PI * x / L) * p_y * cos(a_py * PI * y / L) * cos(a_pz * PI * z / L);

  // Temperature
  ADScalar T = P / RHO / R;

  // Perfect gas energies
  ADScalar E = 1./(Gamma-1.)*P/RHO;
  ADScalar ET = E + .5 * U.dot(U);

  // The shear strain tensor
  NumberVector<NDIM, typename ADScalar::derivatives_type> GradU = gradient(U);

  // The identity tensor I
  NumberVector<NDIM, NumberVector<NDIM, Scalar> > Identity = 
    NumberVector<NDIM, Scalar>::identity();

  // The shear stress tensor
  NumberVector<NDIM, NumberVector<NDIM, ADScalar> > Tau = mu * (GradU + transpose(GradU) - 2./3.*divergence(U)*Identity);

  // Temperature flux
  NumberVector<NDIM, ADScalar> q = -k * T.derivatives();

  // Euler equation residuals
  // Scalar Q_rho = raw_value(divergence(RHO*U));
  // NumberVector<NDIM, Scalar> Q_rho_u = 
  //   raw_value(divergence(RHO*U.outerproduct(U) - Tau) + P.derivatives());

  // energy equation
  Scalar Q_rho_e = raw_value(divergence((RHO*ET+P)*U + q - Tau.dot(U)));

  return Q_rho_e;
}


// public, static method
template <typename Scalar>
Scalar MASA::ad_cns_3d_crossterms<Scalar>::eval_q_rho(Scalar x1, Scalar y1, Scalar z1) const
{
  using std::cos;

  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef SecondDerivType ADScalar;

  const ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  const ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  const ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U;

  // Arbitrary manufactured solution
  U[0] = u_0 + u_x * cos(a_ux * PI * x / L) * u_y * cos(a_uy * PI * y / L) * cos(a_uy * PI * z / L);
  U[1] = v_0 + v_x * cos(a_vx * PI * x / L) * v_y * cos(a_vy * PI * y / L) * cos(a_vy * PI * z / L);
  U[2] = w_0 + w_x * cos(a_wx * PI * x / L) * w_y * cos(a_wy * PI * y / L) * cos(a_wy * PI * z / L);
  ADScalar RHO = rho_0 + rho_x * cos(a_rhox * PI * x / L) * rho_y * cos(a_rhoy * PI * y / L) * cos(a_rhoz * PI * z / L);
  ADScalar P = p_0 + p_x * cos(a_px * PI * x / L) * p_y * cos(a_py * PI * y / L) * cos(a_pz * PI * z / L);

  Scalar Q_rho = raw_value(divergence(RHO*U));

  return Q_rho;
}


// ----------------------------------------
// Analytical Terms
// ----------------------------------------

// example of a public method called from eval_exact_t
template <typename Scalar>
Scalar MASA::ad_cns_3d_crossterms<Scalar>::eval_exact_u(Scalar x, Scalar y, Scalar z)
{
  using std::cos;

  Scalar exact_u;
  exact_u = u_0 + u_x * cos(a_ux * PI * x / L) * u_y * cos(a_uy * PI * y / L) * cos(a_uz * PI * z / L);
  return exact_u;
}

// public method
template <typename Scalar>
Scalar MASA::ad_cns_3d_crossterms<Scalar>::eval_exact_v(Scalar x, Scalar y, Scalar z)
{
  using std::cos;

  Scalar exact_v;
  exact_v = v_0 + v_x * cos(a_vx * PI * x / L) * v_y * cos(a_vy * PI * y / L) * cos(a_vz * PI * z / L);
  return exact_v;
}

// public method
template <typename Scalar>
Scalar MASA::ad_cns_3d_crossterms<Scalar>::eval_exact_w(Scalar x, Scalar y, Scalar z)
{
  using std::cos;

  Scalar exact_v;
  exact_v = v_0 + v_x * cos(a_vx * PI * x / L) * v_y * cos(a_vy * PI * y / L) * cos(a_wz * PI * z / L);
  return exact_v;
}

// public method
template <typename Scalar>
Scalar MASA::ad_cns_3d_crossterms<Scalar>::eval_exact_p(Scalar x, Scalar y, Scalar z)
{
  using std::cos;

  Scalar P = p_0 + p_x * cos(a_px * PI * x / L) * p_y * cos(a_py * PI * y / L) * cos(a_pz * PI * z / L);
  return P;
}

// public method
template <typename Scalar>
Scalar MASA::ad_cns_3d_crossterms<Scalar>::eval_exact_rho(Scalar x, Scalar y, Scalar z)
{
  using std::cos;

  Scalar RHO = rho_0 + rho_x * cos(a_rhox * PI * x / L) * rho_y * cos(a_rhoy * PI * y / L) * cos(a_rhoz * PI * z / L);
  return RHO;
}



// ----------------------------------------
// Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::ad_cns_3d_crossterms);


#endif // HAVE_METAPHYSICL
