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
MASA::variabledens<Scalar>::variabledens()
{
  this->mmsname = "variabledens";
  this->dimension = 3;

  this->register_var("R",&gasConst_);
  this->register_var("Pref",&Pref_);
  this->register_var("Tref",&Tref_);
  this->register_var("Cp",&Cp_);
  this->register_var("k",&k);
  this->register_var("u_0",&u_0);
  this->register_var("v_0",&v_0);
  this->register_var("w_0",&w_0);
  this->register_var("p_0",&p_0);
  this->register_var("a",&a_);
  this->register_var("ah",&ah_);
  this->register_var("mu",&mu);
  this->register_var("MolW",&Mval_);

  this->init_var();

} // done with constructor

template <typename Scalar>
int MASA::variabledens<Scalar>::init_var()
{
  int err = 0;

  err += this->set_var("R",10.0);
  err += this->set_var("Pref",100.0);
  err += this->set_var("Tref",300.0);
  err += this->set_var("k",.0015625);
  err += this->set_var("u_0",1.0);
  err += this->set_var("v_0",1.0);
  err += this->set_var("w_0",1.0);
  err += this->set_var("p_0",1.0);
  err += this->set_var("a",20.0);
  err += this->set_var("ah",10.0);
  err += this->set_var("Cp",0.01);
  err += this->set_var("mu",.00125);
  err += this->set_var("MolW",30.0);

  return err;

} // done with init_var

// ----------------------------------------
// Source Terms
// ----------------------------------------

// public static method, that can be called from eval_q_t
template <typename Scalar>
Scalar MASA::variabledens<Scalar>::eval_q_u(Scalar x1, Scalar y1, Scalar z1)
{
  using std::cos;
  using std::sin;

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
  ADScalar P = -p_0/4.0 * (cos(2.0 * a_ * PI * x) + cos(2.0 * a_ * PI * y) + cos(2.0 * a_ * PI * z));

  // Temperature
  ADScalar T = (cos(ah_ * PI * x) * cos(ah_ * PI * y) * cos(ah_ * PI * z))/Cp_ + Tref_;

  // Density
  ADScalar RHO = Pref_ * Mval_ / (gasConst_ * T);

  // The shear strain tensor
  NumberVector<NDIM, typename ADScalar::derivatives_type> GradU = gradient(U);

  // The identity tensor I
  NumberVector<NDIM, NumberVector<NDIM, Scalar> > Identity = 
    NumberVector<NDIM, Scalar>::identity();

  // The shear stress tensor
  NumberVector<NDIM, NumberVector<NDIM, ADScalar> > Tau = mu * (GradU + transpose(GradU) - 2./3.*divergence(U)*Identity);

  // Temperature flux
  NumberVector<NDIM, ADScalar> q = -k * T.derivatives();

  NumberVector<NDIM, Scalar> Q_rho_u = 
    raw_value(divergence(RHO*U.outerproduct(U) - Tau) + P.derivatives());

  return Q_rho_u[0];
}

// public, static method
template <typename Scalar>
Scalar MASA::variabledens<Scalar>::eval_q_v(Scalar x1, Scalar y1, Scalar z1)
{
  using std::cos;
  using std::sin;

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
  ADScalar P = -p_0/4.0 * (cos(2.0 * a_ * PI * x) + cos(2.0 * a_ * PI * y) + cos(2.0 * a_ * PI * z));

  // Temperature
  ADScalar T = (cos(ah_ * PI * x) * cos(ah_ * PI * y) * cos(ah_ * PI * z))/Cp_ + Tref_;

  // Density
  ADScalar RHO = Pref_ * Mval_ / (gasConst_ * T);

  // The shear strain tensor
  NumberVector<NDIM, typename ADScalar::derivatives_type> GradU = gradient(U);

  // The identity tensor I
  NumberVector<NDIM, NumberVector<NDIM, Scalar> > Identity = 
    NumberVector<NDIM, Scalar>::identity();

  // The shear stress tensor
  NumberVector<NDIM, NumberVector<NDIM, ADScalar> > Tau = mu * (GradU + transpose(GradU) - 2./3.*divergence(U)*Identity);

  // Temperature flux
  NumberVector<NDIM, ADScalar> q = -k * T.derivatives();

  NumberVector<NDIM, Scalar> Q_rho_u = 
    raw_value(divergence(RHO*U.outerproduct(U) - Tau) + P.derivatives());

  return Q_rho_u[1];

}

// public, static method
template <typename Scalar>
Scalar MASA::variabledens<Scalar>::eval_q_w(Scalar x1, Scalar y1, Scalar z1)
{
  using std::cos;
  using std::sin;

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
  ADScalar P = -p_0/4.0 * (cos(2.0 * a_ * PI * x) + cos(2.0 * a_ * PI * y) + cos(2.0 * a_ * PI * z));

  // Temperature
  ADScalar T = (cos(ah_ * PI * x) * cos(ah_ * PI * y) * cos(ah_ * PI * z))/Cp_ + Tref_;

  // Density
  ADScalar RHO = Pref_ * Mval_ / (gasConst_ * T);

  // The shear strain tensor
  NumberVector<NDIM, typename ADScalar::derivatives_type> GradU = gradient(U);

  // The identity tensor I
  NumberVector<NDIM, NumberVector<NDIM, Scalar> > Identity = 
    NumberVector<NDIM, Scalar>::identity();

  // The shear stress tensor
  NumberVector<NDIM, NumberVector<NDIM, ADScalar> > Tau = mu * (GradU + transpose(GradU) - 2./3.*divergence(U)*Identity);

  // Temperature flux
  NumberVector<NDIM, ADScalar> q = -k * T.derivatives();

  NumberVector<NDIM, Scalar> Q_rho_u = 
    raw_value(divergence(RHO*U.outerproduct(U) - Tau) + P.derivatives());

  return Q_rho_u[2];

}

// public, static method
template <typename Scalar>
Scalar MASA::variabledens<Scalar>::eval_q_e(Scalar x1, Scalar y1, Scalar z1)
{
  using std::cos;
  using std::sin;

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
  ADScalar P = -p_0/4.0 * (cos(2.0 * a_ * PI * x) + cos(2.0 * a_ * PI * y) + cos(2.0 * a_ * PI * z));

  // Enthalpy
  ADScalar H = cos(ah_ * PI * x) * cos(ah_ * PI * y) * cos(ah_ * PI * z);

  // Temperature
  ADScalar T = H/Cp_ + Tref_;

  // Density
  ADScalar RHO = Pref_ * Mval_ / (gasConst_ * T);

  // The shear strain tensor
  NumberVector<NDIM, typename ADScalar::derivatives_type> GradU = gradient(U);

  // The identity tensor I
  NumberVector<NDIM, NumberVector<NDIM, Scalar> > Identity = 
    NumberVector<NDIM, Scalar>::identity();

  // The shear stress tensor
  NumberVector<NDIM, NumberVector<NDIM, ADScalar> > Tau = mu * (GradU + transpose(GradU) - 2./3.*divergence(U)*Identity);

  // Temperature flux
  NumberVector<NDIM, ADScalar> q = -k * H.derivatives();

  // Euler equation residuals
  // Scalar Q_rho = raw_value(divergence(RHO*U));
  // NumberVector<NDIM, Scalar> Q_rho_u = 
  //   raw_value(divergence(RHO*U.outerproduct(U) - Tau) + P.derivatives());

  // energy equation
  //Scalar Q_rho_e = raw_value(divergence(RHO*H*U + q - Tau.dot(U)));
  Scalar Q_rho_e = raw_value(divergence(RHO*H*U + q));
  
  return Q_rho_e;
}


// public, static method
template <typename Scalar>
Scalar MASA::variabledens<Scalar>::eval_q_rho(Scalar x1, Scalar y1, Scalar z1)
{
  using std::cos;
  using std::sin;

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
  ADScalar P = -p_0/4.0 * (cos(2.0 * a_ * PI * x) + cos(2.0 * a_ * PI * y) + cos(2.0 * a_ * PI * z));

  // Temperature
  ADScalar T = (cos(ah_ * PI * x) * cos(ah_ * PI * y) * cos(ah_ * PI * z))/Cp_ + Tref_;

  // Density
  ADScalar RHO = Pref_ * Mval_ / (gasConst_ * T);

  Scalar Q_rho = raw_value(divergence(RHO*U));

  return Q_rho;
}


// ----------------------------------------
// Analytical Terms
// ----------------------------------------

// example of a public method called from eval_exact_t
template <typename Scalar>
Scalar MASA::variabledens<Scalar>::eval_exact_u(Scalar x, Scalar y, Scalar z)
{
  using std::cos;
  using std::sin;

  Scalar exact_u;
  exact_u = -u_0 * cos(a_ * PI * x) * sin(a_ * PI * y) * sin(a_ * PI * z);
  return exact_u;
}

// public method
template <typename Scalar>
Scalar MASA::variabledens<Scalar>::eval_exact_v(Scalar x, Scalar y, Scalar z)
{
  using std::cos;
  using std::sin;

  Scalar exact_v;
  exact_v = v_0 * sin(a_ * PI * x) * cos(a_ * PI * y) * sin(a_ * PI * z);
  return exact_v;
}

// public method
template <typename Scalar>
Scalar MASA::variabledens<Scalar>::eval_exact_w(Scalar x, Scalar y, Scalar z)
{
  using std::cos;
  using std::sin;

  Scalar exact_w;
  exact_w = -w_0 * sin(a_ * PI * x) * sin(a_ * PI * y) * cos(a_ * PI * z);
  return exact_w;
}

// public method
template <typename Scalar>
Scalar MASA::variabledens<Scalar>::eval_exact_p(Scalar x, Scalar y, Scalar z)
{
  using std::cos;
  using std::sin;

  Scalar P = -p_0/4.0 * (cos(2.0 * a_ * PI * x) + cos(2.0 * a_ * PI * y) + cos(2.0 * a_ * PI * z));
  return P;
}

// public method
template <typename Scalar>
Scalar MASA::variabledens<Scalar>::eval_exact_rho(Scalar x, Scalar y, Scalar z)
{
  using std::cos;
  using std::sin;

  Scalar T = (cos(ah_ * PI * x) * cos(ah_ * PI * y) * cos(ah_ * PI * z))/Cp_ + Tref_;

  Scalar RHO = Pref_ * Mval_ / (gasConst_ * T);
  return RHO;
}

// public method
template <typename Scalar>
Scalar MASA::variabledens<Scalar>::eval_exact_t(Scalar x, Scalar y, Scalar z)
{
  using std::cos;
  using std::sin;

  Scalar T = (cos(ah_ * PI * x) * cos(ah_ * PI * y) * cos(ah_ * PI * z))/Cp_ + Tref_;

  return T;
}


// ----------------------------------------
// Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::variabledens);


#endif // HAVE_METAPHYSICL
