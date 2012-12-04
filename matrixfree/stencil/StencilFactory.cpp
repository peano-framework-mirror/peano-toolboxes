#include "matrixfree/stencil/StencilFactory.h"

#include "peano/utils/Loop.h"


matrixfree::stencil::Stencil
matrixfree::stencil::StencilFactory::exchangeCoordinates( const Stencil& stencil, int coord0, int coord1 ) {
  Stencil result;

  d3for3(destination)
    tarch::la::Vector<DIMENSIONS,int> source;
    source = destination;
    source(coord0) = destination(coord1);
    source(coord1) = destination(coord0);
    result( peano::utils::dLinearised(destination,3) ) = stencil( peano::utils::dLinearised(source,3) );
  enddforx

  return result;
}


tarch::la::Vector<5,double> matrixfree::stencil::StencilFactory::get1DLinearInterpolationStencil() {
  tarch::la::Vector<5,double> result;
  result(0) =  1.0 / 3.0;
  result(1) =  2.0 / 3.0;
  result(2) =  3.0 / 3.0;
  result(3) =  2.0 / 3.0;
  result(4) =  1.0 / 3.0;
  return result;
}


tarch::la::Vector<3,double> matrixfree::stencil::StencilFactory::get1DMassStencilWithoutHScaling() {
	tarch::la::Vector<3,double> result;
  result(0) =  1.0 / 6.0;
  result(1) =  4.0 / 6.0;
  result(2) =  1.0 / 6.0;
	return result;
}


tarch::la::Vector<3,double> matrixfree::stencil::StencilFactory::get1DIdentityWithoutHScaling() {
  tarch::la::Vector<3,double> result;
  result(0) =  0.0;
  result(1) =  1.0;
  result(2) =  0.0;
  return result;
}



tarch::la::Vector<3,double> matrixfree::stencil::StencilFactory::get1DMeanValueStencil() {
  tarch::la::Vector<3,double> result;
  result(0) =  1.0 / 2.0;
  result(1) =  0.0;
  result(2) =  1.0 / 2.0;
  return result;
}


tarch::la::Vector<3,double> matrixfree::stencil::StencilFactory::get1DMassStencil(double h) {
  return get1DMassStencilWithoutHScaling() * h;
}


tarch::la::Vector<3,double> matrixfree::stencil::StencilFactory::get1DLaplaceStencilWithoutHScaling() {
  tarch::la::Vector<3,double> result;
  result(0) = -1.0;
  result(1) =  2.0;
  result(2) = -1.0;
  return result;
}

tarch::la::Vector<3,double> matrixfree::stencil::StencilFactory::get1DLaplaceStencil(double h) {
  assertion( h>0.0 );
  return get1DLaplaceStencilWithoutHScaling() / h;
}



tarch::la::Vector<3*3,double> matrixfree::stencil::StencilFactory::stencilProduct(
  const tarch::la::Vector<3,double>& a,
  const tarch::la::Vector<3,double>& b
) {
  tarch::la::Vector<3*3,double> result;

  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      result(i+j*3) = a(i) * b(j);
    }
  }

  return result;
}


tarch::la::Vector<5*5,double> matrixfree::stencil::StencilFactory::stencilProduct(
  const tarch::la::Vector<5,double>& a,
  const tarch::la::Vector<5,double>& b
) {
  tarch::la::Vector<5*5,double> result;

  for (int i=0; i<5; i++) {
    for (int j=0; j<5; j++) {
      result(i+j*5) = a(i) * b(j);
    }
  }

  return result;
}


tarch::la::Vector<3*3*3,double> matrixfree::stencil::StencilFactory::stencilProduct(
  const tarch::la::Vector<3,double>& a,
  const tarch::la::Vector<3,double>& b,
  const tarch::la::Vector<3,double>& c
) {
  tarch::la::Vector<3*3*3,double> result;

  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
        result(i+j*3+k*3*3) = a(i) * b(j) * c(k);
      }
    }
  }

  return result;
}


tarch::la::Vector<5*5*5,double> matrixfree::stencil::StencilFactory::stencilProduct(
  const tarch::la::Vector<5,double>& a,
  const tarch::la::Vector<5,double>& b,
  const tarch::la::Vector<5,double>& c
) {
  tarch::la::Vector<5*5*5,double> result;

  for (int i=0; i<5; i++) {
    for (int j=0; j<5; j++) {
      for (int k=0; k<5; k++) {
        result(i+j*5+k*5*5) = a(i) * b(j) * c(k);
      }
    }
  }

  return result;
}


tarch::la::Vector<THREE_POWER_D,double> matrixfree::stencil::StencilFactory::getLaplacian(double scaling, const tarch::la::Vector<DIMENSIONS,double>& h) {
  tarch::la::Vector<THREE_POWER_D,double> result;

  #if defined(Dim2)
  result =
    matrixfree::stencil::StencilFactory::stencilProduct(
      matrixfree::stencil::StencilFactory::get1DLaplaceStencil(h(0)),
      matrixfree::stencil::StencilFactory::get1DMassStencil(h(1))
    ) +
    matrixfree::stencil::StencilFactory::stencilProduct(
      matrixfree::stencil::StencilFactory::get1DMassStencil(h(0)),
      matrixfree::stencil::StencilFactory::get1DLaplaceStencil(h(1))
    );
  #elif defined(Dim3)
  result =
    matrixfree::stencil::StencilFactory::stencilProduct(
      matrixfree::stencil::StencilFactory::get1DLaplaceStencil(h(0)),
      matrixfree::stencil::StencilFactory::get1DMassStencil(h(1)),
      matrixfree::stencil::StencilFactory::get1DMassStencil(h(2))
    ) +
    matrixfree::stencil::StencilFactory::stencilProduct(
      matrixfree::stencil::StencilFactory::get1DMassStencil(h(0)),
      matrixfree::stencil::StencilFactory::get1DLaplaceStencil(h(1)),
      matrixfree::stencil::StencilFactory::get1DMassStencil(h(2))
    ) +
    matrixfree::stencil::StencilFactory::stencilProduct(
      matrixfree::stencil::StencilFactory::get1DMassStencil(h(0)),
      matrixfree::stencil::StencilFactory::get1DMassStencil(h(2)),
      matrixfree::stencil::StencilFactory::get1DLaplaceStencil(h(2))
    );
  #else
  assertionMsg( false, "dimension not supported" );
  #endif

  result *= scaling;

  return result;
}


tarch::la::Vector<THREE_POWER_D,double> matrixfree::stencil::StencilFactory::getLaplacian(const tarch::la::Vector<DIMENSIONS,double>& scaling, const tarch::la::Vector<DIMENSIONS,double>& h) {
  tarch::la::Vector<THREE_POWER_D,double> result;

  #if defined(Dim2)
  result =
    scaling(0) *
    matrixfree::stencil::StencilFactory::stencilProduct(
      matrixfree::stencil::StencilFactory::get1DLaplaceStencil(h(0)),
      matrixfree::stencil::StencilFactory::get1DMassStencil(h(1))
    ) +
    scaling(1) *
    matrixfree::stencil::StencilFactory::stencilProduct(
      matrixfree::stencil::StencilFactory::get1DMassStencil(h(0)),
      matrixfree::stencil::StencilFactory::get1DLaplaceStencil(h(1))
    );
  #elif defined(Dim3)
  result =
    scaling(0) *
    matrixfree::stencil::StencilFactory::stencilProduct(
      matrixfree::stencil::StencilFactory::get1DLaplaceStencil(h(0)),
      matrixfree::stencil::StencilFactory::get1DMassStencil(h(1)),
      matrixfree::stencil::StencilFactory::get1DMassStencil(h(2))
    ) +
    scaling(1) *
    matrixfree::stencil::StencilFactory::stencilProduct(
      matrixfree::stencil::StencilFactory::get1DMassStencil(h(0)),
      matrixfree::stencil::StencilFactory::get1DLaplaceStencil(h(1)),
      matrixfree::stencil::StencilFactory::get1DMassStencil(h(2))
    ) +
    scaling(2) *
    matrixfree::stencil::StencilFactory::stencilProduct(
      matrixfree::stencil::StencilFactory::get1DMassStencil(h(0)),
      matrixfree::stencil::StencilFactory::get1DMassStencil(h(2)),
      matrixfree::stencil::StencilFactory::get1DLaplaceStencil(h(2))
    );
  #else
  assertionMsg( false, "dimension not supported" );
  #endif

  return result;
}


tarch::la::Vector<THREE_POWER_D,double> matrixfree::stencil::StencilFactory::getMassMatrix(const tarch::la::Vector<DIMENSIONS,double>& h) {
  tarch::la::Vector<THREE_POWER_D,double> result;

  #if defined(Dim2)
  result =
    matrixfree::stencil::StencilFactory::stencilProduct(
      matrixfree::stencil::StencilFactory::get1DMassStencil(h(0)),
      matrixfree::stencil::StencilFactory::get1DMassStencil(h(1))
    );
  #elif defined(Dim3)
  result =
    matrixfree::stencil::StencilFactory::stencilProduct(
      matrixfree::stencil::StencilFactory::get1DMassStencil(h(0)),
      matrixfree::stencil::StencilFactory::get1DMassStencil(h(1)),
      matrixfree::stencil::StencilFactory::get1DMassStencil(h(2))
    );
  #else
  assertionMsg( false, "dimension not supported" );
  #endif

  return result;
}


tarch::la::Vector<THREE_POWER_D,double> matrixfree::stencil::StencilFactory::getIdentity(const tarch::la::Vector<DIMENSIONS,double>& h) {
  tarch::la::Vector<THREE_POWER_D,double> result;

  #if defined(Dim2)
  result =
    matrixfree::stencil::StencilFactory::stencilProduct(
      matrixfree::stencil::StencilFactory::get1DIdentityWithoutHScaling(),
      matrixfree::stencil::StencilFactory::get1DIdentityWithoutHScaling()
    )
    * h(0) * h(1);
  #elif defined(Dim3)
  result =
    matrixfree::stencil::StencilFactory::stencilProduct(
      matrixfree::stencil::StencilFactory::get1DIdentityWithoutHScaling(),
      matrixfree::stencil::StencilFactory::get1DIdentityWithoutHScaling(),
      matrixfree::stencil::StencilFactory::get1DIdentityWithoutHScaling()
    )
    * h(0) * h(1) * h(2);
  #else
  assertionMsg( false, "dimension not supported" );
  #endif

  return result;
}


tarch::la::Vector<FIVE_POWER_D,double> matrixfree::stencil::StencilFactory::getDLinearInterpolation() {
  tarch::la::Vector<FIVE_POWER_D,double> result;
  #ifdef Dim2
  result =
    matrixfree::stencil::StencilFactory::stencilProduct(
      matrixfree::stencil::StencilFactory::get1DLinearInterpolationStencil(),
      matrixfree::stencil::StencilFactory::get1DLinearInterpolationStencil()
    );
  #elif Dim3
  result =
    matrixfree::stencil::StencilFactory::stencilProduct(
      matrixfree::stencil::StencilFactory::get1DLinearInterpolationStencil(),
      matrixfree::stencil::StencilFactory::get1DLinearInterpolationStencil(),
      matrixfree::stencil::StencilFactory::get1DLinearInterpolationStencil()
    );
  #else
  assertionMsg( false, "dimension not supported" );
  #endif

  return result;
}
