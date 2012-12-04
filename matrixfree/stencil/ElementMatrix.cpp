#include "matrixfree/stencil/ElementMatrix.h"

#include "peano/utils/Loop.h"

matrixfree::stencil::ElementMatrix::ElementMatrix() {
  _centralElement = 0;
  int base = 1;
  for (int d=0; d<DIMENSIONS; d++ ) {
    _centralElement += base;
    base            *= 3;
  }
}


matrixfree::stencil::ElementWiseAssemblyMatrix
matrixfree::stencil::ElementMatrix::getElementWiseAssemblyMatrix( const matrixfree::stencil::Stencil& stencil ) const {
  /**
   * @todo Die Abbildung sollte man im Konstruktor einmal bauen und dann hier
   * nur noch anwenden. Deshalb ist das Ding ja eine Methode und kein
   * statisches Ding.
   */
  ElementWiseAssemblyMatrix result;

  dfor2(j)
  dfor2(i)
    tarch::la::Vector<DIMENSIONS,int> stencilEntry;
    double    commonFacesPowerTwo = 1.0;
    for (int d=0; d<DIMENSIONS; d++) {
      stencilEntry(d) = i(d)-j(d)+1;
      if (i(d)==j(d)) commonFacesPowerTwo *= 2.0;
    }
    result(jScalar,iScalar) = stencil(peano::utils::dLinearised(stencilEntry,3)) / commonFacesPowerTwo;
  enddforx
  enddforx

  return result;
}


double matrixfree::stencil::ElementMatrix::getDiagonalElement( const ElementWiseAssemblyMatrix& matrix ) const {
  return matrix(0,0) * TWO_POWER_D;
}


double matrixfree::stencil::ElementMatrix::getDiagonalElement( const Stencil& stencil ) const {
  #ifdef Dim2
  assertionEquals1( _centralElement, 4, stencil );
  #endif
  return stencil(_centralElement);
}



matrixfree::stencil::Stencil
matrixfree::stencil::ElementMatrix::reconstructStencil(const ElementWiseAssemblyMatrix& matrix ) const {
  matrixfree::stencil::Stencil result(0.0);

  dfor2(j)
    dfor2(k)
      tarch::la::Vector<DIMENSIONS,int> stencilEntry;
      for (int d=0; d<DIMENSIONS; d++) {
        if (j(d)>k(d)) {
          stencilEntry(d)=2;
        }
        else if (j(d)<k(d)) {
          stencilEntry(d)=0;
        }
        else {
          stencilEntry(d)=1;
        }
      }
      result(
        getStencilEntryLinearisedIndex(stencilEntry)
      ) += matrix(jScalar,kScalar);
    enddforx
  enddforx

  return result;
}


matrixfree::stencil::ElementWiseAssemblyMatrix
matrixfree::stencil::ElementMatrix::getElementWiseAssemblyMatrix( const VectorOfStencils& vectorOfStencils ) const {
  /**
   * @todo Die Abbildung sollte man im Konstruktor einmal bauen und dann hier
   * nur noch anwenden. Deshalb ist das Ding ja eine Methode und kein
   * statisches Ding.
   */
  ElementWiseAssemblyMatrix result;

  dfor2(j)
  dfor2(i)
    const int                         stencilOffset = jScalar * THREE_POWER_D;
    tarch::la::Vector<DIMENSIONS,int> stencilEntry;
    double                            commonFacesPowerTwo = 1.0;
    for (int d=0; d<DIMENSIONS; d++) {
      if (i(d)==j(d)) {
        stencilEntry(d)      = 1;
        commonFacesPowerTwo *= 2.0;
      }
      else if (i(d)<j(d)) {
        stencilEntry(d)      = 0;
      }
      else {
        stencilEntry(d)      = 2;
      }
    }
    const int vectorOfStencilIndex = stencilOffset+peano::utils::dLinearised(stencilEntry,3);
    result(jScalar,iScalar) = vectorOfStencils( vectorOfStencilIndex ) / commonFacesPowerTwo;
  enddforx
  enddforx

  return result;
}
