// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#ifndef _MATRIXFREE_STENCIl_ELEMENT_MATRIX_H_
#define _MATRIXFREE_STENCIl_ELEMENT_MATRIX_H_


#include "matrixfree/stencil/Stencil.h"
#include "peano/grid/VertexEnumerator.h"


namespace matrixfree {
  namespace stencil {
    class ElementMatrix;
  }
}


/**
 * Derive Element-wise Matrices From Stencils
 *
 * Many components in the Peano framework switch from a stencil notation to
 * element-wise assembly matrices frequently. As both representations (stencils
 * and element-wise matrices) are equivalent, it is reasonable to code only one
 * representation, and to derive the alternative representation automatically.
 * This class provides the required operation.
 *
 * - Stencils are stored within vertices (red lines for a 5-point stencil in
 *   2d, e.g.).
 * - However, sometimes we don't have access to neighbouring cells (big blue
 *   cell, e.g., is one coarse cell where the solver automaton descends to
 *   finer grids)
 * - Consequently, we have to convert the @f$ 2^d @f$ vertex stencils within
 *   each cell into a element-wise assembly matrix @f$ \mathbb{R}^{2^d \times 2^d}@f$.
 *   With this matrix, we then can accumulate the residuals.
 * - PeProt automatically generates a routine that takes the @f$ 2^d @f$
 *   vertices and writes their @f$ 3^d @f$ stencils into one big vector of
 *   length @f$ 2^d \cdot 3^d @f$.
 * - This class provides a routine that transforms this big vector into an
 *   element assembly matrix.
 * - There's also a shortcut operaration that takes only one stencil. In this
 *   case, we assume that the stencil is the same for all @f$ 2^d @f$ adjacent
 *   vertices of the cell.
 *
 * @image html ElementWiseMatrix.png
 *
 * !!! Example
 *
 * We have the stencil
 *
 * @f$ \left[ \begin{array}{ccc} -1 & -1 & -1 \\ -1 & 8 & -1 \\ -1 & -1 & -1 \end{array} \right] @f$
 *
 * for (d=2). It is stored as a vector of length nine. Within a cell, we have
 * four vertices. The values u of these vertices can be converted into one
 * vector of length four due to a call of the static operation
 * \code
tarch::la::Vector<TWO_POWER_D,double> u = MyVertices::readU(...);
 \endcode
 * Such a static operation is generated, if you add the statement
 * \code
  read scalar: U
 \endcode
 * to your Peano specification file. Now, if you wanna compute matrix times
 * vector, you basically compute $Au$ with A being a @f$ 4 \times 4 @f$
 * matrix. A is well-defined by the stencil. In this example, it would be
 *
 * @f$ \left[ \begin{array}{cccc} 2 & -0.5 & -0.5 & -1.0 \\ -0.5 & 2 & -1.0 & -0.5 \\ -0.5 & -1.0 & 2 & -0.5 \\ -1.0 & -0.5 & -0.5 & 2 \end{array} \right] @f$
 *
 * The static methods of this class convert your stencil into the element-wise
 * stiffness matrix. There's two variants of conversion operations:
 *
 * - The simple variant takes one stencil and creates the stiffness matrix. It
 *   resembles the example.
 * - A more flexible variant takes @f$ 2^d @f$ stencils and creates the
 *   element-wise system matrix. This is the method of choice if you have
 *   different stencils on the different vertices.
 *
 *
 * @author Tobias Weinzierl
 */
class matrixfree::stencil::ElementMatrix {
  private:
    int _centralElement;
  public:
    ElementMatrix();

    /**
     * Derive element-wise matrix from stencil.
     */
    ElementWiseAssemblyMatrix getElementWiseAssemblyMatrix( const Stencil& stencil ) const;

    /**
     * @todo Docu -> das ist eine Rekonstruktion, die u.U. falsch ist (bei springenden Stencils)
     */
    Stencil reconstructStencil(const ElementWiseAssemblyMatrix& matrix ) const;

    /**
     * Derive element-wise matrix from @f$ 2^d @f$ stencils, i.e. this operation
     * takes one big @f$ 2^d \cdot 3^d @f$ vector as argument and returns an
     * @f$ 2^d \times 2^d @f$ matrix.
     */
    ElementWiseAssemblyMatrix getElementWiseAssemblyMatrix( const VectorOfStencils& vectorOfStencils ) const;

    double getDiagonalElement( const ElementWiseAssemblyMatrix& matrix ) const;
    double getDiagonalElement( const Stencil& stencil ) const;
};


#endif
