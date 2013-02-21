// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#ifndef _MATRIXFREE_SOLVER_MULTIGRID_H_
#define _MATRIXFREE_SOLVER_MULTIGRID_H_

#include "tarch/logging/Log.h"
#include "tarch/la/Vector.h"
#include "tarch/la/Matrix.h"
#include "tarch/la/GramSchmidt.h"
#include "peano/utils/Globals.h"
#include "peano/utils/Loop.h"
#include "peano/grid/VertexEnumerator.h"
#include "matrixfree/stencil/ElementMatrix.h"

namespace matrixfree {
  namespace solver{
     class Multigrid;
  }
}

/**
 * Basis utility class for multigrid
 *
 * @author Marion Weinzierl, Tobias Weinzierl
 */
class matrixfree::solver::Multigrid {
  private:
  	static tarch::logging::Log _log;

	  matrixfree::stencil::ElementMatrix _elementMatrix;

    /**
     * Stencil for d-linear interpolation
     *
     * The class provides a number of standard routines for permanently
     * reoccuring problems such as d-linear interpolation. For this, it
     * internally holds a d-linear interpolation stencil.
     */
    tarch::la::Vector<FIVE_POWER_D,double>  _dLinearInterpolation;

    /**
     * The class keeps track of the number of stencil evaluations.
     */
    int                                     _numberOfStencilEvaluations;

  public:

    Multigrid();
    virtual ~Multigrid();

    void setup();

    tarch::la::Vector<TWO_POWER_D, double> prolongCellValues(
        const tarch::la::Vector<TWO_POWER_D_TIMES_FIVE_POWER_D, double>&	coarseGridVerticesP,
        const tarch::la::Vector<TWO_POWER_D_TIMES_D,int>&                  fineGridVerticesPositions,
    	tarch::la::Vector<TWO_POWER_D, double>								cellValues
    );

    tarch::la::Vector<TWO_POWER_D, double> restrictCellValues(
   	    const tarch::la::Vector<TWO_POWER_D_TIMES_FIVE_POWER_D, double>&	coarseGridVerticesR,
   	    const tarch::la::Vector<TWO_POWER_D_TIMES_D,int>&                  fineGridVerticesPositions,
   		tarch::la::Vector<TWO_POWER_D, double> 								cellValues
    );

    tarch::la::Vector<TWO_POWER_D, double> computeCellResidual(
    		tarch::la::Vector<TWO_POWER_D, double> cellValues,
    		tarch::la::Vector<TWO_POWER_D, double> rhs,
    		tarch::la::Matrix<TWO_POWER_D, TWO_POWER_D, double> elementWiseAssemblyMatrix
    );


    /**
     * !!! Compute Transfer Contribution Weight for one Individual Fine Grid Vertex
     *     
     * If you wanna compute @f$ Pu_\mbox{coarse} @f$ or @f$ Ru_\mbox{fine} @f$ (for
     * example in a a touchVertex FirstTime() event) for one individual fine
     * grid vertex, you can do this in @f$ 2^d @f$ steps: For each coarse grid
     * vertex
     * - you either compute its contribution/weight to the fine grid value and
     *   scale it with the corresponding coarse grid vertex's value
     *   (prolongation) or
     * - you compute its contribution/weight for the fine grid vertex and
     *   scale it with the fine grid vertex's value (restriction)
     * then you add this value to the fine grid vertex or the corresponding
     * coarse grid value, respectively.
     *
     * To split it up into @f$ 2^d @f$ steps does make sense,
     * as the P stencil of each vertex might be different. This operation
     * implements how to get the weight between fine grid vertex and one
     * coarse grid vertex.
     *
     * @image html Multigrid.png
     *
     * !! Usage (compute Pu, i.e. prolongation, for one individual fine grid vertex)
     *
     * - Create temporary value. It is typically a double, and it represents
     *   the new fine grid value.
     * - Run over all @f$ 2^d @f$ (typically a dfor2 loop):
     *   - We will now compute the contribution of one coarse grid vertex to
     *     our fine grid vertex.
     *   - Take the coarse grid inter-grid transfer stencil. It is a
     *     @f$ 5^d @f$ vector.
     *   - Pass the stencil, the @f$ (0,1)^d @f$ position of the current
     *     coarse grid vertex relative to the current fine grid vertex, and the
     *     fine grid vertex position within the @f$ 4^d @f$ patch to this
     *     operation. You get a weight.
     *   - Multiply this weight with the coarse grid vertex's value that is to
     *     be interpolated.
     *   - Add the result to your temporary variable.
     * - The temporary variable now holds the contributions from all
     *   @f$ 2^d @f$ coarse grid values. Set it on the fine grid vertex.
     *
     * Usage:
     * \code
  double interpolatedValue = 0.0;
  dfor2(k) // run over the coarse grid vertices (0,0), (1,0), (0,1),
           // and (1,1) of this one coarse cell if d=2. If d!=2 it
           // does a couple of iterations more. k is an integer vector.
    tarch::la::Vector<FIVE_POWER_D,double> coarseGridStencil = coming from your coarse grid vertex k;
    double                                 coarseGridValue   = coming from your coarse grid vertex k;
    interpolatedValue +=
      coarseGridValue *
      _multigrid.computeContributionWeightOfInterGridTransfer(
        k,
        coarseGridStencil,
        fineGridPositionOfVertexWithinPatch
      );
  enddforx
  fineGridVertex.setXXX( interpolatedValue );
       \endcode
     *
     *
     * !! Usage (compute Ru, i.e. restriction, for one individual fine grid vertex)
     *
     * - Create temporary value. It is typically a double, and it represents
     *   the new fine grid value.
     * - Run over all @f$ 2^d @f$ (typically a dfor2 loop):
     *   - We will now compute the contribution of one coarse grid vertex to
     *     our fine grid vertex.
     *   - Take the coarse grid inter-grid transfer stencil. It is a
     *     @f$ 5^d @f$ vector.
     *   - Pass the stencil, the @f$ (0,1)^d @f$ position of the current
     *     coarse grid vertex relative to the current fine grid vertex, and the
     *     fine grid vertex position within the @f$ 4^d @f$ patch to this
     *     operation. You get a weight.
     *   - Multiply this weight with the fine grid vertex's value that is to
     *     be interpolated.
     *   - Add the result to the current coarse grid value.
     *
     * Usage:
     * \code
  dfor2(k) // run over the coarse grid vertices (0,0), (1,0), (0,1),
           // and (1,1) of this one coarse cell if d=2. If d!=2 it
           // does a couple of iterations more. k is an integer vector.
    tarch::la::Vector<FIVE_POWER_D,double> coarseGridStencil = coming from your coarse grid vertex k;
    double contributionToCurrentCoarseGridVertex =
      fineGridValue *
      _multigrid.computeContributionWeightOfInterGridTransfer(
        k,
        coarseGridStencil,
        fineGridPositionOfVertexWithinPatch
      );
  enddforx
       \endcode
     *
     *
     * The operation is not const as it has to update the internal counters.
     */
    double computeContributionWeightOfInterGridTransfer(
      const tarch::la::Vector<DIMENSIONS,int>&       currentCoarseGridVertex,
      const tarch::la::Vector<FIVE_POWER_D,double>&  currentCoarseGridVertexsInterGridTransferOperator,
      const tarch::la::Vector<DIMENSIONS,int>&       fineGridPositionOfVertex
    );


    /*
     * Get the 4x4 cell prolongation or restriction matrix for a fine grid cell from the 5x5 prolongation operator of the coarse grid vertices.
     * For Galerkin this is without scaling (according to the adjacent cells in the stencil), as the operator A ist already scaled correctly.
     * Otherwise, you have to scale it with 1/TWO_POWER_D.
     */
    tarch::la::Matrix<TWO_POWER_D, TWO_POWER_D, double> calculateCellInterGridTransferOperator(
      const tarch::la::Vector<TWO_POWER_D_TIMES_FIVE_POWER_D, double>&		coarseGridVerticesInterGridTransferOperators,
      const tarch::la::Vector<TWO_POWER_D_TIMES_D,int>&                    fineGridVerticesPositions
    );

    /*
     * Calculate the cell-wise (i.e., four times four) Petrov-Galerkin coarse grid operators for the coarse vertices of one cell.
     * As with the elementwise assembly matrix, the contributions of the adjacent cell-wise operators have to be added in order to get the complete
     * vertex-wise (i.e., 3 times 3) operator.
     *
     * To be called in leaveCell(...)
     *
     * @param fineGridVerticesPositions One d-dimensional vector for each of
     *           the @f$ 2^d @f$ vertices.
     */
    tarch::la::Matrix<TWO_POWER_D, TWO_POWER_D, double> calculatePetrovGalerkinCoarseGridOperator(
        const tarch::la::Vector<TWO_POWER_D_TIMES_FIVE_POWER_D, double>&	coarseGridVerticesP,
        const tarch::la::Vector<TWO_POWER_D_TIMES_FIVE_POWER_D, double>&  	coarseGridVerticesR,
        const tarch::la::Vector<TWO_POWER_D_TIMES_D,int>&                  fineGridVerticesPositions,
        tarch::la::Matrix<TWO_POWER_D, TWO_POWER_D, double>                	elementWiseAssemblyMatrix
    );

    tarch::la::Vector<TWO_POWER_D_TIMES_THREE_POWER_D,double> computeBoxMGIntergridTransferOperator(tarch::la::Vector<THREE_POWER_D_TIMES_FOUR_POWER_D,double> stencils);

    tarch::la::Vector<16, double> solveGammaSystems(
        tarch::la::Vector<3, double> stencil1,
        tarch::la::Vector<3, double> stencil2,
        tarch::la::Vector<3, double> stencil3,
        tarch::la::Vector<3, double> stencil4,
        tarch::la::Vector<3, double> stencil5,
        tarch::la::Vector<3, double> stencil6,
        tarch::la::Vector<3, double> stencil7,
        tarch::la::Vector<3, double> stencil8);

    tarch::la::Vector<16, double> solveIotaSystem(
        tarch::la::Vector<9, double> stencil1,
        tarch::la::Vector<9, double> stencil2,
        tarch::la::Vector<9, double> stencil3,
        tarch::la::Vector<9, double> stencil4,
        tarch::la::Vector<16, double> gammas);

    tarch::la::Vector<2, double> solve2x2System(
        tarch::la::Vector<3, double> stencil1,
        tarch::la::Vector<3, double> stencil2);

    tarch::la::Vector<4, double> solve4x4System(
        tarch::la::Matrix<4, 4, double> a,
        tarch::la::Vector<4, double> rhs);

    int getNumberOfStencilUpdates() const;
    void clearNumberOfStencilUpdates();

    /**
     * Compute position relative to one coarse grid vertex.
     *
     * @image html StencilRelativeToCoarseGrid.png
     *
     * @param coarseGridVertexPosition Index of current coarse grid vertex. It is a
     *        d-dimensional integer vector that holds only 0 and 1. (0,0) in
     *        d=2, e.g., is the bottom left vertex, (1,1) is the top right one.
     * @param fineGridVertexPosition Fine grid position of the fine grid
     *        vertex within a @f$ 4^d @f$ patch. It is an integer vector of
     *        length d that holds entries of @f$ \{ 0,1,2,3 \} @f$.
     * @param entryOfCoarseGridStencil Is an out parameter and afterwards
     *        contains a d-dimensional integer entry out of @f$ \{ 0,1,2,3,4 \} @f$
     *        that tells you for a 5-point stencil which entry affects the fine
     *        grid vertex.
     * @param coarseGridStencilInfluencesFineGridVertex Tells you whether the
     *        5-point stencil on the coarse grid influences the fine grid vertex
     *        at all. If the fine grid vertex, e.g., is (3,0) and the coarse grid
     *        vertex is (0,1), then this flag is false.
     */
    static void getPositionIn5PowDStencilRelativeToKthCoarseVertex(
      const tarch::la::Vector<DIMENSIONS,int>&       coarseGridVertexPosition,
      const tarch::la::Vector<DIMENSIONS,int>&       fineGridVertexPosition,
      tarch::la::Vector<DIMENSIONS,int>&             entryOfCoarseGridStencil,
      int&                                           indexOfCoarseGridStencil,
      bool&                                          coarseGridStencilInfluencesFineGridVertex
    );

    static int getPositionInCellInterGridTransferOperatorVector(
        const int coarseGridVertexNumber,
        const int positionInOperator);

    static int getPositionInCellStencilVector(
        const int coarseGridVertexNumber,
        const int positionInOperator);

    tarch::la::Vector<TWO_POWER_D_TIMES_THREE_POWER_D, double> addUpdateToStencils(
    	tarch::la::Vector<TWO_POWER_D_TIMES_THREE_POWER_D, double>&		verticesStencils,
    	tarch::la::Matrix<TWO_POWER_D, TWO_POWER_D, double> 			cellWiseStencilUpdate
    );

    tarch::la::Vector<TWO_POWER_D_TIMES_FIVE_POWER_D, double> fillInIntergridTransferOperators(
        tarch::la::Vector<TWO_POWER_D_TIMES_THREE_POWER_D,double> operators3x3,
        tarch::la::Vector<TWO_POWER_D_TIMES_FIVE_POWER_D, double> operators5x5
    );

    static tarch::la::Vector<TWO_POWER_D, int> getPositionsInA(
    	const int coarseVertexNumber
    );

    static tarch::la::Vector<THREE_POWER_D, int> getPositionsInIntergridTransferOperator(
      const int coarseVertexNumber
    );

    /**
     * Simple d-linear Interpolation
     *
     * I moved this to the multigrid solver, as it is a standard routine
     * required all the time - for example by most visualisation environments.
     * It requires the @f$ 2^d @f$ coarse grid values and the vertex's position
     * within a @f$ 3^d @f$ patch. Also, the method expects a scaling factor to
     * avoid that someone forgets to scale appropriately.
     *
     * If you interpolate scalar data such as the solution, the scaling
     * typically is 1.0. If you interpolate FEM data such as the residual, the
     * scaling has to reflect the integral of the supports, i.e. if you have
     * the Laplacian, the scaling is 1.0/static_cast<double>(THREE_POWER_D).
     *
     * @param coarseGridValues         The @f$ 2^d @f$ coarse grid values to be interpolated.
     * @param fineGridPositionOfVertex Position of the vertex within the @f$ 3^d @f$ patch, i.e. an integer vector whose entries are somewhere between 0 and 4.
     */
    double getDLinearInterpolatedValue(
      const tarch::la::Vector<TWO_POWER_D,double>&           coarseGridValues,
      const tarch::la::Vector<DIMENSIONS,int>&               fineGridPositionOfVertex,
      double                                                 scaling
    );
};

#endif /* MULTIGRID_H_ */
