/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

#pragma once

/**
 * @file NuConvexSet.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2012/12/06
 *
 * Header file for module NuConvexSet.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(NuConvexSet_RECURSES)
#error Recursive header files inclusion detected in NuConvexSet.h
#else // defined(NuConvexSet_RECURSES)
/** Prevents recursive inclusion of headers. */
#define NuConvexSet_RECURSES

#if !defined NuConvexSet_h
/** Prevents repeated inclusion of headers. */
#define NuConvexSet_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/geometry/surfaces/COBAGenericNaivePlane.h"
#include "BasicHPolytopeND.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class NuConvexSet
  /**
   * Description of template class 'NuConvexSet' <p>
   * \brief Aim:
   */
  template < typename TSpace, typename TVisitor,
	     typename TVertexEmbedder, 
	     typename TInternalInteger = DGtal::int64_t >
  class NuConvexSet
  {
    // ----------------------- public types ------------------------------
  public:
    typedef TSpace Space;
    typedef TVisitor Visitor;
    typedef TInternalInteger InternalInteger;
    typedef TVertexEmbedder VertexEmbedder;
    typedef typename Visitor::Vertex Vertex;

    BOOST_STATIC_ASSERT(( ConceptUtils::SameType
			   < Vertex, typename VertexEmbedder::Argument >::value ));
    typedef COBAGenericNaivePlane< Space, InternalInteger > GenericNaivePlane;
    typedef typename VertexEmbedder::Value RealPoint;
    typedef typename RealPoint::Coordinate Scalar;
    typedef RealPoint RealVector;
    typedef BasicHPolytopeND<RealVector> HPolytope;
    typedef typename Visitor::Node Node;     // vertex
    typedef typename Visitor::Scalar Distance; // distance


    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor.
     */
    ~NuConvexSet();

    /**
       Constructor.
       @param visitor the visitor is duplicated.
       @param embedder the embedder is referenced.
     */
    NuConvexSet( const Visitor & visitor, const VertexEmbedder & embedder );

    /**
       The start vertex is given by the visitor.
       nu = p/q
       @param diameter the diameter is necessary for COBA digital plane algorithm.
    */
    void init( InternalInteger p, InternalInteger q,
	       InternalInteger diameter );

    bool compute( Scalar distanceUpperBound = -1 );


    // ----------------------- Interface --------------------------------------
  public:

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;

    // ------------------------- Protected Datas ------------------------------
  private:

    // ------------------------- Private Datas --------------------------------
  private:

    Visitor myVisitor;
    const VertexEmbedder & myEmbedder;
    GenericNaivePlane myPlane;
    HPolytope myPolytope;

    // ------------------------- Hidden services ------------------------------
  protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    NuConvexSet();

  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    NuConvexSet ( const NuConvexSet & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    NuConvexSet & operator= ( const NuConvexSet & other );

    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class NuConvexSet


  /**
   * Overloads 'operator<<' for displaying objects of class 'NuConvexSet'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'NuConvexSet' to write.
   * @return the output stream after the writing.
   */
  template < typename TSpace, typename TVisitor,
	     typename TVertexEmbedder, 
	     typename TInternalInteger >
  std::ostream&
  operator<< ( std::ostream & out, 
	       const NuConvexSet< TSpace, TVisitor, TVertexEmbedder, TInternalInteger > & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "NuConvexSet.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined NuConvexSet_h

#undef NuConvexSet_RECURSES
#endif // else defined(NuConvexSet_RECURSES)
