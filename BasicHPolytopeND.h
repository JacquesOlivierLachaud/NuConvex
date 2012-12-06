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
 * @file BasicHPolytopeND.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2012/12/06
 *
 * Header file for module BasicHPolytopeND.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(BasicHPolytopeND_RECURSES)
#error Recursive header files inclusion detected in BasicHPolytopeND.h
#else // defined(BasicHPolytopeND_RECURSES)
/** Prevents recursive inclusion of headers. */
#define BasicHPolytopeND_RECURSES

#if !defined BasicHPolytopeND_h
/** Prevents repeated inclusion of headers. */
#define BasicHPolytopeND_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <vector>
#include "DGtal/base/Common.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class BasicHPolytopeND
  /**
   * Description of template class 'BasicHPolytopeND' <p> \brief Aim:
   * This class represents a polytope as an intersection of
   * half-planes, i.e. a H-convex representation. It supports addition
   * of half-planes.
   */
  template <typename TVector>
  class BasicHPolytopeND
  {
    // ----------------------- Standard types ------------------------------
  public:
    typedef TVector Vector;
    typedef typename Vector::Component Component;
    typedef Component Scalar;

    /**
       Defines the space N.x <= u, where \a N is the normal vector and
       \a u is the upper bound.
    */
    struct ClosedHalfSpace {
      Vector normal;
      Scalar upperBound;

      ClosedHalfSpace( const Vector & theNormal, Scalar theUpperBound );

      /**
	 @param v any vector of the space.
	 @return the value N.v - u
      */
      Scalar value( const Vector & v ) const;

      /**
	 @param v any vector of the space.
	 @return 'true' iff N.v <= u
      */
      bool operator()( const Vector & v ) const;
    };

    typedef std::vector<ClosedHalfSpace> Container;
    typedef typename Container::size_type Size;
    typedef typename Container::const_iterator ConstIterator;
    
    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor.
     */
    ~BasicHPolytopeND();

    /**
     * Constructor. The empty H-polytope is the whole space.
     */
    BasicHPolytopeND();

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    BasicHPolytopeND ( const BasicHPolytopeND & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    BasicHPolytopeND & operator= ( const BasicHPolytopeND & other );

    /**
       Adds the given half space to the set of half-spaces defining the polytope.
       @param chs a closed half space.
    */
    void add( const ClosedHalfSpace & chs );

    /**
       Complexity is O(size).
       @param x some vector in space.
       @return 'true' iff x belongs to the polytope.
    */
    bool operator()( const Vector & x ) const;

    /**
       @return 'true' iff the H-polytope is the intersection of 0
       half-planes, i.e. it is the whole set.
    */
    bool empty() const;

    /**
       @return the number of intersected half-planes that defines the H-polytope.
    */
    Size size() const;

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
  protected:
    /**
       The container that stores the half-spaces that defines the polytope.
    */
    Container myHalfSpaces;

    // ------------------------- Private Datas --------------------------------
  private:

    // ------------------------- Hidden services ------------------------------
  protected:


    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class BasicHPolytopeND


  /**
   * Overloads 'operator<<' for displaying objects of class 'BasicHPolytopeND'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'BasicHPolytopeND' to write.
   * @return the output stream after the writing.
   */
  template <typename TVector>
  std::ostream&
  operator<< ( std::ostream & out, const BasicHPolytopeND<TVector> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "BasicHPolytopeND.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined BasicHPolytopeND_h

#undef BasicHPolytopeND_RECURSES
#endif // else defined(BasicHPolytopeND_RECURSES)
