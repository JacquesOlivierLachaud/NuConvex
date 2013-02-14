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
 * @file MaximalPlaneSummary.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2012/12/07
 *
 * Header file for module MaximalPlaneSummary.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(MaximalPlaneSummary_RECURSES)
#error Recursive header files inclusion detected in MaximalPlaneSummary.h
#else // defined(MaximalPlaneSummary_RECURSES)
/** Prevents recursive inclusion of headers. */
#define MaximalPlaneSummary_RECURSES

#if !defined MaximalPlaneSummary_h
/** Prevents repeated inclusion of headers. */
#define MaximalPlaneSummary_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/kernel/CSpace.h"
#include "DGtal/kernel/SimpleMatrix.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class MaximalPlaneSummary
  /**
   * Description of template class 'MaximalPlaneSummary' <p> \brief
   * Aim: Stores in a compact way a lot of information on a maximal
   * plane, like normal vector, upper and lower bounding plane,
   * maximal radius, and the projected areas.
   *
   * Model of boost::DefaultConstructible, boost::CopyConstructible,
   * boost::Assignable.
   *
   * @tparam TSpace the type of digital space in which the maximal
   * plane is defined.
   *
   * @see NuConvexSet
   */
  template <typename TSpace>
  struct MaximalPlaneSummary
  {
  public:
    BOOST_CONCEPT_ASSERT(( CSpace< TSpace > ));
    // BOOST_STATIC_ASSERT(( TSpace::dimension == 3 ));

    typedef TSpace Space;
    typedef typename Space::Point Point;
    typedef typename Space::Vector Vector;
    typedef typename Space::RealPoint RealPoint;
    typedef typename Space::RealVector RealVector;
    typedef typename Space::Size Size;
    typedef typename RealVector::Component Scalar;
    typedef SimpleMatrix< Scalar, 3, 3 > RealMatrix;
    
    // ------------------------- Public Datas ------------------------------
  public:
    RealPoint center;      /**< the coordinates of the point which
                              defines the maximal plane. */
    unsigned int majorAxis;/**< the major axis of this plane. */
    Scalar radius;         /**< the maximal radius of the plane. */
    RealVector normal;     /**< the unit normal direction. */
    Scalar upper;          /**< upper plane offset: normal.x = upper */
    Scalar lower;          /**< lower plane offset: normal.x = lower */
    Scalar projectedArea;  /**< the area of the projected nu-convex set. */
    Size nb;               /**< the number of points. */
    RealVector centroid;   /**< the centroid of the maximal plane. */
    RealVector eigenvalues;/**< the eigenvalues of the covariance matrix. */
    RealMatrix eigenvectors;/**< the eigenvectors of the covariance matrix. */

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor.
     */
    ~MaximalPlaneSummary();

    /**
     * Constructor.
     */
    MaximalPlaneSummary();

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    MaximalPlaneSummary ( const MaximalPlaneSummary & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    MaximalPlaneSummary & operator= ( const MaximalPlaneSummary & other );

    /**
       @param p any point in the space.
       @return the projection of point \a p onto the upper plane of the maximal plane.
    */
    RealPoint projectOntoUpperPlane( const RealPoint & p ) const;

    /**
       @param p any point in the space.
       @return the projection of point \a p onto the lower plane of the maximal plane.
    */
    RealPoint projectOntoLowerPlane( const RealPoint & p ) const;

    /**
       @param p any point in the space.
       @return the projection of point \a p onto the PCA plane.
    */
    RealPoint projectOntoPCAPlane( const RealPoint & p ) const;
    RealPoint projectOntoMedianPlane( const RealPoint & p ) const;

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

    // ------------------------- Hidden services ------------------------------
  protected:


    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class MaximalPlaneSummary


  /**
   * Overloads 'operator<<' for displaying objects of class 'MaximalPlaneSummary'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'MaximalPlaneSummary' to write.
   * @return the output stream after the writing.
   */
  template <typename TSpace>
  std::ostream&
  operator<< ( std::ostream & out, const MaximalPlaneSummary<TSpace> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "MaximalPlaneSummary.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined MaximalPlaneSummary_h

#undef MaximalPlaneSummary_RECURSES
#endif // else defined(MaximalPlaneSummary_RECURSES)
