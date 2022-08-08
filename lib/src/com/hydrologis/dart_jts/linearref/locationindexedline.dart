import '../geom/coordinate.dart';
import '../geom/geometry.dart';
import '../geom/linestring.dart';
import '../geom/multilinestring.dart';
import 'extractlinebylocation.dart';
import 'linearlocation.dart';
import 'locationindexofpoint.dart';
import 'locationindexofline.dart';

/*
 * Copyright (c) 2016 Vivid Solutions.
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse  License v1.0
 * and Eclipse Distribution License v. 1.0 which accompanies this distribution.
 * The Eclipse  License is available at http://www.eclipse.org/legal/epl-v10.html
 * and the Eclipse Distribution License is available at
 *
 * http://www.eclipse.org/org/documents/edl-v10.php.
 */

/**
 * Supports linear referencing
 * along a linear {@link Geometry}
 * using {@link LinearLocation}s as the index.
 */
class LocationIndexedLine {
  Geometry linearGeom;

  /**
   * Constructs an object which allows linear referencing along
   * a given linear {@link Geometry}.
   *
   * @param linearGeom the linear geometry to reference along
   */
  LocationIndexedLine(this.linearGeom) {
    checkGeometryType();
  }

  void checkGeometryType() {
    if (!(linearGeom is LineString || linearGeom is MultiLineString))
      throw new ArgumentError("Input geometry must be linear");
  }

  /**
   * Computes the {@link Coordinate} for the point
   * on the line at the given index.
   * If the index is out of range the first or last point on the
   * line will be returned.
   * The Z-ordinate of the computed point will be interpolated from
   * the Z-ordinates of the line segment containing it, if they exist.
   *
   * @param index the index of the desired point
   * @return the Coordinate at the given index
   */
  Coordinate extractPoint(LinearLocation index) {
    return index.getCoordinate(linearGeom);
  }

  /**
   * Computes the {@link Coordinate} for the point
   * on the line at the given index, offset by the given distance.
   * If the index is out of range the first or last point on the
   * line will be returned.
   * The computed point is offset to the left of the line if the offset distance is
   * positive, to the right if negative.
   * 
   * The Z-ordinate of the computed point will be interpolated from
   * the Z-ordinates of the line segment containing it, if they exist.
   *
   * @param index the index of the desired point
   * @param offsetDistance the distance the point is offset from the segment
   *    (positive is to the left, negative is to the right)
   * @return the Coordinate at the given index
   */
  Coordinate extractPointWithOffset(LinearLocation index, double offsetDistance) {
    LinearLocation indexLow = index.toLowest(linearGeom);
    return indexLow.getSegment(linearGeom).pointAlongOffset(indexLow.getSegmentFraction(), offsetDistance);
  }

  /**
   * Computes the {@link LineString} for the interval
   * on the line between the given indices.
   * If the start location is after the end location,
   * the computed linear geometry has reverse orientation to the input line.
   *
   * @param startIndex the index of the start of the interval
   * @param endIndex the index of the end of the interval
   * @return the linear geometry between the indices
   */
  Geometry extractLine(LinearLocation startIndex, LinearLocation endIndex) {
    return ExtractLineByLocation.extractStatic(linearGeom, startIndex, endIndex);
  }

  /**
   * Computes the index for a given point on the line.
   * <p>
   * The supplied point does not <i>necessarily</i> have to lie precisely
   * on the line, but if it is far from the line the accuracy and
   * performance of this function is not guaranteed.
   * Use {@link #project} to compute a guaranteed result for points
   * which may be far from the line.
   *
   * @param pt a point on the line
   * @return the index of the point
   * @see #project(Coordinate)
   */
  LinearLocation indexOf(Coordinate pt) {
    return LocationIndexOfPoint.indexOfStatic(linearGeom, pt);
  }

  /**
   * Finds the index for a point on the line
   * which is greater than the given index.
   * If no such index exists, returns <tt>minIndex</tt>.
   * This method can be used to determine all indexes for
   * a point which occurs more than once on a non-simple line.
   * It can also be used to disambiguate cases where the given point lies
   * slightly off the line and is equidistant from two different
   * points on the line.
   *
   * The supplied point does not <i>necessarily</i> have to lie precisely
   * on the line, but if it is far from the line the accuracy and
   * performance of this function is not guaranteed.
   * Use {@link #project} to compute a guaranteed result for points
   * which may be far from the line.
   *
   * @param pt a point on the line
   * @param minIndex the value the returned index must be greater than
   * @return the index of the point greater than the given minimum index
   *
   * @see #project(Coordinate)
   */
  LinearLocation indexOfAfter(Coordinate pt, LinearLocation minIndex) {
    return LocationIndexOfPoint.indexOfAfterStatic(linearGeom, pt, minIndex);
  }

  /**
   * Computes the indices for a subline of the line.
   * (The subline must <i>conform</i> to the line; that is,
   * all vertices in the subline (except possibly the first and last)
   * must be vertices of the line and occur in the same order).
   *
   * @param subLine a subLine of the line
   * @return a pair of indices for the start and end of the subline.
   */
  List<LinearLocation> indicesOf(Geometry subLine) {
    return LocationIndexOfLine.indicesOfStatic(linearGeom, subLine);
  }

  /**
   * Computes the index for the closest point on the line to the given point.
   * If more than one point has the closest distance the first one along the line
   * is returned.
   * (The point does not necessarily have to lie precisely on the line.)
   *
   * @param pt a point on the line
   * @return the index of the point
   */
  LinearLocation project(Coordinate pt) {
    return LocationIndexOfPoint.indexOfStatic(linearGeom, pt);
  }

  /**
   * Returns the index of the start of the line
   * @return the location index
   */
  LinearLocation getStartIndex() {
    return new LinearLocation();
  }

  /**
   * Returns the index of the end of the line
   * @return the location index
   */
  LinearLocation getEndIndex() {
    return LinearLocation.getEndLocation(linearGeom);
  }

  /**
   * Tests whether an index is in the valid index range for the line.
   *
   * @param index the index to test
   * @return <code>true</code> if the index is in the valid range
   */
  bool isValidIndex(LinearLocation index) {
    return index.isValid(linearGeom);
  }

  /**
   * Computes a valid index for this line
   * by clamping the given index to the valid range of index values
   *
   * @return a valid index value
   */
  LinearLocation clampIndex(LinearLocation index) {
    LinearLocation loc = index.copy();
    loc.clamp(linearGeom);
    return loc;
  }
}
