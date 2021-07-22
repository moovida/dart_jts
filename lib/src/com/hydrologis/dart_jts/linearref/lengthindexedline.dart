part of dart_jts;

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
 * Supports linear referencing along a linear {@link Geometry}
 * using the length along the line as the index.
 * Negative length values are taken as measured in the reverse direction
 * from the end of the geometry.
 * Out-of-range index values are handled by clamping
 * them to the valid range of values.
 * Non-simple lines (i.e. which loop back to cross or touch
 * themselves) are supported.
 */
class LengthIndexedLine {
  Geometry linearGeom;

  /**
   * Constructs an object which allows a linear {@link Geometry}
   * to be linearly referenced using length as an index.
   *
   * @param linearGeom the linear geometry to reference along
   */
  LengthIndexedLine(this.linearGeom);

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
  Coordinate extractPoint(double index) {
    LinearLocation loc = LengthLocationMap.getLocationStatic(linearGeom, index);
    return loc.getCoordinate(linearGeom);
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
  Coordinate extractPointWithOffset(double index, double offsetDistance) {
    LinearLocation loc = LengthLocationMap.getLocationStatic(linearGeom, index);
    LinearLocation locLow = loc.toLowest(linearGeom);
    return locLow
        .getSegment(linearGeom)
        .pointAlongOffset(locLow.getSegmentFraction(), offsetDistance);
  }

  /**
   * Computes the {@link LineString} for the interval
   * on the line between the given indices.
   * If the endIndex lies before the startIndex,
   * the computed geometry is reversed.
   *
   * @param startIndex the index of the start of the interval
   * @param endIndex the index of the end of the interval
   * @return the linear interval between the indices
   */
  Geometry extractLine(double startIndex, double endIndex) {
    double startIndex2 = clampIndex(startIndex);
    double endIndex2 = clampIndex(endIndex);
    // if extracted line is zero-length, resolve start lower as well to ensure they are equal
    bool resolveStartLower = startIndex2 == endIndex2;
    LinearLocation startLoc =
        locationOfWithResolver(startIndex2, resolveStartLower);
//    LinearLocation endLoc = locationOf(endIndex2, true);
//    LinearLocation startLoc = locationOf(startIndex2);
    LinearLocation endLoc = locationOf(endIndex2);
    return ExtractLineByLocation.extractStatic(linearGeom, startLoc, endLoc);
  }

  LinearLocation locationOf(double index) {
    return LengthLocationMap.getLocationStatic(linearGeom, index);
  }

  LinearLocation locationOfWithResolver(double index, bool resolveLower) {
    return LengthLocationMap.getLocationStaticWithResolver(
        linearGeom, index, resolveLower);
  }

  /**
   * Computes the minimum index for a point on the line.
   * If the line is not simple (i.e. loops back on itself)
   * a single point may have more than one possible index.
   * In this case, the smallest index is returned.
   *
   * The supplied point does not <i>necessarily</i> have to lie precisely
   * on the line, but if it is far from the line the accuracy and
   * performance of this function is not guaranteed.
   * Use {@link #project} to compute a guaranteed result for points
   * which may be far from the line.
   *
   * @param pt a point on the line
   * @return the minimum index of the point
   *
   * @see #project(Coordinate)
   */
  double indexOf(Coordinate pt) {
    return LengthIndexOfPoint.indexOfStatic(linearGeom, pt);
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
  double indexOfAfter(Coordinate pt, double minIndex) {
    return LengthIndexOfPoint.indexOfAfterStatic(linearGeom, pt, minIndex);
  }

  /**
   * Computes the indices for a subline of the line.
   * (The subline must <b>conform</b> to the line; that is,
   * all vertices in the subline (except possibly the first and last)
   * must be vertices of the line and occur in the same order).
   *
   * @param subLine a subLine of the line
   * @return a pair of indices for the start and end of the subline.
   */
  List<double> indicesOf(Geometry subLine) {
    List<LinearLocation> locIndex =
        LocationIndexOfLine.indicesOfStatic(linearGeom, subLine);
    List<double> index = [
      LengthLocationMap.getLengthStatic(linearGeom, locIndex[0]),
      LengthLocationMap.getLengthStatic(linearGeom, locIndex[1])
    ];
    return index;
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
  double project(Coordinate pt) {
    return LengthIndexOfPoint.indexOfStatic(linearGeom, pt);
  }

  /**
   * Returns the index of the start of the line
   * @return the start index
   */
  double getStartIndex() {
    return 0.0;
  }

  /**
   * Returns the index of the end of the line
   * @return the end index
   */
  double getEndIndex() {
    return linearGeom.getLength();
  }

  /**
   * Tests whether an index is in the valid index range for the line.
   *
   * @param index the index to test
   * @return <code>true</code> if the index is in the valid range
   */
  bool isValidIndex(double index) {
    return (index >= getStartIndex() && index <= getEndIndex());
  }

  /**
   * Computes a valid index for this line
   * by clamping the given index to the valid range of index values
   *
   * @return a valid index value
   */
  double clampIndex(double index) {
    double posIndex = positiveIndex(index);
    double startIndex = getStartIndex();
    if (posIndex < startIndex) return startIndex;

    double endIndex = getEndIndex();
    if (posIndex > endIndex) return endIndex;

    return posIndex;
  }

  double positiveIndex(double index) {
    if (index >= 0.0) return index;
    return linearGeom.getLength() + index;
  }
}
