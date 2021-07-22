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
 * Determines the location of a subline along a linear {@link Geometry}.
 * The location is reported as a pair of {@link LinearLocation}s.
 * <p>
 * <b>Note:</b> Currently this algorithm is not guaranteed to
 * return the correct substring in some situations where
 * an endpoint of the test line occurs more than once in the input line.
 * (However, the common case of a ring is always handled correctly).
 */
class LocationIndexOfLine {
  /**
  * MD - this algorithm has been extracted into a class
  * because it is intended to validate that the subline truly is a subline,
  * and also to use the internal vertex information to unambiguously locate the subline.
  */
  static List<LinearLocation> indicesOfStatic(
      Geometry linearGeom, Geometry subLine) {
    LocationIndexOfLine locater = new LocationIndexOfLine(linearGeom);
    return locater.indicesOf(subLine);
  }

  Geometry linearGeom;

  LocationIndexOfLine(this.linearGeom);

  List<LinearLocation> indicesOf(Geometry subLine) {
    Coordinate startPt =
        (subLine.getGeometryN(0) as LineString).getCoordinateN(0);
    LineString lastLine =
        subLine.getGeometryN(subLine.getNumGeometries() - 1) as LineString;
    Coordinate endPt = lastLine.getCoordinateN(lastLine.getNumPoints() - 1);

    LocationIndexOfPoint locPt = new LocationIndexOfPoint(linearGeom);
    List<LinearLocation> subLineLoc = [];
    subLineLoc[0] = locPt.indexOf(startPt);

    // check for case where subline is zero length
    if (subLine.getLength() == 0.0) {
      subLineLoc[1] = subLineLoc[0].copy();
    } else {
      subLineLoc[1] = locPt.indexOfAfter(endPt, subLineLoc[0]);
    }
    return subLineLoc;
  }
}
