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
 * Extracts the subline of a linear {@link Geometry} between
 * two {@link LinearLocation}s on the line.
 */
class ExtractLineByLocation {
  /**
   * Computes the subline of a {@link LineString} between
   * two {@link LinearLocation}s on the line.
   * If the start location is after the end location,
   * the computed linear geometry has reverse orientation to the input line.
   *
   * @param line the line to use as the baseline
   * @param start the start location
   * @param end the end location
   * @return the extracted subline
   */
  static Geometry extractStatic(
      Geometry line, LinearLocation start, LinearLocation end) {
    ExtractLineByLocation ls = new ExtractLineByLocation(line);
    return ls.extract(start, end);
  }

  Geometry line;

  ExtractLineByLocation(this.line);

  /**
   * Extracts a subline of the input.
   * If <code>end < start</code> the linear geometry computed will be reversed.
   *
   * @param start the start location
   * @param end the end location
   * @return a linear geometry
   */
  Geometry extract(LinearLocation start, LinearLocation end) {
    if (end.compareTo(start) < 0) {
      return reverse(computeLinear(end, start));
    }
    return computeLinear(start, end);
  }

  Geometry reverse(Geometry linear) {
    if (linear is Lineal) return linear.reverse();

    Assert.shouldNeverReachHere("non-linear geometry encountered");
    throw ArgumentError("non-linear geometry encountered");
  }

  /**
   * Assumes input is valid (e.g. start <= end)
   *
   * @param start
   * @param end
   * @return a linear geometry
   */
  LineString computeLine(LinearLocation start, LinearLocation end) {
    List<Coordinate> coordinates = line.getCoordinates();
    CoordinateList newCoordinates = new CoordinateList();

    int startSegmentIndex = start.getSegmentIndex();
    if (start.getSegmentFraction() > 0.0) startSegmentIndex += 1;
    int lastSegmentIndex = end.getSegmentIndex();
    if (end.getSegmentFraction() == 1.0) lastSegmentIndex += 1;
    if (lastSegmentIndex >= coordinates.length)
      lastSegmentIndex = coordinates.length - 1;
    // not needed - LinearLocation values should always be correct
    //Assert.isTrue(end.getSegmentFraction() <= 1.0, "invalid segment fraction value");

    if (!start.isVertex()) newCoordinates.add(start.getCoordinate(line));
    for (int i = startSegmentIndex; i <= lastSegmentIndex; i++) {
      newCoordinates.add(coordinates[i]);
    }
    if (!end.isVertex()) newCoordinates.add(end.getCoordinate(line));

    // ensure there is at least one coordinate in the result
    if (newCoordinates._backingList.length <= 0)
      newCoordinates.add(start.getCoordinate(line));

    List<Coordinate> newCoordinateArray = newCoordinates.toCoordinateArray();
    /**
     * Ensure there is enough coordinates to build a valid line.
     * Make a 2-point line with duplicate coordinates, if necessary.
     * There will always be at least one coordinate in the coordList.
     */
    if (newCoordinateArray.length <= 1) {
      newCoordinateArray = <Coordinate>[
        newCoordinateArray[0],
        newCoordinateArray[0]
      ];
    }
    return line.getFactory().createLineString(newCoordinateArray);
  }

  /**
   * Assumes input is valid (e.g. start <= end)
   *
   * @param start
   * @param end
   * @return a linear geometry
   */
  Geometry computeLinear(LinearLocation start, LinearLocation end) {
    LinearGeometryBuilder builder =
        new LinearGeometryBuilder(line.getFactory());
    builder.setFixInvalidLines(true);

    if (!start.isVertex()) builder.add(start.getCoordinate(line));

    for (LinearIterator it = new LinearIterator.withStart(line, start);
        it.hasNext();
        it.next()) {
      if (end.compareLocationValues(
              it.getComponentIndex(), it.getVertexIndex(), 0.0) <
          0) break;

      Coordinate pt = it.getSegmentStart();
      builder.add(pt);
      if (it.isEndOfLine()) builder.endLine();
    }
    if (!end.isVertex()) builder.add(end.getCoordinate(line));

    return builder.getGeometry();
  }

  /**
   * Computes a valid and normalized location
   * compatible with the values in a LinearIterator.
   * (I.e. segmentFractions of 1.0 are converted to the next highest coordinate index)
   */
  /*
   LinearLocation normalize(LinearLocation loc)
  {
    int componentIndex = loc.getComponentIndex();
    int segmentIndex = loc.getSegmentIndex();
    double segmentFraction = loc.getSegmentFraction();

    if (segmentFraction < 0.0) {
      segmentFraction = 0.0;
    }
    if (segmentFraction > 1.0) {
      segmentFraction = 1.0;
    }

    if (componentIndex < 0) {
      componentIndex = 0;
      segmentIndex = 0;
      segmentFraction = 0.0;
    }
    if (segmentIndex < 0) {
      segmentIndex = 0;
      segmentFraction = 0.0;
    }

    if (segmentFraction == 1.0) {
      segmentFraction = 0.0;
      segmentIndex += 1;
    }

    return new LinearLocation(componentIndex, segmentIndex, segmentFraction);
  }
  */
}
