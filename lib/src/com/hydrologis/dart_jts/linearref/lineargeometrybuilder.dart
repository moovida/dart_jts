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
 * Builds a linear geometry ({@link LineString} or {@link MultiLineString})
 * incrementally (point-by-point).
 *
 * @version 1.7
 */
class LinearGeometryBuilder {
  GeometryFactory geomFact;
  List<LineString> lines = [];
  CoordinateList? coordList = null;

  bool ignoreInvalidLines = false;
  bool fixInvalidLines = false;

  Coordinate? lastPt = null;

  LinearGeometryBuilder(this.geomFact);

  /**
   * Allows invalid lines to be ignored rather than causing Exceptions.
   * An invalid line is one which has only one unique point.
   *
   * @param ignoreInvalidLines <code>true</code> if short lines are to be ignored
   */
  void setIgnoreInvalidLines(bool ignoreInvalidLines) {
    this.ignoreInvalidLines = ignoreInvalidLines;
  }

  /**
   * Allows invalid lines to be ignored rather than causing Exceptions.
   * An invalid line is one which has only one unique point.
   *
   * @param fixInvalidLines <code>true</code> if short lines are to be ignored
   */
  void setFixInvalidLines(bool fixInvalidLines) {
    this.fixInvalidLines = fixInvalidLines;
  }

  /**
   * Adds a point to the current line.
   *
   * @param pt the Coordinate to add
   */
  void add(Coordinate pt) {
    addWithRepeated(pt, true);
  }

  /**
   * Adds a point to the current line.
   *
   * @param pt the Coordinate to add
   */
  void addWithRepeated(Coordinate pt, bool allowRepeatedPoints) {
    if (coordList == null) coordList = new CoordinateList();
    coordList!.addCoord(pt, allowRepeatedPoints);
    lastPt = pt;
  }

  Coordinate? getLastCoordinate() {
    return lastPt;
  }

  /**
   * Terminate the current LineString.
   */
  void endLine() {
    if (coordList == null) {
      return;
    }
    if (ignoreInvalidLines && coordList!._backingList.length < 2) {
      coordList = null;
      return;
    }
    List<Coordinate> rawPts = coordList!.toCoordinateArray();
    List<Coordinate> pts = rawPts;
    if (fixInvalidLines) pts = validCoordinateSequence(rawPts);

    coordList = null;
    LineString? line = null;
    try {
      line = geomFact.createLineString(pts);
    } catch (ex) {
      // exception is due to too few points in line.
      // only propagate if not ignoring short lines
      if (!ignoreInvalidLines) throw ex;
    }

    if (line != null) lines.add(line);
  }

  List<Coordinate> validCoordinateSequence(List<Coordinate> pts) {
    if (pts.length >= 2) return pts;
    List<Coordinate> validPts = [pts[0], pts[0]];
    return validPts;
  }

  Geometry getGeometry() {
    // end last line in case it was not done by user
    endLine();
    return geomFact.buildGeometry(lines);
  }
}
