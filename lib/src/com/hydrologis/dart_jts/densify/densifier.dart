import 'package:dart_jts/dart_jts.dart';


/*
 * Copyright (c) 2016 Vivid Solutions.
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License 2.0
 * and Eclipse Distribution License v. 1.0 which accompanies this distribution.
 * The Eclipse Public License is available at http://www.eclipse.org/legal/epl-v20.html
 * and the Eclipse Distribution License is available at
 *
 * http://www.eclipse.org/org/documents/edl-v10.php.
 */

/**
 * Densifies a {@link Geometry} by inserting extra vertices along the line segments
 * contained in the geometry.
 * All segments in the created densified geometry will be <b>no longer</b>
 * than the given distance tolerance
 * (that is, all segments in the output will have length less than or equal to
 * the distance tolerance).
 * The coordinates created during densification respect the input geometry's
 * {@link PrecisionModel}.
 * <p>
 * By default polygonal results are processed to ensure they are valid.
 * This processing is costly, and it is very rare for results to be invalid.
 * Validation processing can be disabled by calling the {@link #setValidate(bool)} method.
 * <p>
 * <b>Note:</b> At some future point this class will
 * offer a variety of densification strategies.
 *
 * @author Martin Davis
 */
class Densifier {
  /**
   * Densifies a geometry using a given distance tolerance,
   * and respecting the input geometry's {@link PrecisionModel}.
   *
   * @param geom the geometry to densify
   * @param distanceTolerance the distance tolerance to densify
   * @return the densified geometry
   */
  static Geometry densify(Geometry geom, double distanceTolerance) {
    Densifier densifier = new Densifier(geom);
    densifier.setDistanceTolerance(distanceTolerance);
    return densifier.getResultGeometry();
  }

  /**
   * Densifies a list of coordinates.
   *
   * @param pts the coordinate list
   * @param distanceTolerance the densify tolerance
   * @return the densified coordinate sequence
   */
  static List<Coordinate> densifyPoints(List<Coordinate> pts,
      double distanceTolerance, PrecisionModel precModel) {
    LineSegment seg = new LineSegment.empty();
    CoordinateList coordList = new CoordinateList();
    for (int i = 0; i < pts.length - 1; i++) {
      seg.p0 = pts[i];
      seg.p1 = pts[i + 1];
      coordList.addCoord(seg.p0, false);
      double len = seg.getLength();

      // check if no densification is required
      if (len <= distanceTolerance) continue;

      // densify the segment
      int densifiedSegCount = ((len / distanceTolerance)).ceil();
      double densifiedSegLen = len / densifiedSegCount;
      for (int j = 1; j < densifiedSegCount; j++) {
        double segFract = (j * densifiedSegLen) / len;
        Coordinate p = seg.pointAlong(segFract);
        precModel.makeCoordinatePrecise(p);
        coordList.addCoord(p, false);
      }
    }
    // this check handles empty sequences
    if (pts.length > 0) coordList.addCoord(pts[pts.length - 1], false);
    return coordList.toCoordinateArray(true);
  }

 late Geometry inputGeom;

  late double distanceTolerance;

  /**
   * Indicates whether areas should be topologically validated.
   */
  bool isValidated = true;

  /**
   * Creates a new densifier instance.
   *
   * @param inputGeom
   */
  Densifier(Geometry inputGeom) {
    this.inputGeom = inputGeom;
  }

  /**
   * Sets the distance tolerance for the densification. All line segments
   * in the densified geometry will be no longer than the distance tolerance.
   * The distance tolerance must be positive.
   *
   * @param distanceTolerance
   *          the densification tolerance to use
   */
  void setDistanceTolerance(double distanceTolerance) {
    if (distanceTolerance <= 0.0)
      throw new Exception("Tolerance must be positive");
    this.distanceTolerance = distanceTolerance;
  }

  /**
   * Sets whether polygonal results are processed to ensure they are valid.
   *
   * @param isValidated true if the results should be validated
   */
  void setValidate(bool isValidated) {
    this.isValidated = isValidated;
  }

  /**
   * Gets the densified geometry.
   *
   * @return the densified geometry
   */
  Geometry getResultGeometry() {
    return (new DensifyTransformer(distanceTolerance, isValidated))
        .transform(inputGeom);
  }
}

class DensifyTransformer extends GeometryTransformer {
  late double distanceTolerance;
  late bool isValidated;

  DensifyTransformer(double distanceTolerance, bool isValidated) {
    this.distanceTolerance = distanceTolerance;
    this.isValidated = isValidated;
  }

  CoordinateSequence transformCoordinates(
      CoordinateSequence coords, Geometry parent) {
    List<Coordinate> inputPts = coords.toCoordinateArray();
    List<Coordinate> newPts = Densifier.densifyPoints(
        inputPts, distanceTolerance, parent.getPrecisionModel());
    // prevent creation of invalid linestrings
    if (parent is LineString && newPts.length == 1) {
      newPts = [Coordinate.empty2D()];
    }
    return factory.getCoordinateSequenceFactory().create(newPts);
  }

   Geometry transformPolygon(Polygon geom, Geometry parent) {
    Geometry roughGeom = super.transformPolygon(geom, parent);
    // don't try and correct if the parent is going to do this
    if (parent is MultiPolygon) {
      return roughGeom;
    }
    return createValidArea(roughGeom);
  }

   Geometry transformMultiPolygon(MultiPolygon geom, Geometry parent) {
    Geometry roughGeom = super.transformMultiPolygon(geom, parent);
    return createValidArea(roughGeom);
  }

  /**
   * Creates a valid area geometry from one that possibly has bad topology
   * (i.e. self-intersections). Since buffer can handle invalid topology, but
   * always returns valid geometry, constructing a 0-width buffer "corrects"
   * the topology. Note this only works for area geometries, since buffer
   * always returns areas. This also may return empty geometries, if the input
   * has no actual area.
   *
   * @param roughAreaGeom
   *          an area geometry possibly containing self-intersections
   * @return a valid area geometry
   */
  Geometry createValidArea(Geometry roughAreaGeom) {
    // if valid no need to process to make valid
    if (!isValidated || roughAreaGeom.isValid()) return roughAreaGeom;
    return roughAreaGeom.buffer(0.0);
  }
}
