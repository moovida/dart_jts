part of dart_jts;

/**
 * Performs an overlay operation on inputs which are both point geometries.
 * <p>
 * Semantics are:
 * <ul>
 * <li>Points are rounded to the precision model if provided
 * <li>Points with identical XY values are merged to a single point
 * <li>Extended ordinate values are preserved in the output,
 * apart from merging
 * <li>An empty result is returned as <code>POINT EMPTY</code>
 * </ul>
 *
 * @author Martin Davis
 *
 */
class OverlayPoints {
  /**
   * Performs an overlay operation on inputs which are both point geometries.
   *
   * @param geom0 the first geometry argument
   * @param geom1 the second geometry argument
   * @param opCode the code for the desired overlay operation
   * @param pm the precision model to use
   * @return the result of the overlay operation
   */
  static Geometry? overlay(
      int opCode, Geometry geom0, Geometry geom1, PrecisionModel pm) {
    OverlayPoints overlay = OverlayPoints(opCode, geom0, geom1, pm);
    return overlay.getResult();
  }

  late int opCode;
  late Geometry geom0;
  late Geometry geom1;
  late PrecisionModel pm;
  late GeometryFactory geometryFactory;
  late List<Point> resultList;

  /**
   * Creates an instance of an overlay operation on inputs which are both point geometries.
   *
   * @param geom0 the first geometry argument
   * @param geom1 the second geometry argument
   * @param opCode the code for the desired overlay operation
   * @param pm the precision model to use
   */
  OverlayPoints(int opCode, Geometry geom0, Geometry geom1, PrecisionModel pm) {
    this.opCode = opCode;
    this.geom0 = geom0;
    this.geom1 = geom1;
    this.pm = pm;
    geometryFactory = geom0.getFactory();
  }

  /**
   * Gets the result of the overlay.
   *
   * @return the overlay result
   */
  Geometry? getResult() {
    Map<Coordinate, Point> map0 = buildPointMap(geom0);
    Map<Coordinate, Point> map1 = buildPointMap(geom1);

    resultList = List.empty(growable: true);
    switch (opCode) {
      case OverlayNG.INTERSECTION:
        computeIntersection(map0, map1, resultList);
        break;
      case OverlayNG.UNION:
        computeUnion(map0, map1, resultList);
        break;
      case OverlayNG.DIFFERENCE:
        computeDifference(map0, map1, resultList);
        break;
      case OverlayNG.SYMDIFFERENCE:
        computeDifference(map0, map1, resultList);
        computeDifference(map1, map0, resultList);
        break;
    }
    if (resultList.isEmpty) {
      return OverlayUtil.createEmptyResult(0, geometryFactory);
    } else {
      return geometryFactory.buildGeometry(resultList);
    }
  }

  void computeIntersection(Map<Coordinate, Point> map0,
      Map<Coordinate, Point> map1, List<Point> resultList) {
    for (var entry in map0.entries) {
      if (map1.containsKey(entry.key)) {
        resultList.add(copyPoint(entry.value));
      }
    }
  }

  void computeDifference(Map<Coordinate, Point> map0,
      Map<Coordinate, Point> map1, List<Point> resultList) {
    for (MapEntry<Coordinate, Point> entry in map0.entries) {
      if (!map1.containsKey(entry.key)) {
        resultList.add(copyPoint(entry.value));
      }
    }
  }

  void computeUnion(Map<Coordinate, Point> map0, Map<Coordinate, Point> map1,
      List<Point> resultList) {
    // copy all A points
    for (Point p in map0.values) {
      resultList.add(copyPoint(p));
    }

    for (MapEntry<Coordinate, Point> entry in map1.entries) {
      if (!map0.containsKey(entry.key)) {
        resultList.add(copyPoint(entry.value));
      }
    }
  }

  Point copyPoint(Point pt) {
    // if pm is floating, the point coordinate is not changed
    if (OverlayUtil.isFloating(pm)) return pt.copy() as Point;

    // pm is fixed.  Round off X&Y ordinates, copy other ordinates unchanged
    CoordinateSequence seq = pt.getCoordinateSequence();
    CoordinateSequence seq2 = seq.copy();
    seq2.setOrdinate(0, CoordinateSequence.X, pm.makePrecise(seq.getX(0)));
    seq2.setOrdinate(0, CoordinateSequence.Y, pm.makePrecise(seq.getY(0)));
    return geometryFactory.createPointSeq(seq2);
  }

  Map<Coordinate, Point> buildPointMap(Geometry geoms) {
    Map<Coordinate, Point> map = {};
    geoms.applyGCF(GeometryComponentFilterImpl(map, pm));
    return map;
  }

  /**
   * Round the key point if precision model is fixed.
   * Note: return value is only copied if rounding is performed.
   *
   * @param pt
   * @return
   */
  static Coordinate roundCoord(Point pt, PrecisionModel pm) {
    Coordinate p = pt.getCoordinates()[0];
    if (OverlayUtil.isFloating(pm)) return p;
    Coordinate p2 = p.copy();
    pm.makeCoordinatePrecise(p2);
    return p2;
  }
}

class GeometryComponentFilterImpl implements GeometryComponentFilter {
  Map<Coordinate, Point> map;
  PrecisionModel pm;
  GeometryComponentFilterImpl(this.map, this.pm);
  @override
  void filter(Geometry geom) {
    if (!(geom is Point)) return;
    if (geom.isEmpty()) return;

    Point pt = geom;
    Coordinate p = OverlayPoints.roundCoord(pt, pm);
    /**
     * Only add first occurrence of a point.
     * This provides the merging semantics of overlay
     */
    if (!map.containsKey(p)) map.putIfAbsent(p, () => pt);
  }
}
