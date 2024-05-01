part of dart_jts;

/**
 * Computes an overlay where one input is Point(s) and one is not.
 * This class supports overlay being used as an efficient way
 * to find points within or outside a polygon.
 * <p>
 * Input semantics are:
 * <ul>
 * <li>Duplicates are removed from Point output
 * <li>Non-point output is rounded and noded using the given precision model
 * </ul>
 * Output semantics are:
 * <ul>
 * <ii>An empty result is an empty atomic geometry
 *     with dimension determined by the inputs and the operation,
 *     as per overlay semantics<li>
 * </ul>
 * For efficiency the following optimizations are used:
 * <ul>
 * <li>Input points are not included in the noding of the non-point input geometry
 * (in particular, they do not participate in snap-rounding if that is used).
 * <li>If the non-point input geometry is not included in the output
 * it is not rounded and noded.  This means that points
 * are compared to the non-rounded geometry.
 * This will be apparent in the result.
 * </ul>
 *
 * @author Martin Davis
 *
 */
class OverlayMixedPoints {
  static Geometry? overlay(
      int opCode, Geometry geom0, Geometry geom1, PrecisionModel pm) {
    OverlayMixedPoints overlay = OverlayMixedPoints(opCode, geom0, geom1, pm);
    return overlay.getResult();
  }

  int opCode;
  PrecisionModel pm;
  late Geometry geomPoint;
  late Geometry geomNonPointInput;
  late GeometryFactory geometryFactory;
  late bool isPointRHS;

  Geometry? geomNonPoint;
  int? geomNonPointDim;
  PointOnGeometryLocator? locator;
  late int resultDim;

  OverlayMixedPoints(this.opCode, Geometry geom0, Geometry geom1, this.pm) {
    geometryFactory = geom0.getFactory();
    resultDim = OverlayUtil.resultDimension(
        opCode, geom0.getDimension(), geom1.getDimension());

    // name the dimensional geometries
    if (geom0.getDimension() == 0) {
      this.geomPoint = geom0;
      this.geomNonPointInput = geom1;
      this.isPointRHS = false;
    } else {
      this.geomPoint = geom1;
      this.geomNonPointInput = geom0;
      this.isPointRHS = true;
    }
  }

  Geometry? getResult() {
    // reduce precision of non-point input, if required
    geomNonPoint = prepareNonPoint(geomNonPointInput);
    geomNonPointDim = geomNonPoint!.getDimension();
    locator = createLocator(geomNonPoint!);

    List<Coordinate> coords = extractCoordinates(geomPoint, pm);

    switch (opCode) {
      case OverlayNG.INTERSECTION:
        return computeIntersection(coords);
      case OverlayNG.UNION:
      case OverlayNG.SYMDIFFERENCE:
        // UNION and SYMDIFFERENCE have same output
        return computeUnion(coords);
      case OverlayNG.DIFFERENCE:
        return computeDifference(coords);
    }
    Assert.shouldNeverReachHere("Unknown overlay op code");
    return null;
  }

  PointOnGeometryLocator createLocator(Geometry geomNonPoint) {
    if (geomNonPointDim == 2) {
      return IndexedPointInAreaLocator(geomNonPoint);
    } else {
      return IndexedPointOnLineLocator(geomNonPoint);
    }
  }

  Geometry? prepareNonPoint(Geometry geomInput) {
    // if non-point not in output no need to node it
    if (resultDim == 0) {
      return geomInput;
    }

    // Node and round the non-point geometry for output
    Geometry? geomPrep = OverlayNG.union(geomNonPointInput, pm);
    return geomPrep;
  }

  Geometry computeIntersection(List<Coordinate> coords) {
    return createPointResult(findPoints(true, coords));
  }

  Geometry computeUnion(List<Coordinate> coords) {
    List<Point> resultPointList = findPoints(false, coords);
    List<LineString> resultLineList = List.empty();
    if (geomNonPointDim == 1) {
      resultLineList = extractLines(geomNonPoint!);
    }
    List<Polygon> resultPolyList = List.empty();
    if (geomNonPointDim == 2) {
      resultPolyList = extractPolygons(geomNonPoint!);
    }

    return OverlayUtil.createResultGeometry(
        resultPolyList, resultLineList, resultPointList, geometryFactory);
  }

  Geometry? computeDifference(List<Coordinate> coords) {
    if (isPointRHS) {
      return copyNonPoint();
    }
    return createPointResult(findPoints(false, coords));
  }

  Geometry createPointResult(List<Point> points) {
    if (points.isEmpty) {
      return geometryFactory.createEmpty(0);
    } else if (points.length == 1) {
      return points[0];
    }
    return geometryFactory.createMultiPoint(points);
  }

  List<Point> findPoints(bool isCovered, List<Coordinate> coords) {
    Set<Coordinate> resultCoords = HashSet<Coordinate>();
    // keep only points contained
    for (Coordinate coord in coords) {
      if (hasLocation(isCovered, coord)) {
        // copy coordinate to avoid aliasing
        resultCoords.add(coord.copy());
      }
    }
    return createPoints(resultCoords);
  }

  List<Point> createPoints(Set<Coordinate> coords) {
    List<Point> points = List.empty(growable: true);
    for (Coordinate coord in coords) {
      Point point = geometryFactory.createPoint(coord);
      points.add(point);
    }
    return points;
  }

  bool hasLocation(bool isCovered, Coordinate coord) {
    bool isExterior = Location.EXTERIOR == locator!.locate(coord);
    if (isCovered) {
      return !isExterior;
    }
    return isExterior;
  }

  /**
   * Copy the non-point input geometry if not
   * already done by precision reduction process.
   *
   * @return a copy of the non-point geometry
   */
  Geometry? copyNonPoint() {
    if (geomNonPointInput != geomNonPoint) return geomNonPoint;
    return geomNonPoint!.copy();
  }

  static List<Coordinate> extractCoordinates(
      Geometry points, PrecisionModel pm) {
    CoordinateList coords = CoordinateList();
    points.applyCF(CoordinateFilterMixedPoints(pm, coords));
    return coords.toCoordinateArray(true);
  }

  static List<Polygon> extractPolygons(Geometry geom) {
    List<Polygon> list = List.empty(growable: true);
    for (int i = 0; i < geom.getNumGeometries(); i++) {
      Polygon poly = geom.getGeometryN(i) as Polygon;
      if (!poly.isEmpty()) {
        list.add(poly);
      }
    }
    return list;
  }

  static List<LineString> extractLines(Geometry geom) {
    List<LineString> list = List.empty(growable: true);
    for (int i = 0; i < geom.getNumGeometries(); i++) {
      LineString line = geom.getGeometryN(i) as LineString;
      if (!line.isEmpty()) {
        list.add(line);
      }
    }
    return list;
  }
}

class CoordinateFilterMixedPoints extends CoordinateFilter {
  CoordinateList coords;
  PrecisionModel pm;
  CoordinateFilterMixedPoints(this.pm, this.coords);
  @override
  void filter(Coordinate? coord) {
    if (coord != null) {
      Coordinate? p = OverlayUtil.roundCoordinate(coord, pm);
      coords.addCoord(p!, false);
    }
  }
}
