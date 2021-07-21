part of dart_jts;

/**
 * Determines the {@link Location} of {@link Coordinate}s relative to
 * an areal geometry, using indexing for efficiency.
 * This algorithm is suitable for use in cases where
 * many points will be tested against a given area.
 * <p>
 * The Location is computed precisely, in that points
 * located on the geometry boundary or segments will
 * return {@link Location.BOUNDARY}.
 * <p>
 * {@link Polygonal} and {@link LinearRing} geometries
 * are supported.
 *
 * Thread-safe and immutable.
 *
 * @author Martin Davis
 *
 */
class IndexedPointInAreaLocator implements PointOnGeometryLocator {
  late IntervalIndexedGeometry index;

  /**
   * Creates a new locator for a given {@link Geometry}.
   * {@link Polygonal} and {@link LinearRing} geometries
   * are supported.
   *
   * @param g the Geometry to locate in
   */
  IndexedPointInAreaLocator(Geometry g) {
    if (!(g is Polygonal || g is LinearRing))
      throw new ArgumentError("Argument must be Polygonal or LinearRing");
    index = new IntervalIndexedGeometry(g);
  }

  /**
   * Determines the {@link Location} of a point in an areal {@link Geometry}.
   *
   * @param p the point to test
   * @return the location of the point in the geometry
   */
  int locate(Coordinate p) {
    RayCrossingCounter rcc = new RayCrossingCounter(p);

    SegmentVisitor visitor = new SegmentVisitor(rcc);
    index.queryWithVisitor(p.y, p.y, visitor);

    /*
     // MD - slightly slower alternative
    List segs = index.query(p.y, p.y);
    countSegs(rcc, segs);
    */

    return rcc.getLocation();
  }
}

class SegmentVisitor implements ItemVisitor {
  RayCrossingCounter counter;

  SegmentVisitor(this.counter);

  void visitItem(Object item) {
    LineSegment seg = item as LineSegment;
    counter.countSegment(seg.getCoordinate(0), seg.getCoordinate(1));
  }
}

class IntervalIndexedGeometry {
  bool isEmpty = false;
  SortedPackedIntervalRTree index = new SortedPackedIntervalRTree();

  IntervalIndexedGeometry(Geometry geom) {
    if (geom.isEmpty())
      isEmpty = true;
    else
      init(geom);
  }

  void init(Geometry geom) {
    List lines = LinearComponentExtracter.getLines(geom);
    for (Iterator i = lines.iterator; i.moveNext();) {
      LineString line = i.current as LineString;
      List<Coordinate> pts = line.getCoordinates();
      addLine(pts);
    }
  }

  void addLine(List<Coordinate> pts) {
    for (int i = 1; i < pts.length; i++) {
      LineSegment seg = new LineSegment.fromCoordinates(pts[i - 1], pts[i]);
      double min = math.min(seg.p0.y, seg.p1.y);
      double max = math.max(seg.p0.y, seg.p1.y);
      index.insert(min, max, seg);
    }
  }

  List query(double min, double max) {
    if (isEmpty) return [];

    ArrayListVisitor visitor = new ArrayListVisitor();
    index.query(min, max, visitor);
    return visitor.getItems();
  }

  void queryWithVisitor(double min, double max, ItemVisitor visitor) {
    if (isEmpty) return;
    index.query(min, max, visitor);
  }
}

/**
 * An interface for classes which determine the {@link Location} of
 * points in a {@link Geometry}.
 *
 * @author Martin Davis
 */
abstract class PointOnGeometryLocator {
  /**
   * Determines the {@link Location} of a point in the {@link Geometry}.
   *
   * @param p the point to test
   * @return the location of the point in the geometry
   */
  int locate(Coordinate p);
}

/**
 * Computes the location of points
 * relative to a {@link Polygonal} {@link Geometry},
 * using a simple <tt>O(n)</tt> algorithm.
 * <p>
 * The algorithm used reports
 * if a point lies in the interior, exterior,
 * or exactly on the boundary of the Geometry.
 * <p>
 * Instance methods are provided to implement
 * the interface {@link PointInAreaLocator}.
 * However, they provide no performance
 * advantage over the class methods.
 * <p>
 * This algorithm is suitable for use in cases where
 * only a few points will be tested.
 * If many points will be tested,
 * {@link IndexedPointInAreaLocator} may provide better performance.
 *
 * @version 1.7
 */
class SimplePointInAreaLocator implements PointOnGeometryLocator {
  /**
   * Determines the {@link Location} of a point in an areal {@link Geometry}.
   * The return value is one of:
   * <ul>
   * <li>{@link Location.INTERIOR} if the point is in the geometry interior
   * <li>{@link Location.BOUNDARY} if the point lies exactly on the boundary
   * <li>{@link Location.EXTERIOR} if the point is outside the geometry
   * </ul>
   *
   * @param p the point to test
   * @param geom the areal geometry to test
   * @return the Location of the point in the geometry
   */
  static int locatePointInGeom(Coordinate p, Geometry geom) {
    if (geom.isEmpty()) return Location.EXTERIOR;
    /**
     * Do a fast check against the geometry envelope first
     */
    if (!geom.getEnvelopeInternal().intersectsCoordinate(p))
      return Location.EXTERIOR;

    return locateInGeometry(p, geom);
  }

  /**
   * Determines whether a point is contained in a {@link Geometry},
   * or lies on its boundary.
   * This is a convenience method for
   * <pre>
   *  Location.EXTERIOR != locate(p, geom)
   * </pre>
   *
   * @param p the point to test
   * @param geom the geometry to test
   * @return true if the point lies in or on the geometry
   */
  static bool isContained(Coordinate p, Geometry geom) {
    return Location.EXTERIOR != locatePointInGeom(p, geom);
  }

  static int locateInGeometry(Coordinate p, Geometry geom) {
    if (geom is Polygon) {
      return locatePointInPolygon(p, geom);
    }

    if (geom is GeometryCollection) {
      Iterator geomi = GeometryCollectionIterator(geom);
      while (geomi.moveNext()) {
        Geometry g2 = geomi.current;
        if (g2 != geom) {
          int loc = locateInGeometry(p, g2);
          if (loc != Location.EXTERIOR) return loc;
        }
      }
    }
    return Location.EXTERIOR;
  }

  /**
   * Determines the {@link Location} of a point in a {@link Polygon}.
   * The return value is one of:
   * <ul>
   * <li>{@link Location.INTERIOR} if the point is in the geometry interior
   * <li>{@link Location.BOUNDARY} if the point lies exactly on the boundary
   * <li>{@link Location.EXTERIOR} if the point is outside the geometry
   * </ul>
   *
   * This method is provided for backwards compatibility only.
   * Use {@link #locate(Coordinate, Geometry)} instead.
   *
   * @param p the point to test
   * @param poly the geometry to test
   * @return the Location of the point in the polygon
   *
   */
  static int locatePointInPolygon(Coordinate p, Polygon poly) {
    if (poly.isEmpty()) return Location.EXTERIOR;
    LinearRing shell = poly.getExteriorRing();
    int shellLoc = locatePointInRing(p, shell);
    if (shellLoc != Location.INTERIOR) return shellLoc;

    // now test if the point lies in or on the holes
    for (int i = 0; i < poly.getNumInteriorRing(); i++) {
      LinearRing hole = poly.getInteriorRingN(i);
      int holeLoc = locatePointInRing(p, hole);
      if (holeLoc == Location.BOUNDARY) return Location.BOUNDARY;
      if (holeLoc == Location.INTERIOR) return Location.EXTERIOR;
      // if in EXTERIOR of this hole keep checking the other ones
    }
    // If not in any hole must be inside polygon
    return Location.INTERIOR;
  }

  /**
   * Determines whether a point lies in a {@link Polygon}.
   * If the point lies on the polygon boundary it is
   * considered to be inside.
   *
   * @param p the point to test
   * @param poly the geometry to test
   * @return true if the point lies in or on the polygon
   */
  static bool containsPointInPolygon(Coordinate p, Polygon poly) {
    return Location.EXTERIOR != locatePointInPolygon(p, poly);
  }

  /**
   * Determines whether a point lies in a LinearRing,
   * using the ring envelope to short-circuit if possible.
   *
   * @param p the point to test
   * @param ring a linear ring
   * @return true if the point lies inside the ring
   */
  static int locatePointInRing(Coordinate p, LinearRing ring) {
    // short-circuit if point is not in ring envelope
    if (!ring.getEnvelopeInternal().intersectsCoordinate(p))
      return Location.EXTERIOR;
    return PointLocation.locateInRing(p, ring.getCoordinates());
  }

  Geometry geom;

  /**
   * Create an instance of a point-in-area locator,
   * using the provided areal geometry.
   *
   * @param geom the areal geometry to locate in
   */
  SimplePointInAreaLocator(this.geom);

  /**
   * Determines the {@link Location} of a point in an areal {@link Geometry}.
   * The return value is one of:
   * <ul>
   * <li>{@link Location.INTERIOR} if the point is in the geometry interior
   * <li>{@link Location.BOUNDARY} if the point lies exactly on the boundary
   * <li>{@link Location.EXTERIOR} if the point is outside the geometry
   * </ul>
   *
   * @param p the point to test
   * @return the Location of the point in the geometry
   */
  int locate(Coordinate p) {
    return SimplePointInAreaLocator.locatePointInGeom(p, geom);
  }
}
