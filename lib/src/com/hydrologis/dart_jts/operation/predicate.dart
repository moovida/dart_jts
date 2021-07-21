part of dart_jts;

/**
 * Optimized implementation of the <tt>contains</tt> spatial predicate
 * for cases where the first {@link Geometry} is a rectangle.
 * This class works for all input geometries, including
 * {@link GeometryCollection}s.
 * <p>
 * As a further optimization,
 * this class can be used to test
 * many geometries against a single
 * rectangle in a slightly more efficient way.
 *
 * @version 1.7
 */
class RectangleContains {
  /**
   * Tests whether a rectangle contains a given geometry.
   *
   * @param rectangle a rectangular Polygon
   * @param b a Geometry of any type
   * @return true if the geometries intersect
   */
  static bool containsStatic(Polygon rectangle, Geometry b) {
    RectangleContains rc = new RectangleContains(rectangle);
    return rc.contains(b);
  }

  late Envelope rectEnv;

  /**
   * Create a new contains computer for two geometries.
   *
   * @param rectangle a rectangular geometry
   */
  RectangleContains(Polygon rectangle) {
    rectEnv = rectangle.getEnvelopeInternal();
  }

  bool contains(Geometry geom) {
    // the test geometry must be wholly contained in the rectangle envelope
    if (!rectEnv.containsEnvelope(geom.getEnvelopeInternal())) return false;

    /**
     * Check that geom is not contained entirely in the rectangle boundary.
     * According to the somewhat odd spec of the SFS, if this
     * is the case the geometry is NOT contained.
     */
    if (isContainedInBoundary(geom)) return false;
    return true;
  }

  bool isContainedInBoundary(Geometry geom) {
    // polygons can never be wholely contained in the boundary
    if (geom is Polygon) return false;
    if (geom is Point) return isPointContainedInBoundary(geom as Point);
    if (geom is LineString)
      return isLineStringContainedInBoundary(geom as LineString);

    for (int i = 0; i < geom.getNumGeometries(); i++) {
      Geometry comp = geom.getGeometryN(i);
      if (!isContainedInBoundary(comp)) return false;
    }
    return true;
  }

  bool isPointContainedInBoundary(Point point) {
    return isPointContainedInBoundaryCoord(point.getCoordinate()!);
  }

  /**
   * Tests if a point is contained in the boundary of the target rectangle.
   *
   * @param pt the point to test
   * @return true if the point is contained in the boundary
   */
  bool isPointContainedInBoundaryCoord(Coordinate pt) {
    /**
     * contains = false iff the point is properly contained in the rectangle.
     *
     * This code assumes that the point lies in the rectangle envelope
     */
    return pt.x == rectEnv.getMinX() ||
        pt.x == rectEnv.getMaxX() ||
        pt.y == rectEnv.getMinY() ||
        pt.y == rectEnv.getMaxY();
  }

  /**
   * Tests if a linestring is completely contained in the boundary of the target rectangle.
   * @param line the linestring to test
   * @return true if the linestring is contained in the boundary
   */
  bool isLineStringContainedInBoundary(LineString line) {
    CoordinateSequence seq = line.getCoordinateSequence();
    Coordinate p0 = new Coordinate.empty2D();
    Coordinate p1 = new Coordinate.empty2D();
    for (int i = 0; i < seq.size() - 1; i++) {
      seq.getCoordinateInto(i, p0);
      seq.getCoordinateInto(i + 1, p1);

      if (!isLineSegmentContainedInBoundary(p0, p1)) return false;
    }
    return true;
  }

  /**
   * Tests if a line segment is contained in the boundary of the target rectangle.
   * @param p0 an endpoint of the segment
   * @param p1 an endpoint of the segment
   * @return true if the line segment is contained in the boundary
   */
  bool isLineSegmentContainedInBoundary(Coordinate p0, Coordinate p1) {
    if (p0.equals(p1)) return isPointContainedInBoundaryCoord(p0);

    // we already know that the segment is contained in the rectangle envelope
    if (p0.x == p1.x) {
      if (p0.x == rectEnv.getMinX() || p0.x == rectEnv.getMaxX()) return true;
    } else if (p0.y == p1.y) {
      if (p0.y == rectEnv.getMinY() || p0.y == rectEnv.getMaxY()) return true;
    }
    /**
     * Either
     *   both x and y values are different
     * or
     *   one of x and y are the same, but the other ordinate is not the same as a boundary ordinate
     *
     * In either case, the segment is not wholely in the boundary
     */
    return false;
  }
}

/**
 * Implementation of the <tt>intersects</tt> spatial predicate
 * optimized for the case where one {@link Geometry} is a rectangle.
 * This class works for all
 * input geometries, including {@link GeometryCollection}s.
 * <p>
 * As a further optimization,
 * this class can be used in batch style
 * to test many geometries
 * against a single rectangle.
 *
 * @version 1.7
 */
class RectangleIntersects {
  /**
   * Tests whether a rectangle intersects a given geometry.
   *
   * @param rectangle
   *          a rectangular Polygon
   * @param b
   *          a Geometry of any type
   * @return true if the geometries intersect
   */
  static bool intersectsStatic(Polygon rectangle, Geometry b) {
    RectangleIntersects rp = new RectangleIntersects(rectangle);
    return rp.intersects(b);
  }

  Polygon rectangle;

  late Envelope rectEnv;

  /**
   * Create a new intersects computer for a rectangle.
   *
   * @param rectangle
   *          a rectangular Polygon
   */
  RectangleIntersects(this.rectangle) {
    rectEnv = rectangle.getEnvelopeInternal();
  }

  /**
   * Tests whether the given Geometry intersects
   * the query rectangle.
   *
   * @param geom the Geometry to test (may be of any type)
   * @return true if the geometry intersects the query rectangle
   */
  bool intersects(Geometry geom) {
    if (!rectEnv.intersectsEnvelope(geom.getEnvelopeInternal())) return false;

    /**
     * Test if rectangle envelope intersects any component envelope.
     * This handles Point components as well
     */
    EnvelopeIntersectsVisitor visitor = new EnvelopeIntersectsVisitor(rectEnv);
    visitor.applyTo(geom);
    if (visitor.intersects()) return true;

    /**
     * Test if any rectangle vertex is contained in the target geometry
     */
    GeometryContainsPointVisitor ecpVisitor =
        new GeometryContainsPointVisitor(rectangle);
    ecpVisitor.applyTo(geom);
    if (ecpVisitor.containsPoint()) return true;

    /**
     * Test if any target geometry line segment intersects the rectangle
     */
    RectangleIntersectsSegmentVisitor riVisitor =
        new RectangleIntersectsSegmentVisitor(rectangle);
    riVisitor.applyTo(geom);
    if (riVisitor.intersects()) return true;

    return false;
  }
}

/**
 * Tests whether it can be concluded that a rectangle intersects a geometry,
 * based on the relationship of the envelope(s) of the geometry.
 *
 * @author Martin Davis
 * @version 1.7
 */
class EnvelopeIntersectsVisitor extends ShortCircuitedGeometryVisitor {
  Envelope rectEnv;

  bool _intersects = false;

  EnvelopeIntersectsVisitor(this.rectEnv);

  /**
   * Reports whether it can be concluded that an intersection occurs,
   * or whether further testing is required.
   *
   * @return true if an intersection must occur
   * or false if no conclusion about intersection can be made
   */
  bool intersects() {
    return _intersects;
  }

  void visit(Geometry element) {
    Envelope elementEnv = element.getEnvelopeInternal();

    // disjoint => no intersection
    if (!rectEnv.intersectsEnvelope(elementEnv)) {
      return;
    }
    // rectangle contains target env => must intersect
    if (rectEnv.containsEnvelope(elementEnv)) {
      _intersects = true;
      return;
    }
    /**
     * Since the envelopes intersect and the test element is connected, if the
     * test envelope is completely bisected by an edge of the rectangle the
     * element and the rectangle must touch (This is basically an application of
     * the Jordan Curve Theorem). The alternative situation is that the test
     * envelope is "on a corner" of the rectangle envelope, i.e. is not
     * completely bisected. In this case it is not possible to make a conclusion
     * about the presence of an intersection.
     */
    if (elementEnv.getMinX() >= rectEnv.getMinX() &&
        elementEnv.getMaxX() <= rectEnv.getMaxX()) {
      _intersects = true;
      return;
    }
    if (elementEnv.getMinY() >= rectEnv.getMinY() &&
        elementEnv.getMaxY() <= rectEnv.getMaxY()) {
      _intersects = true;
      return;
    }
  }

  bool isDone() {
    return intersects == true;
  }
}

/**
 * A visitor which tests whether it can be
 * concluded that a geometry contains a vertex of
 * a query geometry.
 *
 * @author Martin Davis
 * @version 1.7
 */
class GeometryContainsPointVisitor extends ShortCircuitedGeometryVisitor {
  late CoordinateSequence rectSeq;

  late Envelope rectEnv;

  bool _containsPoint = false;

  GeometryContainsPointVisitor(Polygon rectangle) {
    this.rectSeq = rectangle.getExteriorRing().getCoordinateSequence();
    rectEnv = rectangle.getEnvelopeInternal();
  }

  /**
   * Reports whether it can be concluded that a corner point of the rectangle is
   * contained in the geometry, or whether further testing is required.
   *
   * @return true if a corner point is contained
   * or false if no conclusion about intersection can be made
   */
  bool containsPoint() {
    return _containsPoint;
  }

  void visit(Geometry geom) {
    // if test geometry is not polygonal this check is not needed
    if (!(geom is Polygon)) return;

    // skip if envelopes do not intersect
    Envelope elementEnv = geom.getEnvelopeInternal();
    if (!rectEnv.intersectsEnvelope(elementEnv)) return;

    // test each corner of rectangle for inclusion
    Coordinate rectPt = new Coordinate.empty2D();
    for (int i = 0; i < 4; i++) {
      rectSeq.getCoordinateInto(i, rectPt);
      if (!elementEnv.containsCoordinate(rectPt)) continue;
      // check rect point in poly (rect is known not to touch polygon at this
      // point)
      if (SimplePointInAreaLocator.containsPointInPolygon(
          rectPt, geom as Polygon)) {
        _containsPoint = true;
        return;
      }
    }
  }

  bool isDone() {
    return containsPoint == true;
  }
}

/**
 * A visitor to test for intersection between the query
 * rectangle and the line segments of the geometry.
 *
 * @author Martin Davis
 *
 */
class RectangleIntersectsSegmentVisitor extends ShortCircuitedGeometryVisitor {
  late Envelope rectEnv;
  late RectangleLineIntersector rectIntersector;

  bool hasIntersection = false;
  Coordinate p0 = new Coordinate.empty2D();
  Coordinate p1 = new Coordinate.empty2D();

  /**
   * Creates a visitor for checking rectangle intersection
   * with segments
   *
   * @param rectangle the query rectangle
   */
  RectangleIntersectsSegmentVisitor(Polygon rectangle) {
    rectEnv = rectangle.getEnvelopeInternal();
    rectIntersector = new RectangleLineIntersector(rectEnv);
  }

  /**
   * Reports whether any segment intersection exists.
   *
   * @return true if a segment intersection exists
   * or false if no segment intersection exists
   */
  bool intersects() {
    return hasIntersection;
  }

  void visit(Geometry geom) {
    /**
     * It may be the case that the rectangle and the
     * envelope of the geometry component are disjoint,
     * so it is worth checking this simple condition.
     */
    Envelope elementEnv = geom.getEnvelopeInternal();
    if (!rectEnv.intersectsEnvelope(elementEnv)) return;

    // check segment intersections
    // get all lines from geometry component
    // (there may be more than one if it's a multi-ring polygon)
    List lines = LinearComponentExtracter.getLines(geom);
    checkIntersectionWithLineStrings(lines);
  }

  void checkIntersectionWithLineStrings(List lines) {
    for (Iterator i = lines.iterator; i.moveNext();) {
      LineString testLine = i.current as LineString;
      checkIntersectionWithSegments(testLine);
      if (hasIntersection) return;
    }
  }

  void checkIntersectionWithSegments(LineString testLine) {
    CoordinateSequence seq1 = testLine.getCoordinateSequence();
    for (int j = 1; j < seq1.size(); j++) {
      seq1.getCoordinateInto(j - 1, p0);
      seq1.getCoordinateInto(j, p1);

      if (rectIntersector.intersects(p0, p1)) {
        hasIntersection = true;
        return;
      }
    }
  }

  bool isDone() {
    return hasIntersection == true;
  }
}
