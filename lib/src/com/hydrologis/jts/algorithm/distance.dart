part of dart_jts;

/**
 * A ConnectedElementPointFilter extracts a single point
 * from each connected element in a Geometry
 * (e.g. a polygon, linestring or point)
 * and returns them in a list. The elements of the list are
 * {@link org.locationtech.jts.operation.distance.GeometryLocation}s.
 *
 * @version 1.7
 */
class ConnectedElementLocationFilter implements GeometryFilter {
  /**
   * Returns a list containing a point from each Polygon, LineString, and Point
   * found inside the specified geometry. Thus, if the specified geometry is
   * not a GeometryCollection, an empty list will be returned. The elements of the list
   * are {@link org.locationtech.jts.operation.distance.GeometryLocation}s.
   */
  static List getLocations(Geometry geom) {
    List locations = [];
    geom.applyGF(new ConnectedElementLocationFilter(locations));
    return locations;
  }

  List locations;

  ConnectedElementLocationFilter(List locations) {
    this.locations = locations;
  }

  void filter(Geometry geom) {
    if (geom is Point || geom is LineString || geom is Polygon) locations.add(new GeometryLocation(geom, 0, geom.getCoordinate()));
  }
}

/// Functions to compute distance between basic geometric structures.
///
/// @author Martin Davis
///
class Distance {
  /// Computes the distance from a line segment AB to a line segment CD
  ///
  /// Note: NON-ROBUST!
  ///
  /// @param A
  ///          a point of one line
  /// @param B
  ///          the second point of (must be different to A)
  /// @param C
  ///          one point of the line
  /// @param D
  ///          another point of the line (must be different to A)
  static double segmentToSegment(Coordinate A, Coordinate B, Coordinate C, Coordinate D) {
    // check for zero-length segments
    if (A.equals(B)) return Distance.pointToSegment(A, C, D);
    if (C.equals(D)) return Distance.pointToSegment(D, A, B);

    // AB and CD are line segments
    /*
     * from comp.graphics.algo
     * 
     * Solving the above for r and s yields 
     * 
     *     (Ay-Cy)(Dx-Cx)-(Ax-Cx)(Dy-Cy) 
     * r = ----------------------------- (eqn 1) 
     *     (Bx-Ax)(Dy-Cy)-(By-Ay)(Dx-Cx)
     * 
     *     (Ay-Cy)(Bx-Ax)-(Ax-Cx)(By-Ay)  
     * s = ----------------------------- (eqn 2)
     *     (Bx-Ax)(Dy-Cy)-(By-Ay)(Dx-Cx) 
     *     
     * Let P be the position vector of the
     * intersection point, then 
     *   P=A+r(B-A) or 
     *   Px=Ax+r(Bx-Ax) 
     *   Py=Ay+r(By-Ay) 
     * By examining the values of r & s, you can also determine some other limiting
     * conditions: 
     *   If 0<=r<=1 & 0<=s<=1, intersection exists 
     *      r<0 or r>1 or s<0 or s>1 line segments do not intersect 
     *   If the denominator in eqn 1 is zero, AB & CD are parallel 
     *   If the numerator in eqn 1 is also zero, AB & CD are collinear.
     */

    bool noIntersection = false;
    if (!Envelope.intersectsEnvelopeCoords(A, B, C, D)) {
      noIntersection = true;
    } else {
      double denom = (B.x - A.x) * (D.y - C.y) - (B.y - A.y) * (D.x - C.x);

      if (denom == 0) {
        noIntersection = true;
      } else {
        double r_num = (A.y - C.y) * (D.x - C.x) - (A.x - C.x) * (D.y - C.y);
        double s_num = (A.y - C.y) * (B.x - A.x) - (A.x - C.x) * (B.y - A.y);

        double s = s_num / denom;
        double r = r_num / denom;

        if ((r < 0) || (r > 1) || (s < 0) || (s > 1)) {
          noIntersection = true;
        }
      }
    }
    if (noIntersection) {
      return MathUtils.min(
          Distance.pointToSegment(A, C, D), Distance.pointToSegment(B, C, D), Distance.pointToSegment(C, A, B), Distance.pointToSegment(D, A, B));
    }
    // segments intersect
    return 0.0;
  }

  /// Computes the distance from a point to a sequence of line segments.
  ///
  /// @param p
  ///          a point
  /// @param line
  ///          a sequence of contiguous line segments defined by their vertices
  /// @return the minimum distance between the point and the line segments
  static double pointToSegmentString(Coordinate p, List<Coordinate> line) {
    if (line.length == 0) throw new ArgumentError("Line array must contain at least one vertex");
    // this handles the case of length = 1
    double minDistance = p.distance(line[0]);
    for (int i = 0; i < line.length - 1; i++) {
      double dist = Distance.pointToSegment(p, line[i], line[i + 1]);
      if (dist < minDistance) {
        minDistance = dist;
      }
    }
    return minDistance;
  }

  /// Computes the distance from a point p to a line segment AB
  ///
  /// Note: NON-ROBUST!
  ///
  /// @param p
  ///          the point to compute the distance for
  /// @param A
  ///          one point of the line
  /// @param B
  ///          another point of the line (must be different to A)
  /// @return the distance from p to line segment AB
  static double pointToSegment(Coordinate p, Coordinate A, Coordinate B) {
    // if start = end, then just compute distance to one of the endpoints
    if (A.x == B.x && A.y == B.y) return p.distance(A);

    // otherwise use comp.graphics.algorithms Frequently Asked Questions method
    /*
     * (1) r = AC dot AB 
     *         --------- 
     *         ||AB||^2 
     *         
     * r has the following meaning: 
     *   r=0 P = A 
     *   r=1 P = B 
     *   r<0 P is on the backward extension of AB 
     *   r>1 P is on the forward extension of AB 
     *   0<r<1 P is interior to AB
     */

    double len2 = (B.x - A.x) * (B.x - A.x) + (B.y - A.y) * (B.y - A.y);
    double r = ((p.x - A.x) * (B.x - A.x) + (p.y - A.y) * (B.y - A.y)) / len2;

    if (r <= 0.0) return p.distance(A);
    if (r >= 1.0) return p.distance(B);

    /*
     * (2) s = (Ay-Cy)(Bx-Ax)-(Ax-Cx)(By-Ay) 
     *         ----------------------------- 
     *                    L^2
     * 
     * Then the distance from C to P = |s|*L.
     * 
     * This is the same calculation as {@link #distancePointLinePerpendicular}.
     * Unrolled here for performance.
     */
    double s = ((A.y - p.y) * (B.x - A.x) - (A.x - p.x) * (B.y - A.y)) / len2;
    return s.abs() * math.sqrt(len2);
  }

  /// Computes the perpendicular distance from a point p to the (infinite) line
  /// containing the points AB
  ///
  /// @param p
  ///          the point to compute the distance for
  /// @param A
  ///          one point of the line
  /// @param B
  ///          another point of the line (must be different to A)
  /// @return the distance from p to line AB
  static double pointToLinePerpendicular(Coordinate p, Coordinate A, Coordinate B) {
    // use comp.graphics.algorithms Frequently Asked Questions method
    /*
     * (2) s = (Ay-Cy)(Bx-Ax)-(Ax-Cx)(By-Ay) 
     *         ----------------------------- 
     *                    L^2
     * 
     * Then the distance from C to P = |s|*L.
     */
    double len2 = (B.x - A.x) * (B.x - A.x) + (B.y - A.y) * (B.y - A.y);
    double s = ((A.y - p.y) * (B.x - A.x) - (A.x - p.x) * (B.y - A.y)) / len2;

    return s.abs() * math.sqrt(len2);
  }
}

/**
 * Represents the location of a point on a Geometry.
 * Maintains both the actual point location
 * (which may not be exact, if the point is not a vertex)
 * as well as information about the component
 * and segment index where the point occurs.
 * Locations inside area Geometrys will not have an associated segment index,
 * so in this case the segment index will have the sentinel value of
 * {@link #INSIDE_AREA}.
 *
 * @version 1.7
 */
class GeometryLocation {
  /**
   * A special value of segmentIndex used for locations inside area geometries.
   * These locations are not located on a segment,
   * and thus do not have an associated segment index.
   */
  static final int INSIDE_AREA = -1;

  Geometry component = null;
  int segIndex;
  Coordinate pt = null;

  /**
   * Constructs a GeometryLocation specifying a point on a geometry, as well as the
   * segment that the point is on
   * (or {@link #INSIDE_AREA} if the point is not on a segment).
   *
   * @param component the component of the geometry containing the point
   * @param segIndex the segment index of the location, or INSIDE_AREA
   * @param pt the coordinate of the location
   */
  GeometryLocation(Geometry component, int segIndex, Coordinate pt) {
    this.component = component;
    this.segIndex = segIndex;
    this.pt = pt;
  }

  /**
   * Constructs a GeometryLocation specifying a point inside an area geometry.
   *
   * @param component the component of the geometry containing the point
   * @param pt the coordinate of the location
   */
  GeometryLocation.fromPointInArea(Geometry component, Coordinate pt) : this(component, INSIDE_AREA, pt);

  /**
   * Returns the geometry component on (or in) which this location occurs.
   */
  Geometry getGeometryComponent() {
    return component;
  }

  /**
   * Returns the segment index for this location. If the location is inside an
   * area, the index will have the value {@link #INSIDE_AREA};
   *
   * @return the segment index for the location, or INSIDE_AREA
   */
  int getSegmentIndex() {
    return segIndex;
  }

  /**
   * Returns the {@link Coordinate} of this location.
   */
  Coordinate getCoordinate() {
    return pt;
  }

  /**
   * Tests whether this location represents a point inside an area geometry.
   */
  bool isInsideArea() {
    return segIndex == INSIDE_AREA;
  }

  String toString() {
    return component.getGeometryType() + "[$segIndex]-${WKTWriter.toPoint(pt)}";
  }
}

/**
 * Find two points on two {@link Geometry}s which lie
 * within a given distance, or else are the nearest points
 * on the geometries (in which case this also
 * provides the distance between the geometries).
 * <p>
 * The distance computation also finds a pair of points in the input geometries
 * which have the minimum distance between them.
 * If a point lies in the interior of a line segment,
 * the coordinate computed is a close
 * approximation to the exact point.
 * <p>
 * The algorithms used are straightforward O(n^2)
 * comparisons.  This worst-case performance could be improved on
 * by using Voronoi techniques or spatial indexes.
 *
 * @version 1.7
 */
class DistanceOp {
  /**
   * Compute the distance between the nearest points of two geometries.
   * @param g0 a {@link Geometry}
   * @param g1 another {@link Geometry}
   * @return the distance between the geometries
   */
  static double distanceStatic(Geometry g0, Geometry g1) {
    DistanceOp distOp = new DistanceOp(g0, g1);
    return distOp.distance();
  }

  /**
   * Test whether two geometries lie within a given distance of each other.
   * @param g0 a {@link Geometry}
   * @param g1 another {@link Geometry}
   * @param distance the distance to test
   * @return true if g0.distance(g1) &lt;= distance
   */
  static bool isWithinDistanceStatic(Geometry g0, Geometry g1, double distance) {
    // check envelope distance for a short-circuit negative result
    double envDist = g0.getEnvelopeInternal().distance(g1.getEnvelopeInternal());
    if (envDist > distance) return false;

    // MD - could improve this further with a positive short-circuit based on envelope MinMaxDist

    DistanceOp distOp = new DistanceOp.withTerminateDistance(g0, g1, distance);
    return distOp.distance() <= distance;
  }

  /**
   * Compute the the nearest points of two geometries.
   * The points are presented in the same order as the input Geometries.
   *
   * @param g0 a {@link Geometry}
   * @param g1 another {@link Geometry}
   * @return the nearest points in the geometries
   */
  static List<Coordinate> nearestPointsStatic(Geometry g0, Geometry g1) {
    DistanceOp distOp = new DistanceOp(g0, g1);
    return distOp.nearestPoints();
  }

  /**
   * Compute the the closest points of two geometries.
   * The points are presented in the same order as the input Geometries.
   *
   * @param g0 a {@link Geometry}
   * @param g1 another {@link Geometry}
   * @return the closest points in the geometries
   * @deprecated renamed to nearestPoints
   */
  static List<Coordinate> closestPointsStatic(Geometry g0, Geometry g1) {
    DistanceOp distOp = new DistanceOp(g0, g1);
    return distOp.nearestPoints();
  }

  // input
  List<Geometry> geom;
  double terminateDistance = 0.0;

  // working
  PointLocator ptLocator = new PointLocator();
  List<GeometryLocation> minDistanceLocation;
  double minDistance = double.maxFinite;

  /**
   * Constructs a DistanceOp that computes the distance and nearest points between
   * the two specified geometries.
   * @param g0 a Geometry
   * @param g1 a Geometry
   */
  DistanceOp(Geometry g0, Geometry g1) : this.withTerminateDistance(g0, g1, 0.0);

  /**
   * Constructs a DistanceOp that computes the distance and nearest points between
   * the two specified geometries.
   * @param g0 a Geometry
   * @param g1 a Geometry
   * @param terminateDistance the distance on which to terminate the search
   */
  DistanceOp.withTerminateDistance(Geometry g0, Geometry g1, double terminateDistance) {
    this.geom = List(2);
    geom[0] = g0;
    geom[1] = g1;
    this.terminateDistance = terminateDistance;
  }

  /**
   * Report the distance between the nearest points on the input geometries.
   *
   * @return the distance between the geometries
   * or 0 if either input geometry is empty
   * @throws IllegalArgumentException if either input geometry is null
   */
  double distance() {
    if (geom[0] == null || geom[1] == null) throw new ArgumentError("null geometries are not supported");
    if (geom[0].isEmpty() || geom[1].isEmpty()) return 0.0;

    computeMinDistance();
    return minDistance;
  }

  /**
   * Report the coordinates of the nearest points in the input geometries.
   * The points are presented in the same order as the input Geometries.
   *
   * @return a pair of {@link Coordinate}s of the nearest points
   */
  List<Coordinate> nearestPoints() {
    computeMinDistance();
    List<Coordinate> nearestPts = [minDistanceLocation[0].getCoordinate(), minDistanceLocation[1].getCoordinate()];
    return nearestPts;
  }

  /**
   *
   * @return a pair of {@link Coordinate}s of the nearest points
   * @deprecated renamed to nearestPoints
   */
  List<Coordinate> closestPoints() {
    return nearestPoints();
  }

  /**
   * Report the locations of the nearest points in the input geometries.
   * The locations are presented in the same order as the input Geometries.
   *
   * @return a pair of {@link GeometryLocation}s for the nearest points
   */
  List<GeometryLocation> nearestLocations() {
    computeMinDistance();
    return minDistanceLocation;
  }

  /**
   *
   * @return a pair of {@link GeometryLocation}s for the nearest points
   * @deprecated renamed to nearestLocations
   */
  List<GeometryLocation> closestLocations() {
    return nearestLocations();
  }

  void updateMinDistance(List<GeometryLocation> locGeom, bool flip) {
    // if not set then don't update
    if (locGeom[0] == null) return;

    if (flip) {
      minDistanceLocation[0] = locGeom[1];
      minDistanceLocation[1] = locGeom[0];
    } else {
      minDistanceLocation[0] = locGeom[0];
      minDistanceLocation[1] = locGeom[1];
    }
  }

  void computeMinDistance() {
    // only compute once!
    if (minDistanceLocation != null) return;

    minDistanceLocation = List(2);
    computeContainmentDistance();
    if (minDistance <= terminateDistance) return;
    computeFacetDistance();
  }

  void computeContainmentDistance() {
    List<GeometryLocation> locPtPoly = List(2);
    // test if either geometry has a vertex inside the other
    computeContainmentDistance1(0, locPtPoly);
    if (minDistance <= terminateDistance) return;
    computeContainmentDistance1(1, locPtPoly);
  }

  void computeContainmentDistance1(int polyGeomIndex, List<GeometryLocation> locPtPoly) {
    Geometry polyGeom = geom[polyGeomIndex];
    // if no polygon then nothing to do
    if (polyGeom.getDimension() < 2) return;

    int locationsIndex = 1 - polyGeomIndex;
    List polys = PolygonExtracter.getPolygons(polyGeom);
    if (polys.length > 0) {
      List insideLocs = ConnectedElementLocationFilter.getLocations(geom[locationsIndex]);
      computeContainmentDistance2(insideLocs, polys, locPtPoly);
      if (minDistance <= terminateDistance) {
        // this assigment is determined by the order of the args in the computeInside call above
        minDistanceLocation[locationsIndex] = locPtPoly[0];
        minDistanceLocation[polyGeomIndex] = locPtPoly[1];
        return;
      }
    }
  }

  void computeContainmentDistance2(List locs, List polys, List<GeometryLocation> locPtPoly) {
    for (int i = 0; i < locs.length; i++) {
      GeometryLocation loc = locs[i];
      for (int j = 0; j < polys.length; j++) {
        computeContainmentDistance3(loc, polys[j], locPtPoly);
        if (minDistance <= terminateDistance) return;
      }
    }
  }

  void computeContainmentDistance3(GeometryLocation ptLoc, Polygon poly, List<GeometryLocation> locPtPoly) {
    Coordinate pt = ptLoc.getCoordinate();
    // if pt is not in exterior, distance to geom is 0
    if (Location.EXTERIOR != ptLocator.locate(pt, poly)) {
      minDistance = 0.0;
      locPtPoly[0] = ptLoc;
      locPtPoly[1] = new GeometryLocation.fromPointInArea(poly, pt);
      ;
      return;
    }
  }

  /**
   * Computes distance between facets (lines and points)
   * of input geometries.
   *
   */
  void computeFacetDistance() {
    List<GeometryLocation> locGeom = List(2);

    /**
     * Geometries are not wholely inside, so compute distance from lines and points
     * of one to lines and points of the other
     */
    List lines0 = LinearComponentExtracter.getLines(geom[0]);
    List lines1 = LinearComponentExtracter.getLines(geom[1]);

    List pts0 = PointExtracter.getPoints(geom[0]);
    List pts1 = PointExtracter.getPoints(geom[1]);

    // exit whenever minDistance goes LE than terminateDistance
    computeMinDistanceLines(lines0, lines1, locGeom);
    updateMinDistance(locGeom, false);
    if (minDistance <= terminateDistance) return;

    locGeom[0] = null;
    locGeom[1] = null;
    computeMinDistanceLinesPoints(lines0, pts1, locGeom);
    updateMinDistance(locGeom, false);
    if (minDistance <= terminateDistance) return;

    locGeom[0] = null;
    locGeom[1] = null;
    computeMinDistanceLinesPoints(lines1, pts0, locGeom);
    updateMinDistance(locGeom, true);
    if (minDistance <= terminateDistance) return;

    locGeom[0] = null;
    locGeom[1] = null;
    computeMinDistancePoints(pts0, pts1, locGeom);
    updateMinDistance(locGeom, false);
  }

  void computeMinDistanceLines(List lines0, List lines1, List<GeometryLocation> locGeom) {
    for (int i = 0; i < lines0.length; i++) {
      LineString line0 = lines0[i];
      for (int j = 0; j < lines1.length; j++) {
        LineString line1 = lines1[j];
        computeMinDistanceLineLineList(line0, line1, locGeom);
        if (minDistance <= terminateDistance) return;
      }
    }
  }

  void computeMinDistancePoints(List points0, List points1, List<GeometryLocation> locGeom) {
    for (int i = 0; i < points0.length; i++) {
      Point pt0 = points0[i];
      for (int j = 0; j < points1.length; j++) {
        Point pt1 = points1[j];
        double dist = pt0.getCoordinate().distance(pt1.getCoordinate());
        if (dist < minDistance) {
          minDistance = dist;
          locGeom[0] = new GeometryLocation(pt0, 0, pt0.getCoordinate());
          locGeom[1] = new GeometryLocation(pt1, 0, pt1.getCoordinate());
        }
        if (minDistance <= terminateDistance) return;
      }
    }
  }

  void computeMinDistanceLinesPoints(List lines, List points, List<GeometryLocation> locGeom) {
    for (int i = 0; i < lines.length; i++) {
      LineString line = lines[i];
      for (int j = 0; j < points.length; j++) {
        Point pt = points[j];
        computeMinDistanceLinePointList(line, pt, locGeom);
        if (minDistance <= terminateDistance) return;
      }
    }
  }

  void computeMinDistanceLineLineList(LineString line0, LineString line1, List<GeometryLocation> locGeom) {
    if (line0.getEnvelopeInternal().distance(line1.getEnvelopeInternal()) > minDistance) return;
    List<Coordinate> coord0 = line0.getCoordinates();
    List<Coordinate> coord1 = line1.getCoordinates();
    // brute force approach!
    for (int i = 0; i < coord0.length - 1; i++) {
      for (int j = 0; j < coord1.length - 1; j++) {
        double dist = Distance.segmentToSegment(coord0[i], coord0[i + 1], coord1[j], coord1[j + 1]);
        if (dist < minDistance) {
          minDistance = dist;
          LineSegment seg0 = new LineSegment.fromCoordinates(coord0[i], coord0[i + 1]);
          LineSegment seg1 = new LineSegment.fromCoordinates(coord1[j], coord1[j + 1]);
          List<Coordinate> closestPt = seg0.closestPoints(seg1);
          locGeom[0] = new GeometryLocation(line0, i, closestPt[0]);
          locGeom[1] = new GeometryLocation(line1, j, closestPt[1]);
        }
        if (minDistance <= terminateDistance) return;
      }
    }
  }

  void computeMinDistanceLinePointList(LineString line, Point pt, List<GeometryLocation> locGeom) {
    if (line.getEnvelopeInternal().distance(pt.getEnvelopeInternal()) > minDistance) return;
    List<Coordinate> coord0 = line.getCoordinates();
    Coordinate coord = pt.getCoordinate();
    // brute force approach!
    for (int i = 0; i < coord0.length - 1; i++) {
      double dist = Distance.pointToSegment(coord, coord0[i], coord0[i + 1]);
      if (dist < minDistance) {
        minDistance = dist;
        LineSegment seg = new LineSegment.fromCoordinates(coord0[i], coord0[i + 1]);
        Coordinate segClosestPoint = seg.closestPoint(coord);
        locGeom[0] = new GeometryLocation(line, i, segClosestPoint);
        locGeom[1] = new GeometryLocation(pt, 0, coord);
      }
      if (minDistance <= terminateDistance) return;
    }
  }
}
