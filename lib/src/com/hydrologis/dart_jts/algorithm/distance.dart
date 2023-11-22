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

  ConnectedElementLocationFilter(this.locations);

  void filter(Geometry geom) {
    if (geom is Point || geom is LineString || geom is Polygon)
      locations.add(new GeometryLocation(geom, 0, geom.getCoordinate()!));
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
  static double segmentToSegment(
      Coordinate A, Coordinate B, Coordinate C, Coordinate D) {
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
          Distance.pointToSegment(A, C, D),
          Distance.pointToSegment(B, C, D),
          Distance.pointToSegment(C, A, B),
          Distance.pointToSegment(D, A, B));
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
    if (line.length == 0)
      throw new ArgumentError("Line array must contain at least one vertex");
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
  static double pointToLinePerpendicular(
      Coordinate p, Coordinate A, Coordinate B) {
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
 * Contains a pair of points and the distance between them.
 * Provides methods to update with a new point pair with
 * either maximum or minimum distance.
 */
class PointPairDistance {
  List<Coordinate> pt = [new Coordinate.empty2D(), new Coordinate.empty2D()];
  double distance = double.nan;
  bool isNull = true;

  PointPairDistance() {}

  void initialize() {
    isNull = true;
  }

  void initialize2C(Coordinate p0, Coordinate p1) {
    pt[0].setCoordinate(p0);
    pt[1].setCoordinate(p1);
    distance = p0.distance(p1);
    isNull = false;
  }

  /**
   * Initializes the points, avoiding recomputing the distance.
   * @param p0
   * @param p1
   * @param distance the distance between p0 and p1
   */
  void initialize2CD(Coordinate p0, Coordinate p1, double distance) {
    pt[0].setCoordinate(p0);
    pt[1].setCoordinate(p1);
    this.distance = distance;
    isNull = false;
  }

  double getDistance() {
    return distance;
  }

  List<Coordinate> getCoordinates() {
    return pt;
  }

  Coordinate getCoordinate(int i) {
    return pt[i];
  }

  void setMaximum(PointPairDistance ptDist) {
    setMaximum2C(ptDist.pt[0], ptDist.pt[1]);
  }

  void setMaximum2C(Coordinate p0, Coordinate p1) {
    if (isNull) {
      initialize2C(p0, p1);
      return;
    }
    double dist = p0.distance(p1);
    if (dist > distance) initialize2CD(p0, p1, dist);
  }

  void setMinimum(PointPairDistance ptDist) {
    setMinimum2C(ptDist.pt[0], ptDist.pt[1]);
  }

  void setMinimum2C(Coordinate p0, Coordinate p1) {
    if (isNull) {
      initialize2C(p0, p1);
      return;
    }
    double dist = p0.distance(p1);
    if (dist < distance) initialize2CD(p0, p1, dist);
  }

  String toString() {
    return WKTWriter.toLineStringFromCoords(pt[0], pt[1]);
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

  Geometry component;
  int segIndex;
  Coordinate pt;

  /**
   * Constructs a GeometryLocation specifying a point on a geometry, as well as the
   * segment that the point is on
   * (or {@link #INSIDE_AREA} if the point is not on a segment).
   *
   * @param component the component of the geometry containing the point
   * @param segIndex the segment index of the location, or INSIDE_AREA
   * @param pt the coordinate of the location
   */
  GeometryLocation(this.component, this.segIndex, this.pt);

  /**
   * Constructs a GeometryLocation specifying a point inside an area geometry.
   *
   * @param component the component of the geometry containing the point
   * @param pt the coordinate of the location
   */
  GeometryLocation.fromPointInArea(Geometry component, Coordinate pt)
      : this(component, INSIDE_AREA, pt);

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
  static bool isWithinDistanceStatic(
      Geometry g0, Geometry g1, double distance) {
    // check envelope distance for a short-circuit negative result
    double envDist =
        g0.getEnvelopeInternal().distance(g1.getEnvelopeInternal());
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
  late List<Geometry?> geom;
  double terminateDistance = 0.0;

  // working
  PointLocator ptLocator = new PointLocator();
  List<GeometryLocation?>? minDistanceLocation;
  double minDistance = double.maxFinite;

  /**
   * Constructs a DistanceOp that computes the distance and nearest points between
   * the two specified geometries.
   * @param g0 a Geometry
   * @param g1 a Geometry
   */
  DistanceOp(Geometry g0, Geometry g1)
      : this.withTerminateDistance(g0, g1, 0.0);

  /**
   * Constructs a DistanceOp that computes the distance and nearest points between
   * the two specified geometries.
   * @param g0 a Geometry
   * @param g1 a Geometry
   * @param terminateDistance the distance on which to terminate the search
   */
  DistanceOp.withTerminateDistance(
      Geometry g0, Geometry g1, double terminateDistance) {
    this.geom = []..length = 2;
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
    if (geom[0] == null || geom[1] == null)
      throw new ArgumentError("null geometries are not supported");
    if (geom[0]!.isEmpty() || geom[1]!.isEmpty()) return 0.0;

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
    List<Coordinate> nearestPts = [
      minDistanceLocation![0]!.getCoordinate(),
      minDistanceLocation![1]!.getCoordinate()
    ];
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
  List<GeometryLocation?>? nearestLocations() {
    computeMinDistance();
    return minDistanceLocation;
  }

  /**
   *
   * @return a pair of {@link GeometryLocation}s for the nearest points
   * @deprecated renamed to nearestLocations
   */
  List<GeometryLocation?>? closestLocations() {
    return nearestLocations();
  }

  void updateMinDistance(List<GeometryLocation?> locGeom, bool flip) {
    // if not set then don't update
    if (locGeom[0] == null) return;

    if (flip) {
      minDistanceLocation![0] = locGeom[1]!;
      minDistanceLocation![1] = locGeom[0]!;
    } else {
      minDistanceLocation![0] = locGeom[0]!;
      minDistanceLocation![1] = locGeom[1]!;
    }
  }

  void computeMinDistance() {
    // only compute once!
    if (minDistanceLocation != null) return;

    minDistanceLocation = []..length = 2;
    computeContainmentDistance();
    if (minDistance <= terminateDistance) return;
    computeFacetDistance();
  }

  void computeContainmentDistance() {
    List<GeometryLocation?> locPtPoly = []..length = 2;
    // test if either geometry has a vertex inside the other
    computeContainmentDistance1(0, locPtPoly);
    if (minDistance <= terminateDistance) return;
    computeContainmentDistance1(1, locPtPoly);
  }

  void computeContainmentDistance1(
      int polyGeomIndex, List<GeometryLocation?> locPtPoly) {
    Geometry? polyGeom = geom[polyGeomIndex];
    // if no polygon then nothing to do
    if (polyGeom!.getDimension() < 2) return;

    int locationsIndex = 1 - polyGeomIndex;
    List polys = PolygonExtracter.getPolygons(polyGeom);
    if (polys.length > 0) {
      List insideLocs =
          ConnectedElementLocationFilter.getLocations(geom[locationsIndex]!);
      computeContainmentDistance2(insideLocs, polys, locPtPoly);
      if (minDistance <= terminateDistance) {
        // this assigment is determined by the order of the args in the computeInside call above
        minDistanceLocation![locationsIndex] = locPtPoly[0];
        minDistanceLocation![polyGeomIndex] = locPtPoly[1];
        return;
      }
    }
  }

  void computeContainmentDistance2(
      List locs, List polys, List<GeometryLocation?> locPtPoly) {
    for (int i = 0; i < locs.length; i++) {
      GeometryLocation loc = locs[i];
      for (int j = 0; j < polys.length; j++) {
        computeContainmentDistance3(loc, polys[j], locPtPoly);
        if (minDistance <= terminateDistance) return;
      }
    }
  }

  void computeContainmentDistance3(
      GeometryLocation ptLoc, Polygon poly, List<GeometryLocation?> locPtPoly) {
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
    List<GeometryLocation?> locGeom = []..length = 2;

    /**
     * Geometries are not wholely inside, so compute distance from lines and points
     * of one to lines and points of the other
     */
    List lines0 = LinearComponentExtracter.getLines(geom[0]!);
    List lines1 = LinearComponentExtracter.getLines(geom[1]!);

    List pts0 = PointExtracter.getPoints(geom[0]!);
    List pts1 = PointExtracter.getPoints(geom[1]!);

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

  void computeMinDistanceLines(
      List lines0, List lines1, List<GeometryLocation?> locGeom) {
    for (int i = 0; i < lines0.length; i++) {
      LineString line0 = lines0[i];
      for (int j = 0; j < lines1.length; j++) {
        LineString line1 = lines1[j];
        computeMinDistanceLineLineList(line0, line1, locGeom);
        if (minDistance <= terminateDistance) return;
      }
    }
  }

  void computeMinDistancePoints(
      List points0, List points1, List<GeometryLocation?> locGeom) {
    for (int i = 0; i < points0.length; i++) {
      Point pt0 = points0[i];
      for (int j = 0; j < points1.length; j++) {
        Point pt1 = points1[j];
        double dist = pt0.getCoordinate()!.distance(pt1.getCoordinate()!);
        if (dist < minDistance) {
          minDistance = dist;
          locGeom[0] = new GeometryLocation(pt0, 0, pt0.getCoordinate()!);
          locGeom[1] = new GeometryLocation(pt1, 0, pt1.getCoordinate()!);
        }
        if (minDistance <= terminateDistance) return;
      }
    }
  }

  void computeMinDistanceLinesPoints(
      List lines, List points, List<GeometryLocation?> locGeom) {
    for (int i = 0; i < lines.length; i++) {
      LineString line = lines[i];
      for (int j = 0; j < points.length; j++) {
        Point pt = points[j];
        computeMinDistanceLinePointList(line, pt, locGeom);
        if (minDistance <= terminateDistance) return;
      }
    }
  }

  void computeMinDistanceLineLineList(
      LineString line0, LineString line1, List<GeometryLocation?> locGeom) {
    if (line0.getEnvelopeInternal().distance(line1.getEnvelopeInternal()) >
        minDistance) return;
    List<Coordinate> coord0 = line0.getCoordinates();
    List<Coordinate> coord1 = line1.getCoordinates();
    // brute force approach!
    for (int i = 0; i < coord0.length - 1; i++) {
      for (int j = 0; j < coord1.length - 1; j++) {
        double dist = Distance.segmentToSegment(
            coord0[i], coord0[i + 1], coord1[j], coord1[j + 1]);
        if (dist < minDistance) {
          minDistance = dist;
          LineSegment seg0 =
              new LineSegment.fromCoordinates(coord0[i], coord0[i + 1]);
          LineSegment seg1 =
              new LineSegment.fromCoordinates(coord1[j], coord1[j + 1]);
          List<Coordinate> closestPt = seg0.closestPoints(seg1);
          locGeom[0] = new GeometryLocation(line0, i, closestPt[0]);
          locGeom[1] = new GeometryLocation(line1, j, closestPt[1]);
        }
        if (minDistance <= terminateDistance) return;
      }
    }
  }

  void computeMinDistanceLinePointList(
      LineString line, Point pt, List<GeometryLocation?> locGeom) {
    if (line.getEnvelopeInternal().distance(pt.getEnvelopeInternal()) >
        minDistance) return;
    List<Coordinate> coord0 = line.getCoordinates();
    Coordinate coord = pt.getCoordinate()!;
    // brute force approach!
    for (int i = 0; i < coord0.length - 1; i++) {
      double dist = Distance.pointToSegment(coord, coord0[i], coord0[i + 1]);
      if (dist < minDistance) {
        minDistance = dist;
        LineSegment seg =
            new LineSegment.fromCoordinates(coord0[i], coord0[i + 1]);
        Coordinate segClosestPoint = seg.closestPoint(coord);
        locGeom[0] = new GeometryLocation(line, i, segClosestPoint);
        locGeom[1] = new GeometryLocation(pt, 0, coord);
      }
      if (minDistance <= terminateDistance) return;
    }
  }
}

/**
 * An algorithm for computing a distance metric
 * which is an approximation to the Hausdorff Distance
 * based on a discretization of the input {@link Geometry}.
 * The algorithm computes the Hausdorff distance restricted to discrete points
 * for one of the geometries.
 * The points can be either the vertices of the geometries (the default),
 * or the geometries with line segments densified by a given fraction.
 * Also determines two points of the Geometries which are separated by the computed distance.
 * <p>
 * This algorithm is an approximation to the standard Hausdorff distance.
 * Specifically,
 * <pre>
 *    for all geometries a, b:    DHD(a, b) &lt;= HD(a, b)
 * </pre>
 * The approximation can be made as close as needed by densifying the input geometries.
 * In the limit, this value will approach the true Hausdorff distance:
 * <pre>
 *    DHD(A, B, densifyFactor) -&gt; HD(A, B) as densifyFactor -&gt; 0.0
 * </pre>
 * The default approximation is exact or close enough for a large subset of useful cases.
 * Examples of these are:
 * <ul>
 * <li>computing distance between Linestrings that are roughly parallel to each other,
 * and roughly equal in length.  This occurs in matching linear networks.
 * <li>Testing similarity of geometries.
 * </ul>
 * An example where the default approximation is not close is:
 * <pre>
 *   A = LINESTRING (0 0, 100 0, 10 100, 10 100)
 *   B = LINESTRING (0 100, 0 10, 80 10)
 *
 *   DHD(A, B) = 22.360679774997898
 *   HD(A, B) ~= 47.8
 * </pre>
 */
class DiscreteHausdorffDistance {
  static double distanceStatic(Geometry g0, Geometry g1) {
    DiscreteHausdorffDistance dist = new DiscreteHausdorffDistance(g0, g1);
    return dist.distance();
  }

  static double distanceStaticDF(Geometry g0, Geometry g1, double densifyFrac) {
    DiscreteHausdorffDistance dist = new DiscreteHausdorffDistance(g0, g1);
    dist.setDensifyFraction(densifyFrac);
    return dist.distance();
  }

  Geometry g0;
  Geometry g1;
  PointPairDistance ptDist = new PointPairDistance();

  /**
   * Value of 0.0 indicates that no densification should take place
   */
  double densifyFrac = 0.0;

  DiscreteHausdorffDistance(this.g0, this.g1);

  /**
   * Sets the fraction by which to densify each segment.
   * Each segment will be (virtually) split into a number of equal-length
   * subsegments, whose fraction of the total length is closest
   * to the given fraction.
   *
   * @param densifyFrac
   */
  void setDensifyFraction(double densifyFrac) {
    if (densifyFrac > 1.0 || densifyFrac <= 0.0)
      throw ArgumentError("Fraction is not in range (0.0 - 1.0]");

    this.densifyFrac = densifyFrac;
  }

  double distance() {
    compute(g0, g1);
    return ptDist.getDistance();
  }

  double orientedDistance() {
    computeOrientedDistance(g0, g1, ptDist);
    return ptDist.getDistance();
  }

  List<Coordinate> getCoordinates() {
    return ptDist.getCoordinates();
  }

  void compute(Geometry g0, Geometry g1) {
    computeOrientedDistance(g0, g1, ptDist);
    computeOrientedDistance(g1, g0, ptDist);
  }

  void computeOrientedDistance(
      Geometry discreteGeom, Geometry geom, PointPairDistance ptDist) {
    MaxPointDistanceFilter distFilter = new MaxPointDistanceFilter(geom);
    discreteGeom.applyCF(distFilter);
    ptDist.setMaximum(distFilter.getMaxPointDistance());

    if (densifyFrac > 0) {
      MaxDensifiedByFractionDistanceFilter fracFilter =
          new MaxDensifiedByFractionDistanceFilter(geom, densifyFrac);
      discreteGeom.applyCSF(fracFilter);
      ptDist.setMaximum(fracFilter.getMaxPointDistance());
    }
  }
}

class MaxDensifiedByFractionDistanceFilter implements CoordinateSequenceFilter {
  PointPairDistance maxPtDist = new PointPairDistance();
  PointPairDistance minPtDist = new PointPairDistance();
  Geometry geom;
  int numSubSegs = 0;

  MaxDensifiedByFractionDistanceFilter(this.geom, double fraction) {
    this.geom = geom;
    numSubSegs = (1.0 / fraction).round(); // TODO this should be Math.rint(a)
  }

  void filter(CoordinateSequence seq, int index) {
    /**
     * This logic also handles skipping Point geometries
     */
    if (index == 0) return;

    Coordinate p0 = seq.getCoordinate(index - 1);
    Coordinate p1 = seq.getCoordinate(index);

    double delx = (p1.x - p0.x) / numSubSegs;
    double dely = (p1.y - p0.y) / numSubSegs;

    for (int i = 0; i < numSubSegs; i++) {
      double x = p0.x + i * delx;
      double y = p0.y + i * dely;
      Coordinate pt = new Coordinate(x, y);
      minPtDist.initialize();
      DistanceToPoint.computeDistance(geom, pt, minPtDist);
      maxPtDist.setMaximum(minPtDist);
    }
  }

  bool isGeometryChanged() {
    return false;
  }

  bool isDone() {
    return false;
  }

  PointPairDistance getMaxPointDistance() {
    return maxPtDist;
  }
}

class MaxPointDistanceFilter implements CoordinateFilter {
  PointPairDistance maxPtDist = new PointPairDistance();
  PointPairDistance minPtDist = new PointPairDistance();
  DistanceToPoint euclideanDist = new DistanceToPoint();
  Geometry geom;

  MaxPointDistanceFilter(this.geom);

  void filter(Coordinate? pt) {
    minPtDist.initialize();
    DistanceToPoint.computeDistance(geom, pt!, minPtDist);
    maxPtDist.setMaximum(minPtDist);
  }

  PointPairDistance getMaxPointDistance() {
    return maxPtDist;
  }
}

/**
 * Computes the Euclidean distance (L2 metric) from a {@link Coordinate} to a {@link Geometry}.
 * Also computes two points on the geometry which are separated by the distance found.
 */
class DistanceToPoint {
  DistanceToPoint() {}

  static void computeDistance(
      Geometry geom, Coordinate pt, PointPairDistance ptDist) {
    if (geom is LineString) {
      computeDistanceL(geom, pt, ptDist);
    } else if (geom is Polygon) {
      computeDistanceP(geom, pt, ptDist);
    } else if (geom is GeometryCollection) {
      for (int i = 0; i < geom.getNumGeometries(); i++) {
        Geometry g = geom.getGeometryN(i);
        computeDistance(g, pt, ptDist);
      }
    } else {
      // assume geom is Point
      ptDist.setMinimum2C(geom.getCoordinate()!, pt);
    }
  }

  static void computeDistanceL(
      LineString line, Coordinate pt, PointPairDistance ptDist) {
    LineSegment tempSegment = new LineSegment.empty();
    var coords = line.getCoordinates();
    for (int i = 0; i < coords.length - 1; i++) {
      tempSegment.setCoordinates(coords[i], coords[i + 1]);
      // this is somewhat inefficient - could do better
      Coordinate closestPt = tempSegment.closestPoint(pt);
      ptDist.setMinimum2C(closestPt, pt);
    }
  }

  static void computeDistanceLS(
      LineSegment segment, Coordinate pt, PointPairDistance ptDist) {
    Coordinate closestPt = segment.closestPoint(pt);
    ptDist.setMinimum2C(closestPt, pt);
  }

  static void computeDistanceP(
      Polygon poly, Coordinate pt, PointPairDistance ptDist) {
    computeDistance(poly.getExteriorRing(), pt, ptDist);
    for (int i = 0; i < poly.getNumInteriorRing(); i++) {
      computeDistance(poly.getInteriorRingN(i), pt, ptDist);
    }
  }
}

/**
 * Computes the distance between the facets (segments and vertices)
 * of two {@link Geometry}s
 * using a Branch-and-Bound algorithm.
 * The Branch-and-Bound algorithm operates over a
 * traversal of R-trees built
 * on the target and the query geometries.
 * <p>
 * This approach provides the following benefits:
 * <ul>
 * <li>Performance is dramatically improved due to the use of the
 * R-tree index
 * and the pruning due to the Branch-and-Bound approach
 * <li>The spatial index on the target geometry is cached
 * which allow reuse in an repeated query situation.
 * </ul>
 * Using this technique is usually much more performant
 * than using the brute-force {@link Geometry#distance(Geometry)}
 * when one or both input geometries are large,
 * or when evaluating many distance computations against
 * a single geometry.
 * <p>
 * This class is thread-safe.
 *
 * @author Martin Davis
 *
 */
class IndexedFacetDistance {
  static final FacetSequenceDistance FACET_SEQ_DIST =
      new FacetSequenceDistance();

  /**
   * Computes the distance between facets of two geometries.
   * <p>
   * For geometries with many segments or points,
   * this can be faster than using a simple distance
   * algorithm.
   *
   * @param g1 a geometry
   * @param g2 a geometry
   * @return the distance between facets of the geometries
   */

  static double distance(Geometry g1, Geometry g2) {
    IndexedFacetDistance dist = new IndexedFacetDistance(g1);
    return dist.facetDistance(g2);
  }

  /**
   * Tests whether the facets of two geometries lie within a given distance.
   *
   * @param g1 a geometry
   * @param g2 a geometry
   * @param distance the distance limit
   * @return true if two facets lie with the given distance
   */

  static bool isWithinDistance(Geometry g1, Geometry g2, double distance) {
    IndexedFacetDistance dist = new IndexedFacetDistance(g1);
    return dist.isWithinDistanceGeometry(g2, distance);
  }

  /**
   * Computes the nearest points of the facets of two geometries.
   *
   * @param g1 a geometry
   * @param g2 a geometry
   * @return the nearest points on the facets of the geometries
   */
  static List<Coordinate>? nearestPoints(Geometry g1, Geometry g2) {
    IndexedFacetDistance dist = new IndexedFacetDistance(g1);
    return dist.nearestPointsToGeometry(g2);
  }

  STRtree? cachedTree;
  Geometry baseGeometry;

  /**
   * Creates a new distance-finding instance for a given target {@link Geometry}.
   * <p>
   * Distances will be computed to all facets of the input geometry.
   * The facets of the geometry are the discrete segments and points
   * contained in its components.
   * In the case of {@link Lineal} and {@link Puntal} inputs,
   * this is equivalent to computing the conventional distance.
   * In the case of {@link Polygonal} inputs, this is equivalent
   * to computing the distance to the polygon boundaries.
   *
   * @param geom a Geometry, which may be of any type.
   */
  IndexedFacetDistance(this.baseGeometry)
      : cachedTree = FacetSequenceTreeBuilder.build(baseGeometry);

  /**
   * Computes the distance from the base geometry to
   * the given geometry.
   *
   * @param g the geometry to compute the distance to
   *
   * @return the computed distance
   */
  double facetDistance(Geometry g) {
    STRtree tree2 = FacetSequenceTreeBuilder.build(g);
    List<Object>? obj =
        cachedTree!.nearestNeighbourWithTree(tree2, FACET_SEQ_DIST);
    FacetSequence fs1 = obj![0] as FacetSequence;
    FacetSequence fs2 = obj[1] as FacetSequence;
    return fs1.distance(fs2);
  }

  /**
   * Computes the nearest locations on the base geometry
   * and the given geometry.
   *
   * @param g the geometry to compute the nearest location to
   * @return the nearest locations
   */
  List<GeometryLocation> nearestLocations(Geometry g) {
    STRtree tree2 = FacetSequenceTreeBuilder.build(g);
    List<Object>? obj =
        cachedTree!.nearestNeighbourWithTree(tree2, FACET_SEQ_DIST);
    FacetSequence fs1 = obj![0] as FacetSequence;
    FacetSequence fs2 = obj[1] as FacetSequence;
    return fs1.nearestLocations(fs2);
  }

  /**
   * Compute the nearest locations on the target geometry
   * and the given geometry.
   *
   * @param g the geometry to compute the nearest point to
   * @return the nearest points
   */
  List<Coordinate>? nearestPointsToGeometry(Geometry g) {
    List<GeometryLocation> minDistanceLocation = nearestLocations(g);
    List<Coordinate>? nearestPts = toPoints(minDistanceLocation);
    return nearestPts;
  }

  static List<Coordinate>? toPoints(List<GeometryLocation>? locations) {
    if (locations == null) return null;
    List<Coordinate> nearestPts = <Coordinate>[
      locations[0].getCoordinate(),
      locations[1].getCoordinate()
    ];
    return nearestPts;
  }

  /**
   * Tests whether the base geometry lies within
   * a specified distance of the given geometry.
   *
   * @param g the geometry to test
   * @param maxDistance the maximum distance to test
   * @return true if the geometry lies with the specified distance
   */
  bool isWithinDistanceGeometry(Geometry g, double maxDistance) {
    // short-ciruit check
    double envDist =
        baseGeometry.getEnvelopeInternal().distance(g.getEnvelopeInternal());
    if (envDist > maxDistance) return false;

    STRtree tree2 = FacetSequenceTreeBuilder.build(g);
    return cachedTree!
        .isWithinDistanceWithTree(tree2, FACET_SEQ_DIST, maxDistance);
  }
}

class FacetSequenceDistance implements ItemDistance {
  double distance(ItemBoundable item1, ItemBoundable item2) {
    FacetSequence fs1 = item1.getItem() as FacetSequence;
    FacetSequence fs2 = item2.getItem() as FacetSequence;
    return fs1.distance(fs2);
  }
}

class FacetSequenceTreeBuilder {
  // 6 seems to be a good facet sequence size
  static final int FACET_SEQUENCE_SIZE = 6;

  // Seems to be better to use a minimum node capacity
  static final int STR_TREE_NODE_CAPACITY = 4;

  static STRtree build(Geometry g) {
    STRtree tree = STRtree.withCapacity(STR_TREE_NODE_CAPACITY);
    List<FacetSequence> sections = computeFacetSequences(g);
    for (var section in sections) {
      tree.insert(section.getEnvelope(), section);
    }
    tree.build();
    return tree;
  }

  /**
   * Creates facet sequences
   *
   * @param g
   * @return List<GeometryFacetSequence>
   */
  static List<FacetSequence> computeFacetSequences(Geometry g) {
    FacetSequenceTreeFilter filter = FacetSequenceTreeFilter();
    g.applyGCF(filter);
    return filter.sections;
  }
}

class FacetSequenceTreeFilter extends GeometryComponentFilter {
  static final int FACET_SEQUENCE_SIZE = 6;
  List<FacetSequence> sections = List.empty(growable: true);
  List<FacetSequence> filterSequence(Geometry geom) {
    filter(geom);
    return sections;
  }

  void filter(Geometry geom) {
    CoordinateSequence seq;
    if (geom is LineString) {
      seq = geom.getCoordinateSequence();
      addFacetSequences(geom, seq, sections);
    } else if (geom is Point) {
      seq = geom.getCoordinateSequence();
      addFacetSequences(geom, seq, sections);
    }
  }

  static void addFacetSequences(
      Geometry geom, CoordinateSequence pts, List sections) {
    int i = 0;
    int size = pts.size();
    while (i <= size - 1) {
      int end = i + FACET_SEQUENCE_SIZE + 1;
      // if only one point remains after this section, include it in this
      // section
      if (end >= size - 1) end = size;
      FacetSequence sect = new FacetSequence(geom, pts, i, end);
      sections.add(sect);
      i = i + FACET_SEQUENCE_SIZE;
    }
  }
}

class FacetSequence {
  Geometry? geom;
  CoordinateSequence pts;
  int? start;
  int? end;

  /**
   * Creates a new sequence of facets based on a {@link CoordinateSequence}
   * contained in the given {@link Geometry}.
   *
   * @param geom the geometry containing the facets
   * @param pts the sequence containing the facet points
   * @param start the index of the start point
   * @param end the index of the end point + 1
   */
  FacetSequence(this.geom, this.pts, this.start, this.end);

  /**
   * Creates a new sequence of facets based on a {@link CoordinateSequence}.
   *
   * @param pts the sequence containing the facet points
   * @param start the index of the start point
   * @param end the index of the end point + 1
   */
  FacetSequence.fromSequence(this.pts, this.start, this.end);

  /**
   * Creates a new sequence for a single point from a {@link CoordinateSequence}.
   *
   * @param pts the sequence containing the facet point
   * @param start the index of the point
   */
  FacetSequence.fromPoint(this.pts, this.start);

  Envelope getEnvelope() {
    Envelope env = new Envelope.empty();
    if (start != null && end != null) {
      for (int? i = start; i! < end!; i++) {
        env.expandToInclude(pts.getX(i), pts.getY(i));
      }
    }
    return env;
  }

  int? size() {
    return end! - start!;
  }

  Coordinate getCoordinate(int index) {
    return pts.getCoordinate(start! + index);
  }

  bool isPoint() {
    return end! - start! == 1;
  }

  /**
   * Computes the distance between this and another
   * <tt>FacetSequence</tt>.
   *
   * @param facetSeq the sequence to compute the distance to
   * @return the minimum distance between the sequences
   */
  double distance(FacetSequence facetSeq) {
    bool isPoint = this.isPoint();
    bool isPointOther = facetSeq.isPoint();
    double distance;

    if (isPoint && isPointOther) {
      Coordinate pt = pts.getCoordinate(start!);
      Coordinate seqPt = facetSeq.pts.getCoordinate(facetSeq.start!);
      distance = pt.distance(seqPt);
    } else if (isPoint) {
      Coordinate pt = pts.getCoordinate(start!);
      distance = computeDistancePointLine(pt, facetSeq, null);
    } else if (isPointOther) {
      Coordinate seqPt = facetSeq.pts.getCoordinate(facetSeq.start!);
      distance = computeDistancePointLine(seqPt, this, null);
    } else {
      distance = computeDistanceLineLine(facetSeq, null);
    }
    return distance;
  }

  /**
   * Computes the locations of the nearest points between this sequence
   * and another sequence.
   * The locations are presented in the same order as the input sequences.
   *
   * @return a pair of {@link GeometryLocation}s for the nearest points
   */
  List<GeometryLocation> nearestLocations(FacetSequence facetSeq) {
    bool isPoint = this.isPoint();
    bool isPointOther = facetSeq.isPoint();
    List<GeometryLocation>? locs = List.empty(growable: true);

    if (isPoint && isPointOther) {
      Coordinate pt = pts.getCoordinate(start!);
      Coordinate seqPt = facetSeq.pts.getCoordinate(facetSeq.start!);
      locs.add(GeometryLocation(geom!, start!, Coordinate.fromCoordinate(pt)));
      locs.add(GeometryLocation(facetSeq.geom!, facetSeq.start!,
          new Coordinate.fromCoordinate(seqPt)));
    } else if (isPoint) {
      Coordinate pt = pts.getCoordinate(start!);
      computeDistancePointLine(pt, facetSeq, locs);
    } else if (isPointOther) {
      Coordinate seqPt = facetSeq.pts.getCoordinate(facetSeq.start!);
      computeDistancePointLine(seqPt, this, locs);
      // unflip the locations
      GeometryLocation tmp = locs[0];
      locs[0] = locs[1];
      locs[1] = tmp;
    } else {
      computeDistanceLineLine(facetSeq, locs);
    }
    return locs;
  }

  double computeDistanceLineLine(
      FacetSequence facetSeq, List<GeometryLocation>? locs) {
    // both linear - compute minimum segment-segment distance
    double minDistance = double.maxFinite;
    if (start != null && end != null) {
      for (int? i = start!; i! < end! - 1; i++) {
        Coordinate p0 = pts.getCoordinate(i);
        Coordinate p1 = pts.getCoordinate(i + 1);
        if (facetSeq.start != null && facetSeq.end != null) {
          for (int j = facetSeq.start!; j < facetSeq.end! - 1; j++) {
            Coordinate q0 = facetSeq.pts.getCoordinate(j);
            Coordinate q1 = facetSeq.pts.getCoordinate(j + 1);

            double dist = Distance.segmentToSegment(p0, p1, q0, q1);
            if (dist < minDistance) {
              minDistance = dist;
              if (locs != null)
                updateNearestLocationsLineLine(
                    i, p0, p1, facetSeq, j, q0, q1, locs);
              if (minDistance <= 0.0) return minDistance;
            }
          }
        }
      }
    }
    return minDistance;
  }

  void updateNearestLocationsLineLine(
      int i,
      Coordinate p0,
      Coordinate p1,
      FacetSequence facetSeq,
      int j,
      Coordinate q0,
      Coordinate q1,
      List<GeometryLocation> locs) {
    LineSegment seg0 = LineSegment.fromCoordinates(p0, p1);
    LineSegment seg1 = new LineSegment.fromCoordinates(q0, q1);
    List<Coordinate> closestPt = seg0.closestPoints(seg1);
    locs[0] =
        GeometryLocation(geom!, i, Coordinate.fromCoordinate(closestPt[0]));
    locs[1] = GeometryLocation(
        facetSeq.geom!, j, new Coordinate.fromCoordinate(closestPt[1]));
  }

  double computeDistancePointLine(
      Coordinate pt, FacetSequence facetSeq, List<GeometryLocation>? locs) {
    double minDistance = double.maxFinite;

    for (int i = facetSeq.start!; i < facetSeq.end! - 1; i++) {
      Coordinate q0 = facetSeq.pts.getCoordinate(i);
      Coordinate q1 = facetSeq.pts.getCoordinate(i + 1);
      double dist = Distance.pointToSegment(pt, q0, q1);
      if (dist < minDistance) {
        minDistance = dist;
        if (locs != null)
          updateNearestLocationsPointLine(pt, facetSeq, i, q0, q1, locs);
        if (minDistance <= 0.0) return minDistance;
      }
    }
    return minDistance;
  }

  void updateNearestLocationsPointLine(Coordinate pt, FacetSequence facetSeq,
      int i, Coordinate q0, Coordinate q1, List<GeometryLocation> locs) {
    locs[0] = GeometryLocation(geom!, start!, Coordinate.fromCoordinate(pt));
    LineSegment seg = LineSegment.fromCoordinates(q0, q1);
    Coordinate segClosestPoint = seg.closestPoint(pt);
    locs[1] = GeometryLocation(
        facetSeq.geom!, i, Coordinate.fromCoordinate(segClosestPoint));
  }

  String toString() {
    StringBuffer buf = StringBuffer();
    buf.write("LINESTRING ( ");
    for (int i = start!; i < end!; i++) {
      if (i > start!) buf.write(", ");
      Coordinate p = pts.getCoordinate(i);
      buf.write(p.x);
      buf.write(' ');
      buf.write(p.y);
    }
    buf.write(' )');
    return buf.toString();
  }
}
