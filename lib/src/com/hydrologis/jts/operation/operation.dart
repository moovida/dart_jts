part of dart_jts;

/**
 * Tests whether a <code>Geometry</code> is simple.
 * In general, the SFS specification of simplicity
 * follows the rule:
 * <ul>
 *    <li> A Geometry is simple if and only if the only self-intersections are at
 *    boundary points.
 * </ul>
 * <p>
 * Simplicity is defined for each {@link Geometry} type as follows:
 * <ul>
 * <li><b>Polygonal</b> geometries are simple by definition, so
 * <code>isSimple</code> trivially returns true.
 * (Note: this means that <tt>isSimple</tt> cannot be used to test
 * for (invalid) self-intersections in <tt>Polygon</tt>s.
 * In order to check if a <tt>Polygonal</tt> geometry has self-intersections,
 * use {@link Geometry#isValid()}).
 * <li><b>Linear</b> geometries are simple iff they do <i>not</i> self-intersect at interior points
 * (i.e. points other than boundary points).
 * This is equivalent to saying that no two linear components satisfy the SFS {@link Geometry#touches(Geometry)}
 * predicate.
 * <li><b>Zero-dimensional (point)</b> geometries are simple if and only if they have no
 * repeated points.
 * <li><b>Empty</b> geometries are <i>always</i> simple, by definition
 * </ul>
 * For {@link Lineal} geometries the evaluation of simplicity
 * can be customized by supplying a {@link BoundaryNodeRule}
 * to define how boundary points are determined.
 * The default is the SFS-standard {@link BoundaryNodeRule#MOD2_BOUNDARY_RULE}.
 * Note that under the <tt>Mod-2</tt> rule, closed <tt>LineString</tt>s (rings)
 * will never satisfy the <tt>touches</tt> predicate at their endpoints, since these are
 * interior points, not boundary points.
 * If it is required to test whether a set of <code>LineString</code>s touch
 * only at their endpoints, use <code>IsSimpleOp</code> with {@link BoundaryNodeRule#ENDPOINT_BOUNDARY_RULE}.
 * For example, this can be used to validate that a set of lines form a topologically valid
 * linear network.
 *
 * @see BoundaryNodeRule
 *
 * @version 1.7
 */
class IsSimpleOp {
  Geometry inputGeom;
  bool isClosedEndpointsInInterior = true;
  Coordinate nonSimpleLocation = null;

  /**
   * Creates a simplicity checker using the default SFS Mod-2 Boundary Node Rule
   *
   * @deprecated use IsSimpleOp(Geometry)
   */
  IsSimpleOp() {}

  /**
   * Creates a simplicity checker using the default SFS Mod-2 Boundary Node Rule
   *
   * @param geom the geometry to test
   */
  IsSimpleOp.withGeom(Geometry geom) {
    this.inputGeom = geom;
  }

  /**
   * Creates a simplicity checker using a given {@link BoundaryNodeRule}
   *
   * @param geom the geometry to test
   * @param boundaryNodeRule the rule to use.
   */
  IsSimpleOp.withGeomAndRule(Geometry geom, BoundaryNodeRule boundaryNodeRule) {
    this.inputGeom = geom;
    isClosedEndpointsInInterior = !boundaryNodeRule.isInBoundary(2);
  }

  /**
   * Tests whether the geometry is simple.
   *
   * @return true if the geometry is simple
   */
  bool isSimple() {
    nonSimpleLocation = null;
    return computeSimple(inputGeom);
  }

  bool computeSimple(Geometry geom) {
    nonSimpleLocation = null;
    if (geom.isEmpty()) return true;
    if (geom is LineString) return isSimpleLinearGeometry(geom);
    if (geom is MultiLineString) return isSimpleLinearGeometry(geom);
    if (geom is MultiPoint) return isSimpleMultiPoint(geom);
    if (geom is Polygonal) return isSimplePolygonal(geom);
    if (geom is GeometryCollection) return isSimpleGeometryCollection(geom);
    // all other geometry types are simple by definition
    return true;
  }

  /**
   * Gets a coordinate for the location where the geometry
   * fails to be simple.
   * (i.e. where it has a non-boundary self-intersection).
   * {@link #isSimple} must be called before this method is called.
   *
   * @return a coordinate for the location of the non-boundary self-intersection
   * or null if the geometry is simple
   */
  Coordinate getNonSimpleLocation() {
    return nonSimpleLocation;
  }

  /**
   * Reports whether a {@link LineString} is simple.
   *
   * @param geom the lineal geometry to test
   * @return true if the geometry is simple
   * @deprecated use isSimple()
   */
  bool isSimpleLine(LineString geom) {
    return isSimpleLinearGeometry(geom);
  }

  /**
   * Reports whether a {@link MultiLineString} geometry is simple.
   *
   * @param geom the lineal geometry to test
   * @return true if the geometry is simple
   * @deprecated use isSimple()
   */
  bool isSimpleMultiLine(MultiLineString geom) {
    return isSimpleLinearGeometry(geom);
  }

  /**
   * A MultiPoint is simple iff it has no repeated points
   * @deprecated use isSimple()
   */
  bool isSimpleMultiPoint(MultiPoint mp) {
    return isSimpleMultiPoint(mp);
  }

  bool isSimpleMultiPoint_(MultiPoint mp) {
    if (mp.isEmpty()) return true;
    Set points = new SplayTreeSet();
    for (int i = 0; i < mp.getNumGeometries(); i++) {
      Point pt = mp.getGeometryN(i);
      Coordinate p = pt.getCoordinate();
      if (points.contains(p)) {
        nonSimpleLocation = p;
        return false;
      }
      points.add(p);
    }
    return true;
  }

  /**
   * Computes simplicity for polygonal geometries.
   * Polygonal geometries are simple if and only if
   * all of their component rings are simple.
   *
   * @param geom a Polygonal geometry
   * @return true if the geometry is simple
   */
  bool isSimplePolygonal(Geometry geom) {
    List rings = LinearComponentExtracter.getLines(geom);
    for (Iterator i = rings.iterator; i.moveNext();) {
      LinearRing ring = i.current as LinearRing;
      if (!isSimpleLinearGeometry(ring)) return false;
    }
    return true;
  }

  /**
   * Semantics for GeometryCollection is
   * simple iff all components are simple.
   *
   * @param geom
   * @return true if the geometry is simple
   */
  bool isSimpleGeometryCollection(Geometry geom) {
    for (int i = 0; i < geom.getNumGeometries(); i++) {
      Geometry comp = geom.getGeometryN(i);
      if (!computeSimple(comp)) return false;
    }
    return true;
  }

  bool isSimpleLinearGeometry(Geometry geom) {
    if (geom.isEmpty()) return true;
    GeometryGraph graph = new GeometryGraph(0, geom);
    LineIntersector li = new RobustLineIntersector();
    SegmentIntersector si = graph.computeSelfNodes(li, true);
    // if no self-intersection, must be simple
    if (!si.hasIntersection()) return true;
    if (si.hasProperIntersection()) {
      nonSimpleLocation = si.getProperIntersectionPoint();
      return false;
    }
    if (hasNonEndpointIntersection(graph)) return false;
    if (isClosedEndpointsInInterior) {
      if (hasClosedEndpointIntersection(graph)) return false;
    }
    return true;
  }

  /**
   * For all edges, check if there are any intersections which are NOT at an endpoint.
   * The Geometry is not simple if there are intersections not at endpoints.
   */
  bool hasNonEndpointIntersection(GeometryGraph graph) {
    for (Iterator i = graph.getEdgeIterator(); i.moveNext();) {
      Edge e = i.current as Edge;
      int maxSegmentIndex = e.getMaximumSegmentIndex();
      for (Iterator eiIt = e.getEdgeIntersectionList().iterator(); eiIt.moveNext();) {
        EdgeIntersection ei = eiIt.current as EdgeIntersection;
        if (!ei.isEndPoint(maxSegmentIndex)) {
          nonSimpleLocation = ei.getCoordinate();
          return true;
        }
      }
    }
    return false;
  }

  /**
   * Tests that no edge intersection is the endpoint of a closed line.
   * This ensures that closed lines are not touched at their endpoint,
   * which is an interior point according to the Mod-2 rule
   * To check this we compute the degree of each endpoint.
   * The degree of endpoints of closed lines
   * must be exactly 2.
   */
  bool hasClosedEndpointIntersection(GeometryGraph graph) {
    Map endPoints = new SplayTreeMap();
    for (Iterator i = graph.getEdgeIterator(); i.moveNext();) {
      Edge e = i.current as Edge;
      bool isClosed = e.isClosed();
      Coordinate p0 = e.getCoordinateWithIndex(0);
      addEndpoint(endPoints, p0, isClosed);
      Coordinate p1 = e.getCoordinateWithIndex(e.getNumPoints() - 1);
      addEndpoint(endPoints, p1, isClosed);
    }

    for (Iterator i = endPoints.values.iterator; i.moveNext();) {
      EndpointInfo eiInfo = i.current as EndpointInfo;
      if (eiInfo.isClosed && eiInfo.degree != 2) {
        nonSimpleLocation = eiInfo.getCoordinate();
        return true;
      }
    }
    return false;
  }

  /**
   * Add an endpoint to the map, creating an entry for it if none exists
   */
  void addEndpoint(Map endPoints, Coordinate p, bool isClosed) {
    EndpointInfo eiInfo = endPoints[p] as EndpointInfo;
    if (eiInfo == null) {
      eiInfo = new EndpointInfo(p);
      endPoints[p] = eiInfo;
    }
    eiInfo.addEndpoint(isClosed);
  }
}

class EndpointInfo {
  Coordinate pt;
  bool isClosed;
  int degree;

  EndpointInfo(Coordinate pt) {
    this.pt = pt;
    isClosed = false;
    degree = 0;
  }

  Coordinate getCoordinate() {
    return pt;
  }

  void addEndpoint(bool isClosed) {
    degree++;
    this.isClosed |= isClosed;
  }
}

/**
 * Implements the algorithms required to compute the <code>isValid()</code> method
 * for {@link Geometry}s.
 * See the documentation for the various geometry types for a specification of validity.
 *
 * @version 1.7
 */
class IsValidOp {
  /**
   * Tests whether a {@link Geometry} is valid.
   * @param geom the Geometry to test
   * @return true if the geometry is valid
   */
  static bool isValidStatic(Geometry geom) {
    IsValidOp isValidOp = new IsValidOp(geom);
    return isValidOp.isValid();
  }

  /**
   * Checks whether a coordinate is valid for processing.
   * Coordinates are valid iff their x and y ordinates are in the
   * range of the floating point representation.
   *
   * @param coord the coordinate to validate
   * @return <code>true</code> if the coordinate is valid
   */
  static bool isValidStaticCoord(Coordinate coord) {
    if (coord.x.isNaN) return false;
    if (coord.x.isInfinite) return false;
    if (coord.y.isNaN) return false;
    if (coord.y.isInfinite) return false;
    return true;
  }

  /**
   * Find a point from the list of testCoords
   * that is NOT a node in the edge for the list of searchCoords
   *
   * @return the point found, or <code>null</code> if none found
   */
  static Coordinate findPtNotNode(List<Coordinate> testCoords, LinearRing searchRing, GeometryGraph graph) {
    // find edge corresponding to searchRing.
    Edge searchEdge = graph.findEdgeFromLine(searchRing);
    // find a point in the testCoords which is not a node of the searchRing
    EdgeIntersectionList eiList = searchEdge.getEdgeIntersectionList();
    // somewhat inefficient - is there a better way? (Use a node map, for instance?)
    for (int i = 0; i < testCoords.length; i++) {
      Coordinate pt = testCoords[i];
      if (!eiList.isIntersection(pt)) return pt;
    }
    return null;
  }

  Geometry parentGeometry; // the base Geometry to be validated
  /**
   * If the following condition is TRUE JTS will validate inverted shells and exverted holes
   * (the ESRI SDE model)
   */
  bool isSelfTouchingRingFormingHoleValid = false;
  TopologyValidationError validErr;

  IsValidOp(Geometry parentGeometry) {
    this.parentGeometry = parentGeometry;
  }

  /**
   * Sets whether polygons using <b>Self-Touching Rings</b> to form
   * holes are reported as valid.
   * If this flag is set, the following Self-Touching conditions
   * are treated as being valid:
   * <ul>
   * <li>the shell ring self-touches to create a hole touching the shell
   * <li>a hole ring self-touches to create two holes touching at a point
   * </ul>
   * <p>
   * The default (following the OGC SFS standard)
   * is that this condition is <b>not</b> valid (<code>false</code>).
   * <p>
   * This does not affect whether Self-Touching Rings
   * disconnecting the polygon interior are considered valid
   * (these are considered to be <b>invalid</b> under the SFS, and many other
   * spatial models as well).
   * This includes "bow-tie" shells,
   * which self-touch at a single point causing the interior to
   * be disconnected,
   * and "C-shaped" holes which self-touch at a single point causing an island to be formed.
   *
   * @param isValid states whether geometry with this condition is valid
   */
  void setSelfTouchingRingFormingHoleValid(bool isValid) {
    isSelfTouchingRingFormingHoleValid = isValid;
  }

  /**
   * Computes the validity of the geometry,
   * and returns <tt>true</tt> if it is valid.
   *
   * @return true if the geometry is valid
   */
  bool isValid() {
    checkValid(parentGeometry);
    return validErr == null;
  }

  /**
   * Computes the validity of the geometry,
   * and if not valid returns the validation error for the geometry,
   * or null if the geometry is valid.
   *
   * @return the validation error, if the geometry is invalid
   * or null if the geometry is valid
   */
  TopologyValidationError getValidationError() {
    checkValid(parentGeometry);
    return validErr;
  }

  void checkValid(Geometry g) {
    validErr = null;

    // empty geometries are always valid!
    if (g.isEmpty()) return;

    if (g is Point)
      checkValid(g);
    else if (g is MultiPoint)
      checkValid(g);
    // LineString also handles LinearRings
    else if (g is LinearRing)
      checkValid(g);
    else if (g is LineString)
      checkValid(g);
    else if (g is Polygon)
      checkValid(g);
    else if (g is MultiPolygon)
      checkValid(g);
    else if (g is GeometryCollection)
      checkValid(g);
    else
      throw new UnsupportedError(g.runtimeType.toString());
  }

  /**
   * Checks validity of a Point.
   */
  void checkValidP(Point g) {
    checkInvalidCoordinatesList(g.getCoordinates());
  }

  /**
   * Checks validity of a MultiPoint.
   */
  void checkValidMP(MultiPoint g) {
    checkInvalidCoordinatesList(g.getCoordinates());
  }

  /**
   * Checks validity of a LineString.  Almost anything goes for linestrings!
   */
  void checkValidL(LineString g) {
    checkInvalidCoordinatesList(g.getCoordinates());
    if (validErr != null) return;
    GeometryGraph graph = new GeometryGraph(0, g);
    checkTooFewPoints(graph);
  }

  /**
   * Checks validity of a LinearRing.
   */
  void checkValidLR(LinearRing g) {
    checkInvalidCoordinatesList(g.getCoordinates());
    if (validErr != null) return;
    checkClosedRing(g);
    if (validErr != null) return;

    GeometryGraph graph = new GeometryGraph(0, g);
    checkTooFewPoints(graph);
    if (validErr != null) return;

    LineIntersector li = new RobustLineIntersector();
    graph.computeSelfNodes3(li, true, true);
    checkNoSelfIntersectingRings(graph);
  }

  /**
   * Checks the validity of a polygon.
   * Sets the validErr flag.
   */
  void checkValidPol(Polygon g) {
    checkInvalidCoordinates(g);
    if (validErr != null) return;
    checkClosedRings(g);
    if (validErr != null) return;

    GeometryGraph graph = new GeometryGraph(0, g);

    checkTooFewPoints(graph);
    if (validErr != null) return;
    checkConsistentArea(graph);
    if (validErr != null) return;

    if (!isSelfTouchingRingFormingHoleValid) {
      checkNoSelfIntersectingRings(graph);
      if (validErr != null) return;
    }
    checkHolesInShell(g, graph);
    if (validErr != null) return;
    //SLOWcheckHolesNotNested(g);
    checkHolesNotNested(g, graph);
    if (validErr != null) return;
    checkConnectedInteriors(graph);
  }

  void checkValidMPol(MultiPolygon g) {
    for (int i = 0; i < g.getNumGeometries(); i++) {
      Polygon p = g.getGeometryN(i) as Polygon;
      checkInvalidCoordinates(p);
      if (validErr != null) return;
      checkClosedRings(p);
      if (validErr != null) return;
    }

    GeometryGraph graph = new GeometryGraph(0, g);

    checkTooFewPoints(graph);
    if (validErr != null) return;
    checkConsistentArea(graph);
    if (validErr != null) return;
    if (!isSelfTouchingRingFormingHoleValid) {
      checkNoSelfIntersectingRings(graph);
      if (validErr != null) return;
    }
    for (int i = 0; i < g.getNumGeometries(); i++) {
      Polygon p = g.getGeometryN(i) as Polygon;
      checkHolesInShell(p, graph);
      if (validErr != null) return;
    }
    for (int i = 0; i < g.getNumGeometries(); i++) {
      Polygon p = g.getGeometryN(i) as Polygon;
      checkHolesNotNested(p, graph);
      if (validErr != null) return;
    }
    checkShellsNotNested(g, graph);
    if (validErr != null) return;
    checkConnectedInteriors(graph);
  }

  void checkValidGC(GeometryCollection gc) {
    for (int i = 0; i < gc.getNumGeometries(); i++) {
      Geometry g = gc.getGeometryN(i);
      checkValid(g);
      if (validErr != null) return;
    }
  }

  void checkInvalidCoordinatesList(List<Coordinate> coords) {
    for (int i = 0; i < coords.length; i++) {
      if (!isValidStaticCoord(coords[i])) {
        validErr = new TopologyValidationError.withCoordinate(TopologyValidationError.INVALID_COORDINATE, coords[i]);
        return;
      }
    }
  }

  void checkInvalidCoordinates(Polygon poly) {
    checkInvalidCoordinatesList(poly.getExteriorRing().getCoordinates());
    if (validErr != null) return;
    for (int i = 0; i < poly.getNumInteriorRing(); i++) {
      checkInvalidCoordinatesList(poly.getInteriorRingN(i).getCoordinates());
      if (validErr != null) return;
    }
  }

  void checkClosedRings(Polygon poly) {
    checkClosedRing(poly.getExteriorRing());
    if (validErr != null) return;
    for (int i = 0; i < poly.getNumInteriorRing(); i++) {
      checkClosedRing(poly.getInteriorRingN(i));
      if (validErr != null) return;
    }
  }

  void checkClosedRing(LinearRing ring) {
    if (ring.isEmpty()) return;
    if (!ring.isClosed()) {
      Coordinate pt = null;
      if (ring.getNumPoints() >= 1) pt = ring.getCoordinateN(0);
      validErr = new TopologyValidationError.withCoordinate(TopologyValidationError.RING_NOT_CLOSED, pt);
    }
  }

  void checkTooFewPoints(GeometryGraph graph) {
    if (graph.hasTooFewPoints()) {
      validErr = new TopologyValidationError.withCoordinate(TopologyValidationError.TOO_FEW_POINTS, graph.getInvalidPoint());
      return;
    }
  }

  /**
   * Checks that the arrangement of edges in a polygonal geometry graph
   * forms a consistent area.
   *
   * @param graph
   *
   * @see ConsistentAreaTester
   */
  void checkConsistentArea(GeometryGraph graph) {
    ConsistentAreaTester cat = new ConsistentAreaTester(graph);
    bool isValidArea = cat.isNodeConsistentArea();
    if (!isValidArea) {
      validErr = new TopologyValidationError.withCoordinate(TopologyValidationError.SELF_INTERSECTION, cat.getInvalidPoint());
      return;
    }
    if (cat.hasDuplicateRings()) {
      validErr = new TopologyValidationError.withCoordinate(TopologyValidationError.DUPLICATE_RINGS, cat.getInvalidPoint());
    }
  }

  /**
   * Check that there is no ring which self-intersects (except of course at its endpoints).
   * This is required by OGC topology rules (but not by other models
   * such as ESRI SDE, which allow inverted shells and exverted holes).
   *
   * @param graph the topology graph of the geometry
   */
  void checkNoSelfIntersectingRings(GeometryGraph graph) {
    for (Iterator i = graph.getEdgeIterator(); i.moveNext();) {
      Edge e = i.current as Edge;
      checkNoSelfIntersectingRing(e.getEdgeIntersectionList());
      if (validErr != null) return;
    }
  }

  /**
   * Check that a ring does not self-intersect, except at its endpoints.
   * Algorithm is to count the number of times each node along edge occurs.
   * If any occur more than once, that must be a self-intersection.
   */
  void checkNoSelfIntersectingRing(EdgeIntersectionList eiList) {
    Set nodeSet = new SplayTreeSet();
    bool isFirst = true;
    for (Iterator i = eiList.iterator(); i.moveNext();) {
      EdgeIntersection ei = i.current as EdgeIntersection;
      if (isFirst) {
        isFirst = false;
        continue;
      }
      if (nodeSet.contains(ei.coord)) {
        validErr = new TopologyValidationError.withCoordinate(TopologyValidationError.RING_SELF_INTERSECTION, ei.coord);
        return;
      } else {
        nodeSet.add(ei.coord);
      }
    }
  }

  /**
   * Tests that each hole is inside the polygon shell.
   * This routine assumes that the holes have previously been tested
   * to ensure that all vertices lie on the shell or on the same side of it
   * (i.e. that the hole rings do not cross the shell ring).
   * In other words, this test is only correct if the ConsistentArea test is passed first.
   * Given this, a simple point-in-polygon test of a single point in the hole can be used,
   * provided the point is chosen such that it does not lie on the shell.
   *
   * @param p the polygon to be tested for hole inclusion
   * @param graph a GeometryGraph incorporating the polygon
   */
  void checkHolesInShell(Polygon p, GeometryGraph graph) {
    LinearRing shell = p.getExteriorRing();
    bool isShellEmpty = shell.isEmpty();
    //PointInRing pir = new SimplePointInRing(shell);
    //PointInRing pir = new SIRtreePointInRing(shell);
    //PointInRing pir = new MCPointInRing(shell);
    PointOnGeometryLocator pir = new IndexedPointInAreaLocator(shell);

    for (int i = 0; i < p.getNumInteriorRing(); i++) {
      LinearRing hole = p.getInteriorRingN(i);
      Coordinate holePt = null;
      if (hole.isEmpty()) continue;
      holePt = findPtNotNode(hole.getCoordinates(), shell, graph);
      /**
       * If no non-node hole vertex can be found, the hole must
       * split the polygon into disconnected interiors.
       * This will be caught by a subsequent check.
       */
      if (holePt == null) return;

      bool outside = isShellEmpty || (Location.EXTERIOR == pir.locate(holePt));
      if (outside) {
        validErr = new TopologyValidationError.withCoordinate(TopologyValidationError.HOLE_OUTSIDE_SHELL, holePt);
        return;
      }
    }
  }

  /**
   * Tests that no hole is nested inside another hole.
   * This routine assumes that the holes are disjoint.
   * To ensure this, holes have previously been tested
   * to ensure that:
   * <ul>
   * <li>they do not partially overlap
   *      (checked by <code>checkRelateConsistency</code>)
   * <li>they are not identical
   *      (checked by <code>checkRelateConsistency</code>)
   * </ul>
   */
  void checkHolesNotNested(Polygon p, GeometryGraph graph) {
    IndexedNestedRingTester nestedTester = new IndexedNestedRingTester(graph);
    //SimpleNestedRingTester nestedTester = new SimpleNestedRingTester(arg[0]);
    //SweeplineNestedRingTester nestedTester = new SweeplineNestedRingTester(arg[0]);

    for (int i = 0; i < p.getNumInteriorRing(); i++) {
      LinearRing innerHole = p.getInteriorRingN(i);
      if (innerHole.isEmpty()) continue;
      nestedTester.add(innerHole);
    }
    bool isNonNested = nestedTester.isNonNested();
    if (!isNonNested) {
      validErr = new TopologyValidationError.withCoordinate(TopologyValidationError.NESTED_HOLES, nestedTester.getNestedPoint());
    }
  }

  /**
   * Tests that no element polygon is wholly in the interior of another element polygon.
   * <p>
   * Preconditions:
   * <ul>
   * <li>shells do not partially overlap
   * <li>shells do not touch along an edge
   * <li>no duplicate rings exist
   * </ul>
   * This routine relies on the fact that while polygon shells may touch at one or
   * more vertices, they cannot touch at ALL vertices.
   */
  void checkShellsNotNested(MultiPolygon mp, GeometryGraph graph) {
    for (int i = 0; i < mp.getNumGeometries(); i++) {
      Polygon p = mp.getGeometryN(i) as Polygon;
      LinearRing shell = p.getExteriorRing();
      for (int j = 0; j < mp.getNumGeometries(); j++) {
        if (i == j) continue;
        Polygon p2 = mp.getGeometryN(j) as Polygon;
        checkShellNotNested(shell, p2, graph);
        if (validErr != null) return;
      }
    }
  }

  /**
   * Check if a shell is incorrectly nested within a polygon.  This is the case
   * if the shell is inside the polygon shell, but not inside a polygon hole.
   * (If the shell is inside a polygon hole, the nesting is valid.)
   * <p>
   * The algorithm used relies on the fact that the rings must be properly contained.
   * E.g. they cannot partially overlap (this has been previously checked by
   * <code>checkRelateConsistency</code> )
   */
  void checkShellNotNested(LinearRing shell, Polygon p, GeometryGraph graph) {
    List<Coordinate> shellPts = shell.getCoordinates();
    // test if shell is inside polygon shell
    LinearRing polyShell = p.getExteriorRing();
    if (polyShell.isEmpty()) return;
    List<Coordinate> polyPts = polyShell.getCoordinates();
    Coordinate shellPt = findPtNotNode(shellPts, polyShell, graph);
    // if no point could be found, we can assume that the shell is outside the polygon
    if (shellPt == null) return;
    bool insidePolyShell = PointLocation.isInRing(shellPt, polyPts);
    if (!insidePolyShell) return;

    // if no holes, this is an error!
    if (p.getNumInteriorRing() <= 0) {
      validErr = new TopologyValidationError.withCoordinate(TopologyValidationError.NESTED_SHELLS, shellPt);
      return;
    }

    /**
     * Check if the shell is inside one of the holes.
     * This is the case if one of the calls to checkShellInsideHole
     * returns a null coordinate.
     * Otherwise, the shell is not properly contained in a hole, which is an error.
     */
    Coordinate badNestedPt = null;
    for (int i = 0; i < p.getNumInteriorRing(); i++) {
      LinearRing hole = p.getInteriorRingN(i);
      badNestedPt = checkShellInsideHole(shell, hole, graph);
      if (badNestedPt == null) return;
    }
    validErr = new TopologyValidationError.withCoordinate(TopologyValidationError.NESTED_SHELLS, badNestedPt);
  }

  /**
   * This routine checks to see if a shell is properly contained in a hole.
   * It assumes that the edges of the shell and hole do not
   * properly intersect.
   *
   * @return <code>null</code> if the shell is properly contained, or
   *   a Coordinate which is not inside the hole if it is not
   *
   */
  Coordinate checkShellInsideHole(LinearRing shell, LinearRing hole, GeometryGraph graph) {
    List<Coordinate> shellPts = shell.getCoordinates();
    List<Coordinate> holePts = hole.getCoordinates();
    // TODO: improve performance of this - by sorting pointlists for instance?
    Coordinate shellPt = findPtNotNode(shellPts, hole, graph);
    // if point is on shell but not hole, check that the shell is inside the hole
    if (shellPt != null) {
      bool insideHole = PointLocation.isInRing(shellPt, holePts);
      if (!insideHole) {
        return shellPt;
      }
    }
    Coordinate holePt = findPtNotNode(holePts, shell, graph);
    // if point is on hole but not shell, check that the hole is outside the shell
    if (holePt != null) {
      bool insideShell = PointLocation.isInRing(holePt, shellPts);
      if (insideShell) {
        return holePt;
      }
      return null;
    }
    assert(true, "points in shell and hole appear to be equal");
    return null;
  }

  void checkConnectedInteriors(GeometryGraph graph) {
    ConnectedInteriorTester cit = new ConnectedInteriorTester(graph);
    if (!cit.isInteriorsConnected()) validErr = new TopologyValidationError.withCoordinate(TopologyValidationError.DISCONNECTED_INTERIOR, cit.getCoordinate());
  }
}

/**
 * Contains information about the nature and location of a {@link Geometry}
 * validation error
 *
 * @version 1.7
 */
class TopologyValidationError {
  /**
   * Not used
   * @deprecated
   */
  static final int ERROR = 0;

  /**
   * No longer used - repeated points are considered valid as per the SFS
   * @deprecated
   */
  static final int REPEATED_POINT = 1;

  /**
   * Indicates that a hole of a polygon lies partially or completely in the exterior of the shell
   */
  static final int HOLE_OUTSIDE_SHELL = 2;

  /**
   * Indicates that a hole lies in the interior of another hole in the same polygon
   */
  static final int NESTED_HOLES = 3;

  /**
   * Indicates that the interior of a polygon is disjoint
   * (often caused by set of contiguous holes splitting the polygon into two parts)
   */
  static final int DISCONNECTED_INTERIOR = 4;

  /**
   * Indicates that two rings of a polygonal geometry intersect
   */
  static final int SELF_INTERSECTION = 5;

  /**
   * Indicates that a ring self-intersects
   */
  static final int RING_SELF_INTERSECTION = 6;

  /**
   * Indicates that a polygon component of a MultiPolygon lies inside another polygonal component
   */
  static final int NESTED_SHELLS = 7;

  /**
   * Indicates that a polygonal geometry contains two rings which are identical
   */
  static final int DUPLICATE_RINGS = 8;

  /**
   * Indicates that either
   * <ul>
   * <li>a LineString contains a single point
   * <li>a LinearRing contains 2 or 3 points
   * </ul>
   */
  static final int TOO_FEW_POINTS = 9;

  /**
   * Indicates that the <code>X</code> or <code>Y</code> ordinate of
   * a Coordinate is not a valid numeric value (e.g. {@link Double#NaN} )
   */
  static final int INVALID_COORDINATE = 10;

  /**
   * Indicates that a ring is not correctly closed
   * (the first and the last coordinate are different)
   */
  static final int RING_NOT_CLOSED = 11;

  /**
   * Messages corresponding to error codes
   */
  static final List<String> errMsg = [
    "Topology Validation Error",
    "Repeated Point",
    "Hole lies outside shell",
    "Holes are nested",
    "Interior is disconnected",
    "Self-intersection",
    "Ring Self-intersection",
    "Nested shells",
    "Duplicate Rings",
    "Too few distinct points in geometry component",
    "Invalid Coordinate",
    "Ring is not closed"
  ];

  int errorType;
  Coordinate pt;

  /**
   * Creates a validation error with the given type and location
   *
   * @param errorType the type of the error
   * @param pt the location of the error
   */
  TopologyValidationError.withCoordinate(int errorType, Coordinate pt) {
    this.errorType = errorType;
    if (pt != null) this.pt = pt.copy();
  }

  /**
   * Creates a validation error of the given type with a null location
   *
   * @param errorType the type of the error
   *
   */
  TopologyValidationError(int errorType) : this.withCoordinate(errorType, null);

  /**
   * Returns the location of this error (on the {@link Geometry} containing the error).
   *
   * @return a {@link Coordinate} on the input geometry
   */
  Coordinate getCoordinate() {
    return pt;
  }

  /**
   * Gets the type of this error.
   *
   * @return the error type
   */
  int getErrorType() {
    return errorType;
  }

  /**
   * Gets an error message describing this error.
   * The error message does not describe the location of the error.
   *
   * @return the error message
   */
  String getMessage() {
    return errMsg[errorType];
  }

  /**
   * Gets a message describing the type and location of this error.
   * @return the error message
   */
  String toString() {
    String locStr = "";
    if (pt != null) locStr = " at or near point $pt";
    return getMessage() + locStr;
  }
}

/**
 * Tests whether any of a set of {@link LinearRing}s are
 * nested inside another ring in the set, using a spatial
 * index to speed up the comparisons.
 *
 * @version 1.7
 */
class IndexedNestedRingTester {
  GeometryGraph graph; // used to find non-node vertices
  List rings = [];
  Envelope totalEnv = new Envelope.empty();
  SpatialIndex index;
  Coordinate nestedPt;

  IndexedNestedRingTester(GeometryGraph graph) {
    this.graph = graph;
  }

  Coordinate getNestedPoint() {
    return nestedPt;
  }

  void add(LinearRing ring) {
    rings.add(ring);
    totalEnv.expandToIncludeEnvelope(ring.getEnvelopeInternal());
  }

  bool isNonNested() {
    buildIndex();

    for (int i = 0; i < rings.length; i++) {
      LinearRing innerRing = rings[i] as LinearRing;
      List<Coordinate> innerRingPts = innerRing.getCoordinates();

      List results = index.query(innerRing.getEnvelopeInternal());
//System.out.println(results.size());
      for (int j = 0; j < results.length; j++) {
        LinearRing searchRing = results[j] as LinearRing;
        List<Coordinate> searchRingPts = searchRing.getCoordinates();

        if (innerRing == searchRing) continue;

        if (!innerRing.getEnvelopeInternal().intersectsEnvelope(searchRing.getEnvelopeInternal())) continue;

        Coordinate innerRingPt = IsValidOp.findPtNotNode(innerRingPts, searchRing, graph);

        /**
         * If no non-node pts can be found, this means
         * that the searchRing touches ALL of the innerRing vertices.
         * This indicates an invalid polygon, since either
         * the two holes create a disconnected interior,
         * or they touch in an infinite number of points
         * (i.e. along a line segment).
         * Both of these cases are caught by other tests,
         * so it is safe to simply skip this situation here.
         */
        if (innerRingPt == null) continue;

        bool isInside = PointLocation.isInRing(innerRingPt, searchRingPts);
        if (isInside) {
          nestedPt = innerRingPt;
          return false;
        }
      }
    }
    return true;
  }

  void buildIndex() {
    index = new STRtree();

    for (int i = 0; i < rings.length; i++) {
      LinearRing ring = rings[i] as LinearRing;
      Envelope env = ring.getEnvelopeInternal();
      index.insert(env, ring);
    }
  }
}
