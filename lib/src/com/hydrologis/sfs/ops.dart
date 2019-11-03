part of dart_sfs;
/*
 * Copyright (c) 2016 Vivid Solutions.
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * and Eclipse Distribution License v. 1.0 which accompanies this distribution.
 * The Eclipse Public License is available at http://www.eclipse.org/legal/epl-v10.html
 * and the Eclipse Distribution License is available at
 *
 * http://www.eclipse.org/org/documents/edl-v10.php.
 */

/// An interface for rules which determine whether node points
/// which are in boundaries of {@link Lineal} geometry components
/// are in the boundary of the parent geometry collection.
/// The SFS specifies a single kind of boundary node rule,
/// the {@link Mod2BoundaryNodeRule} rule.
/// However, other kinds of Boundary Node Rules are appropriate
/// in specific situations (for instance, linear network topology
/// usually follows the {@link EndPointBoundaryNodeRule}.)
/// Some JTS operations
/// (such as {@link RelateOp}, {@link BoundaryOp} and {@link IsSimpleOp})
/// allow the BoundaryNodeRule to be specified,
/// and respect the supplied rule when computing the results of the operation.
/// <p>
/// An example use case for a non-SFS-standard Boundary Node Rule is
/// that of checking that a set of {@link LineString}s have
/// valid linear network topology, when turn-arounds are represented
/// as closed rings.  In this situation, the entry road to the
/// turn-around is only valid when it touches the turn-around ring
/// at the single (common) endpoint.  This is equivalent
/// to requiring the set of <tt>LineString</tt>s to be
/// <b>simple</b> under the {@link EndPointBoundaryNodeRule}.
/// The SFS-standard {@link Mod2BoundaryNodeRule} is not
/// sufficient to perform this test, since it
/// states that closed rings have <b>no</b> boundary points.
/// <p>
/// This interface and its subclasses follow the <tt>Strategy</tt> design pattern.
///
/// @author Martin Davis
/// @version 1.7
///
/// @see RelateOp
/// @see BoundaryOp
/// @see IsSimpleOp
/// @see PointLocator
abstract class BoundaryNodeRule {

  /// Tests whether a point that lies in <tt>boundaryCount</tt>
  /// geometry component boundaries is considered to form part of the boundary
  /// of the parent geometry.
  ///
  /// @param boundaryCount the number of component boundaries that this point occurs in
  /// @return true if points in this number of boundaries lie in the parent boundary
  bool isInBoundary(int boundaryCount);

  /// The Mod-2 Boundary Node Rule (which is the rule specified in the OGC SFS).
  /// @see Mod2BoundaryNodeRule
  static final BoundaryNodeRule MOD2_BOUNDARY_RULE = new Mod2BoundaryNodeRule();

  /// The Endpoint Boundary Node Rule.
  /// @see EndPointBoundaryNodeRule
  static final BoundaryNodeRule ENDPOINT_BOUNDARY_RULE = new EndPointBoundaryNodeRule();

  /// The MultiValent Endpoint Boundary Node Rule.
  /// @see MultiValentEndPointBoundaryNodeRule
  static final BoundaryNodeRule MULTIVALENT_ENDPOINT_BOUNDARY_RULE = new MultiValentEndPointBoundaryNodeRule();

  /// The Monovalent Endpoint Boundary Node Rule.
  /// @see MonoValentEndPointBoundaryNodeRule
  static final BoundaryNodeRule MONOVALENT_ENDPOINT_BOUNDARY_RULE = new MonoValentEndPointBoundaryNodeRule();

  /// The Boundary Node Rule specified by the OGC Simple Features Specification,
  /// which is the same as the Mod-2 rule.
  /// @see Mod2BoundaryNodeRule
  static final BoundaryNodeRule OGC_SFS_BOUNDARY_RULE = MOD2_BOUNDARY_RULE;

}


/// A {@link BoundaryNodeRule} specifies that points are in the
/// boundary of a lineal geometry iff
/// the point lies on the boundary of an odd number
/// of components.
/// Under this rule {@link LinearRing}s and closed
/// {@link LineString}s have an empty boundary.
/// <p>
/// This is the rule specified by the <i>OGC SFS</i>,
/// and is the default rule used in JTS.
///
/// @author Martin Davis
/// @version 1.7
class Mod2BoundaryNodeRule
    implements BoundaryNodeRule {
  bool isInBoundary(int boundaryCount) {
// the "Mod-2 Rule"
    return boundaryCount % 2 == 1;
  }
}

/// A {@link BoundaryNodeRule} which specifies that any points which are endpoints
/// of lineal components are in the boundary of the
/// parent geometry.
/// This corresponds to the "intuitive" topological definition
/// of boundary.
/// Under this rule {@link LinearRing}s have a non-empty boundary
/// (the common endpoint of the underlying LineString).
/// <p>
/// This rule is useful when dealing with linear networks.
/// For example, it can be used to check
/// whether linear networks are correctly noded.
/// The usual network topology constraint is that linear segments may touch only at endpoints.
/// In the case of a segment touching a closed segment (ring) at one point,
/// the Mod2 rule cannot distinguish between the permitted case of touching at the
/// node point and the invalid case of touching at some other interior (non-node) point.
/// The EndPoint rule does distinguish between these cases,
/// so is more appropriate for use.
///
/// @author Martin Davis
/// @version 1.7
class EndPointBoundaryNodeRule
    implements BoundaryNodeRule {
  bool isInBoundary(int boundaryCount) {
    return boundaryCount > 0;
  }
}

/// A {@link BoundaryNodeRule} which determines that only
/// endpoints with valency greater than 1 are on the boundary.
/// This corresponds to the boundary of a {@link MultiLineString}
/// being all the "attached" endpoints, but not
/// the "unattached" ones.
///
/// @author Martin Davis
/// @version 1.7
class MultiValentEndPointBoundaryNodeRule
    implements BoundaryNodeRule {
  bool isInBoundary(int boundaryCount) {
    return boundaryCount > 1;
  }
}

/// A {@link BoundaryNodeRule} which determines that only
/// endpoints with valency of exactly 1 are on the boundary.
/// This corresponds to the boundary of a {@link MultiLineString}
/// being all the "unattached" endpoints.
///
/// @author Martin Davis
/// @version 1.7
class MonoValentEndPointBoundaryNodeRule
    implements BoundaryNodeRule {
  bool isInBoundary(int boundaryCount) {
    return boundaryCount == 1;
  }
}


/// Tests whether a <code>Geometry</code> is simple.
/// In general, the SFS specification of simplicity
/// follows the rule:
/// <ul>
///    <li> A Geometry is simple if and only if the only self-intersections are at
///    boundary points.
/// </ul>
/// <p>
/// Simplicity is defined for each {@link Geometry} type as follows:
/// <ul>
/// <li><b>Polygonal</b> geometries are simple by definition, so
/// <code>isSimple</code> trivially returns true.
/// (Note: this means that <tt>isSimple</tt> cannot be used to test
/// for (invalid) self-intersections in <tt>Polygon</tt>s.
/// In order to check if a <tt>Polygonal</tt> geometry has self-intersections,
/// use {@link Geometry#isValid()}).
/// <li><b>Linear</b> geometries are simple iff they do <i>not</i> self-intersect at interior points
/// (i.e. points other than boundary points).
/// This is equivalent to saying that no two linear components satisfy the SFS {@link Geometry#touches(Geometry)}
/// predicate.
/// <li><b>Zero-dimensional (point)</b> geometries are simple if and only if they have no
/// repeated points.
/// <li><b>Empty</b> geometries are <i>always</i> simple, by definition
/// </ul>
/// For {@link Lineal} geometries the evaluation of simplicity
/// can be customized by supplying a {@link BoundaryNodeRule}
/// to define how boundary points are determined.
/// The default is the SFS-standard {@link BoundaryNodeRule#MOD2_BOUNDARY_RULE}.
/// Note that under the <tt>Mod-2</tt> rule, closed <tt>LineString</tt>s (rings)
/// will never satisfy the <tt>touches</tt> predicate at their endpoints, since these are
/// interior points, not boundary points.
/// If it is required to test whether a set of <code>LineString</code>s touch
/// only at their endpoints, use <code>IsSimpleOp</code> with {@link BoundaryNodeRule#ENDPOINT_BOUNDARY_RULE}.
/// For example, this can be used to validate that a set of lines form a topologically valid
/// linear network.
///
/// @see BoundaryNodeRule
///
/// @version 1.7
//class IsSimpleOp {
//  Geometry inputGeom;
//  bool isClosedEndpointsInInterior = true;
//  Coordinate nonSimpleLocation = null;
//
//  /// Creates a simplicity checker using the default SFS Mod-2 Boundary Node Rule
//  ///
//  /// @param geom the geometry to test
//  IsSimpleOp(Geometry geom) {
//    this.inputGeom = geom;
//  }
//
//  /// Creates a simplicity checker using a given {@link BoundaryNodeRule}
//  ///
//  /// @param geom the geometry to test
//  /// @param boundaryNodeRule the rule to use.
//  IsSimpleOp.withRule(Geometry geom, BoundaryNodeRule boundaryNodeRule) {
//    this.inputGeom = geom;
//    isClosedEndpointsInInterior = !boundaryNodeRule.isInBoundary(2);
//  }
//
//  /// Tests whether the geometry is simple.
//  ///
//  /// @return true if the geometry is simple
//  bool isSimple() {
//    nonSimpleLocation = null;
//    return computeSimple(inputGeom);
//  }
//
//  bool computeSimple(Geometry geom) {
//    nonSimpleLocation = null;
//    if (geom.isEmpty) return true;
//    if (geom is LineString) return isSimpleLinearGeometry(geom);
//    if (geom is MultiLineString) return isSimpleLinearGeometry(geom);
//    if (geom is MultiPoint) {
//      return isSimpleMultiPoint(geom);
//    }
//    if (geom is Polygon) return isSimplePolygonal(geom);
//    if (geom is MultiPolygon) return isSimplePolygonal(geom);
//    if (geom is GeometryCollection) {
//      return isSimpleGeometryCollection(geom);
//    }
//    // all other geometry types are simple by definition
//    return true;
//  }
//
//  /// Gets a coordinate for the location where the geometry
//  /// fails to be simple.
//  /// (i.e. where it has a non-boundary self-intersection).
//  /// {@link #isSimple} must be called before this method is called.
//  ///
//  /// @return a coordinate for the location of the non-boundary self-intersection
//  /// or null if the geometry is simple
//  Coordinate
//  getNonSimpleLocation() {
//    return nonSimpleLocation;
//  }
//
//  /// A MultiPoint is simple iff it has no repeated points
//  bool  isSimpleMultiPoint(MultiPoint mp) {
//    if (mp.isEmpty) return true;
//    Set points = SplayTreeSet();
//    for (int i = 0; i < mp.getNumGeometries(); i++) {
//      Point pt = mp.getGeometryN(i);
//      Coordinate p = pt.getCoordinate();
//      if (points.contains(p)) {
//        nonSimpleLocation = p;
//        return false;
//      }
//      points.add(p);
//    }
//    return true;
//  }
//
//  /// Computes simplicity for polygonal geometries.
//  /// Polygonal geometries are simple if and only if
//  /// all of their component rings are simple.
//  ///
//  /// @param geom a Polygonal geometry
//  /// @return true if the geometry is simple
//  bool  isSimplePolygonal(Geometry geom) {
//    List<LineString> rings = LinearComponentExtracter.getLines(geom);
//    for (int i = 0; i < rings.length; i++) {
//      if (!isSimpleLinearGeometry(rings[i])) {
//        return false;
//      }
//    }
//    return true;
//  }
//
//  /// Semantics for GeometryCollection is
//  /// simple iff all components are simple.
//  ///
//  /// @param geom
//  /// @return true if the geometry is simple
//  bool
//  isSimpleGeometryCollection(Geometry geom) {
//    for (int i = 0; i < geom.getNumGeometries(); i++) {
//      Geometry comp = geom.getGeometryN(i);
//      if (!computeSimple(comp))
//        return false;
//    }
//    return true;
//  }
//
//  /**
//   * Reports whether a linear geometry is simple.
//   *
//   * @param geom the lineal geometry to test
//   * @return true if the geometry is simple
//   */
//  bool   isSimpleLinearGeometry(Geometry geom) {
//    if (geom.isEmpty) return true;
//    GeometryGraph graph = new GeometryGraph(0, geom);
//    LineIntersector li = new RobustLineIntersector();
//    SegmentIntersector si = graph.computeSelfNodes(li, true);
//    // if no self-intersection, must be simple
//    if (!si.hasIntersection()) return true;
//    if (si.hasProperIntersection()) {
//      nonSimpleLocation = si.getProperIntersectionPoint();
//      return false;
//    }
//    if (hasNonEndpointIntersection(graph)) return false;
//    if (isClosedEndpointsInInterior) {
//      if (hasClosedEndpointIntersection(graph)) return false;
//    }
//    return true;
//  }
//
//  /// For all edges, check if there are any intersections which are NOT at an endpoint.
//  /// The Geometry is not simple if there are intersections not at endpoints.
//  bool  hasNonEndpointIntersection(GeometryGraph graph) {
//    for (Iterator i = graph.getEdgeIterator(); i.hasNext();) {
//      Edge e = (Edge) i.next();
//      int maxSegmentIndex = e.getMaximumSegmentIndex();
//      for (Iterator eiIt = e.getEdgeIntersectionList().iterator(); eiIt
//          .hasNext();) {
//        EdgeIntersection ei = (EdgeIntersection) eiIt.next();
//        if (!ei.isEndPoint(maxSegmentIndex)) {
//          nonSimpleLocation = ei.getCoordinate();
//          return true;
//        }
//      }
//    }
//    return false;
//  }
//
//
//  /**
//   * Tests that no edge intersection is the endpoint of a closed line.
//   * This ensures that closed lines are not touched at their endpoint,
//   * which is an interior point according to the Mod-2 rule
//   * To check this we compute the degree of each endpoint.
//   * The degree of endpoints of closed lines
//   * must be exactly 2.
//   */
//  bool hasClosedEndpointIntersection(GeometryGraph graph)
//  {
//  Map endPoints = new TreeMap();
//  for (Iterator i = graph.getEdgeIterator(); i.hasNext(); ) {
//  Edge e = (Edge) i.next();
//  int maxSegmentIndex = e.getMaximumSegmentIndex();
//  bool isClosed = e.isClosed();
//  Coordinate p0 = e.getCoordinate(0);
//  addEndpoint(endPoints, p0, isClosed);
//  Coordinate p1 = e.getCoordinate(e.getNumPoints() - 1);
//  addEndpoint(endPoints, p1, isClosed);
//  }
//
//  for (Iterator i = endPoints.values().iterator(); i.hasNext(); ) {
//  EndpointInfo eiInfo = (EndpointInfo) i.next();
//  if (eiInfo.isClosed && eiInfo.degree != 2) {
//  nonSimpleLocation = eiInfo.getCoordinate();
//  return true;
//  }
//  }
//  return false;
//  }
//
//  /**
//   * Add an endpoint to the map, creating an entry for it if none exists
//   */
//  void addEndpoint(Map endPoints, Coordinate p, bool isClosed)
//  {
//  EndpointInfo eiInfo = (EndpointInfo) endPoints.get(p);
//  if (eiInfo == null) {
//  eiInfo = new EndpointInfo(p);
//  endPoints.put(p, eiInfo);
//  }
//  eiInfo.addEndpoint(isClosed);
//  }
//
//}

class EndpointInfo {
  Coordinate pt;
  bool isClosed;
  int degree;

   EndpointInfo(Coordinate pt)
  {
    this.pt = pt;
    isClosed = false;
    degree = 0;
  }

   Coordinate getCoordinate() { return pt; }

   void addEndpoint(bool isClosed)
  {
    degree++;
    this.isClosed |= isClosed;
  }
}


/// Computes the boundary of a {@link Geometry}.
/// This operation will always return a {@link Geometry} of the appropriate
/// dimension for the boundary (even if the input geometry is empty).
/// The boundary of zero-dimensional geometries (Points) is
/// always the empty {@link GeometryCollection}.
///
/// @author Martin Davis
/// @author Andrea Antonello - dart port
/// @version 1.7

class BoundaryOp {
  static Geometry getBoundaryFromGeometry(Geometry g) {
    BoundaryOp bop = BoundaryOp(g);
    return bop.getBoundary();
  }

  Geometry geom;
  GeometryFactory geomFact = GeometryFactory.defaultPrecision();
  Mod2BoundaryNodeRule bnRule = Mod2BoundaryNodeRule();

  BoundaryOp(Geometry geom) {
    this.geom = geom;
    this.bnRule = BoundaryNodeRule.MOD2_BOUNDARY_RULE;
  }

  BoundaryOp.withRule(Geometry geom, BoundaryNodeRule bnRule)
  {
    this.geom = geom;
    this.bnRule = bnRule;
  }

  Geometry getBoundary() {
    if (geom is LineString) return boundaryLineString(geom);
    if (geom is MultiLineString) {
      return boundaryMultiLineString(geom);
    }
    return geom.getBoundary();
  }

  MultiPoint getEmptyMultiPoint() {
    return geomFact.createMultiPoint();
  }

  Geometry boundaryMultiLineString(MultiLineString mLine) {
    if (geom.isEmpty()) {
      return getEmptyMultiPoint();
    }

    List<Coordinate> bdyPts = computeBoundaryCoordinates(mLine);

    // return Point or MultiPoint
    if (bdyPts.length == 1) {
      return geomFact.createPoint(bdyPts[0]);
    }
    // this handles 0 points case as well
    return geomFact.createMultiPointFromCoords(bdyPts);
  }

/*
// MD - superseded
  private Coordinate[] computeBoundaryFromGeometryGraph(MultiLineString mLine)
  {
    GeometryGraph g = new GeometryGraph(0, mLine, bnRule);
    Coordinate[] bdyPts = g.getBoundaryPoints();
    return bdyPts;
  }
*/

  Map<Coordinate, Counter> endpointMap;

  List<Coordinate> computeBoundaryCoordinates(MultiLineString mLine) {
    List<Coordinate> bdyPts = [];
    endpointMap = SplayTreeMap();
    for (int i = 0; i < mLine.getNumGeometries(); i++) {
      LineString line = mLine.getGeometryN(i);
      if (line.isEmpty) {
        continue;
      }
      addEndpoint(line.getCoordinateN(0));
      addEndpoint(line.getCoordinateN(line.length - 1));
    }

    endpointMap.forEach((coord, counter) {
      int valence = counter.count;
      if (bnRule.isInBoundary(valence)) {
        bdyPts.add(coord);
      }
    });

    return bdyPts;
  }

  void addEndpoint(Coordinate pt) {
    Counter counter = endpointMap[pt];
    if (counter == null) {
      counter = Counter();
      endpointMap[pt] = counter;
    }
    counter.count++;
  }

  Geometry boundaryLineString(LineString line) {
    if (geom.isEmpty()) {
      return getEmptyMultiPoint();
    }

    if (line.isClosed) {
      // check whether endpoints of valence 2 are on the boundary or not
      bool closedEndpointOnBoundary = bnRule.isInBoundary(2);
      if (closedEndpointOnBoundary) {
        return line.getStartPoint();
      } else {
        return geomFact.createMultiPoint();
      }
    }
    return geomFact
        .createMultiPoint([line.getStartPoint(), line.getEndPoint()]);
  }
}

/// Stores an integer count, for use as a Map entry.
///
/// @author Martin Davis
/// @version 1.7
class Counter {
  /// The value of the count
  int count;
}
