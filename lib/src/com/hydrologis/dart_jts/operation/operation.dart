part of dart_jts;

/**
 * Computes the topological relationship between two Geometries.
 * <p>
 * RelateComputer does not need to build a complete graph structure to compute
 * the IntersectionMatrix.  The relationship between the geometries can
 * be computed by simply examining the labelling of edges incident on each node.
 * <p>
 * RelateComputer does not currently support arbitrary GeometryCollections.
 * This is because GeometryCollections can contain overlapping Polygons.
 * In order to correct compute relate on overlapping Polygons, they
 * would first need to be noded and merged (if not explicitly, at least
 * implicitly).
 *
 * @version 1.7
 */
class RelateComputer {
  LineIntersector li = new RobustLineIntersector();
  PointLocator ptLocator = new PointLocator();
  List<GeometryGraph> arg; // the arg(s) of the operation
  NodeMap nodes = new NodeMap(new RelateNodeFactory());

  // this intersection matrix will hold the results compute for the relate
  IntersectionMatrix? im = null;
  List isolatedEdges = [];

  // the intersection point found (if any)
  Coordinate? invalidPoint;

  RelateComputer(this.arg);

  IntersectionMatrix computeIM() {
    IntersectionMatrix im = new IntersectionMatrix();
    // since Geometries are finite and embedded in a 2-D space, the EE element must always be 2
    im.set(Location.EXTERIOR, Location.EXTERIOR, 2);

    // if the Geometries don't overlap there is nothing to do
    if (!arg[0]
        .getGeometry()!
        .getEnvelopeInternal()
        .intersectsEnvelope(arg[1].getGeometry()!.getEnvelopeInternal())) {
      computeDisjointIM(im);
      return im;
    }
    arg[0].computeSelfNodes(li, false);
    arg[1].computeSelfNodes(li, false);

    // compute intersections between edges of the two input geometries
    SegmentIntersector intersector =
        arg[0].computeEdgeIntersections(arg[1], li, false);
//System.out.println("computeIM: # segment intersection tests: " + intersector.numTests);
    computeIntersectionNodes(0);
    computeIntersectionNodes(1);
    /**
     * Copy the labelling for the nodes in the parent Geometries.  These override
     * any labels determined by intersections between the geometries.
     */
    copyNodesAndLabels(0);
    copyNodesAndLabels(1);

    // complete the labelling for any nodes which only have a label for a single geometry
//Debug.addWatch(nodes.find(new Coordinate(110, 200)));
//Debug.printWatch();
    labelIsolatedNodes();
//Debug.printWatch();

    // If a proper intersection was found, we can set a lower bound on the IM.
    computeProperIntersectionIM(intersector, im);

    /**
     * Now process improper intersections
     * (eg where one or other of the geometries has a vertex at the intersection point)
     * We need to compute the edge graph at all nodes to determine the IM.
     */

    // build EdgeEnds for all intersections
    EdgeEndBuilder eeBuilder = new EdgeEndBuilder();
    List ee0 = eeBuilder.computeEdgeEnds(arg[0].getEdgeIterator());
    insertEdgeEnds(ee0);
    List ee1 = eeBuilder.computeEdgeEnds(arg[1].getEdgeIterator());
    insertEdgeEnds(ee1);

//Debug.println("==== NodeList ===");
//Debug.print(nodes);

    labelNodeEdges();

    /**
     * Compute the labeling for isolated components
     * <br>
     * Isolated components are components that do not touch any other components in the graph.
     * They can be identified by the fact that they will
     * contain labels containing ONLY a single element, the one for their parent geometry.
     * We only need to check components contained in the input graphs, since
     * isolated components will not have been replaced by new components formed by intersections.
     */
//debugPrintln("Graph A isolated edges - ");
    labelIsolatedEdges(0, 1);
//debugPrintln("Graph B isolated edges - ");
    labelIsolatedEdges(1, 0);

    // update the IM from all components
    updateIM(im);
    return im;
  }

  void insertEdgeEnds(List ee) {
    for (Iterator i = ee.iterator; i.moveNext();) {
      EdgeEnd e = i.current as EdgeEnd;
      nodes.add(e);
    }
  }

  void computeProperIntersectionIM(
      SegmentIntersector intersector, IntersectionMatrix im) {
    // If a proper intersection is found, we can set a lower bound on the IM.
    int dimA = arg[0].getGeometry()!.getDimension();
    int dimB = arg[1].getGeometry()!.getDimension();
    bool hasProper = intersector.hasProperIntersection();
    bool hasProperInterior = intersector.hasProperInteriorIntersection();

    // For Geometry's of dim 0 there can never be proper intersections.

    /**
     * If edge segments of Areas properly intersect, the areas must properly overlap.
     */
    if (dimA == 2 && dimB == 2) {
      if (hasProper) im.setAtLeastDimensionSymbols("212101212");
    }
    /**
     * If an Line segment properly intersects an edge segment of an Area,
     * it follows that the Interior of the Line intersects the Boundary of the Area.
     * If the intersection is a proper <i>interior</i> intersection, then
     * there is an Interior-Interior intersection too.
     * Note that it does not follow that the Interior of the Line intersects the Exterior
     * of the Area, since there may be another Area component which contains the rest of the Line.
     */
    else if (dimA == 2 && dimB == 1) {
      if (hasProper) im.setAtLeastDimensionSymbols("FFF0FFFF2");
      if (hasProperInterior) im.setAtLeastDimensionSymbols("1FFFFF1FF");
    } else if (dimA == 1 && dimB == 2) {
      if (hasProper) im.setAtLeastDimensionSymbols("F0FFFFFF2");
      if (hasProperInterior) im.setAtLeastDimensionSymbols("1F1FFFFFF");
    }
    /* If edges of LineStrings properly intersect *in an interior point*, all
        we can deduce is that
        the interiors intersect.  (We can NOT deduce that the exteriors intersect,
        since some other segments in the geometries might cover the points in the
        neighbourhood of the intersection.)
        It is important that the point be known to be an interior point of
        both Geometries, since it is possible in a self-intersecting geometry to
        have a proper intersection on one segment that is also a boundary point of another segment.
    */
    else if (dimA == 1 && dimB == 1) {
      if (hasProperInterior) im.setAtLeastDimensionSymbols("0FFFFFFFF");
    }
  }

  /**
   * Copy all nodes from an arg geometry into this graph.
   * The node label in the arg geometry overrides any previously computed
   * label for that argIndex.
   * (E.g. a node may be an intersection node with
   * a computed label of BOUNDARY,
   * but in the original arg Geometry it is actually
   * in the interior due to the Boundary Determination Rule)
   */
  void copyNodesAndLabels(int argIndex) {
    for (Iterator i = arg[argIndex].getNodeIterator(); i.moveNext();) {
      Node graphNode = i.current as Node;
      Node newNode = nodes.addNodeFromCoordinate(graphNode.getCoordinate());
      newNode.setLabelWithIndex(
          argIndex, graphNode.getLabel()!.getLocation(argIndex));
//node.print(System.out);
    }
  }

  /**
   * Insert nodes for all intersections on the edges of a Geometry.
   * Label the created nodes the same as the edge label if they do not already have a label.
   * This allows nodes created by either self-intersections or
   * mutual intersections to be labelled.
   * Endpoint nodes will already be labelled from when they were inserted.
   */
  void computeIntersectionNodes(int argIndex) {
    for (Iterator i = arg[argIndex].getEdgeIterator(); i.moveNext();) {
      Edge e = i.current as Edge;
      int eLoc = e.getLabel()!.getLocation(argIndex);
      for (Iterator eiIt = e.getEdgeIntersectionList().iterator();
          eiIt.moveNext();) {
        EdgeIntersection ei = eiIt.current as EdgeIntersection;
        RelateNode n = nodes.addNodeFromCoordinate(ei.coord) as RelateNode;
        if (eLoc == Location.BOUNDARY)
          n.setLabelBoundary(argIndex);
        else {
          if (n.getLabel()!.isNull(argIndex))
            n.setLabelWithIndex(argIndex, Location.INTERIOR);
        }
//Debug.println(n);
      }
    }
  }

  /**
   * For all intersections on the edges of a Geometry,
   * label the corresponding node IF it doesn't already have a label.
   * This allows nodes created by either self-intersections or
   * mutual intersections to be labelled.
   * Endpoint nodes will already be labelled from when they were inserted.
   */
  void labelIntersectionNodes(int argIndex) {
    for (Iterator i = arg[argIndex].getEdgeIterator(); i.moveNext();) {
      Edge e = i.current as Edge;
      int eLoc = e.getLabel()!.getLocation(argIndex);
      for (Iterator eiIt = e.getEdgeIntersectionList().iterator();
          eiIt.moveNext();) {
        EdgeIntersection ei = eiIt.current as EdgeIntersection;
        RelateNode n = nodes.find(ei.coord) as RelateNode;
        if (n.getLabel()!.isNull(argIndex)) {
          if (eLoc == Location.BOUNDARY)
            n.setLabelBoundary(argIndex);
          else
            n.setLabelWithIndex(argIndex, Location.INTERIOR);
        }
//n.print(System.out);
      }
    }
  }

  /**
   * If the Geometries are disjoint, we need to enter their dimension and
   * boundary dimension in the Ext rows in the IM
   */
  void computeDisjointIM(IntersectionMatrix im) {
    Geometry ga = arg[0].getGeometry()!;
    if (!ga.isEmpty()) {
      im.set(Location.INTERIOR, Location.EXTERIOR, ga.getDimension());
      im.set(Location.BOUNDARY, Location.EXTERIOR, ga.getBoundaryDimension());
    }
    Geometry gb = arg[1].getGeometry()!;
    if (!gb.isEmpty()) {
      im.set(Location.EXTERIOR, Location.INTERIOR, gb.getDimension());
      im.set(Location.EXTERIOR, Location.BOUNDARY, gb.getBoundaryDimension());
    }
  }

  void labelNodeEdges() {
    for (Iterator ni = nodes.iterator(); ni.moveNext();) {
      RelateNode node = ni.current as RelateNode;
      node.getEdges()!.computeLabelling(arg);
//Debug.print(node.getEdges());
//node.print(System.out);
    }
  }

  /**
   * update the IM with the sum of the IMs for each component
   */
  void updateIM(IntersectionMatrix im) {
//Debug.println(im);
    for (Iterator ei = isolatedEdges.iterator; ei.moveNext();) {
      Edge e = ei.current as Edge;
      e.updateIM(im);
//Debug.println(im);
    }
    for (Iterator ni = nodes.iterator(); ni.moveNext();) {
      RelateNode node = ni.current as RelateNode;
      node.updateIM(im);
//Debug.println(im);
      node.updateIMFromEdges(im);
//Debug.println(im);
//node.print(System.out);
    }
  }

  /**
   * Processes isolated edges by computing their labelling and adding them
   * to the isolated edges list.
   * Isolated edges are guaranteed not to touch the boundary of the target (since if they
   * did, they would have caused an intersection to be computed and hence would
   * not be isolated)
   */
  void labelIsolatedEdges(int thisIndex, int targetIndex) {
    for (Iterator ei = arg[thisIndex].getEdgeIterator(); ei.moveNext();) {
      Edge e = ei.current as Edge;
      if (e.isIsolated()) {
        labelIsolatedEdge(e, targetIndex, arg[targetIndex].getGeometry()!);
        isolatedEdges.add(e);
      }
    }
  }

  /**
   * Label an isolated edge of a graph with its relationship to the target geometry.
   * If the target has dim 2 or 1, the edge can either be in the interior or the exterior.
   * If the target has dim 0, the edge must be in the exterior
   */
  void labelIsolatedEdge(Edge e, int targetIndex, Geometry target) {
    // this won't work for GeometryCollections with both dim 2 and 1 geoms
    if (target.getDimension() > 0) {
      // since edge is not in boundary, may not need the full generality of PointLocator?
      // Possibly should use ptInArea locator instead?  We probably know here
      // that the edge does not touch the bdy of the target Geometry
      int loc = ptLocator.locate(e.getCoordinate()!, target);
      e.getLabel()!.setAllLocations(targetIndex, loc);
    } else {
      e.getLabel()!.setAllLocations(targetIndex, Location.EXTERIOR);
    }
//System.out.println(e.getLabel());
  }

  /**
   * Isolated nodes are nodes whose labels are incomplete
   * (e.g. the location for one Geometry is null).
   * This is the case because nodes in one graph which don't intersect
   * nodes in the other are not completely labelled by the initial process
   * of adding nodes to the nodeList.
   * To complete the labelling we need to check for nodes that lie in the
   * interior of edges, and in the interior of areas.
   */
  void labelIsolatedNodes() {
    for (Iterator ni = nodes.iterator(); ni.moveNext();) {
      Node n = ni.current;
      Label label = n.getLabel()!;
      // isolated nodes should always have at least one geometry in their label
      Assert.isTrue(
          label.getGeometryCount() > 0, "node with empty label found");
      if (n.isIsolated()) {
        if (label.isNull(0))
          labelIsolatedNode(n, 0);
        else
          labelIsolatedNode(n, 1);
      }
    }
  }

  /**
   * Label an isolated node with its relationship to the target geometry.
   */
  void labelIsolatedNode(Node n, int targetIndex) {
    int loc =
        ptLocator.locate(n.getCoordinate(), arg[targetIndex].getGeometry()!);
    n.getLabel()!.setAllLocations(targetIndex, loc);
//debugPrintln(n.getLabel());
  }
}

/**
 * The base class for operations that require {@link GeometryGraph}s.
 *
 * @version 1.7
 */
class GeometryGraphOperation {
  final LineIntersector li = new RobustLineIntersector();
  PrecisionModel? resultPrecisionModel;

  /**
   * The operation args into an array so they can be accessed by index
   */
  late List<GeometryGraph> arg; // the arg(s) of the operation

  GeometryGraphOperation(Geometry g0, Geometry g1)
      : this.withRule(g0, g1, BoundaryNodeRule.OGC_SFS_BOUNDARY_RULE
//         BoundaryNodeRule.ENDPOINT_BOUNDARY_RULE
            );

  GeometryGraphOperation.withRule(
      Geometry g0, Geometry g1, BoundaryNodeRule boundaryNodeRule) {
    // use the most precise model for the result
    if (g0.getPrecisionModel().compareTo(g1.getPrecisionModel()) >= 0)
      setComputationPrecision(g0.getPrecisionModel());
    else
      setComputationPrecision(g1.getPrecisionModel());

    arg = []; //..length = (2);
    arg.add(new GeometryGraph.args3(0, g0, boundaryNodeRule));
    arg.add(new GeometryGraph.args3(1, g1, boundaryNodeRule));
    // arg[0] = new GeometryGraph.args3(0, g0, boundaryNodeRule);
    // arg[1] = new GeometryGraph.args3(1, g1, boundaryNodeRule);
  }

  GeometryGraphOperation.singleGeom(Geometry g0) {
    setComputationPrecision(g0.getPrecisionModel());

    arg = []; //..length = (1);
    arg.add(new GeometryGraph(0, g0));
    // arg[0] = new GeometryGraph(0, g0);
    ;
  }

  Geometry getArgGeometry(int i) {
    return arg[i].getGeometry()!;
  }

  void setComputationPrecision(PrecisionModel pm) {
    resultPrecisionModel = pm;
    li.setPrecisionModel(resultPrecisionModel);
  }
}

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
  Coordinate? nonSimpleLocation = null;

  /**
   * Creates a simplicity checker using the default SFS Mod-2 Boundary Node Rule
   *
   * @param geom the geometry to test
   */
  IsSimpleOp.withGeom(this.inputGeom);

  /**
   * Creates a simplicity checker using a given {@link BoundaryNodeRule}
   *
   * @param geom the geometry to test
   * @param boundaryNodeRule the rule to use.
   */
  IsSimpleOp.withGeomAndRule(
      this.inputGeom, BoundaryNodeRule boundaryNodeRule) {
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
  Coordinate? getNonSimpleLocation() {
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
      Point pt = mp.getGeometryN(i) as Point;
      Coordinate p = pt.getCoordinate()!;
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
      for (Iterator eiIt = e.getEdgeIntersectionList().iterator();
          eiIt.moveNext();) {
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
  late bool isClosed;
  late int degree;

  EndpointInfo(this.pt) {
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
  static Coordinate? findPtNotNode(
      List<Coordinate> testCoords, LinearRing searchRing, GeometryGraph graph) {
    // find edge corresponding to searchRing.
    Edge searchEdge = graph.findEdgeFromLine(searchRing)!;
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
  TopologyValidationError? validErr;

  IsValidOp(this.parentGeometry);

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
  TopologyValidationError? getValidationError() {
    checkValid(parentGeometry);
    return validErr;
  }

  void checkValid(Geometry g) {
    validErr = null;

    // empty geometries are always valid!
    if (g.isEmpty()) return;

    if (g is Point)
      checkValidP(g);
    else if (g is MultiPoint)
      checkValidMP(g);
    // LineString also handles LinearRings
    else if (g is LinearRing)
      checkValidLR(g);
    else if (g is LineString)
      checkValidL(g);
    else if (g is Polygon)
      checkValidPol(g);
    else if (g is MultiPolygon)
      checkValidMPol(g);
    else if (g is GeometryCollection)
      checkValidGC(g);
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
        validErr = new TopologyValidationError.withCoordinate(
            TopologyValidationError.INVALID_COORDINATE, coords[i]);
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
      Coordinate? pt = null;
      if (ring.getNumPoints() >= 1) pt = ring.getCoordinateN(0);
      validErr = new TopologyValidationError.withCoordinate(
          TopologyValidationError.RING_NOT_CLOSED, pt);
    }
  }

  void checkTooFewPoints(GeometryGraph graph) {
    if (graph.hasTooFewPoints()) {
      validErr = new TopologyValidationError.withCoordinate(
          TopologyValidationError.TOO_FEW_POINTS, graph.getInvalidPoint());
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
      validErr = new TopologyValidationError.withCoordinate(
          TopologyValidationError.SELF_INTERSECTION, cat.getInvalidPoint());
      return;
    }
    if (cat.hasDuplicateRings()) {
      validErr = new TopologyValidationError.withCoordinate(
          TopologyValidationError.DUPLICATE_RINGS, cat.getInvalidPoint());
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
        validErr = new TopologyValidationError.withCoordinate(
            TopologyValidationError.RING_SELF_INTERSECTION, ei.coord);
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
      Coordinate? holePt = null;
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
        validErr = new TopologyValidationError.withCoordinate(
            TopologyValidationError.HOLE_OUTSIDE_SHELL, holePt);
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
      validErr = new TopologyValidationError.withCoordinate(
          TopologyValidationError.NESTED_HOLES, nestedTester.getNestedPoint());
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
    Coordinate? shellPt = findPtNotNode(shellPts, polyShell, graph);
    // if no point could be found, we can assume that the shell is outside the polygon
    if (shellPt == null) return;
    bool insidePolyShell = PointLocation.isInRing(shellPt, polyPts);
    if (!insidePolyShell) return;

    // if no holes, this is an error!
    if (p.getNumInteriorRing() <= 0) {
      validErr = new TopologyValidationError.withCoordinate(
          TopologyValidationError.NESTED_SHELLS, shellPt);
      return;
    }

    /**
     * Check if the shell is inside one of the holes.
     * This is the case if one of the calls to checkShellInsideHole
     * returns a null coordinate.
     * Otherwise, the shell is not properly contained in a hole, which is an error.
     */
    Coordinate? badNestedPt = null;
    for (int i = 0; i < p.getNumInteriorRing(); i++) {
      LinearRing hole = p.getInteriorRingN(i);
      badNestedPt = checkShellInsideHole(shell, hole, graph);
      if (badNestedPt == null) return;
    }
    validErr = new TopologyValidationError.withCoordinate(
        TopologyValidationError.NESTED_SHELLS, badNestedPt);
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
  Coordinate? checkShellInsideHole(
      LinearRing shell, LinearRing hole, GeometryGraph graph) {
    List<Coordinate> shellPts = shell.getCoordinates();
    List<Coordinate> holePts = hole.getCoordinates();
    // TODO: improve performance of this - by sorting pointlists for instance?
    Coordinate? shellPt = findPtNotNode(shellPts, hole, graph);
    // if point is on shell but not hole, check that the shell is inside the hole
    if (shellPt != null) {
      bool insideHole = PointLocation.isInRing(shellPt, holePts);
      if (!insideHole) {
        return shellPt;
      }
    }
    Coordinate? holePt = findPtNotNode(holePts, shell, graph);
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
    if (!cit.isInteriorsConnected())
      validErr = new TopologyValidationError.withCoordinate(
          TopologyValidationError.DISCONNECTED_INTERIOR, cit.getCoordinate());
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

  int errorType = 0;
  Coordinate? pt;

  /**
   * Creates a validation error with the given type and location
   *
   * @param errorType the type of the error
   * @param pt the location of the error
   */
  TopologyValidationError.withCoordinate(int errorType, Coordinate? pt) {
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
  Coordinate? getCoordinate() {
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
  SpatialIndex? index;
  Coordinate? nestedPt;

  IndexedNestedRingTester(this.graph);

  Coordinate? getNestedPoint() {
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

      List results = index!.query(innerRing.getEnvelopeInternal());
//System.out.println(results.size());
      for (int j = 0; j < results.length; j++) {
        LinearRing searchRing = results[j] as LinearRing;
        List<Coordinate> searchRingPts = searchRing.getCoordinates();

        if (innerRing == searchRing) continue;

        if (!innerRing
            .getEnvelopeInternal()
            .intersectsEnvelope(searchRing.getEnvelopeInternal())) continue;

        Coordinate? innerRingPt =
            IsValidOp.findPtNotNode(innerRingPts, searchRing, graph);

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
      index!.insert(env, ring);
    }
  }
}

/// Computes the boundary of a {@link Geometry}.
/// This operation will always return a {@link Geometry} of the appropriate
/// dimension for the boundary (even if the input geometry is empty).
/// The boundary of zero-dimensional geometries (Points) is
/// always the empty {@link GeometryCollection}.
///
/// @author Martin Davis
/// @version 1.7
class BoundaryOp {
  static Geometry? getBoundaryFromGeometry(Geometry g) {
    BoundaryOp bop = BoundaryOp(g);
    return bop.getBoundary();
  }

  Geometry geom;
  GeometryFactory geomFact = GeometryFactory.defaultPrecision();
  Mod2BoundaryNodeRule bnRule = Mod2BoundaryNodeRule();

  BoundaryOp(this.geom) {
    this.bnRule = BoundaryNodeRule.MOD2_BOUNDARY_RULE as Mod2BoundaryNodeRule;
  }

  BoundaryOp.withRule(this.geom, this.bnRule);

  Geometry getBoundary() {
    if (geom is LineString) return boundaryLineString(geom as LineString)!;
    if (geom is MultiLineString) {
      return boundaryMultiLineString(geom as MultiLineString);
    }
    return geom.getBoundary();
  }

  MultiPoint getEmptyMultiPoint() {
    return geomFact.createMultiPointEmpty();
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

  Map<Coordinate, Counter>? endpointMap;

  List<Coordinate> computeBoundaryCoordinates(MultiLineString mLine) {
    List<Coordinate> bdyPts = [];
    endpointMap = SplayTreeMap();
    for (int i = 0; i < mLine.getNumGeometries(); i++) {
      LineString line = mLine.getGeometryN(i) as LineString;
      if (line.isEmpty()) {
        continue;
      }
      addEndpoint(line.getCoordinateN(0));
      addEndpoint(line.getCoordinateN(line.getNumPoints() - 1));
    }

    endpointMap!.forEach((coord, counter) {
      int valence = counter.count;
      if (bnRule.isInBoundary(valence)) {
        bdyPts.add(coord);
      }
    });

    return bdyPts;
  }

  void addEndpoint(Coordinate pt) {
    Counter? counter = endpointMap![pt];
    if (counter == null) {
      counter = Counter();
      endpointMap![pt] = counter;
    }
    counter.count++;
  }

  Geometry? boundaryLineString(LineString line) {
    if (geom.isEmpty()) {
      return getEmptyMultiPoint();
    }

    if (line.isClosed()) {
      // check whether endpoints of valence 2 are on the boundary or not
      bool closedEndpointOnBoundary = bnRule.isInBoundary(2);
      if (closedEndpointOnBoundary) {
        return line.getStartPoint();
      } else {
        return geomFact.createMultiPointEmpty();
      }
    }
    return geomFact
        .createMultiPoint([line.getStartPoint()!, line.getEndPoint()!]);
  }
}

/// Stores an integer count, for use as a Map entry.
///
/// @author Martin Davis
/// @version 1.7
class Counter {
  /// The value of the count
  int count = 0;
}
