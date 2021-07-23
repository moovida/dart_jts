part of dart_jts;

/**
 * @version 1.7
 */
class NodeFactory {
  /**
   * The basic node constructor does not allow for incident edges
   */
  Node createNode(Coordinate coord) {
    return new Node(coord, null);
  }
}

/**
 * A map of nodes, indexed by the coordinate of the node
 * @version 1.7
 */
class NodeMap {
  //Map nodeMap = new HashMap();
  Map<Coordinate, Node?> nodeMap = new SplayTreeMap();
  NodeFactory nodeFact;

  NodeMap(this.nodeFact);

  /**
   * Factory function - subclasses can override to create their own types of nodes
   */
  /*
   Node createNode(Coordinate coord)
  {
    return new Node(coord);
  }
  */
  /**
   * This method expects that a node has a coordinate value.
   */
  Node addNodeFromCoordinate(Coordinate coord) {
    Node? node = nodeMap[coord];
    if (node == null) {
      node = nodeFact.createNode(coord);
      nodeMap[coord] = node;
    }
    return node;
  }

  Node addNode(Node n) {
    Node? node = nodeMap[n.getCoordinate()];
    if (node == null) {
      nodeMap[n.getCoordinate()] = n;
      return n;
    }
    node.mergeLabelFromNode(n);
    return node;
  }

  /**
   * Adds a node for the start point of this EdgeEnd
   * (if one does not already exist in this map).
   * Adds the EdgeEnd to the (possibly new) node.
   */
  void add(EdgeEnd e) {
    Coordinate p = e.getCoordinate()!;
    Node n = addNodeFromCoordinate(p);
    n.add(e);
  }

  /**
   * @return the node if found; null otherwise
   */
  Node? find(Coordinate coord) {
    return nodeMap[coord];
  }

  Iterator iterator() {
    return nodeMap.values.iterator;
  }

  List values() {
    return List.from(nodeMap.values);
  }

  List getBoundaryNodes(int geomIndex) {
    List bdyNodes = [];
    for (Iterator i = iterator(); i.moveNext();) {
      Node node = i.current;
      if (node.getLabel()!.getLocation(geomIndex) == Location.BOUNDARY)
        bdyNodes.add(node);
    }
    return bdyNodes;
  }

//   void print(PrintStream out)
//  {
//    for (Iterator it = iterator(); it.hasNext(); )
//    {
//      Node n = (Node) it.next();
//      n.print(out);
//    }
//  }
}

/**
 * The computation of the <code>IntersectionMatrix</code> relies on the use of a structure
 * called a "topology graph".  The topology graph contains nodes and edges
 * corresponding to the nodes and line segments of a <code>Geometry</code>. Each
 * node and edge in the graph is labeled with its topological location relative to
 * the source geometry.
 * <P>
 * Note that there is no requirement that points of self-intersection be a vertex.
 * Thus to obtain a correct topology graph, <code>Geometry</code>s must be
 * self-noded before constructing their graphs.
 * <P>
 * Two fundamental operations are supported by topology graphs:
 * <UL>
 *   <LI>Computing the intersections between all the edges and nodes of a single graph
 *   <LI>Computing the intersections between the edges and nodes of two different graphs
 * </UL>
 *
 * @version 1.7
 */
class PlanarGraph {
  /**
   * For nodes in the Collection, link the DirectedEdges at the node that are in the result.
   * This allows clients to link only a subset of nodes in the graph, for
   * efficiency (because they know that only a subset is of interest).
   */
  static void linkResultDirectedEdgesStatic(List nodes) {
    for (Iterator nodeit = nodes.iterator; nodeit.moveNext();) {
      Node node = nodeit.current;
      (node.getEdges() as DirectedEdgeStar).linkResultDirectedEdges();
    }
  }

  List edges = [];
  late NodeMap nodes;
  List edgeEndList = [];

  PlanarGraph.withFactory(NodeFactory nodeFact) {
    nodes = new NodeMap(nodeFact);
  }

  PlanarGraph() {
    nodes = new NodeMap(new NodeFactory());
  }

  Iterator getEdgeIterator() {
    return edges.iterator;
  }

  List getEdgeEnds() {
    return edgeEndList;
  }

  bool isBoundaryNode(int geomIndex, Coordinate coord) {
    Node? node = nodes.find(coord);
    if (node == null) return false;
    Label? label = node.getLabel();
    if (label != null && label.getLocation(geomIndex) == Location.BOUNDARY)
      return true;
    return false;
  }

  void insertEdge(Edge e) {
    edges.add(e);
  }

  void add(EdgeEnd e) {
    nodes.add(e);
    edgeEndList.add(e);
  }

  Iterator getNodeIterator() {
    return nodes.iterator();
  }

  List getNodes() {
    return nodes.values();
  }

  Node addNode(Node node) {
    return nodes.addNode(node);
  }

  Node addNodeFromCoordinate(Coordinate coord) {
    return nodes.addNodeFromCoordinate(coord);
  }

  /**
   * @return the node if found; null otherwise
   */
  Node? find(Coordinate coord) {
    return nodes.find(coord);
  }

  /**
   * Add a set of edges to the graph.  For each edge two DirectedEdges
   * will be created.  DirectedEdges are NOT linked by this method.
   */
  void addEdges(List edgesToAdd) {
    // create all the nodes for the edges
    for (Iterator it = edgesToAdd.iterator; it.moveNext();) {
      Edge e = it.current as Edge;
      edges.add(e);

      DirectedEdge de1 = new DirectedEdge(e, true);
      DirectedEdge de2 = new DirectedEdge(e, false);
      de1.setSym(de2);
      de2.setSym(de1);

      add(de1);
      add(de2);
    }
  }

  /**
   * Link the DirectedEdges at the nodes of the graph.
   * This allows clients to link only a subset of nodes in the graph, for
   * efficiency (because they know that only a subset is of interest).
   */
  void linkResultDirectedEdges() {
    for (Iterator nodeit = nodes.iterator(); nodeit.moveNext();) {
      Node node = nodeit.current as Node;
      (node.getEdges() as DirectedEdgeStar).linkResultDirectedEdges();
    }
  }

  /**
   * Link the DirectedEdges at the nodes of the graph.
   * This allows clients to link only a subset of nodes in the graph, for
   * efficiency (because they know that only a subset is of interest).
   */
  void linkAllDirectedEdges() {
    for (Iterator nodeit = nodes.iterator(); nodeit.moveNext();) {
      Node node = nodeit.current as Node;
      (node.getEdges() as DirectedEdgeStar).linkAllDirectedEdges();
    }
  }

  /**
   * Returns the EdgeEnd which has edge e as its base edge
   * (MD 18 Feb 2002 - this should return a pair of edges)
   *
   * @return the edge, if found
   *    <code>null</code> if the edge was not found
   */
  EdgeEnd? findEdgeEnd(Edge e) {
    for (Iterator i = getEdgeEnds().iterator; i.moveNext();) {
      EdgeEnd ee = i.current as EdgeEnd;
      if (ee.getEdge() == e) return ee;
    }
    return null;
  }

  /**
   * Returns the edge whose first two coordinates are p0 and p1
   *
   * @return the edge, if found
   *    <code>null</code> if the edge was not found
   */
  Edge? findEdge(Coordinate p0, Coordinate p1) {
    for (int i = 0; i < edges.length; i++) {
      Edge e = edges[i] as Edge;
      List<Coordinate> eCoord = e.getCoordinates();
      if (p0.equals(eCoord[0]) && p1.equals(eCoord[1])) return e;
    }
    return null;
  }

  /**
   * Returns the edge which starts at p0 and whose first segment is
   * parallel to p1
   *
   * @return the edge, if found
   *    <code>null</code> if the edge was not found
   */
  Edge? findEdgeInSameDirection(Coordinate p0, Coordinate p1) {
    for (int i = 0; i < edges.length; i++) {
      Edge e = edges[i];

      List<Coordinate> eCoord = e.getCoordinates();
      if (matchInSameDirection(p0, p1, eCoord[0], eCoord[1])) return e;

      if (matchInSameDirection(
          p0, p1, eCoord[eCoord.length - 1], eCoord[eCoord.length - 2]))
        return e;
    }
    return null;
  }

  /**
   * The coordinate pairs match if they define line segments lying in the same direction.
   * E.g. the segments are parallel and in the same quadrant
   * (as opposed to parallel and opposite!).
   */
  bool matchInSameDirection(
      Coordinate p0, Coordinate p1, Coordinate ep0, Coordinate ep1) {
    if (!p0.equals(ep0)) return false;

    if (Orientation.index(p0, p1, ep1) == Orientation.COLLINEAR &&
        Quadrant.quadrantFromCoords(p0, p1) ==
            Quadrant.quadrantFromCoords(ep0, ep1)) return true;
    return false;
  }

//   void printEdges(PrintStream out)
//  {
//    out.println("Edges:");
//    for (int i = 0; i < edges.size(); i++) {
//      out.println("edge " + i + ":");
//      Edge e = (Edge) edges.get(i);
//      e.print(out);
//      e.eiList.print(out);
//    }
//  }
  void debugPrint(Object o) {
    print(o);
  }

  void debugPrintln(Object o) {
    print(o);
    print("\n");
  }
}

/**
 * A GeometryGraph is a graph that models a given Geometry
 * @version 1.7
 */
class GeometryGraph extends PlanarGraph {
  /**
   * This method implements the Boundary Determination Rule
   * for determining whether
   * a component (node or edge) that appears multiple times in elements
   * of a MultiGeometry is in the boundary or the interior of the Geometry
   * <br>
   * The SFS uses the "Mod-2 Rule", which this function implements
   * <br>
   * An alternative (and possibly more intuitive) rule would be
   * the "At Most One Rule":
   *    isInBoundary = (componentCount == 1)
   */
/*
   static bool isInBoundary(int boundaryCount)
  {
    // the "Mod-2 Rule"
    return boundaryCount % 2 == 1;
  }
   static int determineBoundary(int boundaryCount)
  {
    return isInBoundary(boundaryCount) ? Location.BOUNDARY : Location.INTERIOR;
  }
*/

  static int determineBoundary(
      BoundaryNodeRule boundaryNodeRule, int boundaryCount) {
    return boundaryNodeRule.isInBoundary(boundaryCount)
        ? Location.BOUNDARY
        : Location.INTERIOR;
  }

  Geometry? parentGeom;

  /**
   * The lineEdgeMap is a map of the linestring components of the
   * parentGeometry to the edges which are derived from them.
   * This is used to efficiently perform findEdge queries
   */
  Map<LineString, Edge> lineEdgeMap = {};

  BoundaryNodeRule boundaryNodeRule;

  /**
   * If this flag is true, the Boundary Determination Rule will used when deciding
   * whether nodes are in the boundary or not
   */
  bool useBoundaryDeterminationRule = true;
  int argIndex =
      0; // the index of this geometry as an argument to a spatial function (used for labelling)
  List? boundaryNodes;
  bool _hasTooFewPoints = false;
  Coordinate? invalidPoint;

  PointOnGeometryLocator? areaPtLocator;

  // for use if geometry is not Polygonal
  final PointLocator ptLocator = new PointLocator();

  EdgeSetIntersector createEdgeSetIntersector() {
    // various options for computing intersections, from slowest to fastest

    // EdgeSetIntersector esi = new SimpleEdgeSetIntersector();
    // EdgeSetIntersector esi = new MonotoneChainIntersector();
    // EdgeSetIntersector esi = new NonReversingChainIntersector();
    // EdgeSetIntersector esi = new SimpleSweepLineIntersector();
    // EdgeSetIntersector esi = new MCSweepLineIntersector();

    //return new SimpleEdgeSetIntersector();
    return new SimpleMCSweepLineIntersector();
  }

  GeometryGraph(int argIndex, Geometry parentGeom)
      : this.args3(
            argIndex, parentGeom, BoundaryNodeRule.OGC_SFS_BOUNDARY_RULE);

  GeometryGraph.args3(int argIndex, this.parentGeom, this.boundaryNodeRule) {
    this.argIndex = argIndex;
    if (parentGeom != null) {
//      precisionModel = parentGeom.getPrecisionModel();
//      SRID = parentGeom.getSRID();
      addGeometry(parentGeom!);
    }
  }

  /**
   * This constructor is used by clients that wish to add Edges explicitly,
   * rather than adding a Geometry.  (An example is BufferOp).
   */
  // no longer used
//   GeometryGraph(int argIndex, PrecisionModel precisionModel, int SRID) {
//    this(argIndex, null);
//    this.precisionModel = precisionModel;
//    this.SRID = SRID;
//  }
//   PrecisionModel getPrecisionModel()
//  {
//    return precisionModel;
//  }
//   int getSRID() { return SRID; }

  bool hasTooFewPoints() {
    return _hasTooFewPoints;
  }

  Coordinate? getInvalidPoint() {
    return invalidPoint;
  }

  Geometry? getGeometry() {
    return parentGeom;
  }

  BoundaryNodeRule getBoundaryNodeRule() {
    return boundaryNodeRule;
  }

  List getBoundaryNodes() {
    if (boundaryNodes == null) boundaryNodes = nodes.getBoundaryNodes(argIndex);
    return boundaryNodes!;
  }

  List<Coordinate> getBoundaryPoints() {
    List coll = getBoundaryNodes();
    List<Coordinate> pts = []; //..length = (coll.length);
    int i = 0;
    for (Iterator it = coll.iterator; it.moveNext();) {
      Node node = it.current as Node;
      pts.add(node.getCoordinate().copy());
      // pts[i++] = node.getCoordinate().copy();
    }
    return pts;
  }

  Edge? findEdgeFromLine(LineString line) {
    return lineEdgeMap[line];
  }

  void computeSplitEdges(List edgelist) {
    for (Iterator i = edges.iterator; i.moveNext();) {
      Edge e = i.current as Edge;
      e.eiList.addSplitEdges(edgelist);
    }
  }

  void addGeometry(Geometry g) {
    if (g.isEmpty()) return;

    // check if this Geometry should obey the Boundary Determination Rule
    // all collections except MultiPolygons obey the rule
    if (g is MultiPolygon) useBoundaryDeterminationRule = false;

    if (g is Polygon)
      addPolygon(g);
    // LineString also handles LinearRings
    else if (g is LineString)
      addLineString(g);
    else if (g is Point)
      addPoint(g);
    else if (g is MultiPoint)
      addCollection(g);
    else if (g is MultiLineString)
      addCollection(g);
    else if (g is MultiPolygon)
      addCollection(g);
    else if (g is GeometryCollection)
      addCollection(g);
    else
      throw new UnsupportedError(g.runtimeType.toString());
  }

  void addCollection(GeometryCollection gc) {
    for (int i = 0; i < gc.getNumGeometries(); i++) {
      Geometry g = gc.getGeometryN(i);
      addGeometry(g);
    }
  }

  /**
   * Add a Point to the graph.
   */
  void addPoint(Point p) {
    Coordinate coord = p.getCoordinate()!;
    insertPoint(argIndex, coord, Location.INTERIOR);
  }

  /**
   * Adds a polygon ring to the graph.
   * Empty rings are ignored.
   *
   * The left and right topological location arguments assume that the ring is oriented CW.
   * If the ring is in the opposite orientation,
   * the left and right locations must be interchanged.
   */
  void addPolygonRing(LinearRing lr, int cwLeft, int cwRight) {
    // don't bother adding empty holes
    if (lr.isEmpty()) return;

    List<Coordinate> coord =
        CoordinateArrays.removeRepeatedPoints(lr.getCoordinates());

    if (coord.length < 4) {
      _hasTooFewPoints = true;
      invalidPoint = coord[0];
      return;
    }

    int left = cwLeft;
    int right = cwRight;
    if (Orientation.isCCW(coord)) {
      left = cwRight;
      right = cwLeft;
    }
    Edge e = Edge(coord, Label.args4(argIndex, Location.BOUNDARY, left, right));
    lineEdgeMap[lr] = e;

    insertEdge(e);
    // insert the endpoint as a node, to mark that it is on the boundary
    insertPoint(argIndex, coord[0], Location.BOUNDARY);
  }

  void addPolygon(Polygon p) {
    addPolygonRing(p.getExteriorRing(), Location.EXTERIOR, Location.INTERIOR);

    for (int i = 0; i < p.getNumInteriorRing(); i++) {
      LinearRing hole = p.getInteriorRingN(i);

      // Holes are topologically labelled opposite to the shell, since
      // the interior of the polygon lies on their opposite side
      // (on the left, if the hole is oriented CW)
      addPolygonRing(hole, Location.INTERIOR, Location.EXTERIOR);
    }
  }

  void addLineString(LineString line) {
    List<Coordinate> coord =
        CoordinateArrays.removeRepeatedPoints(line.getCoordinates());

    if (coord.length < 2) {
      _hasTooFewPoints = true;
      invalidPoint = coord[0];
      return;
    }

    // add the edge for the LineString
    // line edges do not have locations for their left and right sides
    Edge e = new Edge(coord, new Label.args2(argIndex, Location.INTERIOR));
    lineEdgeMap[line] = e;
    insertEdge(e);
    /**
     * Add the boundary points of the LineString, if any.
     * Even if the LineString is closed, add both points as if they were endpoints.
     * This allows for the case that the node already exists and is a boundary point.
     */
    assert(coord.length >= 2, "found LineString with single point");
    insertBoundaryPoint(argIndex, coord[0]);
    insertBoundaryPoint(argIndex, coord[coord.length - 1]);
  }

  /**
   * Add an Edge computed externally.  The label on the Edge is assumed
   * to be correct.
   */
  void addEdge(Edge e) {
    insertEdge(e);
    List<Coordinate> coord = e.getCoordinates();
    // insert the endpoint as a node, to mark that it is on the boundary
    insertPoint(argIndex, coord[0], Location.BOUNDARY);
    insertPoint(argIndex, coord[coord.length - 1], Location.BOUNDARY);
  }

  /**
   * Add a point computed externally.  The point is assumed to be a
   * Point Geometry part, which has a location of INTERIOR.
   */
  void addPointFromCoordinate(Coordinate pt) {
    insertPoint(argIndex, pt, Location.INTERIOR);
  }

  /**
   * Compute self-nodes, taking advantage of the Geometry type to
   * minimize the number of intersection tests.  (E.g. rings are
   * not tested for self-intersection, since they are assumed to be valid).
   *
   * @param li the LineIntersector to use
   * @param computeRingSelfNodes if <code>false</code>, intersection checks are optimized to not test rings for self-intersection
   * @return the computed SegmentIntersector containing information about the intersections found
   */
  SegmentIntersector computeSelfNodes(
      LineIntersector li, bool computeRingSelfNodes) {
    return computeSelfNodes3(li, computeRingSelfNodes, false);
  }

  /**
   * Compute self-nodes, taking advantage of the Geometry type to
   * minimize the number of intersection tests.  (E.g. rings are
   * not tested for self-intersection, since they are assumed to be valid).
   *
   * @param li the LineIntersector to use
   * @param computeRingSelfNodes if <code>false</code>, intersection checks are optimized to not test rings for self-intersection
   * @param isDoneIfProperInt short-circuit the intersection computation if a proper intersection is found
   * @return the computed SegmentIntersector containing information about the intersections found
   */
  SegmentIntersector computeSelfNodes3(
      LineIntersector li, bool computeRingSelfNodes, bool isDoneIfProperInt) {
    SegmentIntersector si = new SegmentIntersector(li, true, false);
    si.setIsDoneIfProperInt(isDoneIfProperInt);
    EdgeSetIntersector esi = createEdgeSetIntersector();
    // optimize intersection search for valid Polygons and LinearRings
    bool isRings = parentGeom is LinearRing ||
        parentGeom is Polygon ||
        parentGeom is MultiPolygon;
    bool computeAllSegments = computeRingSelfNodes || !isRings;
    esi.computeIntersections(edges, si, computeAllSegments);

    //System.out.println("SegmentIntersector # tests = " + si.numTests);
    addSelfIntersectionNodes(argIndex);
    return si;
  }

  SegmentIntersector computeEdgeIntersections(
      GeometryGraph g, LineIntersector li, bool includeProper) {
    SegmentIntersector si = new SegmentIntersector(li, includeProper, true);
    si.setBoundaryNodes(this.getBoundaryNodes(), g.getBoundaryNodes());

    EdgeSetIntersector esi = createEdgeSetIntersector();
    esi.computeIntersections3(edges, g.edges, si);
/*
for (Iterator i = g.edges.iterator; i.moveNext();) {
Edge e = (Edge) i.current
Debug.print(e.getEdgeIntersectionList());
}
*/
    return si;
  }

  void insertPoint(int argIndex, Coordinate coord, int onLocation) {
    Node n = nodes.addNodeFromCoordinate(coord);
    Label? lbl = n.getLabel();
    if (lbl == null) {
      n.label = new Label.args2(argIndex, onLocation);
    } else
      lbl.setLocationWithIndex(argIndex, onLocation);
  }

  /**
   * Adds candidate boundary points using the current {@link BoundaryNodeRule}.
   * This is used to add the boundary
   * points of dim-1 geometries (Curves/MultiCurves).
   */
  void insertBoundaryPoint(int argIndex, Coordinate coord) {
    Node n = nodes.addNodeFromCoordinate(coord);
    // nodes always have labels
    Label? lbl = n.getLabel();
    // the new point to insert is on a boundary
    int boundaryCount = 1;
    // determine the current location for the point (if any)
    int loc = Location.NONE;
    loc = lbl!.getLocationWithPosIndex(argIndex, Position.ON);
    if (loc == Location.BOUNDARY) boundaryCount++;

    // determine the boundary status of the point according to the Boundary Determination Rule
    int newLoc = determineBoundary(boundaryNodeRule, boundaryCount);
    lbl.setLocationWithIndex(argIndex, newLoc);
  }

  void addSelfIntersectionNodes(int argIndex) {
    for (Iterator i = edges.iterator; i.moveNext();) {
      Edge e = i.current;
      int eLoc = e.getLabel()!.getLocation(argIndex);
      for (Iterator eiIt = e.eiList.iterator(); eiIt.moveNext();) {
        EdgeIntersection ei = eiIt.current as EdgeIntersection;
        addSelfIntersectionNode(argIndex, ei.coord, eLoc);
      }
    }
  }

  /**
   * Add a node for a self-intersection.
   * If the node is a potential boundary node (e.g. came from an edge which
   * is a boundary) then insert it as a potential boundary node.
   * Otherwise, just add it as a regular node.
   */
  void addSelfIntersectionNode(int argIndex, Coordinate coord, int loc) {
    // if this node is already a boundary node, don't change it
    if (isBoundaryNode(argIndex, coord)) return;
    if (loc == Location.BOUNDARY && useBoundaryDeterminationRule)
      insertBoundaryPoint(argIndex, coord);
    else
      insertPoint(argIndex, coord, loc);
  }

  // MD - experimental for now
  /**
   * Determines the {@link Location} of the given {@link Coordinate}
   * in this geometry.
   *
   * @param pt the point to test
   * @return the location of the point in the geometry
   */
  int locate(Coordinate pt) {
    if (parentGeom is Polygonal && parentGeom!.getNumGeometries() > 50) {
      // lazily init point locator
      if (areaPtLocator == null) {
        areaPtLocator = new IndexedPointInAreaLocator(parentGeom!);
      }
      return areaPtLocator!.locate(pt);
    }
    return ptLocator.locate(pt, parentGeom!);
  }
}

/**
 * A DirectedEdgeStar is an ordered list of <b>outgoing</b> DirectedEdges around a node.
 * It supports labelling the edges as well as linking the edges to form both
 * MaximalEdgeRings and MinimalEdgeRings.
 *
 * @version 1.7
 */
class DirectedEdgeStar extends EdgeEndStar {
  /**
   * A list of all outgoing edges in the result, in CCW order
   */
  List? resultAreaEdgeList;
  Label? label;

  /**
   * Insert a directed edge in the list
   */
  void insert(EdgeEnd ee) {
    DirectedEdge de = ee as DirectedEdge;
    insertEdgeEnd(de, de);
  }

  Label? getLabel() {
    return label;
  }

  int getOutgoingDegree() {
    int degree = 0;
    for (Iterator it = iterator(); it.moveNext();) {
      DirectedEdge de = it.current as DirectedEdge;
      if (de.isInResult()) degree++;
    }
    return degree;
  }

  int getOutgoingDegreeWithRing(EdgeRing er) {
    int degree = 0;
    for (Iterator it = iterator(); it.moveNext();) {
      DirectedEdge de = it.current as DirectedEdge;
      if (de.getEdgeRing() == er) degree++;
    }
    return degree;
  }

  DirectedEdge? getRightmostEdge() {
    List edges = getEdges();
    int size = edges.length;
    if (size < 1) return null;
    DirectedEdge de0 = edges[0] as DirectedEdge;
    if (size == 1) return de0;
    DirectedEdge deLast = edges[size - 1] as DirectedEdge;

    int quad0 = de0.getQuadrant();
    int quad1 = deLast.getQuadrant();
    if (Quadrant.isNorthern(quad0) && Quadrant.isNorthern(quad1))
      return de0;
    else if (!Quadrant.isNorthern(quad0) && !Quadrant.isNorthern(quad1))
      return deLast;
    else {
      // edges are in different hemispheres - make sure we return one that is non-horizontal
      //Assert.isTrue(de0.getDy() != 0, "should never return horizontal edge!");
      DirectedEdge? nonHorizontalEdge = null;
      if (de0.getDy() != 0)
        return de0;
      else if (deLast.getDy() != 0) return deLast;
    }
// TODO   assert("found two horizontal edges incident on node");
    return null;
  }

  /**
   * Compute the labelling for all dirEdges in this star, as well
   * as the overall labelling
   */
  void computeLabelling(List<GeometryGraph> geom) {
//Debug.print(this);
    super.computeLabelling(geom);

    // determine the overall labelling for this DirectedEdgeStar
    // (i.e. for the node it is based at)
    label = new Label(Location.NONE);
    for (Iterator it = iterator(); it.moveNext();) {
      EdgeEnd ee = it.current as EdgeEnd;
      Edge e = ee.getEdge();
      Label eLabel = e.getLabel()!;
      for (int i = 0; i < 2; i++) {
        int eLoc = eLabel.getLocation(i);
        if (eLoc == Location.INTERIOR || eLoc == Location.BOUNDARY)
          label!.setLocationWithIndex(i, Location.INTERIOR);
      }
    }
//Debug.print(this);
  }

  /**
   * For each dirEdge in the star,
   * merge the label from the sym dirEdge into the label
   */
  void mergeSymLabels() {
    for (Iterator it = iterator(); it.moveNext();) {
      DirectedEdge de = it.current as DirectedEdge;
      Label label = de.getLabel()!;
      label.merge(de.getSym().getLabel()!);
    }
  }

  /**
   * Update incomplete dirEdge labels from the labelling for the node
   */
  void updateLabelling(Label nodeLabel) {
    for (Iterator it = iterator(); it.moveNext();) {
      DirectedEdge de = it.current as DirectedEdge;
      Label label = de.getLabel()!;
      label.setAllLocationsIfNullWithIndex(0, nodeLabel.getLocation(0));
      label.setAllLocationsIfNullWithIndex(1, nodeLabel.getLocation(1));
    }
  }

  List getResultAreaEdges() {
//print(System.out);
    if (resultAreaEdgeList != null) return resultAreaEdgeList!;
    resultAreaEdgeList = [];
    for (Iterator it = iterator(); it.moveNext();) {
      DirectedEdge de = it.current as DirectedEdge;
      if (de.isInResult() || de.getSym().isInResult())
        resultAreaEdgeList!.add(de);
    }
    return resultAreaEdgeList!;
  }

  static const SCANNING_FOR_INCOMING = 1;
  static const LINKING_TO_OUTGOING = 2;

  /**
   * Traverse the star of DirectedEdges, linking the included edges together.
   * To link two dirEdges, the <code>next</code> pointer for an incoming dirEdge
   * is set to the next outgoing edge.
   * <p>
   * DirEdges are only linked if:
   * <ul>
   * <li>they belong to an area (i.e. they have sides)
   * <li>they are marked as being in the result
   * </ul>
   * <p>
   * Edges are linked in CCW order (the order they are stored).
   * This means that rings have their face on the Right
   * (in other words,
   * the topological location of the face is given by the RHS label of the DirectedEdge)
   * <p>
   * PRECONDITION: No pair of dirEdges are both marked as being in the result
   */
  void linkResultDirectedEdges() {
    // make sure edges are copied to resultAreaEdges list
    getResultAreaEdges();
    // find first area edge (if any) to start linking at
    DirectedEdge? firstOut = null;
    DirectedEdge? incoming = null;
    int state = SCANNING_FOR_INCOMING;
    // link edges in CCW order
    for (int i = 0; i < resultAreaEdgeList!.length; i++) {
      DirectedEdge nextOut = resultAreaEdgeList![i] as DirectedEdge;
      DirectedEdge nextIn = nextOut.getSym();

      // skip de's that we're not interested in
      if (!nextOut.getLabel()!.isArea()) continue;

      // record first outgoing edge, in order to link the last incoming edge
      if (firstOut == null && nextOut.isInResult()) firstOut = nextOut;
      // assert: sym.isInResult() == false, since pairs of dirEdges should have been removed already

      switch (state) {
        case SCANNING_FOR_INCOMING:
          if (!nextIn.isInResult()) continue;
          incoming = nextIn;
          state = LINKING_TO_OUTGOING;
          break;
        case LINKING_TO_OUTGOING:
          if (!nextOut.isInResult()) continue;
          incoming!.setNext(nextOut);
          state = SCANNING_FOR_INCOMING;
          break;
      }
    }
//Debug.print(this);
    if (state == LINKING_TO_OUTGOING) {
//Debug.print(firstOut == null, this);
      if (firstOut == null)
        throw new TopologyException(
            "no outgoing dirEdge found ${getCoordinate()}");
      //Assert.isTrue(firstOut != null, "no outgoing dirEdge found (at " + getCoordinate() );
      assert(firstOut.isInResult(), "unable to link last incoming dirEdge");
      incoming!.setNext(firstOut);
    }
  }

  void linkMinimalDirectedEdges(EdgeRing er) {
    // find first area edge (if any) to start linking at
    DirectedEdge? firstOut = null;
    DirectedEdge? incoming = null;
    int state = SCANNING_FOR_INCOMING;
    // link edges in CW order
    for (int i = resultAreaEdgeList!.length - 1; i >= 0; i--) {
      DirectedEdge nextOut = resultAreaEdgeList![i] as DirectedEdge;
      DirectedEdge nextIn = nextOut.getSym();

      // record first outgoing edge, in order to link the last incoming edge
      if (firstOut == null && nextOut.getEdgeRing() == er) firstOut = nextOut;

      switch (state) {
        case SCANNING_FOR_INCOMING:
          if (nextIn.getEdgeRing() != er) continue;
          incoming = nextIn;
          state = LINKING_TO_OUTGOING;
          break;
        case LINKING_TO_OUTGOING:
          if (nextOut.getEdgeRing() != er) continue;
          incoming!.setNextMin(nextOut);
          state = SCANNING_FOR_INCOMING;
          break;
      }
    }
//print(System.out);
    if (state == LINKING_TO_OUTGOING) {
      assert(firstOut != null, "found null for first outgoing dirEdge");
      assert(firstOut!.getEdgeRing() == er,
          "unable to link last incoming dirEdge");
      incoming!.setNextMin(firstOut!);
    }
  }

  void linkAllDirectedEdges() {
    getEdges();
    // find first area edge (if any) to start linking at
    DirectedEdge? prevOut = null;
    DirectedEdge? firstIn = null;
    // link edges in CW order
    for (int i = edgeList!.length - 1; i >= 0; i--) {
      DirectedEdge nextOut = edgeList![i] as DirectedEdge;
      DirectedEdge nextIn = nextOut.getSym();
      if (firstIn == null) firstIn = nextIn;
      if (prevOut != null) nextIn.setNext(prevOut);
      // record outgoing edge, in order to link the last incoming edge
      prevOut = nextOut;
    }
    firstIn!.setNext(prevOut!);
//Debug.print(this);
  }

  /**
   * Traverse the star of edges, maintaining the current location in the result
   * area at this node (if any).
   * If any L edges are found in the interior of the result, mark them as covered.
   */
  void findCoveredLineEdges() {
//Debug.print("findCoveredLineEdges");
//Debug.print(this);
    // Since edges are stored in CCW order around the node,
    // as we move around the ring we move from the right to the left side of the edge

    /**
     * Find first DirectedEdge of result area (if any).
     * The interior of the result is on the RHS of the edge,
     * so the start location will be:
     * - INTERIOR if the edge is outgoing
     * - EXTERIOR if the edge is incoming
     */
    int startLoc = Location.NONE;
    for (Iterator it = iterator(); it.moveNext();) {
      DirectedEdge nextOut = it.current as DirectedEdge;
      DirectedEdge nextIn = nextOut.getSym();
      if (!nextOut.isLineEdge()) {
        if (nextOut.isInResult()) {
          startLoc = Location.INTERIOR;
          break;
        }
        if (nextIn.isInResult()) {
          startLoc = Location.EXTERIOR;
          break;
        }
      }
    }
    // no A edges found, so can't determine if L edges are covered or not
    if (startLoc == Location.NONE) return;

    /**
     * move around ring, keeping track of the current location
     * (Interior or Exterior) for the result area.
     * If L edges are found, mark them as covered if they are in the interior
     */
    int currLoc = startLoc;
    for (Iterator it = iterator(); it.moveNext();) {
      DirectedEdge nextOut = it.current as DirectedEdge;
      DirectedEdge nextIn = nextOut.getSym();
      if (nextOut.isLineEdge()) {
        nextOut.getEdge().setCovered(currLoc == Location.INTERIOR);
//Debug.println(nextOut);
      } else {
        // edge is an Area edge
        if (nextOut.isInResult()) currLoc = Location.EXTERIOR;
        if (nextIn.isInResult()) currLoc = Location.INTERIOR;
      }
    }
  }

  void computeDepths(DirectedEdge de) {
    int edgeIndex = findIndex(de);
    int startDepth = de.getDepth(Position.LEFT);
    int targetLastDepth = de.getDepth(Position.RIGHT);
    // compute the depths from this edge up to the end of the edge array
    int nextDepth = computeDepths3(edgeIndex + 1, edgeList!.length, startDepth);
    // compute the depths for the initial part of the array
    int lastDepth = computeDepths3(0, edgeIndex, nextDepth);
//Debug.print(lastDepth != targetLastDepth, this);
//Debug.print(lastDepth != targetLastDepth, "mismatch: " + lastDepth + " / " + targetLastDepth);
    if (lastDepth != targetLastDepth)
      throw new TopologyException("depth mismatch at ${de.getCoordinate()}");
    //Assert.isTrue(lastDepth == targetLastDepth, "depth mismatch at " + de.getCoordinate());
  }

  /**
   * Compute the DirectedEdge depths for a subsequence of the edge array.
   *
   * @return the last depth assigned (from the R side of the last edge visited)
   */
  int computeDepths3(int startIndex, int endIndex, int startDepth) {
    int currDepth = startDepth;
    for (int i = startIndex; i < endIndex; i++) {
      DirectedEdge nextDe = edgeList![i] as DirectedEdge;
      nextDe.setEdgeDepths(Position.RIGHT, currDepth);
      currDepth = nextDe.getDepth(Position.LEFT);
    }
    return currDepth;
  }

//   void print(PrintStream out)
//  {
//    System.out.println("DirectedEdgeStar: " + getCoordinate());
//    for (Iterator it = iterator; it.moveNext(); ) {
//      DirectedEdge de = (DirectedEdge) it.current;
//      out.print("out ");
//      de.print(out);
//      out.println();
//      out.print("in ");
//      de.getSym().print(out);
//      out.println();
//    }
//  }
}

/**
 * Represents a point on an
 * edge which intersects with another edge.
 * <p>
 * The intersection may either be a single point, or a line segment
 * (in which case this point is the start of the line segment)
 * The intersection point must be precise.
 *
 * @version 1.7
 */
class EdgeIntersection implements Comparable {
  late Coordinate coord; // the point of intersection
  int segmentIndex; // the index of the containing line segment in the parent edge
  double
      dist; // the edge distance of this point along the containing line segment

  EdgeIntersection(Coordinate coord, this.segmentIndex, this.dist) {
    this.coord = new Coordinate.fromCoordinate(coord);
  }

  Coordinate getCoordinate() {
    return coord;
  }

  int getSegmentIndex() {
    return segmentIndex;
  }

  double getDistance() {
    return dist;
  }

  int compareTo(dynamic obj) {
    EdgeIntersection other = obj as EdgeIntersection;
    return compare(other.segmentIndex, other.dist);
  }

  /**
   * @return -1 this EdgeIntersection is located before the argument location
   * @return 0 this EdgeIntersection is at the argument location
   * @return 1 this EdgeIntersection is located after the argument location
   */
  int compare(int segmentIndex, double dist) {
    if (this.segmentIndex < segmentIndex) return -1;
    if (this.segmentIndex > segmentIndex) return 1;
    if (this.dist < dist) return -1;
    if (this.dist > dist) return 1;
    return 0;
  }

  bool isEndPoint(int maxSegmentIndex) {
    if (segmentIndex == 0 && dist == 0.0) return true;
    if (segmentIndex == maxSegmentIndex) return true;
    return false;
  }

//   void print(PrintStream out)
//  {
//    out.print(coord);
//    out.print(" seg # = " + segmentIndex);
//    out.println(" dist = " + dist);
//  }
  String toString() {
    return "$coord seg # = $segmentIndex dist = $dist";
  }
}

/**
 * @version 1.7
 */
abstract class EdgeRing {
  DirectedEdge?
      startDe; // the directed edge which starts the list of edges for this EdgeRing
  int maxNodeDegree = -1;
  List edges = []; // the DirectedEdges making up this EdgeRing
  List pts = [];
  Label label = new Label(Location
      .NONE); // label stores the locations of each geometry on the face surrounded by this ring
  LinearRing? ring; // the ring created for this EdgeRing
  bool _isHole = false;
  EdgeRing?
      shell; // if non-null, the ring is a hole and this EdgeRing is its containing shell
  List holes = []; // a list of EdgeRings which are holes in this EdgeRing

  GeometryFactory geometryFactory;

  EdgeRing(DirectedEdge start, this.geometryFactory) {
    computePoints(start);
    computeRing();
  }

  bool isIsolated() {
    return (label.getGeometryCount() == 1);
  }

  bool isHole() {
    //computePoints();
    return _isHole;
  }

  Coordinate getCoordinate(int i) {
    return pts[i];
  }

  LinearRing? getLinearRing() {
    return ring;
  }

  Label getLabel() {
    return label;
  }

  bool isShell() {
    return shell == null;
  }

  EdgeRing? getShell() {
    return shell;
  }

  void setShell(EdgeRing? shell) {
    this.shell = shell;
    if (shell != null) shell.addHole(this);
  }

  void addHole(EdgeRing ring) {
    holes.add(ring);
  }

  Polygon toPolygon(GeometryFactory geometryFactory) {
    List<LinearRing> holeLR = []; //..length = (holes.length);
    for (int i = 0; i < holes.length; i++) {
      holeLR.add((holes[i] as EdgeRing).getLinearRing()!);
      // holeLR[i] = (holes[i] as EdgeRing).getLinearRing()!;
    }
    Polygon poly = geometryFactory.createPolygon(getLinearRing(), holeLR);
    return poly;
  }

  /**
   * Compute a LinearRing from the point list previously collected.
   * Test if the ring is a hole (i.e. if it is CCW) and set the hole flag
   * accordingly.
   */
  void computeRing() {
    if (ring != null) return; // don't compute more than once
    List<Coordinate> coord = []; //..length = (pts.length);
    for (int i = 0; i < pts.length; i++) {
      coord.add(pts[i]);
      // coord[i] = pts[i];
    }
    ring = geometryFactory.createLinearRing(coord);
    _isHole = Orientation.isCCW(ring!.getCoordinates());
//Debug.println( (isHole ? "hole - " : "shell - ") + WKTWriter.toLineString(new CoordinateArraySequence(ring.getCoordinates())));
  }

  DirectedEdge getNext(DirectedEdge de);

  void setEdgeRing(DirectedEdge de, EdgeRing er);

  /**
   * Returns the list of DirectedEdges that make up this EdgeRing
   */
  List getEdges() {
    return edges;
  }

  /**
   * Collect all the points from the DirectedEdges of this ring into a contiguous list
   */
  void computePoints(DirectedEdge start) {
//System.out.println("buildRing");
    startDe = start;
    DirectedEdge de = start;
    bool isFirstEdge = true;
    do {
//      Assert.isTrue(de != null, "found null Directed Edge");
      if (de == null) throw new TopologyException("Found null DirectedEdge");
      if (de.getEdgeRing() == this)
        throw new TopologyException(
            "Directed Edge visited twice during ring-building at ${de.getCoordinate()}");

      edges.add(de);
//Debug.println(de);
//Debug.println(de.getEdge());
      Label label = de.getLabel()!;
      assert(label.isArea());
      mergeLabel(label);
      addPoints(de.getEdge(), de.isForward(), isFirstEdge);
      isFirstEdge = false;
      setEdgeRing(de, this);
      de = getNext(de);
    } while (de != startDe);
  }

  int getMaxNodeDegree() {
    if (maxNodeDegree < 0) computeMaxNodeDegree();
    return maxNodeDegree;
  }

  void computeMaxNodeDegree() {
    maxNodeDegree = 0;
    DirectedEdge de = startDe!;
    do {
      Node node = de.getNode()!;
      int degree =
          (node.getEdges() as DirectedEdgeStar).getOutgoingDegreeWithRing(this);
      if (degree > maxNodeDegree) maxNodeDegree = degree;
      de = getNext(de);
    } while (de != startDe);
    maxNodeDegree *= 2;
  }

  void setInResult() {
    DirectedEdge de = startDe!;
    do {
      de.getEdge().setInResult(true);
      de = de.getNext();
    } while (de != startDe);
  }

  void mergeLabel(Label deLabel) {
    mergeLabelWithIndex(deLabel, 0);
    mergeLabelWithIndex(deLabel, 1);
  }

  /**
   * Merge the RHS label from a DirectedEdge into the label for this EdgeRing.
   * The DirectedEdge label may be null.  This is acceptable - it results
   * from a node which is NOT an intersection node between the Geometries
   * (e.g. the end node of a LinearRing).  In this case the DirectedEdge label
   * does not contribute any information to the overall labelling, and is simply skipped.
   */
  void mergeLabelWithIndex(Label deLabel, int geomIndex) {
    int loc = deLabel.getLocationWithPosIndex(geomIndex, Position.RIGHT);
    // no information to be had from this label
    if (loc == Location.NONE) return;
    // if there is no current RHS value, set it
    if (label.getLocation(geomIndex) == Location.NONE) {
      label.setLocationWithIndex(geomIndex, loc);
      return;
    }
  }

  void addPoints(Edge edge, bool isForward, bool isFirstEdge) {
    List<Coordinate> edgePts = edge.getCoordinates();
    if (isForward) {
      int startIndex = 1;
      if (isFirstEdge) startIndex = 0;
      for (int i = startIndex; i < edgePts.length; i++) {
        pts.add(edgePts[i]);
      }
    } else {
      // is backward
      int startIndex = edgePts.length - 2;
      if (isFirstEdge) startIndex = edgePts.length - 1;
      for (int i = startIndex; i >= 0; i--) {
        pts.add(edgePts[i]);
      }
    }
  }

  /**
   * This method will cause the ring to be computed.
   * It will also check any holes, if they have been assigned.
   */
  bool containsPoint(Coordinate p) {
    LinearRing shell = getLinearRing()!;
    Envelope env = shell.getEnvelopeInternal();
    if (!env.containsCoordinate(p)) return false;
    if (!PointLocation.isInRing(p, shell.getCoordinates())) return false;

    for (Iterator i = holes.iterator; i.moveNext();) {
      EdgeRing hole = i.current as EdgeRing;
      if (hole.containsPoint(p)) return false;
    }
    return true;
  }
}

/**
 * @version 1.7
 */
class DirectedEdge extends EdgeEnd {
  /**
   * Computes the factor for the change in depth when moving from one location to another.
   * E.g. if crossing from the INTERIOR to the EXTERIOR the depth decreases, so the factor is -1
   */
  static int depthFactor(int currLocation, int nextLocation) {
    if (currLocation == Location.EXTERIOR && nextLocation == Location.INTERIOR)
      return 1;
    else if (currLocation == Location.INTERIOR &&
        nextLocation == Location.EXTERIOR) return -1;
    return 0;
  }

  bool _isForward = false;
  bool _isInResult = false;
  bool _isVisited = false;

  late DirectedEdge sym; // the symmetric edge
  late DirectedEdge
      next; // the next edge in the edge ring for the polygon containing this edge
  late DirectedEdge
      nextMin; // the next edge in the MinimalEdgeRing that contains this edge
  EdgeRing? edgeRing; // the EdgeRing that this edge is part of
  EdgeRing? minEdgeRing; // the MinimalEdgeRing that this edge is part of
  /**
   * The depth of each side (position) of this edge.
   * The 0 element of the array is never used.
   */
  List<int> depth = [0, -999, -999];

  DirectedEdge(Edge edge, bool isForward) : super(edge) {
    this._isForward = isForward;
    if (isForward) {
      init(edge.getCoordinateWithIndex(0), edge.getCoordinateWithIndex(1));
    } else {
      int n = edge.getNumPoints() - 1;
      init(edge.getCoordinateWithIndex(n), edge.getCoordinateWithIndex(n - 1));
    }
    computeDirectedLabel();
  }

  Edge getEdge() {
    return edge;
  }

  void setInResult(bool isInResult) {
    this._isInResult = isInResult;
  }

  bool isInResult() {
    return _isInResult;
  }

  bool isVisited() {
    return _isVisited;
  }

  void setVisited(bool isVisited) {
    this._isVisited = isVisited;
  }

  void setEdgeRing(EdgeRing edgeRing) {
    this.edgeRing = edgeRing;
  }

  EdgeRing? getEdgeRing() {
    return edgeRing;
  }

  void setMinEdgeRing(EdgeRing minEdgeRing) {
    this.minEdgeRing = minEdgeRing;
  }

  EdgeRing? getMinEdgeRing() {
    return minEdgeRing;
  }

  int getDepth(int position) {
    return depth[position];
  }

  void setDepth(int position, int depthVal) {
    if (depth[position] != -999) {
//      if (depth[position] != depthVal) {
//        Debug.print(this);
//      }
      if (depth[position] != depthVal) {
        throw new TopologyException(
            "assigned depths do not match ${getCoordinate()}");
      }
      //Assert.isTrue(depth[position] == depthVal, "assigned depths do not match at " + getCoordinate());
    }
    depth[position] = depthVal;
  }

  int getDepthDelta() {
    int depthDelta = edge.getDepthDelta();
    if (!_isForward) depthDelta = -depthDelta;
    return depthDelta;
  }

  /**
   * setVisitedEdge marks both DirectedEdges attached to a given Edge.
   * This is used for edges corresponding to lines, which will only
   * appear oriented in a single direction in the result.
   */
  void setVisitedEdge(bool isVisited) {
    setVisited(isVisited);
    sym.setVisited(isVisited);
  }

  /**
   * Each Edge gives rise to a pair of symmetric DirectedEdges, in opposite
   * directions.
   * @return the DirectedEdge for the same Edge but in the opposite direction
   */
  DirectedEdge getSym() {
    return sym;
  }

  bool isForward() {
    return _isForward;
  }

  void setSym(DirectedEdge de) {
    sym = de;
  }

  DirectedEdge getNext() {
    return next;
  }

  void setNext(DirectedEdge next) {
    this.next = next;
  }

  DirectedEdge getNextMin() {
    return nextMin;
  }

  void setNextMin(DirectedEdge nextMin) {
    this.nextMin = nextMin;
  }

  /**
   * This edge is a line edge if
   * <ul>
   * <li> at least one of the labels is a line label
   * <li> any labels which are not line labels have all Locations = EXTERIOR
   * </ul>
   */
  bool isLineEdge() {
    bool isLine = label!.isLine(0) || label!.isLine(1);
    bool isExteriorIfArea0 = !label!.isAreaWithIndex(0) ||
        label!.allPositionsEqual(0, Location.EXTERIOR);
    bool isExteriorIfArea1 = !label!.isAreaWithIndex(1) ||
        label!.allPositionsEqual(1, Location.EXTERIOR);

    return isLine && isExteriorIfArea0 && isExteriorIfArea1;
  }

  /**
   * This is an interior Area edge if
   * <ul>
   * <li> its label is an Area label for both Geometries
   * <li> and for each Geometry both sides are in the interior.
   * </ul>
   *
   * @return true if this is an interior Area edge
   */
  bool isInteriorAreaEdge() {
    bool isInteriorAreaEdge = true;
    for (int i = 0; i < 2; i++) {
      if (!(label!.isAreaWithIndex(i) &&
          label!.getLocationWithPosIndex(i, Position.LEFT) ==
              Location.INTERIOR &&
          label!.getLocationWithPosIndex(i, Position.RIGHT) ==
              Location.INTERIOR)) {
        isInteriorAreaEdge = false;
      }
    }
    return isInteriorAreaEdge;
  }

  /**
   * Compute the label in the appropriate orientation for this DirEdge
   */
  void computeDirectedLabel() {
    label = new Label.fromLabel(edge.getLabel()!);
    if (!_isForward) label!.flip();
  }

  /**
   * Set both edge depths.  One depth for a given side is provided.  The other is
   * computed depending on the Location transition and the depthDelta of the edge.
   */
  void setEdgeDepths(int position, int depth) {
    // get the depth transition delta from R to L for this directed Edge
    int depthDelta = getEdge().getDepthDelta();
    if (!_isForward) depthDelta = -depthDelta;

    // if moving from L to R instead of R to L must change sign of delta
    int directionFactor = 1;
    if (position == Position.LEFT) directionFactor = -1;

    int oppositePos = Position.opposite(position);
    int delta = depthDelta * directionFactor;
    //TESTINGint delta = depthDelta * DirectedEdge.depthFactor(loc, oppositeLoc);
    int oppositeDepth = depth + delta;
    setDepth(position, depth);
    setDepth(oppositePos, oppositeDepth);
  }

//   void print(PrintStream out)
//  {
//    super.print(out);
//    out.print(" " + depth[Position.LEFT] + "/" + depth[Position.RIGHT]);
//    out.print(" (" + getDepthDelta() + ")");
//    //out.print(" " + this.hashCode());
//    //if (next != null) out.print(" next:" + next.hashCode());
//    if (_isInResult) out.print(" inResult");
//  }
//   void printEdge(PrintStream out)
//  {
//    print(out);
//    out.print(" ");
//    if (_isForward)
//      edge.print(out);
//    else
//      edge.printReverse(out);
//  }

}

/**
 * Utility functions for working with quadrants, which are numbered as follows:
 * <pre>
 * 1 | 0
 * --+--
 * 2 | 3
 * </pre>
 *
 * @version 1.7
 */
class Quadrant {
  static final int NE = 0;
  static final int NW = 1;
  static final int SW = 2;
  static final int SE = 3;

  /**
   * Returns the quadrant of a directed line segment (specified as x and y
   * displacements, which cannot both be 0).
   *
   * @throws IllegalArgumentException if the displacements are both 0
   */
  static int quadrant(double dx, double dy) {
    if (dx == 0.0 && dy == 0.0) {
      throw ArgumentError("Cannot compute the quadrant for point ( $dx, $dy )");
    }
    if (dx >= 0.0) {
      if (dy >= 0.0)
        return NE;
      else
        return SE;
    } else {
      if (dy >= 0.0)
        return NW;
      else
        return SW;
    }
  }

  /**
   * Returns the quadrant of a directed line segment from p0 to p1.
   *
   * @throws IllegalArgumentException if the points are equal
   */
  static int quadrantFromCoords(Coordinate p0, Coordinate p1) {
    if (p1.x == p0.x && p1.y == p0.y)
      throw ArgumentError(
          "Cannot compute the quadrant for two identical points $p0");

    if (p1.x >= p0.x) {
      if (p1.y >= p0.y)
        return NE;
      else
        return SE;
    } else {
      if (p1.y >= p0.y)
        return NW;
      else
        return SW;
    }
  }

  /**
   * Returns true if the quadrants are 1 and 3, or 2 and 4
   */
  static bool isOpposite(int quad1, int quad2) {
    if (quad1 == quad2) return false;
    int diff = (quad1 - quad2 + 4) % 4;
    // if quadrants are not adjacent, they are opposite
    if (diff == 2) return true;
    return false;
  }

  /**
   * Returns the right-hand quadrant of the halfplane defined by the two quadrants,
   * or -1 if the quadrants are opposite, or the quadrant if they are identical.
   */
  static int commonHalfPlane(int quad1, int quad2) {
    // if quadrants are the same they do not determine a unique common halfplane.
    // Simply return one of the two possibilities
    if (quad1 == quad2) return quad1;
    int diff = (quad1 - quad2 + 4) % 4;
    // if quadrants are not adjacent, they do not share a common halfplane
    if (diff == 2) return -1;
    //
    int min = (quad1 < quad2) ? quad1 : quad2;
    int max = (quad1 > quad2) ? quad1 : quad2;
    // for this one case, the righthand plane is NOT the minimum index;
    if (min == 0 && max == 3) return 3;
    // in general, the halfplane index is the minimum of the two adjacent quadrants
    return min;
  }

  /**
   * Returns whether the given quadrant lies within the given halfplane (specified
   * by its right-hand quadrant).
   */
  static bool isInHalfPlane(int quad, int halfPlane) {
    if (halfPlane == SE) {
      return quad == SE || quad == SW;
    }
    return quad == halfPlane || quad == halfPlane + 1;
  }

  /**
   * Returns true if the given quadrant is 0 or 1.
   */
  static bool isNorthern(int quad) {
    return quad == NE || quad == NW;
  }
}

/**
 * Models the end of an edge incident on a node.
 * EdgeEnds have a direction
 * determined by the direction of the ray from the initial
 * point to the next point.
 * EdgeEnds are comparable under the ordering
 * "a has a greater angle with the x-axis than b".
 * This ordering is used to sort EdgeEnds around a node.
 * @version 1.7
 */
class EdgeEnd implements Comparable {
  Edge edge; // the parent edge of this edge end
  Label? label;

  Node? node; // the node this edge end originates at
  Coordinate? p0, p1; // points of initial line segment
  double dx = 0,
      dy = 0; // the direction vector for this edge from its starting point
  int quadrant = 0;

  EdgeEnd(this.edge);

  EdgeEnd.withCoords(Edge edge, Coordinate p0, Coordinate p1)
      : this.withCoordsLabel(edge, p0, p1, null);

  EdgeEnd.withCoordsLabel(
      this.edge, Coordinate p0, Coordinate p1, Label? label) {
    init(p0, p1);
    this.label = label;
  }

  void init(Coordinate p0, Coordinate p1) {
    this.p0 = p0;
    this.p1 = p1;
    dx = p1.x - p0.x;
    dy = p1.y - p0.y;
    quadrant = Quadrant.quadrant(dx, dy);
    // TODO Assert.isTrue(! (dx == 0 && dy == 0), "EdgeEnd with identical endpoints found");
  }

  Edge getEdge() {
    return edge;
  }

  Label? getLabel() {
    return label;
  }

  Coordinate? getCoordinate() {
    return p0;
  }

  Coordinate? getDirectedCoordinate() {
    return p1;
  }

  int getQuadrant() {
    return quadrant;
  }

  double getDx() {
    return dx;
  }

  double getDy() {
    return dy;
  }

  void setNode(Node node) {
    this.node = node;
  }

  Node? getNode() {
    return node;
  }

  int compareTo(dynamic obj) {
    EdgeEnd e = obj;
    return compareDirection(e);
  }

  /**
   * Implements the total order relation:
   * <p>
   *    a has a greater angle with the positive x-axis than b
   * <p>
   * Using the obvious algorithm of simply computing the angle is not robust,
   * since the angle calculation is obviously susceptible to roundoff.
   * A robust algorithm is:
   * - first compare the quadrant.  If the quadrants
   * are different, it it trivial to determine which vector is "greater".
   * - if the vectors lie in the same quadrant, the computeOrientation function
   * can be used to decide the relative orientation of the vectors.
   */
  int compareDirection(EdgeEnd e) {
    if (dx == e.dx && dy == e.dy) return 0;
    // if the rays are in different quadrants, determining the ordering is trivial
    if (quadrant > e.quadrant) return 1;
    if (quadrant < e.quadrant) return -1;
    // vectors are in the same quadrant - check relative orientation of direction vectors
    // this is > e if it is CCW of e
    return Orientation.index(e.p0!, e.p1!, p1!);
  }

  void computeLabel(BoundaryNodeRule boundaryNodeRule) {
    // subclasses should override this if they are using labels
  }

//   void print(PrintStream out)
//  {
//    double angle = Math.atan2(dy, dx);
//    String className = getClass().getName();
//    int lastDotPos = className.lastIndexOf('.');
//    String name = className.substring(lastDotPos + 1);
//    out.print("  " + name + ": " + p0 + " - " + p1 + " " + quadrant + ":" + angle + "   " + label);
//  }
  String toString() {
    double angle = math.atan2(dy, dx);
    String className = runtimeType.toString();
    int lastDotPos = className.lastIndexOf('.');
    String name = className.substring(lastDotPos + 1);
    return "  $name: $p0  - $p1 $quadrant: $angle    $label";
  }
}

/**
 * A EdgeEndStar is an ordered list of EdgeEnds around a node.
 * They are maintained in CCW order (starting with the positive x-axis) around the node
 * for efficient lookup and topology building.
 *
 * @version 1.7
 */
abstract class EdgeEndStar {
  /**
   * A map which maintains the edges in sorted order around the node
   */
  Map edgeMap = new SplayTreeMap();

  /**
   * A list of all outgoing edges in the result, in CCW order
   */
  List? edgeList;

  /**
   * The location of the point for this star in Geometry i Areas
   */
  List<int> ptInAreaLocation = [Location.NONE, Location.NONE];

  EdgeEndStar() {}

  /**
   * Insert a EdgeEnd into this EdgeEndStar
   */
  void insert(EdgeEnd e);

  /**
   * Insert an EdgeEnd into the map, and clear the edgeList cache,
   * since the list of edges has now changed
   */
  void insertEdgeEnd(EdgeEnd e, Object obj) {
    edgeMap[e] = obj;
    edgeList = null; // edge list has changed - clear the cache
  }

  /**
   * @return the coordinate for the node this star is based at
   */
  Coordinate? getCoordinate() {
    Iterator it = iterator();
    if (!it.moveNext()) return null;
    EdgeEnd e = it.current;
    return e.getCoordinate();
  }

  int getDegree() {
    return edgeMap.length;
  }

  /**
   * Iterator access to the ordered list of edges is optimized by
   * copying the map collection to a list.  (This assumes that
   * once an iterator is requested, it is likely that insertion into
   * the map is complete).
   */
  Iterator iterator() {
    return getEdges().iterator;
  }

  List getEdges() {
    if (edgeList == null) {
      edgeList = List.from(edgeMap.values);
    }
    return edgeList!;
  }

  EdgeEnd getNextCW(EdgeEnd ee) {
    getEdges();
    int i = edgeList!.indexOf(ee);
    int iNextCW = i - 1;
    if (i == 0) iNextCW = edgeList!.length - 1;
    return edgeList![iNextCW];
  }

  void computeLabelling(List<GeometryGraph> geomGraph) {
    computeEdgeEndLabels(geomGraph[0].getBoundaryNodeRule());
    // Propagate side labels  around the edges in the star
    // for each parent Geometry
//Debug.print(this);
    propagateSideLabels(0);
//Debug.print(this);
//Debug.printIfWatch(this);
    propagateSideLabels(1);
//Debug.print(this);
//Debug.printIfWatch(this);

    /**
     * If there are edges that still have null labels for a geometry
     * this must be because there are no area edges for that geometry incident on this node.
     * In this case, to label the edge for that geometry we must test whether the
     * edge is in the interior of the geometry.
     * To do this it suffices to determine whether the node for the edge is in the interior of an area.
     * If so, the edge has location INTERIOR for the geometry.
     * In all other cases (e.g. the node is on a line, on a point, or not on the geometry at all) the edge
     * has the location EXTERIOR for the geometry.
     * <p>
     * Note that the edge cannot be on the BOUNDARY of the geometry, since then
     * there would have been a parallel edge from the Geometry at this node also labelled BOUNDARY
     * and this edge would have been labelled in the previous step.
     * <p>
     * This code causes a problem when dimensional collapses are present, since it may try and
     * determine the location of a node where a dimensional collapse has occurred.
     * The point should be considered to be on the EXTERIOR
     * of the polygon, but locate() will return INTERIOR, since it is passed
     * the original Geometry, not the collapsed version.
     *
     * If there are incident edges which are Line edges labelled BOUNDARY,
     * then they must be edges resulting from dimensional collapses.
     * In this case the other edges can be labelled EXTERIOR for this Geometry.
     *
     * MD 8/11/01 - NOT TRUE!  The collapsed edges may in fact be in the interior of the Geometry,
     * which means the other edges should be labelled INTERIOR for this Geometry.
     * Not sure how solve this...  Possibly labelling needs to be split into several phases:
     * area label propagation, symLabel merging, then finally null label resolution.
     */
    List<bool> hasDimensionalCollapseEdge = [false, false];
    for (Iterator it = iterator(); it.moveNext();) {
      EdgeEnd e = it.current;
      Label label = e.getLabel()!;
      for (int geomi = 0; geomi < 2; geomi++) {
        if (label.isLine(geomi) &&
            label.getLocation(geomi) == Location.BOUNDARY)
          hasDimensionalCollapseEdge[geomi] = true;
      }
    }
//Debug.print(this);
    for (Iterator it = iterator(); it.moveNext();) {
      EdgeEnd e = it.current;
      Label label = e.getLabel()!;
//Debug.println(e);
      for (int geomi = 0; geomi < 2; geomi++) {
        if (label.isAnyNull(geomi)) {
          int loc = Location.NONE;
          if (hasDimensionalCollapseEdge[geomi]) {
            loc = Location.EXTERIOR;
          } else {
            Coordinate p = e.getCoordinate()!;
            loc = getLocation(geomi, p, geomGraph);
          }
          label.setAllLocationsIfNullWithIndex(geomi, loc);
        }
      }
//Debug.println(e);
    }
//Debug.print(this);
//Debug.printIfWatch(this);
  }

  void computeEdgeEndLabels(BoundaryNodeRule boundaryNodeRule) {
    // Compute edge label for each EdgeEnd
    for (Iterator it = iterator(); it.moveNext();) {
      EdgeEnd ee = it.current;
      ee.computeLabel(boundaryNodeRule);
    }
  }

  int getLocation(int geomIndex, Coordinate p, List<GeometryGraph> geom) {
    // compute location only on demand
    if (ptInAreaLocation[geomIndex] == Location.NONE) {
      ptInAreaLocation[geomIndex] = SimplePointInAreaLocator.locatePointInGeom(
          p, geom[geomIndex].getGeometry()!);
    }
    return ptInAreaLocation[geomIndex];
  }

  bool isAreaLabelsConsistent(GeometryGraph geomGraph) {
    computeEdgeEndLabels(geomGraph.getBoundaryNodeRule());
    return checkAreaLabelsConsistent(0);
  }

  bool checkAreaLabelsConsistent(int geomIndex) {
    // Since edges are stored in CCW order around the node,
    // As we move around the ring we move from the right to the left side of the edge
    List edges = getEdges();
    // if no edges, trivially consistent
    if (edges.length <= 0) return true;
    // initialize startLoc to location of last L side (if any)
    int lastEdgeIndex = edges.length - 1;
    Label startLabel = (edges[lastEdgeIndex] as EdgeEnd).getLabel()!;
    int startLoc = startLabel.getLocationWithPosIndex(geomIndex, Position.LEFT);
    // TODO Assert.isTrue(startLoc != Location.NONE, "Found unlabelled area edge");

    int currLoc = startLoc;
    for (Iterator it = iterator(); it.moveNext();) {
      EdgeEnd e = it.current;
      Label label = e.getLabel()!;
      // we assume that we are only checking a area
      // TODO Assert.isTrue(label.isArea(geomIndex), "Found non-area edge");
      int leftLoc = label.getLocationWithPosIndex(geomIndex, Position.LEFT);
      int rightLoc = label.getLocationWithPosIndex(geomIndex, Position.RIGHT);
//System.out.println(leftLoc + " " + rightLoc);
//Debug.print(this);
      // check that edge is really a boundary between inside and outside!
      if (leftLoc == rightLoc) {
        return false;
      }
      // check side location conflict
      //Assert.isTrue(rightLoc == currLoc, "side location conflict " + locStr);
      if (rightLoc != currLoc) {
//Debug.print(this);
        return false;
      }
      currLoc = leftLoc;
    }
    return true;
  }

  void propagateSideLabels(int geomIndex) {
    // Since edges are stored in CCW order around the node,
    // As we move around the ring we move from the right to the left side of the edge
    int startLoc = Location.NONE;

    // initialize loc to location of last L side (if any)
//System.out.println("finding start location");
    for (Iterator it = iterator(); it.moveNext();) {
      EdgeEnd e = it.current;
      Label label = e.getLabel()!;
      if (label.isAreaWithIndex(geomIndex) &&
          label.getLocationWithPosIndex(geomIndex, Position.LEFT) !=
              Location.NONE)
        startLoc = label.getLocationWithPosIndex(geomIndex, Position.LEFT);
    }

    // no labelled sides found, so no labels to propagate
    if (startLoc == Location.NONE) return;

    int currLoc = startLoc;
    for (Iterator it = iterator(); it.moveNext();) {
      EdgeEnd e = it.current;
      Label label = e.getLabel()!;
      // set null ON values to be in current location
      if (label.getLocationWithPosIndex(geomIndex, Position.ON) ==
          Location.NONE) label.setLocation(geomIndex, Position.ON, currLoc);
      // set side labels (if any)
      if (label.isAreaWithIndex(geomIndex)) {
        int leftLoc = label.getLocationWithPosIndex(geomIndex, Position.LEFT);
        int rightLoc = label.getLocationWithPosIndex(geomIndex, Position.RIGHT);
        // if there is a right location, that is the next location to propagate
        if (rightLoc != Location.NONE) {
//Debug.print(rightLoc != currLoc, this);
          if (rightLoc != currLoc)
            throw new TopologyException(
                "side location conflict ${e.getCoordinate()}");
          if (leftLoc == Location.NONE) {
            Assert.shouldNeverReachHere(
                "found single null side (at ${e.getCoordinate()})");
          }
          currLoc = leftLoc;
        } else {
          /** RHS is null - LHS must be null too.
           *  This must be an edge from the other geometry, which has no location
           *  labelling for this geometry.  This edge must lie wholly inside or outside
           *  the other geometry (which is determined by the current location).
           *  Assign both sides to be the current location.
           */
          assert(
              label.getLocationWithPosIndex(geomIndex, Position.LEFT) ==
                  Location.NONE,
              "found single null side");
          label.setLocation(geomIndex, Position.RIGHT, currLoc);
          label.setLocation(geomIndex, Position.LEFT, currLoc);
        }
      }
    }
  }

  int findIndex(EdgeEnd eSearch) {
    iterator; // force edgelist to be computed
    for (int i = 0; i < edgeList!.length; i++) {
      EdgeEnd e = edgeList![i];
      if (e == eSearch) return i;
    }
    return -1;
  }

//   void print(PrintStream out)
//  {
//    System.out.println("EdgeEndStar:   " + getCoordinate());
//    for (Iterator it = iterator; it.moveNext(); ) {
//      EdgeEnd e = (EdgeEnd) it.current;
//      e.print(out);
//    }
//  }

  String toString() {
    StringBuffer buf = new StringBuffer();
    buf.write("EdgeEndStar:   ${getCoordinate()}");
    buf.write("\n");
    for (Iterator it = iterator(); it.moveNext();) {
      EdgeEnd e = it.current;
      buf.write(e);
      buf.write("\n");
    }
    return buf.toString();
  }
}

/**
 * A list of edge intersections along an {@link Edge}.
 * Implements splitting an edge with intersections
 * into multiple resultant edges.
 *
 * @version 1.7
 */
class EdgeIntersectionList {
  // a Map <EdgeIntersection, EdgeIntersection>
  Map nodeMap = new SplayTreeMap();
  Edge edge; // the parent edge

  EdgeIntersectionList(this.edge);

  /**
   * Adds an intersection into the list, if it isn't already there.
   * The input segmentIndex and dist are expected to be normalized.
   * @return the EdgeIntersection found or added
   */
  EdgeIntersection add(Coordinate intPt, int segmentIndex, double dist) {
    EdgeIntersection eiNew = new EdgeIntersection(intPt, segmentIndex, dist);
    EdgeIntersection? ei = nodeMap[eiNew];
    if (ei != null) {
      return ei;
    }
    nodeMap[eiNew] = eiNew;
    return eiNew;
  }

  /**
   * Returns an iterator of {@link EdgeIntersection}s
   *
   * @return an Iterator of EdgeIntersections
   */
  Iterator iterator() {
    return nodeMap.values.iterator;
  }

  /**
   * Tests if the given point is an edge intersection
   *
   * @param pt the point to test
   * @return true if the point is an intersection
   */
  bool isIntersection(Coordinate pt) {
    for (Iterator it = iterator(); it.moveNext();) {
      EdgeIntersection ei = it.current;
      if (ei.coord.equals(pt)) return true;
    }
    return false;
  }

  /**
   * Adds entries for the first and last points of the edge to the list
   */
  void addEndpoints() {
    int maxSegIndex = edge.pts.length - 1;
    add(edge.pts[0], 0, 0.0);
    add(edge.pts[maxSegIndex], maxSegIndex, 0.0);
  }

  /**
   * Creates new edges for all the edges that the intersections in this
   * list split the parent edge into.
   * Adds the edges to the input list (this is so a single list
   * can be used to accumulate all split edges for a Geometry).
   *
   * @param edgeList a list of EdgeIntersections
   */
  void addSplitEdges(List edgeList) {
    // ensure that the list has entries for the first and last point of the edge
    addEndpoints();

    Iterator it = iterator();
    // there should always be at least two entries in the list
    it.moveNext();
    EdgeIntersection eiPrev = it.current;
    while (it.moveNext()) {
      EdgeIntersection ei = it.current;
      Edge newEdge = createSplitEdge(eiPrev, ei);
      edgeList.add(newEdge);

      eiPrev = ei;
    }
  }

  /**
   * Create a new "split edge" with the section of points between
   * (and including) the two intersections.
   * The label for the new edge is the same as the label for the parent edge.
   */
  Edge createSplitEdge(EdgeIntersection ei0, EdgeIntersection ei1) {
//Debug.print("\ncreateSplitEdge"); Debug.print(ei0); Debug.print(ei1);
    int npts = ei1.segmentIndex - ei0.segmentIndex + 2;

    Coordinate lastSegStartPt = edge.pts[ei1.segmentIndex];
    // if the last intersection point is not equal to the its segment start pt,
    // add it to the points list as well.
    // (This check is needed because the distance metric is not totally reliable!)
    // The check for point equality is 2D only - Z values are ignored
    bool useIntPt1 = ei1.dist > 0.0 || !ei1.coord.equals2D(lastSegStartPt);
    if (!useIntPt1) {
      npts--;
    }

    List<Coordinate> pts = []; //..length = (npts);
    // int ipt = 0;
    pts.add(Coordinate.fromCoordinate(ei0.coord));
    // pts[ipt++] = new Coordinate.fromCoordinate(ei0.coord);
    for (int i = ei0.segmentIndex + 1; i <= ei1.segmentIndex; i++) {
      pts.add(edge.pts[i]);
      // pts[ipt++] = edge.pts[i];
    }
    if (useIntPt1) pts.add(ei1.coord);
    // if (useIntPt1) pts[ipt] = ei1.coord;
    return new Edge(pts, new Label.fromLabel(edge.label!));
  }

//   void print(PrintStream out)
//  {
//    out.println("Intersections:");
//    for (Iterator it = iterator; it.moveNext(); ) {
//      EdgeIntersection ei = (EdgeIntersection) it.current;
//      ei.print(out);
//    }
//  }
}

/**
 * A Depth object records the topological depth of the sides
 * of an Edge for up to two Geometries.
 * @version 1.7
 */
class Depth {
  static final int NULL_VALUE = -1;

  static int depthAtLocation(int location) {
    if (location == Location.EXTERIOR) return 0;
    if (location == Location.INTERIOR) return 1;
    return NULL_VALUE;
  }

  List<List<int>> depth = MatrixUtils.createIntMatrixWithDefault(2, 3, 0);

  Depth() {
    // initialize depth array to a sentinel value
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 3; j++) {
        depth[i][j] = NULL_VALUE;
      }
    }
  }

  int getDepth(int geomIndex, int posIndex) {
    return depth[geomIndex][posIndex];
  }

  void setDepth(int geomIndex, int posIndex, int depthValue) {
    depth[geomIndex][posIndex] = depthValue;
  }

  int getLocation(int geomIndex, int posIndex) {
    if (depth[geomIndex][posIndex] <= 0) return Location.EXTERIOR;
    return Location.INTERIOR;
  }

  void add3(int geomIndex, int posIndex, int location) {
    if (location == Location.INTERIOR) depth[geomIndex][posIndex]++;
  }

  /**
   * A Depth object is null (has never been initialized) if all depths are null.
   */
  bool isNull() {
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 3; j++) {
        if (depth[i][j] != NULL_VALUE) return false;
      }
    }
    return true;
  }

  bool isNull1(int geomIndex) {
    return depth[geomIndex][1] == NULL_VALUE;
  }

  bool isNull2(int geomIndex, int posIndex) {
    return depth[geomIndex][posIndex] == NULL_VALUE;
  }

  void add(Label lbl) {
    for (int i = 0; i < 2; i++) {
      for (int j = 1; j < 3; j++) {
        int loc = lbl.getLocationWithPosIndex(i, j);
        if (loc == Location.EXTERIOR || loc == Location.INTERIOR) {
          // initialize depth if it is null, otherwise add this location value
          if (isNull2(i, j)) {
            depth[i][j] = depthAtLocation(loc);
          } else
            depth[i][j] += depthAtLocation(loc);
        }
      }
    }
  }

  int getDelta(int geomIndex) {
    return depth[geomIndex][Position.RIGHT] - depth[geomIndex][Position.LEFT];
  }

  /**
   * Normalize the depths for each geometry, if they are non-null.
   * A normalized depth
   * has depth values in the set { 0, 1 }.
   * Normalizing the depths
   * involves reducing the depths by the same amount so that at least
   * one of them is 0.  If the remaining value is &gt; 0, it is set to 1.
   */
  void normalize() {
    for (int i = 0; i < 2; i++) {
      if (!isNull1(i)) {
        int minDepth = depth[i][1];
        if (depth[i][2] < minDepth) minDepth = depth[i][2];

        if (minDepth < 0) minDepth = 0;
        for (int j = 1; j < 3; j++) {
          int newValue = 0;
          if (depth[i][j] > minDepth) newValue = 1;
          depth[i][j] = newValue;
        }
      }
    }
  }

  String toString() {
    return "A: ${depth[0][1]} ,${depth[0][2]} B: ${depth[1][1]},${depth[1][2]}";
  }
}

/**
 * @version 1.7
 */
class Node extends GraphComponent {
  Coordinate coord; // only non-null if this node is precise
  EdgeEndStar? edges;

  Node(this.coord, EdgeEndStar? edges) {
    this.edges = edges;
    label = new Label.args2(0, Location.NONE);
  }

  Coordinate getCoordinate() {
    return coord;
  }

  EdgeEndStar? getEdges() {
    return edges;
  }

  /**
   * Tests whether any incident edge is flagged as
   * being in the result.
   * This test can be used to determine if the node is in the result,
   * since if any incident edge is in the result, the node must be in the result as well.
   *
   * @return <code>true</code> if any incident edge in the in the result
   */
  bool isIncidentEdgeInResult() {
    for (Iterator it = getEdges()!.getEdges().iterator; it.moveNext();) {
      DirectedEdge de = it.current;
      if (de.getEdge().isInResult()) return true;
    }
    return false;
  }

  bool isIsolated() {
    return (label!.getGeometryCount() == 1);
  }

  /**
   * Basic nodes do not compute IMs
   */
  void computeIM(IntersectionMatrix im) {}

  /**
   * Add the edge to the list of edges at this node
   */
  void add(EdgeEnd e) {
    // Assert: start pt of e is equal to node point
    edges!.insert(e);
    e.setNode(this);
  }

  void mergeLabelFromNode(Node n) {
    mergeLabel(n.label!);
  }

  /**
   * To merge labels for two nodes,
   * the merged location for each LabelElement is computed.
   * The location for the corresponding node LabelElement is set to the result,
   * as long as the location is non-null.
   */

  void mergeLabel(Label label2) {
    for (int i = 0; i < 2; i++) {
      int loc = computeMergedLocation(label2, i);
      int thisLoc = label!.getLocation(i);
      if (thisLoc == Location.NONE) label!.setLocationWithIndex(i, loc);
    }
  }

  void setLabelWithIndex(int argIndex, int onLocation) {
    if (label == null) {
      label = new Label.args2(argIndex, onLocation);
    } else
      label!.setLocationWithIndex(argIndex, onLocation);
  }

  /**
   * Updates the label of a node to BOUNDARY,
   * obeying the mod-2 boundaryDetermination rule.
   */
  void setLabelBoundary(int argIndex) {
    if (label == null) return;

    // determine the current location for the point (if any)
    int loc = Location.NONE;
    if (label != null) loc = label!.getLocation(argIndex);
    // flip the loc
    int newLoc;
    switch (loc) {
      case Location.BOUNDARY:
        newLoc = Location.INTERIOR;
        break;
      case Location.INTERIOR:
        newLoc = Location.BOUNDARY;
        break;
      default:
        newLoc = Location.BOUNDARY;
        break;
    }
    label!.setLocationWithIndex(argIndex, newLoc);
  }

  /**
   * The location for a given eltIndex for a node will be one
   * of { null, INTERIOR, BOUNDARY }.
   * A node may be on both the boundary and the interior of a geometry;
   * in this case, the rule is that the node is considered to be in the boundary.
   * The merged location is the maximum of the two input values.
   */
  int computeMergedLocation(Label label2, int eltIndex) {
    int loc = Location.NONE;
    loc = label!.getLocation(eltIndex);
    if (!label2.isNull(eltIndex)) {
      int nLoc = label2.getLocation(eltIndex);
      if (loc != Location.BOUNDARY) loc = nLoc;
    }
    return loc;
  }

//   void print(PrintStream out)
//  {
//    out.println("node " + coord + " lbl: " + label);
//  }
}

/**
 * @version 1.7
 */
class Edge extends GraphComponent {
  /**
   * Updates an IM from the label for an edge.
   * Handles edges from both L and A geometries.
   */
  static void updateIMStatic(Label label, IntersectionMatrix im) {
    im.setAtLeastIfValid(label.getLocationWithPosIndex(0, Position.ON),
        label.getLocationWithPosIndex(1, Position.ON), 1);
    if (label.isArea()) {
      im.setAtLeastIfValid(label.getLocationWithPosIndex(0, Position.LEFT),
          label.getLocationWithPosIndex(1, Position.LEFT), 2);
      im.setAtLeastIfValid(label.getLocationWithPosIndex(0, Position.RIGHT),
          label.getLocationWithPosIndex(1, Position.RIGHT), 2);
    }
  }

  List<Coordinate> pts;
  Envelope? env;
  late EdgeIntersectionList eiList;
  String? name;
  MonotoneChainEdge? mce;
  bool _isIsolated = true;
  Depth depth = new Depth();
  int depthDelta =
      0; // the change in area depth from the R to L side of this edge

  Edge(this.pts, Label? label) {
    eiList = new EdgeIntersectionList(this);
    this.label = label;
  }

  Edge.fromList(List<Coordinate> pts) : this(pts, null);

  int getNumPoints() {
    return pts.length;
  }

  void setName(String name) {
    this.name = name;
  }

  List<Coordinate> getCoordinates() {
    return pts;
  }

  Coordinate getCoordinateWithIndex(int i) {
    return pts[i];
  }

  Coordinate? getCoordinate() {
    if (pts.length > 0) return pts[0];
    return null;
  }

  Envelope getEnvelope() {
    // compute envelope lazily
    if (env == null) {
      env = new Envelope.empty();
      for (int i = 0; i < pts.length; i++) {
        env!.expandToIncludeCoordinate(pts[i]);
      }
    }
    return env!;
  }

  Depth getDepth() {
    return depth;
  }

  /**
   * The depthDelta is the change in depth as an edge is crossed from R to L
   * @return the change in depth as the edge is crossed from R to L
   */
  int getDepthDelta() {
    return depthDelta;
  }

  void setDepthDelta(int depthDelta) {
    this.depthDelta = depthDelta;
  }

  int getMaximumSegmentIndex() {
    return pts.length - 1;
  }

  EdgeIntersectionList getEdgeIntersectionList() {
    return eiList;
  }

  MonotoneChainEdge getMonotoneChainEdge() {
    if (mce == null) mce = new MonotoneChainEdge(this);
    return mce!;
  }

  bool isClosed() {
    return pts[0].equals(pts[pts.length - 1]);
  }

  /**
   * An Edge is collapsed if it is an Area edge and it consists of
   * two segments which are equal and opposite (eg a zero-width V).
   */
  bool isCollapsed() {
    if (!label!.isArea()) return false;
    if (pts.length != 3) return false;
    if (pts[0].equals(pts[2])) return true;
    return false;
  }

  Edge getCollapsedEdge() {
    List<Coordinate> newPts = [];
    newPts[0] = pts[0];
    newPts[1] = pts[1];
    Edge newe = new Edge(newPts, Label.toLineLabel(label!));
    return newe;
  }

  void setIsolated(bool isIsolated) {
    this._isIsolated = isIsolated;
  }

  bool isIsolated() {
    return _isIsolated;
  }

  /**
   * Adds EdgeIntersections for one or both
   * intersections found for a segment of an edge to the edge intersection list.
   */
  void addIntersections(LineIntersector li, int segmentIndex, int geomIndex) {
    for (int i = 0; i < li.getIntersectionNum(); i++) {
      addIntersection(li, segmentIndex, geomIndex, i);
    }
  }

  /**
   * Add an EdgeIntersection for intersection intIndex.
   * An intersection that falls exactly on a vertex of the edge is normalized
   * to use the higher of the two possible segmentIndexes
   */
  void addIntersection(
      LineIntersector li, int segmentIndex, int geomIndex, int intIndex) {
    Coordinate intPt =
        new Coordinate.fromCoordinate(li.getIntersection(intIndex));
    int normalizedSegmentIndex = segmentIndex;
    double dist = li.getEdgeDistance(geomIndex, intIndex);
//Debug.println("edge intpt: " + intPt + " dist: " + dist);
    // normalize the intersection point location
    int nextSegIndex = normalizedSegmentIndex + 1;
    if (nextSegIndex < pts.length) {
      Coordinate nextPt = pts[nextSegIndex];
//Debug.println("next pt: " + nextPt);

      // Normalize segment index if intPt falls on vertex
      // The check for point equality is 2D only - Z values are ignored
      if (intPt.equals2D(nextPt)) {
//Debug.println("normalized distance");
        normalizedSegmentIndex = nextSegIndex;
        dist = 0.0;
      }
    }
    /**
     * Add the intersection point to edge intersection list.
     */
    EdgeIntersection ei = eiList.add(intPt, normalizedSegmentIndex, dist);
//ei.print(System.out);
  }

  /**
   * Update the IM with the contribution for this component.
   * A component only contributes if it has a labelling for both parent geometries
   */
  void computeIM(IntersectionMatrix im) {
    updateIMStatic(label!, im);
  }

  /**
   * equals is defined to be:
   * <p>
   * e1 equals e2
   * <b>iff</b>
   * the coordinates of e1 are the same or the reverse of the coordinates in e2
   */
  bool equals(Object o) {
    if (!(o is Edge)) return false;
    Edge e = o;

    if (pts.length != e.pts.length) return false;

    bool isEqualForward = true;
    bool isEqualReverse = true;
    int iRev = pts.length;
    for (int i = 0; i < pts.length; i++) {
      if (!pts[i].equals2D(e.pts[i])) {
        isEqualForward = false;
      }
      if (!pts[i].equals2D(e.pts[--iRev])) {
        isEqualReverse = false;
      }
      if (!isEqualForward && !isEqualReverse) return false;
    }
    return true;
  }

  /**
   * @return true if the coordinate sequences of the Edges are identical
   */
  bool isPointwiseEqual(Edge e) {
    if (pts.length != e.pts.length) return false;

    for (int i = 0; i < pts.length; i++) {
      if (!pts[i].equals2D(e.pts[i])) {
        return false;
      }
    }
    return true;
  }

  String toString() {
    StringBuffer builder = new StringBuffer();
    builder.write("edge " + (name != null ? name! : "") + ": ");
    builder.write("LINESTRING (");
    for (int i = 0; i < pts.length; i++) {
      if (i > 0) builder.write(",");
      builder.write(pts[i].x);
      builder.write(" ");
      builder.write(pts[i].y);
    }
    builder.write(")  $label $depthDelta");
    return builder.toString();
  }

//   void print(PrintStream out)
//  {
//    out.print("edge " + name + ": ");
//    out.print("LINESTRING (");
//    for (int i = 0; i < pts.length; i++) {
//      if (i > 0) out.print(",");
//      out.print(pts[i].x + " " + pts[i].y);
//    }
//    out.print(")  " + label + " " + depthDelta);
//  }
//   void printReverse(PrintStream out)
//  {
//    out.print("edge " + name + ": ");
//    for (int i = pts.length - 1; i >= 0; i--) {
//      out.print(pts[i] + " ");
//    }
//    out.println("");
//  }

}

/**
 * A TopologyLocation is the labelling of a
 * GraphComponent's topological relationship to a single Geometry.
 * <p>
 * If the parent component is an area edge, each side and the edge itself
 * have a topological location.  These locations are named
 * <ul>
 * <li> ON: on the edge
 * <li> LEFT: left-hand side of the edge
 * <li> RIGHT: right-hand side
 * </ul>
 * If the parent component is a line edge or node, there is a single
 * topological relationship attribute, ON.
 * <p>
 * The possible values of a topological location are
 * {Location.NONE, Location.EXTERIOR, Location.BOUNDARY, Location.INTERIOR}
 * <p>
 * The labelling is stored in an array location[j] where
 * where j has the values ON, LEFT, RIGHT
 * @version 1.7
 */
class TopologyLocation {
  late List<int> location;

  TopologyLocation.fromList(List<int> location) {
    init(location.length);
  }

  /**
   * Constructs a TopologyLocation specifying how points on, to the left of, and to the
   * right of some GraphComponent relate to some Geometry. Possible values for the
   * parameters are Location.NULL, Location.EXTERIOR, Location.BOUNDARY,
   * and Location.INTERIOR.
   * @see Location
   */
  TopologyLocation(int on, int left, int right) {
    init(3);
    location[Position.ON] = on;
    location[Position.LEFT] = left;
    location[Position.RIGHT] = right;
  }

  TopologyLocation.fromOn(int on) {
    init(1);
    location[Position.ON] = on;
  }

  TopologyLocation.fromTL(TopologyLocation gl) {
    init(gl.location.length);
    // if (gl != null) {
    for (int i = 0; i < location.length; i++) {
      location[i] = gl.location[i];
    }
    // }
  }

  void init(int size) {
    location = List.filled(size, 0);
    setAllLocations(Location.NONE);
  }

  int get(int posIndex) {
    if (posIndex < location.length) return location[posIndex];
    return Location.NONE;
  }

  /**
   * @return true if all locations are NULL
   */
  bool isNull() {
    for (int i = 0; i < location.length; i++) {
      if (location[i] != Location.NONE) return false;
    }
    return true;
  }

  /**
   * @return true if any locations are NULL
   */
  bool isAnyNull() {
    for (int i = 0; i < location.length; i++) {
      if (location[i] == Location.NONE) return true;
    }
    return false;
  }

  bool isEqualOnSide(TopologyLocation le, int locIndex) {
    return location[locIndex] == le.location[locIndex];
  }

  bool isArea() {
    return location.length > 1;
  }

  bool isLine() {
    return location.length == 1;
  }

  void flip() {
    if (location.length <= 1) return;
    int temp = location[Position.LEFT];
    location[Position.LEFT] = location[Position.RIGHT];
    location[Position.RIGHT] = temp;
  }

  void setAllLocations(int locValue) {
    for (int i = 0; i < location.length; i++) {
      location[i] = locValue;
    }
  }

  void setAllLocationsIfNull(int locValue) {
    for (int i = 0; i < location.length; i++) {
      if (location[i] == Location.NONE) location[i] = locValue;
    }
  }

  void setLocationWithIndex(int locIndex, int locValue) {
    location[locIndex] = locValue;
  }

  void setLocation(int locValue) {
    setLocationWithIndex(Position.ON, locValue);
  }

  List<int?> getLocations() {
    return location;
  }

  void setLocations(int on, int left, int right) {
    location[Position.ON] = on;
    location[Position.LEFT] = left;
    location[Position.RIGHT] = right;
  }

  bool allPositionsEqual(int loc) {
    for (int i = 0; i < location.length; i++) {
      if (location[i] != loc) return false;
    }
    return true;
  }

  /**
   * merge updates only the NULL attributes of this object
   * with the attributes of another.
   */
  void merge(TopologyLocation gl) {
    // if the src is an Area label & and the dest is not, increase the dest to be an Area
    if (gl.location.length > location.length) {
      List<int> newLoc = List.filled(3, 0);
      newLoc[Position.ON] = location[Position.ON];
      newLoc[Position.LEFT] = Location.NONE;
      newLoc[Position.RIGHT] = Location.NONE;
      location = newLoc;
    }
    for (int i = 0; i < location.length; i++) {
      if (location[i] == Location.NONE && i < gl.location.length)
        location[i] = gl.location[i];
    }
  }

  String toString() {
    StringBuffer buf = new StringBuffer();
    if (location.length > 1)
      buf.write(Location.toLocationSymbol(location[Position.LEFT]));
    buf.write(Location.toLocationSymbol(location[Position.ON]));
    if (location.length > 1)
      buf.write(Location.toLocationSymbol(location[Position.RIGHT]));
    return buf.toString();
  }
}

/**
 * A Position indicates the position of a Location relative to a graph component
 * (Node, Edge, or Area).
 * @version 1.7
 */
class Position {
  /** An indicator that a Location is <i>on</i> a GraphComponent */
  static final int ON = 0;

  /** An indicator that a Location is to the <i>left</i> of a GraphComponent */
  static final int LEFT = 1;

  /** An indicator that a Location is to the <i>right</i> of a GraphComponent */
  static final int RIGHT = 2;

  /**
   * Returns LEFT if the position is RIGHT, RIGHT if the position is LEFT, or the position
   * otherwise.
   */
  static int opposite(int position) {
    if (position == LEFT) return RIGHT;
    if (position == RIGHT) return LEFT;
    return position;
  }
}

/**
 * A <code>Label</code> indicates the topological relationship of a component
 * of a topology graph to a given <code>Geometry</code>.
 * This class supports labels for relationships to two <code>Geometry</code>s,
 * which is sufficient for algorithms for binary operations.
 * <P>
 * Topology graphs support the concept of labeling nodes and edges in the graph.
 * The label of a node or edge specifies its topological relationship to one or
 * more geometries.  (In fact, since JTS operations have only two arguments labels
 * are required for only two geometries).  A label for a node or edge has one or
 * two elements, depending on whether the node or edge occurs in one or both of the
 * input <code>Geometry</code>s.  Elements contain attributes which categorize the
 * topological location of the node or edge relative to the parent
 * <code>Geometry</code>; that is, whether the node or edge is in the interior,
 * boundary or exterior of the <code>Geometry</code>.  Attributes have a value
 * from the set <code>{Interior, Boundary, Exterior}</code>.  In a node each
 * element has  a single attribute <code>&lt;On&gt;</code>.  For an edge each element has a
 * triplet of attributes <code>&lt;Left, On, Right&gt;</code>.
 * <P>
 * It is up to the client code to associate the 0 and 1 <code>TopologyLocation</code>s
 * with specific geometries.
 * @version 1.7
 *
 */
class Label {
  // converts a Label to a Line label (that is, one with no side Locations)
  static Label toLineLabel(Label label) {
    Label lineLabel = new Label(Location.NONE);
    for (int i = 0; i < 2; i++) {
      lineLabel.setLocationWithIndex(i, label.getLocation(i));
    }
    return lineLabel;
  }

  List<TopologyLocation?> elt = []..length = (2);

  /**
   * Construct a Label with a single location for both Geometries.
   * Initialize the locations to Null
   */
  Label(int onLoc) {
    elt[0] = new TopologyLocation.fromOn(onLoc);
    elt[1] = new TopologyLocation.fromOn(onLoc);
  }

  /**
   * Construct a Label with a single location for both Geometries.
   * Initialize the location for the Geometry index.
   */
  Label.args2(int geomIndex, int onLoc) {
    elt[0] = new TopologyLocation.fromOn(Location.NONE);
    elt[1] = new TopologyLocation.fromOn(Location.NONE);
    elt[geomIndex]!.setLocation(onLoc);
  }

  /**
   * Construct a Label with On, Left and Right locations for both Geometries.
   * Initialize the locations for both Geometries to the given values.
   */
  Label.args3(int onLoc, int leftLoc, int rightLoc) {
    elt[0] = new TopologyLocation(onLoc, leftLoc, rightLoc);
    elt[1] = new TopologyLocation(onLoc, leftLoc, rightLoc);
  }

  /**
   * Construct a Label with On, Left and Right locations for both Geometries.
   * Initialize the locations for the given Geometry index.
   */
  Label.args4(int geomIndex, int onLoc, int leftLoc, int rightLoc) {
    elt[0] = new TopologyLocation(Location.NONE, Location.NONE, Location.NONE);
    elt[1] = new TopologyLocation(Location.NONE, Location.NONE, Location.NONE);
    elt[geomIndex]!.setLocations(onLoc, leftLoc, rightLoc);
  }

  /**
   * Construct a Label with the same values as the argument Label.
   */
  Label.fromLabel(Label lbl) {
    elt[0] = new TopologyLocation.fromTL(lbl.elt[0]!);
    elt[1] = new TopologyLocation.fromTL(lbl.elt[1]!);
  }

  void flip() {
    elt[0]!.flip();
    elt[1]!.flip();
  }

  int getLocationWithPosIndex(int geomIndex, int posIndex) {
    return elt[geomIndex]!.get(posIndex);
  }

  int getLocation(int geomIndex) {
    return elt[geomIndex]!.get(Position.ON);
  }

  void setLocation(int geomIndex, int posIndex, int location) {
    elt[geomIndex]!.setLocationWithIndex(posIndex, location);
  }

  void setLocationWithIndex(int geomIndex, int location) {
    elt[geomIndex]!.setLocationWithIndex(Position.ON, location);
  }

  void setAllLocations(int geomIndex, int location) {
    elt[geomIndex]!.setAllLocations(location);
  }

  void setAllLocationsIfNullWithIndex(int geomIndex, int location) {
    elt[geomIndex]!.setAllLocationsIfNull(location);
  }

  void setAllLocationsIfNull(int location) {
    setAllLocationsIfNullWithIndex(0, location);
    setAllLocationsIfNullWithIndex(1, location);
  }

  /**
   * Merge this label with another one.
   * Merging updates any null attributes of this label with the attributes from lbl
   */
  void merge(Label lbl) {
    for (int i = 0; i < 2; i++) {
      if (elt[i] == null && lbl.elt[i] != null) {
        elt[i] = new TopologyLocation.fromTL(lbl.elt[i]!);
      } else {
        elt[i]!.merge(lbl.elt[i]!);
      }
    }
  }

  int getGeometryCount() {
    int count = 0;
    if (!elt[0]!.isNull()) count++;
    if (!elt[1]!.isNull()) count++;
    return count;
  }

  bool isNull(int geomIndex) {
    return elt[geomIndex]!.isNull();
  }

  bool isAnyNull(int geomIndex) {
    return elt[geomIndex]!.isAnyNull();
  }

  bool isArea() {
    return elt[0]!.isArea() || elt[1]!.isArea();
  }

  bool isAreaWithIndex(int geomIndex) {
    /*  Testing
  	if (elt[0].getLocations().length != elt[1].getLocations().length) {
  		System.out.println(this);
  	}
  		*/
    return elt[geomIndex]!.isArea();
  }

  bool isLine(int geomIndex) {
    return elt[geomIndex]!.isLine();
  }

  bool isEqualOnSide(Label lbl, int side) {
    return this.elt[0]!.isEqualOnSide(lbl.elt[0]!, side) &&
        this.elt[1]!.isEqualOnSide(lbl.elt[1]!, side);
  }

  bool allPositionsEqual(int geomIndex, int loc) {
    return elt[geomIndex]!.allPositionsEqual(loc);
  }

  /**
   * Converts one GeometryLocation to a Line location
   */
  void toLine(int geomIndex) {
    if (elt[geomIndex]!.isArea())
      elt[geomIndex] = new TopologyLocation.fromOn(elt[geomIndex]!.location[0]);
  }

  String toString() {
    StringBuffer buf = new StringBuffer();
    if (elt[0] != null) {
      buf.write("A:");
      buf.write(elt[0].toString());
    }
    if (elt[1] != null) {
      buf.write(" B:");
      buf.write(elt[1].toString());
    }
    return buf.toString();
  }
}

/**
 * A GraphComponent is the parent class for the objects'
 * that form a graph.  Each GraphComponent can carry a
 * Label.
 * @version 1.7
 */
abstract class GraphComponent {
  Label? label;

  /**
   * isInResult indicates if this component has already been included in the result
   */
  bool _isInResult = false;
  bool _isCovered = false;
  bool _isCoveredSet = false;
  bool _isVisited = false;

  GraphComponent() {}

  GraphComponent.fromLabel(this.label);

  Label? getLabel() {
    return label;
  }

  void setLabel(Label label) {
    this.label = label;
  }

  void setInResult(bool isInResult) {
    this._isInResult = isInResult;
  }

  bool isInResult() {
    return _isInResult;
  }

  void setCovered(bool isCovered) {
    this._isCovered = isCovered;
    this._isCoveredSet = true;
  }

  bool isCovered() {
    return _isCovered;
  }

  bool isCoveredSet() {
    return _isCoveredSet;
  }

  bool isVisited() {
    return _isVisited;
  }

  void setVisited(bool isVisited) {
    this._isVisited = isVisited;
  }

  /**
   * @return a coordinate in this component (or null, if there are none)
   */
  Coordinate? getCoordinate();

  /**
   * compute the contribution to an IM for this component
   */
  void computeIM(IntersectionMatrix im);

  /**
   * An isolated component is one that does not intersect or touch any other
   * component.  This is the case if the label has valid locations for
   * only a single Geometry.
   *
   * @return true if this component is isolated
   */
  bool isIsolated();

  /**
   * Update the IM with the contribution for this component.
   * A component only contributes if it has a labelling for both parent geometries
   */
  void updateIM(IntersectionMatrix im) {
//  TODO  Assert.isTrue(label.getGeometryCount() >= 2, "found partial label");
    computeIM(im);
  }
}

/**
 * A EdgeList is a list of Edges.  It supports locating edges
 * that are pointwise equals to a target edge.
 * @version 1.7
 */
class EdgeList {
  List edges = [];

  /**
   * An index of the edges, for fast lookup.
   *
   */
  Map ocaMap = new SplayTreeMap();

  EdgeList() {}

  /**
   * Insert an edge unless it is already in the list
   */
  void add(Edge e) {
    edges.add(e);
    OrientedCoordinateArray oca =
        new OrientedCoordinateArray(e.getCoordinates());
    ocaMap[oca] = e;
  }

  void addAll(List<Edge> edgeColl) {
    for (var e in edgeColl) {
      add(e);
    }
  }

  List getEdges() {
    return edges;
  }

  /**
   * If there is an edge equal to e already in the list, return it.
   * Otherwise return null.
   * @return  equal edge, if there is one already in the list
   *          null otherwise
   */
  Edge? findEqualEdge(Edge e) {
    OrientedCoordinateArray oca =
        new OrientedCoordinateArray(e.getCoordinates());
    // will return null if no edge matches
    Edge? matchEdge = ocaMap[oca];
    return matchEdge;
  }

  Iterator iterator() {
    return edges.iterator;
  }

  Edge get(int i) {
    return edges[i];
  }

  /**
   * If the edge e is already in the list, return its index.
   * @return  index, if e is already in the list
   *          -1 otherwise
   */
  int findEdgeIndex(Edge e) {
    return edges.indexOf(e);
  }

//   void print(PrintStream out)
//  {
//    out.print("MULTILINESTRING ( ");
//    for (int j = 0; j < edges.size(); j++) {
//      Edge e = (Edge) edges.get(j);
//      if (j > 0) out.print(",");
//      out.print("(");
//      Coordinate[] pts = e.getCoordinates();
//      for (int i = 0; i < pts.length; i++) {
//        if (i > 0) out.print(",");
//        out.print(pts[i].x + " " + pts[i].y);
//      }
//      out.println(")");
//    }
//    out.print(")  ");
//  }

}
