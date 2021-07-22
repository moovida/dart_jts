part of dart_jts;

/**
 * Implements the SFS <tt>relate()</tt> generalized spatial predicate on two {@link Geometry}s.
 * <p>
 * The class supports specifying a custom {@link BoundaryNodeRule}
 * to be used during the relate computation.
 * <p>
 * If named spatial predicates are used on the result {@link IntersectionMatrix}
 * of the RelateOp, the result may or not be affected by the
 * choice of <tt>BoundaryNodeRule</tt>, depending on the exact nature of the pattern.
 * For instance, {@link IntersectionMatrix#isIntersects()} is insensitive
 * to the choice of <tt>BoundaryNodeRule</tt>,
 * whereas {@link IntersectionMatrix#isTouches(int, int)} is affected by the rule chosen.
 * <p>
 * <b>Note:</b> custom Boundary Node Rules do not (currently)
 * affect the results of other {@link Geometry} methods (such
 * as {@link Geometry#getBoundary}.  The results of
 * these methods may not be consistent with the relationship computed by
 * a custom Boundary Node Rule.
 *
 * @version 1.7
 */
class RelateOp extends GeometryGraphOperation {
  /**
   * Computes the {@link IntersectionMatrix} for the spatial relationship
   * between two {@link Geometry}s, using the default (OGC SFS) Boundary Node Rule
   *
   * @param a a Geometry to test
   * @param b a Geometry to test
   * @return the IntersectionMatrix for the spatial relationship between the geometries
   */
  static IntersectionMatrix relateStatic(Geometry a, Geometry b) {
    RelateOp relOp = new RelateOp(a, b);
    IntersectionMatrix im = relOp.getIntersectionMatrix();
    return im;
  }

  /**
   * Computes the {@link IntersectionMatrix} for the spatial relationship
   * between two {@link Geometry}s using a specified Boundary Node Rule.
   *
   * @param a a Geometry to test
   * @param b a Geometry to test
   * @param boundaryNodeRule the Boundary Node Rule to use
   * @return the IntersectionMatrix for the spatial relationship between the input geometries
   */
  static IntersectionMatrix relateStaticWithRule(
      Geometry a, Geometry b, BoundaryNodeRule boundaryNodeRule) {
    RelateOp relOp = new RelateOp.withRule(a, b, boundaryNodeRule);
    IntersectionMatrix im = relOp.getIntersectionMatrix();
    return im;
  }

  late RelateComputer relate;

  /**
   * Creates a new Relate operation, using the default (OGC SFS) Boundary Node Rule.
   *
   * @param g0 a Geometry to relate
   * @param g1 another Geometry to relate
   */
  RelateOp(Geometry g0, Geometry g1) : super(g0, g1) {
    relate = new RelateComputer(arg);
  }

  /**
   * Creates a new Relate operation with a specified Boundary Node Rule.
   *
   * @param g0 a Geometry to relate
   * @param g1 another Geometry to relate
   * @param boundaryNodeRule the Boundary Node Rule to use
   */
  RelateOp.withRule(Geometry g0, Geometry g1, BoundaryNodeRule boundaryNodeRule)
      : super.withRule(g0, g1, boundaryNodeRule) {
    relate = new RelateComputer(arg);
  }

  /**
   * Gets the IntersectionMatrix for the spatial relationship
   * between the input geometries.
   *
   * @return the IntersectionMatrix for the spatial relationship between the input geometries
   */
  IntersectionMatrix getIntersectionMatrix() {
    return relate.computeIM();
  }
}

/**
 * Computes the {@link EdgeEnd}s which arise from a noded {@link Edge}.
 *
 * @version 1.7
 */
class EdgeEndBuilder {
  EdgeEndBuilder() {}

  List computeEdgeEnds(Iterator edges) {
    List l = [];
    for (Iterator i = edges; i.moveNext();) {
      Edge e = i.current as Edge;
      computeEdgeEndsWithEdge(e, l);
    }
    return l;
  }

  /**
   * Creates stub edges for all the intersections in this
   * Edge (if any) and inserts them into the graph.
   */
  void computeEdgeEndsWithEdge(Edge edge, List l) {
    EdgeIntersectionList eiList = edge.getEdgeIntersectionList();
//Debug.print(eiList);
    // ensure that the list has entries for the first and last point of the edge
    eiList.addEndpoints();

    Iterator it = eiList.iterator();
    EdgeIntersection? eiPrev = null;
    EdgeIntersection? eiCurr = null;
    // no intersections, so there is nothing to do
    if (!it.moveNext()) return;
    EdgeIntersection? eiNext = it.current as EdgeIntersection;
    do {
      eiPrev = eiCurr;
      eiCurr = eiNext;
      eiNext = null;
      if (it.moveNext()) eiNext = it.current as EdgeIntersection;

      if (eiCurr != null) {
        createEdgeEndForPrev(edge, l, eiCurr, eiPrev);
        createEdgeEndForNext(edge, l, eiCurr, eiNext);
      }
    } while (eiCurr != null);
  }

  /**
   * Create a EdgeStub for the edge before the intersection eiCurr.
   * The previous intersection is provided
   * in case it is the endpoint for the stub edge.
   * Otherwise, the previous point from the parent edge will be the endpoint.
   * <br>
   * eiCurr will always be an EdgeIntersection, but eiPrev may be null.
   */
  void createEdgeEndForPrev(
      Edge edge, List l, EdgeIntersection eiCurr, EdgeIntersection? eiPrev) {
    int iPrev = eiCurr.segmentIndex;
    if (eiCurr.dist == 0.0) {
      // if at the start of the edge there is no previous edge
      if (iPrev == 0) return;
      iPrev--;
    }
    Coordinate pPrev = edge.getCoordinateWithIndex(iPrev);
    // if prev intersection is past the previous vertex, use it instead
    if (eiPrev != null && eiPrev.segmentIndex >= iPrev) pPrev = eiPrev.coord;

    Label label = new Label.fromLabel(edge.getLabel()!);
    // since edgeStub is oriented opposite to it's parent edge, have to flip sides for edge label
    label.flip();
    EdgeEnd e = new EdgeEnd.withCoordsLabel(edge, eiCurr.coord, pPrev, label);
//e.print(System.out);  System.out.println();
    l.add(e);
  }

  /**
   * Create a StubEdge for the edge after the intersection eiCurr.
   * The next intersection is provided
   * in case it is the endpoint for the stub edge.
   * Otherwise, the next point from the parent edge will be the endpoint.
   * <br>
   * eiCurr will always be an EdgeIntersection, but eiNext may be null.
   */
  void createEdgeEndForNext(
      Edge edge, List l, EdgeIntersection eiCurr, EdgeIntersection? eiNext) {
    int iNext = eiCurr.segmentIndex + 1;
    // if there is no next edge there is nothing to do
    if (iNext >= edge.getNumPoints() && eiNext == null) return;

    Coordinate pNext = edge.getCoordinateWithIndex(iNext);

    // if the next intersection is in the same segment as the current, use it as the endpoint
    if (eiNext != null && eiNext.segmentIndex == eiCurr.segmentIndex)
      pNext = eiNext.coord;

    EdgeEnd e = new EdgeEnd.withCoordsLabel(
        edge, eiCurr.coord, pNext, new Label.fromLabel(edge.getLabel()!));
//Debug.println(e);
    l.add(e);
  }
}

/**
 * A collection of {@link EdgeEnd}s which obey the following invariant:
 * They originate at the same node and have the same direction.
 *
 * @version 1.7
 */
class EdgeEndBundle extends EdgeEnd {
//   BoundaryNodeRule boundaryNodeRule;
  List edgeEnds = [];

  EdgeEndBundle.withRule(BoundaryNodeRule? boundaryNodeRule, EdgeEnd e)
      : super.withCoordsLabel(e.getEdge(), e.getCoordinate()!,
            e.getDirectedCoordinate()!, new Label.fromLabel(e.getLabel()!)) {
    insert(e);
    /*
    if (boundaryNodeRule != null)
      this.boundaryNodeRule = boundaryNodeRule;
    else
      boundaryNodeRule = BoundaryNodeRule.OGC_SFS_BOUNDARY_RULE;
    */
  }

  EdgeEndBundle.withEdgeEnd(EdgeEnd e) : this.withRule(null, e);

  Label? getLabel() {
    return label;
  }

  Iterator iterator() {
    return edgeEnds.iterator;
  }

  List getEdgeEnds() {
    return edgeEnds;
  }

  void insert(EdgeEnd e) {
    // Assert: start point is the same
    // Assert: direction is the same
    edgeEnds.add(e);
  }

  /**
   * This computes the overall edge label for the set of
   * edges in this EdgeStubBundle.  It essentially merges
   * the ON and side labels for each edge.  These labels must be compatible
   */
  void computeLabel(BoundaryNodeRule boundaryNodeRule) {
    // create the label.  If any of the edges belong to areas,
    // the label must be an area label
    bool isArea = false;
    for (Iterator it = iterator(); it.moveNext();) {
      EdgeEnd e = it.current as EdgeEnd;
      if (e.getLabel()!.isArea()) isArea = true;
    }
    if (isArea)
      label = new Label.args3(Location.NONE, Location.NONE, Location.NONE);
    else
      label = new Label(Location.NONE);

    // compute the On label, and the side labels if present
    for (int i = 0; i < 2; i++) {
      computeLabelOn(i, boundaryNodeRule);
      if (isArea) computeLabelSides(i);
    }
  }

  /**
   * Compute the overall ON location for the list of EdgeStubs.
   * (This is essentially equivalent to computing the self-overlay of a single Geometry)
   * edgeStubs can be either on the boundary (e.g. Polygon edge)
   * OR in the interior (e.g. segment of a LineString)
   * of their parent Geometry.
   * In addition, GeometryCollections use a {@link BoundaryNodeRule} to determine
   * whether a segment is on the boundary or not.
   * Finally, in GeometryCollections it can occur that an edge is both
   * on the boundary and in the interior (e.g. a LineString segment lying on
   * top of a Polygon edge.) In this case the Boundary is given precedence.
   * <br>
   * These observations result in the following rules for computing the ON location:
   * <ul>
   * <li> if there are an odd number of Bdy edges, the attribute is Bdy
   * <li> if there are an even number >= 2 of Bdy edges, the attribute is Int
   * <li> if there are any Int edges, the attribute is Int
   * <li> otherwise, the attribute is NULL.
   * </ul>
   */
  void computeLabelOn(int geomIndex, BoundaryNodeRule boundaryNodeRule) {
    // compute the ON location value
    int boundaryCount = 0;
    bool foundInterior = false;

    for (Iterator it = iterator(); it.moveNext();) {
      EdgeEnd e = it.current as EdgeEnd;
      int loc = e.getLabel()!.getLocation(geomIndex);
      if (loc == Location.BOUNDARY) boundaryCount++;
      if (loc == Location.INTERIOR) foundInterior = true;
    }
    int loc = Location.NONE;
    if (foundInterior) loc = Location.INTERIOR;
    if (boundaryCount > 0) {
      loc = GeometryGraph.determineBoundary(boundaryNodeRule, boundaryCount);
    }
    label!.setLocationWithIndex(geomIndex, loc);
  }

  /**
   * Compute the labelling for each side
   */
  void computeLabelSides(int geomIndex) {
    computeLabelSide(geomIndex, Position.LEFT);
    computeLabelSide(geomIndex, Position.RIGHT);
  }

  /**
   * To compute the summary label for a side, the algorithm is:
   *   FOR all edges
   *     IF any edge's location is INTERIOR for the side, side location = INTERIOR
   *     ELSE IF there is at least one EXTERIOR attribute, side location = EXTERIOR
   *     ELSE  side location = NULL
   *  <br>
   *  Note that it is possible for two sides to have apparently contradictory information
   *  i.e. one edge side may indicate that it is in the interior of a geometry, while
   *  another edge side may indicate the exterior of the same geometry.  This is
   *  not an incompatibility - GeometryCollections may contain two Polygons that touch
   *  along an edge.  This is the reason for Interior-primacy rule above - it
   *  results in the summary label having the Geometry interior on <b>both</b> sides.
   */
  void computeLabelSide(int geomIndex, int side) {
    for (Iterator it = iterator(); it.moveNext();) {
      EdgeEnd e = it.current as EdgeEnd;
      if (e.getLabel()!.isArea()) {
        int loc = e.getLabel()!.getLocationWithPosIndex(geomIndex, side);
        if (loc == Location.INTERIOR) {
          label!.setLocation(geomIndex, side, Location.INTERIOR);
          return;
        } else if (loc == Location.EXTERIOR)
          label!.setLocation(geomIndex, side, Location.EXTERIOR);
      }
    }
  }

  /**
   * Update the IM with the contribution for the computed label for the EdgeStubs.
   */
  void updateIM(IntersectionMatrix im) {
    Edge.updateIMStatic(label!, im);
  }
//   void print(PrintStream out)
//  {
//    out.println("EdgeEndBundle--> Label: " + label);
//    for (Iterator it = iterator(); it.hasNext(); ) {
//      EdgeEnd ee = (EdgeEnd) it.next();
//      ee.print(out);
//      out.println();
//    }
//  }
}

/**
 * An ordered list of {@link EdgeEndBundle}s around a {@link RelateNode}.
 * They are maintained in CCW order (starting with the positive x-axis) around the node
 * for efficient lookup and topology building.
 *
 * @version 1.7
 */
class EdgeEndBundleStar extends EdgeEndStar {
  /**
   * Creates a new empty EdgeEndBundleStar
   */
  EdgeEndBundleStar() {}

  /**
   * Insert a EdgeEnd in order in the list.
   * If there is an existing EdgeStubBundle which is parallel, the EdgeEnd is
   * added to the bundle.  Otherwise, a new EdgeEndBundle is created
   * to contain the EdgeEnd.
   * <br>
   */
  void insert(EdgeEnd e) {
    EdgeEndBundle? eb = edgeMap[e];
    if (eb == null) {
      eb = new EdgeEndBundle.withEdgeEnd(e);
      insertEdgeEnd(e, eb);
    } else {
      eb.insert(e);
    }
  }

  /**
   * Update the IM with the contribution for the EdgeStubs around the node.
   */
  void updateIM(IntersectionMatrix im) {
    for (Iterator it = iterator(); it.moveNext();) {
      EdgeEndBundle esb = it.current as EdgeEndBundle;
      esb.updateIM(im);
    }
  }
}

/**
 * Represents a node in the topological graph used to compute spatial relationships.
 *
 * @version 1.7
 */
class RelateNode extends Node {
  RelateNode(Coordinate coord, EdgeEndStar edges) : super(coord, edges);

  /**
   * Update the IM with the contribution for this component.
   * A component only contributes if it has a labelling for both parent geometries
   */
  void computeIM(IntersectionMatrix im) {
    im.setAtLeastIfValid(label!.getLocation(0), label!.getLocation(1), 0);
  }

  /**
   * Update the IM with the contribution for the EdgeEnds incident on this node.
   */
  void updateIMFromEdges(IntersectionMatrix im) {
    (edges as EdgeEndBundleStar).updateIM(im);
  }
}

/**
 * Used by the {@link NodeMap} in a {@link RelateNodeGraph} to create {@link RelateNode}s.
 *
 * @version 1.7
 */
class RelateNodeFactory extends NodeFactory {
  Node createNode(Coordinate coord) {
    return new RelateNode(coord, new EdgeEndBundleStar());
  }
}

/**
 * Implements the simple graph of Nodes and EdgeEnd which is all that is
 * required to determine topological relationships between Geometries.
 * Also supports building a topological graph of a single Geometry, to
 * allow verification of valid topology.
 * <p>
 * It is <b>not</b> necessary to create a fully linked
 * PlanarGraph to determine relationships, since it is sufficient
 * to know how the Geometries interact locally around the nodes.
 * In fact, this is not even feasible, since it is not possible to compute
 * exact intersection points, and hence the topology around those nodes
 * cannot be computed robustly.
 * The only Nodes that are created are for improper intersections;
 * that is, nodes which occur at existing vertices of the Geometries.
 * Proper intersections (e.g. ones which occur between the interior of line segments)
 * have their topology determined implicitly, without creating a Node object
 * to represent them.
 *
 * @version 1.7
 */
class RelateNodeGraph {
  NodeMap nodes = new NodeMap(new RelateNodeFactory());

  RelateNodeGraph() {}

  Iterator getNodeIterator() {
    return nodes.iterator();
  }

  void build(GeometryGraph geomGraph) {
    // compute nodes for intersections between previously noded edges
    computeIntersectionNodes(geomGraph, 0);
    /**
     * Copy the labelling for the nodes in the parent Geometry.  These override
     * any labels determined by intersections.
     */
    copyNodesAndLabels(geomGraph, 0);

    /**
     * Build EdgeEnds for all intersections.
     */
    EdgeEndBuilder eeBuilder = new EdgeEndBuilder();
    List eeList = eeBuilder.computeEdgeEnds(geomGraph.getEdgeIterator());
    insertEdgeEnds(eeList);

//Debug.println("==== NodeList ===");
//Debug.print(nodes);
  }

  /**
   * Insert nodes for all intersections on the edges of a Geometry.
   * Label the created nodes the same as the edge label if they do not already have a label.
   * This allows nodes created by either self-intersections or
   * mutual intersections to be labelled.
   * Endpoint nodes will already be labelled from when they were inserted.
   * <p>
   * Precondition: edge intersections have been computed.
   */
  void computeIntersectionNodes(GeometryGraph geomGraph, int argIndex) {
    for (Iterator edgeIt = geomGraph.getEdgeIterator(); edgeIt.moveNext();) {
      Edge e = edgeIt.current as Edge;
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
   * Copy all nodes from an arg geometry into this graph.
   * The node label in the arg geometry overrides any previously computed
   * label for that argIndex.
   * (E.g. a node may be an intersection node with
   * a computed label of BOUNDARY,
   * but in the original arg Geometry it is actually
   * in the interior due to the Boundary Determination Rule)
   */
  void copyNodesAndLabels(GeometryGraph geomGraph, int argIndex) {
    for (Iterator nodeIt = geomGraph.getNodeIterator(); nodeIt.moveNext();) {
      Node graphNode = nodeIt.current as Node;
      Node newNode = nodes.addNodeFromCoordinate(graphNode.getCoordinate());
      newNode.setLabelWithIndex(
          argIndex, graphNode.getLabel()!.getLocation(argIndex));
//node.print(System.out);
    }
  }

  void insertEdgeEnds(List ee) {
    for (Iterator i = ee.iterator; i.moveNext();) {
      EdgeEnd e = i.current as EdgeEnd;
      nodes.add(e);
    }
  }
}
