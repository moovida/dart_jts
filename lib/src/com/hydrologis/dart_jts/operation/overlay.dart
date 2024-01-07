part of dart_jts;

/**
 * A ring of {@link DirectedEdge}s which may contain nodes of degree &gt; 2.
 * A <tt>MaximalEdgeRing</tt> may represent two different spatial entities:
 * <ul>
 * <li>a single polygon possibly containing inversions (if the ring is oriented CW)
 * <li>a single hole possibly containing exversions (if the ring is oriented CCW)
 * </ul>
 * If the MaximalEdgeRing represents a polygon,
 * the interior of the polygon is strongly connected.
 * <p>
 * These are the form of rings used to define polygons under some spatial data models.
 * However, under the OGC SFS model, {@link MinimalEdgeRing}s are required.
 * A MaximalEdgeRing can be converted to a list of MinimalEdgeRings using the
 * {@link #buildMinimalRings() } method.
 *
 * @version 1.7
 * @see org.locationtech.jts.operation.overlay.MinimalEdgeRing
 */
class MaximalEdgeRing extends EdgeRing {
  MaximalEdgeRing(DirectedEdge start, GeometryFactory geometryFactory)
      : super(start, geometryFactory);

  DirectedEdge getNext(DirectedEdge de) {
    return de.getNext();
  }

  void setEdgeRing(DirectedEdge de, EdgeRing er) {
    de.setEdgeRing(er);
  }

  /**
   * For all nodes in this EdgeRing,
   * link the DirectedEdges at the node to form minimalEdgeRings
   */
  void linkDirectedEdgesForMinimalEdgeRings() {
    DirectedEdge de = startDe!;
    do {
      Node node = de.getNode()!;
      (node.getEdges() as DirectedEdgeStar).linkMinimalDirectedEdges(this);
      de = de.getNext();
    } while (de != startDe);
  }

  List buildMinimalRings() {
    List minEdgeRings = [];
    DirectedEdge de = startDe!;
    do {
      if (de.getMinEdgeRing() == null) {
        EdgeRing minEr = new MinimalEdgeRing(de, geometryFactory);
        minEdgeRings.add(minEr);
      }
      de = de.getNext();
    } while (de != startDe);
    return minEdgeRings;
  }
}

/**
 * A ring of {@link Edge}s with the property that no node
 * has degree greater than 2.  These are the form of rings required
 * to represent polygons under the OGC SFS spatial data model.
 *
 * @version 1.7
 * @see org.locationtech.jts.operation.overlay.MaximalEdgeRing
 */
class MinimalEdgeRing extends EdgeRing {
  MinimalEdgeRing(DirectedEdge start, GeometryFactory geometryFactory)
      : super(start, geometryFactory);

  DirectedEdge getNext(DirectedEdge de) {
    return de.getNextMin();
  }

  void setEdgeRing(DirectedEdge de, EdgeRing er) {
    de.setMinEdgeRing(er);
  }
}

/**
 * Forms {@link Polygon}s out of a graph of {@link DirectedEdge}s.
 * The edges to use are marked as being in the result Area.
 * <p>
 *
 * @version 1.7
 */
class PolygonBuilder {
  GeometryFactory geometryFactory;
  List shellList = [];

  PolygonBuilder(this.geometryFactory);

  /**
   * Add a complete graph.
   * The graph is assumed to contain one or more polygons,
   * possibly with holes.
   */
  void addGraph(PlanarGraph graph) {
    add(graph.getEdgeEnds(), graph.getNodes());
  }

  /**
   * Add a set of edges and nodes, which form a graph.
   * The graph is assumed to contain one or more polygons,
   * possibly with holes.
   */
  void add(List dirEdges, List nodes) {
    PlanarGraph.linkResultDirectedEdgesStatic(nodes);
    List maxEdgeRings = buildMaximalEdgeRings(dirEdges);
    List freeHoleList = [];
    List edgeRings =
        buildMinimalEdgeRings(maxEdgeRings, shellList, freeHoleList);
    sortShellsAndHoles(edgeRings, shellList, freeHoleList);
    placeFreeHoles(shellList, freeHoleList);
    //Assert: every hole on freeHoleList has a shell assigned to it
  }

  List<Polygon> getPolygons() {
    List<Polygon> resultPolyList = computePolygons(shellList);
    return resultPolyList;
  }

  /**
   * for all DirectedEdges in result, form them into MaximalEdgeRings
   */
  List buildMaximalEdgeRings(List dirEdges) {
    List maxEdgeRings = [];
    for (DirectedEdge de in dirEdges) {
      if (de.isInResult() && de.getLabel()!.isArea()) {
        // if this edge has not yet been processed
        if (de.getEdgeRing() == null) {
          MaximalEdgeRing er = new MaximalEdgeRing(de, geometryFactory);
          maxEdgeRings.add(er);
          er.setInResult();
//System.out.println("max node degree = " + er.getMaxDegree());
        }
      }
    }
    return maxEdgeRings;
  }

  List buildMinimalEdgeRings(
      List maxEdgeRings, List shellList, List freeHoleList) {
    List edgeRings = [];
    for (MaximalEdgeRing er in maxEdgeRings) {
      if (er.getMaxNodeDegree() > 2) {
        er.linkDirectedEdgesForMinimalEdgeRings();
        List minEdgeRings = er.buildMinimalRings();
        // at this point we can go ahead and attempt to place holes, if this EdgeRing is a polygon
        EdgeRing? shell = findShell(minEdgeRings);
        if (shell != null) {
          placePolygonHoles(shell, minEdgeRings);
          shellList.add(shell);
        } else {
          freeHoleList.addAll(minEdgeRings);
        }
      } else {
        edgeRings.add(er);
      }
    }
    return edgeRings;
  }

  /**
   * This method takes a list of MinimalEdgeRings derived from a MaximalEdgeRing,
   * and tests whether they form a Polygon.  This is the case if there is a single shell
   * in the list.  In this case the shell is returned.
   * The other possibility is that they are a series of connected holes, in which case
   * no shell is returned.
   *
   * @return the shell EdgeRing, if there is one
   * or null, if all the rings are holes
   */
  EdgeRing? findShell(List minEdgeRings) {
    int shellCount = 0;
    EdgeRing? shell;
    for (EdgeRing er in minEdgeRings) {
      if (!er.isHole()) {
        shell = er;
        shellCount++;
      }
    }
    Assert.isTrue(shellCount <= 1, "found two shells in MinimalEdgeRing list");
    return shell;
  }

  /**
   * This method assigns the holes for a Polygon (formed from a list of
   * MinimalEdgeRings) to its shell.
   * Determining the holes for a MinimalEdgeRing polygon serves two purposes:
   * <ul>
   * <li>it is faster than using a point-in-polygon check later on.
   * <li>it ensures correctness, since if the PIP test was used the point
   * chosen might lie on the shell, which might return an incorrect result from the
   * PIP test
   * </ul>
   */
  void placePolygonHoles(EdgeRing shell, List minEdgeRings) {
    for (MinimalEdgeRing er in minEdgeRings) {
      if (er.isHole()) {
        er.setShell(shell);
      }
    }
  }

  /**
   * For all rings in the input list,
   * determine whether the ring is a shell or a hole
   * and add it to the appropriate list.
   * Due to the way the DirectedEdges were linked,
   * a ring is a shell if it is oriented CW, a hole otherwise.
   */
  void sortShellsAndHoles(List edgeRings, List shellList, List freeHoleList) {
    for (EdgeRing er in edgeRings) {
//      er.setInResult();
      if (er.isHole()) {
        freeHoleList.add(er);
      } else {
        shellList.add(er);
      }
    }
  }

  /**
   * This method determines finds a containing shell for all holes
   * which have not yet been assigned to a shell.
   * These "free" holes should
   * all be <b>properly</b> contained in their parent shells, so it is safe to use the
   * <code>findEdgeRingContaining</code> method.
   * (This is the case because any holes which are NOT
   * properly contained (i.e. are connected to their
   * parent shell) would have formed part of a MaximalEdgeRing
   * and been handled in a previous step).
   *
   * @throws TopologyException if a hole cannot be assigned to a shell
   */
  void placeFreeHoles(List shellList, List freeHoleList) {
    for (EdgeRing hole in freeHoleList) {
      // only place this hole if it doesn't yet have a shell
      if (hole.getShell() == null) {
        EdgeRing? shell = findEdgeRingContaining(hole, shellList);
        if (shell == null)
          throw new TopologyException.withCoord(
              "unable to assign hole to a shell", hole.getCoordinate(0));
//        Assert.isTrue(shell != null, "unable to assign hole to a shell");
        hole.setShell(shell);
      }
    }
  }

  /**
   * Find the innermost enclosing shell EdgeRing containing the argument EdgeRing, if any.
   * The innermost enclosing ring is the <i>smallest</i> enclosing ring.
   * The algorithm used depends on the fact that:
   * <br>
   *  ring A contains ring B iff envelope(ring A) contains envelope(ring B)
   * <br>
   * This routine is only safe to use if the chosen point of the hole
   * is known to be properly contained in a shell
   * (which is guaranteed to be the case if the hole does not touch its shell)
   *
   * @return containing EdgeRing, if there is one
   * or null if no containing EdgeRing is found
   */
  static EdgeRing? findEdgeRingContaining(EdgeRing testEr, List shellList) {
    LinearRing testRing = testEr.getLinearRing()!;
    Envelope testEnv = testRing.getEnvelopeInternal();
    Coordinate? testPt = testRing.getCoordinateN(0);

    EdgeRing? minShell = null;
    Envelope? minShellEnv = null;
    for (EdgeRing tryShell in shellList) {
      LinearRing tryShellRing = tryShell.getLinearRing()!;
      Envelope tryShellEnv = tryShellRing.getEnvelopeInternal();
      // the hole envelope cannot equal the shell envelope
      // (also guards against testing rings against themselves)
      if (tryShellEnv == testEnv) continue;
      // hole must be contained in shell
      if (!tryShellEnv.containsEnvelope(testEnv)) continue;

      testPt = CoordinateArrays.ptNotInList(
          testRing.getCoordinates(), tryShellRing.getCoordinates());
      bool isContained = false;
      if (PointLocation.isInRing(testPt!, tryShellRing.getCoordinates()))
        isContained = true;

      // check if this new containing ring is smaller than the current minimum ring
      if (isContained) {
        if (minShell == null || minShellEnv!.containsEnvelope(tryShellEnv)) {
          minShell = tryShell;
          minShellEnv = minShell.getLinearRing()!.getEnvelopeInternal();
        }
      }
    }
    return minShell;
  }

  List<Polygon> computePolygons(List shellList) {
    List<Polygon> resultPolyList = [];
    // add Polygons for all shells
    for (EdgeRing er in shellList) {
      Polygon poly = er.toPolygon(geometryFactory);
      resultPolyList.add(poly);
    }
    return resultPolyList;
  }
}

/**
 * Functions to reduce the precision of a geometry
 * by rounding it to a given precision model.
 * <p>
 * This class handles only polygonal and linear inputs.
 * For full functionality see {@link org.locationtech.jts.precision.GeometryPrecisionReducer}.
 *
 * @see org.locationtech.jts.precision.GeometryPrecisionReducer
 * @author Martin Davis
 */
class PrecisionReducer {
  /**
   * Reduces the precision of a geometry by rounding and snapping it to the
   * supplied {@link PrecisionModel}.
   * The input geometry must be polygonal or linear.
   * <p>
   * The output is always a valid geometry.  This implies that input components
   * may be merged if they are closer than the grid precision.
   * if merging is not desired, then the individual geometry components
   * should be processed separately.
   * <p>
   * The output is fully noded
   * (i.e. coincident lines are merged and noded).
   * This provides an effective way to node / snap-round a collection of {@link LineString}s.
   *
   * @param geom the geometry to reduce
   * @param pm the precision model to use
   * @return the precision-reduced geometry
   *
   * @throws IllegalArgumentException if the reduction fails due to invalid input geometry is invalid
   */
  static Geometry? reducePrecision(Geometry geom, PrecisionModel pm) {
    OverlayNG ov = OverlayNG.fromSingle(geom, pm);
    /**
     * Ensure reducing a area only produces polygonal result.
     * (I.e. collapse lines are not output)
     */
    if (geom.getDimension() == 2) {
      ov.setAreaResultOnly(true);
    }
    try {
      Geometry? reduced = ov.getResult();
      return reduced;
    } on Error catch (ex,stacktrace){
      print(stacktrace);
        print(ex);
      throw Exception("Reduction failed, possible invalid input");
    }
  }

  PrecisionReducer() {
    // no instantiation for now
  }
}

/**
 * Computes the geometric overlay of two {@link Geometry}s.  The overlay
 * can be used to determine any boolean combination of the geometries.
 *
 * @version 1.7
 */
class OverlayOp extends GeometryGraphOperation {
  /**
   * The spatial functions supported by this class.
   * These operations implement various boolean combinations of the resultants of the overlay.
   */

  /**
   * The code for the Intersection overlay operation.
   */
  static const INTERSECTION = 1;

  /**
   * The code for the Union overlay operation.
   */
  static const UNION = 2;

  /**
   *  The code for the Difference overlay operation.
   */
  static const DIFFERENCE = 3;

  /**
   *  The code for the Symmetric Difference overlay operation.
   */
  static const SYMDIFFERENCE = 4;

  /**
   * Computes an overlay operation for
   * the given geometry arguments.
   *
   * @param geom0 the first geometry argument
   * @param geom1 the second geometry argument
   * @param opCode the code for the desired overlay operation
   * @return the result of the overlay operation
   * @throws TopologyException if a robustness problem is encountered
   */
  static Geometry? overlayOp(Geometry geom0, Geometry geom1, int opCode) {
    OverlayOp gov = new OverlayOp(geom0, geom1);
    Geometry? geomOv = gov.getResultGeometry(opCode);
    return geomOv;
  }

  /**
   * Tests whether a point with a given topological {@link Label}
   * relative to two geometries is contained in
   * the result of overlaying the geometries using
   * a given overlay operation.
   * <p>
   * The method handles arguments of {@link Location#NONE} correctly
   *
   * @param label the topological label of the point
   * @param opCode the code for the overlay operation to test
   * @return true if the label locations correspond to the overlayOpCode
   */
  static bool isResultOfOpLabel(Label label, int opCode) {
    int loc0 = label.getLocation(0);
    int loc1 = label.getLocation(1);
    return isResultOfOp(loc0, loc1, opCode);
  }

  /**
   * Tests whether a point with given {@link Location}s
   * relative to two geometries is contained in
   * the result of overlaying the geometries using
   * a given overlay operation.
   * <p>
   * The method handles arguments of {@link Location#NONE} correctly
   *
   * @param loc0 the code for the location in the first geometry
   * @param loc1 the code for the location in the second geometry
   * @param overlayOpCode the code for the overlay operation to test
   * @return true if the locations correspond to the overlayOpCode
   */
  static bool isResultOfOp(int loc0, int loc1, int overlayOpCode) {
    if (loc0 == Location.BOUNDARY) loc0 = Location.INTERIOR;
    if (loc1 == Location.BOUNDARY) loc1 = Location.INTERIOR;
    switch (overlayOpCode) {
      case INTERSECTION:
        return loc0 == Location.INTERIOR && loc1 == Location.INTERIOR;
      case UNION:
        return loc0 == Location.INTERIOR || loc1 == Location.INTERIOR;
      case DIFFERENCE:
        return loc0 == Location.INTERIOR && loc1 != Location.INTERIOR;
      case SYMDIFFERENCE:
        return (loc0 == Location.INTERIOR && loc1 != Location.INTERIOR) ||
            (loc0 != Location.INTERIOR && loc1 == Location.INTERIOR);
    }
    return false;
  }

  final PointLocator ptLocator = new PointLocator();
  late GeometryFactory geomFact;
  Geometry? resultGeom;

  late PlanarGraph graph;
  EdgeList edgeList = EdgeList();

  List<Polygon> resultPolyList = List.empty(growable: true);
  List<LineString> resultLineList = List.empty(growable: true);
  List<Point> resultPointList = List.empty(growable: true);

  /**
   * Constructs an instance to compute a single overlay operation
   * for the given geometries.
   *
   * @param g0 the first geometry argument
   * @param g1 the second geometry argument
   */
  OverlayOp(Geometry g0, Geometry g1) : super(g0, g1) {
    graph = PlanarGraph.withFactory(OverlayNodeFactory());
    /**
     * Use factory of primary geometry.
     * Note that this does NOT handle mixed-precision arguments
     * where the second arg has greater precision than the first.
     */
    geomFact = g0.getFactory();
  }

  /**
   * Gets the result of the overlay for a given overlay operation.
   * <p>
   * Note: this method can be called once only.
   *
   * @param overlayOpCode the overlay operation to perform
   * @return the compute result geometry
   * @throws TopologyException if a robustness problem is encountered
   */
  Geometry? getResultGeometry(int overlayOpCode) {
    computeOverlay(overlayOpCode);
    return resultGeom;
  }

  /**
   * Gets the graph constructed to compute the overlay.
   *
   * @return the overlay graph
   */
  PlanarGraph getGraph() {
    return graph;
  }

  void computeOverlay(int opCode) {
    // copy points from input Geometries.
    // This ensures that any Point geometries
    // in the input are considered for inclusion in the result set
    copyPoints(0);
    copyPoints(1);

    // node the input Geometries
    arg[0].computeSelfNodes(li, false);
    arg[1].computeSelfNodes(li, false);

    // compute intersections between edges of the two input geometries
    arg[0].computeEdgeIntersections(arg[1], li, true);

    List<Edge> baseSplitEdges = List.empty(growable: true);
    arg[0].computeSplitEdges(baseSplitEdges);
    arg[1].computeSplitEdges(baseSplitEdges);
    // add the noded edges to this result graph
    insertUniqueEdges(baseSplitEdges);

    computeLabelsFromDepths();
    replaceCollapsedEdges();

    /**
     * Check that the noding completed correctly.
     *
     * This test is slow, but necessary in order to catch robustness failure
     * situations.
     * If an exception is thrown because of a noding failure,
     * then snapping will be performed, which will hopefully avoid the problem.
     * In the future hopefully a faster check can be developed.
     *
     */
    EdgeNodingValidator.checkValidStatic(edgeList.getEdges());

    graph.addEdges(edgeList.getEdges());
    computeLabelling();
//Debug.printWatch();
    labelIncompleteNodes();
//Debug.printWatch();
//nodeMap.print(System.out);

    /**
     * The ordering of building the result Geometries is important.
     * Areas must be built before lines, which must be built before points.
     * This is so that lines which are covered by areas are not included
     * explicitly, and similarly for points.
     */
    findResultAreaEdges(opCode);
    cancelDuplicateResultEdges();

    PolygonBuilder polyBuilder = new PolygonBuilder(geomFact);
    polyBuilder.addGraph(graph);
    resultPolyList = polyBuilder.getPolygons();

    LineBuilder lineBuilder = LineBuilder(this, geomFact, ptLocator);
    resultLineList = lineBuilder.build(opCode);

    PointBuilder pointBuilder = new PointBuilder(this, geomFact);
    resultPointList = pointBuilder.build(opCode);

    // gather the results from all calculations into a single Geometry for the result set
    resultGeom = computeGeometry(
        resultPointList, resultLineList, resultPolyList, opCode);
  }

  void insertUniqueEdges(List edges) {
    for (var e in edges) {
      insertUniqueEdge(e);
    }
  }

  /**
   * Insert an edge from one of the noded input graphs.
   * Checks edges that are inserted to see if an
   * identical edge already exists.
   * If so, the edge is not inserted, but its label is merged
   * with the existing edge.
   */
  void insertUniqueEdge(Edge e) {
//<FIX> MD 8 Oct 03  speed up identical edge lookup
    // fast lookup
    Edge? existingEdge = edgeList.findEqualEdge(e);

    // If an identical edge already exists, simply update its label
    if (existingEdge != null) {
      Label? existingLabel = existingEdge.getLabel();

      Label? labelToMerge = e.getLabel();
      // check if new edge is in reverse direction to existing edge
      // if so, must flip the label before merging it
      if (!existingEdge.isPointwiseEqual(e)) {
        labelToMerge = Label.fromLabel(e.getLabel()!);
        labelToMerge.flip();
      }
      Depth depth = existingEdge.getDepth();
      // if this is the first duplicate found for this edge, initialize the depths
      ///*
      if (depth.isNull()) {
        depth.add(existingLabel!);
      }
      //*/
      depth.add(labelToMerge!);
      existingLabel!.merge(labelToMerge);
    } else {
      // no matching existing edge was found
      // add this new edge to the list of edges in this graph
      //e.setName(name + edges.size());
      //e.getDepth().add(e.getLabel());
      edgeList.add(e);
    }
  }

  /**
   * If either of the GeometryLocations for the existing label is
   * exactly opposite to the one in the labelToMerge,
   * this indicates a dimensional collapse has happened.
   * In this case, convert the label for that Geometry to a Line label
   */
  /**
   * Update the labels for edges according to their depths.
   * For each edge, the depths are first normalized.
   * Then, if the depths for the edge are equal,
   * this edge must have collapsed into a line edge.
   * If the depths are not equal, update the label
   * with the locations corresponding to the depths
   * (i.e. a depth of 0 corresponds to a Location of EXTERIOR,
   * a depth of 1 corresponds to INTERIOR)
   */
  void computeLabelsFromDepths() {
    for (var e in edgeList.edges) {
      Label? lbl = e.getLabel();
      Depth depth = e.getDepth();
      /**
       * Only check edges for which there were duplicates,
       * since these are the only ones which might
       * be the result of dimensional collapses.
       */
      if (!depth.isNull()) {
        depth.normalize();
        for (int i = 0; i < 2; i++) {
          if (!lbl!.isNull(i) && lbl.isArea() && !depth.isNull1(i)) {
            /**
             * if the depths are equal, this edge is the result of
             * the dimensional collapse of two or more edges.
             * It has the same location on both sides of the edge,
             * so it has collapsed to a line.
             */
            if (depth.getDelta(i) == 0) {
              lbl.toLine(i);
            } else {
              /**
               * This edge may be the result of a dimensional collapse,
               * but it still has different locations on both sides.  The
               * label of the edge must be updated to reflect the resultant
               * side locations indicated by the depth values.
               */
              Assert.isTrue(!depth.isNull2(i, Position.LEFT),
                  "depth of LEFT side has not been initialized");
              lbl.setLocation(
                  i, Position.LEFT, depth.getLocation(i, Position.LEFT));
              Assert.isTrue(!depth.isNull2(i, Position.RIGHT),
                  "depth of RIGHT side has not been initialized");
              lbl.setLocation(
                  i, Position.RIGHT, depth.getLocation(i, Position.RIGHT));
            }
          }
        }
      }
    }
  }

  /**
   * If edges which have undergone dimensional collapse are found,
   * replace them with a new edge which is a L edge
   */
  void replaceCollapsedEdges() {
    List<Edge> newEdges = List.empty(growable: true);
    var toRemove = [];
    for (Iterator it = edgeList.iterator(); it.moveNext();) {
      Edge e = it.current;
      if (e.isCollapsed()) {
//        it.remove();
        toRemove.add(e);
        newEdges.add(e.getCollapsedEdge());
      }
    }

    edgeList.addAll(newEdges);
  }

  /**
   * Copy all nodes from an arg geometry into this graph.
   * The node label in the arg geometry overrides any previously computed
   * label for that argIndex.
   * (E.g. a node may be an intersection node with
   * a previously computed label of BOUNDARY,
   * but in the original arg Geometry it is actually
   * in the interior due to the Boundary Determination Rule)
   */
  void copyPoints(int argIndex) {
    for (Iterator i = arg[argIndex].getNodeIterator(); i.moveNext();) {
      Node graphNode = i.current;
      Node newNode = graph.addNodeFromCoordinate(graphNode.getCoordinate());
      if (graphNode.getLabel() != null) {
        newNode.setLabelWithIndex(
            argIndex, graphNode.getLabel()!.getLocation(argIndex));
      }
    }
  }

  /**
   * Compute initial labelling for all DirectedEdges at each node.
   * In this step, DirectedEdges will acquire a complete labelling
   * (i.e. one with labels for both Geometries)
   * only if they
   * are incident on a node which has edges for both Geometries
   */
  void computeLabelling() {
    for (var node in graph.getNodes()) {
      node.getEdges().computeLabelling(arg);
    }
    mergeSymLabels();
    updateNodeLabelling();
  }

  /**
   * For nodes which have edges from only one Geometry incident on them,
   * the previous step will have left their dirEdges with no labelling for the other
   * Geometry.  However, the sym dirEdge may have a labelling for the other
   * Geometry, so merge the two labels.
   */
  void mergeSymLabels() {
    for (var node in graph.getNodes()) {
      node.getEdges().mergeSymLabels();
    }
  }

  void updateNodeLabelling() {
    // update the labels for nodes
    // The label for a node is updated from the edges incident on it
    // (Note that a node may have already been labelled
    // because it is a point in one of the input geometries)
    for (var node in graph.getNodes()) {
      Label lbl = node.getEdges().getLabel();
      node.getLabel().merge(lbl);
    }
  }

  /**
   * Incomplete nodes are nodes whose labels are incomplete.
   * (e.g. the location for one Geometry is null).
   * These are either isolated nodes,
   * or nodes which have edges from only a single Geometry incident on them.
   *
   * Isolated nodes are found because nodes in one graph which don't intersect
   * nodes in the other are not completely labelled by the initial process
   * of adding nodes to the nodeList.
   * To complete the labelling we need to check for nodes that lie in the
   * interior of edges, and in the interior of areas.
   * <p>
   * When each node labelling is completed, the labelling of the incident
   * edges is updated, to complete their labelling as well.
   */
  void labelIncompleteNodes() {
    // int nodeCount = 0;
    for (var n in graph.getNodes()) {
      Label label = n.getLabel();
      if (n.isIsolated()) {
        // nodeCount++;
        if (label.isNull(0))
          labelIncompleteNode(n, 0);
        else
          labelIncompleteNode(n, 1);
      }
      // now update the labelling for the DirectedEdges incident on this node
      n.getEdges().updateLabelling(label);
//n.print(System.out);
    }
  }

  /**
   * Label an isolated node with its relationship to the target geometry.
   */
  void labelIncompleteNode(Node n, int targetIndex) {
    int loc =
        ptLocator.locate(n.getCoordinate(), arg[targetIndex].getGeometry()!);
    if (n.getLabel() != null) {
      n.getLabel()!.setLocationWithIndex(targetIndex, loc);
    }
  }

  /**
   * Find all edges whose label indicates that they are in the result area(s),
   * according to the operation being performed.  Since we want polygon shells to be
   * oriented CW, choose dirEdges with the interior of the result on the RHS.
   * Mark them as being in the result.
   * Interior Area edges are the result of dimensional collapses.
   * They do not form part of the result area boundary.
   */
  void findResultAreaEdges(int opCode) {
    for (var e in graph.getEdgeEnds()) {
      var de = e as DirectedEdge;
      Label? label = de.getLabel();
      if (label!.isArea() &&
          !de.isInteriorAreaEdge() &&
          isResultOfOp(label.getLocationWithPosIndex(0, Position.RIGHT),
              label.getLocationWithPosIndex(1, Position.RIGHT), opCode)) {
        de.setInResult(true);
      }
    }
  }

  /**
   * If both a dirEdge and its sym are marked as being in the result, cancel
   * them out.
   */
  void cancelDuplicateResultEdges() {
    // remove any dirEdges whose sym is also included
    // (they "cancel each other out")
    for (var e in graph.getEdgeEnds()) {
      var de = e as DirectedEdge;
      DirectedEdge sym = de.getSym();
      if (de.isInResult() && sym.isInResult()) {
        de.setInResult(false);
        sym.setInResult(false);
      }
    }
  }

  /**
   * Tests if a point node should be included in the result or not.
   *
   * @param coord the point coordinate
   * @return true if the coordinate point is covered by a result Line or Area geometry
   */
  bool isCoveredByLA(Coordinate coord) {
    if (isCovered(coord, resultLineList)) return true;
    if (isCovered(coord, resultPolyList)) return true;
    return false;
  }

  /**
   * Tests if an L edge should be included in the result or not.
   *
   * @param coord the point coordinate
   * @return true if the coordinate point is covered by a result Area geometry
   */
  bool isCoveredByA(Coordinate coord) {
    if (isCovered(coord, resultPolyList)) return true;
    return false;
  }

  /**
   * @return true if the coord is located in the interior or boundary of
   * a geometry in the list.
   */
  bool isCovered(Coordinate coord, List geomList) {
    for (var geom in geomList) {
      int loc = ptLocator.locate(coord, geom);
      if (loc != Location.EXTERIOR) return true;
    }
    return false;
  }

  Geometry computeGeometry(
      List<Point> resultPointList,
      List<LineString> resultLineList,
      List<Polygon> resultPolyList,
      int opcode) {
    List<Geometry> geomList = List.empty(growable: true);
    // element geometries of the result are always in the order P,L,A
    geomList.addAll(resultPointList);
    geomList.addAll(resultLineList);
    geomList.addAll(resultPolyList);

    //*
    if (geomList.isEmpty)
      return createEmptyResult(
          opcode, arg[0].getGeometry()!, arg[1].getGeometry()!, geomFact);
    //*/

    // build the most specific geometry possible
    return geomFact.buildGeometry(geomList);
  }

  /**
   * Creates an empty result geometry of the appropriate dimension,
   * based on the given overlay operation and the dimensions of the inputs.
   * The created geometry is always an atomic geometry,
   * not a collection.
   * <p>
   * The empty result is constructed using the following rules:
   * <ul>
   * <li>{@link #INTERSECTION} - result has the dimension of the lowest input dimension
   * <li>{@link #UNION} - result has the dimension of the highest input dimension
   * <li>{@link #DIFFERENCE} - result has the dimension of the left-hand input
   * <li>{@link #SYMDIFFERENCE} - result has the dimension of the highest input dimension
   * (since the symmetric Difference is the union of the differences).
   * </ul>
   *
   * @param overlayOpCode the code for the overlay operation being performed
   * @param a an input geometry
   * @param b an input geometry
   * @param geomFact the geometry factory being used for the operation
   * @return an empty atomic geometry of the appropriate dimension
   */
  static Geometry createEmptyResult(
      int overlayOpCode, Geometry a, Geometry b, GeometryFactory geomFact) {
    int resultDim = resultDimension(overlayOpCode, a, b);

    /**
     * Handles resultSDim = -1, although should not happen
     */
    return geomFact.createEmpty(resultDim);
  }

  static int resultDimension(int opCode, Geometry g0, Geometry g1) {
    int dim0 = g0.getDimension();
    int dim1 = g1.getDimension();

    int resultDimension = -1;
    switch (opCode) {
      case INTERSECTION:
        resultDimension = math.min(dim0, dim1);
        break;
      case UNION:
        resultDimension = math.max(dim0, dim1);
        break;
      case DIFFERENCE:
        resultDimension = dim0;
        break;
      case SYMDIFFERENCE:
        /**
       * This result is chosen because
       * <pre>
       * SymDiff = Union(Diff(A, B), Diff(B, A)
       * </pre>
       * and Union has the dimension of the highest-dimension argument.
       */
        resultDimension = math.max(dim0, dim1);
        break;
    }
    return resultDimension;
  }
}

/**
 * Constructs {@link Point}s from the nodes of an overlay graph.
 * @version 1.7
 */
class PointBuilder {
  OverlayOp _op;
  GeometryFactory _geometryFactory;
  List<Point> resultPointList = List.empty(growable: true);

  PointBuilder(this._op, this._geometryFactory);

  /**
   * Computes the Point geometries which will appear in the result,
   * given the specified overlay operation.
   *
   * @return a list of the Points objects in the result
   */
  List<Point> build(int opCode) {
    _extractNonCoveredResultNodes(opCode);
    /**
     * It can happen that connected result nodes are still covered by
     * result geometries, so must perform this filter.
     * (For instance, this can happen during topology collapse).
     */
    return resultPointList;
  }

  /**
   * Determines nodes which are in the result, and creates {@link Point}s for them.
   *
   * This method determines nodes which are candidates for the result via their
   * labelling and their graph topology.
   *
   * @param opCode the overlay operation
   */
  void _extractNonCoveredResultNodes(int opCode) {
    // testing only
    //if (true) return resultNodeList;

    for (var n in _op.getGraph().getNodes()) {
      // filter out nodes which are known to be in the result
      if (n.isInResult()) continue;
      // if an incident edge is in the result, then the node coordinate is included already
      if (n.isIncidentEdgeInResult()) continue;
      if (n.getEdges().getDegree() == 0 || opCode == OverlayOp.INTERSECTION) {
        /**
         * For nodes on edges, only INTERSECTION can result in edge nodes being included even
         * if none of their incident edges are included
         */
        Label label = n.getLabel();
        if (OverlayOp.isResultOfOpLabel(label, opCode)) {
          _filterCoveredNodeToPoint(n);
        }
      }
    }
    //System.out.println("connectedResultNodes collected = " + connectedResultNodes.size());
  }

  /**
   * Converts non-covered nodes to Point objects and adds them to the result.
   *
   * A node is covered if it is contained in another element Geometry
   * with higher dimension (e.g. a node point might be contained in a polygon,
   * in which case the point can be eliminated from the result).
   *
   * @param n the node to test
   */
  void _filterCoveredNodeToPoint(Node n) {
    Coordinate coord = n.getCoordinate();
    if (!_op.isCoveredByLA(coord)) {
      Point pt = _geometryFactory.createPoint(coord);
      resultPointList.add(pt);
    }
  }
}

/**
 * Forms JTS LineStrings out of a the graph of {@link DirectedEdge}s
 * created by an {@link OverlayOp}.
 *
 * @version 1.7
 */
class LineBuilder {
  OverlayOp _op;
  GeometryFactory _geometryFactory;
  PointLocator _ptLocator;

  List<Edge> _lineEdgesList = List.empty(growable: true);
  List<LineString> _resultLineList = List.empty(growable: true);

  LineBuilder(this._op, this._geometryFactory, this._ptLocator);
  /**
   * @return a list of the LineStrings in the result of the specified overlay operation
   */
  List<LineString> build(int opCode) {
    _findCoveredLineEdges();
    _collectLines(opCode);
    buildLines(opCode);
    return _resultLineList;
  }

  /**
   * Find and mark L edges which are "covered" by the result area (if any).
   * L edges at nodes which also have A edges can be checked by checking
   * their depth at that node.
   * L edges at nodes which do not have A edges can be checked by doing a
   * point-in-polygon test with the previously computed result areas.
   */
  void _findCoveredLineEdges() {
    // first set covered for all L edges at nodes which have A edges too
    for (var node in _op.getGraph().getNodes()) {
      (node.getEdges()).findCoveredLineEdges();
    }

    /**
   * For all L edges which weren't handled by the above,
   * use a point-in-poly test to determine whether they are covered
   */
    for (var ee in _op.getGraph().getEdgeEnds()) {
      var de = ee as DirectedEdge;
      Edge e = de.getEdge();
      if (de.isLineEdge() && !e.isCoveredSet()) {
        bool isCovered = _op.isCoveredByA(de.getCoordinate()!);
        e.setCovered(isCovered);
      }
    }
  }

  void _collectLines(int opCode) {
    for (var e in _op.getGraph().getEdgeEnds()) {
      var de = e as DirectedEdge;
      _collectLineEdge(de, opCode, _lineEdgesList);
      _collectBoundaryTouchEdge(de, opCode, _lineEdgesList);
    }
  }

  /**
   * Collect line edges which are in the result.
   * Line edges are in the result if they are not part of
   * an area boundary, if they are in the result of the overlay operation,
   * and if they are not covered by a result area.
   *
   * @param de the directed edge to test
   * @param opCode the overlap operation
   * @param edges the list of included line edges
   */
  void _collectLineEdge(DirectedEdge de, int opCode, List edges) {
    Label? label = de.getLabel();
    Edge e = de.getEdge();
    // include L edges which are in the result
    if (de.isLineEdge()) {
      if (!de.isVisited() &&
          OverlayOp.isResultOfOpLabel(label!, opCode) &&
          !e.isCovered()) {
        edges.add(e);
        de.setVisitedEdge(true);
      }
    }
  }

  /**
   * Collect edges from Area inputs which should be in the result but
   * which have not been included in a result area.
   * This happens ONLY:
   * <ul>
   * <li>during an intersection when the boundaries of two
   * areas touch in a line segment
   * <li> OR as a result of a dimensional collapse.
   * </ul>
   */
  void _collectBoundaryTouchEdge(DirectedEdge de, int opCode, List edges) {
    Label? label = de.getLabel();
    if (de.isLineEdge()) return; // only interested in area edges
    if (de.isVisited()) return; // already processed
    if (de.isInteriorAreaEdge())
      return; // added to handle dimensional collapses
    if (de.getEdge().isInResult())
      return; // if the edge linework is already included, don't include it again

    // sanity check for labelling of result edgerings
    Assert.isTrue(!(de.isInResult() || de.getSym().isInResult()) ||
        !de.getEdge().isInResult());

    // include the linework if it's in the result of the operation
    if (OverlayOp.isResultOfOpLabel(label!, opCode) &&
        opCode == OverlayOp.INTERSECTION) {
      edges.add(de.getEdge());
      de.setVisitedEdge(true);
    }
  }

  void buildLines(int opCode) {
    for (var e in _lineEdgesList) {
      LineString line = _geometryFactory.createLineString(e.getCoordinates());
      _resultLineList.add(line);
      e.setInResult(true);
    }
  }

  void labelIsolatedLines(List edgesList) {
    for (var e in edgesList) {
      Label label = e.getLabel();
      if (e.isIsolated()) {
        if (label.isNull(0))
          _labelIsolatedLine(e, 0);
        else
          _labelIsolatedLine(e, 1);
      }
    }
  }

  /**
   * Label an isolated node with its relationship to the target geometry.
   */
  void _labelIsolatedLine(Edge e, int targetIndex) {
    int loc =
        _ptLocator.locate(e.getCoordinate()!, _op.getArgGeometry(targetIndex));
    e.getLabel()?.setLocationWithIndex(targetIndex, loc);
  }
}
