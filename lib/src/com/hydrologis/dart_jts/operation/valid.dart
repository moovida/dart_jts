part of dart_jts;

/**
 * Checks that a {@link GeometryGraph} representing an area
 * (a {@link Polygon} or {@link MultiPolygon} )
 * has consistent semantics for area geometries.
 * This check is required for any reasonable polygonal model
 * (including the OGC-SFS model, as well as models which allow ring self-intersection at single points)
 * <p>
 * Checks include:
 * <ul>
 * <li>test for rings which properly intersect
 * (but not for ring self-intersection, or intersections at vertices)
 * <li>test for consistent labelling at all node points
 * (this detects vertex intersections with invalid topology,
 * i.e. where the exterior side of an edge lies in the interior of the area)
 * <li>test for duplicate rings
 * </ul>
 * If an inconsistency is found the location of the problem
 * is recorded and is available to the caller.
 *
 * @version 1.7
 */
class ConsistentAreaTester {
  final LineIntersector li = new RobustLineIntersector();
  GeometryGraph geomGraph;
  RelateNodeGraph nodeGraph = new RelateNodeGraph();

  // the intersection point found (if any)
  Coordinate? invalidPoint;

  /**
   * Creates a new tester for consistent areas.
   *
   * @param geomGraph the topology graph of the area geometry
   */
  ConsistentAreaTester(this.geomGraph);

  /**
   * @return the intersection point, or <code>null</code> if none was found
   */
  Coordinate? getInvalidPoint() {
    return invalidPoint;
  }

  /**
   * Check all nodes to see if their labels are consistent with area topology.
   *
   * @return <code>true</code> if this area has a consistent node labelling
   */
  bool isNodeConsistentArea() {
    /**
     * To fully check validity, it is necessary to
     * compute ALL intersections, including self-intersections within a single edge.
     */
    SegmentIntersector intersector =
        geomGraph.computeSelfNodes3(li, true, true);
    /**
     * A proper intersection means that the area is not consistent.
     */
    if (intersector.hasProperIntersection()) {
      invalidPoint = intersector.getProperIntersectionPoint();
      return false;
    }

    nodeGraph.build(geomGraph);

    return isNodeEdgeAreaLabelsConsistent();
  }

  /**
   * Check all nodes to see if their labels are consistent.
   * If any are not, return false
   *
   * @return <code>true</code> if the edge area labels are consistent at this node
   */
  bool isNodeEdgeAreaLabelsConsistent() {
    for (Iterator nodeIt = nodeGraph.getNodeIterator(); nodeIt.moveNext();) {
      RelateNode node = nodeIt.current as RelateNode;
      if (!node.getEdges()!.isAreaLabelsConsistent(geomGraph)) {
        invalidPoint = node.getCoordinate().copy();
        return false;
      }
    }
    return true;
  }

  /**
   * Checks for two duplicate rings in an area.
   * Duplicate rings are rings that are topologically equal
   * (that is, which have the same sequence of points up to point order).
   * If the area is topologically consistent (determined by calling the
   * <code>isNodeConsistentArea</code>,
   * duplicate rings can be found by checking for EdgeBundles which contain
   * more than one EdgeEnd.
   * (This is because topologically consistent areas cannot have two rings sharing
   * the same line segment, unless the rings are equal).
   * The start point of one of the equal rings will be placed in
   * invalidPoint.
   *
   * @return true if this area Geometry is topologically consistent but has two duplicate rings
   */
  bool hasDuplicateRings() {
    for (Iterator nodeIt = nodeGraph.getNodeIterator(); nodeIt.moveNext();) {
      RelateNode node = nodeIt.current as RelateNode;
      for (Iterator i = node.getEdges()!.iterator(); i.moveNext();) {
        EdgeEndBundle eeb = i.current as EdgeEndBundle;
        if (eeb.getEdgeEnds().length > 1) {
          invalidPoint = eeb.getEdge().getCoordinateWithIndex(0);
          return true;
        }
      }
    }
    return false;
  }
}

/**
 * This class tests that the interior of an area {@link Geometry}
 * ( {@link Polygon}  or {@link MultiPolygon} )
 * is connected.
 * This can happen if:
 * <ul>
 * <li>a shell self-intersects
 * <li>one or more holes form a connected chain touching a shell at two different points
 * <li>one or more holes form a ring around a subset of the interior
 * </ul>
 * If a disconnected situation is found the location of the problem is recorded.
 *
 * @version 1.7
 */
class ConnectedInteriorTester {
  static Coordinate? findDifferentPoint(List<Coordinate> coord, Coordinate pt) {
    for (int i = 0; i < coord.length; i++) {
      if (!coord[i].equals(pt)) return coord[i];
    }
    return null;
  }

  GeometryFactory geometryFactory = new GeometryFactory.defaultPrecision();

  GeometryGraph geomGraph;

  // save a coordinate for any disconnected interior found
  // the coordinate will be somewhere on the ring surrounding the disconnected interior
  Coordinate? disconnectedRingcoord;

  ConnectedInteriorTester(this.geomGraph);

  Coordinate? getCoordinate() {
    return disconnectedRingcoord;
  }

  bool isInteriorsConnected() {
    // node the edges, in case holes touch the shell
    List splitEdges = [];
    geomGraph.computeSplitEdges(splitEdges);

    // form the edges into rings
    PlanarGraph graph = new PlanarGraph.withFactory(new OverlayNodeFactory());
    graph.addEdges(splitEdges);
    setInteriorEdgesInResult(graph);
    graph.linkResultDirectedEdges();
    List edgeRings = buildEdgeRings(graph.getEdgeEnds());

    /**
     * Mark all the edges for the edgeRings corresponding to the shells
     * of the input polygons.  Note only ONE ring gets marked for each shell.
     */
    visitShellInteriors(geomGraph.getGeometry(), graph);

    /**
     * If there are any unvisited shell edges
     * (i.e. a ring which is not a hole and which has the interior
     * of the parent area on the RHS)
     * this means that one or more holes must have split the interior of the
     * polygon into at least two pieces.  The polygon is thus invalid.
     */
    return !hasUnvisitedShellEdge(edgeRings);
  }

  void setInteriorEdgesInResult(PlanarGraph graph) {
    for (Iterator it = graph.getEdgeEnds().iterator; it.moveNext();) {
      DirectedEdge de = it.current as DirectedEdge;
      if (de.getLabel()!.getLocationWithPosIndex(0, Position.RIGHT) ==
          Location.INTERIOR) {
        de.setInResult(true);
      }
    }
  }

  /**
   * Form DirectedEdges in graph into Minimal EdgeRings.
   * (Minimal Edgerings must be used, because only they are guaranteed to provide
   * a correct isHole computation)
   */
  List buildEdgeRings(List dirEdges) {
    List edgeRings = [];
    for (Iterator it = dirEdges.iterator; it.moveNext();) {
      DirectedEdge de = it.current as DirectedEdge;
      // if this edge has not yet been processed
      if (de.isInResult() && de.getEdgeRing() == null) {
        MaximalEdgeRing er = new MaximalEdgeRing(de, geometryFactory);

        er.linkDirectedEdgesForMinimalEdgeRings();
        List minEdgeRings = er.buildMinimalRings();
        edgeRings.addAll(minEdgeRings);
      }
    }
    return edgeRings;
  }

  /**
   * Mark all the edges for the edgeRings corresponding to the shells
   * of the input polygons.
   * Only ONE ring gets marked for each shell - if there are others which remain unmarked
   * this indicates a disconnected interior.
   */
  void visitShellInteriors(Geometry? g, PlanarGraph graph) {
    if (g is Polygon) {
      visitInteriorRing(g.getExteriorRing(), graph);
    }
    if (g is MultiPolygon) {
      for (int i = 0; i < g.getNumGeometries(); i++) {
        Polygon p = g.getGeometryN(i) as Polygon;
        visitInteriorRing(p.getExteriorRing(), graph);
      }
    }
  }

  void visitInteriorRing(LineString ring, PlanarGraph graph) {
    if (ring.isEmpty()) return;
    List<Coordinate> pts = ring.getCoordinates();
    Coordinate pt0 = pts[0];
    /**
     * Find first point in coord list different to initial point.
     * Need special check since the first point may be repeated.
     */
    Coordinate pt1 = findDifferentPoint(pts, pt0)!;
    Edge e = graph.findEdgeInSameDirection(pt0, pt1)!;
    DirectedEdge de = graph.findEdgeEnd(e) as DirectedEdge;
    DirectedEdge? intDe = null;
    if (de.getLabel()!.getLocationWithPosIndex(0, Position.RIGHT) ==
        Location.INTERIOR) {
      intDe = de;
    } else if (de
            .getSym()
            .getLabel()!
            .getLocationWithPosIndex(0, Position.RIGHT) ==
        Location.INTERIOR) {
      intDe = de.getSym();
    }
    Assert.isTrue(intDe != null, "unable to find dirEdge with Interior on RHS");

    visitLinkedDirectedEdges(intDe!);
  }

  void visitLinkedDirectedEdges(DirectedEdge? start) {
    DirectedEdge? startDe = start;
    DirectedEdge? de = start;
    do {
      Assert.isTrue(de != null, "found null Directed Edge");
      de!.setVisited(true);
      de = de.getNext();
    } while (de != startDe);
  }

  /**
   * Check if any shell ring has an unvisited edge.
   * A shell ring is a ring which is not a hole and which has the interior
   * of the parent area on the RHS.
   * (Note that there may be non-hole rings with the interior on the LHS,
   * since the interior of holes will also be polygonized into CW rings
   * by the linkAllDirectedEdges() step)
   *
   * @return true if there is an unvisited edge in a non-hole ring
   */
  bool hasUnvisitedShellEdge(List edgeRings) {
    for (int i = 0; i < edgeRings.length; i++) {
      EdgeRing er = edgeRings[i] as EdgeRing;
      // don't check hole rings
      if (er.isHole()) continue;
      List edges = er.getEdges();
      DirectedEdge de = edges[0] as DirectedEdge;
      // don't check CW rings which are holes
      // (MD - this check may now be irrelevant)
      if (de.getLabel()!.getLocationWithPosIndex(0, Position.RIGHT) !=
          Location.INTERIOR) continue;

      /**
       * the edgeRing is CW ring which surrounds the INT of the area, so check all
       * edges have been visited.  If any are unvisited, this is a disconnected part of the interior
       */
      for (int j = 0; j < edges.length; j++) {
        de = edges[j] as DirectedEdge;
//Debug.print("visted? "); Debug.println(de);
        if (!de.isVisited()) {
//Debug.print("not visited "); Debug.println(de);
          disconnectedRingcoord = de.getCoordinate();
          return true;
        }
      }
    }
    return false;
  }
}

/**
 * Creates nodes for use in the {@link PlanarGraph}s constructed during
 * overlay operations.
 *
 * @version 1.7
 */
class OverlayNodeFactory extends NodeFactory {
  Node createNode(Coordinate coord) {
    return new Node(coord, new DirectedEdgeStar());
  }
}
