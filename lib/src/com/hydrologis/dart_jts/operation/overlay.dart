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
        EdgeRing shell = findEdgeRingContaining(hole, shellList);
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
  static EdgeRing findEdgeRingContaining(EdgeRing testEr, List shellList) {
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
    return minShell!;
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
