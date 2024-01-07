part of dart_jts;

class MaximalEdgeRingNG {
  static const int STATE_FIND_INCOMING = 1;
  static const int STATE_LINK_OUTGOING = 2;

  /**
   * Traverses the star of edges originating at a node
   * and links consecutive result edges together
   * into <b>maximal</b> edge rings.
   * To link two edges the <code>resultNextMax</code> pointer
   * for an <b>incoming</b> result edge
   * is set to the next <b>outgoing</b> result edge.
   * <p>
   * Edges are linked when:
   * <ul>
   * <li>they belong to an area (i.e. they have sides)
   * <li>they are marked as being in the result
   * </ul>
   * <p>
   * Edges are linked in CCW order
   * (which is the order they are linked in the underlying graph).
   * This means that rings have their face on the Right
   * (in other words,
   * the topological location of the face is given by the RHS label of the DirectedEdge).
   * This produces rings with CW orientation.
   * <p>
   * PRECONDITIONS:
   * - This edge is in the result
   * - This edge is not yet linked
   * - The edge and its sym are NOT both marked as being in the result
   */
  static void linkResultAreaMaxRingAtNode(OverlayEdge nodeEdge) {
    Assert.isTrue(nodeEdge.isInResultArea, "Attempt to link non-result edge");
    // assertion is only valid if building a polygonal geometry (ie not a coverage)
    //Assert.isTrue(! nodeEdge.symOE().isInResultArea(), "Found both half-edges in result");

    /**
     * Since the node edge is an out-edge,
     * make it the last edge to be linked
     * by starting at the next edge.
     * The node edge cannot be an in-edge as well,
     * but the next one may be the first in-edge.
     */
    OverlayEdge? endOut = nodeEdge.oNextOE();
    OverlayEdge? currOut = endOut;
    int state = STATE_FIND_INCOMING;
    OverlayEdge? currResultIn = null;
    do {
      /**
       * If an edge is linked this node has already been processed
       * so can skip further processing
       */
      if (currResultIn != null && currResultIn.isResultMaxLinked()) return;

      switch (state) {
        case STATE_FIND_INCOMING:
          OverlayEdge currIn = currOut!.symOE();
          if (!currIn.isInResultArea) break;
          currResultIn = currIn;
          state = STATE_LINK_OUTGOING;
          //Debug.println("Found result in-edge:  " + currResultIn);
          break;
        case STATE_LINK_OUTGOING:
          if (!currOut!.isInResultArea) break;
          // link the in edge to the out edge
          currResultIn!.setNextResultMax(currOut);
          state = STATE_FIND_INCOMING;
          //Debug.println("Linked Max edge:  " + currResultIn + " -> " + currOut);
          break;
      }
      currOut = currOut!.oNextOE();
    } while (currOut != endOut);
    //Debug.println("AFTER: " + toString(nodeEdge));
    if (state == STATE_LINK_OUTGOING) {
//Debug.print(firstOut == null, this);
      throw TopologyException("no outgoing edge found");
    }
  }

  OverlayEdge startEdge;

  MaximalEdgeRingNG(this.startEdge) {
    _attachEdges(startEdge);
  }

  void _attachEdges(OverlayEdge startEdge) {
    OverlayEdge? edge = startEdge;
    do {
      if (edge == null) throw TopologyException("Ring edge is null");
      if (edge.getEdgeRingMax() == this)
        throw TopologyException(
            "Ring edge visited twice at " + edge.getCoordinate().toString());
      if (edge.nextResultMax() == null) {
        throw TopologyException("Ring edge missing at");
      }
      edge.setEdgeRingMax(this);
      edge = edge.nextResultMax();
    } while (edge != startEdge);
  }

  List<OverlayEdgeRing> buildMinimalRings(GeometryFactory geometryFactory) {
    _linkMinimalRings();

    List<OverlayEdgeRing> minEdgeRings = List.empty(growable: true);
    OverlayEdge? e = startEdge;
    do {
      if (e!.getEdgeRing() == null) {
        OverlayEdgeRing minEr = OverlayEdgeRing(e, geometryFactory);
        minEdgeRings.add(minEr);
      }
      e = e.nextResultMax();
    } while (e != startEdge);
    return minEdgeRings;
  }

  void _linkMinimalRings() {
    OverlayEdge? e = startEdge;
    do {
      linkMinRingEdgesAtNode(e!, this);
      e = e.nextResultMax();
    } while (e != startEdge);
  }

  /**
   * Links the edges of a {@link MaximalEdgeRingNG} around this node
   * into minimal edge rings ({@link OverlayEdgeRing}s).
   * Minimal ring edges are linked in the opposite orientation (CW)
   * to the maximal ring.
   * This changes self-touching rings into a two or more separate rings,
   * as per the OGC SFS polygon topology semantics.
   * This relinking must be done to each max ring separately,
   * rather than all the node result edges, since there may be
   * more than one max ring incident at the node.
   *
   * @param nodeEdge an edge originating at this node
   * @param maxRing the maximal ring to link
   */
  static void linkMinRingEdgesAtNode(
      OverlayEdge nodeEdge, MaximalEdgeRingNG maxRing) {
    //Assert.isTrue(nodeEdge.isInResult(), "Attempt to link non-result edge");

    /**
     * The node edge is an out-edge,
     * so it is the first edge linked
     * with the next CCW in-edge
     */
    OverlayEdge endOut = nodeEdge;
    OverlayEdge? currMaxRingOut = endOut;
    OverlayEdge? currOut = endOut.oNextOE();
    do {
      if (_isAlreadyLinked(currOut!.symOE(), maxRing)) return;

      if (currMaxRingOut == null) {
        currMaxRingOut = selectMaxOutEdge(currOut, maxRing);
      } else {
        currMaxRingOut = linkMaxInEdge(currOut, currMaxRingOut, maxRing);
      }
      currOut = currOut.oNextOE();
    } while (currOut != endOut);
    //Debug.println("AFTER: " + toString(nodeEdge));
    if (currMaxRingOut != null) {
      throw TopologyException("Unmatched edge found during min-ring linking");
    }
  }

  /**
   * Tests if an edge of the maximal edge ring is already linked into
   * a minimal {@link OverlayEdgeRing}.  If so, this node has already been processed
   * earlier in the maximal edgering linking scan.
   *
   * @param edge an edge of a maximal edgering
   * @param maxRing the maximal edgering
   * @return true if the edge has already been linked into a minimal edgering.
   */
  static bool _isAlreadyLinked(OverlayEdge edge, MaximalEdgeRingNG maxRing) {
    bool isLinked = edge.getEdgeRingMax() == maxRing && edge.isResultLinked();
    return isLinked;
  }

  static OverlayEdge? selectMaxOutEdge(
      OverlayEdge currOut, MaximalEdgeRingNG maxEdgeRing) {
    // select if currOut edge is part of this max ring
    if (currOut.getEdgeRingMax() == maxEdgeRing) return currOut;
    // otherwise skip this edge
    return null;
  }

  static OverlayEdge? linkMaxInEdge(OverlayEdge currOut,
      OverlayEdge currMaxRingOut, MaximalEdgeRingNG maxEdgeRing) {
    OverlayEdge currIn = currOut.symOE();
    // currIn is not in this max-edgering, so keep looking
    if (currIn.getEdgeRingMax() != maxEdgeRing) return currMaxRingOut;

    //Debug.println("Found result in-edge:  " + currIn);

    currIn.setNextResult(currMaxRingOut);
    //Debug.println("Linked Min Edge:  " + currIn + " -> " + currMaxRingOut);
    // return null to indicate to scan for the next max-ring out-edge
    return null;
  }

  String toString() {
    List<Coordinate> pts = _getCoordinates();
    return WKTWriter.toLineString(pts);
  }

  List<Coordinate> _getCoordinates() {
    CoordinateList coords = CoordinateList();
    OverlayEdge? edge = startEdge;
    do {
      coords.add(edge!.orig());
      if (edge.nextResultMax() == null) {
        break;
      }
      edge = edge.nextResultMax();
    } while (edge != startEdge);
    // add last coordinate
    coords.add(edge!.dest());
    return coords.toCoordinateArray(true);
  }
}
