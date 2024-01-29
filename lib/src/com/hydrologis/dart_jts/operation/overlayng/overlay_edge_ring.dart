part of dart_jts;

class OverlayEdgeRing {
  OverlayEdge _startEdge;
  late LinearRing _ring;
  bool isHoleFlag = false;
  late List<Coordinate> _ringPts;
  IndexedPointInAreaLocator? _locator;
  OverlayEdgeRing? _shell;
  List<OverlayEdgeRing> holes = List.empty(
      growable: true); // a list of EdgeRings which are holes in this EdgeRing

  OverlayEdgeRing(this._startEdge, GeometryFactory geometryFactory) {
    _ringPts = _computeRingPts(_startEdge);
    _computeRing(_ringPts, geometryFactory);
  }

  LinearRing getRing() {
    return _ring;
  }

  Envelope _getEnvelope() {
    return _ring.getEnvelopeInternal();
  }

  /**
   * Tests whether this ring is a hole.
   * @return <code>true</code> if this ring is a hole
   */
  bool isHole() {
    return isHoleFlag;
  }

  /**
   * Sets the containing shell ring of a ring that has been determined to be a hole.
   *
   * @param shell the shell ring
   */
  void setShell(OverlayEdgeRing shell) {
    _shell = shell;
    if (_shell != null) shell.addHole(this);
  }

  /**
   * Tests whether this ring has a shell assigned to it.
   *
   * @return true if the ring has a shell
   */
  bool hasShell() {
    return _shell != null;
  }

  /**
   * Gets the shell for this ring.  The shell is the ring itself if it is not a hole, otherwise its parent shell.
   *
   * @return the shell for this ring
   */
  OverlayEdgeRing? getShell() {
    if (isHole()) return _shell;
    return this;
  }

  void addHole(OverlayEdgeRing ring) {
    holes.add(ring);
  }

  List<Coordinate> _computeRingPts(OverlayEdge start) {
    OverlayEdge edge = start;
    CoordinateList pts = CoordinateList();
    do {
      if (edge.getEdgeRing() == this)
        throw TopologyException("Edge visited twice during ring-building at " +
            edge.getCoordinate().toString());

      edge.addCoordinates(pts);
      edge.setEdgeRing(this);
      if (edge.nextResult() == null)
        throw TopologyException("Found null edge in ring");
      if (edge.nextResult() != null) {
        edge = edge.nextResult()!;
      }
    } while (edge != start);
    pts.closeRing();
    return pts.toCoordinateArray(true);
  }

  void _computeRing(List<Coordinate> ringPts, GeometryFactory geometryFactory) {
    // don't compute more than once
    _ring = geometryFactory.createLinearRing(ringPts);
    isHoleFlag = Orientation.isCCW(_ring.getCoordinates());
  }

  /**
   * Computes the list of coordinates which are contained in this ring.
   * The coordinates are computed once only and cached.
   *
   * @return an array of the {@link Coordinate}s in this ring
   */
  List<Coordinate> getCoordinates() {
    return _ringPts;
  }

  /**
   * Finds the innermost enclosing shell OverlayEdgeRing
   * containing this OverlayEdgeRing, if any.
   * The innermost enclosing ring is the <i>smallest</i> enclosing ring.
   * The algorithm used depends on the fact that:
   * <br>
   *  ring A contains ring B if envelope(ring A) contains envelope(ring B)
   * <br>
   * This routine is only safe to use if the chosen point of the hole
   * is known to be properly contained in a shell
   * (which is guaranteed to be the case if the hole does not touch its shell)
   * <p>
   * To improve performance of this function the caller should
   * make the passed shellList as small as possible (e.g.
   * by using a spatial index filter beforehand).
   *
   * @return containing EdgeRing or null if no containing EdgeRing is found
   */
  OverlayEdgeRing? findEdgeRingContaining(List<OverlayEdgeRing> erList) {
    OverlayEdgeRing? minContainingRing = null;

    for (OverlayEdgeRing edgeRing in erList) {
      if (edgeRing._contains(this)) {
        if (minContainingRing == null ||
            minContainingRing
                ._getEnvelope()
                .containsEnvelope(edgeRing._getEnvelope())) {
          minContainingRing = edgeRing;
        }
      }
    }
    return minContainingRing;
  }

  PointOnGeometryLocator _getLocator() {
    if (_locator == null) {
      _locator = IndexedPointInAreaLocator(getRing());
    }
    return _locator!;
  }

  int locate(Coordinate pt) {
    /**
   * Use an indexed point-in-polygon for performance
   */
    return _getLocator().locate(pt);
  }

  /**
   * Tests if an edgeRing is properly contained in this ring.
   * Relies on property that edgeRings never overlap (although they may
   * touch at single vertices).
   *
   * @param ring ring to test
   * @return true if ring is properly contained
   */
  bool _contains(OverlayEdgeRing ring) {
    // the test envelope must be properly contained
    // (guards against testing rings against themselves)
    Envelope env = _getEnvelope();
    Envelope testEnv = ring._getEnvelope();
    if (!env.containsEnvelope(testEnv)) return false;
    return _isPointInOrOut(ring);
  }

  bool _isPointInOrOut(OverlayEdgeRing ring) {
    // in most cases only one or two points will be checked
    for (Coordinate pt in ring.getCoordinates()) {
      int loc = locate(pt);
      if (loc == Location.INTERIOR) {
        return true;
      }
      if (loc == Location.EXTERIOR) {
        return false;
      }
      // pt is on BOUNDARY, so keep checking for a determining location
    }
    return false;
  }

  Coordinate getCoordinate() {
    return _ringPts[0];
  }

  /**
   * Computes the {@link Polygon} formed by this ring and any contained holes.
   *
   * @return the {@link Polygon} formed by this ring and its holes.
   */
  Polygon toPolygon(GeometryFactory factory) {
    List<LinearRing>? holeLR = List.empty(growable: true);
    for (int i = 0; i < holes.length; i++) {
      holeLR.add(holes[i].getRing());
    }
    Polygon poly = factory.createPolygon(_ring, holeLR);
    return poly;
  }

  OverlayEdge getEdge() {
    return _startEdge;
  }
}
