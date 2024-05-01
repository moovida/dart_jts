part of dart_jts;

class OverlayEdge extends HalfEdge {
  /**
   * Creates a single OverlayEdge.
   *
   * @param pts
   * @param lbl
   * @param direction
   *
   * @return a new edge based on the given coordinates and direction
   */
  static OverlayEdge createEdge(
      List<Coordinate> pts, OverlayLabel lbl, bool direction) {
    Coordinate origin;
    Coordinate dirPt;
    if (direction) {
      origin = pts[0];
      dirPt = pts[1];
    } else {
      int ilast = pts.length - 1;
      origin = pts[ilast];
      dirPt = pts[ilast - 1];
    }
    return OverlayEdge(origin, dirPt, direction, lbl, pts);
  }

  static OverlayEdge createEdgePair(List<Coordinate> pts, OverlayLabel lbl) {
    OverlayEdge e0 = OverlayEdge.createEdge(pts, lbl, true);
    OverlayEdge e1 = OverlayEdge.createEdge(pts, lbl, false);
    e0.link(e1);
    return e0;
  }

  List<Coordinate> pts;

  /**
   * <code>true</code> indicates direction is forward along segString
   * <code>false</code> is reverse direction
   * The label must be interpreted accordingly.
   */
  bool direction;
  Coordinate dirPt;
  OverlayLabel label;

  bool isInResultArea = false;
  bool isInResultLine = false;
  bool isVisited = false;

  /**
   * Link to next edge in the result ring.
   * The origin of the edge is the dest of this edge.
   */
  OverlayEdge? nextResultEdge;

  OverlayEdgeRing? edgeRing;

  MaximalEdgeRingNG? maxEdgeRing;

  OverlayEdge? nextResultMaxEdge;

  OverlayEdge(Coordinate orig, this.dirPt, this.direction, this.label, this.pts)
      : super(orig);

  bool isForward() {
    return direction;
  }

  Coordinate directionPt() {
    return dirPt;
  }

  OverlayLabel getLabel() {
    return label;
  }

  int getLocation(int index, int position) {
    return label.getLocation(index, position, direction);
  }

  Coordinate getCoordinate() {
    return orig();
  }

  List<Coordinate> getCoordinates() {
    return pts;
  }

  List<Coordinate> getCoordinatesOriented() {
    if (direction) {
      return pts;
    }
    List<Coordinate> copy = List.empty(growable: true);
    copy.addAll(pts);
    CoordinateArrays.reverse(copy);
    return copy;
  }

  /**
   * Adds the coordinates of this edge to the given list,
   * in the direction of the edge.
   * Duplicate coordinates are removed
   * (which means that this is safe to use for a path
   * of connected edges in the topology graph).
   *
   * @param coords the coordinate list to add to
   */
  void addCoordinates(CoordinateList coords) {
    bool isFirstEdge = coords.size() > 0;
    if (direction) {
      int startIndex = 1;
      if (isFirstEdge) startIndex = 0;
      for (int i = startIndex; i < pts.length; i++) {
        coords.addCoord(pts[i], false);
      }
    } else {
      // is backward
      int startIndex = pts.length - 2;
      if (isFirstEdge) startIndex = pts.length - 1;
      for (int i = startIndex; i >= 0; i--) {
        coords.addCoord(pts[i], false);
      }
    }
  }

  /**
   * Gets the symmetric pair edge of this edge.
   *
   * @return the symmetric pair edge
   */
  OverlayEdge symOE() {
    return super.getSymEdge() as OverlayEdge;
  }

  /**
   * Gets the next edge CCW around the origin of this edge,
   * with the same origin.
   * If the origin vertex has degree 1 then this is the edge itself.
   *
   * @return the next edge around the origin
   */
  OverlayEdge? oNextOE() {
    return oNext() as OverlayEdge;
  }

  bool isInResultAreaBoth() {
    return isInResultArea && symOE().isInResultArea;
  }

  void unmarkFromResultAreaBoth() {
    isInResultArea = false;
    symOE().isInResultArea = false;
  }

  void markInResultArea() {
    isInResultArea = true;
  }

  void markInResultAreaBoth() {
    isInResultArea = true;
    symOE().isInResultArea = true;
  }

  void markInResultLine() {
    isInResultLine = true;
    symOE().isInResultLine = true;
  }

  bool isInResult() {
    return isInResultArea || isInResultLine;
  }

  bool isInResultEither() {
    return isInResult() || symOE().isInResult();
  }

  void setNextResult(OverlayEdge e) {
    // Assert: e.orig() == this.dest();
    nextResultEdge = e;
  }

  OverlayEdge? nextResult() {
    return nextResultEdge;
  }

  bool isResultLinked() {
    return nextResultEdge != null;
  }

  void setNextResultMax(OverlayEdge e) {
    // Assert: e.orig() == this.dest();
    nextResultMaxEdge = e;
  }

  OverlayEdge? nextResultMax() {
    return nextResultMaxEdge;
  }

  bool isResultMaxLinked() {
    return nextResultMaxEdge != null;
  }

  void markVisited() {
    isVisited = true;
  }

  void markVisitedBoth() {
    markVisited();
    symOE().markVisited();
  }

  void setEdgeRing(OverlayEdgeRing edgeRing) {
    this.edgeRing = edgeRing;
  }

  OverlayEdgeRing? getEdgeRing() {
    return edgeRing;
  }

  MaximalEdgeRingNG? getEdgeRingMax() {
    return maxEdgeRing;
  }

  void setEdgeRingMax(MaximalEdgeRingNG maximalEdgeRing) {
    maxEdgeRing = maximalEdgeRing;
  }

  String toString() {
    Coordinate orig = super.orig();
    Coordinate dest = super.dest();
    String dirPtStr = (pts.length > 2) ? ", " + directionPt().toString() : "";

    return "OE( " +
        orig.toString() +
        dirPtStr +
        " .. " +
        dest.toString() +
        " ) " +
        label.toStringLocal(direction) +
        resultSymbol() +
        " / Sym: " +
        symOE().getLabel().toStringLocal(symOE().direction) +
        symOE().resultSymbol();
  }

  String resultSymbol() {
    if (isInResultArea) return " resA";
    if (isInResultLine) return " resL";
    return "";
  }
}
