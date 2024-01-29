part of dart_jts;

/**
 * Extracts Point resultants from an overlay graph
 * created by an Intersection operation
 * between non-Point inputs.
 * Points may be created during intersection
 * if lines or areas touch one another at single points.
 * Intersection is the only overlay operation which can
 * result in Points from non-Point inputs.
 * <p>
 * Overlay operations where one or more inputs
 * are Points are handled via a different code path.
 *
 *
 * @author Martin Davis
 *
 * @see OverlayPoints
 *
 */
class IntersectionPointBuilder {
  GeometryFactory _geometryFactory;
  OverlayGraph _graph;
  List<Point> points = List.empty(growable: true);

  /**
   * Controls whether lines created by area topology collapses
   * to participate in the result computation.
   * True provides the original JTS semantics.
   */
  bool _isAllowCollapseLines = !OverlayNG.STRICT_MODE_DEFAULT;

  IntersectionPointBuilder(this._graph, this._geometryFactory);

  void setStrictMode(bool isStrictMode) {
    _isAllowCollapseLines = !isStrictMode;
  }

  List<Point> getPoints() {
    _addResultPoints();
    return points;
  }

  void _addResultPoints() {
    for (OverlayEdge nodeEdge in _graph.getNodeEdges()) {
      if (_isResultPoint(nodeEdge)) {
        Point pt =
            _geometryFactory.createPoint(nodeEdge.getCoordinate().copy());
        points.add(pt);
      }
    }
  }

  /**
   * Tests if a node is a result point.
   * This is the case if the node is incident on edges from both
   * inputs, and none of the edges are themselves in the result.
   *
   * @param nodeEdge an edge originating at the node
   * @return true if this node is a result point
   */
  bool _isResultPoint(OverlayEdge nodeEdge) {
    bool isEdgeOfA = false;
    bool isEdgeOfB = false;

    OverlayEdge edge = nodeEdge;
    do {
      if (edge.isInResult()) return false;
      OverlayLabel label = edge.getLabel();
      isEdgeOfA |= _isEdgeOf(label, 0);
      isEdgeOfB |= _isEdgeOf(label, 1);
      edge = edge.oNext() as OverlayEdge;
    } while (edge != nodeEdge);
    bool isNodeInBoth = isEdgeOfA && isEdgeOfB;
    return isNodeInBoth;
  }

  bool _isEdgeOf(OverlayLabel label, int i) {
    if (!_isAllowCollapseLines && label.isBoundaryCollapse()) return false;
    return label.isBoundary(i) || label.isLine();
  }
}
