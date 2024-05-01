part of dart_jts;

/**
 * A planar graph of edges, representing
 * the topology resulting from an overlay operation.
 * Each source edge is represented
 * by a pair of {@link OverlayEdge}s, with opposite (symmetric) orientation.
 * The pair of OverlayEdges share the edge coordinates
 * and a single {@link OverlayLabel}.
 *
 * @author Martin Davis
 *
 */
class OverlayGraph {
  List<OverlayEdge> edges = List.empty(growable: true);
  Map<Coordinate, OverlayEdge> nodeMap = HashMap<Coordinate, OverlayEdge>();

  /**
   * Creates an empty graph.
   */
  OverlayGraph();

  /**
   * Gets the set of edges in this graph.
   * Only one of each symmetric pair of OverlayEdges is included.
   * The opposing edge can be found by using {@link OverlayEdge#sym()}.
   *
   * @return the collection of representative edges in this graph
   */
  List<OverlayEdge> getEdges() {
    return edges;
  }

  /**
   * Gets the collection of edges representing the nodes in this graph.
   * For each star of edges originating at a node
   * a single representative edge is included.
   * The other edges around the node can be found by following the next and prev links.
   *
   * @return the collection of representative node edges
   */
  List<OverlayEdge> getNodeEdges() {
    return nodeMap.values.toList();
  }

  /**
   * Gets an edge originating at the given node point.
   *
   * @param nodePt the node coordinate to query
   * @return an edge originating at the point, or null if none exists
   */
  OverlayEdge? getNodeEdge(Coordinate nodePt) {
    return nodeMap[nodePt];
  }

  /**
   * Gets the representative edges marked as being in the result area.
   *
   * @return the result area edges
   */
  List<OverlayEdge> getResultAreaEdges() {
    List<OverlayEdge> resultEdges = List.empty(growable: true);
    for (OverlayEdge edge in getEdges()) {
      if (edge.isInResultArea) {
        resultEdges.add(edge);
      }
    }
    return resultEdges;
  }

  /**
   * Adds a new edge to this graph, for the given linework and topology information.
   * A pair of {@link OverlayEdge}s with opposite (symmetric) orientation is added,
   * sharing the same {@link OverlayLabel}.
   *
   * @param pts the edge vertices
   * @param label the edge topology information
   * @return the created graph edge with same orientation as the linework
   */
  OverlayEdge addEdge(List<Coordinate> pts, OverlayLabel label) {
    //if (! isValidEdge(orig, dest)) return null;
    OverlayEdge e = OverlayEdge.createEdgePair(pts, label);
    //Debug.println("added edge: " + e);
    insert(e);
    insert(e.symOE());
    return e;
  }

  /**
   * Inserts a single half-edge into the graph.
   * The sym edge must also be inserted.
   *
   * @param e the half-edge to insert
   */
  void insert(OverlayEdge e) {
    edges.add(e);

    /**
     * If the edge origin node is already in the graph,
     * insert the edge into the star of edges around the node.
     * Otherwise, add a new node for the origin.
     */
    OverlayEdge? nodeEdge = nodeMap[e.orig()];
    if (nodeEdge != null) {
      nodeEdge.insert(e);
    } else {
      nodeMap.putIfAbsent(e.orig(), () => e);
    }
  }
}
