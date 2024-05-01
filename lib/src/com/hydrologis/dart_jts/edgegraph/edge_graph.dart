part of dart_jts;

/**
 * A graph comprised of {@link HalfEdge}s.
 * It supports tracking the vertices in the graph
 * via edges incident on them,
 * to allow efficient lookup of edges and vertices.
 * <p>
 * This class may be subclassed to use a
 * different subclass of HalfEdge,
 * by overriding {@link #createEdge(Coordinate)}.
 * If additional logic is required to initialize
 * edges then {@link EdgeGraph#addEdge(Coordinate, Coordinate)}
 * can be overridden as well.
 *
 * @author Martin Davis
 *
 */
class EdgeGraph {
  Map<Coordinate, HalfEdge> vertexMap = {};

  EdgeGraph() {}

  /**
   * Creates a single HalfEdge.
   * Override to use a different HalfEdge subclass.
   *
   * @param orig the origin location
   * @return a new HalfEdge with the given origin
   */
  HalfEdge createEdge(Coordinate orig) {
    return new HalfEdge(orig);
  }

  /**
   * Creates a HalfEge pair, using the HalfEdge type of the graph subclass.
   *
   * @param p0
   * @param p1
   * @return
   */
  HalfEdge create(Coordinate p0, Coordinate p1) {
    HalfEdge e0 = createEdge(p0);
    HalfEdge e1 = createEdge(p1);
    e0.link(e1);
    return e0;
  }

  /**
   * Adds an edge between the coordinates orig and dest
   * to this graph.
   * Only valid edges can be added (in particular, zero-length segments cannot be added)
   *
   * @param orig the edge origin location
   * @param dest the edge destination location.
   * @return the created edge
   * @return null if the edge was invalid and not added
   *
   * @see #isValidEdge(Coordinate, Coordinate)
   */
  HalfEdge? addEdge(Coordinate orig, Coordinate dest) {
    if (!isValidEdge(orig, dest)) return null;

    /**
     * Attempt to find the edge already in the graph.
     * Return it if found.
     * Otherwise, use a found edge with same origin (if any) to construct new edge.
     */
    HalfEdge? eAdj = vertexMap[orig];
    HalfEdge? eSame = null;
    if (eAdj != null) {
      eSame = eAdj.find(dest);
    }
    if (eSame != null) {
      return eSame;
    }

    HalfEdge e = insert(orig, dest, eAdj!);
    return e;
  }

  /**
   * Tests if the given coordinates form a valid edge (with non-zero length).
   *
   * @param orig the start coordinate
   * @param dest the end coordinate
   * @return true if the edge formed is valid
   */
  static bool isValidEdge(Coordinate orig, Coordinate dest) {
    int cmp = dest.compareTo(orig);
    return cmp != 0;
  }

  /**
   * Inserts an edge not already present into the graph.
   *
   * @param orig the edge origin location
   * @param dest the edge destination location
   * @param eAdj an existing edge with same orig (if any)
   * @return the created edge
   */
  HalfEdge insert(Coordinate orig, Coordinate dest, HalfEdge? eAdj) {
    // edge does not exist, so create it and insert in graph
    HalfEdge e = create(orig, dest);
    if (eAdj != null) {
      eAdj.insert(e);
    } else {
      vertexMap.putIfAbsent(orig, () => e);
    }

    HalfEdge? eAdjDest = vertexMap[dest];
    if (eAdjDest != null) {
      eAdjDest.insert(e.symEdge);
    } else {
      vertexMap.putIfAbsent(dest, () => e.symEdge);
    }
    return e;
  }

  /**
   * Gets all {@link HalfEdge}s in the graph.
   * Both edges of edge pairs are included.
   *
   * @return a collection of the graph edges
   */
  List<HalfEdge> getVertexEdges() {
    return vertexMap.values.toList();
  }

  /**
   * Finds an edge in this graph with the given origin
   * and destination, if one exists.
   *
   * @param orig the origin location
   * @param dest the destination location.
   * @return an edge with the given orig and dest, or null if none exists
   */
  HalfEdge? findEdge(Coordinate orig, Coordinate dest) {
    HalfEdge? e = vertexMap[orig];
    if (e == null) return null;
    return e.find(dest);
  }
}
