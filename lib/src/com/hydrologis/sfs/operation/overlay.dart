part of dart_sfs;

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
  MaximalEdgeRing(DirectedEdge start, GeometryFactory geometryFactory) : super(start, geometryFactory);

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
    DirectedEdge de = startDe;
    do {
      Node node = de.getNode();
      (node.getEdges() as DirectedEdgeStar).linkMinimalDirectedEdges(this);
      de = de.getNext();
    } while (de != startDe);
  }

  List buildMinimalRings() {
    List minEdgeRings = [];
    DirectedEdge de = startDe;
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
  MinimalEdgeRing(DirectedEdge start, GeometryFactory geometryFactory) : super(start, geometryFactory);

  DirectedEdge getNext(DirectedEdge de) {
    return de.getNextMin();
  }

  void setEdgeRing(DirectedEdge de, EdgeRing er) {
    de.setMinEdgeRing(er);
  }
}
