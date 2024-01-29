part of dart_jts;

/**
 * A {@link HalfEdge} which supports
 * marking edges with a bool flag.
 * Useful for algorithms which perform graph traversals.
 *
 * @author Martin Davis
 *
 */
class MarkHalfEdge extends HalfEdge {
  /**
   * Tests whether the given edge is marked.
   *
   * @param e the edge to test
   * @return true if the edge is marked
   */
  static bool isHalfEdgeMarked(HalfEdge e) {
    return (e as MarkHalfEdge).isMarked;
  }

  /**
   * Marks the given edge.
   *
   * @param e the edge to mark
   */
  static void markHalfEdge(HalfEdge e) {
    (e as MarkHalfEdge).mark();
  }

  /**
   * Sets the mark for the given edge to a bool value.
   *
   * @param e the edge to set
   * @param isMarked the mark value
   */
  static void setMarkHalfEdge(HalfEdge e, bool isMarked) {
    (e as MarkHalfEdge).setMark(isMarked);
  }

  /**
   * Sets the mark for the given edge pair to a bool value.
   *
   * @param e an edge of the pair to update
   * @param isMarked the mark value to set
   */
  static void setMarkBoth(HalfEdge e, bool isMarked) {
    final ee = e as MarkHalfEdge;
    ee.setMark(isMarked);
    final eee = ee.symEdge as MarkHalfEdge;
    eee.setMark(isMarked);
  }

  /**
   * Marks the edges in a pair.
   *
   * @param e an edge of the pair to mark
   */
  static void markBoth(HalfEdge e) {
    final ee = e as MarkHalfEdge;
    ee.mark();
    final eee = ee.symEdge as MarkHalfEdge;
    eee.mark();
  }

  bool isMarked = false;

  /**
   * Creates a new marked edge.
   *
   * @param orig the coordinate of the edge origin
   */
  MarkHalfEdge(Coordinate orig) : super(orig);

  /**
   * Marks this edge.
   *
   */
  void mark() {
    isMarked = true;
  }

  /**
   * Sets the value of the mark on this edge.
   *
   * @param isMarked the mark value to set
   */
  void setMark(bool isMarked) {
    this.isMarked = isMarked;
  }
}
