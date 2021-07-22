part of dart_jts;

/**
 * The action for the internal iterator for performing
 * overlap queries on a MonotoneChain
 *
 * @version 1.7
 */
class MonotoneChainOverlapAction {
  LineSegment overlapSeg1 = LineSegment.empty();
  LineSegment overlapSeg2 = LineSegment.empty();

  /**
   * This function can be overridden if the original chains are needed
   *
   * @param start1 the index of the start of the overlapping segment from mc1
   * @param start2 the index of the start of the overlapping segment from mc2
   */
  void overlap(MonotoneChainI mc1, int start1, MonotoneChainI mc2, int start2) {
    mc1.getLineSegment(start1, overlapSeg1);
    mc2.getLineSegment(start2, overlapSeg2);
    overlapLS(overlapSeg1, overlapSeg2);
  }

  /**
   * This is a convenience function which can be overridden to obtain the actual
   * line segments which overlap
   * @param seg1
   * @param seg2
   */
  void overlapLS(LineSegment seg1, LineSegment seg2) {}
}

/**
 * The action for the internal iterator for performing
 * envelope select queries on a MonotoneChain
 *
 * @version 1.7
 */
class MonotoneChainSelectAction {
  // these envelopes are used during the MonotoneChain search process
  //Envelope tempEnv1 = new Envelope();

  LineSegment selectedSegment = LineSegment.empty();

  /**
   * This method is overridden
   * to process a segment
   * in the context of the parent chain.
   *
   * @param mc the parent chain
   * @param startIndex the index of the start vertex of the segment being processed
   */
  void select(MonotoneChainI mc, int startIndex) {
    mc.getLineSegment(startIndex, selectedSegment);
    // call this routine in case select(segmenet) was overridden
    selectLS(selectedSegment);
  }

  /**
   * This is a convenience method which can be overridden to obtain the actual
   * line segment which is selected.
   *
   * @param seg
   */
  void selectLS(LineSegment seg) {}
}

/**
 * Monotone Chains are a way of partitioning the segments of a linestring to
 * allow for fast searching of intersections.
 * They have the following properties:
 * <ol>
 * <li>the segments within a monotone chain never intersect each other
 * <li>the envelope of any contiguous subset of the segments in a monotone chain
 * is equal to the envelope of the endpoints of the subset.
 * </ol>
 * Property 1 means that there is no need to test pairs of segments from within
 * the same monotone chain for intersection.
 * <p>
 * Property 2 allows
 * an efficient binary search to be used to find the intersection points of two monotone chains.
 * For many types of real-world data, these properties eliminate a large number of
 * segment comparisons, producing substantial speed gains.
 * <p>
 * One of the goals of this implementation of MonotoneChains is to be
 * as space and time efficient as possible. One design choice that aids this
 * is that a MonotoneChain is based on a subarray of a list of points.
 * This means that new arrays of points (potentially very large) do not
 * have to be allocated.
 * <p>
 *
 * MonotoneChains support the following kinds of queries:
 * <ul>
 * <li>Envelope select: determine all the segments in the chain which
 * intersect a given envelope
 * <li>Overlap: determine all the pairs of segments in two chains whose
 * envelopes overlap
 * </ul>
 *
 * This implementation of MonotoneChains uses the concept of internal iterators
 * ({@link MonotoneChainSelectAction} and {@link MonotoneChainOverlapAction})
 * to return the results for queries.
 * This has time and space advantages, since it
 * is not necessary to build lists of instantiated objects to represent the segments
 * returned by the query.
 * Queries made in this manner are thread-safe.
 *
 * MonotoneChains support being assigned an integer id value
 * to provide a total ordering for a set of chains.
 * This can be used during some kinds of processing to
 * avoid redundant comparisons
 * (i.e. by comparing only chains where the first id is less than the second).
 *
 * @version 1.7
 */
class MonotoneChainI {
  List<Coordinate> pts;
  int start;
  int end;
  Envelope? env = null;
  Object? context; // user-defined information
  int id = 0; // useful for optimizing chain comparisons

  /**
   * Creates a new MonotoneChain based on the given array of points.
   * @param pts the points containing the chain
   * @param start the index of the first coordinate in the chain
   * @param end the index of the last coordinate in the chain
   * @param context a user-defined data object
   */
  MonotoneChainI(this.pts, this.start, this.end, this.context);

  /**
   * Sets the id of this chain.
   * Useful for assigning an ordering to a set of
   * chains, which can be used to avoid redundant processing.
   *
   * @param id an id value
   */
  void setId(int id) {
    this.id = id;
  }

  /**
   * Gets the id of this chain.
   *
   * @return the id value
   */
  int getId() {
    return id;
  }

  /**
   * Gets the user-defined context data value.
   *
   * @return a data value
   */
  Object? getContext() {
    return context;
  }

  /**
   * Gets the envelope of the chain.
   *
   * @return the envelope of the chain
   */
  Envelope getEnvelope() {
    if (env == null) {
      /**
       * The monotonicity property allows fast envelope determination
       */
      Coordinate p0 = pts[start];
      Coordinate p1 = pts[end];
      env = new Envelope.fromCoordinates(p0, p1);
    }
    return env!;
  }

  /**
   * Gets the index of the start of the monotone chain
   * in the underlying array of points.
   *
   * @return the start index of the chain
   */
  int getStartIndex() {
    return start;
  }

  /**
   * Gets the index of the end of the monotone chain
   * in the underlying array of points.
   *
   * @return the end index of the chain
   */
  int getEndIndex() {
    return end;
  }

  /**
   * Gets the line segment starting at <code>index</code>
   *
   * @param index index of segment
   * @param ls line segment to extract into
   */
  void getLineSegment(int index, LineSegment ls) {
    ls.p0 = pts[index];
    ls.p1 = pts[index + 1];
  }

  /**
   * Return the subsequence of coordinates forming this chain.
   * Allocates a new array to hold the Coordinates
   */
  List<Coordinate> getCoordinates() {
    List<Coordinate> coord = []; //..length = end - start + 1;
    // int index = 0;
    for (int i = start; i <= end; i++) {
      coord.add(pts[i]);
      // coord[index++] = pts[i];
    }
    return coord;
  }

  /**
   * Determine all the line segments in the chain whose envelopes overlap
   * the searchEnvelope, and process them.
   * <p>
   * The monotone chain search algorithm attempts to optimize
   * performance by not calling the select action on chain segments
   * which it can determine are not in the search envelope.
   * However, it *may* call the select action on segments
   * which do not intersect the search envelope.
   * This saves on the overhead of checking envelope intersection
   * each time, since clients may be able to do this more efficiently.
   *
   * @param searchEnv the search envelope
   * @param mcs the select action to execute on selected segments
   */
  void select(Envelope searchEnv, MonotoneChainSelectAction mcs) {
    computeSelect(searchEnv, start, end, mcs);
  }

  void computeSelect(
      Envelope searchEnv, int start0, int end0, MonotoneChainSelectAction mcs) {
    Coordinate p0 = pts[start0];
    Coordinate p1 = pts[end0];

//Debug.println("trying:" + p0 + p1 + " [ " + start0 + ", " + end0 + " ]");
    // terminating condition for the recursion
    if (end0 - start0 == 1) {
      //Debug.println("computeSelect:" + p0 + p1);
      mcs.select(this, start0);
      return;
    }
    // nothing to do if the envelopes don't overlap
    if (!searchEnv.intersectsEnvelopeCoordinates(p0, p1)) return;

    // the chains overlap, so split each in half and iterate  (binary search)
    int mid = ((start0 + end0) / 2).toInt();

    // Assert: mid != start or end (since we checked above for end - start <= 1)
    // check terminating conditions before recursing
    if (start0 < mid) {
      computeSelect(searchEnv, start0, mid, mcs);
    }
    if (mid < end0) {
      computeSelect(searchEnv, mid, end0, mcs);
    }
  }

  /**
   * Determine all the line segments in two chains which may overlap, and process them.
   * <p>
   * The monotone chain search algorithm attempts to optimize
   * performance by not calling the overlap action on chain segments
   * which it can determine do not overlap.
   * However, it *may* call the overlap action on segments
   * which do not actually interact.
   * This saves on the overhead of checking intersection
   * each time, since clients may be able to do this more efficiently.
   *
   * @param searchEnv the search envelope
   * @param mco the overlap action to execute on selected segments
   */
  void computeOverlaps3(MonotoneChainI mc, MonotoneChainOverlapAction mco) {
    computeOverlaps6(start, end, mc, mc.start, mc.end, mco);
  }

  /**
   * Uses an efficient mutual binary search strategy
   * to determine which pairs of chain segments
   * may overlap, and calls the given overlap action on them.
   *
   * @param start0 the start index of this chain section
   * @param end0 the end index of this chain section
   * @param mc the target monotone chain
   * @param start1 the start index of the target chain section
   * @param end1 the end index of the target chain section
   * @param mco the overlap action to execute on selected segments
   */
  void computeOverlaps6(int start0, int end0, MonotoneChainI mc, int start1,
      int end1, MonotoneChainOverlapAction mco) {
//Debug.println("computeIntersectsForChain:" + p00 + p01 + p10 + p11);
    // terminating condition for the recursion
    if (end0 - start0 == 1 && end1 - start1 == 1) {
      mco.overlap(this, start0, mc, start1);
      return;
    }
    // nothing to do if the envelopes of these subchains don't overlap
    if (!overlaps(start0, end0, mc, start1, end1)) return;

    // the chains overlap, so split each in half and iterate  (binary search)
    int mid0 = ((start0 + end0) / 2).toInt();
    int mid1 = ((start1 + end1) / 2).toInt();

    // Assert: mid != start or end (since we checked above for end - start <= 1)
    // check terminating conditions before recursing
    if (start0 < mid0) {
      if (start1 < mid1) computeOverlaps6(start0, mid0, mc, start1, mid1, mco);
      if (mid1 < end1) computeOverlaps6(start0, mid0, mc, mid1, end1, mco);
    }
    if (mid0 < end0) {
      if (start1 < mid1) computeOverlaps6(mid0, end0, mc, start1, mid1, mco);
      if (mid1 < end1) computeOverlaps6(mid0, end0, mc, mid1, end1, mco);
    }
  }

  /**
   * Tests whether the envelope of a section of the chain
   * overlaps (intersects) the envelope of a section of another target chain.
   * This test is efficient due to the monotonicity property
   * of the sections (i.e. the envelopes can be are determined
   * from the section endpoints
   * rather than a full scan).
   *
   * @param start0 the start index of this chain section
   * @param end0 the end index of this chain section
   * @param mc the target monotone chain
   * @param start1 the start index of the target chain section
   * @param end1 the end index of the target chain section
   * @return true if the section envelopes overlap
   */
  bool overlaps(int start0, int end0, MonotoneChainI mc, int start1, int end1) {
    return Envelope.intersectsEnvelopeCoords(
        pts[start0], pts[end0], mc.pts[start1], mc.pts[end1]);
  }
}

/**
 * Constructs {@link MonotoneChain}s
 * for sequences of {@link Coordinate}s.
 *
 * @version 1.7
 */
class MonotoneChainBuilder {
  /**
   * Computes a list of the {@link MonotoneChain}s
   * for a list of coordinates.
   *
   * @param pts the list of points to compute chains for
   * @return a list of the monotone chains for the points
   */
  static List getChains(List<Coordinate> pts) {
    return getChainsWithContext(pts, null);
  }

  /**
   * Computes a list of the {@link MonotoneChain}s
   * for a list of coordinates,
   * attaching a context data object to each.
   *
   * @param pts the list of points to compute chains for
   * @param context a data object to attach to each chain
   * @return a list of the monotone chains for the points
   */
  static List getChainsWithContext(List<Coordinate> pts, Object? context) {
    List mcList = [];
    int chainStart = 0;
    do {
      int chainEnd = findChainEnd(pts, chainStart);
      MonotoneChainI mc =
          new MonotoneChainI(pts, chainStart, chainEnd, context);
      mcList.add(mc);
      chainStart = chainEnd;
    } while (chainStart < pts.length - 1);
    return mcList;
  }

  /**
   * Finds the index of the last point in a monotone chain
   * starting at a given point.
   * Repeated points (0-length segments) are included
   * in the monotone chain returned.
   *
   * @param pts the points to scan
   * @param start the index of the start of this chain
   * @return the index of the last point in the monotone chain
   * starting at <code>start</code>.
   */
  static int findChainEnd(List<Coordinate> pts, int start) {
    int safeStart = start;
    // skip any zero-length segments at the start of the sequence
    // (since they cannot be used to establish a quadrant)
    while (safeStart < pts.length - 1 &&
        pts[safeStart].equals2D(pts[safeStart + 1])) {
      safeStart++;
    }
    // check if there are NO non-zero-length segments
    if (safeStart >= pts.length - 1) {
      return pts.length - 1;
    }
    // determine overall quadrant for chain (which is the starting quadrant)
    int chainQuad =
        Quadrant.quadrantFromCoords(pts[safeStart], pts[safeStart + 1]);
    int last = start + 1;
    while (last < pts.length) {
      // skip zero-length segments, but include them in the chain
      if (!pts[last - 1].equals2D(pts[last])) {
        // compute quadrant for next possible segment in chain
        int quad = Quadrant.quadrantFromCoords(pts[last - 1], pts[last]);
        if (quad != chainQuad) break;
      }
      last++;
    }
    return last - 1;
  }
}
