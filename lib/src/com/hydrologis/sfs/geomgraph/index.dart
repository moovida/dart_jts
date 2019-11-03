part of dart_sfs;


/**
 * Computes the intersection of line segments,
 * and adds the intersection to the edges containing the segments.
 *
 * @version 1.7
 */
 class SegmentIntersector
{

   static bool isAdjacentSegments(int i1, int i2)
  {
    return (i1 - i2).abs() == 1;
  }

  /**
   * These variables keep track of what types of intersections were
   * found during ALL edges that have been intersected.
   */
   bool _hasIntersection = false;
   bool _hasProper = false;
   bool _hasProperInterior = false;
  // the proper intersection point found
   Coordinate _properIntersectionPoint = null;

   LineIntersector _li;
   bool _includeProper;
   bool _recordIsolated;
   bool _isSelfIntersection;
  // bool intersectionFound;
   int _numIntersections = 0;

  // testing only
   int _numTests = 0;

   List<List> bdyNodes;
   bool _isDone = false;
   bool _isDoneWhenProperInt = false;


   SegmentIntersector(LineIntersector li,  bool includeProper, bool recordIsolated)
  {
    this._li = li;
    this._includeProper = includeProper;
    this._recordIsolated = recordIsolated;
  }

   void setBoundaryNodes( List bdyNodes0,
       List bdyNodes1)
  {
    bdyNodes = List(2);
    bdyNodes[0] = bdyNodes0;
    bdyNodes[1] = bdyNodes1;
  }

   void setIsDoneIfProperInt(bool isDoneWhenProperInt) {
    this._isDoneWhenProperInt = isDoneWhenProperInt;
  }

   bool isDone() {
    return _isDone;
  }
  /**
   * @return the proper intersection point, or <code>null</code> if none was found
   */
   Coordinate getProperIntersectionPoint()  {    return _properIntersectionPoint;  }

   bool hasIntersection() { return _hasIntersection; }
  /**
   * A proper intersection is an intersection which is interior to at least two
   * line segments.  Note that a proper intersection is not necessarily
   * in the interior of the entire Geometry, since another edge may have
   * an endpoint equal to the intersection, which according to SFS semantics
   * can result in the point being on the Boundary of the Geometry.
   */
   bool hasProperIntersection() { return _hasProper; }
  /**
   * A proper interior intersection is a proper intersection which is <b>not</b>
   * contained in the set of boundary nodes set for this SegmentIntersector.
   */
   bool hasProperInteriorIntersection() { return _hasProperInterior; }


  /**
   * A trivial intersection is an apparent self-intersection which in fact
   * is simply the point shared by adjacent line segments.
   * Note that closed edges require a special check for the point shared by the beginning
   * and end segments.
   */
   bool isTrivialIntersection(Edge e0, int segIndex0, Edge e1, int segIndex1)
  {
    if (e0 == e1) {
      if (_li.getIntersectionNum() == 1) {
        if (isAdjacentSegments(segIndex0, segIndex1))
          return true;
        if (e0.isClosed()) {
          int maxSegIndex = e0.getNumPoints() - 1;
          if (    (segIndex0 == 0 && segIndex1 == maxSegIndex)
              ||  (segIndex1 == 0 && segIndex0 == maxSegIndex) ) {
            return true;
          }
        }
      }
    }
    return false;
  }

  /**
   * This method is called by clients of the EdgeIntersector class to test for and add
   * intersections for two segments of the edges being intersected.
   * Note that clients (such as MonotoneChainEdges) may choose not to intersect
   * certain pairs of segments for efficiency reasons.
   */
   void addIntersections(
      Edge e0,  int segIndex0,
      Edge e1,  int segIndex1
      )
  {
    if (e0 == e1 && segIndex0 == segIndex1) return;
    _numTests++;
    Coordinate p00 = e0.getCoordinates()[segIndex0];
    Coordinate p01 = e0.getCoordinates()[segIndex0 + 1];
    Coordinate p10 = e1.getCoordinates()[segIndex1];
    Coordinate p11 = e1.getCoordinates()[segIndex1 + 1];

    _li.computeIntersection(p00, p01, p10, p11);
//if (li.hasIntersection() && li.isProper()) Debug.println(li);
    /**
     *  Always record any non-proper intersections.
     *  If includeProper is true, record any proper intersections as well.
     */
    if (_li.hasIntersection()) {
      if (_recordIsolated) {
        e0.setIsolated(false);
        e1.setIsolated(false);
      }
      //intersectionFound = true;
      _numIntersections++;
      // if the segments are adjacent they have at least one trivial intersection,
      // the shared endpoint.  Don't bother adding it if it is the
      // only intersection.
      if (! isTrivialIntersection(e0, segIndex0, e1, segIndex1)) {
        _hasIntersection = true;
        if (_includeProper || ! _li.isProper() ) {
//Debug.println(li);
          e0.addIntersections(_li, segIndex0, 0);
          e1.addIntersections(_li, segIndex1, 1);
        }
        if (_li.isProper()) {
          _properIntersectionPoint = _li.getIntersection(0).copy();
          _hasProper = true;
          if (_isDoneWhenProperInt) {
            _isDone = true;
          }
          if (! isBoundaryPoint(_li, bdyNodes))
            _hasProperInterior = true;
        }
        //if (li.isCollinear())
        //hasCollinear = true;
      }
    }
  }

   bool isBoundaryPoint(LineIntersector li, List<List> bdyNodes)
  {
  if (bdyNodes == null) return false;
  if (isBoundaryPointInternal(li, bdyNodes[0])) return true;
  if (isBoundaryPointInternal(li, bdyNodes[1])) return true;
  return false;
  }

   bool isBoundaryPointInternal(LineIntersector li, List<Node> bdyNodes)
  {

    for(int i = 0; i < bdyNodes.length; i++){
      Node node = bdyNodes[i];
      Coordinate pt = node.getCoordinate();
      if (li.isIntersection(pt)) return true;
    }
    return false;
  }

}



/**
 * MonotoneChains are a way of partitioning the segments of an edge to
 * allow for fast searching of intersections.
 * They have the following properties:
 * <ol>
 * <li>the segments within a monotone chain will never intersect each other
 * <li>the envelope of any contiguous subset of the segments in a monotone chain
 * is simply the envelope of the endpoints of the subset.
 * </ol>
 * Property 1 means that there is no need to test pairs of segments from within
 * the same monotone chain for intersection.
 * Property 2 allows
 * binary search to be used to find the intersection points of two monotone chains.
 * For many types of real-world data, these properties eliminate a large number of
 * segment comparisons, producing substantial speed gains.
 * @version 1.7
 */
class MonotoneChainEdge {

  Edge e;
  List<Coordinate> pts; // cache a reference to the coord array, for efficiency
  // the lists of start/end indexes of the monotone chains.
  // Includes the end point of the edge as a sentinel
  List<int> startIndex;

  MonotoneChainEdge(Edge e) {
    this.e = e;
    pts = e.getCoordinates();
    MonotoneChainIndexer mcb = new MonotoneChainIndexer();
    startIndex = mcb.getChainStartIndices(pts);
  }

  List<Coordinate> getCoordinates() { return pts; }
  List<int> getStartIndexes() { return startIndex; }

  double getMinX(int chainIndex)
  {
    double x1 = pts[startIndex[chainIndex]].x;
    double x2 = pts[startIndex[chainIndex + 1]].x;
    return x1 < x2 ? x1 : x2;
  }
  double getMaxX(int chainIndex)
  {
    double x1 = pts[startIndex[chainIndex]].x;
    double x2 = pts[startIndex[chainIndex + 1]].x;
    return x1 > x2 ? x1 : x2;
  }

  void computeIntersects(MonotoneChainEdge mce, SegmentIntersector si)
  {
    for (int i = 0; i < startIndex.length - 1; i++) {
      for (int j = 0; j < mce.startIndex.length - 1; j++) {
        computeIntersectsForChain(  i,
            mce,  j,
            si );
      }
    }
  }
  void computeIntersectsForChain(
      int chainIndex0,
      MonotoneChainEdge mce,
      int chainIndex1,
      SegmentIntersector si)
  {
    computeIntersectsForChain(startIndex[chainIndex0], startIndex[chainIndex0 + 1],
        mce,
        mce.startIndex[chainIndex1], mce.startIndex[chainIndex1 + 1],
        si );
  }

  private void computeIntersectsForChain(
      int start0, int end0,
      MonotoneChainEdge mce,
      int start1, int end1,
      SegmentIntersector ei)
  {
//Debug.println("computeIntersectsForChain:" + p00 + p01 + p10 + p11);

    // terminating condition for the recursion
    if (end0 - start0 == 1 && end1 - start1 == 1) {
      ei.addIntersections(e, start0, mce.e, start1);
      return;
    }
    // nothing to do if the envelopes of these chains don't overlap
    if (! overlaps(start0, end0, mce, start1, end1)) return;

    // the chains overlap, so split each in half and iterate  (binary search)
    int mid0 = (start0 + end0) / 2;
    int mid1 = (start1 + end1) / 2;

    // Assert: mid != start or end (since we checked above for end - start <= 1)
    // check terminating conditions before recursing
    if (start0 < mid0) {
      if (start1 < mid1) computeIntersectsForChain(start0, mid0, mce, start1,  mid1, ei);
      if (mid1 < end1)   computeIntersectsForChain(start0, mid0, mce, mid1,    end1, ei);
    }
    if (mid0 < end0) {
      if (start1 < mid1) computeIntersectsForChain(mid0,   end0, mce, start1,  mid1, ei);
      if (mid1 < end1)   computeIntersectsForChain(mid0,   end0, mce, mid1,    end1, ei);
    }
  }

  /**
   * Tests whether the envelopes of two chain sections overlap (intersect).
   *
   * @param start0
   * @param end0
   * @param mce
   * @param start1
   * @param end1
   * @return true if the section envelopes overlap
   */
  private boolean overlaps(
      int start0, int end0,
      MonotoneChainEdge mce,
      int start1, int end1)
  {
    return Envelope.intersects(pts[start0], pts[end0], mce.pts[start1], mce.pts[end1]);
  }

}
