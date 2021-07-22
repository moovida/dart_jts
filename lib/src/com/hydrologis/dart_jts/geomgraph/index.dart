part of dart_jts;

class IntervalRTreeBranchNode extends IntervalRTreeNode {
  IntervalRTreeNode node1;
  IntervalRTreeNode node2;

  IntervalRTreeBranchNode(this.node1, this.node2) {
    buildExtent(node1, node2);
  }

  void buildExtent(IntervalRTreeNode n1, IntervalRTreeNode n2) {
    min = math.min(n1.min, n2.min);
    max = math.max(n1.max, n2.max);
  }

  void query(double queryMin, double queryMax, ItemVisitor visitor) {
    if (!intersects(queryMin, queryMax)) {
//			System.out.println("Does NOT Overlap branch: " + this);
      return;
    }
//		System.out.println("Overlaps branch: " + this);
    if (node1 != null) node1.query(queryMin, queryMax, visitor);
    if (node2 != null) node2.query(queryMin, queryMax, visitor);
  }
}

class IntervalRTreeLeafNode extends IntervalRTreeNode {
  Object item;

  IntervalRTreeLeafNode(double min, double max, this.item) {
    this.min = min;
    this.max = max;
  }

  void query(double queryMin, double queryMax, ItemVisitor visitor) {
    if (!intersects(queryMin, queryMax)) return;

    visitor.visitItem(item);
  }
}

/**
 * A static index on a set of 1-dimensional intervals,
 * using an R-Tree packed based on the order of the interval midpoints.
 * It supports range searching,
 * where the range is an interval of the real line (which may be a single point).
 * A common use is to index 1-dimensional intervals which
 * are the projection of 2-D objects onto an axis of the coordinate system.
 * <p>
 * This index structure is <i>static</i>
 * - items cannot be added or removed once the first query has been made.
 * The advantage of this characteristic is that the index performance
 * can be optimized based on a fixed set of items.
 *
 * @author Martin Davis
 */
class SortedPackedIntervalRTree {
  List leaves = [];

  /**
   * If root is null that indicates
   * that the tree has not yet been built,
   * OR nothing has been added to the tree.
   * In both cases, the tree is still open for insertions.
   */
  IntervalRTreeNode? root = null;

  SortedPackedIntervalRTree() {}

  /**
   * Adds an item to the index which is associated with the given interval
   *
   * @param min the lower bound of the item interval
   * @param max the upper bound of the item interval
   * @param item the item to insert
   *
   * @throws IllegalStateException if the index has already been queried
   */
  void insert(double min, double max, Object item) {
    if (root != null)
      throw new StateError("Index cannot be added to once it has been queried");
    leaves.add(new IntervalRTreeLeafNode(min, max, item));
  }

  void init() {
    // already built
    if (root != null) return;

    /**
     * if leaves is empty then nothing has been inserted.
     * In this case it is safe to leave the tree in an open state
     */
    if (leaves.isEmpty) return;

    buildRoot();
  }

  void buildRoot() // TODO check how to make this methos synchronized
  {
    if (root != null) return;
    root = buildTree();
  }

  IntervalRTreeNode buildTree() {
    // sort the leaf nodes
    leaves.sort(nodeComparator);

    // now group nodes into blocks of two and build tree up recursively
    List src = leaves;
    List? temp = null;
    List dest = [];

    while (true) {
      buildLevel(src, dest);
      if (dest.length == 1) return dest[0];

      temp = src;
      src = dest;
      dest = temp;
    }
  }

  int level = 0;

  void buildLevel(List src, List dest) {
    level++;
    dest.clear();
    for (int i = 0; i < src.length; i += 2) {
      IntervalRTreeNode n1 = src[i] as IntervalRTreeNode;
      IntervalRTreeNode? n2 =
          (i + 1 < src.length) ? src[i] as IntervalRTreeNode : null;
      if (n2 == null) {
        dest.add(n1);
      } else {
        IntervalRTreeNode node =
            new IntervalRTreeBranchNode(src[i], src[i + 1]);
//        printNode(node);
//				System.out.println(node);
        dest.add(node);
      }
    }
  }

  void printNode(IntervalRTreeNode node) {
    print(WKTWriter.toLineStringFromCoords(
        new Coordinate(node.min, level.toDouble()),
        new Coordinate(node.max, level.toDouble())));
  }

  /**
   * Search for intervals in the index which intersect the given closed interval
   * and apply the visitor to them.
   *
   * @param min the lower bound of the query interval
   * @param max the upper bound of the query interval
   * @param visitor the visitor to pass any matched items to
   */
  void query(double min, double max, ItemVisitor visitor) {
    init();

    // if root is null tree must be empty
    if (root == null) return;

    root!.query(min, max, visitor);
  }
}

abstract class IntervalRTreeNode {
  double min = double.infinity;
  double max = double.negativeInfinity;

  double getMin() {
    return min;
  }

  double getMax() {
    return max;
  }

  void query(double queryMin, double queryMax, ItemVisitor visitor);

  bool intersects(double queryMin, double queryMax) {
    if (min > queryMax || max < queryMin) return false;
    return true;
  }

  String toString() {
    return WKTWriter.toLineStringFromCoords(
        new Coordinate(min, 0), new Coordinate(max, 0));
  }
}

Comparator nodeComparator = (o1, o2) {
  IntervalRTreeNode n1 = o1 as IntervalRTreeNode;
  IntervalRTreeNode n2 = o2 as IntervalRTreeNode;
  double mid1 = (n1.min + n1.max) / 2;
  double mid2 = (n2.min + n2.max) / 2;
  if (mid1 < mid2) return -1;
  if (mid1 > mid2) return 1;
  return 0;
};

/**
 * A visitor for items in a {@link SpatialIndex}.
 *
 * @version 1.7
 */

abstract class ItemVisitor {
  /**
   * Visits an item in the index.
   *
   * @param item the index item to be visited
   */
  void visitItem(Object item);
}

/**
 * Builds an array of all visited items.
 *
 * @version 1.7
 */
class ArrayListVisitor implements ItemVisitor {
  List items = [];

  /**
   * Creates a new instance.
   */
  ArrayListVisitor() {}

  /**
   * Visits an item.
   *
   * @param item the item to visit
   */
  void visitItem(Object item) {
    items.add(item);
  }

  /**
   * Gets the array of visited items.
   *
   * @return the array of items
   */
  List getItems() {
    return items;
  }
}

/**
 * @version 1.7
 */
class SweepLineEvent implements Comparable {
  static final int INSERT = 1;
  static final int DELETE = 2;

  Object? label; // used for red-blue intersection detection
  double xValue;
  int eventType = 0;
  SweepLineEvent? insertEvent = null; // null if this is an INSERT event
  int deleteEventIndex = 0;
  Object? obj;

  /**
   * Creates an INSERT event.
   *
   * @param label the edge set label for this object
   * @param x the event location
   * @param obj the object being inserted
   */
  SweepLineEvent(this.label, this.xValue, this.obj) {
    this.eventType = INSERT;
  }

  /**
   * Creates a DELETE event.
   *
   * @param x the event location
   * @param insertEvent the corresponding INSERT event
   */
  SweepLineEvent.withEvent(this.xValue, this.insertEvent) {
    eventType = DELETE;
  }

  bool isInsert() {
    return eventType == INSERT;
  }

  bool isDelete() {
    return eventType == DELETE;
  }

  SweepLineEvent? getInsertEvent() {
    return insertEvent;
  }

  int getDeleteEventIndex() {
    return deleteEventIndex;
  }

  void setDeleteEventIndex(int deleteEventIndex) {
    this.deleteEventIndex = deleteEventIndex;
  }

  Object? getObject() {
    return obj;
  }

  bool isSameLabel(SweepLineEvent ev) {
    // no label set indicates single group
    if (label == null) return false;
    return label == ev.label;
  }

  /**
   * Events are ordered first by their x-value, and then by their eventType.
   * Insert events are sorted before Delete events, so that
   * items whose Insert and Delete events occur at the same x-value will be
   * correctly handled.
   */
  int compareTo(dynamic o) {
    SweepLineEvent pe = o;
    if (xValue < pe.xValue) return -1;
    if (xValue > pe.xValue) return 1;
    if (eventType < pe.eventType) return -1;
    if (eventType > pe.eventType) return 1;
    return 0;
  }
}

/**
 * An EdgeSetIntersector computes all the intersections between the
 * edges in the set.  It adds the computed intersections to each edge
 * they are found on.  It may be used in two scenarios:
 * <ul>
 * <li>determining the internal intersections between a single set of edges
 * <li>determining the mutual intersections between two different sets of edges
 * </ul>
 * It uses a {@link SegmentIntersector} to compute the intersections between
 * segments and to record statistics about what kinds of intersections were found.
 *
 * @version 1.7
 */
abstract class EdgeSetIntersector {
  EdgeSetIntersector() {}

  /**
   * Computes all self-intersections between edges in a set of edges,
   * allowing client to choose whether self-intersections are computed.
   *
   * @param edges a list of edges to test for intersections
   * @param si the SegmentIntersector to use
   * @param testAllSegments true if self-intersections are to be tested as well
   */
  void computeIntersections(
      List edges, SegmentIntersector si, bool testAllSegments);

  /**
   * Computes all mutual intersections between two sets of edges.
   */
  void computeIntersections3(List edges0, List edges1, SegmentIntersector si);
}

/**
 * Finds all intersections in one or two sets of edges,
 * using an x-axis sweepline algorithm in conjunction with Monotone Chains.
 * While still O(n^2) in the worst case, this algorithm
 * drastically improves the average-case time.
 * The use of MonotoneChains as the items in the index
 * seems to offer an improvement in performance over a sweep-line alone.
 *
 * @version 1.7
 */
class SimpleMCSweepLineIntersector extends EdgeSetIntersector {
  List events = [];

  // statistics information
  int nOverlaps = 0;

  /**
   * A SimpleMCSweepLineIntersector creates monotone chains from the edges
   * and compares them using a simple sweep-line along the x-axis.
   */
  SimpleMCSweepLineIntersector() {}

  void computeIntersections(
      List edges, SegmentIntersector si, bool testAllSegments) {
    if (testAllSegments)
      addEdgesWithSet(edges, null);
    else
      addEdges(edges);
    computeIntersections1(si);
  }

  void computeIntersections3(List edges0, List edges1, SegmentIntersector si) {
    addEdgesWithSet(edges0, edges0);
    addEdgesWithSet(edges1, edges1);
    computeIntersections1(si);
  }

  void addEdges(List edges) {
    for (Iterator i = edges.iterator; i.moveNext();) {
      Edge edge = i.current;
      // edge is its own group
      addEdge(edge, edge);
    }
  }

  void addEdgesWithSet(List edges, Object? edgeSet) {
    for (Iterator i = edges.iterator; i.moveNext();) {
      Edge edge = i.current;
      addEdge(edge, edgeSet);
    }
  }

  void addEdge(Edge edge, Object? edgeSet) {
    MonotoneChainEdge mce = edge.getMonotoneChainEdge();
    List<int> startIndex = mce.getStartIndexes();
    for (int i = 0; i < startIndex.length - 1; i++) {
      MonotoneChain mc = new MonotoneChain(mce, i);
      SweepLineEvent insertEvent =
          new SweepLineEvent(edgeSet, mce.getMinX(i), mc);
      events.add(insertEvent);
      events.add(new SweepLineEvent.withEvent(mce.getMaxX(i), insertEvent));
    }
  }

  /**
   * Because Delete Events have a link to their corresponding Insert event,
   * it is possible to compute exactly the range of events which must be
   * compared to a given Insert event object.
   */
  void prepareEvents() {
    events.sort();
    // set DELETE event indexes
    for (int i = 0; i < events.length; i++) {
      SweepLineEvent ev = events[i];
      if (ev.isDelete()) {
        ev.getInsertEvent()!.setDeleteEventIndex(i);
      }
    }
  }

  void computeIntersections1(SegmentIntersector si) {
    nOverlaps = 0;
    prepareEvents();

    for (int i = 0; i < events.length; i++) {
      SweepLineEvent ev = events[i];
      if (ev.isInsert()) {
        processOverlaps(i, ev.getDeleteEventIndex(), ev, si);
      }
      if (si.isDone()) {
        break;
      }
    }
  }

  void processOverlaps(
      int start, int end, SweepLineEvent ev0, SegmentIntersector si) {
    MonotoneChain mc0 = ev0.getObject() as MonotoneChain;
    /**
     * Since we might need to test for self-intersections,
     * include current INSERT event object in list of event objects to test.
     * Last index can be skipped, because it must be a Delete event.
     */
    for (int i = start; i < end; i++) {
      SweepLineEvent ev1 = events[i];
      if (ev1.isInsert()) {
        MonotoneChain mc1 = ev1.getObject() as MonotoneChain;
        // don't compare edges in same group, if labels are present
        if (!ev0.isSameLabel(ev1)) {
          mc0.computeIntersections(mc1, si);
          nOverlaps++;
        }
      }
    }
  }
}

/**
 * Computes the intersection of line segments,
 * and adds the intersection to the edges containing the segments.
 *
 * @version 1.7
 */
class SegmentIntersector {
  static bool isAdjacentSegments(int i1, int i2) {
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
  Coordinate? _properIntersectionPoint = null;

  late LineIntersector _li;
  bool _includeProper = false;
  bool _recordIsolated = false;
  bool _isSelfIntersection = false;

  // bool intersectionFound;
  int _numIntersections = 0;

  // testing only
  int _numTests = 0;

  List<List>? bdyNodes;
  bool _isDone = false;
  bool _isDoneWhenProperInt = false;

  SegmentIntersector(
      LineIntersector li, bool includeProper, bool recordIsolated) {
    this._li = li;
    this._includeProper = includeProper;
    this._recordIsolated = recordIsolated;
  }

  void setBoundaryNodes(List bdyNodes0, List bdyNodes1) {
    bdyNodes = []; //..length = (2);
    bdyNodes!.add(bdyNodes0);
    bdyNodes!.add(bdyNodes1);
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
  Coordinate? getProperIntersectionPoint() {
    return _properIntersectionPoint;
  }

  bool hasIntersection() {
    return _hasIntersection;
  }

  /**
   * A proper intersection is an intersection which is interior to at least two
   * line segments.  Note that a proper intersection is not necessarily
   * in the interior of the entire Geometry, since another edge may have
   * an endpoint equal to the intersection, which according to SFS semantics
   * can result in the point being on the Boundary of the Geometry.
   */
  bool hasProperIntersection() {
    return _hasProper;
  }

  /**
   * A proper interior intersection is a proper intersection which is <b>not</b>
   * contained in the set of boundary nodes set for this SegmentIntersector.
   */
  bool hasProperInteriorIntersection() {
    return _hasProperInterior;
  }

  /**
   * A trivial intersection is an apparent self-intersection which in fact
   * is simply the point shared by adjacent line segments.
   * Note that closed edges require a special check for the point shared by the beginning
   * and end segments.
   */
  bool isTrivialIntersection(Edge e0, int segIndex0, Edge e1, int segIndex1) {
    if (e0 == e1) {
      if (_li.getIntersectionNum() == 1) {
        if (isAdjacentSegments(segIndex0, segIndex1)) return true;
        if (e0.isClosed()) {
          int maxSegIndex = e0.getNumPoints() - 1;
          if ((segIndex0 == 0 && segIndex1 == maxSegIndex) ||
              (segIndex1 == 0 && segIndex0 == maxSegIndex)) {
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
  void addIntersections(Edge e0, int segIndex0, Edge e1, int segIndex1) {
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
      if (!isTrivialIntersection(e0, segIndex0, e1, segIndex1)) {
        _hasIntersection = true;
        if (_includeProper || !_li.isProper()) {
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
          if (!isBoundaryPoint(_li, bdyNodes)) _hasProperInterior = true;
        }
        //if (li.isCollinear())
        //hasCollinear = true;
      }
    }
  }

  bool isBoundaryPoint(LineIntersector li, List<List>? bdyNodes) {
    if (bdyNodes == null) return false;
    if (isBoundaryPointInternal(li, bdyNodes[0].cast<Node>())) return true;
    if (isBoundaryPointInternal(li, bdyNodes[1].cast<Node>())) return true;
    return false;
  }

  bool isBoundaryPointInternal(LineIntersector li, List<Node> bdyNodes) {
    for (int i = 0; i < bdyNodes.length; i++) {
      Node node = bdyNodes[i];
      Coordinate pt = node.getCoordinate();
      if (li.isIntersection(pt)) return true;
    }
    return false;
  }
}

/**
 * @version 1.7
 */
class MonotoneChain {
  MonotoneChainEdge mce;
  int chainIndex;

  MonotoneChain(this.mce, this.chainIndex);

  void computeIntersections(MonotoneChain mc, SegmentIntersector si) {
    this.mce.computeIntersectsForChain(chainIndex, mc.mce, mc.chainIndex, si);
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
  late List<Coordinate>
      pts; // cache a reference to the coord array, for efficiency
  // the lists of start/end indexes of the monotone chains.
  // Includes the end point of the edge as a sentinel
  late List<int> startIndex;

  MonotoneChainEdge(this.e) {
    pts = e.getCoordinates();
    MonotoneChainIndexer mcb = new MonotoneChainIndexer();
    startIndex = mcb.getChainStartIndices(pts);
  }

  List<Coordinate> getCoordinates() {
    return pts;
  }

  List<int> getStartIndexes() {
    return startIndex;
  }

  double getMinX(int chainIndex) {
    double x1 = pts[startIndex[chainIndex]].x;
    double x2 = pts[startIndex[chainIndex + 1]].x;
    return x1 < x2 ? x1 : x2;
  }

  double getMaxX(int chainIndex) {
    double x1 = pts[startIndex[chainIndex]].x;
    double x2 = pts[startIndex[chainIndex + 1]].x;
    return x1 > x2 ? x1 : x2;
  }

  void computeIntersects(MonotoneChainEdge mce, SegmentIntersector si) {
    for (int i = 0; i < startIndex.length - 1; i++) {
      for (int j = 0; j < mce.startIndex.length - 1; j++) {
        computeIntersectsForChain(i, mce, j, si);
      }
    }
  }

  void computeIntersectsForChain(int chainIndex0, MonotoneChainEdge mce,
      int chainIndex1, SegmentIntersector si) {
    computeIntersectsForChain6(
        startIndex[chainIndex0],
        startIndex[chainIndex0 + 1],
        mce,
        mce.startIndex[chainIndex1],
        mce.startIndex[chainIndex1 + 1],
        si);
  }

  void computeIntersectsForChain6(int start0, int end0, MonotoneChainEdge mce,
      int start1, int end1, SegmentIntersector ei) {
//Debug.println("computeIntersectsForChain:" + p00 + p01 + p10 + p11);

    // terminating condition for the recursion
    if (end0 - start0 == 1 && end1 - start1 == 1) {
      ei.addIntersections(e, start0, mce.e, start1);
      return;
    }
    // nothing to do if the envelopes of these chains don't overlap
    if (!overlaps(start0, end0, mce, start1, end1)) return;

    // the chains overlap, so split each in half and iterate  (binary search)
    int mid0 = (start0 + end0) ~/ 2;
    int mid1 = (start1 + end1) ~/ 2;

    // Assert: mid != start or end (since we checked above for end - start <= 1)
    // check terminating conditions before recursing
    if (start0 < mid0) {
      if (start1 < mid1)
        computeIntersectsForChain6(start0, mid0, mce, start1, mid1, ei);
      if (mid1 < end1)
        computeIntersectsForChain6(start0, mid0, mce, mid1, end1, ei);
    }
    if (mid0 < end0) {
      if (start1 < mid1)
        computeIntersectsForChain6(mid0, end0, mce, start1, mid1, ei);
      if (mid1 < end1)
        computeIntersectsForChain6(mid0, end0, mce, mid1, end1, ei);
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
  bool overlaps(
      int start0, int end0, MonotoneChainEdge mce, int start1, int end1) {
    return Envelope.intersectsEnvelopeCoords(
        pts[start0], pts[end0], mce.pts[start1], mce.pts[end1]);
  }
}

/**
 * MonotoneChains are a way of partitioning the segments of an edge to
 * allow for fast searching of intersections.
 * Specifically, a sequence of contiguous line segments
 * is a monotone chain iff all the vectors defined by the oriented segments
 * lies in the same quadrant.
 * <p>
 * Monotone Chains have the following useful properties:
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
 * <p>
 * Note that due to the efficient intersection test, there is no need to limit the size
 * of chains to obtain fast performance.
 *
 * @version 1.7
 */
class MonotoneChainIndexer {
  static List<int> toIntArray(List list) {
    List<int> array = []; //..length = (list.length);
    for (int i = 0; i < list.length; i++) {
      array.add(list[i] as int);
    }
    return array;
  }

  MonotoneChainIndexer() {}

  List<int> getChainStartIndices(List<Coordinate> pts) {
    // find the startpoint (and endpoints) of all monotone chains in this edge
    int start = 0;
    List<int> startIndexList = []; //List(pts.length ~/ 2, );
    // use heuristic to size initial array
    //startIndexList.ensureCapacity(pts.length / 4);
    startIndexList.add(start);
    do {
      int last = findChainEnd(pts, start);
      startIndexList.add(last);
      start = last;
    } while (start < pts.length - 1);
    // copy list to an array of ints, for efficiency
    return startIndexList;
  }

  List<int> OLDgetChainStartIndices(List<Coordinate> pts) {
    // find the startpoint (and endpoints) of all monotone chains in this edge
    int start = 0;
    List startIndexList = [];
    startIndexList.add(start);
    do {
      int last = findChainEnd(pts, start);
      startIndexList.add(last);
      start = last;
    } while (start < pts.length - 1);
    // copy list to an array of ints, for efficiency
    List<int> startIndex = toIntArray(startIndexList);
    return startIndex;
  }

  /**
   * @return the index of the last point in the monotone chain
   */
  int findChainEnd(List<Coordinate> pts, int start) {
    // determine quadrant for chain
    int chainQuad = Quadrant.quadrantFromCoords(pts[start], pts[start + 1]);
    int last = start + 1;
    while (last < pts.length) {
      //if (last - start > 100) break;
      // compute quadrant for next possible segment in chain
      int quad = Quadrant.quadrantFromCoords(pts[last - 1], pts[last]);
      if (quad != chainQuad) break;
      last++;
    }
    return last - 1;
  }
}
