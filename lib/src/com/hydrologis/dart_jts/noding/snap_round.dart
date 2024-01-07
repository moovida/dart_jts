part of dart_jts;

/**
 * Uses Snap Rounding to compute a rounded,
 * fully noded arrangement from a set of {@link SegmentString}s,
 * in a performant way, and avoiding unnecessary noding.
 * <p>
 * Implements the Snap Rounding technique described in
 * the papers by Hobby, Guibas &amp; Marimont, and Goodrich et al.
 * Snap Rounding enforces that all output vertices lie on a uniform grid,
 * which is determined by the provided {@link PrecisionModel}.
 * <p>
 * Input vertices do not have to be rounded to the grid beforehand;
 * this is done during the snap-rounding process.
 * In fact, rounding cannot be done a priori,
 * since rounding vertices by themselves can distort the rounded topology
 * of the arrangement (i.e. by moving segments away from hot pixels
 * that would otherwise intersect them, or by moving vertices
 * across segments).
 * <p>
 * To minimize the number of introduced nodes,
 * the Snap-Rounding Noder avoids creating nodes
 * at edge vertices if there is no intersection or snap at that location.
 * However, if two different input edges contain identical segments,
 * each of the segment vertices will be noded.
 * This still provides fully-noded output.
 * This is the same behaviour provided by other noders,
 * such as {@link MCIndexNoder} and {@link org.locationtech.jts.noding.snap.SnappingNoder}.
 *
 * @version 1.7
 */
class SnapRoundingNoder implements Noder {
  /**
   * The division factor used to determine
   * nearness distance tolerance for intersection detection.
   */
  static final int NEARNESS_FACTOR = 100;

  PrecisionModel pm;
  late HotPixelIndex pixelIndex;

  List<NodedSegmentString> snappedResult = List.empty(growable: true);

  SnapRoundingNoder(this.pm) {
    this.pm = pm;
    pixelIndex = new HotPixelIndex(pm);
  }

  /**
   * @return a List of NodedSegmentStrings representing the substrings
   *
   */
  List getNodedSubstrings() {
    return NodedSegmentString.getNodedSubstrings(snappedResult);
  }

  /**
   * Computes the nodes in the snap-rounding line arrangement.
   * The nodes are added to the {@link NodedSegmentString}s provided as the input.
   *
   * @param inputSegmentStrings a List of NodedSegmentStrings
   */
  void computeNodes(List inputSegmentStrings) {
    snappedResult = _snapRound(inputSegmentStrings as List<NodedSegmentString>);
  }

  List<NodedSegmentString> _snapRound(List<NodedSegmentString> segStrings) {
    /**
     * Determine hot pixels for intersections and vertices.
     * This is done BEFORE the input lines are rounded,
     * to avoid distorting the line arrangement
     * (rounding can cause vertices to move across edges).
     */
    _addIntersectionPixels(segStrings);
    _addVertexPixels(segStrings);

    List<NodedSegmentString> snapped = _computeSnaps(segStrings);
    return snapped;
  }

  /**
   * Detects interior intersections in the collection of {@link SegmentString}s,
   * and adds nodes for them to the segment strings.
   * Also creates HotPixel nodes for the intersection points.
   *
   * @param segStrings the input NodedSegmentStrings
   */
  void _addIntersectionPixels(List<NodedSegmentString> segStrings) {
    /**
     * nearness tolerance is a small fraction of the grid size.
     */
    double snapGridSize = 1.0 / pm.getScale();
    double nearnessTol = snapGridSize / NEARNESS_FACTOR;

    SnapRoundingIntersectionAdder intAdder =
        SnapRoundingIntersectionAdder(nearnessTol);
    MCIndexNoder noder = MCIndexNoder.withTolerance(intAdder,nearnessTol);
    noder.computeNodes(segStrings);
    List<Coordinate> intPts = intAdder.getIntersections();
    pixelIndex.addNodes(intPts);
  }

  /**
   * Creates HotPixels for each vertex in the input segStrings.
   * The HotPixels are not marked as nodes, since they will
   * only be nodes in the final line arrangement
   * if they interact with other segments (or they are already
   * created as intersection nodes).
   *
   * @param segStrings the input NodedSegmentStrings
   */
  void _addVertexPixels(List<NodedSegmentString> segStrings) {
    for (SegmentString nss in segStrings) {
      List<Coordinate> pts = nss.getCoordinates();
      pixelIndex.addNodes(pts);
    }
  }

  Coordinate _roundCoordinate(Coordinate pt) {
    Coordinate p2 = pt.copy();
    pm.makeCoordinatePrecise(p2);
    return p2;
  }

  /**
   * Gets a list of the rounded coordinates.
   * Duplicate (collapsed) coordinates are removed.
   *
   * @param pts the coordinates to round
   * @return array of rounded coordinates
   */
  List<Coordinate> _roundCoordList(List<Coordinate> pts) {
    CoordinateList roundPts = new CoordinateList();

    for (int i = 0; i < pts.length; i++) {
      roundPts.addCoord(_roundCoordinate(pts[i]), false);
    }
    return roundPts.toCoordinateArray(true);
  }

  /**
   * Computes new segment strings which are rounded and contain
   * intersections added as a result of snapping segments to snap points (hot pixels).
   *
   * @param segStrings segments to snap
   * @return the snapped segment strings
   */
  List<NodedSegmentString> _computeSnaps(List<NodedSegmentString> segStrings) {
    List<NodedSegmentString> snapped = List.empty(growable: true);
    for (NodedSegmentString ss in segStrings) {
      NodedSegmentString? snappedSS = _computeSegmentSnaps(ss);
      if(snappedSS != null) {
        snapped.add(snappedSS);
      }
    }
    /**
     * Some intersection hot pixels may have been marked as nodes in the previous
     * loop, so add nodes for them.
     */
    for (NodedSegmentString ss in snapped) {
      _addVertexNodeSnaps(ss);
    }
    return snapped;
  }

  /**
   * Add snapped vertices to a segment string.
   * If the segment string collapses completely due to rounding,
   * null is returned.
   *
   * @param ss the segment string to snap
   * @return the snapped segment string, or null if it collapses completely
   */
  NodedSegmentString? _computeSegmentSnaps(NodedSegmentString ss) {
    //List<Coordinate> pts = ss.getCoordinates();
    /**
   * Get edge coordinates, including added intersection nodes.
   * The coordinates are now rounded to the grid,
   * in preparation for snapping to the Hot Pixels
   */
    List<Coordinate> pts = ss.getNodedCoordinates();
    List<Coordinate> ptsRound = _roundCoordList(pts);

    // if complete collapse this edge can be eliminated
    if (ptsRound.length <= 1) {
      return null;
    }

    // Create new nodedSS to allow adding any hot pixel nodes
    NodedSegmentString snapSS = new NodedSegmentString(ptsRound, ss.getData());

    int snapSSindex = 0;
    for (int i = 0; i < pts.length - 1; i++) {
      Coordinate currSnap = snapSS.getCoordinate(snapSSindex);

      /**
   * If the segment has collapsed completely, skip it
   */
      Coordinate p1 = pts[i + 1];
      Coordinate p1Round = _roundCoordinate(p1);
      if (p1Round.equals2D(currSnap)) continue;

      Coordinate p0 = pts[i];

      /**
   * Add any Hot Pixel intersections with *original* segment to rounded segment.
   * (It is important to check original segment because rounding can
   * move it enough to intersect other hot pixels not intersecting original segment)
   */
      _snapSegment(p0, p1, snapSS, snapSSindex);
      snapSSindex++;
    }
    return snapSS;
  }

  /**
   * Snaps a segment in a segmentString to HotPixels that it intersects.
   *
   * @param p0 the segment start coordinate
   * @param p1 the segment end coordinate
   * @param ss the segment string to add intersections to
   * @param segIndex the index of the segment
   */
  void _snapSegment(
      Coordinate p0, Coordinate p1, NodedSegmentString ss, int segIndex) {
    pixelIndex.query(
        p0, p1, new KdNodeVisitorSnapSegment(ss, p0, p1, segIndex));
  }

  /**
   * Add nodes for any vertices in hot pixels that were
   * added as nodes during segment noding.
   *
   * @param ss a noded segment string
   */
  void _addVertexNodeSnaps(NodedSegmentString ss) {
    List<Coordinate> pts = ss.getCoordinates();
    for (int i = 1; i < pts.length - 1; i++) {
      Coordinate p0 = pts[i];
      _snapVertexNode(p0, ss, i);
    }
  }

  void _snapVertexNode(Coordinate p0, NodedSegmentString ss, int segIndex) {
    pixelIndex.query(p0, p0, KdNodeVisitorVertex(ss, p0, segIndex));
  }
}

class KdNodeVisitorVertex extends KdNodeVisitor {
  NodedSegmentString ss;
  Coordinate p0;
  int segIndex;

  KdNodeVisitorVertex(this.ss, this.p0, this.segIndex);

  @override
  void visit(KdNode node) {
    HotPixel hp = node.getData() as HotPixel;
    if (hp.isNode && hp.getCoordinate().equals2D(p0)) {
      ss.addIntersection(p0, segIndex);
    }
  }
}

class KdNodeVisitorSnapSegment extends KdNodeVisitor {
  NodedSegmentString ss;
  Coordinate p0;
  Coordinate p1;
  int segIndex;
  KdNodeVisitorSnapSegment(this.ss, this.p0, this.p1, this.segIndex);
  @override
  void visit(KdNode node) {
    HotPixel? hp = node.getData() as HotPixel?;

    /**
     * If the hot pixel is not a node, and it contains one of the segment vertices,
     * then that vertex is the source for the hot pixel.
     * To avoid over-noding a node is not added at this point.
     * The hot pixel may be subsequently marked as a node,
     * in which case the intersection will be added during the final vertex noding phase.
     */
    if (!hp!.isNode) {
      if (hp.intersects(p0) || hp.intersects(p1)) return;
    }
    /**
     * Add a node if the segment intersects the pixel.
     * Mark the HotPixel as a node (since it may not have been one before).
     * This ensures the vertex for it is added as a node during the final vertex noding phase.
     */
    if (hp.intersectsSegment(p0, p1)) {
      //System.out.println("Added intersection: " + hp.getCoordinate());
      ss.addIntersection(hp.getCoordinate(), segIndex);
      hp.setToNode();
    }
  }
}

/**
 * An index which creates unique {@link HotPixel}s for provided points,
 * and performs range queries on them.
 * The points passed to the index do not needed to be
 * rounded to the specified scale factor; this is done internally
 * when creating the HotPixels for them.
 *
 * @author mdavis
 *
 */
class HotPixelIndex {
  PrecisionModel precModel;
  late double scaleFactor;

  /**
   * Use a kd-tree to index the pixel centers for optimum performance.
   * Since HotPixels have an extent, range queries to the
   * index must enlarge the query range by a suitable value
   * (using the pixel width is safest).
   */
  KdTree index = new KdTree.withZeroTolerance();

  HotPixelIndex(this.precModel) {
    scaleFactor = precModel.getScale();
  }
  /**
   * Adds a list of points as non-node pixels.
   *
   * @param pts the points to add
   */
  void add(List<Coordinate> pts) {
    /**
   * Shuffle the points before adding.
   * This avoids having long monontic runs of points
   * causing an unbalanced KD-tree, which would create
   * performance and robustness issues.
   */
    Iterator<Coordinate> it = CoordinateShuffler(pts);
    while (it.moveNext()) {
      addCoordinate(it.current);
    }
  }

  /**
   * Adds a list of points as node pixels.
   *
   * @param pts the points to add
   */
  void addNodes(List<Coordinate> pts) {
    /**
   * Node points are not shuffled, since they are
   * added after the vertex points, and hence the KD-tree should
   * be reasonably balanced already.
   */
    for (Coordinate pt in pts) {
      HotPixel hp = addCoordinate(pt);
      hp.setToNode();
    }
  }

  /**
   * Adds a point as a Hot Pixel.
   * If the point has been added already, it is marked as a node.
   *
   * @param p the point to add
   * @return the HotPixel for the point
   */
  HotPixel addCoordinate(Coordinate p) {
    // TODO: is there a faster way of doing this?
    Coordinate pRound = _round(p);

    HotPixel? hp = _find(pRound);
    /**
   * Hot Pixels which are added more than once
   * must have more than one vertex in them
   * and thus must be nodes.
   */
    if (hp != null) {
      hp.setToNode();
      return hp;
    }

    /**
   * A pixel containing the point was not found, so create a new one.
   * It is initially set to NOT be a node
   * (but may become one later on).
   */
    hp = new HotPixel(pRound, scaleFactor);
    index.insert(hp.getCoordinate(), hp);
    return hp;
  }

  HotPixel? _find(Coordinate pixelPt) {
    KdNode? kdNode = index.queryCoordinate(pixelPt);
    if (kdNode == null) return null;
    return kdNode.getData() as HotPixel;
  }

  Coordinate _round(Coordinate pt) {
    Coordinate p2 = pt.copy();
    precModel.makeCoordinatePrecise(p2);
    return p2;
  }

  /**
   * Visits all the hot pixels which may intersect a segment (p0-p1).
   * The visitor must determine whether each hot pixel actually intersects
   * the segment.
   *
   * @param p0 the segment start point
   * @param p1 the segment end point
   * @param visitor the visitor to apply
   */
  void query(Coordinate p0, Coordinate p1, KdNodeVisitor visitor) {
    Envelope queryEnv = Envelope.fromCoordinates(p0, p1);
    // expand query range to account for HotPixel extent
    // expand by full width of one pixel to be safe
    queryEnv.expandByDistance(1.0 / scaleFactor);
    index.query(queryEnv, visitor);
  }

  /**
   * Adds a list of points as non-node pixels.
   *
   * @param pts the points to add
   */
  void addList(List<Coordinate> pts) {
    /**
   * Shuffle the points before adding.
   * This avoids having long monontic runs of points
   * causing an unbalanced KD-tree, which would create
   * performance and robustness issues.
   */
    Iterator<Coordinate> it = CoordinateShuffler(pts).iterator;
    while (it.moveNext()) {
      addCoordinate(it.current);
    }
  }

  /**
   * Adds a list of points as node pixels.
   *
   * @param pts the points to add
   */
  void addNodesList(List<Coordinate> pts) {
    /**
   * Node points are not shuffled, since they are
   * added after the vertex points, and hence the KD-tree should
   * be reasonably balanced already.
   */
    for (Coordinate pt in pts) {
      HotPixel hp = addCoordinate(pt);
      hp.setToNode();
    }
  }
}

/**
 * Utility class to shuffle an array of {@link Coordinate}s using
 * the Fisher-Yates shuffle algorithm
 *
 * @see <a href="https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle">Fihser-Yates shuffle</a>
 */
class CoordinateShuffler implements Iterator<Coordinate> {
  final math.Random rnd = new math.Random(13);
  late List<Coordinate> coordinates;
  final List<int> indices = List.empty(growable: true);
  late int index;
  late Coordinate currentPt;

  /**
   * Creates an instance of this class
   * @param pts An array of {@link Coordinate}s.
   */
  CoordinateShuffler(List<Coordinate> pts) {
    coordinates = pts;
    for (int i = 0; i < pts.length; i++) indices.add(i);
    index = pts.length - 1;
  }

  @override
  bool moveNext() {
    int j = rnd.nextInt(index + 1);
    currentPt = coordinates[indices[j]];
    indices[j] = indices[index--];
    return index < 0;
  }

  Iterator<Coordinate> get iterator => coordinates.iterator;

  @override
  // TODO: implement current
  Coordinate get current => currentPt;
}

/**
 * Finds intersections between line segments which will be snap-rounded,
 * and adds them as nodes to the segments.
 * <p>
 * Intersections are detected and computed using full precision.
 * Snapping takes place in a subsequent phase.
 * <p>
 * The intersection points are recorded, so that HotPixels can be created for them.
 * <p>
 * To avoid robustness issues with vertices which lie very close to line segments
 * a heuristic is used:
 * nodes are created if a vertex lies within a tolerance distance
 * of the interior of a segment.
 * The tolerance distance is chosen to be significantly below the snap-rounding grid size.
 * This has empirically proven to eliminate noding failures.
 *
 * @version 1.17
 */
class SnapRoundingIntersectionAdder extends SegmentIntersectorN {
  late LineIntersector li;
  List<Coordinate> intersections = List.empty(growable: true);
  double nearnessTol;

  /**
   * Creates an intersector which finds all snapped interior intersections,
   * and adds them as nodes.
   *
   * @param nearnessTol the intersection distance tolerance
   */
  SnapRoundingIntersectionAdder(this.nearnessTol) {
    li = new RobustLineIntersector();
  }

  /**
   * Gets the created intersection nodes,
   * so they can be processed as hot pixels.
   *
   * @return a list of the intersection points
   */
  List<Coordinate> getIntersections() {
    return intersections;
  }

  /**
   * This method is called by clients
   * of the {@link SegmentIntersector} class to process
   * intersections for two segments of the {@link SegmentString}s being intersected.
   * Note that some clients (such as <code>MonotoneChain</code>s) may optimize away
   * this call for segment pairs which they have determined do not intersect
   * (e.g. by an disjoint envelope test).
   */
  void processIntersections(
      SegmentString e0, int segIndex0, SegmentString e1, int segIndex1) {
    // don't bother intersecting a segment with itself
    if (e0 == e1 && segIndex0 == segIndex1) return;

    Coordinate p00 = e0.getCoordinate(segIndex0);
    Coordinate p01 = e0.getCoordinate(segIndex0 + 1);
    Coordinate p10 = e1.getCoordinate(segIndex1);
    Coordinate p11 = e1.getCoordinate(segIndex1 + 1);

    li.computeIntersection(p00, p01, p10, p11);
//if (li.hasIntersection() && li.isProper()) Debug.println(li);

    if (li.hasIntersection()) {
      if (li.isInteriorIntersection()) {
        for (int intIndex = 0; intIndex < li.getIntersectionNum(); intIndex++) {
          intersections.add(li.getIntersection(intIndex));
        }
        (e0 as NodedSegmentString).addIntersections(li, segIndex0, 0);
        (e1 as NodedSegmentString).addIntersections(li, segIndex1, 1);
        return;
      }
    }

    /**
     * Segments did not actually intersect, within the limits of orientation index robustness.
     *
     * To avoid certain robustness issues in snap-rounding,
     * also treat very near vertex-segment situations as intersections.
     */
    _processNearVertex(p00, e1, segIndex1, p10, p11);
    _processNearVertex(p01, e1, segIndex1, p10, p11);
    _processNearVertex(p10, e0, segIndex0, p00, p01);
    _processNearVertex(p11, e0, segIndex0, p00, p01);
  }

  /**
   * If an endpoint of one segment is near
   * the <i>interior</i> of the other segment, add it as an intersection.
   * EXCEPT if the endpoint is also close to a segment endpoint
   * (since this can introduce "zigs" in the linework).
   * <p>
   * This resolves situations where
   * a segment A endpoint is extremely close to another segment B,
   * but is not quite crossing.  Due to robustness issues
   * in orientation detection, this can
   * result in the snapped segment A crossing segment B
   * without a node being introduced.
   *
   * @param p
   * @param edge
   * @param segIndex
   * @param p0
   * @param p1
   */
  void _processNearVertex(Coordinate p, SegmentString edge, int segIndex,
      Coordinate p0, Coordinate p1) {
    /**
     * Don't add intersection if candidate vertex is near endpoints of segment.
     * This avoids creating "zig-zag" linework
     * (since the vertex could actually be outside the segment envelope).
     */
    if (p.distance(p0) < nearnessTol) return;
    if (p.distance(p1) < nearnessTol) return;

    double distSeg = Distance.pointToSegment(p, p0, p1);
    if (distSeg < nearnessTol) {
      intersections.add(p);
      (edge as NodedSegmentString).addIntersection(p, segIndex);
    }
  }

  /**
   * Always process all intersections
   *
   * @return false always
   */
  bool isDone() {
    return false;
  }
}
