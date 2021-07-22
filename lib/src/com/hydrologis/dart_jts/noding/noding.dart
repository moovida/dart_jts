part of dart_jts;

/**
 * Allows comparing {@link Coordinate} arrays
 * in an orientation-independent way.
 *
 * @author Martin Davis
 * @version 1.7
 */
class OrientedCoordinateArray implements Comparable {
  List<Coordinate> pts;
  late bool _orientation;

  /**
   * Creates a new {@link OrientedCoordinateArray}
   * for the given {@link Coordinate} array.
   *
   * @param pts the coordinates to orient
   */
  OrientedCoordinateArray(this.pts) {
    _orientation = orientation(pts);
  }

  /**
   * Computes the canonical orientation for a coordinate array.
   *
   * @param pts the array to test
   * @return <code>true</code> if the points are oriented forwards
   * or <code>false</code if the points are oriented in reverse
   */
  static bool orientation(List<Coordinate> pts) {
    return CoordinateArrays.increasingDirection(pts) == 1;
  }

  /**
   * Compares two {@link OrientedCoordinateArray}s for their relative order
   *
   * @return -1 this one is smaller;
   * 0 the two objects are equal;
   * 1 this one is greater
   */

  int compareTo(var o1) {
    OrientedCoordinateArray oca = o1;
    int comp = compareOriented(pts, _orientation, oca.pts, oca._orientation);
/*
    // MD - testing only
    int oldComp = SegmentStringDissolver.ptsComp.compare(pts, oca.pts);
    if ((oldComp == 0 || comp == 0) && oldComp != comp) {
      System.out.println("bidir mismatch");

      bool orient1 = orientation(pts);
      bool orient2 = orientation(oca.pts);
      int comp2 = compareOriented(pts, orientation,
                               oca.pts, oca.orientation);
      int oldComp2 = SegmentStringDissolver.ptsComp.compare(pts, oca.pts);
    }
    */
    return comp;
  }

  static int compareOriented(List<Coordinate> pts1, bool orientation1,
      List<Coordinate> pts2, bool orientation2) {
    int dir1 = orientation1 ? 1 : -1;
    int dir2 = orientation2 ? 1 : -1;
    int limit1 = orientation1 ? pts1.length : -1;
    int limit2 = orientation2 ? pts2.length : -1;

    int i1 = orientation1 ? 0 : pts1.length - 1;
    int i2 = orientation2 ? 0 : pts2.length - 1;
    while (true) {
      int compPt = pts1[i1].compareTo(pts2[i2]);
      if (compPt != 0) return compPt;
      i1 += dir1;
      i2 += dir2;
      bool done1 = i1 == limit1;
      bool done2 = i2 == limit2;
      if (done1 && !done2) return -1;
      if (!done1 && done2) return 1;
      if (done1 && done2) return 0;
    }
  }
}

/**
 * Uses Snap Rounding to compute a rounded,
 * fully noded arrangement from a set of {@link SegmentString}s.
 * Implements the Snap Rounding technique described in
 * papers by Hobby, Guibas &amp; Marimont, and Goodrich et al.
 * Snap Rounding assumes that all vertices lie on a uniform grid;
 * hence the precision model of the input must be fixed precision,
 * and all the input vertices must be rounded to that precision.
 * <p>
 * This implementation uses a monotone chains and a spatial index to
 * speed up the intersection tests.
 * <p>
 * This implementation appears to be fully robust using an integer precision model.
 * It will function with non-integer precision models, but the
 * results are not 100% guaranteed to be correctly noded.
 *
 * @version 1.7
 */
class MCIndexSnapRounder implements Noder {
  PrecisionModel pm;
  late LineIntersector li;
  late double scaleFactor;
  late MCIndexNoder noder;
  late MCIndexPointSnapper pointSnapper;
  late List nodedSegStrings;

  MCIndexSnapRounder(this.pm) {
    li = new RobustLineIntersector();
    li.setPrecisionModel(pm);
    scaleFactor = pm.getScale();
  }

  List getNodedSubstrings() {
    return NodedSegmentString.getNodedSubstrings(nodedSegStrings);
  }

  void computeNodes(List inputSegmentStrings) {
    this.nodedSegStrings = inputSegmentStrings;
    noder = new MCIndexNoder.empty();
    pointSnapper = new MCIndexPointSnapper(noder.getIndex() as STRtree);
    snapRound(inputSegmentStrings, li);

    // testing purposes only - remove in final version
    //checkCorrectness(inputSegmentStrings);
  }

  void checkCorrectness(List inputSegmentStrings) {
    List resultSegStrings =
        NodedSegmentString.getNodedSubstrings(inputSegmentStrings);
    NodingValidator nv = new NodingValidator(resultSegStrings);
    try {
      nv.checkValid();
    } catch (ex) {
      print(ex);
    }
  }

  void snapRound(List segStrings, LineIntersector li) {
    List intersections = findInteriorIntersections(segStrings, li);
    computeIntersectionSnaps(intersections);
    computeVertexSnaps(segStrings);
  }

  /**
   * Computes all interior intersections in the collection of {@link SegmentString}s,
   * and returns their {@link Coordinate}s.
   *
   * Does NOT node the segStrings.
   *
   * @return a list of Coordinates for the intersections
   */
  List findInteriorIntersections(List segStrings, LineIntersector li) {
    InteriorIntersectionFinderAdder intFinderAdder =
        new InteriorIntersectionFinderAdder(li);
    noder.setSegmentIntersector(intFinderAdder);
    noder.computeNodes(segStrings);
    return intFinderAdder.getInteriorIntersections();
  }

  /**
   * Snaps segments to nodes created by segment intersections.
   */
  void computeIntersectionSnaps(List snapPts) {
    for (Coordinate snapPt in snapPts) {
      HotPixel hotPixel = new HotPixel(snapPt, scaleFactor, li);
      pointSnapper.snapPX(hotPixel);
    }
  }

  /**
   * Snaps segments to all vertices.
   *
   * @param edges the list of segment strings to snap together
   */
  void computeVertexSnaps(List edges) {
    for (NodedSegmentString edge0 in edges) {
      computeVertexSnapsNSS(edge0);
    }
  }

  /**
   * Snaps segments to the vertices of a Segment String.
   */
  void computeVertexSnapsNSS(NodedSegmentString e) {
    List<Coordinate> pts0 = e.getCoordinates();
    for (int i = 0; i < pts0.length; i++) {
      HotPixel hotPixel = new HotPixel(pts0[i], scaleFactor, li);
      bool isNodeAdded = pointSnapper.snap(hotPixel, e, i);
      // if a node is created for a vertex, that vertex must be noded too
      if (isNodeAdded) {
        e.addIntersection(pts0[i], i);
      }
    }
  }
}

/**
 * Nodes a set of {@link SegmentString}s using a index based
 * on {@link MonotoneChain}s and a {@link SpatialIndex}.
 * The {@link SpatialIndex} used should be something that supports
 * envelope (range) queries efficiently (such as a <code>Quadtree</code>}
 * or {@link STRtree} (which is the default index provided).
 *
 * @version 1.7
 */
class MCIndexNoder extends SinglePassNoder {
  List monoChains = [];
  SpatialIndex index = new STRtree();
  int idCounter = 0;
  late List nodedSegStrings;

  // statistics
  int nOverlaps = 0;

  MCIndexNoder.empty() : super.empty();

  MCIndexNoder(SegmentIntersectorN si) : super(si);

  List getMonotoneChains() {
    return monoChains;
  }

  SpatialIndex getIndex() {
    return index;
  }

  List getNodedSubstrings() {
    return NodedSegmentString.getNodedSubstrings(nodedSegStrings);
  }

  void computeNodes(List inputSegStrings) {
    this.nodedSegStrings = inputSegStrings;

    for (SegmentString ss in inputSegStrings) {
      add(ss);
    }
    intersectChains();
//System.out.println("MCIndexNoder: # chain overlaps = " + nOverlaps);
  }

  void intersectChains() {
    MonotoneChainOverlapAction overlapAction =
        new SegmentOverlapAction(segInt!);

    for (MonotoneChainI queryChain in monoChains) {
      List overlapChains = index.query(queryChain.getEnvelope());
      for (MonotoneChainI testChain in overlapChains) {
        /**
         * following test makes sure we only compare each pair of chains once
         * and that we don't compare a chain to itself
         */
        if (testChain.getId() > queryChain.getId()) {
          queryChain.computeOverlaps3(testChain, overlapAction);
          nOverlaps++;
        }
        // short-circuit if possible
        if (segInt!.isDone()) return;
      }
    }
  }

  void add(SegmentString segStr) {
    List segChains = MonotoneChainBuilder.getChainsWithContext(
        segStr.getCoordinates(), segStr);
    for (MonotoneChainI mc in segChains) {
      mc.setId(idCounter++);
      index.insert(mc.getEnvelope(), mc);
      monoChains.add(mc);
    }
  }
}

class SegmentOverlapAction extends MonotoneChainOverlapAction {
  SegmentIntersectorN si;

  SegmentOverlapAction(this.si);

  void overlap(MonotoneChainI mc1, int start1, MonotoneChainI mc2, int start2) {
    SegmentString ss1 = mc1.getContext() as SegmentString;
    SegmentString ss2 = mc2.getContext() as SegmentString;
    si.processIntersections(ss1, start1, ss2, start2);
  }
}

/**
 * Processes possible intersections detected by a {@link Noder}.
 * The {@link SegmentIntersector} is passed to a {@link Noder}.
 * The {@link SegmentIntersector#processIntersections(SegmentString, int, SegmentString, int)} method is called whenever the {@link Noder}
 * detects that two SegmentStrings <i>might</i> intersect.
 * This class may be used either to find all intersections, or
 * to detect the presence of an intersection.  In the latter case,
 * Noders may choose to short-circuit their computation by calling the
 * {@link #isDone()} method.
 * This class is an example of the <i>Strategy</i> pattern.
 *
 * @version 1.7
 */
abstract class SegmentIntersectorN {
  /**
   * This method is called by clients
   * of the {@link SegmentIntersector} interface to process
   * intersections for two segments of the {@link SegmentString}s being intersected.
   */
  void processIntersections(
      SegmentString e0, int segIndex0, SegmentString e1, int segIndex1);

  /**
   * Reports whether the client of this class
   * needs to continue testing all intersections in an arrangement.
   *
   * @return true if there is no need to continue testing segments
   */
  bool isDone();
}

/**
 * An interface for classes which represent a sequence of contiguous line segments.
 * SegmentStrings can carry a context object, which is useful
 * for preserving topological or parentage information.
 *
 * @version 1.7
 */
abstract class SegmentString {
  /**
   * Gets the user-defined data for this segment string.
   *
   * @return the user-defined data
   */
  Object getData();

  /**
   * Sets the user-defined data for this segment string.
   *
   * @param data an Object containing user-defined data
   */
  void setData(Object data);

  int size();

  Coordinate getCoordinate(int i);

  List<Coordinate> getCoordinates();

  bool isClosed();
}

/**
 * Base class for {@link Noder}s which make a single
 * pass to find intersections.
 * This allows using a custom {@link SegmentIntersector}
 * (which for instance may simply identify intersections, rather than
 * insert them).
 *
 * @version 1.7
 */
abstract class SinglePassNoder implements Noder {
  SegmentIntersectorN? segInt;

  SinglePassNoder.empty() {}

  SinglePassNoder(SegmentIntersectorN segInt) {
    setSegmentIntersector(segInt);
  }

  /**
   * Sets the SegmentIntersector to use with this noder.
   * A SegmentIntersector will normally add intersection nodes
   * to the input segment strings, but it may not - it may
   * simply record the presence of intersections.
   * However, some Noders may require that intersections be added.
   *
   * @param segInt
   */
  void setSegmentIntersector(SegmentIntersectorN segInt) {
    this.segInt = segInt;
  }

  /**
   * Computes the noding for a collection of {@link SegmentString}s.
   * Some Noders may add all these nodes to the input SegmentStrings;
   * others may only add some or none at all.
   *
   * @param segStrings a collection of {@link SegmentString}s to node
   */
  void computeNodes(List segStrings);

  /**
   * Returns a {@link List} of fully noded {@link SegmentString}s.
   * The SegmentStrings have the same context as their parent.
   *
   * @return a List of SegmentStrings
   */
  List getNodedSubstrings();
}

/**
 * Computes all intersections between segments in a set of {@link SegmentString}s.
 * Intersections found are represented as {@link SegmentNode}s and added to the
 * {@link SegmentString}s in which they occur.
 * As a final step in the noding a new set of segment strings split
 * at the nodes may be returned.
 *
 * @version 1.7
 */
abstract class Noder {
  /**
   * Computes the noding for a collection of {@link SegmentString}s.
   * Some Noders may add all these nodes to the input SegmentStrings;
   * others may only add some or none at all.
   *
   * @param segStrings a collection of {@link SegmentString}s to node
   */
  void computeNodes(List segStrings);

  /**
   * Returns a {@link List} of fully noded {@link SegmentString}s.
   * The SegmentStrings have the same context as their parent.
   *
   * @return a List of SegmentStrings
   */
  List getNodedSubstrings();
}

/**
 * "Snaps" all {@link SegmentString}s in a {@link SpatialIndex} containing
 * {@link MonotoneChain}s to a given {@link HotPixel}.
 *
 * @version 1.7
 */
class MCIndexPointSnapper {
  // static final int nSnaps = 0;

  STRtree index;

  MCIndexPointSnapper(this.index);

  /**
   * Snaps (nodes) all interacting segments to this hot pixel.
   * The hot pixel may represent a vertex of an edge,
   * in which case this routine uses the optimization
   * of not noding the vertex itself
   *
   * @param hotPixel the hot pixel to snap to
   * @param parentEdge the edge containing the vertex, if applicable, or <code>null</code>
   * @param hotPixelVertexIndex the index of the hotPixel vertex, if applicable, or -1
   * @return <code>true</code> if a node was added for this pixel
   */
  bool snap(
      HotPixel hotPixel, SegmentString? parentEdge, int hotPixelVertexIndex) {
    final Envelope pixelEnv = hotPixel.getSafeEnvelope();
    final HotPixelSnapAction hotPixelSnapAction =
        new HotPixelSnapAction(hotPixel, parentEdge, hotPixelVertexIndex);

    index.queryWithVisitor(
        pixelEnv, MonotoneChainIVisitor(pixelEnv, hotPixelSnapAction));
    return hotPixelSnapAction.isNodeAdded();
  }

  bool snapPX(HotPixel hotPixel) {
    return snap(hotPixel, null, -1);
  }
}

class MonotoneChainIVisitor implements ItemVisitor {
  final Envelope pixelEnv;
  final HotPixelSnapAction hotPixelSnapAction;

  MonotoneChainIVisitor(this.pixelEnv, this.hotPixelSnapAction);

  /**
   * Visits an item.
   *
   * @param item the item to visit
   */
  void visitItem(Object item) {
    MonotoneChainI testChain = item as MonotoneChainI;
    testChain.select(pixelEnv, hotPixelSnapAction);
  }
}

class HotPixelSnapAction extends MonotoneChainSelectAction {
  HotPixel hotPixel;
  SegmentString? parentEdge;

// is -1 if hotPixel is not a vertex
  int hotPixelVertexIndex;
  bool _isNodeAdded = false;

  HotPixelSnapAction(this.hotPixel, this.parentEdge, this.hotPixelVertexIndex);

  /**
   * Reports whether the HotPixel caused a node to be added in any target
   * segmentString (including its own). If so, the HotPixel must be added as a
   * node as well.
   *
   * @return true if a node was added in any target segmentString.
   */
  bool isNodeAdded() {
    return _isNodeAdded;
  }

  /**
   * Check if a segment of the monotone chain intersects
   * the hot pixel vertex and introduce a snap node if so.
   * Optimized to avoid noding segments which
   * contain the vertex (which otherwise
   * would cause every vertex to be noded).
   */
  void select(MonotoneChainI mc, int startIndex) {
    NodedSegmentString ss = mc.getContext() as NodedSegmentString;
    /**
     * Check to avoid snapping a hotPixel vertex to the same vertex.
     * This method is called for segments which intersects the
     * hot pixel,
     * so need to check if either end of the segment is equal to the hot pixel
     * and if so, do not snap.
     *
     * Sep 22 2012 - MD - currently do need to snap to every vertex,
     * since otherwise the testCollapse1 test in SnapRoundingTest fails.
     */
    if (parentEdge == ss) {
// exit if hotpixel is equal to endpoint of target segment
      if (startIndex == hotPixelVertexIndex ||
          startIndex + 1 == hotPixelVertexIndex) return;
    }
// snap and record if a node was created
    _isNodeAdded |= hotPixel.addSnappedNode(ss, startIndex);
  }
}

/**
 * Implements a "hot pixel" as used in the Snap Rounding algorithm.
 * A hot pixel contains the interior of the tolerance square and
 * the boundary
 * <b>minus</b> the top and right segments.
 * <p>
 * The hot pixel operations are all computed in the integer domain
 * to avoid rounding problems.
 *
 * @version 1.7
 */
class HotPixel {
  // testing only
//   static int nTests = 0;

  LineIntersector li;

  Coordinate pt;
  late Coordinate originalPt;
  Coordinate? ptScaled;

  Coordinate? p0Scaled;
  Coordinate? p1Scaled;

  double scaleFactor;

  double minx = 0.0;
  double maxx = 0.0;
  double miny = 0.0;
  double maxy = 0.0;

  /**
   * The corners of the hot pixel, in the order:
   *  10
   *  23
   */
  List<Coordinate> corner = []; //..length = (4);

  Envelope? safeEnv = null;

  /**
   * Creates a new hot pixel, using a given scale factor.
   * The scale factor must be strictly positive (non-zero).
   *
   * @param pt the coordinate at the centre of the pixel
   * @param scaleFactor the scaleFactor determining the pixel size.  Must be &gt; 0
   * @param li the intersector to use for testing intersection with line segments
   *
   */
  HotPixel(this.pt, this.scaleFactor, this.li) {
    originalPt = pt;
    this.scaleFactor = scaleFactor;
    this.li = li;
    //tolerance = 0.5;
    if (scaleFactor <= 0) throw ArgumentError("Scale factor must be non-zero");
    if (scaleFactor != 1.0) {
      this.pt = new Coordinate(scale(pt.x), scale(pt.y));
      p0Scaled = new Coordinate.empty2D();
      p1Scaled = new Coordinate.empty2D();
    }
    initCorners(this.pt);
  }

  /**
   * Gets the coordinate this hot pixel is based at.
   *
   * @return the coordinate of the pixel
   */
  Coordinate getCoordinate() {
    return originalPt;
  }

  static final double SAFE_ENV_EXPANSION_FACTOR = 0.75;

  /**
   * Returns a "safe" envelope that is guaranteed to contain the hot pixel.
   * The envelope returned will be larger than the exact envelope of the
   * pixel.
   *
   * @return an envelope which contains the hot pixel
   */
  Envelope getSafeEnvelope() {
    if (safeEnv == null) {
      double safeTolerance = SAFE_ENV_EXPANSION_FACTOR / scaleFactor;
      safeEnv = new Envelope(
          originalPt.x - safeTolerance,
          originalPt.x + safeTolerance,
          originalPt.y - safeTolerance,
          originalPt.y + safeTolerance);
    }
    return safeEnv!;
  }

  void initCorners(Coordinate pt) {
    double tolerance = 0.5;
    minx = pt.x - tolerance;
    maxx = pt.x + tolerance;
    miny = pt.y - tolerance;
    maxy = pt.y + tolerance;

    corner.add(new Coordinate(maxx, maxy));
    corner.add(new Coordinate(minx, maxy));
    corner.add(new Coordinate(minx, miny));
    corner.add(new Coordinate(maxx, miny));
    // corner[0] = new Coordinate(maxx, maxy);
    // corner[1] = new Coordinate(minx, maxy);
    // corner[2] = new Coordinate(minx, miny);
    // corner[3] = new Coordinate(maxx, miny);
  }

  double scale(double val) {
    return (val * scaleFactor).roundToDouble();
  }

  /**
   * Tests whether the line segment (p0-p1)
   * intersects this hot pixel.
   *
   * @param p0 the first coordinate of the line segment to test
   * @param p1 the second coordinate of the line segment to test
   * @return true if the line segment intersects this hot pixel
   */
  bool intersects(Coordinate p0, Coordinate p1) {
    if (scaleFactor == 1.0) return intersectsScaled(p0, p1);

    copyScaled(p0, p0Scaled!);
    copyScaled(p1, p1Scaled!);
    return intersectsScaled(p0Scaled!, p1Scaled!);
  }

  void copyScaled(Coordinate p, Coordinate pScaled) {
    pScaled.x = scale(p.x);
    pScaled.y = scale(p.y);
  }

  bool intersectsScaled(Coordinate p0, Coordinate p1) {
    double segMinx = math.min(p0.x, p1.x);
    double segMaxx = math.max(p0.x, p1.x);
    double segMiny = math.min(p0.y, p1.y);
    double segMaxy = math.max(p0.y, p1.y);

    bool isOutsidePixelEnv =
        maxx < segMinx || minx > segMaxx || maxy < segMiny || miny > segMaxy;
    if (isOutsidePixelEnv) return false;
    bool intersects = intersectsToleranceSquare(p0, p1);
//    bool intersectsPixelClosure = intersectsPixelClosure(p0, p1);

//    if (intersectsPixel != intersects) {
//      Debug.println("Found hot pixel intersection mismatch at " + pt);
//      Debug.println("Test segment: " + p0 + " " + p1);
//    }

/*
    if (scaleFactor != 1.0) {
      bool intersectsScaled = intersectsScaledTest(p0, p1);
      if (intersectsScaled != intersects) {
        intersectsScaledTest(p0, p1);
//        Debug.println("Found hot pixel scaled intersection mismatch at " + pt);
//        Debug.println("Test segment: " + p0 + " " + p1);
      }
      return intersectsScaled;
    }
*/

    Assert.isTrue(
        !(isOutsidePixelEnv && intersects), "Found bad envelope test");
//    if (isOutsideEnv && intersects) {
//      Debug.println("Found bad envelope test");
//    }

    return intersects;
    //return intersectsPixelClosure;
  }

  /**
   * Tests whether the segment p0-p1 intersects the hot pixel tolerance square.
   * Because the tolerance square point set is partially open (along the
   * top and right) the test needs to be more sophisticated than
   * simply checking for any intersection.
   * However, it can take advantage of the fact that the hot pixel edges
   * do not lie on the coordinate grid.
   * It is sufficient to check if any of the following occur:
   * <ul>
   * <li>a proper intersection between the segment and any hot pixel edge
   * <li>an intersection between the segment and <b>both</b> the left and bottom hot pixel edges
   * (which detects the case where the segment intersects the bottom left hot pixel corner)
   * <li>an intersection between a segment endpoint and the hot pixel coordinate
   * </ul>
   *
   * @param p0
   * @param p1
   * @return
   */
  bool intersectsToleranceSquare(Coordinate p0, Coordinate p1) {
    bool intersectsLeft = false;
    bool intersectsBottom = false;
    //System.out.println("Hot Pixel: " + WKTWriter.toLineString(corner));
    //System.out.println("Line: " + WKTWriter.toLineString(p0, p1));

    li.computeIntersection(p0, p1, corner[0], corner[1]);
    if (li.isProper()) return true;

    li.computeIntersection(p0, p1, corner[1], corner[2]);
    if (li.isProper()) return true;
    if (li.hasIntersection()) intersectsLeft = true;

    li.computeIntersection(p0, p1, corner[2], corner[3]);
    if (li.isProper()) return true;
    if (li.hasIntersection()) intersectsBottom = true;

    li.computeIntersection(p0, p1, corner[3], corner[0]);
    if (li.isProper()) return true;

    if (intersectsLeft && intersectsBottom) return true;

    if (p0.equals(pt)) return true;
    if (p1.equals(pt)) return true;

    return false;
  }

  /**
   * Test whether the given segment intersects
   * the closure of this hot pixel.
   * This is NOT the test used in the standard snap-rounding
   * algorithm, which uses the partially closed tolerance square
   * instead.
   * This routine is provided for testing purposes only.
   *
   * @param p0 the start point of a line segment
   * @param p1 the end point of a line segment
   * @return <code>true</code> if the segment intersects the closure of the pixel's tolerance square
   */
  bool intersectsPixelClosure(Coordinate p0, Coordinate p1) {
    li.computeIntersection(p0, p1, corner[0], corner[1]);
    if (li.hasIntersection()) return true;
    li.computeIntersection(p0, p1, corner[1], corner[2]);
    if (li.hasIntersection()) return true;
    li.computeIntersection(p0, p1, corner[2], corner[3]);
    if (li.hasIntersection()) return true;
    li.computeIntersection(p0, p1, corner[3], corner[0]);
    if (li.hasIntersection()) return true;

    return false;
  }

  /**
   * Adds a new node (equal to the snap pt) to the specified segment
   * if the segment passes through the hot pixel
   *
   * @param segStr
   * @param segIndex
   * @return true if a node was added to the segment
   */
  bool addSnappedNode(NodedSegmentString segStr, int segIndex) {
    Coordinate p0 = segStr.getCoordinate(segIndex);
    Coordinate p1 = segStr.getCoordinate(segIndex + 1);

    if (intersects(p0, p1)) {
      //System.out.println("snapped: " + snapPt);
      //System.out.println("POINT (" + snapPt.x + " " + snapPt.y + ")");
      segStr.addIntersection(getCoordinate(), segIndex);

      return true;
    }
    return false;
  }
}

/**
 * Represents a list of contiguous line segments,
 * and supports noding the segments.
 * The line segments are represented by an array of {@link Coordinate}s.
 * Intended to optimize the noding of contiguous segments by
 * reducing the number of allocated objects.
 * SegmentStrings can carry a context object, which is useful
 * for preserving topological or parentage information.
 * All noded substrings are initialized with the same context object.
 *
 * @version 1.7
 */
class NodedSegmentString implements NodableSegmentString {
  /**
   * Gets the {@link SegmentString}s which result from splitting this string at node points.
   *
   * @param segStrings a List of NodedSegmentStrings
   * @return a List of NodedSegmentStrings representing the substrings
   */
  static List getNodedSubstrings(List segStrings) {
    List resultEdgelist = [];
    getNodedSubstrings2(segStrings, resultEdgelist);
    return resultEdgelist;
  }

  /**
   * Adds the noded {@link SegmentString}s which result from splitting this string at node points.
   *
   * @param segStrings a List of NodedSegmentStrings
   * @param resultEdgelist a List which will collect the NodedSegmentStrings representing the substrings
   */
  static void getNodedSubstrings2(List segStrings, List resultEdgelist) {
    for (NodedSegmentString ss in segStrings) {
      ss.getNodeList().addSplitEdges(resultEdgelist);
    }
  }

  late SegmentNodeList nodeList;
  List<Coordinate> pts;
  Object data;

  /**
   * Creates a new segment string from a list of vertices.
   *
   * @param pts the vertices of the segment string
   * @param data the user-defined data of this segment string (may be null)
   */
  NodedSegmentString(this.pts, this.data) {
    nodeList = new SegmentNodeList(this);
  }

  /**
   * Gets the user-defined data for this segment string.
   *
   * @return the user-defined data
   */
  Object getData() {
    return data;
  }

  /**
   * Sets the user-defined data for this segment string.
   *
   * @param data an Object containing user-defined data
   */
  void setData(Object data) {
    this.data = data;
  }

  SegmentNodeList getNodeList() {
    return nodeList;
  }

  int size() {
    return pts.length;
  }

  Coordinate getCoordinate(int i) {
    return pts[i];
  }

  List<Coordinate> getCoordinates() {
    return pts;
  }

  bool isClosed() {
    return pts[0].equals(pts[pts.length - 1]);
  }

  /**
   * Gets the octant of the segment starting at vertex <code>index</code>.
   *
   * @param index the index of the vertex starting the segment.  Must not be
   * the last index in the vertex list
   * @return the octant of the segment at the vertex
   */
  int getSegmentOctant(int index) {
    if (index == pts.length - 1) return -1;
    return safeOctant(getCoordinate(index), getCoordinate(index + 1));
//    return Octant.octant(getCoordinate(index), getCoordinate(index + 1));
  }

  int safeOctant(Coordinate p0, Coordinate p1) {
    if (p0.equals2D(p1)) return 0;
    return Octant.octantCoords(p0, p1);
  }

  /**
   * Adds EdgeIntersections for one or both
   * intersections found for a segment of an edge to the edge intersection list.
   */
  void addIntersections(LineIntersector li, int segmentIndex, int geomIndex) {
    for (int i = 0; i < li.getIntersectionNum(); i++) {
      addIntersectionLI(li, segmentIndex, geomIndex, i);
    }
  }

  /**
   * Add an SegmentNode for intersection intIndex.
   * An intersection that falls exactly on a vertex
   * of the SegmentString is normalized
   * to use the higher of the two possible segmentIndexes
   */
  void addIntersectionLI(
      LineIntersector li, int segmentIndex, int geomIndex, int intIndex) {
    Coordinate intPt =
        new Coordinate.fromCoordinate(li.getIntersection(intIndex));
    addIntersection(intPt, segmentIndex);
  }

  /**
   * Adds an intersection node for a given point and segment to this segment string.
   *
   * @param intPt the location of the intersection
   * @param segmentIndex the index of the segment containing the intersection
   */
  void addIntersection(Coordinate intPt, int segmentIndex) {
    addIntersectionNode(intPt, segmentIndex);
  }

  /**
   * Adds an intersection node for a given point and segment to this segment string.
   * If an intersection already exists for this exact location, the existing
   * node will be returned.
   *
   * @param intPt the location of the intersection
   * @param segmentIndex the index of the segment containing the intersection
   * @return the intersection node for the point
   */
  SegmentNode addIntersectionNode(Coordinate intPt, int segmentIndex) {
    int normalizedSegmentIndex = segmentIndex;
    //Debug.println("edge intpt: " + intPt + " dist: " + dist);
    // normalize the intersection point location
    int nextSegIndex = normalizedSegmentIndex + 1;
    if (nextSegIndex < pts.length) {
      Coordinate nextPt = pts[nextSegIndex];
      //Debug.println("next pt: " + nextPt);

      // Normalize segment index if intPt falls on vertex
      // The check for point equality is 2D only - Z values are ignored
      if (intPt.equals2D(nextPt)) {
        //Debug.println("normalized distance");
        normalizedSegmentIndex = nextSegIndex;
      }
    }
    /**
     * Add the intersection point to edge intersection list.
     */
    SegmentNode ei = nodeList.add(intPt, normalizedSegmentIndex);
    return ei;
  }

  String toString() {
    return WKTWriter.toLineStringFromSequence(CoordinateArraySequence(pts));
  }
}

/**
 * An interface for classes which support adding nodes to
 * a segment string.
 *
 * @author Martin Davis
 */
abstract class NodableSegmentString extends SegmentString {
  /**
   * Adds an intersection node for a given point and segment to this segment string.
   *
   * @param intPt the location of the intersection
   * @param segmentIndex the index of the segment containing the intersection
   */
  void addIntersection(Coordinate intPt, int segmentIndex);
}

/**
 * Represents an intersection point between two {@link SegmentString}s.
 *
 * @version 1.7
 */
class SegmentNode implements Comparable<SegmentNode> {
  NodedSegmentString segString;
  late Coordinate coord; // the point of intersection
  int segmentIndex; // the index of the containing line segment in the parent edge
  int segmentOctant;
  bool _isInterior = false;

  SegmentNode(
      this.segString, Coordinate coord, this.segmentIndex, this.segmentOctant) {
    this.coord = new Coordinate.fromCoordinate(coord);
    _isInterior = !coord.equals2D(segString.getCoordinate(segmentIndex));
  }

  /**
   * Gets the {@link Coordinate} giving the location of this node.
   *
   * @return the coordinate of the node
   */
  Coordinate getCoordinate() {
    return coord;
  }

  bool isInterior() {
    return _isInterior;
  }

  bool isEndPoint(int maxSegmentIndex) {
    if (segmentIndex == 0 && !_isInterior) return true;
    if (segmentIndex == maxSegmentIndex) return true;
    return false;
  }

  /**
   * @return -1 this SegmentNode is located before the argument location;
   * 0 this SegmentNode is at the argument location;
   * 1 this SegmentNode is located after the argument location
   */
  int compareTo(SegmentNode obj) {
    SegmentNode other = obj;

    if (segmentIndex < other.segmentIndex) return -1;
    if (segmentIndex > other.segmentIndex) return 1;

    if (coord.equals2D(other.coord)) return 0;

    // an exterior node is the segment start point, so always sorts first
    // this guards against a robustness problem where the octants are not reliable
    if (!_isInterior) return -1;
    if (!other._isInterior) return 1;

    return SegmentPointComparator.compare(segmentOctant, coord, other.coord);
    //return segment.compareNodePosition(this, other);
  }

  String toString() {
    return "$segmentIndex:" + coord.toString();
  }
}

/**
 * Implements a robust method of comparing the relative position of two
 * points along the same segment.
 * The coordinates are assumed to lie "near" the segment.
 * This means that this algorithm will only return correct results
 * if the input coordinates
 * have the same precision and correspond to rounded values
 * of exact coordinates lying on the segment.
 *
 * @version 1.7
 */
class SegmentPointComparator {
  /**
   * Compares two {@link Coordinate}s for their relative position along a segment
   * lying in the specified {@link Octant}.
   *
   * @return -1 node0 occurs first;
   * 0 the two nodes are equal;
   * 1 node1 occurs first
   */
  static int compare(int octant, Coordinate p0, Coordinate p1) {
    // nodes can only be equal if their coordinates are equal
    if (p0.equals2D(p1)) return 0;

    int xSign = relativeSign(p0.x, p1.x);
    int ySign = relativeSign(p0.y, p1.y);

    switch (octant) {
      case 0:
        return compareValue(xSign, ySign);
      case 1:
        return compareValue(ySign, xSign);
      case 2:
        return compareValue(ySign, -xSign);
      case 3:
        return compareValue(-xSign, ySign);
      case 4:
        return compareValue(-xSign, -ySign);
      case 5:
        return compareValue(-ySign, -xSign);
      case 6:
        return compareValue(-ySign, xSign);
      case 7:
        return compareValue(xSign, -ySign);
    }
    Assert.shouldNeverReachHere("invalid octant value");
    return 0;
  }

  static int relativeSign(double x0, double x1) {
    if (x0 < x1) return -1;
    if (x0 > x1) return 1;
    return 0;
  }

  static int compareValue(int compareSign0, int compareSign1) {
    if (compareSign0 < 0) return -1;
    if (compareSign0 > 0) return 1;
    if (compareSign1 < 0) return -1;
    if (compareSign1 > 0) return 1;
    return 0;
  }
}

/**
 * A list of the {@link SegmentNode}s present along a noded {@link SegmentString}.
 *
 * @version 1.7
 */
class SegmentNodeList {
  Map nodeMap = new SplayTreeMap();
  NodedSegmentString edge; // the parent edge

  SegmentNodeList(this.edge);

  NodedSegmentString getEdge() {
    return edge;
  }

  /**
   * Adds an intersection into the list, if it isn't already there.
   * The input segmentIndex and dist are expected to be normalized.
   *
   * @return the SegmentIntersection found or added
   */
  SegmentNode add(Coordinate intPt, int segmentIndex) {
    SegmentNode eiNew = new SegmentNode(
        edge, intPt, segmentIndex, edge.getSegmentOctant(segmentIndex));
    SegmentNode? ei = nodeMap[eiNew];
    if (ei != null) {
      // debugging sanity check
      Assert.isTrue(ei.coord.equals2D(intPt),
          "Found equal nodes with different coordinates");
//      if (! ei.coord.equals2D(intPt))
//        Debug.println("Found equal nodes with different coordinates");

      return ei;
    }
    // node does not exist, so create it
    nodeMap[eiNew] = eiNew;
    return eiNew;
  }

  /**
   * returns an iterator of SegmentNodes
   */
  Iterator iterator() {
    return nodeMap.values.iterator;
  }

  /**
   * Adds nodes for the first and last points of the edge
   */
  void addEndpoints() {
    int maxSegIndex = edge.size() - 1;
    add(edge.getCoordinate(0), 0);
    add(edge.getCoordinate(maxSegIndex), maxSegIndex);
  }

  /**
   * Adds nodes for any collapsed edge pairs.
   * Collapsed edge pairs can be caused by inserted nodes, or they can be
   * pre-existing in the edge vertex list.
   * In order to provide the correct fully noded semantics,
   * the vertex at the base of a collapsed pair must also be added as a node.
   */
  void addCollapsedNodes() {
    List<int> collapsedVertexIndexes = [];

    findCollapsesFromInsertedNodes(collapsedVertexIndexes);
    findCollapsesFromExistingVertices(collapsedVertexIndexes);

    // node the collapses
    for (Iterator it = collapsedVertexIndexes.iterator; it.moveNext();) {
      int vertexIndex = it.current;
      add(edge.getCoordinate(vertexIndex), vertexIndex);
    }
  }

  /**
   * Adds nodes for any collapsed edge pairs
   * which are pre-existing in the vertex list.
   */
  void findCollapsesFromExistingVertices(List<int> collapsedVertexIndexes) {
    for (int i = 0; i < edge.size() - 2; i++) {
      Coordinate p0 = edge.getCoordinate(i);
      Coordinate p1 = edge.getCoordinate(i + 1);
      Coordinate p2 = edge.getCoordinate(i + 2);
      if (p0.equals2D(p2)) {
        // add base of collapse as node
        collapsedVertexIndexes.add(i + 1);
      }
    }
  }

  /**
   * Adds nodes for any collapsed edge pairs caused by inserted nodes
   * Collapsed edge pairs occur when the same coordinate is inserted as a node
   * both before and after an existing edge vertex.
   * To provide the correct fully noded semantics,
   * the vertex must be added as a node as well.
   */
  void findCollapsesFromInsertedNodes(List<int> collapsedVertexIndexes) {
    List<int> collapsedVertexIndex = List.filled(1, 0);
    Iterator it = iterator();
    // there should always be at least two entries in the list, since the endpoints are nodes
    it.moveNext();
    SegmentNode eiPrev = it.current;
    while (it.moveNext()) {
      SegmentNode ei = it.current;
      bool isCollapsed = findCollapseIndex(eiPrev, ei, collapsedVertexIndex);
      if (isCollapsed) collapsedVertexIndexes.add(collapsedVertexIndex[0]);

      eiPrev = ei;
    }
  }

  bool findCollapseIndex(
      SegmentNode ei0, SegmentNode ei1, List<int> collapsedVertexIndex) {
    // only looking for equal nodes
    if (!ei0.coord.equals2D(ei1.coord)) return false;

    int numVerticesBetween = ei1.segmentIndex - ei0.segmentIndex;
    if (!ei1.isInterior()) {
      numVerticesBetween--;
    }

    // if there is a single vertex between the two equal nodes, this is a collapse
    if (numVerticesBetween == 1) {
      collapsedVertexIndex[0] = ei0.segmentIndex + 1;
      return true;
    }
    return false;
  }

  /**
   * Creates new edges for all the edges that the intersections in this
   * list split the parent edge into.
   * Adds the edges to the provided argument list
   * (this is so a single list can be used to accumulate all split edges
   * for a set of {@link SegmentString}s).
   */
  void addSplitEdges(List edgeList) {
    // ensure that the list has entries for the first and last point of the edge
    addEndpoints();
    addCollapsedNodes();

    Iterator it = iterator();
    // there should always be at least two entries in the list, since the endpoints are nodes
    it.moveNext();
    SegmentNode eiPrev = it.current;
    while (it.moveNext()) {
      SegmentNode ei = it.current;
      SegmentString newEdge = createSplitEdge(eiPrev, ei);
      /*
      if (newEdge.size() < 2)
        throw new RuntimeException("created single point edge: " + newEdge.toString());
      */
      edgeList.add(newEdge);
      eiPrev = ei;
    }
    //checkSplitEdgesCorrectness(testingSplitEdges);
  }

  /**
   * Checks the correctness of the set of split edges corresponding to this edge.
   *
   * @param splitEdges the split edges for this edge (in order)
   */
  void checkSplitEdgesCorrectness(List splitEdges) {
    List<Coordinate> edgePts = edge.getCoordinates();

    // check that first and last points of split edges are same as endpoints of edge
    SegmentString split0 = splitEdges[0];
    Coordinate pt0 = split0.getCoordinate(0);
    if (!pt0.equals2D(edgePts[0]))
      throw new RuntimeException(
          "bad split edge start point at " + pt0.toString());

    SegmentString splitn = splitEdges[splitEdges.length - 1];
    List<Coordinate> splitnPts = splitn.getCoordinates();
    Coordinate ptn = splitnPts[splitnPts.length - 1];
    if (!ptn.equals2D(edgePts[edgePts.length - 1]))
      throw new RuntimeException(
          "bad split edge end point at " + ptn.toString());
  }

  /**
   * Create a new "split edge" with the section of points between
   * (and including) the two intersections.
   * The label for the new edge is the same as the label for the parent edge.
   */
  SegmentString createSplitEdge(SegmentNode ei0, SegmentNode ei1) {
    List<Coordinate> pts = createSplitEdgePts(ei0, ei1);
    return new NodedSegmentString(pts, edge.getData());
  }

  /**
   * Extracts the points for a split edge running between two nodes.
   * The extracted points should contain no duplicate points.
   * There should always be at least two points extracted
   * (which will be the given nodes).
   *
   * @param ei0 the start node of the split edge
   * @param ei1 the end node of the split edge
   * @return the points for the split edge
   */
  List<Coordinate> createSplitEdgePts(SegmentNode ei0, SegmentNode ei1) {
//Debug.println("\ncreateSplitEdge"); Debug.print(ei0); Debug.print(ei1);
    int npts = ei1.segmentIndex - ei0.segmentIndex + 2;

    // if only two points in split edge they must be the node points
    if (npts == 2)
      return [
        new Coordinate.fromCoordinate(ei0.coord),
        new Coordinate.fromCoordinate(ei1.coord)
      ];

    Coordinate lastSegStartPt = edge.getCoordinate(ei1.segmentIndex);
    /**
     * If the last intersection point is not equal to the its segment start pt,
     * add it to the points list as well.
     * This check is needed because the distance metric is not totally reliable!
     *
     * Also ensure that the created edge always has at least 2 points.
     *
     * The check for point equality is 2D only - Z values are ignored
     */
    bool useIntPt1 = ei1.isInterior() || !ei1.coord.equals2D(lastSegStartPt);
    if (!useIntPt1) {
      npts--;
    }

    List<Coordinate> pts = []; //..length = (npts);
    // int ipt = 0;
    pts.add(new Coordinate.fromCoordinate(ei0.coord));
    // pts[ipt++] = new Coordinate.fromCoordinate(ei0.coord);
    for (int i = ei0.segmentIndex + 1; i <= ei1.segmentIndex; i++) {
      pts.add(edge.getCoordinate(i));
      // pts[ipt++] = edge.getCoordinate(i);
    }
    if (useIntPt1) pts.add(new Coordinate.fromCoordinate(ei1.coord));
    // if (useIntPt1) pts[ipt] = new Coordinate.fromCoordinate(ei1.coord);
    return pts;
  }

  /**
   * Gets the list of coordinates for the fully noded segment string,
   * including all original segment string vertices and vertices
   * introduced by nodes in this list.
   * Repeated coordinates are collapsed.
   *
   * @return an array of Coordinates
   *
   */
  List<Coordinate> getSplitCoordinates() {
    CoordinateList coordList = new CoordinateList();
    // ensure that the list has entries for the first and last point of the edge
    addEndpoints();

    Iterator it = iterator();
    // there should always be at least two entries in the list, since the endpoints are nodes
    SegmentNode eiPrev = it.current;
    while (it.moveNext()) {
      SegmentNode ei = it.current;
      addEdgeCoordinates(eiPrev, ei, coordList);
      eiPrev = ei;
    }
    return coordList.toCoordinateArray();
  }

  void addEdgeCoordinates(
      SegmentNode ei0, SegmentNode ei1, CoordinateList coordList) {
    List<Coordinate> pts = createSplitEdgePts(ei0, ei1);
    coordList.addList(pts, false);
  }

//   void print(PrintStream out)
//  {
//    out.println("Intersections:");
//    for (Iterator it = iterator(); it.hasNext(); ) {
//      SegmentNode ei = (SegmentNode) it.next();
//      ei.print(out);
//    }
//  }
}

// INCOMPLETE!
class NodeVertexIterator // implements Iterator
{
  SegmentNodeList nodeList;
  late NodedSegmentString edge;
  late Iterator nodeIt;
  SegmentNode? currNode = null;
  SegmentNode? nextNode = null;
  int currSegIndex = 0;

  NodeVertexIterator(this.nodeList) {
    edge = nodeList.getEdge();
    nodeIt = nodeList.iterator();
    readNextNode();
  }

  bool hasNext() {
    if (nextNode == null) return false;
    return true;
  }

  Object? next() {
    if (currNode == null) {
      currNode = nextNode;
      currSegIndex = currNode!.segmentIndex;
      readNextNode();
      return currNode;
    }
    // check for trying to read too far
    if (nextNode == null) return null;

    if (nextNode!.segmentIndex == currNode!.segmentIndex) {
      currNode = nextNode;
      currSegIndex = currNode!.segmentIndex;
      readNextNode();
      return currNode;
    }

    if (nextNode!.segmentIndex > currNode!.segmentIndex) {}
    return null;
  }

  void readNextNode() {
    if (nodeIt.moveNext())
      nextNode = nodeIt.current;
    else
      nextNode = null;
  }

  /**
   *  Not implemented.
   *
   *@throws  UnsupportedOperationException  This method is not implemented.
   */
  void remove() {
    throw UnimplementedError("remove");
  }
}

/**
 * Methods for computing and working with octants of the Cartesian plane
 * Octants are numbered as follows:
 * <pre>
 *  \2|1/
 * 3 \|/ 0
 * ---+--
 * 4 /|\ 7
 *  /5|6\
 * </pre>
 * If line segments lie along a coordinate axis, the octant is the lower of the two
 * possible values.
 *
 * @version 1.7
 */
class Octant {
  /**
   * Returns the octant of a directed line segment (specified as x and y
   * displacements, which cannot both be 0).
   */
  static int octant(double dx, double dy) {
    if (dx == 0.0 && dy == 0.0)
      throw ArgumentError("Cannot compute the octant for point ( $dx, $dy )");

    double adx = dx.abs();
    double ady = dy.abs();

    if (dx >= 0) {
      if (dy >= 0) {
        if (adx >= ady)
          return 0;
        else
          return 1;
      } else {
        // dy < 0
        if (adx >= ady)
          return 7;
        else
          return 6;
      }
    } else {
      // dx < 0
      if (dy >= 0) {
        if (adx >= ady)
          return 3;
        else
          return 2;
      } else {
        // dy < 0
        if (adx >= ady)
          return 4;
        else
          return 5;
      }
    }
  }

  /**
   * Returns the octant of a directed line segment from p0 to p1.
   */
  static int octantCoords(Coordinate p0, Coordinate p1) {
    double dx = p1.x - p0.x;
    double dy = p1.y - p0.y;
    if (dx == 0.0 && dy == 0.0)
      throw ArgumentError(
          "Cannot compute the octant for two identical points $p0");
    return octant(dx, dy);
  }

  Octant() {}
}

/**
 * Validates that a collection of {@link SegmentString}s is correctly noded.
 * Throws an appropriate exception if an noding error is found.
 *
 * @version 1.7
 */
class NodingValidator {
  LineIntersector li = new RobustLineIntersector();

  List segStrings;

  static final GeometryFactory fact = new GeometryFactory.defaultPrecision();

  NodingValidator(this.segStrings);

  void checkValid() {
    // MD - is this call required?  Or could it be done in the Interior Intersection code?
    checkEndPtVertexIntersections();
    checkInteriorIntersections();
    checkCollapses();
  }

  /**
   * Checks if a segment string contains a segment pattern a-b-a (which implies a self-intersection)
   */
  void checkCollapses() {
    for (SegmentString ss in segStrings) {
      checkCollapsesSS(ss);
    }
  }

  void checkCollapsesSS(SegmentString ss) {
    List<Coordinate> pts = ss.getCoordinates();
    for (int i = 0; i < pts.length - 2; i++) {
      checkCollapse(pts[i], pts[i + 1], pts[i + 2]);
    }
  }

  void checkCollapse(Coordinate p0, Coordinate p1, Coordinate p2) {
    if (p0.equals(p2))
      throw new RuntimeException("found non-noded collapse at " +
          fact.createLineString([p0, p1, p2]).toText());
  }

  /**
   * Checks all pairs of segments for intersections at an interior point of a segment
   */
  void checkInteriorIntersections() {
    for (SegmentString ss0 in segStrings) {
      for (SegmentString ss1 in segStrings) {
        checkInteriorIntersectionsSS(ss0, ss1);
      }
    }
  }

  void checkInteriorIntersectionsSS(SegmentString ss0, SegmentString ss1) {
    List<Coordinate> pts0 = ss0.getCoordinates();
    List<Coordinate> pts1 = ss1.getCoordinates();
    for (int i0 = 0; i0 < pts0.length - 1; i0++) {
      for (int i1 = 0; i1 < pts1.length - 1; i1++) {
        checkInteriorIntersectionsWithIndex(ss0, i0, ss1, i1);
      }
    }
  }

  void checkInteriorIntersectionsWithIndex(
      SegmentString e0, int segIndex0, SegmentString e1, int segIndex1) {
    if (e0 == e1 && segIndex0 == segIndex1) return;
//numTests++;
    Coordinate p00 = e0.getCoordinates()[segIndex0];
    Coordinate p01 = e0.getCoordinates()[segIndex0 + 1];
    Coordinate p10 = e1.getCoordinates()[segIndex1];
    Coordinate p11 = e1.getCoordinates()[segIndex1 + 1];

    li.computeIntersection(p00, p01, p10, p11);
    if (li.hasIntersection()) {
      if (li.isProper() ||
          hasInteriorIntersection(li, p00, p01) ||
          hasInteriorIntersection(li, p10, p11)) {
        throw new RuntimeException("found non-noded intersection at " +
            p00.toString() +
            "-" +
            p01.toString() +
            " and " +
            p10.toString() +
            "-" +
            p11.toString());
      }
    }
  }

  /**
   *@return true if there is an intersection point which is not an endpoint of the segment p0-p1
   */
  bool hasInteriorIntersection(
      LineIntersector li, Coordinate p0, Coordinate p1) {
    for (int i = 0; i < li.getIntersectionNum(); i++) {
      Coordinate intPt = li.getIntersection(i);
      if (!(intPt.equals(p0) || intPt.equals(p1))) return true;
    }
    return false;
  }

  /**
   * Checks for intersections between an endpoint of a segment string
   * and an interior vertex of another segment string
   */
  void checkEndPtVertexIntersections() {
    for (SegmentString ss in segStrings) {
      List<Coordinate> pts = ss.getCoordinates();
      checkEndPtVertexIntersections2(pts[0], segStrings);
      checkEndPtVertexIntersections2(pts[pts.length - 1], segStrings);
    }
  }

  void checkEndPtVertexIntersections2(Coordinate testPt, List segStrings) {
    for (SegmentString ss in segStrings) {
      List<Coordinate> pts = ss.getCoordinates();
      for (int j = 1; j < pts.length - 1; j++) {
        if (pts[j].equals(testPt))
          throw new RuntimeException(
              "found endpt/interior pt intersection at index $j :pt " +
                  testPt.toString());
      }
    }
  }
}

/**
 * Finds <b>interior</b> intersections between line segments in {@link NodedSegmentString}s,
 * and adds them as nodes
 * using {@link NodedSegmentString#addIntersection(LineIntersector, int, int, int)}.
 * <p>
 * This class is used primarily for Snap-Rounding.
 * For general-purpose noding, use {@link IntersectionAdder}.
 *
 * @version 1.7
 * @see IntersectionAdder
 */
class InteriorIntersectionFinderAdder implements SegmentIntersectorN {
  LineIntersector li;
  late List interiorIntersections;

  /**
   * Creates an intersection finder which finds all proper intersections
   *
   * @param li the LineIntersector to use
   */
  InteriorIntersectionFinderAdder(this.li) {
    interiorIntersections = [];
  }

  List getInteriorIntersections() {
    return interiorIntersections;
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

    Coordinate p00 = e0.getCoordinates()[segIndex0];
    Coordinate p01 = e0.getCoordinates()[segIndex0 + 1];
    Coordinate p10 = e1.getCoordinates()[segIndex1];
    Coordinate p11 = e1.getCoordinates()[segIndex1 + 1];

    li.computeIntersection(p00, p01, p10, p11);
//if (li.hasIntersection() && li.isProper()) Debug.println(li);

    if (li.hasIntersection()) {
      if (li.isInteriorIntersection()) {
        for (int intIndex = 0; intIndex < li.getIntersectionNum(); intIndex++) {
          interiorIntersections.add(li.getIntersection(intIndex));
        }
        (e0 as NodedSegmentString).addIntersections(li, segIndex0, 0);
        (e1 as NodedSegmentString).addIntersections(li, segIndex1, 1);
      }
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

/**
 * Wraps a {@link Noder} and transforms its input
 * into the integer domain.
 * This is intended for use with Snap-Rounding noders,
 * which typically are only intended to work in the integer domain.
 * Offsets can be provided to increase the number of digits of available precision.
 * <p>
 * Clients should be aware that rescaling can involve loss of precision,
 * which can cause zero-length line segments to be created.
 * These in turn can cause problems when used to build a planar graph.
 * This situation should be checked for and collapsed segments removed if necessary.
 *
 * @version 1.7
 */
class ScaledNoder implements Noder {
  Noder noder;
  double scaleFactor = 0;
  double offsetX = 0;
  double offsetY = 0;
  bool isScaled = false;

  ScaledNoder(Noder noder, double scaleFactor)
      : this.withOffests(noder, scaleFactor, 0, 0);

  ScaledNoder.withOffests(
      this.noder, this.scaleFactor, double offsetX, double offsetY) {
    // no need to scale if input precision is already integral
    isScaled = !isIntegerPrecision();
  }

  bool isIntegerPrecision() {
    return scaleFactor == 1.0;
  }

  List getNodedSubstrings() {
    List splitSS = noder.getNodedSubstrings();
    if (isScaled) rescale(splitSS);
    return splitSS;
  }

  void computeNodes(List inputSegStrings) {
    List intSegStrings = inputSegStrings;
    if (isScaled) intSegStrings = scale(inputSegStrings);
    noder.computeNodes(intSegStrings);
  }

  List scale(List segStrings) {
    List<NodedSegmentString> nodedSegmentStrings = [];
    for (SegmentString ss in segStrings) {
      nodedSegmentStrings.add(new NodedSegmentString(
          scaleCoords(ss.getCoordinates()), ss.getData()));
    }
    return nodedSegmentStrings;
  }

  List<Coordinate> scaleCoords(List<Coordinate> pts) {
    List<Coordinate> roundPts = []; //..length = (pts.length);
    for (int i = 0; i < pts.length; i++) {
      // roundPts[i] = Coordinate.fromXYZ(
      roundPts.add(Coordinate.fromXYZ(
        ((pts[i].x - offsetX) * scaleFactor).roundToDouble(),
        ((pts[i].y - offsetY) * scaleFactor).roundToDouble(),
        pts[i].getZ(),
      ));
    }
    List<Coordinate> roundPtsNoDup =
        CoordinateArrays.removeRepeatedPoints(roundPts);
    return roundPtsNoDup;
  }

  // double scale(double val) { return (double) Math.round(val * scaleFactor); }

  void rescale(List segStrings) {
    for (SegmentString ss in segStrings) {
      rescaleCoords(ss.getCoordinates());
    }
  }

  void rescaleCoords(List<Coordinate> pts) {
    for (int i = 0; i < pts.length; i++) {
      pts[i].x = pts[i].x / scaleFactor + offsetX;
      pts[i].y = pts[i].y / scaleFactor + offsetY;
    }

//    if (pts.length == 2 && pts[0].equals2D(pts[1])) {
//      System.out.println(pts);
//    }
  }

// double rescale(double val) { return val / scaleFactor; }
}

/**
 * Computes the possible intersections between two line segments in {@link NodedSegmentString}s
 * and adds them to each string
 * using {@link NodedSegmentString#addIntersection(LineIntersector, int, int, int)}.
 *
 * @version 1.7
 */
class IntersectionAdder implements SegmentIntersectorN {
  static bool isAdjacentSegments(int i1, int i2) {
    return (i1 - i2).abs() == 1;
  }

  /**
   * These variables keep track of what types of intersections were
   * found during ALL edges that have been intersected.
   */
  bool _hasIntersection = false;
  bool hasProper = false;
  bool hasProperInterior = false;
  bool hasInterior = false;

  // the proper intersection point found
  Coordinate? properIntersectionPoint = null;

  LineIntersector li;
  bool isSelfIntersection = false;

  // bool intersectionFound;
  int numIntersections = 0;
  int numInteriorIntersections = 0;
  int numProperIntersections = 0;

  // testing only
  int numTests = 0;

  IntersectionAdder(this.li);

  LineIntersector getLineIntersector() {
    return li;
  }

  /**
   * @return the proper intersection point, or <code>null</code> if none was found
   */
  Coordinate? getProperIntersectionPoint() {
    return properIntersectionPoint;
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
    return hasProper;
  }

  /**
   * A proper interior intersection is a proper intersection which is <b>not</b>
   * contained in the set of boundary nodes set for this SegmentIntersector.
   */
  bool hasProperInteriorIntersection() {
    return hasProperInterior;
  }

  /**
   * An interior intersection is an intersection which is
   * in the interior of some segment.
   */
  bool hasInteriorIntersection() {
    return hasInterior;
  }

  /**
   * A trivial intersection is an apparent self-intersection which in fact
   * is simply the point shared by adjacent line segments.
   * Note that closed edges require a special check for the point shared by the beginning
   * and end segments.
   */
  bool isTrivialIntersection(
      SegmentString e0, int segIndex0, SegmentString e1, int segIndex1) {
    if (e0 == e1) {
      if (li.getIntersectionNum() == 1) {
        if (isAdjacentSegments(segIndex0, segIndex1)) return true;
        if (e0.isClosed()) {
          int maxSegIndex = e0.size() - 1;
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
   * This method is called by clients
   * of the {@link SegmentIntersector} class to process
   * intersections for two segments of the {@link SegmentString}s being intersected.
   * Note that some clients (such as <code>MonotoneChain</code>s) may optimize away
   * this call for segment pairs which they have determined do not intersect
   * (e.g. by an disjoint envelope test).
   */
  void processIntersections(
      SegmentString e0, int segIndex0, SegmentString e1, int segIndex1) {
    if (e0 == e1 && segIndex0 == segIndex1) return;
    numTests++;
    Coordinate p00 = e0.getCoordinates()[segIndex0];
    Coordinate p01 = e0.getCoordinates()[segIndex0 + 1];
    Coordinate p10 = e1.getCoordinates()[segIndex1];
    Coordinate p11 = e1.getCoordinates()[segIndex1 + 1];

    li.computeIntersection(p00, p01, p10, p11);
//if (li.hasIntersection() && li.isProper()) Debug.println(li);
    if (li.hasIntersection()) {
      //intersectionFound = true;
      numIntersections++;
      if (li.isInteriorIntersection()) {
        numInteriorIntersections++;
        hasInterior = true;
//System.out.println(li);
      }
      // if the segments are adjacent they have at least one trivial intersection,
      // the shared endpoint.  Don't bother adding it if it is the
      // only intersection.
      if (!isTrivialIntersection(e0, segIndex0, e1, segIndex1)) {
        _hasIntersection = true;
        (e0 as NodedSegmentString).addIntersections(li, segIndex0, 0);
        (e1 as NodedSegmentString).addIntersections(li, segIndex1, 1);
        if (li.isProper()) {
          numProperIntersections++;
//Debug.println(li.toString());  Debug.println(li.getIntersection(0));
          //properIntersectionPoint = (Coordinate) li.getIntersection(0).clone();
          hasProper = true;
          hasProperInterior = true;
        }
      }
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
