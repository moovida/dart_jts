part of dart_jts;

/**
 * Builds a set of noded, unique, labelled Edges from
 * the edges of the two input geometries.
 * <p>
 * It performs the following steps:
 * <ul>
 * <li>Extracts input edges, and attaches topological information
 * <li>if clipping is enabled, handles clipping or limiting input geometry
 * <li>chooses a {@link Noder} based on provided precision model, unless a custom one is supplied
 * <li>calls the chosen Noder, with precision model
 * <li>removes any fully collapsed noded edges
 * <li>builds {@link EdgeNG}s and merges them
 * </ul>
 *
 * @author mdavis
 *
 */
class EdgeNodingBuilder {
  /**
   * Limiting is skipped for Lines with few vertices,
   * to avoid additional copying.
   */
  static const int MIN_LIMIT_PTS = 20;

  /**
   * Indicates whether floating precision noder output is validated.
   */
  static const bool IS_NODING_VALIDATED = true;

  static Noder _createFixedPrecisionNoder(PrecisionModel pm) {
    Noder noder = SnapRoundingNoder(pm);
    return noder;
  }

  static Noder _createFloatingPrecisionNoder(bool doValidation) {
    MCIndexNoder mcNoder = MCIndexNoder.empty();
    LineIntersector li = RobustLineIntersector();
    mcNoder.setSegmentIntersector(IntersectionAdder(li));

    Noder noder = mcNoder;
    if (doValidation) {
      noder = ValidatingNoder(mcNoder);
    }
    return noder;
  }

  PrecisionModel _pm;
  List<NodedSegmentString> _inputEdges = List.empty(growable: true);
  Noder? _customNoder;

  Envelope? _clipEnv = null;
  RingClipper? _clipper;
  LineLimiter? _limiter;

  List<bool> hasEdges = List.filled(2, false);

  /**
   * Creates a new builder, with an optional custom noder.
   * If the noder is not provided, a suitable one will
   * be used based on the supplied precision model.
   *
   * @param pm the precision model to use
   * @param noder an optional custom noder to use (may be null)
   */
  EdgeNodingBuilder(this._pm, this._customNoder);

  /**
   * Gets a noder appropriate for the precision model supplied.
   * This is one of:
   * <ul>
   * <li>Fixed precision: a snap-rounding noder (which should be fully robust)
   * <li>Floating precision: a conventional nodel (which may be non-robust).
   * In this case, a validation step is applied to the output from the noder.
   * </ul>
   *
   * @return
   */
  Noder _getNoder() {
    if (_customNoder != null) return _customNoder!;
    if (OverlayUtil.isFloating(_pm))
      return _createFloatingPrecisionNoder(IS_NODING_VALIDATED);
    return _createFixedPrecisionNoder(_pm);
  }

  void setClipEnvelope(Envelope clipEnv) {
    _clipEnv = clipEnv;
    _clipper = RingClipper(clipEnv);
    _limiter = LineLimiter(clipEnv);
  }

  /**
   * Reports whether there are noded edges
   * for the given input geometry.
   * If there are none, this indicates that either
   * the geometry was empty, or has completely collapsed
   * (because it is smaller than the noding precision).
   *
   * @param geomIndex index of input geometry
   * @return true if there are edges for the geometry
   */
  bool hasEdgesFor(int geomIndex) {
    return hasEdges[geomIndex];
  }

  /**
   * Creates a set of labelled {EdgeNG}s.
   * representing the fully noded edges of the input geometries.
   * Coincident edges (from the same or both geometries)
   * are merged along with their labels
   * into a single unique, fully labelled edge.
   *
   * @param geom0 the first geometry
   * @param geom1 the second geometry
   * @return the noded, merged, labelled edges
   */
  List<EdgeNG> build(Geometry? geom0, Geometry? geom1) {
    _add(geom0, 0);
    _add(geom1, 1);
    List<EdgeNG> nodedEdges = _node(_inputEdges);

    /**
     * Merge the noded edges to eliminate duplicates.
     * Labels are combined.
     */
    List<EdgeNG> mergedEdges = EdgeMerger.merge(nodedEdges);
    return mergedEdges;
  }

  /**
   * Nodes a set of segment strings and creates {@link EdgeNG}s from the result.
   * The input segment strings each carry a {@link EdgeSourceInfo} object,
   * which is used to provide source topology info to the constructed Edges
   * (and then is discarded).
   *
   * @param segStrings
   * @return
   */
  List<EdgeNG> _node(List<NodedSegmentString> segStrings) {
    Noder noder = _getNoder();
    noder.computeNodes(segStrings);
    var list = noder.getNodedSubstrings();
    List<SegmentString> nodedSS = List.empty(growable: true);
    for (var e in list) {
      nodedSS.add(e as SegmentString);
    }
    List<EdgeNG> edges = _createEdges(nodedSS);
    return edges;
  }

  List<EdgeNG> _createEdges(List<SegmentString> segStrings) {
    List<EdgeNG> edges = List.empty(growable: true);
    for (SegmentString ss in segStrings) {
      List<Coordinate> pts = ss.getCoordinates();

      //-- don't create edges from collapsed lines
      if (EdgeNG.isCollapsed(pts)) continue;

      EdgeSourceInfo info = ss.getData() as EdgeSourceInfo;
      //-- Record that a non-collapsed edge exists for the parent geometry
      hasEdges[info.getIndex()] = true;
      edges.add(EdgeNG(ss.getCoordinates(), info));
    }
    return edges;
  }

  void _add(Geometry? g, int geomIndex) {
    if (g == null || g.isEmpty()) return;

    if (_isClippedCompletely(g.getEnvelopeInternal())) return;

    if (g is Polygon)
      _addPolygon(g, geomIndex);
    // LineString also handles LinearRings
    else if (g is LineString)
      _addLine(g, geomIndex);
    else if (g is MultiLineString)
      _addCollection(g, geomIndex);
    else if (g is MultiPolygon)
      _addCollection(g, geomIndex);
    else if (g is GeometryCollection)
      _addGeometryCollection(g, geomIndex, g.getDimension());
    // ignore Point geometries - they are handled elsewhere
  }

  void _addCollection(GeometryCollection gc, int geomIndex) {
    for (int i = 0; i < gc.getNumGeometries(); i++) {
      Geometry g = gc.getGeometryN(i);
      _add(g, geomIndex);
    }
  }

  void _addGeometryCollection(
      GeometryCollection gc, int geomIndex, int expectedDim) {
    for (int i = 0; i < gc.getNumGeometries(); i++) {
      Geometry g = gc.getGeometryN(i);
      // check for mixed-dimension input, which is not supported
      if (g.getDimension() != expectedDim) {
        throw Exception("Overlay input is mixed-dimension");
      }
      _add(g, geomIndex);
    }
  }

  void _addPolygon(Polygon poly, int geomIndex) {
    LinearRing shell = poly.getExteriorRing();
    _addPolygonRing(shell, false, geomIndex);

    for (int i = 0; i < poly.getNumInteriorRing(); i++) {
      LinearRing hole = poly.getInteriorRingN(i);

      // Holes are topologically labelled opposite to the shell, since
      // the interior of the polygon lies on their opposite side
      // (on the left, if the hole is oriented CW)
      _addPolygonRing(hole, true, geomIndex);
    }
  }

  /**
   * Adds a polygon ring to the graph.
   * Empty rings are ignored.
   */
  void _addPolygonRing(LinearRing ring, bool isHole, int index) {
    // don't add empty rings
    if (ring.isEmpty()) return;

    if (_isClippedCompletely(ring.getEnvelopeInternal())) return;

    List<Coordinate> pts = _clip(ring);

    /**
     * Don't add edges that collapse to a point
     */
    if (pts.length < 2) {
      return;
    }

    //if (pts.length < ring.getNumPoints()) System.out.println("Ring clipped: " + ring.getNumPoints() + " => " + pts.length);

    int depthDelta = _computeDepthDelta(ring, isHole);
    EdgeSourceInfo info = EdgeSourceInfo(index, depthDelta, isHole);
    _addEdge(pts, info);
  }

  /**
   * Tests whether a geometry (represented by its envelope)
   * lies completely outside the clip extent(if any).
   *
   * @param env the geometry envelope
   * @return true if the geometry envelope is outside the clip extent.
   */
  bool _isClippedCompletely(Envelope env) {
    if (_clipEnv == null) return false;
    return _clipEnv!.disjoint(env);
  }

  /**
   * If a clipper is present,
   * clip the line to the clip extent.
   * Otherwise, remove duplicate points from the ring.
   * <p>
   * If clipping is enabled, then every ring MUST
   * be clipped, to ensure that holes are clipped to
   * be inside the shell.
   * This means it is not possible to skip
   * clipping for rings with few vertices.
   *
   * @param ring the line to clip
   * @return the points in the clipped line
   */
  List<Coordinate> _clip(LinearRing ring) {
    List<Coordinate> pts = ring.getCoordinates();
    Envelope env = ring.getEnvelopeInternal();

    /**
     * If no clipper or ring is completely contained then no need to clip.
     * But repeated points must be removed to ensure correct noding.
     */
    if (_clipper == null || _clipEnv!.coversEnvelope(env)) {
      return _removeRepeatedPoints(ring);
    }

    return _clipper!.clip(pts);
  }

  /**
   * Removes any repeated points from a linear component.
   * This is required so that noding can be computed correctly.
   *
   * @param line the line to process
   * @return the points of the line with repeated points removed
   */
  static List<Coordinate> _removeRepeatedPoints(LineString line) {
    List<Coordinate> pts = line.getCoordinates();
    return CoordinateArrays.removeRepeatedPoints(pts);
  }

  static int _computeDepthDelta(LinearRing ring, bool isHole) {
    /**
     * Compute the orientation of the ring, to
     * allow assigning side interior/exterior labels correctly.
     * JTS canonical orientation is that shells are CW, holes are CCW.
     *
     * It is important to compute orientation on the original ring,
     * since topology collapse can make the orientation computation give the wrong answer.
     */
    bool isCCW =
        Orientation.isCCW(ring.getCoordinateSequence().toCoordinateArray());
    /**
     * Compute whether ring is in canonical orientation or not.
     * Canonical orientation for the overlay process is
     * Shells : CW, Holes: CCW
     */
    bool isOriented = true;
    if (!isHole)
      isOriented = !isCCW;
    else {
      isOriented = isCCW;
    }
    /**
     * Depth delta can now be computed.
     * Canonical depth delta is 1 (Exterior on L, Interior on R).
     * It is flipped to -1 if the ring is oppositely oriented.
     */
    int depthDelta = isOriented ? 1 : -1;
    return depthDelta;
  }

  /**
   * Adds a line geometry, limiting it if enabled,
   * and otherwise removing repeated points.
   *
   * @param line the line to add
   * @param geomIndex the index of the parent geometry
   */
  void _addLine(LineString line, int geomIndex) {
    // don't add empty lines
    if (line.isEmpty()) return;

    if (_isClippedCompletely(line.getEnvelopeInternal())) return;

    if (_isToBeLimited(line)) {
      List<List<Coordinate>> sections = _limit(line);
      for (List<Coordinate> pts in sections) {
        _addLineList(pts, geomIndex);
      }
    } else {
      List<Coordinate> ptsNoRepeat = _removeRepeatedPoints(line);
      _addLineList(ptsNoRepeat, geomIndex);
    }
  }

  void _addLineList(List<Coordinate> pts, int geomIndex) {
    /**
     * Don't add edges that collapse to a point
     */
    if (pts.length < 2) {
      return;
    }

    EdgeSourceInfo info = EdgeSourceInfo.fromIndex(geomIndex);
    _addEdge(pts, info);
  }

  void _addEdge(List<Coordinate> pts, EdgeSourceInfo info) {
    NodedSegmentString ss = NodedSegmentString(pts, info);
    _inputEdges.add(ss);
  }

  /**
   * Tests whether it is worth limiting a line.
   * Lines that have few vertices or are covered
   * by the clip extent do not need to be limited.
   *
   * @param line line to test
   * @return true if the line should be limited
   */
  bool _isToBeLimited(LineString line) {
    List<Coordinate> pts = line.getCoordinates();
    if (_limiter == null || pts.length <= MIN_LIMIT_PTS) {
      return false;
    }
    Envelope env = line.getEnvelopeInternal();
    /**
     * If line is completely contained then no need to limit
     */
    if (_clipEnv!.coversEnvelope(env)) {
      return false;
    }
    return true;
  }

  /**
   * If limiter is provided,
   * limit the line to the clip envelope.
   *
   * @param line the line to clip
   * @return the point sections in the clipped line
   */
  List<List<Coordinate>> _limit(LineString line) {
    List<Coordinate> pts = line.getCoordinates();
    return _limiter!.limit(pts);
  }
}
