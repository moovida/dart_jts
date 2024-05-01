part of dart_jts;

/**
 * Computes the geometric overlay of two {@link Geometry}s,
 * using an explicit precision model to allow robust computation.
 * <p>
 * The overlay can be used to determine any of the
 * following set-theoretic operations (bool combinations) of the geometries:</p>
 * <ul>
 * <li>{@link #INTERSECTION} - all points which lie in both geometries</li>
 * <li>{@link #UNION} - all points which lie in at least one geometry</li>
 * <li>{@link #DIFFERENCE} - all points which lie in the first geometry but not the second</li>
 * <li>{@link #SYMDIFFERENCE} - all points which lie in one geometry but not both</li>
 * </ul>
 * Input geometries may have different dimension.
 * Input collections must be homogeneous (all elements must have the same dimension).
 * Inputs may be <b>simple</b> {@link GeometryCollection}s.
 * A GeometryCollection is simple if it can be flattened into a valid Multi-geometry;
 * i.e. it is homogeneous and does not contain any overlapping Polygons.
 * <p>
 * The precision model used for the computation can be supplied
 * independent of the precision model of the input geometry.
 * The main use for this is to allow using a fixed precision
 * for geometry with a floating precision model.
 * This does two things: ensures robust computation;
 * and forces the output to be validly rounded to the precision model.</p>
 * <p>
 * For fixed precision models noding is performed using a {@link SnapRoundingNoder}.
 * This provides robust computation (as long as precision is limited to
 * around 13 decimal digits).</p>
 * <p>
 * For floating precision an {@link MCIndexNoder} is used.
 * This is not fully robust, so can sometimes result in
 * {@link TopologyException}s being thrown.
 * For robust full-precision overlay see {@link OverlayNGRobust}.</p>
 * <p>
 * A custom {@link Noder} can be supplied.
 * This allows using a more performant noding strategy in specific cases,
 * for instance in {@link CoverageUnion}.</p>
 * <p>
 * <b>Note:</b If a {@link SnappingNoder} is used
 * it is best to specify a fairly small snap tolerance,
 * since the intersection clipping optimization can
 * interact with the snapping to alter the result.</p>
 * <p>
 * Optionally the overlay computation can process using strict mode
 * (via {@link #setStrictMode(bool)}.
 * In strict mode result semantics are:</p>
 * <ul>
 * <li>Lines and Points resulting from topology collapses are not included in the result</li>
 * <li>Result geometry is homogeneous
 *     for the {@link #INTERSECTION} and {@link #DIFFERENCE} operations.</li>
 * <li>Result geometry is homogeneous
 *     for the {@link #UNION} and {@link #SYMDIFFERENCE} operations
 *     if the inputs have the same dimension</li>
 * </ul>
 * <p>
 * Strict mode has the following benefits:</p>
 * <ul>
 * <li>Results are simpler</li>
 * <li>Overlay operations are chainable
 *     without needing to remove lower-dimension elements</li>
 * </ul>
 * <p>
 * The original JTS overlay semantics corresponds to non-strict mode.</p>
 * <p>
 * If a robustness error occurs, a {@link TopologyException} is thrown.
 * These are usually caused by numerical rounding causing the noding output
 * to not be fully noded.
 * For robust computation with full-precision {@link OverlayNGRobust} can be used.</p>
 *
 * @author mdavis
 *
 * @see OverlayNGRobust
 *
 */
class OverlayNG {
  /**
   * The code for the Intersection overlay operation.
   */
  static const INTERSECTION = OverlayOp.INTERSECTION;

  /**
   * The code for the Union overlay operation.
   */
  static const UNION = OverlayOp.UNION;

  /**
   *  The code for the Difference overlay operation.
   */
  static const DIFFERENCE = OverlayOp.DIFFERENCE;

  /**
   *  The code for the Symmetric Difference overlay operation.
   */
  static const SYMDIFFERENCE = OverlayOp.SYMDIFFERENCE;

  /**
   * The default setting for Strict Mode.
   *
   * The original JTS overlay semantics used non-strict result
   * semantics, including;
   * - An Intersection result can be mixed-dimension,
   *   due to inclusion of intersection components of all dimensions
   * - Results can include lines caused by Area topology collapse
   */
  static const STRICT_MODE_DEFAULT = false;

  /**
   * Tests whether a point with a given topological {@link Label}
   * relative to two geometries is contained in
   * the result of overlaying the geometries using
   * a given overlay operation.
   * <p>
   * The method handles arguments of {@link Location#NONE} correctly
   *
   * @param label the topological label of the point
   * @param opCode the code for the overlay operation to test
   * @return true if the label locations correspond to the overlayOpCode
   */
  static bool isResultOfOpPoint(OverlayLabel label, int opCode) {
    int loc0 = label.getLineLocation(0);
    int loc1 = label.getLineLocation(1);
    return isResultOfOp(opCode, loc0, loc1);
  }

  /**
   * Tests whether a point with given {@link Location}s
   * relative to two geometries would be contained in
   * the result of overlaying the geometries using
   * a given overlay operation.
   * This is used to determine whether components
   * computed during the overlay process should be
   * included in the result geometry.
   * <p>
   * The method handles arguments of {@link Location#NONE} correctly.
   *
   * @param overlayOpCode the code for the overlay operation to test
   * @param loc0 the code for the location in the first geometry
   * @param loc1 the code for the location in the second geometry
   *
   * @return true if a point with given locations is in the result of the overlay operation
   */
  static bool isResultOfOp(int overlayOpCode, int loc0, int loc1) {
    if (loc0 == Location.BOUNDARY) loc0 = Location.INTERIOR;
    if (loc1 == Location.BOUNDARY) loc1 = Location.INTERIOR;
    switch (overlayOpCode) {
      case OverlayOp.INTERSECTION:
        return loc0 == Location.INTERIOR && loc1 == Location.INTERIOR;
      case UNION:
        return loc0 == Location.INTERIOR || loc1 == Location.INTERIOR;
      case DIFFERENCE:
        return loc0 == Location.INTERIOR && loc1 != Location.INTERIOR;
      case SYMDIFFERENCE:
        return (loc0 == Location.INTERIOR && loc1 != Location.INTERIOR) ||
            (loc0 != Location.INTERIOR && loc1 == Location.INTERIOR);
    }
    return false;
  }

  /**
   * Computes an overlay operation for
   * the given geometry operands, with the
   * noding strategy determined by the precision model.
   *
   * @param geom0 the first geometry argument
   * @param geom1 the second geometry argument
   * @param opCode the code for the desired overlay operation
   * @param pm the precision model to use
   * @return the result of the overlay operation
   */
  static Geometry? overlayWithPM(
      Geometry geom0, Geometry geom1, int opCode, PrecisionModel pm) {
    OverlayNG ov = OverlayNG.withPM(geom0, geom1, pm, opCode);
    Geometry? geomOv = ov.getResult();
    return geomOv;
  }

  /**
   * Computes an overlay operation on the given geometry operands,
   * using a supplied {@link Noder}.
   *
   * @param geom0 the first geometry argument
   * @param geom1 the second geometry argument
   * @param opCode the code for the desired overlay operation
   * @param pm the precision model to use (which may be null if the noder does not use one)
   * @param noder the noder to use
   * @return the result of the overlay operation
   */
  static Geometry? overlayWithNoder(Geometry geom0, Geometry geom1, int opCode,
      PrecisionModel pm, Noder noder) {
    OverlayNG ov = OverlayNG.withPM(geom0, geom1, pm, opCode);
    ov.setNoder(noder);
    Geometry? geomOv = ov.getResult();
    return geomOv;
  }

  /**
   * Computes an overlay operation on the given geometry operands,
   * using a supplied {@link Noder}.
   *
   * @param geom0 the first geometry argument
   * @param geom1 the second geometry argument
   * @param opCode the code for the desired overlay operation
   * @param noder the noder to use
   * @return the result of the overlay operation
   */
  static Geometry? overlayWithPMAndNoder(
      Geometry geom0, Geometry geom1, int opCode, Noder noder) {
    OverlayNG ov = OverlayNG(geom0, geom1, opCode);
    ov.setNoder(noder);
    Geometry? geomOv = ov.getResult();
    return geomOv;
  }

  /**
   * Computes an overlay operation on
   * the given geometry operands,
   * using the precision model of the geometry.
   * and an appropriate noder.
   * <p>
   * The noder is chosen according to the precision model specified.
   * <ul>
   * <li>For {@link PrecisionModel#FIXED}
   * a snap-rounding noder is used, and the computation is robust.
   * <li>For {@link PrecisionModel#FLOATING}
   * a non-snapping noder is used,
   * and this computation may not be robust.
   * If errors occur a {@link TopologyException} is thrown.
   * </ul>
   *
   *
   * @param geom0 the first argument geometry
   * @param geom1 the second argument geometry
   * @param opCode the code for the desired overlay operation
   * @return the result of the overlay operation
   */
  static Geometry? overlay(Geometry geom0, Geometry geom1, int opCode) {
    OverlayNG ov = OverlayNG(geom0, geom1, opCode);
    return ov.getResult();
  }

  /**
   * Computes a union operation on
   * the given geometry, with the supplied precision model.
   * <p>
   * The input must be a valid geometry.
   * Collections must be homogeneous.
   * <p>
   * To union an overlapping set of polygons in a more performant way use {@link UnaryUnionNG}.
   * To union a polyonal coverage or linear network in a more performant way,
   * use {@link CoverageUnion}.
   *
   * @param geom0 the geometry
   * @param pm the precision model to use
   * @return the result of the union operation
   *
   * @see OverlayMixedPoints
   */
  static Geometry? union(Geometry geom, PrecisionModel pm) {
    OverlayNG ov = OverlayNG.fromSingle(geom, pm);
    Geometry? geomOv = ov.getResult();
    return geomOv;
  }

  /**
   * Computes a union of a single geometry using a custom noder.
   * <p>
   * The primary use of this is to support coverage union.
   * Because of this the overlay is performed using strict mode.
   *
   * @param geom the geometry to union
   * @param pm the precision model to use (maybe be null)
   * @param noder the noder to use
   * @return the result geometry
   *
   * @see CoverageUnion
   */
  static Geometry? unionWithNoder(
      Geometry geom, PrecisionModel pm, Noder noder) {
    OverlayNG? ov = OverlayNG.fromSingle(geom, pm);
    ov.setNoder(noder);
    ov.setStrictMode(true);
    Geometry? geomOv = ov.getResult();
    return geomOv;
  }

  int? opCode;
  late InputGeometry inputGeom;
  late GeometryFactory geomFact;
  PrecisionModel? pm;
  Noder? noder;
  bool isStrictMode = STRICT_MODE_DEFAULT;
  bool isOptimized = true;
  bool isAreaResultOnly = false;
  bool isOutputEdges = false;
  bool isOutputResultEdges = false;
  bool isOutputNodedEdges = false;

  /**
   * Creates an overlay operation on the given geometries,
   * with a defined precision model.
   * The noding strategy is determined by the precision model.
   *
   * @param geom0 the A operand geometry
   * @param geom1 the B operand geometry (may be null)
   * @param pm the precision model to use
   * @param opCode the overlay opcode
   */
  OverlayNG.withPM(Geometry geom0, Geometry? geom1, this.pm, this.opCode) {
    geomFact = geom0.getFactory();
    inputGeom = InputGeometry(geom0, geom1);
  }

  /**
   * Creates an overlay operation on the given geometries
   * using the precision model of the geometries.
   * <p>
   * The noder is chosen according to the precision model specified.
   * <ul>
   * <li>For {@link PrecisionModel#FIXED}
   * a snap-rounding noder is used, and the computation is robust.
   * <li>For {@link PrecisionModel#FLOATING}
   * a non-snapping noder is used,
   * and this computation may not be robust.
   * If errors occur a {@link TopologyException} is thrown.
   * </ul>
   *
   * @param geom0 the A operand geometry
   * @param geom1 the B operand geometry (may be null)
   * @param opCode the overlay opcode
   */
  factory OverlayNG(Geometry geom0, Geometry? geom1, int opCode) {
    return OverlayNG.withPM(
        geom0, geom1!, geom0.getFactory().getPrecisionModel(), opCode);
  }

  /**
   * Creates a union of a single geometry with a given precision model.
   *
   * @param geom the geometry
   * @param pm the precision model to use
   */
  factory OverlayNG.fromSingle(Geometry geom, PrecisionModel pm) {
    return OverlayNG.withPM(geom, null, pm, UNION);
  }

  /**
   * Sets whether the overlay results are computed according to strict mode
   * semantics.
   * <ul>
   * <li>Lines resulting from topology collapse are not included
   * <li>Result geometry is homogeneous
   *     for the {@link #INTERSECTION} and {@link #DIFFERENCE} operations.
   * <li>Result geometry is homogeneous
   *     for the {@link #UNION} and {@link #SYMDIFFERENCE} operations
   *     if the inputs have the same dimension
   * </ul>
   *
   * @param isStrictMode true if strict mode is to be used
   */
  void setStrictMode(bool isStrictMode) {
    this.isStrictMode = isStrictMode;
  }

  /**
   * Sets whether overlay processing optimizations are enabled.
   * It may be useful to disable optimizations
   * for testing purposes.
   * Default is TRUE (optimization enabled).
   *
   * @param isOptimized whether to optimize processing
   */
  void setOptimized(bool isOptimized) {
    this.isOptimized = isOptimized;
  }

  /**
   * Sets whether the result can contain only {@link Polygon} components.
   * This is used if it is known that the result must be an (possibly empty) area.
   *
   * @param isAreaResultOnly true if the result should contain only area components
   */
  void setAreaResultOnly(bool isAreaResultOnly) {
    this.isAreaResultOnly = isAreaResultOnly;
  }

  //------ Testing options -------

  /**
   *
   * @param isOutputEdges
   */
  void setOutputEdges(bool isOutputEdges) {
    this.isOutputEdges = isOutputEdges;
  }

  void setOutputNodedEdges(bool isOutputNodedEdges) {
    this.isOutputEdges = true;
    this.isOutputNodedEdges = isOutputNodedEdges;
  }

  void setOutputResultEdges(bool isOutputResultEdges) {
    this.isOutputResultEdges = isOutputResultEdges;
  }
  //---------------------------------

  void setNoder(Noder noder) {
    this.noder = noder;
  }

  /**
   * Gets the result of the overlay operation.
   *
   * @return the result of the overlay operation.
   *
   * @throws IllegalArgumentException if the input is not supported (e.g. a mixed-dimension geometry)
   * @throws TopologyException if a robustness error occurs
   */
  Geometry? getResult() {
    // handle empty inputs which determine result
    if (OverlayUtil.isEmptyResult(
        opCode!, inputGeom.getGeometry(0), inputGeom.getGeometry(1), pm!)) {
      return createEmptyResult();
    }

    /**
     * The elevation model is only computed if the input geometries have Z values.
     */
    ElevationModel elevModel = ElevationModel.create(
        inputGeom.getGeometry(0)!, inputGeom.getGeometry(1));
    Geometry? result;
    if (inputGeom.isAllPoints()) {
      // handle Point-Point inputs
      result = OverlayPoints.overlay(
          opCode!, inputGeom.getGeometry(0)!, inputGeom.getGeometry(1)!, pm!);
    } else if (!inputGeom.isSingle() && inputGeom.hasPoints()) {
      // handle Point-nonPoint inputs
      result = OverlayMixedPoints.overlay(
          opCode!, inputGeom.getGeometry(0)!, inputGeom.getGeometry(1)!, pm!);
    } else {
      // handle case where both inputs are formed of edges (Lines and Polygons)
      result = computeEdgeOverlay();
    }
    /**
     * This is a no-op if the elevation model was not computed due to Z not present
     */
    elevModel.populateZ(result!);
    return result;
  }

  Geometry? computeEdgeOverlay() {
    List<EdgeNG> edges = nodeEdges();

    OverlayGraph graph = buildGraph(edges);

    if (isOutputNodedEdges) {
      return OverlayUtil.toLines(graph, isOutputEdges, geomFact);
    }

    labelGraph(graph);
    //for (OverlayEdge e : graph.getEdges()) {  Debug.println(e);  }

    if (isOutputEdges || isOutputResultEdges) {
      return OverlayUtil.toLines(graph, isOutputEdges, geomFact);
    }

    Geometry? result = extractResult(opCode!, graph);

    /**
     * Heuristic check on result area.
     * Catches cases where noding causes vertex to move
     * and make topology graph area "invert".
     */
    if (OverlayUtil.isFloating(pm)) {
      bool isAreaConsistent = OverlayUtil.isResultAreaConsistent(
          inputGeom.getGeometry(0), inputGeom.getGeometry(1), opCode!, result!);
      if (!isAreaConsistent)
        throw TopologyException(
            "Result area inconsistent with overlay operation");
    }
    return result;
  }

  List<EdgeNG> nodeEdges() {
    /**
     * Node the edges, using whatever noder is being used
     */
    EdgeNodingBuilder nodingBuilder = EdgeNodingBuilder(pm!, noder);

    /**
     * Optimize Intersection and Difference by clipping to the
     * result extent, if enabled.
     */
    if (isOptimized) {
      Envelope? clipEnv = OverlayUtil.clippingEnvelope(opCode!, inputGeom, pm!);
      if (clipEnv != null) nodingBuilder.setClipEnvelope(clipEnv);
    }

    List<EdgeNG> mergedEdges =
        nodingBuilder.build(inputGeom.getGeometry(0), inputGeom.getGeometry(1));

    /**
     * Record if an input geometry has collapsed.
     * This is used to avoid trying to locate disconnected edges
     * against a geometry which has collapsed completely.
     */
    inputGeom.setCollapsed(0, !nodingBuilder.hasEdgesFor(0));
    inputGeom.setCollapsed(1, !nodingBuilder.hasEdgesFor(1));

    return mergedEdges;
  }

  OverlayGraph buildGraph(List<EdgeNG> edges) {
    OverlayGraph graph = OverlayGraph();
    for (EdgeNG e in edges) {
      graph.addEdge(e.getCoordinates(), e.createLabel());
    }
    return graph;
  }

  void labelGraph(OverlayGraph graph) {
    OverlayLabeller labeller = OverlayLabeller(graph, inputGeom);
    labeller.computeLabelling();
    labeller.markResultAreaEdges(opCode!);
    labeller.unmarkDuplicateEdgesFromResultArea();
  }

  /**
   * Extracts the result geometry components from the fully labelled topology graph.
   * <p>
   * This method implements the semantic that the result of an
   * intersection operation is homogeneous with highest dimension.
   * In other words,
   * if an intersection has components of a given dimension
   * no lower-dimension components are output.
   * For example, if two polygons intersect in an area,
   * no linestrings or points are included in the result,
   * even if portions of the input do meet in lines or points.
   * This semantic choice makes more sense for typical usage,
   * in which only the highest dimension components are of interest.
   *
   * @param opCode the overlay operation
   * @param graph the topology graph
   * @return the result geometry
   */
  Geometry? extractResult(int opCode, OverlayGraph graph) {
    bool isAllowMixedIntResult = !isStrictMode;

    //--- Build polygons
    List<OverlayEdge> resultAreaEdges = graph.getResultAreaEdges();
    PolygonBuilderNG polyBuilder = PolygonBuilderNG(resultAreaEdges, geomFact);
    List<Polygon> resultPolyList = polyBuilder.getPolygons();
    bool hasResultAreaComponents = resultPolyList.length > 0;

    List<LineString> resultLineList = List.empty(growable: true);
    List<Point> resultPointList = List.empty(growable: true);

    if (!isAreaResultOnly) {
      //--- Build lines
      bool allowResultLines = !hasResultAreaComponents ||
          isAllowMixedIntResult ||
          opCode == SYMDIFFERENCE ||
          opCode == UNION;
      if (allowResultLines) {
        LineBuilderNG lineBuilder = LineBuilderNG(
            inputGeom, graph, hasResultAreaComponents, opCode, geomFact);
        lineBuilder.setStrictMode(isStrictMode);
        resultLineList = lineBuilder.getLines();
      }
      /**
       * Operations with point inputs are handled elsewhere.
       * Only an Intersection op can produce point results
       * from non-point inputs.
       */
      bool hasResultComponents =
          hasResultAreaComponents || resultLineList.length > 0;
      bool allowResultPoints = !hasResultComponents || isAllowMixedIntResult;
      if (opCode == INTERSECTION && allowResultPoints) {
        IntersectionPointBuilder pointBuilder =
            IntersectionPointBuilder(graph, geomFact);
        pointBuilder.setStrictMode(isStrictMode);
        resultPointList = pointBuilder.getPoints();
      }
    }

    if (isEmpty(resultPolyList) &&
        isEmpty(resultLineList) &&
        isEmpty(resultPointList)) return createEmptyResult();

    Geometry resultGeom = OverlayUtil.createResultGeometry(
        resultPolyList, resultLineList, resultPointList, geomFact);
    return resultGeom;
  }

  static bool isEmpty(List? list) {
    return list == null || list.isEmpty;
  }

  Geometry? createEmptyResult() {
    return OverlayUtil.createEmptyResult(
        OverlayUtil.resultDimension(
            opCode!, inputGeom.getDimension(0), inputGeom.getDimension(1)),
        geomFact);
  }
}

/**
 * Computes a robust clipping envelope for a pair of polygonal geometries.
 * The envelope is computed to be large enough to include the full
 * length of all geometry line segments which intersect
 * a given target envelope.
 * This ensures that line segments which might intersect are
 * not perturbed when clipped using {@link RingClipper}.
 *
 * @author Martin Davis
 *
 */
class RobustClipEnvelopeComputer {
  static Envelope getEnvelopeExpanded(
      Geometry a, Geometry b, Envelope targetEnv) {
    RobustClipEnvelopeComputer cec = RobustClipEnvelopeComputer(targetEnv);
    cec.add(a);
    cec.add(b);
    return cec.getEnvelope();
  }

  Envelope targetEnv;
  late Envelope clipEnv;

  RobustClipEnvelopeComputer(this.targetEnv) {
    this.targetEnv = targetEnv;
    clipEnv = targetEnv.copy();
  }

  Envelope getEnvelope() {
    return clipEnv;
  }

  void add(Geometry? g) {
    if (g == null || g.isEmpty()) return;

    if (g is Polygon)
      addPolygon(g);
    else if (g is GeometryCollection) addCollection(g);
  }

  void addCollection(GeometryCollection gc) {
    for (int i = 0; i < gc.getNumGeometries(); i++) {
      Geometry g = gc.getGeometryN(i);
      add(g);
    }
  }

  void addPolygon(Polygon poly) {
    LinearRing shell = poly.getExteriorRing();
    addPolygonRing(shell);

    for (int i = 0; i < poly.getNumInteriorRing(); i++) {
      LinearRing hole = poly.getInteriorRingN(i);
      addPolygonRing(hole);
    }
  }

  /**
   * Adds a polygon ring to the graph. Empty rings are ignored.
   */
  void addPolygonRing(LinearRing ring) {
    // don't add empty lines
    if (ring.isEmpty()) return;

    CoordinateSequence seq = ring.getCoordinateSequence();
    for (int i = 1; i < seq.size(); i++) {
      addSegment(seq.getCoordinate(i - 1), seq.getCoordinate(i));
    }
  }

  void addSegment(Coordinate p1, Coordinate p2) {
    if (intersectsSegment(targetEnv, p1, p2)) {
      clipEnv.expandToIncludeCoordinate(p1);
      clipEnv.expandToIncludeCoordinate(p2);
    }
  }

  static bool intersectsSegment(Envelope env, Coordinate p1, Coordinate p2) {
    /**
     * This is a crude test of whether segment intersects envelope.
     * It could be refined by checking exact intersection.
     * This could be based on the algorithm in the HotPixel.intersectsScaled method.
     */
    return env.intersectsEnvelopeCoordinates(p1, p2);
  }
}

class ElevationCell {
  int numZ = 0;
  double sumZ = 0.0;
  late double avgZ;

  void add(double z) {
    numZ++;
    sumZ += z;
  }

  void compute() {
    avgZ = double.nan;
    if (numZ > 0) avgZ = sumZ / numZ;
  }

  double getZ() {
    return avgZ;
  }
}

/**
 * Implements the logic to compute the full labeling
 * for the edges in an {@link OverlayGraph}.
 *
 * @author mdavis
 *
 */
class OverlayLabeller {
  OverlayGraph _graph;
  InputGeometry _inputGeometry;
  late List<OverlayEdge> _edges;

  OverlayLabeller(this._graph, this._inputGeometry) {
    _edges = _graph.getEdges();
  }

  /**
   * Computes the topological labelling for the edges in the graph.
   *
   */
  void computeLabelling() {
    List<OverlayEdge> nodes = _graph.getNodeEdges();

    _labelAreaNodeEdges(nodes);
    _labelConnectedLinearEdges();

    //TODO: is there a way to avoid scanning all edges in these steps?
    /**
     * At this point collapsed edges labeled with location UNKNOWN
     * must be disconnected from the area edges of the parent.
     * This can occur with a collapsed hole or shell.
     * The edges can be labeled based on their parent ring role (shell or hole).
     */
    _labelCollapsedEdges();
    _labelConnectedLinearEdges();

    _labelDisconnectedEdges();
  }

  /**
   * Labels edges around nodes based on the arrangement
   * of incident area boundary edges.
   * Also propagates the labeling to connected linear edges.
   *
   * @param nodes the nodes to label
   */
  void _labelAreaNodeEdges(List<OverlayEdge> nodes) {
    for (OverlayEdge nodeEdge in nodes) {
      propagateAreaLocations(nodeEdge, 0);
      if (_inputGeometry.hasEdges(1)) {
        propagateAreaLocations(nodeEdge, 1);
      }
    }
  }

  /**
   * Scans around a node CCW, propagating the side labels
   * for a given area geometry to all edges (and their sym)
   * with unknown locations for that geometry.
   * @param e2
   *
   * @param geomIndex the geometry to propagate locations for
   */
  void propagateAreaLocations(OverlayEdge nodeEdge, int geomIndex) {
    /**
     * Only propagate for area geometries
     */
    if (!_inputGeometry.isArea(geomIndex)) return;
    /**
     * No need to propagate if node has only one edge.
     * This handles dangling edges created by overlap limiting.
     */
    if (nodeEdge.degree() == 1) return;

    OverlayEdge? eStart = _findPropagationStartEdge(nodeEdge, geomIndex);
    // no labelled edge found, so nothing to propagate
    if (eStart == null) {
      return;
    }
    // initialize currLoc to location of L side
    int currLoc = eStart.getLocation(geomIndex, Position.LEFT);
    OverlayEdge? e = eStart.oNextOE();

    //Debug.println("\npropagateSideLabels geomIndex = " + geomIndex + " : " + eStart);
    //Debug.print("BEFORE: " + toString(eStart));

    do {
      OverlayLabel label = e!.getLabel();
      if (!label.isBoundary(geomIndex)) {
        /**
         * If this is not a Boundary edge for this input area,
         * its location is now known relative to this input area
         */
        label.setLocationLine(geomIndex, currLoc);
      } else {
        // must be a boundary edge
        Assert.isTrue(label.hasSides(geomIndex));
        /**
         *  This is a boundary edge for the input area geom.
         *  Update the current location from its labels.
         *  Also check for topological consistency.
         */
        int locRight = e.getLocation(geomIndex, Position.RIGHT);
        if (locRight != currLoc) {
          throw TopologyException(
              "side location conflict: arg " + geomIndex.toString());
        }
        int locLeft = e.getLocation(geomIndex, Position.LEFT);
        if (locLeft == Location.NONE) {
          Assert.shouldNeverReachHere(
              "found single null side at " + e.toString());
        }
        currLoc = locLeft;
      }
      e = e.oNextOE();
    } while (e != eStart);
    //Debug.print("AFTER: " + toString(eStart));
  }

  /**
   * Finds a boundary edge for this geom originating at the given
   * node, if one exists.
   * A boundary edge should exist if this is a node on the boundary
   * of the parent area geometry.
   *
   * @param nodeEdge an edge for this node
   * @param geomIndex the parent geometry index
   * @return a boundary edge, or null if no boundary edge exists
   */
  static OverlayEdge? _findPropagationStartEdge(
      OverlayEdge nodeEdge, int geomIndex) {
    OverlayEdge eStart = nodeEdge;
    do {
      OverlayLabel label = eStart.getLabel();
      if (label.isBoundary(geomIndex)) {
        Assert.isTrue(label.hasSides(geomIndex));
        return eStart;
      }
      eStart = eStart.oNext() as OverlayEdge;
    } while (eStart != nodeEdge);
    return null;
  }

  /**
   * At this point collapsed edges with unknown location
   * must be disconnected from the boundary edges of the parent
   * (because otherwise the location would have
   * been propagated from them).
   * They can be now located based on their parent ring role (shell or hole).
   * (This cannot be done earlier, because the location
   * based on the boundary edges must take precedence.
   * There are situations where a collapsed edge has a location
   * which is different to its ring role -
   * e.g. a narrow gore in a polygon, which is in
   * the interior of the reduced polygon, but whose
   * ring role would imply the location EXTERIOR.)
   *
   * Note that collapsed edges can NOT have location determined via a PIP location check,
   * because that is done against the unreduced input geometry,
   * which may give an invalid result due to topology collapse.
   *
   * The labeling is propagated to other connected linear edges,
   * since there may be NOT_PART edges which are connected,
   * and they can be labeled in the same way.
   * (These would get labeled anyway during subsequent disconnected labeling pass,
   * but may be more efficient and accurate to do it here.)
   */
  void _labelCollapsedEdges() {
    for (OverlayEdge edge in _edges) {
      if (edge.getLabel().isLineLocationUnknown(0)) {
        _labelCollapsedEdge(edge, 0);
      }
      if (edge.getLabel().isLineLocationUnknown(1)) {
        _labelCollapsedEdge(edge, 1);
      }
    }
  }

  void _labelCollapsedEdge(OverlayEdge edge, int geomIndex) {
    //Debug.println("\n------  labelCollapsedEdge - geomIndex= " + geomIndex);
    //Debug.print("BEFORE: " + edge.toStringNode());
    OverlayLabel label = edge.getLabel();
    if (!label.isCollapse(geomIndex)) return;
    /**
       * This must be a collapsed edge which is disconnected
       * from any area edges (e.g. a fully collapsed shell or hole).
       * It can be labeled according to its parent source ring role.
       */
    label.setLocationCollapse(geomIndex);
    //Debug.print("AFTER: " + edge.toStringNode());
  }

  /**
   * There can be edges which have unknown location
   * but are connected to a linear edge with known location.
   * In this case linear location is propagated to the connected edges.
   */
  void _labelConnectedLinearEdges() {
    //TODO: can these be merged to avoid two scans?
    _propagateLinearLocations(0);
    if (_inputGeometry.hasEdges(1)) {
      _propagateLinearLocations(1);
    }
  }

  /**
   * Performs a breadth-first graph traversal to find and label
   * connected linear edges.
   *
   * @param geomIndex the index of the input geometry to label
   */
  void _propagateLinearLocations(int geomIndex) {
    // find located linear edges
    List<OverlayEdge> linearEdges =
        _findLinearEdgesWithLocation(_edges, geomIndex);
    if (linearEdges.isEmpty) return;

    Queue<OverlayEdge> edgeStack = Queue();
    for (var e in linearEdges) {
      edgeStack.addLast(e);
    }
    bool isInputLine = _inputGeometry.isLine(geomIndex);
    // traverse connected linear edges, labeling unknown ones
    while (!edgeStack.isEmpty) {
      OverlayEdge lineEdge = edgeStack.removeFirst();
      // assert: lineEdge.getLabel().isLine(geomIndex);

      // for any edges around origin with unknown location for this geomIndex,
      // add those edges to stack to continue traversal
      _propagateLinearLocationAtNode(
          lineEdge, geomIndex, isInputLine, edgeStack);
    }
  }

  static void _propagateLinearLocationAtNode(OverlayEdge eNode, int geomIndex,
      bool isInputLine, Queue<OverlayEdge> edgeStack) {
    int lineLoc = eNode.getLabel().getLineLocation(geomIndex);
    /**
     * If the parent geom is a Line
     * then only propagate EXTERIOR locations.
     */
    if (isInputLine && lineLoc != Location.EXTERIOR) return;

    //Debug.println("propagateLinearLocationAtNode ----- using location for " + geomIndex + " from: " + eNode);
    OverlayEdge? e = eNode.oNextOE();
    do {
      OverlayLabel label = e!.getLabel();
      //Debug.println("check " + geomIndex + ": " + e);
      if (label.isLineLocationUnknown(geomIndex)) {
        /**
         * If edge is not a boundary edge,
         * its location is now known for this area
         */
        label.setLocationLine(geomIndex, lineLoc);
        //Debug.println("*** setting "+ geomIndex + ": " + e);

        /**
         * Add sym edge to stack for graph traversal
         * (Don't add e itself, since e origin node has now been scanned)
         */
        edgeStack.addFirst(e.symOE());
      }
      e = e.oNextOE();
    } while (e != eNode);
  }

  /**
   * Finds all OverlayEdges which are linear
   * (i.e. line or collapsed) and have a known location
   * for the given input geometry.
   *
   * @param geomIndex the index of the input geometry
   * @return list of linear edges with known location
   */
  static List<OverlayEdge> _findLinearEdgesWithLocation(
      List<OverlayEdge> edges, int geomIndex) {
    List<OverlayEdge> linearEdges = List.empty(growable: true);
    for (OverlayEdge edge in edges) {
      OverlayLabel lbl = edge.getLabel();
      // keep if linear with known location
      if (lbl.isLinear(geomIndex) && !lbl.isLineLocationUnknown(geomIndex)) {
        linearEdges.add(edge);
      }
    }
    return linearEdges;
  }

  /**
   * At this point there may still be edges which have unknown location
   * relative to an input geometry.
   * This must be because they are NOT_PART edges for that geometry,
   * and are disconnected from any edges of that geometry.
   * An example of this is rings of one geometry wholly contained
   * in another geometry.
   * The location must be fully determined to compute a
   * correct result for all overlay operations.
   *
   * If the input geometry is an Area the edge location can
   * be determined via a PIP test.
   * If the input is not an Area the location is EXTERIOR.
   */
  void _labelDisconnectedEdges() {
    for (OverlayEdge edge in _edges) {
      //Debug.println("\n------  checking for Disconnected edge " + edge);
      if (edge.getLabel().isLineLocationUnknown(0)) {
        _labelDisconnectedEdge(edge, 0);
      }
      if (edge.getLabel().isLineLocationUnknown(1)) {
        _labelDisconnectedEdge(edge, 1);
      }
    }
  }

  /**
   * Determines the location of an edge relative to a target input geometry.
   * The edge has no location information
   * because it is disconnected from other
   * edges that would provide that information.
   * The location is determined by checking
   * if the edge lies inside the target geometry area (if any).
   *
   * @param edge the edge to label
   * @param geomIndex the input geometry to label against
   */
  void _labelDisconnectedEdge(OverlayEdge edge, int geomIndex) {
    OverlayLabel label = edge.getLabel();
    //Assert.isTrue(label.isNotPart(geomIndex));

    /**
     * if target geom is not an area then
     * edge must be EXTERIOR, since to be
     * INTERIOR it would have been labelled
     * when it was created.
     */
    if (!_inputGeometry.isArea(geomIndex)) {
      label.setLocationAll(geomIndex, Location.EXTERIOR);
      return;
    }
    ;

    //Debug.println("\n------  labelDisconnectedEdge - geomIndex= " + geomIndex);
    //Debug.print("BEFORE: " + edge.toStringNode());
    /**
     * Locate edge in input area using a Point-In-Poly check.
     * This should be safe even with precision reduction,
     * because since the edge has remained disconnected
     * its interior-exterior relationship
     * can be determined relative to the original input geometry.
     */
    //int edgeLoc = locateEdge(geomIndex, edge);
    int edgeLoc = _locateEdgeBothEnds(geomIndex, edge);
    label.setLocationAll(geomIndex, edgeLoc);
    //Debug.print("AFTER: " + edge.toStringNode());
  }

  /**
   * Determines the {@link Location} for an edge within an Area geometry
   * via point-in-polygon location.
   * <p>
   * NOTE this is only safe to use for disconnected edges,
   * since the test is carried out against the original input geometry,
   * and precision reduction may cause incorrect results for edges
   * which are close enough to a boundary to become connected.
   *
   * @param geomIndex the parent geometry index
   * @param edge the edge to locate
   * @return the location of the edge
   */
  int _locateEdge(int geomIndex, OverlayEdge edge) {
    int loc = _inputGeometry.locatePointInArea(geomIndex, edge.orig());
    int edgeLoc =
        loc != Location.EXTERIOR ? Location.INTERIOR : Location.EXTERIOR;
    return edgeLoc;
  }

  /**
   * Determines the {@link Location} for an edge within an Area geometry
   * via point-in-polygon location,
   * by checking that both endpoints are interior to the target geometry.
   * Checking both endpoints ensures correct results in the presence of topology collapse.
   * <p>
   * NOTE this is only safe to use for disconnected edges,
   * since the test is carried out against the original input geometry,
   * and precision reduction may cause incorrect results for edges
   * which are close enough to a boundary to become connected.
   *
   * @param geomIndex the parent geometry index
   * @param edge the edge to locate
   * @return the location of the edge
   */
  int _locateEdgeBothEnds(int geomIndex, OverlayEdge edge) {
    /*
     * To improve the robustness of the point location,
     * check both ends of the edge.
     * EdgeNG is only labelled INTERIOR if both ends are.
     */
    int locOrig = _inputGeometry.locatePointInArea(geomIndex, edge.orig());
    int locDest = _inputGeometry.locatePointInArea(geomIndex, edge.dest());
    bool isInt = locOrig != Location.EXTERIOR && locDest != Location.EXTERIOR;
    int edgeLoc = isInt ? Location.INTERIOR : Location.EXTERIOR;
    return edgeLoc;
  }

  void markResultAreaEdges(int overlayOpCode) {
    for (OverlayEdge edge in _edges) {
      markInResultArea(edge, overlayOpCode);
    }
  }

  /**
   * Marks an edge which forms part of the boundary of the result area.
   * This is determined by the overlay operation being executed,
   * and the location of the edge.
   * The relevant location is either the right side of a boundary edge,
   * or the line location of a non-boundary edge.
   *
   * @param e the edge to mark
   * @param overlayOpCode the overlay operation
   */
  void markInResultArea(OverlayEdge e, int overlayOpCode) {
    OverlayLabel label = e.getLabel();
    if (label.isBoundaryEither() &&
        OverlayNG.isResultOfOp(
            overlayOpCode,
            label.getLocationBoundaryOrLine(0, Position.RIGHT, e.isForward()),
            label.getLocationBoundaryOrLine(
                1, Position.RIGHT, e.isForward()))) {
      e.markInResultArea();
    }
    //Debug.println("markInResultArea: " + e);
  }

  /**
   * Unmarks result area edges where the sym edge
   * is also marked as in the result.
   * This has the effect of merging edge-adjacent result areas,
   * as required by polygon validity rules.
   */
  void unmarkDuplicateEdgesFromResultArea() {
    for (OverlayEdge edge in _edges) {
      if (edge.isInResultAreaBoth()) {
        edge.unmarkFromResultAreaBoth();
      }
    }
  }

  static String toStringEdge(OverlayEdge? nodeEdge) {
    Coordinate orig = nodeEdge!.orig();
    StringBuffer sb = StringBuffer();
    sb.write("Node( " + orig.toString() + " )" + "\n");
    OverlayEdge? e = nodeEdge;
    do {
      sb.write("  -> " + e!.toString());
      if (e.isResultLinked()) {
        sb.write(" Link: ");
        sb.write(e.nextResult());
      }
      sb.write("\n");
      e = e.oNextOE();
    } while (e != nodeEdge);
    return sb.toString();
  }
}

/**
 * Represents the linework for edges in the topology
 * derived from (up to) two parent geometries.
 * An edge may be the result of the merging of
 * two or more edges which have the same linework
 * (although possibly different orientations).
 * In this case the topology information is
 * derived from the merging of the information in the
 * source edges.
 * Merged edges can occur in the following situations
 * <ul>
 * <li>Due to coincident edges of polygonal or linear geometries.
 * <li>Due to topology collapse caused by snapping or rounding
 * of polygonal geometries.
 * </ul>
 * The source edges may have the same parent geometry,
 * or different ones, or a mix of the two.
 *
 * @author mdavis
 *
 */
class EdgeNG {
  /**
   * Tests if the given point sequence
   * is a collapsed line.
   * A collapsed edge has fewer than two distinct points.
   *
   * @param pts the point sequence to check
   * @return true if the points form a collapsed line
   */
  static bool isCollapsed(List<Coordinate> pts) {
    if (pts.length < 2) return true;
    // zero-length line
    if (pts[0].equals2D(pts[1])) return true;
    // TODO: is pts > 2 with equal points ever expected?
    if (pts.length > 2) {
      if (pts[pts.length - 1].equals2D(pts[pts.length - 2])) return true;
    }
    return false;
  }

  List<Coordinate> _pts;

  int _aDim = OverlayLabel.DIM_UNKNOWN;
  int _aDepthDelta = 0;
  bool _aIsHole = false;

  int _bDim = OverlayLabel.DIM_UNKNOWN;
  int _bDepthDelta = 0;
  bool _bIsHole = false;

  EdgeNG(this._pts, EdgeSourceInfo info) {
    _copyInfo(info);
  }

  List<Coordinate> getCoordinates() {
    return _pts;
  }

  Coordinate getCoordinate(int index) {
    return _pts[index];
  }

  int size() {
    return _pts.length;
  }

  bool direction() {
    List<Coordinate> pts = getCoordinates();
    if (pts.length < 2) {
      throw Exception("EdgeNG must have >= 2 points");
    }
    Coordinate p0 = pts[0];
    Coordinate p1 = pts[1];

    Coordinate pn0 = pts[pts.length - 1];
    Coordinate pn1 = pts[pts.length - 2];

    int cmp = 0;
    int cmp0 = p0.compareTo(pn0);
    if (cmp0 != 0) cmp = cmp0;

    if (cmp == 0) {
      int cmp1 = p1.compareTo(pn1);
      if (cmp1 != 0) cmp = cmp1;
    }

    if (cmp == 0) {
      throw Exception(
          "EdgeNG direction cannot be determined because endpoints are equal");
    }

    return cmp == -1;
  }

  /**
   * Compares two coincident edges to determine
   * whether they have the same or opposite direction.
   *
   * @param edge1 an edge
   * @param edge2 an edge
   * @return true if the edges have the same direction, false if not
   */
  bool relativeDirection(EdgeNG edge2) {
    // assert: the edges match (have the same coordinates up to direction)
    if (!getCoordinate(0).equals2D(edge2.getCoordinate(0))) return false;
    if (!getCoordinate(1).equals2D(edge2.getCoordinate(1))) return false;
    return true;
  }

  OverlayLabel createLabel() {
    OverlayLabel lbl = OverlayLabel.empty();
    _initLabel(lbl, 0, _aDim, _aDepthDelta, _aIsHole);
    _initLabel(lbl, 1, _bDim, _bDepthDelta, _bIsHole);
    return lbl;
  }

  /**
   * Populates the label for an edge resulting from an input geometry.
   *
   * <ul>
   * <li>If the edge is not part of the input, the label is left as NOT_PART
   * <li>If input is an Area and the edge is on the boundary
   * (which may include some collapses),
   * edge is marked as an AREA edge and side locations are assigned
   * <li>If input is an Area and the edge is collapsed
   * (depth delta = 0),
   * the label is set to COLLAPSE.
   * The location will be determined later
   * by evaluating the final graph topology.
   * <li>If input is a Line edge is set to a LINE edge.
   * For line edges the line location is not significant
   * (since there is no parent area for which to determine location).
   * </ul>
   *
   * @param lbl
   * @param geomIndex
   * @param dim
   * @param depthDelta
   */
  static void _initLabel(
      OverlayLabel lbl, int geomIndex, int dim, int depthDelta, bool isHole) {
    int dimLabel = _labelDim(dim, depthDelta);

    switch (dimLabel) {
      case OverlayLabel.DIM_NOT_PART:
        lbl.initNotPart(geomIndex);
        break;
      case OverlayLabel.DIM_BOUNDARY:
        lbl.initBoundary(geomIndex, _locationLeft(depthDelta),
            _locationRight(depthDelta), isHole);
        break;
      case OverlayLabel.DIM_COLLAPSE:
        lbl.initCollapse(geomIndex, isHole);
        break;
      case OverlayLabel.DIM_LINE:
        lbl.initLine(geomIndex);
        break;
    }
  }

  static int _labelDim(int dim, int depthDelta) {
    if (dim == Dimension.FALSE) return OverlayLabel.DIM_NOT_PART;

    if (dim == Dimension.L) return OverlayLabel.DIM_LINE;

    // assert: dim is A
    bool isCollapse = depthDelta == 0;
    if (isCollapse) return OverlayLabel.DIM_COLLAPSE;

    return OverlayLabel.DIM_BOUNDARY;
  }

  /**
   * Tests whether the edge is part of a shell in the given geometry.
   * This is only the case if the edge is a boundary.
   *
   * @param geomIndex the index of the geometry
   * @return true if this edge is a boundary and part of a shell
   */
  bool _isShell(int geomIndex) {
    if (geomIndex == 0) {
      return _aDim == OverlayLabel.DIM_BOUNDARY && !_aIsHole;
    }
    return _bDim == OverlayLabel.DIM_BOUNDARY && !_bIsHole;
  }

  static int _locationRight(int depthDelta) {
    int delSignI = _delSign(depthDelta);
    switch (delSignI) {
      case 0:
        return OverlayLabel.LOC_UNKNOWN;
      case 1:
        return Location.INTERIOR;
      case -1:
        return Location.EXTERIOR;
    }
    return OverlayLabel.LOC_UNKNOWN;
  }

  static int _locationLeft(int depthDelta) {
    // TODO: is it always safe to ignore larger depth deltas?
    int delSignI = _delSign(depthDelta);
    switch (delSignI) {
      case 0:
        return OverlayLabel.LOC_UNKNOWN;
      case 1:
        return Location.EXTERIOR;
      case -1:
        return Location.INTERIOR;
    }
    return OverlayLabel.LOC_UNKNOWN;
  }

  static int _delSign(int depthDel) {
    if (depthDel > 0) return 1;
    if (depthDel < 0) return -1;
    return 0;
  }

  void _copyInfo(EdgeSourceInfo info) {
    if (info.getIndex() == 0) {
      _aDim = info.getDimension();
      _aIsHole = info.isHole;
      _aDepthDelta = info.getDepthDelta();
    } else {
      _bDim = info.getDimension();
      _bIsHole = info.isHole;
      _bDepthDelta = info.getDepthDelta();
    }
  }

  /**
   * Merges an edge into this edge,
   * updating the topology info accordingly.
   *
   * @param edge
   */
  void merge(EdgeNG edge) {
    /**
   * Marks this
   * as a shell edge if any contributing edge is a shell.
   * Update hole status first, since it depends on edge dim
   */
    _aIsHole = _isHoleMerged(0, this, edge);
    _bIsHole = _isHoleMerged(1, this, edge);

    if (edge._aDim > _aDim) _aDim = edge._aDim;
    if (edge._bDim > _bDim) _bDim = edge._bDim;

    bool relDir = relativeDirection(edge);
    int flipFactor = relDir ? 1 : -1;
    _aDepthDelta += flipFactor * edge._aDepthDelta;
    _bDepthDelta += flipFactor * edge._bDepthDelta;
    /*
    if (aDepthDelta > 1) {
      Debug.println(this);
    }
    */
  }

  static bool _isHoleMerged(int geomIndex, EdgeNG edge1, EdgeNG edge2) {
    // TOD: this might be clearer with tri-state logic for isHole?
    bool isShell1 = edge1._isShell(geomIndex);
    bool isShell2 = edge2._isShell(geomIndex);
    bool isShellMerged = isShell1 || isShell2;
    // flip since isHole is stored
    return !isShellMerged;
  }

  String toString() {
    String ptsStr = _oStringPts(_pts);

    String aInfo = infoString(0, _aDim, _aIsHole, _aDepthDelta);
    String bInfo = infoString(1, _bDim, _bIsHole, _bDepthDelta);

    return "EdgeNG( " + ptsStr + " ) " + aInfo + "/" + bInfo;
  }

  String toLineString() {
    return WKTWriter.toLineString(_pts);
  }

  static String _oStringPts(List<Coordinate> pts) {
    Coordinate orig = pts[0];
    Coordinate dest = pts[pts.length - 1];
    String dirPtStr = (pts.length > 2) ? ", " + pts[1].toString() : "";
    String ptsStr = orig.toString() + dirPtStr + " .. " + dest.toString();
    return ptsStr;
  }

  static String infoString(int index, int dim, bool isHole, int depthDelta) {
    return (index == 0 ? "A:" : "B:") +
        OverlayLabel.dimensionSymbol(dim) +
        _ringRoleSymbol(dim, isHole) +
        depthDelta.toString(); // force to string
  }

  static String _ringRoleSymbol(int dim, bool isHole) {
    if (_hasAreaParent(dim)) return "" + OverlayLabel.ringRoleSymbol(isHole);
    return "";
  }

  static bool _hasAreaParent(int dim) {
    return dim == OverlayLabel.DIM_BOUNDARY || dim == OverlayLabel.DIM_COLLAPSE;
  }
}

class PolygonBuilderNG {
  GeometryFactory _geometryFactory;
  List<OverlayEdgeRing> _shellList = List.empty(growable: true);
  List<OverlayEdgeRing> _freeHoleList = List.empty(growable: true);
  late bool _isEnforcePolygonal;

  factory PolygonBuilderNG(
      List<OverlayEdge> resultAreaEdges, GeometryFactory gf) {
    return PolygonBuilderNG.enforcedPolygonal(resultAreaEdges, gf, true);
  }

  PolygonBuilderNG.enforcedPolygonal(List<OverlayEdge> resultAreaEdges,
      this._geometryFactory, this._isEnforcePolygonal) {
    _buildRings(resultAreaEdges);
  }

  List<Polygon> getPolygons() {
    return _computePolygons(_shellList);
  }

  List<OverlayEdgeRing> getShellRings() {
    return _shellList;
  }

  List<Polygon> _computePolygons(List<OverlayEdgeRing> shellList) {
    List<Polygon> resultPolyList = List.empty(growable: true);
    // add Polygons for all shells
    for (OverlayEdgeRing er in shellList) {
      Polygon poly = er.toPolygon(_geometryFactory);
      resultPolyList.add(poly);
    }
    return resultPolyList;
  }

  void _buildRings(List<OverlayEdge> resultAreaEdges) {
    _linkResultAreaEdgesMax(resultAreaEdges);
    List<MaximalEdgeRingNG> maxRings = _buildMaximalRings(resultAreaEdges);
    _buildMinimalRings(maxRings);
    _placeFreeHoles(_shellList, _freeHoleList);
    //Assert: every hole on freeHoleList has a shell assigned to it
  }

  void _linkResultAreaEdgesMax(List<OverlayEdge> resultEdges) {
    for (OverlayEdge edge in resultEdges) {
      //Assert.isTrue(edge.isInResult());
      // TODO: find some way to skip nodes which are already linked
      MaximalEdgeRingNG.linkResultAreaMaxRingAtNode(edge);
    }
  }

  /**
   * For all OverlayEdges in result, form them into MaximalEdgeRings
   */
  static List<MaximalEdgeRingNG> _buildMaximalRings(List<OverlayEdge> edges) {
    List<MaximalEdgeRingNG> edgeRings = List.empty(growable: true);
    for (OverlayEdge e in edges) {
      if (e.isInResultArea && e.getLabel().isBoundaryEither()) {
        // if this edge has not yet been processed
        if (e.getEdgeRingMax() == null) {
          MaximalEdgeRingNG er = MaximalEdgeRingNG(e);
          edgeRings.add(er);
        }
      }
    }
    return edgeRings;
  }

  void _buildMinimalRings(List<MaximalEdgeRingNG> maxRings) {
    for (MaximalEdgeRingNG erMax in maxRings) {
      List<OverlayEdgeRing> minRings =
          erMax.buildMinimalRings(_geometryFactory);
      _assignShellsAndHoles(minRings);
    }
  }

  void _assignShellsAndHoles(List<OverlayEdgeRing> minRings) {
    /**
     * Two situations may occur:
     * - the rings are a shell and some holes
     * - rings are a set of holes
     * This code identifies the situation
     * and places the rings appropriately
     */
    OverlayEdgeRing? shell = findSingleShell(minRings);
    if (shell != null) {
      _assignHoles(shell, minRings);
      _shellList.add(shell);
    } else {
      // all rings are holes; their shell will be found later
      _freeHoleList.addAll(minRings);
    }
  }

  /**
   * Finds the single shell, if any, out of
   * a list of minimal rings derived from a maximal ring.
   * The other possibility is that they are a set of (connected) holes,
   * in which case no shell will be found.
   *
   * @return the shell ring, if there is one
   * or null, if all rings are holes
   */
  OverlayEdgeRing? findSingleShell(List<OverlayEdgeRing> edgeRings) {
    int shellCount = 0;
    OverlayEdgeRing? shell = null;
    for (OverlayEdgeRing er in edgeRings) {
      if (!er.isHole()) {
        shell = er;
        shellCount++;
      }
    }
    Assert.isTrue(shellCount <= 1, "found two shells in EdgeRing list");
    return shell;
  }

  /**
   * For the set of minimal rings comprising a maximal ring,
   * assigns the holes to the shell known to contain them.
   * Assigning the holes directly to the shell serves two purposes:
   * <ul>
   * <li>it is faster than using a point-in-polygon check later on.
   * <li>it ensures correctness, since if the PIP test was used the point
   * chosen might lie on the shell, which might return an incorrect result from the
   * PIP test
   * </ul>
   */
  static void _assignHoles(
      OverlayEdgeRing shell, List<OverlayEdgeRing> edgeRings) {
    for (OverlayEdgeRing er in edgeRings) {
      if (er.isHole()) {
        er.setShell(shell);
      }
    }
  }

  /**
   * Place holes have not yet been assigned to a shell.
   * These "free" holes should
   * all be <b>properly</b> contained in their parent shells, so it is safe to use the
   * <code>findEdgeRingContaining</code> method.
   * (This is the case because any holes which are NOT
   * properly contained (i.e. are connected to their
   * parent shell) would have formed part of a MaximalEdgeRingNG
   * and been handled in a previous step).
   *
   * @throws TopologyException if a hole cannot be assigned to a shell
   */
  void _placeFreeHoles(
      List<OverlayEdgeRing> shellList, List<OverlayEdgeRing> freeHoleList) {
    // TODO: use a spatial index to improve performance
    for (OverlayEdgeRing hole in freeHoleList) {
      // only place this hole if it doesn't yet have a shell
      if (hole.getShell() == null) {
        OverlayEdgeRing? shell = hole.findEdgeRingContaining(shellList);
        // only when building a polygon-valid result
        if (_isEnforcePolygonal && shell == null) {
          throw TopologyException("unable to assign free hole to a shell");
        }
        hole.setShell(shell!);
      }
    }
  }
}
