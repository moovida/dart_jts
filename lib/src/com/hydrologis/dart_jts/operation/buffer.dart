part of dart_jts;

/**
 * A value class containing the parameters which
 * specify how a buffer should be constructed.
 * <p>
 * The parameters allow control over:
 * <ul>
 * <li>Quadrant segments (accuracy of approximation for circular arcs)
 * <li>End Cap style
 * <li>Join style
 * <li>Mitre limit
 * <li>whether the buffer is single-sided
 * </ul>
 *
 * @author Martin Davis
 *
 */
class BufferParameters {
  /**
   * Specifies a round line buffer end cap style.
   */
  static const int CAP_ROUND = 1;

  /**
   * Specifies a flat line buffer end cap style.
   */
  static const int CAP_FLAT = 2;

  /**
   * Specifies a square line buffer end cap style.
   */
  static const int CAP_SQUARE = 3;

  /**
   * Specifies a round join style.
   */
  static const int JOIN_ROUND = 1;

  /**
   * Specifies a mitre join style.
   */
  static const int JOIN_MITRE = 2;

  /**
   * Specifies a bevel join style.
   */
  static const int JOIN_BEVEL = 3;

  /**
   * The default number of facets into which to divide a fillet of 90 degrees.
   * A value of 8 gives less than 2% max error in the buffer distance.
   * For a max error of &lt; 1%, use QS = 12.
   * For a max error of &lt; 0.1%, use QS = 18.
   */
  static const int DEFAULT_QUADRANT_SEGMENTS = 8;

  /**
   * The default mitre limit
   * Allows fairly pointy mitres.
   */
  static const double DEFAULT_MITRE_LIMIT = 5.0;

  /**
   * The default simplify factor
   * Provides an accuracy of about 1%, which matches the accuracy of the default Quadrant Segments parameter.
   */
  static const double DEFAULT_SIMPLIFY_FACTOR = 0.01;

  int quadrantSegments = DEFAULT_QUADRANT_SEGMENTS;
  int endCapStyle = CAP_ROUND;
  int joinStyle = JOIN_ROUND;
  double mitreLimit = DEFAULT_MITRE_LIMIT;
  bool _isSingleSided = false;
  double simplifyFactor = DEFAULT_SIMPLIFY_FACTOR;

  /**
   * Creates a default set of parameters
   *
   */
  BufferParameters() {}

  /**
   * Creates a set of parameters with the
   * given quadrantSegments value.
   *
   * @param quadrantSegments the number of quadrant segments to use
   */
  BufferParameters.withQS(int quadrantSegments) {
    setQuadrantSegments(quadrantSegments);
  }

  /**
   * Creates a set of parameters with the
   * given quadrantSegments and endCapStyle values.
   *
   * @param quadrantSegments the number of quadrant segments to use
   * @param endCapStyle the end cap style to use
   */
  BufferParameters.withQS_ECS(int quadrantSegments, int endCapStyle) {
    setQuadrantSegments(quadrantSegments);
    setEndCapStyle(endCapStyle);
  }

  /**
   * Creates a set of parameters with the
   * given parameter values.
   *
   * @param quadrantSegments the number of quadrant segments to use
   * @param endCapStyle the end cap style to use
   * @param joinStyle the join style to use
   * @param mitreLimit the mitre limit to use
   */
  BufferParameters.withQS_ECS_JS_ML(
      int quadrantSegments, int endCapStyle, int joinStyle, double mitreLimit) {
    setQuadrantSegments(quadrantSegments);
    setEndCapStyle(endCapStyle);
    setJoinStyle(joinStyle);
    setMitreLimit(mitreLimit);
  }

  /**
   * Gets the number of quadrant segments which will be used
   *
   * @return the number of quadrant segments
   */
  int getQuadrantSegments() {
    return quadrantSegments;
  }

  /**
   * Sets the number of line segments used to approximate an angle fillet.
   * <ul>
   * <li>If <tt>quadSegs</tt> &gt;= 1, joins are round, and <tt>quadSegs</tt> indicates the number of
   * segments to use to approximate a quarter-circle.
   * <li>If <tt>quadSegs</tt> = 0, joins are bevelled (flat)
   * <li>If <tt>quadSegs</tt> &lt; 0, joins are mitred, and the value of qs
   * indicates the mitre ration limit as
   * <pre>
   * mitreLimit = |<tt>quadSegs</tt>|
   * </pre>
   * </ul>
   * For round joins, <tt>quadSegs</tt> determines the maximum
   * error in the approximation to the true buffer curve.
   * The default value of 8 gives less than 2% max error in the buffer distance.
   * For a max error of &lt; 1%, use QS = 12.
   * For a max error of &lt; 0.1%, use QS = 18.
   * The error is always less than the buffer distance
   * (in other words, the computed buffer curve is always inside the true
   * curve).
   *
   * @param quadSegs the number of segments in a fillet for a quadrant
   */
  void setQuadrantSegments(int quadSegs) {
    quadrantSegments = quadSegs;

    /**
     * Indicates how to construct fillets.
     * If qs >= 1, fillet is round, and qs indicates number of
     * segments to use to approximate a quarter-circle.
     * If qs = 0, fillet is bevelled flat (i.e. no filleting is performed)
     * If qs < 0, fillet is mitred, and absolute value of qs
     * indicates maximum length of mitre according to
     *
     * mitreLimit = |qs|
     */
    if (quadrantSegments == 0) joinStyle = JOIN_BEVEL;
    if (quadrantSegments < 0) {
      joinStyle = JOIN_MITRE;
      mitreLimit = quadrantSegments.abs().toDouble();
    }

    if (quadSegs <= 0) {
      quadrantSegments = 1;
    }

    /**
     * If join style was set by the quadSegs value,
     * use the default for the actual quadrantSegments value.
     */
    if (joinStyle != JOIN_ROUND) {
      quadrantSegments = DEFAULT_QUADRANT_SEGMENTS;
    }
  }

  /**
   * Computes the maximum distance error due to a given level
   * of approximation to a true arc.
   *
   * @param quadSegs the number of segments used to approximate a quarter-circle
   * @return the error of approximation
   */
  static double bufferDistanceError(int quadSegs) {
    double alpha = math.pi / 2.0 / quadSegs;
    return 1 - math.cos(alpha / 2.0);
  }

  /**
   * Gets the end cap style.
   *
   * @return the end cap style
   */
  int getEndCapStyle() {
    return endCapStyle;
  }

  /**
   * Specifies the end cap style of the generated buffer.
   * The styles supported are {@link #CAP_ROUND}, {@link #CAP_FLAT}, and {@link #CAP_SQUARE}.
   * The default is CAP_ROUND.
   *
   * @param endCapStyle the end cap style to specify
   */
  void setEndCapStyle(int endCapStyle) {
    this.endCapStyle = endCapStyle;
  }

  /**
   * Gets the join style
   *
   * @return the join style code
   */
  int getJoinStyle() {
    return joinStyle;
  }

  /**
   * Sets the join style for outside (reflex) corners between line segments.
   * Allowable values are {@link #JOIN_ROUND} (which is the default),
   * {@link #JOIN_MITRE} and {link JOIN_BEVEL}.
   *
   * @param joinStyle the code for the join style
   */
  void setJoinStyle(int joinStyle) {
    this.joinStyle = joinStyle;
  }

  /**
   * Gets the mitre ratio limit.
   *
   * @return the limit value
   */
  double getMitreLimit() {
    return mitreLimit;
  }

  /**
   * Sets the limit on the mitre ratio used for very sharp corners.
   * The mitre ratio is the ratio of the distance from the corner
   * to the end of the mitred offset corner.
   * When two line segments meet at a sharp angle,
   * a miter join will extend far beyond the original geometry.
   * (and in the extreme case will be infinitely far.)
   * To prevent unreasonable geometry, the mitre limit
   * allows controlling the maximum length of the join corner.
   * Corners with a ratio which exceed the limit will be beveled.
   *
   * @param mitreLimit the mitre ratio limit
   */
  void setMitreLimit(double mitreLimit) {
    this.mitreLimit = mitreLimit;
  }

  /**
   * Sets whether the computed buffer should be single-sided.
   * A single-sided buffer is constructed on only one side of each input line.
   * <p>
   * The side used is determined by the sign of the buffer distance:
   * <ul>
   * <li>a positive distance indicates the left-hand side
   * <li>a negative distance indicates the right-hand side
   * </ul>
   * The single-sided buffer of point geometries is
   * the same as the regular buffer.
   * <p>
   * The End Cap Style for single-sided buffers is
   * always ignored,
   * and forced to the equivalent of <tt>CAP_FLAT</tt>.
   *
   * @param isSingleSided true if a single-sided buffer should be constructed
   */
  void setSingleSided(bool isSingleSided) {
    this._isSingleSided = isSingleSided;
  }

  /**
   * Tests whether the buffer is to be generated on a single side only.
   *
   * @return true if the generated buffer is to be single-sided
   */
  bool get isSingleSided => _isSingleSided;

  /**
   * Gets the simplify factor.
   *
   * @return the simplify factor
   */
  double getSimplifyFactor() {
    return simplifyFactor;
  }

  /**
   * Sets the factor used to determine the simplify distance tolerance
   * for input simplification.
   * Simplifying can increase the performance of computing buffers.
   * Generally the simplify factor should be greater than 0.
   * Values between 0.01 and .1 produce relatively good accuracy for the generate buffer.
   * Larger values sacrifice accuracy in return for performance.
   *
   * @param simplifyFactor a value greater than or equal to zero.
   */
  void setSimplifyFactor(double simplifyFactor) {
    this.simplifyFactor = simplifyFactor < 0 ? 0 : simplifyFactor;
  }
}

/**
 * Computes the buffer of a geometry, for both positive and negative buffer distances.
 * <p>
 * In GIS, the positive (or negative) buffer of a geometry is defined as
 * the Minkowski sum (or difference) of the geometry
 * with a circle of radius equal to the absolute value of the buffer distance.
 * In the CAD/CAM world buffers are known as <i>offset curves</i>.
 * In morphological analysis the
 * operation of positive and negative buffering
 * is referred to as <i>erosion</i> and <i>dilation</i>
 * <p>
 * The buffer operation always returns a polygonal result.
 * The negative or zero-distance buffer of lines and points is always an empty {@link Polygon}.
 * <p>
 * Since true buffer curves may contain circular arcs,
 * computed buffer polygons are only approximations to the true geometry.
 * The user can control the accuracy of the approximation by specifying
 * the number of linear segments used to approximate arcs.
 * This is specified via {@link BufferParameters#setQuadrantSegments(int)} or {@link #setQuadrantSegments(int)}.
 * <p>
 * The <b>end cap style</b> of a linear buffer may be {@link BufferParameters#setEndCapStyle(int) specified}. The
 * following end cap styles are supported:
 * <ul>
 * <li>{@link BufferParameters#CAP_ROUND} - the usual round end caps
 * <li>{@link BufferParameters#CAP_FLAT} - end caps are truncated flat at the line ends
 * <li>{@link BufferParameters#CAP_SQUARE} - end caps are squared off at the buffer distance beyond the line ends
 * </ul>
 * <p>
 * The <b>join style</b> of the corners in a buffer may be {@link BufferParameters#setJoinStyle(int) specified}. The
 * following join styles are supported:
 * <ul>
 * <li>{@link BufferParameters#JOIN_ROUND} - the usual round join
 * <li>{@link BufferParameters#JOIN_MITRE} - corners are "sharp" (up to a {@link BufferParameters#getMitreLimit() distance limit})
 * <li>{@link BufferParameters#JOIN_BEVEL} - corners are beveled (clipped off).
 * </ul>
 * <p>
 * The buffer algorithm can perform simplification on the input to increase performance.
 * The simplification is performed a way that always increases the buffer area
 * (so that the simplified input covers the original input).
 * The degree of simplification can be {@link BufferParameters#setSimplifyFactor(double) specified},
 * with a {@link BufferParameters#DEFAULT_SIMPLIFY_FACTOR default} used otherwise.
 * Note that if the buffer distance is zero then so is the computed simplify tolerance,
 * no matter what the simplify factor.
 *
 * @version 1.7
 */
class BufferOp {
  /**
   * Specifies a round line buffer end cap style.
   * @deprecated use BufferParameters
   */
  static final int CAP_ROUND = BufferParameters.CAP_ROUND;

  /**
   * Specifies a butt (or flat) line buffer end cap style.
   * @deprecated use BufferParameters
   */
  static final int CAP_BUTT = BufferParameters.CAP_FLAT;

  /**
   * Specifies a butt (or flat) line buffer end cap style.
   * @deprecated use BufferParameters
   */
  static final int CAP_FLAT = BufferParameters.CAP_FLAT;

  /**
   * Specifies a square line buffer end cap style.
   * @deprecated use BufferParameters
   */
  static final int CAP_SQUARE = BufferParameters.CAP_SQUARE;

  /**
   * A number of digits of precision which leaves some computational "headroom"
   * for floating point operations.
   *
   * This value should be less than the decimal precision of double-precision values (16).
   */
  static int MAX_PRECISION_DIGITS = 12;

  /**
   * Compute a scale factor to limit the precision of
   * a given combination of Geometry and buffer distance.
   * The scale factor is determined by
   * the number of digits of precision in the (geometry + buffer distance),
   * limited by the supplied <code>maxPrecisionDigits</code> value.
   * <p>
   * The scale factor is based on the absolute magnitude of the (geometry + buffer distance).
   * since this determines the number of digits of precision which must be handled.
   *
   * @param g the Geometry being buffered
   * @param distance the buffer distance
   * @param maxPrecisionDigits the max # of digits that should be allowed by
   *          the precision determined by the computed scale factor
   *
   * @return a scale factor for the buffer computation
   */
  static double precisionScaleFactor(
      Geometry g, double distance, int maxPrecisionDigits) {
    Envelope env = g.getEnvelopeInternal();
    double envMax = MathUtils.max4(env.getMaxX().abs(), env.getMaxY().abs(),
        env.getMinX().abs(), env.getMinY().abs());

    double expandByDistance = distance > 0.0 ? distance : 0.0;
    double bufEnvMax = envMax + 2 * expandByDistance;

    // the smallest power of 10 greater than the buffer envelope
    int bufEnvPrecisionDigits =
        (math.log(bufEnvMax) / math.log(10) + 1.0).toInt();
    int minUnitLog10 = maxPrecisionDigits - bufEnvPrecisionDigits;

    double scaleFactor = math.pow(10.0, minUnitLog10).toDouble();
    return scaleFactor;
  }

  /*
   static double OLDprecisionScaleFactor(Geometry g,
      double distance,
    int maxPrecisionDigits)
  {
    Envelope env = g.getEnvelopeInternal();
    double envSize = Math.max(env.getHeight(), env.getWidth());
    double expandByDistance = distance > 0.0 ? distance : 0.0;
    double bufEnvSize = envSize + 2 * expandByDistance;

    // the smallest power of 10 greater than the buffer envelope
    int bufEnvLog10 = (int) (Math.log(bufEnvSize) / Math.log(10) + 1.0);
    int minUnitLog10 = bufEnvLog10 - maxPrecisionDigits;
    // scale factor is inverse of min Unit size, so flip sign of exponent
    double scaleFactor = Math.pow(10.0, -minUnitLog10);
    return scaleFactor;
  }
  */

  /**
   * Computes the buffer of a geometry for a given buffer distance.
   *
   * @param g the geometry to buffer
   * @param distance the buffer distance
   * @return the buffer of the input geometry
   */
  static Geometry bufferOp(Geometry g, double distance) {
    BufferOp gBuf = new BufferOp(g);
    Geometry geomBuf = gBuf.getResultGeometry(distance);
//BufferDebug.saveBuffer(geomBuf);
    //BufferDebug.runCount++;
    return geomBuf;
  }

  /**
   * Computes the buffer for a geometry for a given buffer distance
   * and accuracy of approximation.
   *
   * @param g the geometry to buffer
   * @param distance the buffer distance
   * @param params the buffer parameters to use
   * @return the buffer of the input geometry
   *
   */
  static Geometry bufferOpWithParams(
      Geometry g, double distance, BufferParameters params) {
    BufferOp bufOp = new BufferOp.withParams(g, params);
    Geometry geomBuf = bufOp.getResultGeometry(distance);
    return geomBuf;
  }

  /**
   * Computes the buffer for a geometry for a given buffer distance
   * and accuracy of approximation.
   *
   * @param g the geometry to buffer
   * @param distance the buffer distance
   * @param quadrantSegments the number of segments used to approximate a quarter circle
   * @return the buffer of the input geometry
   *
   */
  static Geometry bufferOp3(Geometry g, double distance, int quadrantSegments) {
    BufferOp bufOp = new BufferOp(g);
    bufOp.setQuadrantSegments(quadrantSegments);
    Geometry geomBuf = bufOp.getResultGeometry(distance);
    return geomBuf;
  }

  /**
   * Computes the buffer for a geometry for a given buffer distance
   * and accuracy of approximation.
   *
   * @param g the geometry to buffer
   * @param distance the buffer distance
   * @param quadrantSegments the number of segments used to approximate a quarter circle
   * @param endCapStyle the end cap style to use
   * @return the buffer of the input geometry
   *
   */
  static Geometry bufferOp4(
      Geometry g, double distance, int quadrantSegments, int endCapStyle) {
    BufferOp bufOp = new BufferOp(g);
    bufOp.setQuadrantSegments(quadrantSegments);
    bufOp.setEndCapStyle(endCapStyle);
    Geometry geomBuf = bufOp.getResultGeometry(distance);
    return geomBuf;
  }

  Geometry argGeom;
  double distance = 0.0;

  BufferParameters bufParams = new BufferParameters();

  Geometry? resultGeometry = null;
  // Exception? saveException; // debugging only

  /**
   * Initializes a buffer computation for the given geometry
   *
   * @param g the geometry to buffer
   */
  BufferOp(this.argGeom);

  /**
   * Initializes a buffer computation for the given geometry
   * with the given set of parameters
   *
   * @param g the geometry to buffer
   * @param bufParams the buffer parameters to use
   */
  BufferOp.withParams(this.argGeom, this.bufParams);

  /**
   * Specifies the end cap style of the generated buffer.
   * The styles supported are {@link BufferParameters#CAP_ROUND}, {@link BufferParameters#CAP_FLAT}, and {@link BufferParameters#CAP_SQUARE}.
   * The default is CAP_ROUND.
   *
   * @param endCapStyle the end cap style to specify
   */
  void setEndCapStyle(int endCapStyle) {
    bufParams.setEndCapStyle(endCapStyle);
  }

  /**
   * Sets the number of segments used to approximate a angle fillet
   *
   * @param quadrantSegments the number of segments in a fillet for a quadrant
   */
  void setQuadrantSegments(int quadrantSegments) {
    bufParams.setQuadrantSegments(quadrantSegments);
  }

  /**
   * Returns the buffer computed for a geometry for a given buffer distance.
   *
   * @param distance the buffer distance
   * @return the buffer of the input geometry
   */
  Geometry getResultGeometry(double distance) {
    this.distance = distance;
    computeGeometry();
    return resultGeometry!;
  }

  void computeGeometry() {
    bufferOriginalPrecision();
    if (resultGeometry != null) return;

    PrecisionModel argPM = argGeom.getFactory().getPrecisionModel();
    if (argPM.getType() == PrecisionModel.FIXED)
      bufferFixedPrecision(argPM);
    else
      bufferReducedPrecision();
  }

  void bufferReducedPrecision() {
    // try and compute with decreasing precision
    for (int precDigits = MAX_PRECISION_DIGITS; precDigits >= 0; precDigits--) {
      try {
        bufferReducedPrecisionWithDigits(precDigits);
      } catch (ex) {
        // update the saved exception to reflect the new input geometry
        // saveException = ex as Exception;
        print(ex);
        // don't propagate the exception - it will be detected by fact that resultGeometry is null
      }
      if (resultGeometry != null) return;
    }

    // tried everything - have to bail
    throw StateError("bufferReducedPrecision error");
  }

  void bufferOriginalPrecision() {
    try {
      // use fast noding by default
      BufferBuilder bufBuilder = new BufferBuilder(bufParams);
      resultGeometry = bufBuilder.buffer(argGeom, distance);
    } catch (ex) {
      // saveException = ex as Exception;
      // don't propagate the exception - it will be detected by fact that resultGeometry is null

      // testing ONLY - propagate exception
      // throw ex;
    }
  }

  void bufferReducedPrecisionWithDigits(int precisionDigits) {
    double sizeBasedScaleFactor =
        precisionScaleFactor(argGeom, distance, precisionDigits);
//    System.out.println("recomputing with precision scale factor = " + sizeBasedScaleFactor);

    PrecisionModel fixedPM =
        PrecisionModel.fixedPrecision(sizeBasedScaleFactor);
    bufferFixedPrecision(fixedPM);
  }

  void bufferFixedPrecision(PrecisionModel fixedPM) {
    Noder noder = new ScaledNoder(
        new MCIndexSnapRounder(PrecisionModel.fixedPrecision(1.0)),
        fixedPM.getScale());

    BufferBuilder bufBuilder = new BufferBuilder(bufParams);
    bufBuilder.setWorkingPrecisionModel(fixedPM);
    bufBuilder.setNoder(noder);
    // this may throw an exception, if robustness errors are encountered
    resultGeometry = bufBuilder.buffer(argGeom, distance);
  }
}

/**
 * Builds the buffer geometry for a given input geometry and precision model.
 * Allows setting the level of approximation for circular arcs,
 * and the precision model in which to carry out the computation.
 * <p>
 * When computing buffers in floating point double-precision
 * it can happen that the process of iterated noding can fail to converge (terminate).
 * In this case a {@link TopologyException} will be thrown.
 * Retrying the computation in a fixed precision
 * can produce more robust results.
 *
 * @version 1.7
 */
class BufferBuilder {
  /**
   * Compute the change in depth as an edge is crossed from R to L
   */
  static int depthDelta(Label label) {
    int lLoc = label.getLocationWithPosIndex(0, Position.LEFT);
    int rLoc = label.getLocationWithPosIndex(0, Position.RIGHT);
    if (lLoc == Location.INTERIOR && rLoc == Location.EXTERIOR)
      return 1;
    else if (lLoc == Location.EXTERIOR && rLoc == Location.INTERIOR) return -1;
    return 0;
  }

  BufferParameters bufParams;

  PrecisionModel? workingPrecisionModel;
  Noder? workingNoder;
  GeometryFactory? geomFact;
  PlanarGraph? graph;
  EdgeList edgeList = new EdgeList();

  /**
   * Creates a new BufferBuilder,
   * using the given parameters.
   *
   * @param bufParams the buffer parameters to use
   */
  BufferBuilder(this.bufParams);

  /**
   * Sets the precision model to use during the curve computation and noding,
   * if it is different to the precision model of the Geometry.
   * If the precision model is less than the precision of the Geometry precision model,
   * the Geometry must have previously been rounded to that precision.
   *
   * @param pm the precision model to use
   */
  void setWorkingPrecisionModel(PrecisionModel pm) {
    workingPrecisionModel = pm;
  }

  /**
   * Sets the {@link Noder} to use during noding.
   * This allows choosing fast but non-robust noding, or slower
   * but robust noding.
   *
   * @param noder the noder to use
   */
  void setNoder(Noder noder) {
    workingNoder = noder;
  }

  Geometry buffer(Geometry g, double distance) {
    PrecisionModel? precisionModel = workingPrecisionModel;
    if (precisionModel == null) precisionModel = g.getPrecisionModel();

    // factory must be the same as the one used by the input
    geomFact = g.getFactory();

    OffsetCurveBuilder curveBuilder =
        new OffsetCurveBuilder(precisionModel, bufParams);

    OffsetCurveSetBuilder curveSetBuilder =
        new OffsetCurveSetBuilder(g, distance, curveBuilder);

    List bufferSegStrList = curveSetBuilder.getCurves();

    // short-circuit test
    if (bufferSegStrList.isEmpty) {
      return createEmptyResultGeometry();
    }

//BufferDebug.runCount++;
//String filename = "run" + BufferDebug.runCount + "_curves";
//System.out.println("saving " + filename);
//BufferDebug.saveEdges(bufferEdgeList, filename);
// DEBUGGING ONLY
//WKTWriter wktWriter = new WKTWriter();
//Debug.println("Rings: " + wktWriter.write(convertSegStrings(bufferSegStrList.iterator())));
//wktWriter.setMaxCoordinatesPerLine(10);
//System.out.println(wktWriter.writeFormatted(convertSegStrings(bufferSegStrList.iterator())));

    computeNodedEdges(bufferSegStrList, precisionModel);
    graph = new PlanarGraph.withFactory(OverlayNodeFactory());
    graph!.addEdges(edgeList.getEdges());

    List subgraphList = createSubgraphs(graph!);
    PolygonBuilder polyBuilder = new PolygonBuilder(geomFact!);
    buildSubgraphs(subgraphList, polyBuilder);
    List<Polygon> resultPolyList = polyBuilder.getPolygons();

    // just in case...
    if (resultPolyList.isEmpty) {
      return createEmptyResultGeometry();
    }

    Geometry resultGeom = geomFact!.buildGeometry(resultPolyList);
    return resultGeom;
  }

  Noder getNoder(PrecisionModel precisionModel) {
    if (workingNoder != null) return workingNoder!;

    // otherwise use a fast (but non-robust) noder
    MCIndexNoder noder = MCIndexNoder.empty();
    LineIntersector li = RobustLineIntersector();
    li.setPrecisionModel(precisionModel);
    noder.setSegmentIntersector(IntersectionAdder(li));
//    Noder noder = new IteratedNoder(precisionModel);
    return noder;
//    Noder noder = new SimpleSnapRounder(precisionModel);
//    Noder noder = new MCIndexSnapRounder(precisionModel);
//    Noder noder = new ScaledNoder(new MCIndexSnapRounder(new PrecisionModel(1.0)),
//                                  precisionModel.getScale());
  }

  void computeNodedEdges(List bufferSegStrList, PrecisionModel precisionModel) {
    Noder noder = getNoder(precisionModel);
    noder.computeNodes(bufferSegStrList);
    List nodedSegStrings = noder.getNodedSubstrings();
// DEBUGGING ONLY
//BufferDebug.saveEdges(nodedEdges, "run" + BufferDebug.runCount + "_nodedEdges");

    for (SegmentString segStr in nodedSegStrings) {
      /**
       * Discard edges which have zero length,
       * since they carry no information and cause problems with topology building
       */
      List<Coordinate> pts = segStr.getCoordinates();
      if (pts.length == 2 && pts[0].equals2D(pts[1])) continue;

      Label oldLabel = segStr.getData() as Label;
      Edge edge =
          new Edge(segStr.getCoordinates(), new Label.fromLabel(oldLabel));
      insertUniqueEdge(edge);
    }
    //saveEdges(edgeList.getEdges(), "run" + runCount + "_collapsedEdges");
  }

  /**
   * Inserted edges are checked to see if an identical edge already exists.
   * If so, the edge is not inserted, but its label is merged
   * with the existing edge.
   */
  void insertUniqueEdge(Edge e) {
//<FIX> MD 8 Oct 03  speed up identical edge lookup
    // fast lookup
    Edge? existingEdge = edgeList.findEqualEdge(e);

    // If an identical edge already exists, simply update its label
    if (existingEdge != null) {
      Label existingLabel = existingEdge.getLabel()!;

      Label labelToMerge = e.getLabel()!;
      // check if new edge is in reverse direction to existing edge
      // if so, must flip the label before merging it
      if (!existingEdge.isPointwiseEqual(e)) {
        labelToMerge = new Label.fromLabel(e.getLabel()!);
        labelToMerge.flip();
      }
      existingLabel.merge(labelToMerge);

      // compute new depth delta of sum of edges
      int mergeDelta = depthDelta(labelToMerge);
      int existingDelta = existingEdge.getDepthDelta();
      int newDelta = existingDelta + mergeDelta;
      existingEdge.setDepthDelta(newDelta);
    } else {
      // no matching existing edge was found
      // add this new edge to the list of edges in this graph
      //e.setName(name + edges.size());
      edgeList.add(e);
      e.setDepthDelta(depthDelta(e.getLabel()!));
    }
  }

  List createSubgraphs(PlanarGraph graph) {
    List subgraphList = [];
    for (Node node in graph.getNodes()) {
      if (!node.isVisited()) {
        BufferSubgraph subgraph = new BufferSubgraph();
        subgraph.create(node);
        subgraphList.add(subgraph);
      }
    }

    /**
     * Sort the subgraphs in descending order of their rightmost coordinate.
     * This ensures that when the Polygons for the subgraphs are built,
     * subgraphs for shells will have been built before the subgraphs for
     * any holes they contain.
     */
    subgraphList.sort(
        (o1, o2) => (o1 as BufferSubgraph).compareTo((o2 as BufferSubgraph)));
    var list = subgraphList.reversed.toList();
    return list;
  }

  /**
   * Completes the building of the input subgraphs by depth-labelling them,
   * and adds them to the PolygonBuilder.
   * The subgraph list must be sorted in rightmost-coordinate order.
   *
   * @param subgraphList the subgraphs to build
   * @param polyBuilder the PolygonBuilder which will build the final polygons
   */
  void buildSubgraphs(List subgraphList, PolygonBuilder polyBuilder) {
    List processedGraphs = [];
    for (BufferSubgraph subgraph in subgraphList) {
      Coordinate p = subgraph.getRightmostCoordinate()!;
//      int outsideDepth = 0;
//      if (polyBuilder.containsPoint(p))
//        outsideDepth = 1;
      SubgraphDepthLocater locater = new SubgraphDepthLocater(processedGraphs);
      int outsideDepth = locater.getDepth(p);
//      try {
      subgraph.computeDepth(outsideDepth);
//      }
//      catch (RuntimeException ex) {
//        // debugging only
//        //subgraph.saveDirEdges();
//        throw ex;
//      }
      subgraph.findResultEdges();
      processedGraphs.add(subgraph);
      polyBuilder.add(subgraph.getDirectedEdges(), subgraph.getNodes());
    }
  }

  static Geometry convertSegStrings(Iterator it) {
    GeometryFactory fact = new GeometryFactory.defaultPrecision();
    List<Geometry> lines = [];
    while (it.moveNext()) {
      SegmentString ss = it.current;
      LineString line = fact.createLineString(ss.getCoordinates());
      lines.add(line);
    }
    return fact.buildGeometry(lines);
  }

  /**
   * Gets the standard result for an empty buffer.
   * Since buffer always returns a polygonal result,
   * this is chosen to be an empty polygon.
   *
   * @return the empty result geometry
   */
  Geometry createEmptyResultGeometry() {
    Geometry emptyGeom = geomFact!.createPolygonEmpty();
    return emptyGeom;
  }
}

/**
 * A connected subset of the graph of
 * {@link DirectedEdge}s and {@link Node}s.
 * Its edges will generate either
 * <ul>
 * <li> a single polygon in the complete buffer, with zero or more holes, or
 * <li> one or more connected holes
 * </ul>
 *
 *
 * @version 1.7
 */
class BufferSubgraph implements Comparable {
  late RightmostEdgeFinder finder;
  List dirEdgeList = [];
  List nodes = [];
  Coordinate? rightMostCoord = null;
  Envelope? env = null;

  BufferSubgraph() {
    finder = new RightmostEdgeFinder();
  }

  List getDirectedEdges() {
    return dirEdgeList;
  }

  List getNodes() {
    return nodes;
  }

  /**
   * Computes the envelope of the edges in the subgraph.
   * The envelope is cached after being computed.
   *
   * @return the envelope of the graph.
   */
  Envelope getEnvelope() {
    if (env == null) {
      Envelope edgeEnv = new Envelope.empty();
      for (DirectedEdge dirEdge in dirEdgeList) {
        List<Coordinate> pts = dirEdge.getEdge().getCoordinates();
        for (int i = 0; i < pts.length - 1; i++) {
          edgeEnv.expandToIncludeCoordinate(pts[i]);
        }
      }
      env = edgeEnv;
    }
    return env!;
  }

  /**
   * Gets the rightmost coordinate in the edges of the subgraph
   */
  Coordinate? getRightmostCoordinate() {
    return rightMostCoord;
  }

  /**
   * Creates the subgraph consisting of all edges reachable from this node.
   * Finds the edges in the graph and the rightmost coordinate.
   *
   * @param node a node to start the graph traversal from
   */
  void create(Node node) {
    addReachable(node);
    finder.findEdge(dirEdgeList);
    rightMostCoord = finder.getCoordinate();
  }

  /**
   * Adds all nodes and edges reachable from this node to the subgraph.
   * Uses an explicit stack to avoid a large depth of recursion.
   *
   * @param node a node known to be in the subgraph
   */
  void addReachable(Node startNode) {
    Queue nodeStack = new Queue();
    nodeStack.add(startNode);
    while (!nodeStack.isEmpty) {
      Node node = nodeStack.removeLast();
      add(node, nodeStack);
    }
  }

  /**
   * Adds the argument node and all its out edges to the subgraph
   * @param node the node to add
   * @param nodeStack the current set of nodes being traversed
   */
  void add(Node node, Queue nodeStack) {
    node.setVisited(true);
    nodes.add(node);
    for (Iterator i = (node.getEdges() as DirectedEdgeStar).iterator();
        i.moveNext();) {
      DirectedEdge de = i.current;
      dirEdgeList.add(de);
      DirectedEdge sym = de.getSym();
      Node symNode = sym.getNode()!;
      /**
       * NOTE: this is a depth-first traversal of the graph.
       * This will cause a large depth of recursion.
       * It might be better to do a breadth-first traversal.
       */
      if (!symNode.isVisited()) nodeStack.addLast(symNode);
    }
  }

  void clearVisitedEdges() {
    for (DirectedEdge de in dirEdgeList) {
      de.setVisited(false);
    }
  }

  void computeDepth(int outsideDepth) {
    clearVisitedEdges();
    // find an outside edge to assign depth to
    DirectedEdge de = finder.getEdge()!;
    Node n = de.getNode()!;
    Label label = de.getLabel()!;
    // right side of line returned by finder is on the outside
    de.setEdgeDepths(Position.RIGHT, outsideDepth);
    copySymDepths(de);

    //computeNodeDepth(n, de);
    computeDepths(de);
  }

  /**
   * Compute depths for all dirEdges via breadth-first traversal of nodes in graph
   * @param startEdge edge to start processing with
   */
  // <FIX> MD - use iteration & queue rather than recursion, for speed and robustness
  void computeDepths(DirectedEdge startEdge) {
    Set nodesVisited = new HashSet();
    Queue nodeQueue = new Queue();

    Node startNode = startEdge.getNode()!;
    nodeQueue.addLast(startNode);
    nodesVisited.add(startNode);
    startEdge.setVisited(true);

    var tmp = 0;
    while (!nodeQueue.isEmpty) {
      // print("count =  ${tmp++}");
//System.out.println(nodes.size() + " queue: " + nodeQueue.size());
      Node n = nodeQueue.removeFirst();
      nodesVisited.add(n);
      // compute depths around node, starting at this edge since it has depths assigned
      computeNodeDepth(n);

      // add all adjacent nodes to process queue,
      // unless the node has been visited already
      for (Iterator i = (n.getEdges() as DirectedEdgeStar).iterator();
          i.moveNext();) {
        DirectedEdge de = i.current;
        DirectedEdge sym = de.getSym();
        if (sym.isVisited()) continue;
        Node adjNode = sym.getNode()!;
        if (!(nodesVisited.contains(adjNode))) {
          nodeQueue.addLast(adjNode);
          nodesVisited.add(adjNode);
        }
      }
    }
  }

  void computeNodeDepth(Node n) {
    // find a visited dirEdge to start at
    DirectedEdge? startEdge = null;
    for (Iterator i = (n.getEdges() as DirectedEdgeStar).iterator();
        i.moveNext();) {
      DirectedEdge de = i.current;
      if (de.isVisited() || de.getSym().isVisited()) {
        startEdge = de;
        break;
      }
    }
    // MD - testing  Result: breaks algorithm
    //if (startEdge == null) return;

    // only compute string append if assertion would fail
    if (startEdge == null)
      throw new TopologyException("unable to find edge to compute depths at " +
          n.getCoordinate().toString());

    (n.getEdges() as DirectedEdgeStar).computeDepths(startEdge);

    // copy depths to sym edges
    for (Iterator i = (n.getEdges() as DirectedEdgeStar).iterator();
        i.moveNext();) {
      DirectedEdge de = i.current;
      de.setVisited(true);
      copySymDepths(de);
    }
  }

  void copySymDepths(DirectedEdge de) {
    DirectedEdge sym = de.getSym();
    sym.setDepth(Position.LEFT, de.getDepth(Position.RIGHT));
    sym.setDepth(Position.RIGHT, de.getDepth(Position.LEFT));
  }

  /**
   * Find all edges whose depths indicates that they are in the result area(s).
   * Since we want polygon shells to be
   * oriented CW, choose dirEdges with the interior of the result on the RHS.
   * Mark them as being in the result.
   * Interior Area edges are the result of dimensional collapses.
   * They do not form part of the result area boundary.
   */
  void findResultEdges() {
    for (DirectedEdge de in dirEdgeList) {
      /**
       * Select edges which have an interior depth on the RHS
       * and an exterior depth on the LHS.
       * Note that because of weird rounding effects there may be
       * edges which have negative depths!  Negative depths
       * count as "outside".
       */
      // <FIX> - handle negative depths
      if (de.getDepth(Position.RIGHT) >= 1 &&
          de.getDepth(Position.LEFT) <= 0 &&
          !de.isInteriorAreaEdge()) {
        de.setInResult(true);
//Debug.print("in result "); Debug.println(de);
      }
    }
  }

  /**
   * BufferSubgraphs are compared on the x-value of their rightmost Coordinate.
   * This defines a partial ordering on the graphs such that:
   * <p>
   * g1 >= g2 <==> Ring(g2) does not contain Ring(g1)
   * <p>
   * where Polygon(g) is the buffer polygon that is built from g.
   * <p>
   * This relationship is used to sort the BufferSubgraphs so that shells are guaranteed to
   * be built before holes.
   */
  int compareTo(dynamic o) {
    BufferSubgraph graph = o as BufferSubgraph;
    if (this.rightMostCoord!.x < graph.rightMostCoord!.x) {
      return -1;
    }
    if (this.rightMostCoord!.x > graph.rightMostCoord!.x) {
      return 1;
    }
    return 0;
  }

/*
// DEBUGGING only - comment out
   static final String SAVE_DIREDGES = "saveDirEdges";
   static int saveCount = 0;
   void saveDirEdges()
  {
    GeometryFactory fact = new GeometryFactory();
    for (Iterator it = dirEdgeList.iterator(); it.hasNext(); ) {
      DirectedEdge de = (DirectedEdge) it.next();
      double dx = de.getDx();
      double dy = de.getDy();
      Coordinate p0 = de.getCoordinate();
      double ang = Math.atan2(dy, dx);
      Coordinate p1 = new Coordinate(
          p0.x + .4 * Math.cos(ang),
          p0.y + .4 * Math.sin(ang));
//      DebugFeature.add(SAVE_DIREDGES,
//                       fact.createLineString(new List<Coordinate> { p0, p1 } ),
//                       de.getDepth(Position.LEFT) + "/" + de.getDepth(Position.RIGHT)
//                       );
    }
  String filepath = "x:\\jts\\testBuffer\\dirEdges" + saveCount++ + ".jml";
    DebugFeature.saveFeatures(SAVE_DIREDGES, filepath);
  }
  */
}

/**
 * A RightmostEdgeFinder find the DirectedEdge in a list which has the highest coordinate,
 * and which is oriented L to R at that point. (I.e. the right side is on the RHS of the edge.)
 *
 * @version 1.7
 */
class RightmostEdgeFinder {
  // Coordinate extremeCoord;
  int minIndex = -1;
  Coordinate? minCoord = null;
  DirectedEdge? minDe = null;
  DirectedEdge? orientedDe = null;

  /**
   * A RightmostEdgeFinder finds the DirectedEdge with the rightmost coordinate.
   * The DirectedEdge returned is guaranteed to have the R of the world on its RHS.
   */
  RightmostEdgeFinder() {}

  DirectedEdge? getEdge() {
    return orientedDe;
  }

  Coordinate? getCoordinate() {
    return minCoord;
  }

  void findEdge(List dirEdgeList) {
    /**
     * Check all forward DirectedEdges only.  This is still general,
     * because each edge has a forward DirectedEdge.
     */
    for (DirectedEdge de in dirEdgeList) {
      if (!de.isForward()) continue;
      checkForRightmostCoordinate(de);
    }

    /**
     * If the rightmost point is a node, we need to identify which of
     * the incident edges is rightmost.
     */
    Assert.isTrue(minIndex != 0 || minCoord!.equals(minDe!.getCoordinate()!),
        "inconsistency in rightmost processing");
    if (minIndex == 0) {
      findRightmostEdgeAtNode();
    } else {
      findRightmostEdgeAtVertex();
    }
    /**
     * now check that the extreme side is the R side.
     * If not, use the sym instead.
     */
    orientedDe = minDe;
    int rightmostSide = getRightmostSide(minDe!, minIndex);
    if (rightmostSide == Position.LEFT) {
      orientedDe = minDe!.getSym();
    }
  }

  void findRightmostEdgeAtNode() {
    Node node = minDe!.getNode()!;
    DirectedEdgeStar star = node.getEdges() as DirectedEdgeStar;
    minDe = star.getRightmostEdge();
    // the DirectedEdge returned by the previous call is not
    // necessarily in the forward direction. Use the sym edge if it isn't.
    if (!minDe!.isForward()) {
      minDe = minDe!.getSym();
      minIndex = minDe!.getEdge().getCoordinates().length - 1;
    }
  }

  void findRightmostEdgeAtVertex() {
    /**
     * The rightmost point is an interior vertex, so it has a segment on either side of it.
     * If these segments are both above or below the rightmost point, we need to
     * determine their relative orientation to decide which is rightmost.
     */
    List<Coordinate> pts = minDe!.getEdge().getCoordinates();
    Assert.isTrue(minIndex > 0 && minIndex < pts.length,
        "rightmost point expected to be interior vertex of edge");
    Coordinate pPrev = pts[minIndex - 1];
    Coordinate pNext = pts[minIndex + 1];
    int orientation = Orientation.index(minCoord!, pNext, pPrev);
    bool usePrev = false;
    // both segments are below min point
    if (pPrev.y < minCoord!.y &&
        pNext.y < minCoord!.y &&
        orientation == Orientation.COUNTERCLOCKWISE) {
      usePrev = true;
    } else if (pPrev.y > minCoord!.y &&
        pNext.y > minCoord!.y &&
        orientation == Orientation.CLOCKWISE) {
      usePrev = true;
    }
    // if both segments are on the same side, do nothing - either is safe
    // to select as a rightmost segment
    if (usePrev) {
      minIndex = minIndex - 1;
    }
  }

  void checkForRightmostCoordinate(DirectedEdge de) {
    List<Coordinate> coord = de.getEdge().getCoordinates();
    for (int i = 0; i < coord.length - 1; i++) {
      // only check vertices which are the start or end point of a non-horizontal segment
      // <FIX> MD 19 Sep 03 - NO!  we can test all vertices, since the rightmost must have a non-horiz segment adjacent to it
      if (minCoord == null || coord[i].x > minCoord!.x) {
        minDe = de;
        minIndex = i;
        minCoord = coord[i];
      }
      //}
    }
  }

  int getRightmostSide(DirectedEdge de, int index) {
    int side = getRightmostSideOfSegment(de, index);
    if (side < 0) side = getRightmostSideOfSegment(de, index - 1);
    if (side < 0) {
      // reaching here can indicate that segment is horizontal
      //Assert.shouldNeverReachHere("problem with finding rightmost side of segment at " + de.getCoordinate());
      // testing only
      minCoord = null;
      checkForRightmostCoordinate(de);
    }
    return side;
  }

  int getRightmostSideOfSegment(DirectedEdge de, int i) {
    Edge e = de.getEdge();
    var coord = e.getCoordinates();

    if (i < 0 || i + 1 >= coord.length) return -1;
    if (coord[i].y == coord[i + 1].y)
      return -1; // indicates edge is parallel to x-axis

    int pos = Position.LEFT;
    if (coord[i].y < coord[i + 1].y) pos = Position.RIGHT;
    return pos;
  }
}

/**
 * Locates a subgraph inside a set of subgraphs,
 * in order to determine the outside depth of the subgraph.
 * The input subgraphs are assumed to have had depths
 * already calculated for their edges.
 *
 * @version 1.7
 */
class SubgraphDepthLocater {
  List subgraphs;
  LineSegment seg = LineSegment.empty();

  SubgraphDepthLocater(this.subgraphs);

  int getDepth(Coordinate p) {
    List stabbedSegments = findStabbedSegments(p);
    // if no segments on stabbing line subgraph must be outside all others.
    if (stabbedSegments.isEmpty) return 0;

    var sorted = List.from(stabbedSegments);
    sorted.sort((o1, o2) => o1.compareTo(o2));

    DepthSegment ds = sorted.first; // min is requested
    return ds.leftDepth;
  }

  /**
   * Finds all non-horizontal segments intersecting the stabbing line.
   * The stabbing line is the ray to the right of stabbingRayLeftPt.
   *
   * @param stabbingRayLeftPt the left-hand origin of the stabbing line
   * @return a List of {@link DepthSegments} intersecting the stabbing line
   */
  List findStabbedSegments(Coordinate stabbingRayLeftPt) {
    List stabbedSegments = [];
    for (BufferSubgraph bsg in subgraphs) {
      // optimization - don't bother checking subgraphs which the ray does not intersect
      Envelope env = bsg.getEnvelope();
      if (stabbingRayLeftPt.y < env.getMinY() ||
          stabbingRayLeftPt.y > env.getMaxY()) continue;

      findStabbedSegments3(
          stabbingRayLeftPt, bsg.getDirectedEdges(), stabbedSegments);
    }
    return stabbedSegments;
  }

  /**
   * Finds all non-horizontal segments intersecting the stabbing line
   * in the list of dirEdges.
   * The stabbing line is the ray to the right of stabbingRayLeftPt.
   *
   * @param stabbingRayLeftPt the left-hand origin of the stabbing line
   * @param stabbedSegments the current list of {@link DepthSegments} intersecting the stabbing line
   */
  void findStabbedSegments3(
      Coordinate stabbingRayLeftPt, List dirEdges, List stabbedSegments) {
    /**
     * Check all forward DirectedEdges only.  This is still general,
     * because each Edge has a forward DirectedEdge.
     */
    for (DirectedEdge de in dirEdges) {
      if (!de.isForward()) continue;
      findStabbedSegmentsDE(stabbingRayLeftPt, de, stabbedSegments);
    }
  }

  /**
   * Finds all non-horizontal segments intersecting the stabbing line
   * in the input dirEdge.
   * The stabbing line is the ray to the right of stabbingRayLeftPt.
   *
   * @param stabbingRayLeftPt the left-hand origin of the stabbing line
   * @param stabbedSegments the current list of {@link DepthSegments} intersecting the stabbing line
   */
  void findStabbedSegmentsDE(Coordinate stabbingRayLeftPt, DirectedEdge dirEdge,
      List stabbedSegments) {
    List<Coordinate> pts = dirEdge.getEdge().getCoordinates();
    for (int i = 0; i < pts.length - 1; i++) {
      seg.p0 = pts[i];
      seg.p1 = pts[i + 1];
      // ensure segment always points upwards
      if (seg.p0.y > seg.p1.y) seg.reverse();

      // skip segment if it is left of the stabbing line
      double maxx = math.max(seg.p0.x, seg.p1.x);
      if (maxx < stabbingRayLeftPt.x) continue;

      // skip horizontal segments (there will be a non-horizontal one carrying the same depth info
      if (seg.isHorizontal()) continue;

      // skip if segment is above or below stabbing line
      if (stabbingRayLeftPt.y < seg.p0.y || stabbingRayLeftPt.y > seg.p1.y)
        continue;

      // skip if stabbing ray is right of the segment
      if (Orientation.index(seg.p0, seg.p1, stabbingRayLeftPt) ==
          Orientation.RIGHT) continue;

      // stabbing line cuts this segment, so record it
      int depth = dirEdge.getDepth(Position.LEFT);
      // if segment direction was flipped, use RHS depth instead
      if (!seg.p0.equals(pts[i])) depth = dirEdge.getDepth(Position.RIGHT);
      DepthSegment ds = new DepthSegment(seg, depth);
      stabbedSegments.add(ds);
    }
  }
}

/**
 * A segment from a directed edge which has been assigned a depth value
 * for its sides.
 */
class DepthSegment implements Comparable {
  late LineSegment upwardSeg;
  int leftDepth = 0;

  DepthSegment(LineSegment seg, int depth) {
// input seg is assumed to be normalized
    upwardSeg = new LineSegment.fromSegment(seg);
//upwardSeg.normalize();
    this.leftDepth = depth;
  }

  /**
   * Defines a comparison operation on DepthSegments
   * which orders them left to right.
   * Assumes the segments are normalized.
   * <p>
   * The definition of the ordering is:
   * <ul>
   * <li>-1 : if DS1.seg is left of or below DS2.seg (DS1 < DS2)
   * <li>1 : if   DS1.seg is right of or above DS2.seg (DS1 > DS2)
   * <li>0 : if the segments are identical
   * </ul>
   *
   * KNOWN BUGS:
   * <ul>
   * <li>The logic does not obey the {@link Comparator.compareTo} contract.
   * This is acceptable for the intended usage, but may cause problems if used with some
   * utilities in the Java standard library (e.g. {@link Lists.sort()}.
   * </ul>
   *
   * @param obj a DepthSegment
   * @return the comparison value
   */
  int compareTo(dynamic obj) {
    DepthSegment other = obj as DepthSegment;

// fast check if segments are trivially ordered along X
    if (upwardSeg.minX() >= other.upwardSeg.maxX()) return 1;
    if (upwardSeg.maxX() <= other.upwardSeg.minX()) return -1;

    /**
     * try and compute a determinate orientation for the segments.
     * Test returns 1 if other is left of this (i.e. this > other)
     */
    int orientIndex = upwardSeg.orientationIndex(other.upwardSeg);
    if (orientIndex != 0) return orientIndex;

    /**
     * If comparison between this and other is indeterminate,
     * try the opposite call order.
     * The sign of the result needs to be flipped.
     */
    orientIndex = -1 * other.upwardSeg.orientationIndex(upwardSeg);
    if (orientIndex != 0) return orientIndex;

// otherwise, use standard lexocographic segment ordering
    return upwardSeg.compareTo(other.upwardSeg);
  }

  /**
   * Compare two collinear segments for left-most ordering.
   * If segs are vertical, use vertical ordering for comparison.
   * If segs are equal, return 0.
   * Segments are assumed to be directed so that the second coordinate is >= to the first
   * (e.g. up and to the right).
   *
   * @param seg0 a segment to compare
   * @param seg1 a segment to compare
   * @return
   */
  int compareX(LineSegment seg0, LineSegment seg1) {
    int compare0 = seg0.p0.compareTo(seg1.p0);
    if (compare0 != 0) return compare0;
    return seg0.p1.compareTo(seg1.p1);
  }

  String toString() {
    return upwardSeg.toString();
  }
}

/**
 * Computes the raw offset curve for a
 * single {@link Geometry} component (ring, line or point).
 * A raw offset curve line is not noded -
 * it may contain self-intersections (and usually will).
 * The final buffer polygon is computed by forming a topological graph
 * of all the noded raw curves and tracing outside contours.
 * The points in the raw curve are rounded
 * to a given {@link PrecisionModel}.
 *
 * @version 1.7
 */
class OffsetCurveBuilder {
  double distance = 0.0;
  PrecisionModel precisionModel;
  BufferParameters bufParams;

  OffsetCurveBuilder(this.precisionModel, this.bufParams);

  /**
   * Gets the buffer parameters being used to generate the curve.
   *
   * @return the buffer parameters being used
   */
  BufferParameters getBufferParameters() {
    return bufParams;
  }

  /**
   * This method handles single points as well as LineStrings.
   * LineStrings are assumed <b>not</b> to be closed (the function will not
   * fail for closed lines, but will generate superfluous line caps).
   *
   * @param inputPts the vertices of the line to offset
   * @param distance the offset distance
   *
   * @return a Coordinate array representing the curve
   * or null if the curve is empty
   */
  List<Coordinate>? getLineCurve(List<Coordinate> inputPts, double distance) {
    this.distance = distance;

    // a zero or negative width buffer of a line/point is empty
    if (distance < 0.0 && !bufParams.isSingleSided) return null;
    if (distance == 0.0) return null;

    double posDistance = distance.abs();
    OffsetSegmentGenerator segGen = getSegGen(posDistance);
    if (inputPts.length <= 1) {
      computePointCurve(inputPts[0], segGen);
    } else {
      if (bufParams.isSingleSided) {
        bool isRightSide = distance < 0.0;
        computeSingleSidedBufferCurve(inputPts, isRightSide, segGen);
      } else
        computeLineBufferCurve(inputPts, segGen);
    }

    List<Coordinate> lineCoord = segGen.getCoordinates();
    return lineCoord;
  }

  /**
   * This method handles the degenerate cases of single points and lines,
   * as well as rings.
   *
   * @return a Coordinate array representing the curve
   * or null if the curve is empty
   */
  List<Coordinate>? getRingCurve(
      List<Coordinate> inputPts, int side, double distance) {
    this.distance = distance;
    if (inputPts.length <= 2) return getLineCurve(inputPts, distance);

    // optimize creating ring for for zero distance
    if (distance == 0.0) {
      return copyCoordinates(inputPts);
    }
    OffsetSegmentGenerator segGen = getSegGen(distance);
    computeRingBufferCurve(inputPts, side, segGen);
    return segGen.getCoordinates();
  }

  List<Coordinate>? getOffsetCurve(List<Coordinate> inputPts, double distance) {
    this.distance = distance;

    // a zero width offset curve is empty
    if (distance == 0.0) return null;

    bool isRightSide = distance < 0.0;
    double posDistance = distance.abs();
    OffsetSegmentGenerator segGen = getSegGen(posDistance);
    if (inputPts.length <= 1) {
      computePointCurve(inputPts[0], segGen);
    } else {
      computeOffsetCurve(inputPts, isRightSide, segGen);
    }
    List<Coordinate> curvePts = segGen.getCoordinates();
    // for right side line is traversed in reverse direction, so have to reverse generated line
    if (isRightSide) CoordinateArrays.reverse(curvePts);
    return curvePts;
  }

  static List<Coordinate> copyCoordinates(List<Coordinate> pts) {
    List<Coordinate> copy = []; //..length = (pts.length);
    for (int i = 0; i < pts.length; i++) {
      copy.add(Coordinate.fromCoordinate(pts[i]));
      // copy[i] = Coordinate.fromCoordinate(pts[i]);
    }
    return copy;
  }

  OffsetSegmentGenerator getSegGen(double distance) {
    return new OffsetSegmentGenerator(precisionModel, bufParams, distance);
  }

  /**
   * Computes the distance tolerance to use during input
   * line simplification.
   *
   * @param distance the buffer distance
   * @return the simplification tolerance
   */
  double simplifyTolerance(double bufDistance) {
    return bufDistance * bufParams.getSimplifyFactor();
  }

  void computePointCurve(Coordinate pt, OffsetSegmentGenerator segGen) {
    switch (bufParams.getEndCapStyle()) {
      case BufferParameters.CAP_ROUND:
        segGen.createCircle(pt);
        break;
      case BufferParameters.CAP_SQUARE:
        segGen.createSquare(pt);
        break;
      // otherwise curve is empty (e.g. for a butt cap);
    }
  }

  void computeLineBufferCurve(
      List<Coordinate> inputPts, OffsetSegmentGenerator segGen) {
    double distTol = simplifyTolerance(distance);

    //--------- compute points for left side of line
    // Simplify the appropriate side of the line before generating
    List<Coordinate> simp1 =
        BufferInputLineSimplifier.simplify(inputPts, distTol);
    // MD - used for testing only (to eliminate simplification)
//    List<Coordinate> simp1 = inputPts;

    int n1 = simp1.length - 1;
    segGen.initSideSegments(simp1[0], simp1[1], Position.LEFT);
    for (int i = 2; i <= n1; i++) {
      segGen.addNextSegment(simp1[i], true);
    }
    segGen.addLastSegment();
    // add line cap for end of line
    segGen.addLineEndCap(simp1[n1 - 1], simp1[n1]);

    //---------- compute points for right side of line
    // Simplify the appropriate side of the line before generating
    List<Coordinate> simp2 =
        BufferInputLineSimplifier.simplify(inputPts, -distTol);
    // MD - used for testing only (to eliminate simplification)
//    List<Coordinate> simp2 = inputPts;
    int n2 = simp2.length - 1;

    // since we are traversing line in opposite order, offset position is still LEFT
    segGen.initSideSegments(simp2[n2], simp2[n2 - 1], Position.LEFT);
    for (int i = n2 - 2; i >= 0; i--) {
      segGen.addNextSegment(simp2[i], true);
    }
    segGen.addLastSegment();
    // add line cap for start of line
    segGen.addLineEndCap(simp2[1], simp2[0]);

    segGen.closeRing();
  }

  /*
   void OLDcomputeLineBufferCurve(List<Coordinate> inputPts)
  {
    int n = inputPts.length - 1;
    
    // compute points for left side of line
    initSideSegments(inputPts[0], inputPts[1], Position.LEFT);
    for (int i = 2; i <= n; i++) {
      addNextSegment(inputPts[i], true);
    }
    addLastSegment();
    // add line cap for end of line
    addLineEndCap(inputPts[n - 1], inputPts[n]);

    // compute points for right side of line
    initSideSegments(inputPts[n], inputPts[n - 1], Position.LEFT);
    for (int i = n - 2; i >= 0; i--) {
      addNextSegment(inputPts[i], true);
    }
    addLastSegment();
    // add line cap for start of line
    addLineEndCap(inputPts[1], inputPts[0]);

    vertexList.closeRing();
  }
  */

  void computeSingleSidedBufferCurve(List<Coordinate> inputPts,
      bool isRightSide, OffsetSegmentGenerator segGen) {
    double distTol = simplifyTolerance(distance);

    if (isRightSide) {
      // add original line
      segGen.addSegments(inputPts, true);

      //---------- compute points for right side of line
      // Simplify the appropriate side of the line before generating
      List<Coordinate> simp2 =
          BufferInputLineSimplifier.simplify(inputPts, -distTol);
      // MD - used for testing only (to eliminate simplification)
      //    List<Coordinate> simp2 = inputPts;
      int n2 = simp2.length - 1;

      // since we are traversing line in opposite order, offset position is still LEFT
      segGen.initSideSegments(simp2[n2], simp2[n2 - 1], Position.LEFT);
      segGen.addFirstSegment();
      for (int i = n2 - 2; i >= 0; i--) {
        segGen.addNextSegment(simp2[i], true);
      }
    } else {
      // add original line
      segGen.addSegments(inputPts, false);

      //--------- compute points for left side of line
      // Simplify the appropriate side of the line before generating
      List<Coordinate> simp1 =
          BufferInputLineSimplifier.simplify(inputPts, distTol);
      // MD - used for testing only (to eliminate simplification)
//      List<Coordinate> simp1 = inputPts;

      int n1 = simp1.length - 1;
      segGen.initSideSegments(simp1[0], simp1[1], Position.LEFT);
      segGen.addFirstSegment();
      for (int i = 2; i <= n1; i++) {
        segGen.addNextSegment(simp1[i], true);
      }
    }
    segGen.addLastSegment();
    segGen.closeRing();
  }

  void computeOffsetCurve(List<Coordinate> inputPts, bool isRightSide,
      OffsetSegmentGenerator segGen) {
    double distTol = simplifyTolerance(distance);

    if (isRightSide) {
      //---------- compute points for right side of line
      // Simplify the appropriate side of the line before generating
      List<Coordinate> simp2 =
          BufferInputLineSimplifier.simplify(inputPts, -distTol);
      // MD - used for testing only (to eliminate simplification)
      //    List<Coordinate> simp2 = inputPts;
      int n2 = simp2.length - 1;

      // since we are traversing line in opposite order, offset position is still LEFT
      segGen.initSideSegments(simp2[n2], simp2[n2 - 1], Position.LEFT);
      segGen.addFirstSegment();
      for (int i = n2 - 2; i >= 0; i--) {
        segGen.addNextSegment(simp2[i], true);
      }
    } else {
      //--------- compute points for left side of line
      // Simplify the appropriate side of the line before generating
      List<Coordinate> simp1 =
          BufferInputLineSimplifier.simplify(inputPts, distTol);
      // MD - used for testing only (to eliminate simplification)
//      List<Coordinate> simp1 = inputPts;

      int n1 = simp1.length - 1;
      segGen.initSideSegments(simp1[0], simp1[1], Position.LEFT);
      segGen.addFirstSegment();
      for (int i = 2; i <= n1; i++) {
        segGen.addNextSegment(simp1[i], true);
      }
    }
    segGen.addLastSegment();
  }

  void computeRingBufferCurve(
      List<Coordinate> inputPts, int side, OffsetSegmentGenerator segGen) {
    // simplify input line to improve performance
    double distTol = simplifyTolerance(distance);
    // ensure that correct side is simplified
    if (side == Position.RIGHT) distTol = -distTol;
    List<Coordinate> simp =
        BufferInputLineSimplifier.simplify(inputPts, distTol);
//    List<Coordinate> simp = inputPts;

    int n = simp.length - 1;
    segGen.initSideSegments(simp[n - 1], simp[0], side);
    for (int i = 1; i <= n; i++) {
      bool addStartPoint = i != 1;
      segGen.addNextSegment(simp[i], addStartPoint);
    }
    segGen.closeRing();
  }
}

/**
 * Generates segments which form an offset curve.
 * Supports all end cap and join options
 * provided for buffering.
 * This algorithm implements various heuristics to
 * produce smoother, simpler curves which are
 * still within a reasonable tolerance of the
 * true curve.
 *
 * @author Martin Davis
 *
 */
class OffsetSegmentGenerator {
  /**
   * Factor which controls how close offset segments can be to
   * skip adding a filler or mitre.
   */
  static final double OFFSET_SEGMENT_SEPARATION_FACTOR = 1.0E-3;

  /**
   * Factor which controls how close curve vertices on inside turns can be to be snapped
   */
  static final double INSIDE_TURN_VERTEX_SNAP_DISTANCE_FACTOR = 1.0E-3;

  /**
   * Factor which controls how close curve vertices can be to be snapped
   */
  static final double CURVE_VERTEX_SNAP_DISTANCE_FACTOR = 1.0E-6;

  /**
   * Factor which determines how short closing segs can be for round buffers
   */
  static final int MAX_CLOSING_SEG_LEN_FACTOR = 80;

  /**
   * the max error of approximation (distance) between a quad segment and the true fillet curve
   */
  double maxCurveSegmentError = 0.0;

  /**
   * The angle quantum with which to approximate a fillet curve
   * (based on the input # of quadrant segments)
   */
  double filletAngleQuantum = 0.0;

  /**
   * The Closing Segment Length Factor controls how long
   * "closing segments" are.  Closing segments are added
   * at the middle of inside corners to ensure a smoother
   * boundary for the buffer offset curve.
   * In some cases (particularly for round joins with default-or-better
   * quantization) the closing segments can be made quite short.
   * This substantially improves performance (due to fewer intersections being created).
   *
   * A closingSegFactor of 0 results in lines to the corner vertex
   * A closingSegFactor of 1 results in lines halfway to the corner vertex
   * A closingSegFactor of 80 results in lines 1/81 of the way to the corner vertex
   * (this option is reasonable for the very common default situation of round joins
   * and quadrantSegs >= 8)
   */
  int closingSegLengthFactor = 1;

  late OffsetSegmentString segList;
  double distance = 0.0;
  PrecisionModel precisionModel;
  BufferParameters bufParams;
  LineIntersector? li;

  late Coordinate s0, s1, s2;
  LineSegment seg0 = new LineSegment.empty();
  LineSegment seg1 = new LineSegment.empty();
  LineSegment offset0 = new LineSegment.empty();
  LineSegment offset1 = new LineSegment.empty();
  int side = 0;
  bool _hasNarrowConcaveAngle = false;

  OffsetSegmentGenerator(this.precisionModel, this.bufParams, double distance) {
    // compute intersections in full precision, to provide accuracy
    // the points are rounded as they are inserted into the curve line
    li = new RobustLineIntersector();
    filletAngleQuantum = math.pi / 2.0 / bufParams.getQuadrantSegments();

    /**
     * Non-round joins cause issues with short closing segments, so don't use
     * them. In any case, non-round joins only really make sense for relatively
     * small buffer distances.
     */
    if (bufParams.getQuadrantSegments() >= 8 &&
        bufParams.getJoinStyle() == BufferParameters.JOIN_ROUND)
      closingSegLengthFactor = MAX_CLOSING_SEG_LEN_FACTOR;
    init(distance);
  }

  /**
   * Tests whether the input has a narrow concave angle
   * (relative to the offset distance).
   * In this case the generated offset curve will contain self-intersections
   * and heuristic closing segments.
   * This is expected behaviour in the case of Buffer curves.
   * For pure Offset Curves,
   * the output needs to be further treated
   * before it can be used.
   *
   * @return true if the input has a narrow concave angle
   */
  bool hasNarrowConcaveAngleWW() {
    return _hasNarrowConcaveAngle;
  }

  void init(double distance) {
    this.distance = distance;
    maxCurveSegmentError = distance * (1 - math.cos(filletAngleQuantum / 2.0));
    segList = new OffsetSegmentString();
    segList.setPrecisionModel(precisionModel);
    /**
     * Choose the min vertex separation as a small fraction of the offset distance.
     */
    segList
        .setMinimumVertexDistance(distance * CURVE_VERTEX_SNAP_DISTANCE_FACTOR);
  }

  void initSideSegments(Coordinate s1, Coordinate s2, int side) {
    this.s1 = s1;
    this.s2 = s2;
    this.side = side;
    seg1.setCoordinates(s1, s2);
    computeOffsetSegment(seg1, side, distance, offset1);
  }

  List<Coordinate> getCoordinates() {
    List<Coordinate> pts = segList.getCoordinates();
    return pts;
  }

  void closeRing() {
    segList.closeRing();
  }

  void addSegments(List<Coordinate> pt, bool isForward) {
    segList.addPts(pt, isForward);
  }

  void addFirstSegment() {
    segList.addPt(offset1.p0);
  }

  /**
   * Add last offset point
   */
  void addLastSegment() {
    segList.addPt(offset1.p1);
  }

  // static double MAX_CLOSING_SEG_LEN = 3.0;

  void addNextSegment(Coordinate p, bool addStartPoint) {
    // s0-s1-s2 are the coordinates of the previous segment and the current one
    s0 = s1;
    s1 = s2;
    s2 = p;
    seg0.setCoordinates(s0, s1);
    computeOffsetSegment(seg0, side, distance, offset0);
    seg1.setCoordinates(s1, s2);
    computeOffsetSegment(seg1, side, distance, offset1);

    // do nothing if points are equal
    if (s1.equals(s2)) return;

    int orientation = Orientation.index(s0, s1, s2);
    bool outsideTurn = (orientation == Orientation.CLOCKWISE &&
            side == Position.LEFT) ||
        (orientation == Orientation.COUNTERCLOCKWISE && side == Position.RIGHT);

    if (orientation == 0) {
      // lines are collinear
      addCollinear(addStartPoint);
    } else if (outsideTurn) {
      addOutsideTurn(orientation, addStartPoint);
    } else {
      // inside turn
      addInsideTurn(orientation, addStartPoint);
    }
  }

  void addCollinear(bool addStartPoint) {
    /**
     * This test could probably be done more efficiently,
     * but the situation of exact collinearity should be fairly rare.
     */
    li!.computeIntersection(s0, s1, s1, s2);
    int numInt = li!.getIntersectionNum();
    /**
     * if numInt is < 2, the lines are parallel and in the same direction. In
     * this case the point can be ignored, since the offset lines will also be
     * parallel.
     */
    if (numInt >= 2) {
      /**
       * segments are collinear but reversing.
       * Add an "end-cap" fillet
       * all the way around to other direction This case should ONLY happen
       * for LineStrings, so the orientation is always CW. (Polygons can never
       * have two consecutive segments which are parallel but reversed,
       * because that would be a self intersection.
       *
       */
      if (bufParams.getJoinStyle() == BufferParameters.JOIN_BEVEL ||
          bufParams.getJoinStyle() == BufferParameters.JOIN_MITRE) {
        if (addStartPoint) segList.addPt(offset0.p1);
        segList.addPt(offset1.p0);
      } else {
        addCornerFillet(
            s1, offset0.p1, offset1.p0, Orientation.CLOCKWISE, distance);
      }
    }
  }

  /**
   * Adds the offset points for an outside (convex) turn
   *
   * @param orientation
   * @param addStartPoint
   */
  void addOutsideTurn(int orientation, bool addStartPoint) {
    /**
     * Heuristic: If offset endpoints are very close together,
     * just use one of them as the corner vertex.
     * This avoids problems with computing mitre corners in the case
     * where the two segments are almost parallel
     * (which is hard to compute a robust intersection for).
     */
    if (offset0.p1.distance(offset1.p0) <
        distance * OFFSET_SEGMENT_SEPARATION_FACTOR) {
      segList.addPt(offset0.p1);
      return;
    }

    if (bufParams.getJoinStyle() == BufferParameters.JOIN_MITRE) {
      addMitreJoin(s1, offset0, offset1, distance);
    } else if (bufParams.getJoinStyle() == BufferParameters.JOIN_BEVEL) {
      addBevelJoin(offset0, offset1);
    } else {
      // add a circular fillet connecting the endpoints of the offset segments
      if (addStartPoint) segList.addPt(offset0.p1);
      // TESTING - comment out to produce beveled joins
      addCornerFillet(s1, offset0.p1, offset1.p0, orientation, distance);
      segList.addPt(offset1.p0);
    }
  }

  /**
   * Adds the offset points for an inside (concave) turn.
   *
   * @param orientation
   * @param addStartPoint
   */
  void addInsideTurn(int orientation, bool addStartPoint) {
    /**
     * add intersection point of offset segments (if any)
     */
    li!.computeIntersection(offset0.p0, offset0.p1, offset1.p0, offset1.p1);
    if (li!.hasIntersection()) {
      segList.addPt(li!.getIntersection(0));
    } else {
      /**
       * If no intersection is detected,
       * it means the angle is so small and/or the offset so
       * large that the offsets segments don't intersect.
       * In this case we must
       * add a "closing segment" to make sure the buffer curve is continuous,
       * fairly smooth (e.g. no sharp reversals in direction)
       * and tracks the buffer correctly around the corner. The curve connects
       * the endpoints of the segment offsets to points
       * which lie toward the centre point of the corner.
       * The joining curve will not appear in the final buffer outline, since it
       * is completely internal to the buffer polygon.
       *
       * In complex buffer cases the closing segment may cut across many other
       * segments in the generated offset curve.  In order to improve the
       * performance of the noding, the closing segment should be kept as short as possible.
       * (But not too short, since that would defeat its purpose).
       * This is the purpose of the closingSegFactor heuristic value.
       */

      /**
       * The intersection test above is vulnerable to robustness errors; i.e. it
       * may be that the offsets should intersect very close to their endpoints,
       * but aren't reported as such due to rounding. To handle this situation
       * appropriately, we use the following test: If the offset points are very
       * close, don't add closing segments but simply use one of the offset
       * points
       */
      _hasNarrowConcaveAngle = true;
      //System.out.println("NARROW ANGLE - distance = " + distance);
      if (offset0.p1.distance(offset1.p0) <
          distance * INSIDE_TURN_VERTEX_SNAP_DISTANCE_FACTOR) {
        segList.addPt(offset0.p1);
      } else {
        // add endpoint of this segment offset
        segList.addPt(offset0.p1);

        /**
         * Add "closing segment" of required length.
         */
        if (closingSegLengthFactor > 0) {
          Coordinate mid0 = new Coordinate(
              (closingSegLengthFactor * offset0.p1.x + s1.x) /
                  (closingSegLengthFactor + 1),
              (closingSegLengthFactor * offset0.p1.y + s1.y) /
                  (closingSegLengthFactor + 1));
          segList.addPt(mid0);
          Coordinate mid1 = new Coordinate(
              (closingSegLengthFactor * offset1.p0.x + s1.x) /
                  (closingSegLengthFactor + 1),
              (closingSegLengthFactor * offset1.p0.y + s1.y) /
                  (closingSegLengthFactor + 1));
          segList.addPt(mid1);
        } else {
          /**
           * This branch is not expected to be used except for testing purposes.
           * It is equivalent to the JTS 1.9 logic for closing segments
           * (which results in very poor performance for large buffer distances)
           */
          segList.addPt(s1);
        }

        //*/
        // add start point of next segment offset
        segList.addPt(offset1.p0);
      }
    }
  }

  /**
   * Compute an offset segment for an input segment on a given side and at a given distance.
   * The offset points are computed in full double precision, for accuracy.
   *
   * @param seg the segment to offset
   * @param side the side of the segment ({@link Position}) the offset lies on
   * @param distance the offset distance
   * @param offset the points computed for the offset segment
   */
  void computeOffsetSegment(
      LineSegment seg, int side, double distance, LineSegment offset) {
    int sideSign = side == Position.LEFT ? 1 : -1;
    double dx = seg.p1.x - seg.p0.x;
    double dy = seg.p1.y - seg.p0.y;
    double len = math.sqrt(dx * dx + dy * dy);
    // u is the vector that is the length of the offset, in the direction of the segment
    double ux = sideSign * distance * dx / len;
    double uy = sideSign * distance * dy / len;
    offset.p0.x = seg.p0.x - uy;
    offset.p0.y = seg.p0.y + ux;
    offset.p1.x = seg.p1.x - uy;
    offset.p1.y = seg.p1.y + ux;
  }

  /**
   * Add an end cap around point p1, terminating a line segment coming from p0
   */
  void addLineEndCap(Coordinate p0, Coordinate p1) {
    LineSegment seg = new LineSegment.fromCoordinates(p0, p1);

    LineSegment offsetL = new LineSegment.empty();
    computeOffsetSegment(seg, Position.LEFT, distance, offsetL);
    LineSegment offsetR = new LineSegment.empty();
    computeOffsetSegment(seg, Position.RIGHT, distance, offsetR);

    double dx = p1.x - p0.x;
    double dy = p1.y - p0.y;
    double angle = math.atan2(dy, dx);

    switch (bufParams.getEndCapStyle()) {
      case BufferParameters.CAP_ROUND:
        // add offset seg points with a fillet between them
        segList.addPt(offsetL.p1);
        addDirectedFillet(p1, angle + math.pi / 2, angle - math.pi / 2,
            Orientation.CLOCKWISE, distance);
        segList.addPt(offsetR.p1);
        break;
      case BufferParameters.CAP_FLAT:
        // only offset segment points are added
        segList.addPt(offsetL.p1);
        segList.addPt(offsetR.p1);
        break;
      case BufferParameters.CAP_SQUARE:
        // add a square defined by extensions of the offset segment endpoints
        Coordinate squareCapSideOffset = new Coordinate.empty2D();
        squareCapSideOffset.x = distance.abs() * math.cos(angle);
        squareCapSideOffset.y = distance.abs() * math.sin(angle);

        Coordinate squareCapLOffset = new Coordinate(
            offsetL.p1.x + squareCapSideOffset.x,
            offsetL.p1.y + squareCapSideOffset.y);
        Coordinate squareCapROffset = new Coordinate(
            offsetR.p1.x + squareCapSideOffset.x,
            offsetR.p1.y + squareCapSideOffset.y);
        segList.addPt(squareCapLOffset);
        segList.addPt(squareCapROffset);
        break;
    }
  }

  /**
   * Adds a mitre join connecting the two reflex offset segments.
   * The mitre will be beveled if it exceeds the mitre ratio limit.
   *
   * @param offset0 the first offset segment
   * @param offset1 the second offset segment
   * @param distance the offset distance
   */
  void addMitreJoin(
      Coordinate p, LineSegment offset0, LineSegment offset1, double distance) {
    /**
     * This computation is unstable if the offset segments are nearly collinear.
     * However, this situation should have been eliminated earlier by the check
     * for whether the offset segment endpoints are almost coincident
     */
    Coordinate? intPt = Intersection.intersection(
        offset0.p0, offset0.p1, offset1.p0, offset1.p1);
    if (intPt != null) {
      double mitreRatio =
          distance <= 0.0 ? 1.0 : intPt.distance(p) / distance.abs();
      if (mitreRatio <= bufParams.getMitreLimit()) {
        segList.addPt(intPt);
        return;
      }
    }
    // at this point either intersection failed or mitre limit was exceeded
    addLimitedMitreJoin(offset0, offset1, distance, bufParams.getMitreLimit());
//      addBevelJoin(offset0, offset1);
  }

  /**
   * Adds a limited mitre join connecting the two reflex offset segments.
   * A limited mitre is a mitre which is beveled at the distance
   * determined by the mitre ratio limit.
   *
   * @param offset0 the first offset segment
   * @param offset1 the second offset segment
   * @param distance the offset distance
   * @param mitreLimit the mitre limit ratio
   */
  void addLimitedMitreJoin(LineSegment offset0, LineSegment offset1,
      double distance, double mitreLimit) {
    Coordinate basePt = seg0.p1;

    double ang0 = Angle.angle2C(basePt, seg0.p0);

    // oriented angle between segments
    double angDiff = Angle.angleBetweenOriented(seg0.p0, basePt, seg1.p1);
    // half of the interior angle
    double angDiffHalf = angDiff / 2;

    // angle for bisector of the interior angle between the segments
    double midAng = Angle.normalize(ang0 + angDiffHalf);
    // rotating this by PI gives the bisector of the reflex angle
    double mitreMidAng = Angle.normalize(midAng + math.pi);

    // the miterLimit determines the distance to the mitre bevel
    double mitreDist = mitreLimit * distance;
    // the bevel delta is the difference between the buffer distance
    // and half of the length of the bevel segment
    double bevelDelta = mitreDist * math.sin(angDiffHalf).abs();
    double bevelHalfLen = distance - bevelDelta;

    // compute the midpoint of the bevel segment
    double bevelMidX = basePt.x + mitreDist * math.cos(mitreMidAng);
    double bevelMidY = basePt.y + mitreDist * math.sin(mitreMidAng);
    Coordinate bevelMidPt = new Coordinate(bevelMidX, bevelMidY);

    // compute the mitre midline segment from the corner point to the bevel segment midpoint
    LineSegment mitreMidLine =
        new LineSegment.fromCoordinates(basePt, bevelMidPt);

    // finally the bevel segment endpoints are computed as offsets from
    // the mitre midline
    Coordinate bevelEndLeft = mitreMidLine.pointAlongOffset(1.0, bevelHalfLen);
    Coordinate bevelEndRight =
        mitreMidLine.pointAlongOffset(1.0, -bevelHalfLen);

    if (side == Position.LEFT) {
      segList.addPt(bevelEndLeft);
      segList.addPt(bevelEndRight);
    } else {
      segList.addPt(bevelEndRight);
      segList.addPt(bevelEndLeft);
    }
  }

  /**
   * Adds a bevel join connecting the two offset segments
   * around a reflex corner.
   *
   * @param offset0 the first offset segment
   * @param offset1 the second offset segment
   */
  void addBevelJoin(LineSegment offset0, LineSegment offset1) {
    segList.addPt(offset0.p1);
    segList.addPt(offset1.p0);
  }

  /**
   * Add points for a circular fillet around a reflex corner.
   * Adds the start and end points
   *
   * @param p base point of curve
   * @param p0 start point of fillet curve
   * @param p1 endpoint of fillet curve
   * @param direction the orientation of the fillet
   * @param radius the radius of the fillet
   */
  void addCornerFillet(Coordinate p, Coordinate p0, Coordinate p1,
      int direction, double radius) {
    double dx0 = p0.x - p.x;
    double dy0 = p0.y - p.y;
    double startAngle = math.atan2(dy0, dx0);
    double dx1 = p1.x - p.x;
    double dy1 = p1.y - p.y;
    double endAngle = math.atan2(dy1, dx1);

    if (direction == Orientation.CLOCKWISE) {
      if (startAngle <= endAngle) startAngle += 2.0 * math.pi;
    } else {
      // direction == COUNTERCLOCKWISE
      if (startAngle >= endAngle) startAngle -= 2.0 * math.pi;
    }
    segList.addPt(p0);
    addDirectedFillet(p, startAngle, endAngle, direction, radius);
    segList.addPt(p1);
  }

  /**
   * Adds points for a circular fillet arc
   * between two specified angles.
   * The start and end point for the fillet are not added -
   * the caller must add them if required.
   *
   * @param direction is -1 for a CW angle, 1 for a CCW angle
   * @param radius the radius of the fillet
   */
  void addDirectedFillet(Coordinate p, double startAngle, double endAngle,
      int direction, double radius) {
    int directionFactor = direction == Orientation.CLOCKWISE ? -1 : 1;

    double totalAngle = (startAngle - endAngle).abs();
    int nSegs = (totalAngle / filletAngleQuantum + 0.5).toInt();

    if (nSegs < 1)
      return; // no segments because angle is less than increment - nothing to do!

    double initAngle, currAngleInc;

    // choose angle increment so that each segment has equal length
    initAngle = 0.0;
    currAngleInc = totalAngle / nSegs;

    double currAngle = initAngle;
    Coordinate pt = new Coordinate.empty2D();
    while (currAngle < totalAngle) {
      double angle = startAngle + directionFactor * currAngle;
      pt.x = p.x + radius * math.cos(angle);
      pt.y = p.y + radius * math.sin(angle);
      segList.addPt(pt);
      currAngle += currAngleInc;
    }
  }

  /**
   * Creates a CW circle around a point
   */
  void createCircle(Coordinate p) {
    // add start point
    Coordinate pt = new Coordinate(p.x + distance, p.y);
    segList.addPt(pt);
    addDirectedFillet(p, 0.0, 2.0 * math.pi, -1, distance);
    segList.closeRing();
  }

  /**
   * Creates a CW square around a point
   */
  void createSquare(Coordinate p) {
    segList.addPt(new Coordinate(p.x + distance, p.y + distance));
    segList.addPt(new Coordinate(p.x + distance, p.y - distance));
    segList.addPt(new Coordinate(p.x - distance, p.y - distance));
    segList.addPt(new Coordinate(p.x - distance, p.y + distance));
    segList.closeRing();
  }
}

/**
 * A dynamic list of the vertices in a constructed offset curve.
 * Automatically removes adjacent vertices
 * which are closer than a given tolerance.
 *
 * @author Martin Davis
 *
 */
class OffsetSegmentString {
  late List ptList;
  PrecisionModel? precisionModel = null;

  /**
   * The distance below which two adjacent points on the curve
   * are considered to be coincident.
   * This is chosen to be a small fraction of the offset distance.
   */
  double minimimVertexDistance = 0.0;

  OffsetSegmentString() {
    ptList = [];
  }

  void setPrecisionModel(PrecisionModel precisionModel) {
    this.precisionModel = precisionModel;
  }

  void setMinimumVertexDistance(double minimimVertexDistance) {
    this.minimimVertexDistance = minimimVertexDistance;
  }

  void addPt(Coordinate pt) {
    Coordinate bufPt = new Coordinate.fromCoordinate(pt);
    precisionModel!.makeCoordinatePrecise(bufPt);
    // don't add duplicate (or near-duplicate) points
    if (isRedundant(bufPt)) return;
    ptList.add(bufPt);
//System.out.println(bufPt);
  }

  void addPts(List<Coordinate> pt, bool isForward) {
    if (isForward) {
      for (int i = 0; i < pt.length; i++) {
        addPt(pt[i]);
      }
    } else {
      for (int i = pt.length - 1; i >= 0; i--) {
        addPt(pt[i]);
      }
    }
  }

  /**
   * Tests whether the given point is redundant
   * relative to the previous
   * point in the list (up to tolerance).
   *
   * @param pt
   * @return true if the point is redundant
   */
  bool isRedundant(Coordinate pt) {
    if (ptList.length < 1) return false;
    Coordinate lastPt = ptList.last;
    double ptDist = pt.distance(lastPt);
    if (ptDist < minimimVertexDistance) return true;
    return false;
  }

  void closeRing() {
    if (ptList.length < 1) return;
    Coordinate startPt = new Coordinate.fromCoordinate(ptList.first);
    Coordinate lastPt = ptList.last;
    if (startPt.equals(lastPt)) return;
    ptList.add(startPt);
  }

  void reverse() {}

  List<Coordinate> getCoordinates() {
    /*
     // check that points are a ring - add the startpoint again if they are not
   if (ptList.size() > 1) {
      Coordinate start  = (Coordinate) ptList.get(0);
      Coordinate end    = (Coordinate) ptList.get(ptList.size() - 1);
      if (! start.equals(end) ) addPt(start);
    }
    */
    List<Coordinate> coord = List.from(ptList);
    return coord;
  }

  String toString() {
    GeometryFactory fact = new GeometryFactory.defaultPrecision();
    LineString line = fact.createLineString(getCoordinates());
    return line.toString();
  }
}

/**
 * Simplifies a buffer input line to
 * remove concavities with shallow depth.
 * <p>
 * The most important benefit of doing this
 * is to reduce the number of points and the complexity of
 * shape which will be buffered.
 * It also reduces the risk of gores created by
 * the quantized fillet arcs (although this issue
 * should be eliminated in any case by the
 * offset curve generation logic).
 * <p>
 * A key aspect of the simplification is that it
 * affects inside (concave or inward) corners only.
 * Convex (outward) corners are preserved, since they
 * are required to ensure that the generated buffer curve
 * lies at the correct distance from the input geometry.
 * <p>
 * Another important heuristic used is that the end segments
 * of the input are never simplified.  This ensures that
 * the client buffer code is able to generate end caps faithfully.
 * <p>
 * No attempt is made to avoid self-intersections in the output.
 * This is acceptable for use for generating a buffer offset curve,
 * since the buffer algorithm is insensitive to invalid polygonal
 * geometry.  However,
 * this means that this algorithm
 * cannot be used as a general-purpose polygon simplification technique.
 *
 * @author Martin Davis
 *
 */
class BufferInputLineSimplifier {
  /**
   * Simplify the input coordinate list.
   * If the distance tolerance is positive,
   * concavities on the LEFT side of the line are simplified.
   * If the supplied distance tolerance is negative,
   * concavities on the RIGHT side of the line are simplified.
   *
   * @param inputLine the coordinate list to simplify
   * @param distanceTol simplification distance tolerance to use
   * @return the simplified coordinate list
   */
  static List<Coordinate> simplify(
      List<Coordinate> inputLine, double distanceTol) {
    BufferInputLineSimplifier simp = new BufferInputLineSimplifier(inputLine);
    return simp.simplifyWithTol(distanceTol);
  }

  static final int INIT = 0;
  static final int DELETE = 1;
  static final int KEEP = 1;

  List<Coordinate> inputLine;
  double distanceTol = 0.0;
  late List<int> isDeleted;
  int angleOrientation = Orientation.COUNTERCLOCKWISE;

  BufferInputLineSimplifier(this.inputLine);

  /**
   * Simplify the input coordinate list.
   * If the distance tolerance is positive,
   * concavities on the LEFT side of the line are simplified.
   * If the supplied distance tolerance is negative,
   * concavities on the RIGHT side of the line are simplified.
   *
   * @param distanceTol simplification distance tolerance to use
   * @return the simplified coordinate list
   */
  List<Coordinate> simplifyWithTol(double distanceTol) {
    this.distanceTol = distanceTol.abs();
    if (distanceTol < 0) angleOrientation = Orientation.CLOCKWISE;

    // rely on fact that bool array is filled with false value
    isDeleted = List.filled(inputLine.length, 0);

    bool isChanged = false;
    do {
      isChanged = deleteShallowConcavities();
    } while (isChanged);

    return collapseLine();
  }

  /**
   * Uses a sliding window containing 3 vertices to detect shallow angles
   * in which the middle vertex can be deleted, since it does not
   * affect the shape of the resulting buffer in a significant way.
   * @return
   */
  bool deleteShallowConcavities() {
    /**
     * Do not simplify end line segments of the line string.
     * This ensures that end caps are generated consistently.
     */
    int index = 1;

    int midIndex = findNextNonDeletedIndex(index);
    int lastIndex = findNextNonDeletedIndex(midIndex);

    bool isChanged = false;
    while (lastIndex < inputLine.length) {
      // test triple for shallow concavity
      bool isMiddleVertexDeleted = false;
      if (isDeletable(index, midIndex, lastIndex, distanceTol)) {
        isDeleted[midIndex] = DELETE;
        isMiddleVertexDeleted = true;
        isChanged = true;
      }
      // move simplification window forward
      if (isMiddleVertexDeleted)
        index = lastIndex;
      else
        index = midIndex;

      midIndex = findNextNonDeletedIndex(index);
      lastIndex = findNextNonDeletedIndex(midIndex);
    }
    return isChanged;
  }

  /**
   * Finds the next non-deleted index, or the end of the point array if none
   * @param index
   * @return the next non-deleted index, if any
   * or inputLine.length if there are no more non-deleted indices
   */
  int findNextNonDeletedIndex(int index) {
    int next = index + 1;
    while (next < inputLine.length && isDeleted[next] == DELETE) next++;
    return next;
  }

  List<Coordinate> collapseLine() {
    CoordinateList coordList = new CoordinateList();
    for (int i = 0; i < inputLine.length; i++) {
      if (isDeleted[i] != DELETE) coordList.add(inputLine[i]);
    }
//    if (coordList.size() < inputLine.length)      System.out.println("Simplified " + (inputLine.length - coordList.size()) + " pts");
    return coordList.toCoordinateArray();
  }

  bool isDeletable(int i0, int i1, int i2, double distanceTol) {
    Coordinate p0 = inputLine[i0];
    Coordinate p1 = inputLine[i1];
    Coordinate p2 = inputLine[i2];

    if (!isConcave(p0, p1, p2)) return false;
    if (!isShallow(p0, p1, p2, distanceTol)) return false;

    // MD - don't use this heuristic - it's too restricting
//  	if (p0.distance(p2) > distanceTol) return false;

    return isShallowSampled(p0, p1, i0, i2, distanceTol);
  }

  bool isShallowConcavity(
      Coordinate p0, Coordinate p1, Coordinate p2, double distanceTol) {
    int orientation = Orientation.index(p0, p1, p2);
    bool isAngleToSimplify = (orientation == angleOrientation);
    if (!isAngleToSimplify) return false;

    double dist = Distance.pointToSegment(p1, p0, p2);
    return dist < distanceTol;
  }

  static final int NUM_PTS_TO_CHECK = 10;

  /**
   * Checks for shallowness over a sample of points in the given section.
   * This helps prevents the simplification from incrementally
   * "skipping" over points which are in fact non-shallow.
   *
   * @param p0 start coordinate of section
   * @param p2 end coordinate of section
   * @param i0 start index of section
   * @param i2 end index of section
   * @param distanceTol distance tolerance
   * @return
   */
  bool isShallowSampled(
      Coordinate p0, Coordinate p2, int i0, int i2, double distanceTol) {
    // check every n'th point to see if it is within tolerance
    int inc = ((i2 - i0) / NUM_PTS_TO_CHECK).toInt();
    if (inc <= 0) inc = 1;

    for (int i = i0; i < i2; i += inc) {
      if (!isShallow(p0, p2, inputLine[i], distanceTol)) return false;
    }
    return true;
  }

  bool isShallow(
      Coordinate p0, Coordinate p1, Coordinate p2, double distanceTol) {
    double dist = Distance.pointToSegment(p1, p0, p2);
    return dist < distanceTol;
  }

  bool isConcave(Coordinate p0, Coordinate p1, Coordinate p2) {
    int orientation = Orientation.index(p0, p1, p2);
    bool isConcave = (orientation == angleOrientation);
    return isConcave;
  }
}

/**
 * Creates all the raw offset curves for a buffer of a {@link Geometry}.
 * Raw curves need to be noded together and polygonized to form the final buffer area.
 *
 * @version 1.7
 */
class OffsetCurveSetBuilder {
  Geometry inputGeom;
  double distance;
  OffsetCurveBuilder curveBuilder;

  List curveList = [];

  OffsetCurveSetBuilder(this.inputGeom, this.distance, this.curveBuilder);

  /**
   * Computes the set of raw offset curves for the buffer.
   * Each offset curve has an attached {@link Label} indicating
   * its left and right location.
   *
   * @return a Collection of SegmentStrings representing the raw buffer curves
   */
  List getCurves() {
    add(inputGeom);
    return curveList;
  }

  /**
   * Creates a {@link SegmentString} for a coordinate list which is a raw offset curve,
   * and adds it to the list of buffer curves.
   * The SegmentString is tagged with a Label giving the topology of the curve.
   * The curve may be oriented in either direction.
   * If the curve is oriented CW, the locations will be:
   * <br>Left: Location.EXTERIOR
   * <br>Right: Location.INTERIOR
   */
  void addCurve(List<Coordinate>? coord, int leftLoc, int rightLoc) {
    // don't add null or trivial curves
    if (coord == null || coord.length < 2) return;
    // add the edge for a coordinate list which is a raw offset curve
    SegmentString e = new NodedSegmentString(
        coord, new Label.args4(0, Location.BOUNDARY, leftLoc, rightLoc));
    curveList.add(e);
  }

  void add(Geometry g) {
    if (g.isEmpty()) return;

    if (g is Polygon)
      addPolygon(g);
    // LineString also handles LinearRings
    else if (g is LineString)
      addLineString(g);
    else if (g is Point)
      addPoint(g);
    else if (g is MultiPoint)
      addCollection(g);
    else if (g is MultiLineString)
      addCollection(g);
    else if (g is MultiPolygon)
      addCollection(g);
    else if (g is GeometryCollection)
      addCollection(g);
    else
      throw new UnsupportedError(g.runtimeType.toString());
  }

  void addCollection(GeometryCollection gc) {
    for (int i = 0; i < gc.getNumGeometries(); i++) {
      Geometry g = gc.getGeometryN(i);
      add(g);
    }
  }

  /**
   * Add a Point to the graph.
   */
  void addPoint(Point p) {
    // a zero or negative width buffer of a line/point is empty
    if (distance <= 0.0) return;
    List<Coordinate> coord = p.getCoordinates();
    List<Coordinate>? curve = curveBuilder.getLineCurve(coord, distance);
    addCurve(curve, Location.EXTERIOR, Location.INTERIOR);
  }

  void addLineString(LineString line) {
    // a zero or negative width buffer of a line/point is empty
    if (distance <= 0.0 && !curveBuilder.getBufferParameters().isSingleSided)
      return;
    List<Coordinate> coord =
        CoordinateArrays.removeRepeatedPoints(line.getCoordinates());

    /**
     * Rings (closed lines) are generated with a continuous curve, 
     * with no end arcs. This produces better quality linework, 
     * and avoids noding issues with arcs around almost-parallel end segments.
     * See JTS #523 and #518.
     * 
     * Singled-sided buffers currently treat rings as if they are lines.
     */
    if (CoordinateArrays.isRing(coord) &&
        !curveBuilder.getBufferParameters().isSingleSided) {
      addRingBothSides(coord, distance);
    } else {
      List<Coordinate>? curve = curveBuilder.getLineCurve(coord, distance);
      addCurve(curve, Location.EXTERIOR, Location.INTERIOR);
    }

    // TESTING
    //List<Coordinate> curveTrim = BufferCurveLoopPruner.prune(curve);
    //addCurve(curveTrim, Location.EXTERIOR, Location.INTERIOR);
  }

  void addPolygon(Polygon p) {
    double offsetDistance = distance;
    int offsetSide = Position.LEFT;
    if (distance < 0.0) {
      offsetDistance = -distance;
      offsetSide = Position.RIGHT;
    }

    LinearRing shell = p.getExteriorRing();
    List<Coordinate> shellCoord =
        CoordinateArrays.removeRepeatedPoints(shell.getCoordinates());
    // optimization - don't bother computing buffer
    // if the polygon would be completely eroded
    if (distance < 0.0 && isErodedCompletely(shell, distance)) return;
    // don't attempt to buffer a polygon with too few distinct vertices
    if (distance <= 0.0 && shellCoord.length < 3) return;

    addRingSide(shellCoord, offsetDistance, offsetSide, Location.EXTERIOR,
        Location.INTERIOR);

    for (int i = 0; i < p.getNumInteriorRing(); i++) {
      LinearRing hole = p.getInteriorRingN(i);
      List<Coordinate> holeCoord =
          CoordinateArrays.removeRepeatedPoints(hole.getCoordinates());

      // optimization - don't bother computing buffer for this hole
      // if the hole would be completely covered
      if (distance > 0.0 && isErodedCompletely(hole, -distance)) continue;

      // Holes are topologically labelled opposite to the shell, since
      // the interior of the polygon lies on their opposite side
      // (on the left, if the hole is oriented CCW)
      addRingSide(holeCoord, offsetDistance, Position.opposite(offsetSide),
          Location.INTERIOR, Location.EXTERIOR);
    }
  }

  void addRingBothSides(List<Coordinate> coord, double distance) {
    addRingSide(
        coord, distance, Position.LEFT, Location.EXTERIOR, Location.INTERIOR);
    /* Add the opposite side of the ring
    */
    addRingSide(
        coord, distance, Position.RIGHT, Location.INTERIOR, Location.EXTERIOR);
  }

  /**
   * Adds an offset curve for one side of a ring.
   * The side and left and right topological location arguments
   * are provided as if the ring is oriented CW.
   * (If the ring is in the opposite orientation,
   * this is detected and 
   * the left and right locations are interchanged and the side is flipped.)
   *
   * @param coord the coordinates of the ring (must not contain repeated points)
   * @param offsetDistance the positive distance at which to create the buffer
   * @param side the side {@link Position} of the ring on which to construct the buffer line
   * @param cwLeftLoc the location on the L side of the ring (if it is CW)
   * @param cwRightLoc the location on the R side of the ring (if it is CW)
   */
  void addRingSide(List<Coordinate> coord, double offsetDistance, int side,
      int cwLeftLoc, int cwRightLoc) {
    // don't bother adding ring if it is "flat" and will disappear in the output
    if (offsetDistance == 0.0 && coord.length < LinearRing.MINIMUM_VALID_SIZE)
      return;

    int leftLoc = cwLeftLoc;
    int rightLoc = cwRightLoc;
    if (coord.length >= LinearRing.MINIMUM_VALID_SIZE &&
        Orientation.isCCW(coord)) {
      leftLoc = cwRightLoc;
      rightLoc = cwLeftLoc;
      side = Position.opposite(side);
    }
    List<Coordinate>? curve =
        curveBuilder.getRingCurve(coord, side, offsetDistance);
    addCurve(curve, leftLoc, rightLoc);
  }

  /**
   * Adds an offset curve for a polygon ring.
   * The side and left and right topological location arguments
   * assume that the ring is oriented CW.
   * If the ring is in the opposite orientation,
   * the left and right locations must be interchanged and the side flipped.
   *
   * @param coord the coordinates of the ring (must not contain repeated points)
   * @param offsetDistance the distance at which to create the buffer
   * @param side the side of the ring on which to construct the buffer line
   * @param cwLeftLoc the location on the L side of the ring (if it is CW)
   * @param cwRightLoc the location on the R side of the ring (if it is CW)
   */
  void addPolygonRing(List<Coordinate> coord, double offsetDistance, int side,
      int cwLeftLoc, int cwRightLoc) {
    // don't bother adding ring if it is "flat" and will disappear in the output
    if (offsetDistance == 0.0 && coord.length < LinearRing.MINIMUM_VALID_SIZE)
      return;

    int leftLoc = cwLeftLoc;
    int rightLoc = cwRightLoc;
    if (coord.length >= LinearRing.MINIMUM_VALID_SIZE &&
        Orientation.isCCW(coord)) {
      leftLoc = cwRightLoc;
      rightLoc = cwLeftLoc;
      side = Position.opposite(side);
    }
    List<Coordinate>? curve =
        curveBuilder.getRingCurve(coord, side, offsetDistance);
    addCurve(curve, leftLoc, rightLoc);
  }

  /**
   * The ringCoord is assumed to contain no repeated points.
   * It may be degenerate (i.e. contain only 1, 2, or 3 points).
   * In this case it has no area, and hence has a minimum diameter of 0.
   *
   * @param ringCoord
   * @param offsetDistance
   * @return
   */
  bool isErodedCompletely(LinearRing ring, double bufferDistance) {
    List<Coordinate> ringCoord = ring.getCoordinates();
    // degenerate ring has no area
    if (ringCoord.length < 4) return bufferDistance < 0;

    // important test to eliminate inverted triangle bug
    // also optimizes erosion test for triangles
    if (ringCoord.length == 4)
      return isTriangleErodedCompletely(ringCoord, bufferDistance);

    // if envelope is narrower than twice the buffer distance, ring is eroded
    Envelope env = ring.getEnvelopeInternal();
    double envMinDimension = math.min(env.getHeight(), env.getWidth());
    if (bufferDistance < 0.0 && 2 * bufferDistance.abs() > envMinDimension)
      return true;

    return false;
    /**
     * The following is a heuristic test to determine whether an
     * inside buffer will be eroded completely.
     * It is based on the fact that the minimum diameter of the ring pointset
     * provides an upper bound on the buffer distance which would erode the
     * ring.
     * If the buffer distance is less than the minimum diameter, the ring
     * may still be eroded, but this will be determined by
     * a full topological computation.
     *
     */
//System.out.println(ring);
/* MD  7 Feb 2005 - there's an unknown bug in the MD code, so disable this for now
    MinimumDiameter md = new MinimumDiameter(ring);
    minDiam = md.getLength();
    //System.out.println(md.getDiameter());
    return minDiam < 2 * Math.abs(bufferDistance);
    */
  }

  /**
   * Tests whether a triangular ring would be eroded completely by the given
   * buffer distance.
   * This is a precise test.  It uses the fact that the inner buffer of a
   * triangle converges on the inCentre of the triangle (the point
   * equidistant from all sides).  If the buffer distance is greater than the
   * distance of the inCentre from a side, the triangle will be eroded completely.
   *
   * This test is important, since it removes a problematic case where
   * the buffer distance is slightly larger than the inCentre distance.
   * In this case the triangle buffer curve "inverts" with incorrect topology,
   * producing an incorrect hole in the buffer.
   *
   * @param triangleCoord
   * @param bufferDistance
   * @return
   */
  bool isTriangleErodedCompletely(
      List<Coordinate> triangleCoord, double bufferDistance) {
    Triangle tri =
        new Triangle(triangleCoord[0], triangleCoord[1], triangleCoord[2]);
    Coordinate inCentre = tri.inCentre();
    double distToCentre = Distance.pointToSegment(inCentre, tri.p0, tri.p1);
    return distToCentre < bufferDistance.abs();
  }
}
