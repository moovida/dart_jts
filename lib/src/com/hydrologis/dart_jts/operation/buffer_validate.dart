part of dart_jts;

/**
 * Validates that the result of a buffer operation
 * is geometrically correct, within a computed tolerance.
 * <p>
 * This is a heuristic test, and may return false positive results
 * (I.e. it may fail to detect an invalid result.)
 * It should never return a false negative result, however
 * (I.e. it should never report a valid result as invalid.)
 * <p>
 * This test may be (much) more expensive than the original
 * buffer computation.
 *
 * @author Martin Davis
 */
class BufferResultValidator {
  static bool VERBOSE = false;

  /**
   * Maximum allowable fraction of buffer distance the
   * actual distance can differ by.
   * 1% sometimes causes an error - 1.2% should be safe.
   */
  static final double MAX_ENV_DIFF_FRAC = .012;

  static bool isValidGDG(Geometry g, double distance, Geometry result) {
    BufferResultValidator validator =
        new BufferResultValidator(g, distance, result);
    if (validator.isValid()) return true;
    return false;
  }

  /**
   * Checks whether the geometry buffer is valid,
   * and returns an error message if not.
   *
   * @param g
   * @param distance
   * @param result
   * @return an appropriate error message
   * or null if the buffer is valid
   */
  static String? isValidMsg(Geometry g, double distance, Geometry result) {
    BufferResultValidator validator =
        new BufferResultValidator(g, distance, result);
    if (!validator.isValid()) return validator.getErrorMessage();
    return null;
  }

  Geometry input;
  double distance = 0.0;
  Geometry result;
  bool _isValid = true;
  String? errorMsg = null;
  Coordinate? errorLocation = null;
  Geometry? errorIndicator = null;

  BufferResultValidator(this.input, this.distance, this.result);

  bool isValid() {
    checkPolygonal();
    if (!_isValid) return _isValid;
    checkExpectedEmpty();
    if (!_isValid) return _isValid;
    checkEnvelope();
    if (!_isValid) return _isValid;
    checkArea();
    if (!_isValid) return _isValid;
    checkDistance();
    return _isValid;
  }

  String? getErrorMessage() {
    return errorMsg;
  }

  Coordinate? getErrorLocation() {
    return errorLocation;
  }

  /**
   * Gets a geometry which indicates the location and nature of a validation failure.
   * <p>
   * If the failure is due to the buffer curve being too far or too close
   * to the input, the indicator is a line segment showing the location and size
   * of the discrepancy.
   *
   * @return a geometric error indicator
   * or null if no error was found
   */
  Geometry? getErrorIndicator() {
    return errorIndicator;
  }

  void report(String checkName) {
    if (!VERBOSE) return;
    print("Check " + checkName + ": " + (_isValid ? "passed" : "FAILED"));
  }

  void checkPolygonal() {
    if (!(result is Polygon || result is MultiPolygon)) _isValid = false;
    errorMsg = "Result is not polygonal";
    errorIndicator = result;
    report("Polygonal");
  }

  void checkExpectedEmpty() {
    // can't check areal features
    if (input.getDimension() >= 2) return;
    // can't check positive distances
    if (distance > 0.0) return;

    // at this point can expect an empty result
    if (!result.isEmpty()) {
      _isValid = false;
      errorMsg = "Result is non-empty";
      errorIndicator = result;
    }
    report("ExpectedEmpty");
  }

  void checkEnvelope() {
    if (distance < 0.0) return;

    double padding = distance * MAX_ENV_DIFF_FRAC;
    if (padding == 0.0) padding = 0.001;

    Envelope expectedEnv =
        new Envelope.fromEnvelope(input.getEnvelopeInternal());
    expectedEnv.expandByDistance(distance);

    Envelope bufEnv = new Envelope.fromEnvelope(result.getEnvelopeInternal());
    bufEnv.expandByDistance(padding);

    if (!bufEnv.containsEnvelope(expectedEnv)) {
      _isValid = false;
      errorMsg = "Buffer envelope is incorrect";
      errorIndicator = input.getFactory().toGeometry(bufEnv);
    }
    report("Envelope");
  }

  void checkArea() {
    double inputArea = input.getArea();
    double resultArea = result.getArea();

    if (distance > 0.0 && inputArea > resultArea) {
      _isValid = false;
      errorMsg = "Area of positive buffer is smaller than input";
      errorIndicator = result;
    }
    if (distance < 0.0 && inputArea < resultArea) {
      _isValid = false;
      errorMsg = "Area of negative buffer is larger than input";
      errorIndicator = result;
    }
    report("Area");
  }

  void checkDistance() {
    BufferDistanceValidator distValid =
        new BufferDistanceValidator(input, distance, result);
    if (!distValid.isValid()) {
      _isValid = false;
      errorMsg = distValid.getErrorMessage();
      errorLocation = distValid.getErrorLocation();
      errorIndicator = distValid.getErrorIndicator();
    }
    report("Distance");
  }
}

/**
 * Validates that a given buffer curve lies an appropriate distance
 * from the input generating it.
 * Useful only for round buffers (cap and join).
 * Can be used for either positive or negative distances.
 * <p>
 * This is a heuristic test, and may return false positive results
 * (I.e. it may fail to detect an invalid result.)
 * It should never return a false negative result, however
 * (I.e. it should never report a valid result as invalid.)
 *
 * @author mbdavis
 *
 */
class BufferDistanceValidator {
  static bool VERBOSE = false;

  /**
   * Maximum allowable fraction of buffer distance the
   * actual distance can differ by.
   * 1% sometimes causes an error - 1.2% should be safe.
   */
  static final double MAX_DISTANCE_DIFF_FRAC = .012;

  Geometry input;
  double bufDistance;
  Geometry result;

  double minValidDistance = 0.0;
  double maxValidDistance = 0.0;

  double minDistanceFound = 0.0;
  double maxDistanceFound = 0.0;

  bool _isValid = true;
  String? errMsg = null;
  Coordinate? errorLocation = null;
  Geometry? errorIndicator = null;

  BufferDistanceValidator(this.input, this.bufDistance, this.result);

  bool isValid() {
    double posDistance = bufDistance.abs();
    double distDelta = MAX_DISTANCE_DIFF_FRAC * posDistance;
    minValidDistance = posDistance - distDelta;
    maxValidDistance = posDistance + distDelta;

    // can't use this test if either is empty
    if (input.isEmpty() || result.isEmpty()) return true;

    if (bufDistance > 0.0) {
      checkPositiveValid();
    } else {
      checkNegativeValid();
    }
    if (VERBOSE) {
      print(
          "Min Dist= $minDistanceFound  err= ${(1.0 - minDistanceFound / bufDistance)}  Max Dist= $maxDistanceFound  err= ${(maxDistanceFound / bufDistance - 1.0)}");
    }
    return _isValid;
  }

  String? getErrorMessage() {
    return errMsg;
  }

  Coordinate? getErrorLocation() {
    return errorLocation;
  }

  /**
   * Gets a geometry which indicates the location and nature of a validation failure.
   * <p>
   * The indicator is a line segment showing the location and size
   * of the distance discrepancy.
   *
   * @return a geometric error indicator
   * or null if no error was found
   */
  Geometry? getErrorIndicator() {
    return errorIndicator;
  }

  void checkPositiveValid() {
    Geometry bufCurve = result.getBoundary();
    checkMinimumDistance(input, bufCurve, minValidDistance);
    if (!_isValid) return;

    checkMaximumDistance(input, bufCurve, maxValidDistance);
  }

  void checkNegativeValid() {
    // Assert: only polygonal inputs can be checked for negative buffers

    // MD - could generalize this to handle GCs too
    if (!(input is Polygon ||
        input is MultiPolygon ||
        input is GeometryCollection)) {
      return;
    }
    Geometry inputCurve = getPolygonLines(input);
    checkMinimumDistance(inputCurve, result, minValidDistance);
    if (!_isValid) return;

    checkMaximumDistance(inputCurve, result, maxValidDistance);
  }

  Geometry getPolygonLines(Geometry g) {
    List<LineString> lines = [];
    LinearComponentExtracter lineExtracter =
        new LinearComponentExtracter(lines);
    List polys = PolygonExtracter.getPolygons(g);
    for (Polygon poly in polys) {
      poly.applyGCF(lineExtracter);
    }
    return g.getFactory().buildGeometry(lines);
  }

  /**
   * Checks that two geometries are at least a minimum distance apart.
   *
   * @param g1 a geometry
   * @param g2 a geometry
   * @param minDist the minimum distance the geometries should be separated by
   */
  void checkMinimumDistance(Geometry g1, Geometry g2, double minDist) {
    DistanceOp distOp = new DistanceOp.withTerminateDistance(g1, g2, minDist);
    minDistanceFound = distOp.distance();

    if (minDistanceFound < minDist) {
      _isValid = false;
      List<Coordinate> pts = distOp.nearestPoints();
      errorLocation = distOp.nearestPoints()[1];
      errorIndicator = g1.getFactory().createLineString(pts);
      errMsg = "Distance between buffer curve and input is too small " +
          "($minDistanceFound at " +
          WKTWriter.toLineStringFromCoords(pts[0], pts[1]) +
          " )";
    }
  }

  /**
   * Checks that the furthest distance from the buffer curve to the input
   * is less than the given maximum distance.
   * This uses the Oriented Hausdorff distance metric.
   * It corresponds to finding
   * the point on the buffer curve which is furthest from <i>some</i> point on the input.
   *
   * @param input a geometry
   * @param bufCurve a geometry
   * @param maxDist the maximum distance that a buffer result can be from the input
   */
  void checkMaximumDistance(Geometry input, Geometry bufCurve, double maxDist) {
//    BufferCurveMaximumDistanceFinder maxDistFinder = new BufferCurveMaximumDistanceFinder(input);
//    maxDistanceFound = maxDistFinder.findDistance(bufCurve);

    DiscreteHausdorffDistance haus =
        new DiscreteHausdorffDistance(bufCurve, input);
    haus.setDensifyFraction(0.25);
    maxDistanceFound = haus.orientedDistance();

    if (maxDistanceFound > maxDist) {
      _isValid = false;
      List<Coordinate> pts = haus.getCoordinates();
      errorLocation = pts[1];
      errorIndicator = input.getFactory().createLineString(pts);
      errMsg = "Distance between buffer curve and input is too large " +
          "($maxDistanceFound at " +
          WKTWriter.toLineStringFromCoords(pts[0], pts[1]) +
          ")";
    }
  }

/*
   void OLDcheckMaximumDistance(Geometry input, Geometry bufCurve, double maxDist)
  {
    BufferCurveMaximumDistanceFinder maxDistFinder = new BufferCurveMaximumDistanceFinder(input);
    maxDistanceFound = maxDistFinder.findDistance(bufCurve);
    
    
    if (maxDistanceFound > maxDist) {
      isValid = false;
      PointPairDistance ptPairDist = maxDistFinder.getDistancePoints();
      errorLocation = ptPairDist.getCoordinate(1);
      errMsg = "Distance between buffer curve and input is too large "
        + "(" + ptPairDist.getDistance()
        + " at " + ptPairDist.toString() +")";
    }
  }
  */

}
