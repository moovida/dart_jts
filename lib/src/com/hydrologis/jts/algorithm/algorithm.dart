part of dart_jts;

/// Implements basic computational geometry algorithms using {@link DD} arithmetic.
///
/// @author Martin Davis
///
class CGAlgorithmsDD {
  CGAlgorithmsDD() {}

  /// Returns the index of the direction of the point <code>q</code> relative to
  /// a vector specified by <code>p1-p2</code>.
  ///
  /// @param p1 the origin point of the vector
  /// @param p2 the final point of the vector
  /// @param q the point to compute the direction to
  ///
  /// @return 1 if q is counter-clockwise (left) from p1-p2
  /// @return -1 if q is clockwise (right) from p1-p2
  /// @return 0 if q is collinear with p1-p2
  static int orientationIndex(Coordinate p1, Coordinate p2, Coordinate q) {
    // fast filter for orientation index
    // avoids use of slow extended-precision arithmetic in many cases
    int index = orientationIndexFilter(p1, p2, q);
    if (index <= 1) return index;

    // normalize coordinates
    DD dx1 = DD.valueOf(p2.x).selfAdd(-p1.x);
    DD dy1 = DD.valueOf(p2.y).selfAdd(-p1.y);
    DD dx2 = DD.valueOf(q.x).selfAdd(-p2.x);
    DD dy2 = DD.valueOf(q.y).selfAdd(-p2.y);

    // sign of determinant - unrolled for performance
    return dx1.selfMultiplyDD(dy2).selfSubtractDD(dy1.selfMultiplyDD(dx2)).signum();
  }

  /// Computes the sign of the determinant of the 2x2 matrix
  /// with the given entries.
  ///
  /// @return -1 if the determinant is negative,
  /// @return  1 if the determinant is positive,
  /// @return  0 if the determinant is 0.
  static int signOfDet2x2DD(DD x1, DD y1, DD x2, DD y2) {
    DD det = x1.multiplyDD(y2).selfSubtractDD(y1.multiplyDD(x2));
    return det.signum();
  }

  /// Computes the sign of the determinant of the 2x2 matrix
  /// with the given entries.
  ///
  /// @return -1 if the determinant is negative,
  /// @return  1 if the determinant is positive,
  /// @return  0 if the determinant is 0.
  static int signOfDet2x2(double dx1, double dy1, double dx2, double dy2) {
    DD x1 = DD.valueOf(dx1);
    DD y1 = DD.valueOf(dy1);
    DD x2 = DD.valueOf(dx2);
    DD y2 = DD.valueOf(dy2);

    DD det = x1.multiplyDD(y2).selfSubtractDD(y1.multiplyDD(x2));
    return det.signum();
  }

  /// A value which is safely greater than the
  /// relative round-off error in double-precision numbers
  static final double DP_SAFE_EPSILON = 1e-15;

  /// A filter for computing the orientation index of three coordinates.
  /// <p>
  /// If the orientation can be computed safely using standard DP
  /// arithmetic, this routine returns the orientation index.
  /// Otherwise, a value i > 1 is returned.
  /// In this case the orientation index must
  /// be computed using some other more robust method.
  /// The filter is fast to compute, so can be used to
  /// avoid the use of slower robust methods except when they are really needed,
  /// thus providing better average performance.
  /// <p>
  /// Uses an approach due to Jonathan Shewchuk, which is in the  domain.
  ///
  /// @param pa a coordinate
  /// @param pb a coordinate
  /// @param pc a coordinate
  /// @return the orientation index if it can be computed safely
  /// @return i > 1 if the orientation index cannot be computed safely
  static int orientationIndexFilter(Coordinate pa, Coordinate pb, Coordinate pc) {
    double detsum;

    double detleft = (pa.x - pc.x) * (pb.y - pc.y);
    double detright = (pa.y - pc.y) * (pb.x - pc.x);
    double det = detleft - detright;

    if (detleft > 0.0) {
      if (detright <= 0.0) {
        return signum(det);
      } else {
        detsum = detleft + detright;
      }
    } else if (detleft < 0.0) {
      if (detright >= 0.0) {
        return signum(det);
      } else {
        detsum = -detleft - detright;
      }
    } else {
      return signum(det);
    }

    double errbound = DP_SAFE_EPSILON * detsum;
    if ((det >= errbound) || (-det >= errbound)) {
      return signum(det);
    }

    return 2;
  }

  static int signum(double x) {
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
  }

  /// Computes an intersection point between two lines
  /// using DD arithmetic.
  /// If the lines are parallel (either identical
  /// or separate) a null value is returned.
  ///
  /// @param p1 an endpoint of line segment 1
  /// @param p2 an endpoint of line segment 1
  /// @param q1 an endpoint of line segment 2
  /// @param q2 an endpoint of line segment 2
  /// @return an intersection point if one exists, or null if the lines are parallel
  static Coordinate intersection(Coordinate p1, Coordinate p2, Coordinate q1, Coordinate q2) {
    DD px = new DD(p1.y).selfSubtract(p2.y);
    DD py = new DD(p2.x).selfSubtract(p1.x);
    DD pw = new DD(p1.x).selfMultiply(p2.y).selfSubtractDD(DD(p2.x).selfMultiply(p1.y));

    DD qx = new DD(q1.y).selfSubtract(q2.y);
    DD qy = new DD(q2.x).selfSubtract(q1.x);
    DD qw = new DD(q1.x).selfMultiply(q2.y).selfSubtractDD(new DD(q2.x).selfMultiply(q1.y));

    DD x = py.multiplyDD(qw).selfSubtractDD(qy.multiplyDD(pw));
    DD y = qx.multiplyDD(pw).selfSubtractDD(px.multiplyDD(qw));
    DD w = px.multiplyDD(qy).selfSubtractDD(qx.multiplyDD(py));

    double xInt = x.selfDivideDD(w).doubleValue();
    double yInt = y.selfDivideDD(w).doubleValue();

    if ((xInt.isNaN) || (xInt.isInfinite || yInt.isNaN) || (yInt.isInfinite)) {
      return null;
    }

    return new Coordinate(xInt, yInt);
  }
}

/// Functions to compute the orientation of basic geometric structures
/// including point triplets (triangles) and rings.
/// Orientation is a fundamental property of planar geometries
/// (and more generally geometry on two-dimensional manifolds).
/// <p>
/// Orientation is notoriously subject to numerical precision errors
/// in the case of collinear or nearly collinear points.
/// JTS uses extended-precision arithmetic to increase
/// the robustness of the computation.
///
/// @author Martin Davis
///
class Orientation {
  /// A value that indicates an orientation of clockwise, or a right turn.
  static final int CLOCKWISE = -1;

  /// A value that indicates an orientation of clockwise, or a right turn.
  static final int RIGHT = CLOCKWISE;

  /// A value that indicates an orientation of counterclockwise, or a left turn.
  static final int COUNTERCLOCKWISE = 1;

  /**
   * A value that indicates an orientation of counterclockwise, or a left turn.
   */
  static final int LEFT = COUNTERCLOCKWISE;

  /**
   * A value that indicates an orientation of collinear, or no turn (straight).
   */
  static final int COLLINEAR = 0;

  /**
   * A value that indicates an orientation of collinear, or no turn (straight).
   */
  static final int STRAIGHT = COLLINEAR;

  /**
   * Returns the orientation index of the direction of the point <code>q</code> relative to
   * a directed infinite line specified by <code>p1-p2</code>.
   * The index indicates whether the point lies to the {@link #LEFT} or {@link #RIGHT}
   * of the line, or lies on it {@link #COLLINEAR}.
   * The index also indicates the orientation of the triangle formed by the three points
   * ( {@link #COUNTERCLOCKWISE}, {@link #CLOCKWISE}, or {@link #STRAIGHT} )
   *
   * @param p1 the origin point of the line vector
   * @param p2 the final point of the line vector
   * @param q the point to compute the direction to
   *
   * @return -1 ( {@link #CLOCKWISE} or {@link #RIGHT} ) if q is clockwise (right) from p1-p2;
   *         1 ( {@link #COUNTERCLOCKWISE} or {@link #LEFT} ) if q is counter-clockwise (left) from p1-p2;
   *         0 ( {@link #COLLINEAR} or {@link #STRAIGHT} ) if q is collinear with p1-p2
   */
  static int index(Coordinate p1, Coordinate p2, Coordinate q) {
    /*
     * MD - 9 Aug 2010 It seems that the basic algorithm is slightly orientation
     * dependent, when computing the orientation of a point very close to a
     * line. This is possibly due to the arithmetic in the translation to the
     * origin.
     * 
     * For instance, the following situation produces identical results in spite
     * of the inverse orientation of the line segment:
     * 
     * Coordinate p0 = new Coordinate(219.3649559090992, 140.84159161824724);
     * Coordinate p1 = new Coordinate(168.9018919682399, -5.713787599646864);
     * 
     * Coordinate p = new Coordinate(186.80814046338352, 46.28973405831556); int
     * orient = orientationIndex(p0, p1, p); int orientInv =
     * orientationIndex(p1, p0, p);
     * 
     * A way to force consistent results is to normalize the orientation of the
     * vector using the following code. However, this may make the results of
     * orientationIndex inconsistent through the triangle of points, so it's not
     * clear this is an appropriate patch.
     * 
     */
    return CGAlgorithmsDD.orientationIndex(p1, p2, q);

    // testing only
    //return ShewchuksDeterminant.orientationIndex(p1, p2, q);
    // previous implementation - not quite fully robust
    //return RobustDeterminant.orientationIndex(p1, p2, q);
  }

  /**
   * Computes whether a ring defined by an array of {@link Coordinate}s is
   * oriented counter-clockwise.
   * <ul>
   * <li>The list of points is assumed to have the first and last points equal.
   * <li>This will handle coordinate lists which contain repeated points.
   * </ul>
   * This algorithm is <b>only</b> guaranteed to work with valid rings. If the
   * ring is invalid (e.g. self-crosses or touches), the computed result may not
   * be correct.
   *
   * @param ring
   *          an array of Coordinates forming a ring
   * @return true if the ring is oriented counter-clockwise.
   * @throws IllegalArgumentException
   *           if there are too few points to determine orientation (&lt; 4)
   */
  static bool isCCW(List<Coordinate> ring) {
    // # of points without closing endpoint
    int nPts = ring.length - 1;
    // sanity check
    if (nPts < 3) throw new ArgumentError("Ring has fewer than 4 points, so orientation cannot be determined");

    // find highest point
    Coordinate hiPt = ring[0];
    int hiIndex = 0;
    for (int i = 1; i <= nPts; i++) {
      Coordinate p = ring[i];
      if (p.y > hiPt.y) {
        hiPt = p;
        hiIndex = i;
      }
    }

    // find distinct point before highest point
    int iPrev = hiIndex;
    do {
      iPrev = iPrev - 1;
      if (iPrev < 0) iPrev = nPts;
    } while (ring[iPrev].equals2D(hiPt) && iPrev != hiIndex);

    // find distinct point after highest point
    int iNext = hiIndex;
    do {
      iNext = (iNext + 1) % nPts;
    } while (ring[iNext].equals2D(hiPt) && iNext != hiIndex);

    Coordinate prev = ring[iPrev];
    Coordinate next = ring[iNext];

    /*
     * This check catches cases where the ring contains an A-B-A configuration
     * of points. This can happen if the ring does not contain 3 distinct points
     * (including the case where the input array has fewer than 4 elements), or
     * it contains coincident line segments.
     */
    if (prev.equals2D(hiPt) || next.equals2D(hiPt) || prev.equals2D(next)) return false;

    int disc = Orientation.index(prev, hiPt, next);

    /*
     * If disc is exactly 0, lines are collinear. There are two possible cases:
     * (1) the lines lie along the x axis in opposite directions (2) the lines
     * lie on top of one another
     * 
     * (1) is handled by checking if next is left of prev ==> CCW (2) will never
     * happen if the ring is valid, so don't check for it (Might want to assert
     * this)
     */
    bool isCCW;
    if (disc == 0) {
      // poly is CCW if prev x is right of next x
      isCCW = (prev.x > next.x);
    } else {
      // if area is positive, points are ordered CCW
      isCCW = (disc > 0);
    }
    return isCCW;
  }

  /**
   * Computes whether a ring defined by an {@link CoordinateSequence} is
   * oriented counter-clockwise.
   * <ul>
   * <li>The list of points is assumed to have the first and last points equal.
   * <li>This will handle coordinate lists which contain repeated points.
   * </ul>
   * This algorithm is <b>only</b> guaranteed to work with valid rings. If the
   * ring is invalid (e.g. self-crosses or touches), the computed result may not
   * be correct.
   *
   * @param ring
   *          a CoordinateSequence forming a ring
   * @return true if the ring is oriented counter-clockwise.
   * @throws IllegalArgumentException
   *           if there are too few points to determine orientation (&lt; 4)
   */
  static bool isCCWFromSeq(CoordinateSequence ring) {
    // # of points without closing endpoint
    int nPts = ring.size() - 1;
    // sanity check
    if (nPts < 3) throw new ArgumentError("Ring has fewer than 4 points, so orientation cannot be determined");

    // find highest point
    Coordinate hiPt = ring.getCoordinate(0);
    int hiIndex = 0;
    for (int i = 1; i <= nPts; i++) {
      Coordinate p = ring.getCoordinate(i);
      if (p.y > hiPt.y) {
        hiPt = p;
        hiIndex = i;
      }
    }

    // find distinct point before highest point
    Coordinate prev;
    int iPrev = hiIndex;
    do {
      iPrev = iPrev - 1;
      if (iPrev < 0) iPrev = nPts;
      prev = ring.getCoordinate(iPrev);
    } while (prev.equals2D(hiPt) && iPrev != hiIndex);

    // find distinct point after highest point
    Coordinate next;
    int iNext = hiIndex;
    do {
      iNext = (iNext + 1) % nPts;
      next = ring.getCoordinate(iNext);
    } while (next.equals2D(hiPt) && iNext != hiIndex);

    /*
     * This check catches cases where the ring contains an A-B-A configuration
     * of points. This can happen if the ring does not contain 3 distinct points
     * (including the case where the input array has fewer than 4 elements), or
     * it contains coincident line segments.
     */
    if (prev.equals2D(hiPt) || next.equals2D(hiPt) || prev.equals2D(next)) return false;

    int disc = Orientation.index(prev, hiPt, next);

    /*
     * If disc is exactly 0, lines are collinear. There are two possible cases:
     * (1) the lines lie along the x axis in opposite directions (2) the lines
     * lie on top of one another
     *
     * (1) is handled by checking if next is left of prev ==> CCW (2) will never
     * happen if the ring is valid, so don't check for it (Might want to assert
     * this)
     */
    bool isCCW;
    if (disc == 0) {
      // poly is CCW if prev x is right of next x
      isCCW = (prev.x > next.x);
    } else {
      // if area is positive, points are ordered CCW
      isCCW = (disc > 0);
    }
    return isCCW;
  }
}

/// Contains functions to compute intersections between lines.
///
/// @author Martin Davis
///
class Intersection {
  /// Computes the intersection point of two lines.
  /// If the lines are parallel or collinear this case is detected
  /// and <code>null</code> is returned.
  /// <p>
  /// In general it is not possible to accurately compute
  /// the intersection point of two lines, due to
  /// numerical roundoff.
  /// This is particularly true when the input lines are nearly parallel.
  /// This routine uses numerical conditioning on the input values
  /// to ensure that the computed value should be very close to the correct value.
  ///
  /// @param p1 an endpoint of line 1
  /// @param p2 an endpoint of line 1
  /// @param q1 an endpoint of line 2
  /// @param q2 an endpoint of line 2
  /// @return the intersection point between the lines, if there is one,
  /// or null if the lines are parallel or collinear
  ///
  /// @see CGAlgorithmsDD#intersection(Coordinate, Coordinate, Coordinate, Coordinate)
  static Coordinate intersection(Coordinate p1, Coordinate p2, Coordinate q1, Coordinate q2) {
    // compute midpoint of "kernel envelope"
    double minX0 = p1.x < p2.x ? p1.x : p2.x;
    double minY0 = p1.y < p2.y ? p1.y : p2.y;
    double maxX0 = p1.x > p2.x ? p1.x : p2.x;
    double maxY0 = p1.y > p2.y ? p1.y : p2.y;

    double minX1 = q1.x < q2.x ? q1.x : q2.x;
    double minY1 = q1.y < q2.y ? q1.y : q2.y;
    double maxX1 = q1.x > q2.x ? q1.x : q2.x;
    double maxY1 = q1.y > q2.y ? q1.y : q2.y;

    double intMinX = minX0 > minX1 ? minX0 : minX1;
    double intMaxX = maxX0 < maxX1 ? maxX0 : maxX1;
    double intMinY = minY0 > minY1 ? minY0 : minY1;
    double intMaxY = maxY0 < maxY1 ? maxY0 : maxY1;

    double midx = (intMinX + intMaxX) / 2.0;
    double midy = (intMinY + intMaxY) / 2.0;

    // condition ordinate values by subtracting midpoint
    double p1x = p1.x - midx;
    double p1y = p1.y - midy;
    double p2x = p2.x - midx;
    double p2y = p2.y - midy;
    double q1x = q1.x - midx;
    double q1y = q1.y - midy;
    double q2x = q2.x - midx;
    double q2y = q2.y - midy;

    // unrolled computation using homogeneous coordinates eqn
    double px = p1y - p2y;
    double py = p2x - p1x;
    double pw = p1x * p2y - p2x * p1y;

    double qx = q1y - q2y;
    double qy = q2x - q1x;
    double qw = q1x * q2y - q2x * q1y;

    double x = py * qw - qy * pw;
    double y = qx * pw - px * qw;
    double w = px * qy - qx * py;

    double xInt = x / w;
    double yInt = y / w;

    // check for parallel lines
    if ((xInt.isNaN) || (xInt.isInfinite || yInt.isNaN) || (yInt.isInfinite)) {
      return null;
    }
    // de-condition intersection point
    return new Coordinate(xInt + midx, yInt + midy);
  }
}

/// A <code>LineIntersector</code> is an algorithm that can both test whether
/// two line segments intersect and compute the intersection point(s)
/// if they do.
/// <p>
/// There are three possible outcomes when determining whether two line segments intersect:
/// <ul>
/// <li>{@link #NO_INTERSECTION} - the segments do not intersect
/// <li>{@link #POINT_INTERSECTION} - the segments intersect in a single point
/// <li>{@link #COLLINEAR_INTERSECTION} - the segments are collinear and they intersect in a line segment
/// </ul>
/// For segments which intersect in a single point, the point may be either an endpoint
/// or in the interior of each segment.
/// If the point lies in the interior of both segments,
/// this is termed a <i>proper intersection</i>.
/// The method {@link #isProper()} test for this situation.
/// <p>
/// The intersection point(s) may be computed in a precise or non-precise manner.
/// Computing an intersection point precisely involves rounding it
/// via a supplied {@link PrecisionModel}.
/// <p>
/// LineIntersectors do not perform an initial envelope intersection test
/// to determine if the segments are disjoint.
/// This is because this class is likely to be used in a context where
/// envelope overlap is already known to occur (or be likely).
///
/// @version 1.7
abstract class LineIntersector {
  /// These are deprecated, due to ambiguous naming
  static final int DONT_INTERSECT = 0;
  static final int DO_INTERSECT = 1;
  static final int COLLINEAR = 2;

  /// Indicates that line segments do not intersect
  static final int NO_INTERSECTION = 0;

  /// Indicates that line segments intersect in a single point
  static final int POINT_INTERSECTION = 1;

  /// Indicates that line segments intersect in a line segment
  static final int COLLINEAR_INTERSECTION = 2;

  /// Computes the "edge distance" of an intersection point p along a segment.
  /// The edge distance is a metric of the point along the edge.
  /// The metric used is a robust and easy to compute metric function.
  /// It is <b>not</b> equivalent to the usual Euclidean metric.
  /// It relies on the fact that either the x or the y ordinates of the
  /// points in the edge are unique, depending on whether the edge is longer in
  /// the horizontal or vertical direction.
  /// <p>
  /// NOTE: This function may produce incorrect distances
  ///  for inputs where p is not precisely on p1-p2
  /// (E.g. p = (139,9) p1 = (139,10), p2 = (280,1) produces distance 0.0, which is incorrect.
  /// <p>
  /// My hypothesis is that the function is safe to use for points which are the
  /// result of <b>rounding</b> points which lie on the line,
  /// but not safe to use for <b>truncated</b> points.
  static double computeEdgeDistance(Coordinate p, Coordinate p0, Coordinate p1) {
    double dx = (p1.x - p0.x).abs();
    double dy = (p1.y - p0.y).abs();

    double dist = -1.0; // sentinel value
    if (p.equals(p0)) {
      dist = 0.0;
    } else if (p.equals(p1)) {
      if (dx > dy)
        dist = dx;
      else
        dist = dy;
    } else {
      double pdx = (p.x - p0.x).abs();
      double pdy = (p.y - p0.y).abs();
      if (dx > dy)
        dist = pdx;
      else
        dist = pdy;
      // <FIX>
      // hack to ensure that non-endpoints always have a non-zero distance
      if (dist == 0.0 && !p.equals(p0)) {
        dist = math.max(pdx, pdy);
      }
    }
    //TODO Assert.isTrue(! (dist == 0.0 && ! p.equals(p0)), "Bad distance calculation");
    return dist;
  }

  /**
   * This function is non-robust, since it may compute the square of large numbers.
   * Currently not sure how to improve this.
   */
  static double nonRobustComputeEdgeDistance(Coordinate p, Coordinate p1, Coordinate p2) {
    double dx = p.x - p1.x;
    double dy = p.y - p1.y;
    double dist = math.sqrt(dx * dx + dy * dy); // dummy value
    //TODO Assert.isTrue(! (dist == 0.0 && ! p.equals(p1)), "Invalid distance calculation");
    return dist;
  }

  int result;
  List<List<Coordinate>> inputLines = MatrixUtils.createMatrix(2, 2, null as Coordinate);
  List<Coordinate> intPt = List(2);

  /// The indexes of the endpoints of the intersection lines, in order along
  /// the corresponding line
  List<List<int>> intLineIndex;
  bool _isProper;
  Coordinate pa;
  Coordinate pb;

  /// If makePrecise is true, computed intersection coordinates will be made precise
  /// using Coordinate#makePrecise
  PrecisionModel precisionModel = null;

// int numIntersects = 0;

  LineIntersector() {
    intPt[0] = Coordinate.empty2D();
    intPt[1] = Coordinate.empty2D();
    // alias the intersection points for ease of reference
    pa = intPt[0];
    pb = intPt[1];
    result = 0;
  }

  /// Force computed intersection to be rounded to a given precision model
  /// @param precisionModel
  /// @deprecated use <code>setPrecisionModel</code> instead
  void setMakePrecise(PrecisionModel precisionModel) {
    this.precisionModel = precisionModel;
  }

  /// Force computed intersection to be rounded to a given precision model.
  /// No getter is provided, because the precision model is not required to be specified.
  /// @param precisionModel
  void setPrecisionModel(PrecisionModel precisionModel) {
    this.precisionModel = precisionModel;
  }

  /// Gets an endpoint of an input segment.
  ///
  /// @param segmentIndex the index of the input segment (0 or 1)
  /// @param ptIndex the index of the endpoint (0 or 1)
  /// @return the specified endpoint
  Coordinate getEndpoint(int segmentIndex, int ptIndex) {
    return inputLines[segmentIndex][ptIndex];
  }

  /// Compute the intersection of a point p and the line p1-p2.
  /// This function computes the bool value of the hasIntersection test.
  /// The actual value of the intersection (if there is one)
  /// is equal to the value of <code>p</code>.
  void computeIntersectionPointLine(Coordinate p, Coordinate p1, Coordinate p2);

  bool isCollinear() {
    return result == COLLINEAR_INTERSECTION;
  }

  /// Computes the intersection of the lines p1-p2 and p3-p4.
  /// This function computes both the bool value of the hasIntersection test
  /// and the (approximate) value of the intersection point itself (if there is one).
  void computeIntersection(Coordinate p1, Coordinate p2, Coordinate p3, Coordinate p4) {
    inputLines[0][0] = p1;
    inputLines[0][1] = p2;
    inputLines[1][0] = p3;
    inputLines[1][1] = p4;
    result = computeIntersect(p1, p2, p3, p4);
//numIntersects++;
  }

  int computeIntersect(Coordinate p1, Coordinate p2, Coordinate q1, Coordinate q2);

/*
   String toString() {
    String str = inputLines[0][0] + "-"
         + inputLines[0][1] + " "
         + inputLines[1][0] + "-"
         + inputLines[1][1] + " : "
               + getTopologySummary();
    return str;
  }
*/

  String toString() {
    return WKTWriter.toLineStringFromCoords(inputLines[0][0], inputLines[0][1]) +
        " - " +
        WKTWriter.toLineStringFromCoords(inputLines[1][0], inputLines[1][1]) +
        getTopologySummary();
  }

  String getTopologySummary() {
    StringBuffer catBuilder = new StringBuffer();
    if (isEndPoint()) catBuilder.write(" endpoint");
    if (_isProper) catBuilder.write(" proper");
    if (isCollinear()) catBuilder.write(" collinear");
    return catBuilder.toString();
  }

  bool isEndPoint() {
    return hasIntersection() && !_isProper;
  }

  /// Tests whether the input geometries intersect.
  ///
  /// @return true if the input geometries intersect
  bool hasIntersection() {
    return result != NO_INTERSECTION;
  }

  /// Returns the number of intersection points found.  This will be either 0, 1 or 2.
  ///
  /// @return the number of intersection points found (0, 1, or 2)
  int getIntersectionNum() {
    return result;
  }

  /// Returns the intIndex'th intersection point
  ///
  /// @param intIndex is 0 or 1
  ///
  /// @return the intIndex'th intersection point
  Coordinate getIntersection(int intIndex) {
    return intPt[intIndex];
  }

  void computeIntLineIndex() {
    if (intLineIndex == null) {
      intLineIndex = MatrixUtils.createMatrix(2, 2, 0 as int);
      computeIntLineIndexWithIndex(0);
      computeIntLineIndexWithIndex(1);
    }
  }

  void computeIntLineIndexWithIndex(int segmentIndex) {
    double dist0 = getEdgeDistance(segmentIndex, 0);
    double dist1 = getEdgeDistance(segmentIndex, 1);
    if (dist0 > dist1) {
      intLineIndex[segmentIndex][0] = 0;
      intLineIndex[segmentIndex][1] = 1;
    } else {
      intLineIndex[segmentIndex][0] = 1;
      intLineIndex[segmentIndex][1] = 0;
    }
  }

  /// Test whether a point is a intersection point of two line segments.
  /// Note that if the intersection is a line segment, this method only tests for
  /// equality with the endpoints of the intersection segment.
  /// It does <b>not</b> return true if
  /// the input point is internal to the intersection segment.
  ///
  /// @return true if the input point is one of the intersection points.
  bool isIntersection(Coordinate pt) {
    for (int i = 0; i < result; i++) {
      if (intPt[i].equals2D(pt)) {
        return true;
      }
    }
    return false;
  }

  /// Tests whether either intersection point is an interior point of one of the input segments.
  ///
  /// @return <code>true</code> if either intersection point is in the interior of one of the input segments
  bool isInteriorIntersection() {
    if (isInteriorIntersectionWithIndex(0)) return true;
    if (isInteriorIntersectionWithIndex(1)) return true;
    return false;
  }

  /// Tests whether either intersection point is an interior point of the specified input segment.
  ///
  /// @return <code>true</code> if either intersection point is in the interior of the input segment
  bool isInteriorIntersectionWithIndex(int inputLineIndex) {
    for (int i = 0; i < result; i++) {
      if (!(intPt[i].equals2D(inputLines[inputLineIndex][0]) || intPt[i].equals2D(inputLines[inputLineIndex][1]))) {
        return true;
      }
    }
    return false;
  }

  /// Tests whether an intersection is proper.
  /// <br>
  /// The intersection between two line segments is considered proper if
  /// they intersect in a single point in the interior of both segments
  /// (e.g. the intersection is a single point and is not equal to any of the
  /// endpoints).
  /// <p>
  /// The intersection between a point and a line segment is considered proper
  /// if the point lies in the interior of the segment (e.g. is not equal to
  /// either of the endpoints).
  ///
  /// @return true if the intersection is proper
  bool isProper() {
    return hasIntersection() && _isProper;
  }

  /// Computes the intIndex'th intersection point in the direction of
  /// a specified input line segment
  ///
  /// @param segmentIndex is 0 or 1
  /// @param intIndex is 0 or 1
  ///
  /// @return the intIndex'th intersection point in the direction of the specified input line segment
  Coordinate getIntersectionAlongSegment(int segmentIndex, int intIndex) {
    // lazily compute int line array
    computeIntLineIndex();
    return intPt[intLineIndex[segmentIndex][intIndex]];
  }

  /// Computes the index (order) of the intIndex'th intersection point in the direction of
  /// a specified input line segment
  ///
  /// @param segmentIndex is 0 or 1
  /// @param intIndex is 0 or 1
  ///
  /// @return the index of the intersection point along the input segment (0 or 1)
  int getIndexAlongSegment(int segmentIndex, int intIndex) {
    computeIntLineIndex();
    return intLineIndex[segmentIndex][intIndex];
  }

  /// Computes the "edge distance" of an intersection point along the specified input line segment.
  ///
  /// @param segmentIndex is 0 or 1
  /// @param intIndex is 0 or 1
  ///
  /// @return the edge distance of the intersection point
  double getEdgeDistance(int segmentIndex, int intIndex) {
    double dist = computeEdgeDistance(intPt[intIndex], inputLines[segmentIndex][0], inputLines[segmentIndex][1]);
    return dist;
  }
}

/**
 * A robust version of {@link LineIntersector}.
 *
 * @version 1.7
 */
class RobustLineIntersector extends LineIntersector {
  RobustLineIntersector() {}

  void computeIntersectionPointLine(Coordinate p, Coordinate p1, Coordinate p2) {
    _isProper = false;
    // do between check first, since it is faster than the orientation test
    if (Envelope.intersectsPoint(p1, p2, p)) {
      if ((Orientation.index(p1, p2, p) == 0) && (Orientation.index(p2, p1, p) == 0)) {
        _isProper = true;
        if (p.equals(p1) || p.equals(p2)) {
          _isProper = false;
        }
        result = LineIntersector.POINT_INTERSECTION;
        return;
      }
    }
    result = LineIntersector.NO_INTERSECTION;
  }

  int computeIntersect(Coordinate p1, Coordinate p2, Coordinate q1, Coordinate q2) {
    _isProper = false;

    // first try a fast test to see if the envelopes of the lines intersect
    if (!Envelope.intersectsEnvelopeCoords(p1, p2, q1, q2)) return LineIntersector.NO_INTERSECTION;

    // for each endpoint, compute which side of the other segment it lies
    // if both endpoints lie on the same side of the other segment,
    // the segments do not intersect
    int Pq1 = Orientation.index(p1, p2, q1);
    int Pq2 = Orientation.index(p1, p2, q2);

    if ((Pq1 > 0 && Pq2 > 0) || (Pq1 < 0 && Pq2 < 0)) {
      return LineIntersector.NO_INTERSECTION;
    }

    int Qp1 = Orientation.index(q1, q2, p1);
    int Qp2 = Orientation.index(q1, q2, p2);

    if ((Qp1 > 0 && Qp2 > 0) || (Qp1 < 0 && Qp2 < 0)) {
      return LineIntersector.NO_INTERSECTION;
    }

    bool collinear = Pq1 == 0 && Pq2 == 0 && Qp1 == 0 && Qp2 == 0;
    if (collinear) {
      return computeCollinearIntersection(p1, p2, q1, q2);
    }

    /**
     * At this point we know that there is a single intersection point
     * (since the lines are not collinear).
     */

    /**
     *  Check if the intersection is an endpoint. If it is, copy the endpoint as
     *  the intersection point. Copying the point rather than computing it
     *  ensures the point has the exact value, which is important for
     *  robustness. It is sufficient to simply check for an endpoint which is on
     *  the other line, since at this point we know that the inputLines must
     *  intersect.
     */
    if (Pq1 == 0 || Pq2 == 0 || Qp1 == 0 || Qp2 == 0) {
      _isProper = false;

      /**
       * Check for two equal endpoints.
       * This is done explicitly rather than by the orientation tests
       * below in order to improve robustness.
       *
       * [An example where the orientation tests fail to be consistent is
       * the following (where the true intersection is at the shared endpoint
       * POINT (19.850257749638203 46.29709338043669)
       *
       * LINESTRING ( 19.850257749638203 46.29709338043669, 20.31970698357233 46.76654261437082 )
       * and
       * LINESTRING ( -48.51001596420236 -22.063180333403878, 19.850257749638203 46.29709338043669 )
       *
       * which used to produce the INCORRECT result: (20.31970698357233, 46.76654261437082, NaN)
       *
       */
      if (p1.equals2D(q1) || p1.equals2D(q2)) {
        intPt[0] = p1;
      } else if (p2.equals2D(q1) || p2.equals2D(q2)) {
        intPt[0] = p2;
      }

      /**
       * Now check to see if any endpoint lies on the interior of the other segment.
       */
      else if (Pq1 == 0) {
        intPt[0] = new Coordinate.fromCoordinate(q1);
      } else if (Pq2 == 0) {
        intPt[0] = new Coordinate.fromCoordinate(q2);
      } else if (Qp1 == 0) {
        intPt[0] = new Coordinate.fromCoordinate(p1);
      } else if (Qp2 == 0) {
        intPt[0] = new Coordinate.fromCoordinate(p2);
      }
    } else {
      _isProper = true;
      intPt[0] = intersection(p1, p2, q1, q2);
    }
    return LineIntersector.POINT_INTERSECTION;
  }

  int computeCollinearIntersection(Coordinate p1, Coordinate p2, Coordinate q1, Coordinate q2) {
    bool p1q1p2 = Envelope.intersectsPoint(p1, p2, q1);
    bool p1q2p2 = Envelope.intersectsPoint(p1, p2, q2);
    bool q1p1q2 = Envelope.intersectsPoint(q1, q2, p1);
    bool q1p2q2 = Envelope.intersectsPoint(q1, q2, p2);

    if (p1q1p2 && p1q2p2) {
      intPt[0] = q1;
      intPt[1] = q2;
      return LineIntersector.COLLINEAR_INTERSECTION;
    }
    if (q1p1q2 && q1p2q2) {
      intPt[0] = p1;
      intPt[1] = p2;
      return LineIntersector.COLLINEAR_INTERSECTION;
    }
    if (p1q1p2 && q1p1q2) {
      intPt[0] = q1;
      intPt[1] = p1;
      return q1.equals(p1) && !p1q2p2 && !q1p2q2 ? LineIntersector.POINT_INTERSECTION : LineIntersector.COLLINEAR_INTERSECTION;
    }
    if (p1q1p2 && q1p2q2) {
      intPt[0] = q1;
      intPt[1] = p2;
      return q1.equals(p2) && !p1q2p2 && !q1p1q2 ? LineIntersector.POINT_INTERSECTION : LineIntersector.COLLINEAR_INTERSECTION;
    }
    if (p1q2p2 && q1p1q2) {
      intPt[0] = q2;
      intPt[1] = p1;
      return q2.equals(p1) && !p1q1p2 && !q1p2q2 ? LineIntersector.POINT_INTERSECTION : LineIntersector.COLLINEAR_INTERSECTION;
    }
    if (p1q2p2 && q1p2q2) {
      intPt[0] = q2;
      intPt[1] = p2;
      return q2.equals(p2) && !p1q1p2 && !q1p1q2 ? LineIntersector.POINT_INTERSECTION : LineIntersector.COLLINEAR_INTERSECTION;
    }
    return LineIntersector.NO_INTERSECTION;
  }

  /**
   * This method computes the actual value of the intersection point.
   * To obtain the maximum precision from the intersection calculation,
   * the coordinates are normalized by subtracting the minimum
   * ordinate values (in absolute value).  This has the effect of
   * removing common significant digits from the calculation to
   * maintain more bits of precision.
   */
  Coordinate intersection(Coordinate p1, Coordinate p2, Coordinate q1, Coordinate q2) {
    Coordinate intPt = intersectionSafe(p1, p2, q1, q2);

    /*
    // TESTING ONLY
    Coordinate intPtDD = CGAlgorithmsDD.intersection(p1, p2, q1, q2);
    double dist = intPt.distance(intPtDD);
    System.out.println(intPt + " - " + intPtDD + " dist = " + dist);
    //intPt = intPtDD;
    */

    /**
     * Due to rounding it can happen that the computed intersection is
     * outside the envelopes of the input segments.  Clearly this
     * is inconsistent.
     * This code checks this condition and forces a more reasonable answer
     *
     * MD - May 4 2005 - This is still a problem.  Here is a failure case:
     *
     * LINESTRING (2089426.5233462777 1180182.3877339689, 2085646.6891757075 1195618.7333999649)
     * LINESTRING (1889281.8148903656 1997547.0560044837, 2259977.3672235999 483675.17050843034)
     * int point = (2097408.2633752143,1144595.8008114607)
     *
     * MD - Dec 14 2006 - This does not seem to be a failure case any longer
     */
    if (!isInSegmentEnvelopes(intPt)) {
//      System.out.println("Intersection outside segment envelopes: " + intPt);

      // compute a safer result
      // copy the coordinate, since it may be rounded later
      intPt = new Coordinate.fromCoordinate(nearestEndpoint(p1, p2, q1, q2));
//    intPt = CentralEndpointIntersector.getIntersection(p1, p2, q1, q2);

//      System.out.println("Segments: " + this);
//      System.out.println("Snapped to " + intPt);
//      checkDD(p1, p2, q1, q2, intPt);
    }
    if (precisionModel != null) {
      precisionModel.makeCoordinatePrecise(intPt);
    }
    return intPt;
  }

  void checkDD(Coordinate p1, Coordinate p2, Coordinate q1, Coordinate q2, Coordinate intPt) {
    Coordinate intPtDD = CGAlgorithmsDD.intersection(p1, p2, q1, q2);
    bool isIn = isInSegmentEnvelopes(intPtDD);
    print("DD in env = $isIn --------------------- $intPtDD");
    if (intPt.distance(intPtDD) > 0.0001) {
      print("Distance = ${intPt.distance(intPtDD)}");
    }
  }

  /**
   * Computes a segment intersection using homogeneous coordinates.
   * Round-off error can cause the raw computation to fail,
   * (usually due to the segments being approximately parallel).
   * If this happens, a reasonable approximation is computed instead.
   *
   * @param p1 a segment endpoint
   * @param p2 a segment endpoint
   * @param q1 a segment endpoint
   * @param q2 a segment endpoint
   * @return the computed intersection point
   */
  Coordinate intersectionSafe(Coordinate p1, Coordinate p2, Coordinate q1, Coordinate q2) {
    Coordinate intPt = Intersection.intersection(p1, p2, q1, q2);
    if (intPt == null) intPt = nearestEndpoint(p1, p2, q1, q2);
    //     System.out.println("Snapped to " + intPt);
    return intPt;
  }

  /**
   * Tests whether a point lies in the envelopes of both input segments.
   * A correctly computed intersection point should return <code>true</code>
   * for this test.
   * Since this test is for debugging purposes only, no attempt is
   * made to optimize the envelope test.
   *
   * @return <code>true</code> if the input point lies within both input segment envelopes
   */
  bool isInSegmentEnvelopes(Coordinate intPt) {
    Envelope env0 = new Envelope.fromCoordinates(inputLines[0][0], inputLines[0][1]);
    Envelope env1 = new Envelope.fromCoordinates(inputLines[1][0], inputLines[1][1]);
    return env0.containsCoordinate(intPt) && env1.containsCoordinate(intPt);
  }

  /**
   * Finds the endpoint of the segments P and Q which
   * is closest to the other segment.
   * This is a reasonable surrogate for the true
   * intersection points in ill-conditioned cases
   * (e.g. where two segments are nearly coincident,
   * or where the endpoint of one segment lies almost on the other segment).
   * <p>
   * This replaces the older CentralEndpoint heuristic,
   * which chose the wrong endpoint in some cases
   * where the segments had very distinct slopes
   * and one endpoint lay almost on the other segment.
   *
   * @param p1 an endpoint of segment P
   * @param p2 an endpoint of segment P
   * @param q1 an endpoint of segment Q
   * @param q2 an endpoint of segment Q
   * @return the nearest endpoint to the other segment
   */
  static Coordinate nearestEndpoint(Coordinate p1, Coordinate p2, Coordinate q1, Coordinate q2) {
    Coordinate nearestPt = p1;
    double minDist = Distance.pointToSegment(p1, q1, q2);

    double dist = Distance.pointToSegment(p2, q1, q2);
    if (dist < minDist) {
      minDist = dist;
      nearestPt = p2;
    }
    dist = Distance.pointToSegment(q1, p1, p2);
    if (dist < minDist) {
      minDist = dist;
      nearestPt = q1;
    }
    dist = Distance.pointToSegment(q2, p1, p2);
    if (dist < minDist) {
      minDist = dist;
      nearestPt = q2;
    }
    return nearestPt;
  }
}

/// Counts the number of segments crossed by a horizontal ray extending to the right
/// from a given point, in an incremental fashion.
/// This can be used to determine whether a point lies in a {@link Polygonal} geometry.
/// The class determines the situation where the point lies exactly on a segment.
/// When being used for Point-In-Polygon determination, this case allows short-circuiting
/// the evaluation.
/// <p>
/// This class handles polygonal geometries with any number of shells and holes.
/// The orientation of the shell and hole rings is unimportant.
/// In order to compute a correct location for a given polygonal geometry,
/// it is essential that <b>all</b> segments are counted which
/// <ul>
/// <li>touch the ray
/// <li>lie in in any ring which may contain the point
/// </ul>
/// The only exception is when the point-on-segment situation is detected, in which
/// case no further processing is required.
/// The implication of the above rule is that segments
/// which can be a priori determined to <i>not</i> touch the ray
/// (i.e. by a test of their bounding box or Y-extent)
/// do not need to be counted.  This allows for optimization by indexing.
/// <p>
/// This implementation uses the extended-precision orientation test,
/// to provide maximum robustness and consistency within
/// other algorithms.
///
/// @author Martin Davis
///
class RayCrossingCounter {
  /// Determines the {@link Location} of a point in a ring.
  /// This method is an exemplar of how to use this class.
  ///
  /// @param p the point to test
  /// @param ring an array of Coordinates forming a ring
  /// @return the location of the point in the ring
  static int locatePointInRingList(Coordinate p, List<Coordinate> ring) {
    RayCrossingCounter counter = new RayCrossingCounter(p);

    for (int i = 1; i < ring.length; i++) {
      Coordinate p1 = ring[i];
      Coordinate p2 = ring[i - 1];
      counter.countSegment(p1, p2);
      if (counter.isOnSegment()) return counter.getLocation();
    }
    return counter.getLocation();
  }

  /**
   * Determines the {@link Location} of a point in a ring.
   *
   * @param p
   *            the point to test
   * @param ring
   *            a coordinate sequence forming a ring
   * @return the location of the point in the ring
   */
  static int locatePointInRing(Coordinate p, CoordinateSequence ring) {
    RayCrossingCounter counter = new RayCrossingCounter(p);

    Coordinate p1 = new Coordinate.empty2D();
    Coordinate p2 = new Coordinate.empty2D();
    for (int i = 1; i < ring.size(); i++) {
      ring.getCoordinateInto(i, p1);
      ring.getCoordinateInto(i - 1, p2);
      counter.countSegment(p1, p2);
      if (counter.isOnSegment()) return counter.getLocation();
    }
    return counter.getLocation();
  }

  Coordinate p;
  int crossingCount = 0;

  // true if the test point lies on an input segment
  bool isPointOnSegment = false;

  RayCrossingCounter(Coordinate p) {
    this.p = p;
  }

  /**
   * Counts a segment
   *
   * @param p1 an endpoint of the segment
   * @param p2 another endpoint of the segment
   */
  void countSegment(Coordinate p1, Coordinate p2) {
    /**
     * For each segment, check if it crosses
     * a horizontal ray running from the test point in the positive x direction.
     */

    // check if the segment is strictly to the left of the test point
    if (p1.x < p.x && p2.x < p.x) return;

    // check if the point is equal to the current ring vertex
    if (p.x == p2.x && p.y == p2.y) {
      isPointOnSegment = true;
      return;
    }
    /**
     * For horizontal segments, check if the point is on the segment.
     * Otherwise, horizontal segments are not counted.
     */
    if (p1.y == p.y && p2.y == p.y) {
      double minx = p1.x;
      double maxx = p2.x;
      if (minx > maxx) {
        minx = p2.x;
        maxx = p1.x;
      }
      if (p.x >= minx && p.x <= maxx) {
        isPointOnSegment = true;
      }
      return;
    }
    /**
     * Evaluate all non-horizontal segments which cross a horizontal ray to the
     * right of the test pt. To avoid double-counting shared vertices, we use the
     * convention that
     * <ul>
     * <li>an upward edge includes its starting endpoint, and excludes its
     * final endpoint
     * <li>a downward edge excludes its starting endpoint, and includes its
     * final endpoint
     * </ul>
     */
    if (((p1.y > p.y) && (p2.y <= p.y)) || ((p2.y > p.y) && (p1.y <= p.y))) {
      int orient = Orientation.index(p1, p2, p);
      if (orient == Orientation.COLLINEAR) {
        isPointOnSegment = true;
        return;
      }
      // Re-orient the result if needed to ensure effective segment direction is upwards
      if (p2.y < p1.y) {
        orient = -orient;
      }
      // The upward segment crosses the ray if the test point lies to the left (CCW) of the segment.
      if (orient == Orientation.LEFT) {
        crossingCount++;
      }
    }
  }

  /**
   * Reports whether the point lies exactly on one of the supplied segments.
   * This method may be called at any time as segments are processed.
   * If the result of this method is <tt>true</tt>,
   * no further segments need be supplied, since the result
   * will never change again.
   *
   * @return true if the point lies exactly on a segment
   */
  bool isOnSegment() {
    return isPointOnSegment;
  }

  /**
   * Gets the {@link Location} of the point relative to
   * the ring, polygon
   * or multipolygon from which the processed segments were provided.
   * <p>
   * This method only determines the correct location
   * if <b>all</b> relevant segments must have been processed.
   *
   * @return the Location of the point
   */
  int getLocation() {
    if (isPointOnSegment) return Location.BOUNDARY;

    // The point is in the interior of the ring if the number of X-crossings is
    // odd.
    if ((crossingCount % 2) == 1) {
      return Location.INTERIOR;
    }
    return Location.EXTERIOR;
  }

  /**
   * Tests whether the point lies in or on
   * the ring, polygon
   * or multipolygon from which the processed segments were provided.
   * <p>
   * This method only determines the correct location
   * if <b>all</b> relevant segments must have been processed.
   *
   * @return true if the point lies in or on the supplied polygon
   */
  bool isPointInPolygon() {
    return getLocation() != Location.EXTERIOR;
  }
}

/// Functions for locating points within basic geometric
/// structures such as lines and rings.
///
/// @author Martin Davis
///
class PointLocation {
  /// Tests whether a point lies on the line defined by a list of
  /// coordinates.
  ///
  /// @param p the point to test
  /// @param line the line coordinates
  /// @return true if the point is a vertex of the line or lies in the interior
  ///         of a line segment in the line
  static bool isOnLineList(Coordinate p, List<Coordinate> line) {
    LineIntersector lineIntersector = new RobustLineIntersector();
    for (int i = 1; i < line.length; i++) {
      Coordinate p0 = line[i - 1];
      Coordinate p1 = line[i];
      lineIntersector.computeIntersectionPointLine(p, p0, p1);
      if (lineIntersector.hasIntersection()) {
        return true;
      }
    }
    return false;
  }

  /// Tests whether a point lies on the line defined by a
  /// {@link CoordinateSequence}.
  ///
  /// @param p the point to test
  /// @param line the line coordinates
  /// @return true if the point is a vertex of the line or lies in the interior
  ///         of a line segment in the line
  static bool isOnLine(Coordinate p, CoordinateSequence line) {
    LineIntersector lineIntersector = new RobustLineIntersector();
    Coordinate p0 = new Coordinate.empty2D();
    Coordinate p1 = new Coordinate.empty2D();
    int n = line.size();
    for (int i = 1; i < n; i++) {
      line.getCoordinateInto(i - 1, p0);
      line.getCoordinateInto(i, p1);
      lineIntersector.computeIntersectionPointLine(p, p0, p1);
      if (lineIntersector.hasIntersection()) {
        return true;
      }
    }
    return false;
  }

  /// Tests whether a point lies inside or on a ring. The ring may be oriented in
  /// either direction. A point lying exactly on the ring boundary is considered
  /// to be inside the ring.
  /// <p>
  /// This method does <i>not</i> first check the point against the envelope of
  /// the ring.
  ///
  /// @param p
  ///          point to check for ring inclusion
  /// @param ring
  ///          an array of coordinates representing the ring (which must have
  ///          first point identical to last point)
  /// @return true if p is inside ring
  ///
  /// @see locatePointInRing
  static bool isInRing(Coordinate p, List<Coordinate> ring) {
    return PointLocation.locateInRing(p, ring) != Location.EXTERIOR;
  }

  /// Determines whether a point lies in the interior, on the boundary, or in the
  /// exterior of a ring. The ring may be oriented in either direction.
  /// <p>
  /// This method does <i>not</i> first check the point against the envelope of
  /// the ring.
  ///
  /// @param p
  ///          point to check for ring inclusion
  /// @param ring
  ///          an array of coordinates representing the ring (which must have
  ///          first point identical to last point)
  /// @return the {@link Location} of p relative to the ring
  static int locateInRing(Coordinate p, List<Coordinate> ring) {
    return RayCrossingCounter.locatePointInRingList(p, ring);
  }
}

/// Computes the topological ({@link Location})
/// of a single point to a {@link Geometry}.
/// A {@link BoundaryNodeRule} may be specified
/// to control the evaluation of whether the point lies on the boundary or not
/// The default rule is to use the the <i>SFS Boundary Determination Rule</i>
/// <p>
/// Notes:
/// <ul>
/// <li>{@link LinearRing}s do not enclose any area - points inside the ring are still in the EXTERIOR of the ring.
/// </ul>
/// Instances of this class are not reentrant.
///
/// @version 1.7
class PointLocator {
  // default is to use OGC SFS rule
  BoundaryNodeRule boundaryRule =
      //BoundaryNodeRule.ENDPOINT_BOUNDARY_RULE;
      BoundaryNodeRule.OGC_SFS_BOUNDARY_RULE;

  bool isIn; // true if the point lies in or on any Geometry element
  int numBoundaries; // the number of sub-elements whose boundaries the point lies in

  PointLocator() {}

  PointLocator.fromRule(BoundaryNodeRule boundaryRule) {
    if (boundaryRule == null) throw new ArgumentError("Rule must be non-null");
    this.boundaryRule = boundaryRule;
  }

  /// Convenience method to test a point for intersection with
  /// a Geometry
  /// @param p the coordinate to test
  /// @param geom the Geometry to test
  /// @return <code>true</code> if the point is in the interior or boundary of the Geometry
  bool intersects(Coordinate p, Geometry geom) {
    return locate(p, geom) != Location.EXTERIOR;
  }

  /**
   * Computes the topological relationship ({@link Location}) of a single point
   * to a Geometry.
   * It handles both single-element
   * and multi-element Geometries.
   * The algorithm for multi-part Geometries
   * takes into account the SFS Boundary Determination Rule.
   *
   * @return the {@link Location} of the point relative to the input Geometry
   */
  int locate(Coordinate p, Geometry geom) {
    if (geom.isEmpty()) return Location.EXTERIOR;

    if (geom is LineString) {
      return locateOnLineString(p, geom);
    } else if (geom is Polygon) {
      return locateInPolygon(p, geom);
    }

    isIn = false;
    numBoundaries = 0;
    computeLocation(p, geom);
    if (boundaryRule.isInBoundary(numBoundaries)) return Location.BOUNDARY;
    if (numBoundaries > 0 || isIn) return Location.INTERIOR;

    return Location.EXTERIOR;
  }

  void computeLocation(Coordinate p, Geometry geom) {
    if (geom is Point) {
      updateLocationInfo(locateOnPoint(p, geom));
    }
    if (geom is LineString) {
      updateLocationInfo(locateOnLineString(p, geom));
    } else if (geom is Polygon) {
      updateLocationInfo(locateInPolygon(p, geom));
    } else if (geom is MultiLineString) {
      MultiLineString ml = geom;
      for (int i = 0; i < ml.getNumGeometries(); i++) {
        LineString l = ml.getGeometryN(i);
        updateLocationInfo(locateOnLineString(p, l));
      }
    } else if (geom is MultiPolygon) {
      MultiPolygon mpoly = geom;
      for (int i = 0; i < mpoly.getNumGeometries(); i++) {
        Polygon poly = mpoly.getGeometryN(i);
        updateLocationInfo(locateInPolygon(p, poly));
      }
    } else if (geom is GeometryCollection) {
      GeometryCollectionIterator geomi = new GeometryCollectionIterator(geom);
      while (geomi.hasNext()) {
        Geometry g2 = geomi.next();
        if (g2 != geom) computeLocation(p, g2);
      }
    }
  }

  void updateLocationInfo(int loc) {
    if (loc == Location.INTERIOR) isIn = true;
    if (loc == Location.BOUNDARY) numBoundaries++;
  }

  int locateOnPoint(Coordinate p, Point pt) {
    // no point in doing envelope test, since equality test is just as fast

    Coordinate ptCoord = pt.getCoordinate();
    if (ptCoord.equals2D(p)) return Location.INTERIOR;
    return Location.EXTERIOR;
  }

  int locateOnLineString(Coordinate p, LineString l) {
    // bounding-box check
    if (!l.getEnvelopeInternal().intersectsCoordinate(p)) return Location.EXTERIOR;

    CoordinateSequence seq = l.getCoordinateSequence();
    if (!l.isClosed()) {
      if (p.equals(seq.getCoordinate(0)) || p.equals(seq.getCoordinate(seq.size() - 1))) {
        return Location.BOUNDARY;
      }
    }
    if (PointLocation.isOnLine(p, seq)) {
      return Location.INTERIOR;
    }
    return Location.EXTERIOR;
  }

  int locateInPolygonRing(Coordinate p, LinearRing ring) {
    // bounding-box check
    if (!ring.getEnvelopeInternal().intersectsCoordinate(p)) return Location.EXTERIOR;

    return PointLocation.locateInRing(p, ring.getCoordinates());
  }

  int locateInPolygon(Coordinate p, Polygon poly) {
    if (poly.isEmpty()) return Location.EXTERIOR;

    LinearRing shell = poly.getExteriorRing();

    int shellLoc = locateInPolygonRing(p, shell);
    if (shellLoc == Location.EXTERIOR) return Location.EXTERIOR;
    if (shellLoc == Location.BOUNDARY) return Location.BOUNDARY;
    // now test if the point lies in or on the holes
    for (int i = 0; i < poly.getNumInteriorRing(); i++) {
      LinearRing hole = poly.getInteriorRingN(i);
      int holeLoc = locateInPolygonRing(p, hole);
      if (holeLoc == Location.INTERIOR) return Location.EXTERIOR;
      if (holeLoc == Location.BOUNDARY) return Location.BOUNDARY;
    }
    return Location.INTERIOR;
  }
}

/**
 * Computes whether a rectangle intersects line segments.
 * <p>
 * Rectangles contain a large amount of inherent symmetry
 * (or to put it another way, although they contain four
 * coordinates they only actually contain 4 ordinates
 * worth of information).
 * The algorithm used takes advantage of the symmetry of
 * the geometric situation
 * to optimize performance by minimizing the number
 * of line intersection tests.
 *
 * @author Martin Davis
 *
 */
class RectangleLineIntersector {
  // for intersection testing, don't need to set precision model
  LineIntersector li = new RobustLineIntersector();

  Envelope rectEnv;

  Coordinate diagUp0;
  Coordinate diagUp1;
  Coordinate diagDown0;
  Coordinate diagDown1;

  /**
   * Creates a new intersector for the given query rectangle,
   * specified as an {@link Envelope}.
   *
   *
   * @param rectEnv the query rectangle, specified as an Envelope
   */
  RectangleLineIntersector(Envelope rectEnv) {
    this.rectEnv = rectEnv;

    /**
     * Up and Down are the diagonal orientations
     * relative to the Left side of the rectangle.
     * Index 0 is the left side, 1 is the right side.
     */
    diagUp0 = new Coordinate(rectEnv.getMinX(), rectEnv.getMinY());
    diagUp1 = new Coordinate(rectEnv.getMaxX(), rectEnv.getMaxY());
    diagDown0 = new Coordinate(rectEnv.getMinX(), rectEnv.getMaxY());
    diagDown1 = new Coordinate(rectEnv.getMaxX(), rectEnv.getMinY());
  }

  /**
   * Tests whether the query rectangle intersects a
   * given line segment.
   *
   * @param p0 the first endpoint of the segment
   * @param p1 the second endpoint of the segment
   * @return true if the rectangle intersects the segment
   */
  bool intersects(Coordinate p0, Coordinate p1) {
    // TODO: confirm that checking envelopes first is faster

    /**
     * If the segment envelope is disjoint from the
     * rectangle envelope, there is no intersection
     */
    Envelope segEnv = new Envelope.fromCoordinates(p0, p1);
    if (!rectEnv.intersectsEnvelope(segEnv)) return false;

    /**
     * If either segment endpoint lies in the rectangle,
     * there is an intersection.
     */
    if (rectEnv.intersectsCoordinate(p0)) return true;
    if (rectEnv.intersectsCoordinate(p1)) return true;

    /**
     * Normalize segment.
     * This makes p0 less than p1,
     * so that the segment runs to the right,
     * or vertically upwards.
     */
    if (p0.compareTo(p1) > 0) {
      Coordinate tmp = p0;
      p0 = p1;
      p1 = tmp;
    }
    /**
     * Compute angle of segment.
     * Since the segment is normalized to run left to right,
     * it is sufficient to simply test the Y ordinate.
     * "Upwards" means relative to the left end of the segment.
     */
    bool isSegUpwards = false;
    if (p1.y > p0.y) isSegUpwards = true;

    /**
     * Since we now know that neither segment endpoint
     * lies in the rectangle, there are two possible
     * situations:
     * 1) the segment is disjoint to the rectangle
     * 2) the segment crosses the rectangle completely.
     *
     * In the case of a crossing, the segment must intersect
     * a diagonal of the rectangle.
     *
     * To distinguish these two cases, it is sufficient
     * to test intersection with
     * a single diagonal of the rectangle,
     * namely the one with slope "opposite" to the slope
     * of the segment.
     * (Note that if the segment is axis-parallel,
     * it must intersect both diagonals, so this is
     * still sufficient.)
     */
    if (isSegUpwards) {
      li.computeIntersection(p0, p1, diagDown0, diagDown1);
    } else {
      li.computeIntersection(p0, p1, diagUp0, diagUp1);
    }
    if (li.hasIntersection()) return true;
    return false;
  }
}

/**
 * Functions for computing length.
 *
 * @author Martin Davis
 *
 */
class Length {
  /**
   * Computes the length of a linestring specified by a sequence of points.
   *
   * @param pts the points specifying the linestring
   * @return the length of the linestring
   */
  static double ofLine(CoordinateSequence pts) {
    // optimized for processing CoordinateSequences
    int n = pts.size();
    if (n <= 1) return 0.0;

    double len = 0.0;

    Coordinate p = new Coordinate.empty2D();
    pts.getCoordinateInto(0, p);
    double x0 = p.x;
    double y0 = p.y;

    for (int i = 1; i < n; i++) {
      pts.getCoordinateInto(i, p);
      double x1 = p.x;
      double y1 = p.y;
      double dx = x1 - x0;
      double dy = y1 - y0;

      len += math.sqrt(dx * dx + dy * dy);

      x0 = x1;
      y0 = y1;
    }
    return len;
  }
}

/**
 * Functions for computing area.
 *
 * @author Martin Davis
 *
 */
class Area {
  /**
   * Computes the area for a ring.
   *
   * @param ring the coordinates forming the ring
   * @return the area of the ring
   */
  static double ofRing(List<Coordinate> ring) {
    return ofRingSigned(ring).abs();
  }

  /**
   * Computes the area for a ring.
   *
   * @param ring the coordinates forming the ring
   * @return the area of the ring
   */
  static double ofRingSeq(CoordinateSequence ring) {
    return ofRingSignedSeq(ring).abs();
  }

  /**
   * Computes the signed area for a ring. The signed area is positive if the
   * ring is oriented CW, negative if the ring is oriented CCW, and zero if the
   * ring is degenerate or flat.
   *
   * @param ring
   *          the coordinates forming the ring
   * @return the signed area of the ring
   */
  static double ofRingSigned(List<Coordinate> ring) {
    if (ring.length < 3) return 0.0;
    double sum = 0.0;
    /**
     * Based on the Shoelace formula.
     * http://en.wikipedia.org/wiki/Shoelace_formula
     */
    double x0 = ring[0].x;
    for (int i = 1; i < ring.length - 1; i++) {
      double x = ring[i].x - x0;
      double y1 = ring[i + 1].y;
      double y2 = ring[i - 1].y;
      sum += x * (y2 - y1);
    }
    return sum / 2.0;
  }

  /**
   * Computes the signed area for a ring. The signed area is:
   * <ul>
   * <li>positive if the ring is oriented CW
   * <li>negative if the ring is oriented CCW
   * <li>zero if the ring is degenerate or flat
   * </ul>
   *
   * @param ring
   *          the coordinates forming the ring
   * @return the signed area of the ring
   */
  static double ofRingSignedSeq(CoordinateSequence ring) {
    int n = ring.size();
    if (n < 3) return 0.0;
    /**
     * Based on the Shoelace formula.
     * http://en.wikipedia.org/wiki/Shoelace_formula
     */
    Coordinate p0 = new Coordinate.empty2D();
    Coordinate p1 = new Coordinate.empty2D();
    Coordinate p2 = new Coordinate.empty2D();
    ring.getCoordinateInto(0, p1);
    ring.getCoordinateInto(1, p2);
    double x0 = p1.x;
    p2.x -= x0;
    double sum = 0.0;
    for (int i = 1; i < n - 1; i++) {
      p0.y = p1.y;
      p1.x = p2.x;
      p1.y = p2.y;
      ring.getCoordinateInto(i + 1, p2);
      p2.x -= x0;
      sum += p1.x * (p0.y - p2.y);
    }
    return sum / 2.0;
  }
}
