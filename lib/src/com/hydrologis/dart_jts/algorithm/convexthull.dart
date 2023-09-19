part of dart_jts;

/*
 * Copyright (c) 2016 Vivid Solutions.
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * and Eclipse Distribution License v. 1.0 which accompanies this distribution.
 * The Eclipse Public License is available at http://www.eclipse.org/legal/epl-v10.html
 * and the Eclipse Distribution License is available at
 *
 * http://www.eclipse.org/org/documents/edl-v10.php.
 */

/**
 * Computes the convex hull of a {@link Geometry}.
 * The convex hull is the smallest convex Geometry that contains all the
 * points in the input Geometry.
 * <p>
 * Uses the Graham Scan algorithm.
 *
 *@version 1.7
 */
class ConvexHull {
  late GeometryFactory geomFactory;
  late List<Coordinate> inputPts;

  /**
   * Create a new convex hull construction for the input {@link Geometry}.
   */
  ConvexHull(Geometry geometry)
      : this.fromPoints(extractCoordinates(geometry), geometry.getFactory());

  /**
   * Create a new convex hull construction for the input {@link Coordinate} array.
   */
  ConvexHull.fromPoints(List<Coordinate> pts, GeometryFactory geomFactory) {
    inputPts = UniqueCoordinateArrayFilter.filterCoordinates(pts);
    //inputPts = pts;
    this.geomFactory = geomFactory;
  }

  static List<Coordinate> extractCoordinates(Geometry geom) {
    UniqueCoordinateArrayFilter filter = new UniqueCoordinateArrayFilter();
    geom.applyCF(filter);
    return filter.getCoordinates();
  }

  /**
   * Returns a {@link Geometry} that represents the convex hull of the input
   * geometry.
   * The returned geometry contains the minimal number of points needed to
   * represent the convex hull.  In particular, no more than two consecutive
   * points will be collinear.
   *
   * @return if the convex hull contains 3 or more points, a {@link Polygon};
   * 2 points, a {@link LineString};
   * 1 point, a {@link Point};
   * 0 points, an empty {@link GeometryCollection}.
   */
  Geometry getConvexHull() {
    if (inputPts.length == 0) {
      return geomFactory.createGeometryCollectionEmpty();
    }
    if (inputPts.length == 1) {
      return geomFactory.createPoint(inputPts[0]);
    }
    if (inputPts.length == 2) {
      return geomFactory.createLineString(inputPts);
    }

    List<Coordinate> reducedPts = inputPts;
    // use heuristic to reduce points, if large
    if (inputPts.length > 50) {
      reducedPts = reduce(inputPts);
    }
    // sort points for Graham scan.
    List<Coordinate> sortedPts = preSort(reducedPts);

    // Use Graham scan to find convex hull.
    List<Coordinate> cHS = grahamScan(sortedPts);

    // Convert stack to an array.
    List<Coordinate> cH = toCoordinateArray(cHS);

    // Convert array to appropriate output geometry.
    return lineOrPolygon(cH);
  }

  /**
   * An alternative to Stack.toArray, which is not present in earlier versions
   * of Java.
   */
  List<Coordinate> toCoordinateArray(List<Coordinate> stack) {
    List<Coordinate> coordinates = [];
    for (int i = 0; i < stack.length; i++) {
      Coordinate coordinate = stack[i];
      coordinates.add(coordinate);
    }
    return coordinates;
  }

  /**
   * Uses a heuristic to reduce the number of points scanned
   * to compute the hull.
   * The heuristic is to find a polygon guaranteed to
   * be in (or on) the hull, and eliminate all points inside it.
   * A quadrilateral defined by the extremal points
   * in the four orthogonal directions
   * can be used, but even more inclusive is
   * to use an octilateral defined by the points in the 8 cardinal directions.
   * <p>
   * Note that even if the method used to determine the polygon vertices
   * is not 100% robust, this does not affect the robustness of the convex hull.
   * <p>
   * To satisfy the requirements of the Graham Scan algorithm, 
   * the returned array has at least 3 entries.
   *
   * @param pts the points to reduce
   * @return the reduced list of points (at least 3)
   */
  List<Coordinate> reduce(List<Coordinate> inputPts) {
    //List<Coordinate> polyPts = computeQuad(inputPts);
    List<Coordinate>? polyPts = computeOctRing(inputPts);
    //List<Coordinate> polyPts = null;

    // unable to compute interior polygon for some reason
    if (polyPts == null) return inputPts;

//    LinearRing ring = geomFactory.createLinearRing(polyPts);
//    System.out.println(ring);

    // add points defining polygon
    var reducedSet = SplayTreeSet<Coordinate>();
    // TreeSet reducedSet = new TreeSet();
    for (int i = 0; i < polyPts.length; i++) {
      reducedSet.add(polyPts[i]);
    }
    /**
     * Add all unique points not in the interior poly.
     * CGAlgorithms.isPointInRing is not defined for points actually on the ring,
     * but this doesn't matter since the points of the interior polygon
     * are forced to be in the reduced set.
     */
    for (int i = 0; i < inputPts.length; i++) {
      if (!PointLocation.isInRing(inputPts[i], polyPts)) {
        reducedSet.add(inputPts[i]);
      }
    }
    List<Coordinate> reducedPts = List.from(
        reducedSet); // CoordinateArrays.toCoordinateArray(reducedSet);

    // ensure that computed array has at least 3 points (not necessarily unique)
    if (reducedPts.length < 3) return padArray3(reducedPts);
    return reducedPts;
  }

  List<Coordinate> padArray3(List<Coordinate> pts) {
    List<Coordinate> pad =
        List.filled(3, Coordinate.empty2D()); // new Coordinate[3];
    for (int i = 0; i < pad.length; i++) {
      if (i < pts.length) {
        pad[i] = pts[i];
      } else
        pad[i] = pts[0];
    }
    return pad;
  }

  List<Coordinate> preSort(List<Coordinate> pts) {
    Coordinate t;

    // find the lowest point in the set. If two or more points have
    // the same minimum y coordinate choose the one with the minimu x.
    // This focal point is put in array location pts[0].
    for (int i = 1; i < pts.length; i++) {
      var pt0 = pts[0];
      if ((pts[i].y < pt0.y) || ((pts[i].y == pt0.y) && (pts[i].x < pt0.x))) {
        t = pt0;
        pts[0] = pts[i];
        pts[i] = t;
      }
    }

    // sort the points radially around the focal point.
    var radialComparator = RadialComparator(pts[0]);
    pts.sublist(1).sort((o1, o2) {
      return radialComparator.compare(o1, o2);
    });

    // Arrays.sort(pts, 1, pts.length, new RadialComparator(pts[0]));

    //radialSort(pts);
    return pts;
  }

  /**
   * Uses the Graham Scan algorithm to compute the convex hull vertices.
   * 
   * @param c a list of points, with at least 3 entries
   * @return a Stack containing the ordered points of the convex hull ring
   */
  List<Coordinate> grahamScan(List<Coordinate> c) {
    Coordinate p;
    List<Coordinate> ps = [];
    ps.add(c[0]);
    ps.add(c[1]);
    ps.add(c[2]);
    for (int i = 3; i < c.length; i++) {
      p = ps.removeLast(); // pop();
      // check for empty stack to guard against robustness problems
      while (!ps.isEmpty && Orientation.index(ps.last, p, c[i]) > 0) {
        p = ps.removeLast();
      }
      ps.add(p);
      ps.add(c[i]);
    }
    ps.add(c[0]);
    return ps;
  }

  /**
   *@return    whether the three coordinates are collinear and c2 lies between
   *      c1 and c3 inclusive
   */
  bool isBetween(Coordinate c1, Coordinate c2, Coordinate c3) {
    if (Orientation.index(c1, c2, c3) != 0) {
      return false;
    }
    if (c1.x != c3.x) {
      if (c1.x <= c2.x && c2.x <= c3.x) {
        return true;
      }
      if (c3.x <= c2.x && c2.x <= c1.x) {
        return true;
      }
    }
    if (c1.y != c3.y) {
      if (c1.y <= c2.y && c2.y <= c3.y) {
        return true;
      }
      if (c3.y <= c2.y && c2.y <= c1.y) {
        return true;
      }
    }
    return false;
  }

  List<Coordinate>? computeOctRing(List<Coordinate> inputPts) {
    List<Coordinate> octPts = computeOctPts(inputPts);
    CoordinateList coordList = new CoordinateList();
    coordList.addList(octPts, false);

    // points must all lie in a line
    if (coordList.size() < 3) {
      return null;
    }
    coordList.closeRing();
    return coordList.toCoordinateArray();
  }

  List<Coordinate> computeOctPts(List<Coordinate> inputPts) {
    List<Coordinate> pts = List.filled(8, Coordinate.empty2D());
    for (int j = 0; j < pts.length; j++) {
      pts[j] = inputPts[0];
    }
    for (int i = 1; i < inputPts.length; i++) {
      if (inputPts[i].x < pts[0].x) {
        pts[0] = inputPts[i];
      }
      if (inputPts[i].x - inputPts[i].y < pts[1].x - pts[1].y) {
        pts[1] = inputPts[i];
      }
      if (inputPts[i].y > pts[2].y) {
        pts[2] = inputPts[i];
      }
      if (inputPts[i].x + inputPts[i].y > pts[3].x + pts[3].y) {
        pts[3] = inputPts[i];
      }
      if (inputPts[i].x > pts[4].x) {
        pts[4] = inputPts[i];
      }
      if (inputPts[i].x - inputPts[i].y > pts[5].x - pts[5].y) {
        pts[5] = inputPts[i];
      }
      if (inputPts[i].y < pts[6].y) {
        pts[6] = inputPts[i];
      }
      if (inputPts[i].x + inputPts[i].y < pts[7].x + pts[7].y) {
        pts[7] = inputPts[i];
      }
    }
    return pts;
  }

  /**
   *@param  vertices  the vertices of a linear ring, which may or may not be
   *      flattened (i.e. vertices collinear)
   *@return           a 2-vertex <code>LineString</code> if the vertices are
   *      collinear; otherwise, a <code>Polygon</code> with unnecessary
   *      (collinear) vertices removed
   */
  Geometry lineOrPolygon(List<Coordinate> coordinates) {
    coordinates = cleanRing(coordinates);
    if (coordinates.length == 3) {
      return geomFactory
          .createLineString(<Coordinate>[coordinates[0], coordinates[1]]);
//      return new LineString(new List<Coordinate>{coordinates[0], coordinates[1]},
//          geometry.getPrecisionModel(), geometry.getSRID());
    }
    LinearRing linearRing = geomFactory.createLinearRing(coordinates);
    return geomFactory.createPolygon(linearRing, null);
  }

  /**
   *@param  vertices  the vertices of a linear ring, which may or may not be
   *      flattened (i.e. vertices collinear)
   *@return           the coordinates with unnecessary (collinear) vertices
   *      removed
   */
  List<Coordinate> cleanRing(List<Coordinate> original) {
    Assert.equals(original[0], original[original.length - 1]);
    List<Coordinate> cleanedRing = [];
    Coordinate? previousDistinctCoordinate = null;
    for (int i = 0; i <= original.length - 2; i++) {
      Coordinate currentCoordinate = original[i];
      Coordinate nextCoordinate = original[i + 1];
      if (currentCoordinate.equals(nextCoordinate)) {
        continue;
      }
      if (previousDistinctCoordinate != null &&
          isBetween(
              previousDistinctCoordinate, currentCoordinate, nextCoordinate)) {
        continue;
      }
      cleanedRing.add(currentCoordinate);
      previousDistinctCoordinate = currentCoordinate;
    }
    cleanedRing.add(original[original.length - 1]);
    return cleanedRing;
    // List<Coordinate> cleanedRingCoordinates = new Coordinate[cleanedRing.size()];
    // return (List<Coordinate>) cleanedRing.toArray(cleanedRingCoordinates);
  }
}

//   Comparator radialomparator = (o1, o2) {
//     IntervalRTreeNode n1 = o1 as IntervalRTreeNode;
//     IntervalRTreeNode n2 = o2 as IntervalRTreeNode;
//     double mid1 = (n1.min + n1.max) / 2;
//     double mid2 = (n2.min + n2.max) / 2;
//     if (mid1 < mid2) return -1;
//     if (mid1 > mid2) return 1;
//     return 0;
// };

/**
   * Compares {@link Coordinate}s for their angle and distance
   * relative to an origin.
   *
   * @author Martin Davis
   * @version 1.7
   */
class RadialComparator {
  Coordinate origin;

  RadialComparator(this.origin);

  int compare(Object o1, Object o2) {
    Coordinate p1 = o1 as Coordinate;
    Coordinate p2 = o2 as Coordinate;
    return polarCompare(origin, p1, p2);
  }

  /**
     * Given two points p and q compare them with respect to their radial
     * ordering about point o.  First checks radial ordering.
     * If points are collinear, the comparison is based
     * on their distance to the origin.
     * <p>
     * p < q iff
     * <ul>
     * <li>ang(o-p) < ang(o-q) (e.g. o-p-q is CCW)
     * <li>or ang(o-p) == ang(o-q) && dist(o,p) < dist(o,q)
     * </ul>
     *
     * @param o the origin
     * @param p a point
     * @param q another point
     * @return -1, 0 or 1 depending on whether p is less than,
     * equal to or greater than q
     */
  int polarCompare(Coordinate o, Coordinate p, Coordinate q) {
    double dxp = p.x - o.x;
    double dyp = p.y - o.y;
    double dxq = q.x - o.x;
    double dyq = q.y - o.y;

    int orient = Orientation.index(o, p, q);

    if (orient == Orientation.COUNTERCLOCKWISE) return 1;
    if (orient == Orientation.CLOCKWISE) return -1;

    // points are collinear - check distance
    double op = dxp * dxp + dyp * dyp;
    double oq = dxq * dxq + dyq * dyq;
    if (op < oq) {
      return -1;
    }
    if (op > oq) {
      return 1;
    }
    return 0;
  }
}
