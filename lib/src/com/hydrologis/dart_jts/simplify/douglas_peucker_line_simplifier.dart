part of dart_jts;

class DouglasPeuckerLineSimplifier {
  List<Coordinate> pts;
  List<bool> usePt = List.empty(growable:true);
  double _distanceTolerance = 1;
  LineSegment? seg;

  DouglasPeuckerLineSimplifier(this.pts);
  ///
  /// Simplifies a linestring (sequence of points) using
  /// the standard Douglas-Peucker algorithm.
  ///
  /// @version 1.7
  ///
  ///
  /// Sets the distance tolerance for the simplification.
  /// All vertices in the simplified linestring will be within this
  /// distance of the original linestring.
  ///
  /// @param distanceTolerance the approximation tolerance to use
  ///
  ///
  void setDistanceTolerance(double distanceTolerance) {
    _distanceTolerance = distanceTolerance;
  }
  List<Coordinate> simplify()
  {
    usePt = List.filled(pts.length, true);
    _simplifySection(0, pts.length - 1);
    CoordinateList coordList = new CoordinateList();
    for (int i = 0; i < pts.length; i++) {
      if (usePt[i])
        coordList.add(Coordinate(pts[i].x,pts[i].y));
    }
    return coordList.toCoordinateArray(true);
  }


  void _simplifySection(int i, int j)
  {
    if((i+1) == j) {
      return;
    }
    seg = LineSegment(pts[i].x,pts[i].y,pts[j].x,pts[j].y);
    double maxDistance = -1.0;
    int maxIndex = i;
    for (int k = i + 1; k < j; k++) {
      double distance = seg!.distanceCoord(pts[k]);
      if (distance > maxDistance) {
        maxDistance = distance;
        maxIndex = k;
      }
    }
    if (maxDistance <= _distanceTolerance) {
      for(int k = i + 1; k < j; k++) {
        usePt[k] = false;
      }
    }
    else {
      _simplifySection(i, maxIndex);
      _simplifySection(maxIndex, j);
    }
  }

}