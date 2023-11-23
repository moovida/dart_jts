part of dart_jts;

class TaggedLineString {
  ///
  /// Represents a {@link LineString} which can be modified to a simplified shape.
  /// This class provides an attribute which specifies the minimum allowable length
  /// for the modified result.
  ///
  /// @version 1.7
  ///
  LineString? parentLine;
  List<TaggedLineSegment> segs = List.empty(growable: true);
  List<LineSegment> resultSegs = List.empty(growable: true);
  late int minimumSize;
  TaggedLineString(this.parentLine, this.minimumSize) {
    init();
  }
  TaggedLineString.fromLineString(LineString parent) {
    parentLine = parent;
    minimumSize = 2;
  }
  int getMinimumSize() {
    return minimumSize;
  }

  LineString? getParent() {
    return parentLine;
  }

  List<Coordinate> getParentCoordinates() {
    return parentLine!.getCoordinates();
  }

  List<Coordinate> getResultCoordinates() {
    return extractCoordinates(resultSegs);
  }

  int getResultSize() {
    int resultSegsSize = resultSegs.length;
    return resultSegsSize == 0 ? 0 : resultSegsSize + 1;
  }

  TaggedLineSegment getSegment(int i) {
    return segs[i];
  }

  void init() {
    List<Coordinate> pts = parentLine!.getCoordinates();
    segs = List.empty(growable: true);
    for (int i = 0; i < pts.length - 1; i++) {
      TaggedLineSegment seg =
          new TaggedLineSegment(pts[i], pts[i + 1], parentLine, i);
      segs.add(seg);
    }
  }

  List<TaggedLineSegment> getSegments() {
    return segs;
  }

  void addToResult(LineSegment seg) {
    resultSegs.add(seg);
  }

  LineString asLineString() {
    return parentLine!
        .getFactory()
        .createLineString(extractCoordinates(resultSegs));
  }

  LinearRing asLinearRing() {
    return parentLine!
        .getFactory()
        .createLinearRing(extractCoordinates(resultSegs));
  }

  List<Coordinate> extractCoordinates(List segs) {
    List<Coordinate> pts = List.empty(growable: true);
    LineSegment? seg = null;
    for (int i = 0; i < segs.length; i++) {
      seg = segs[i] as LineSegment;
      pts.add(seg.p0);
    }
    // add last point
    pts.add(seg!.p1);
    return pts;
  }
}
