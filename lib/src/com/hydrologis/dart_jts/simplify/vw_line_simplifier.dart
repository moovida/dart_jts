part of dart_jts;

///
/// Simplifies a linestring (sequence of points) using the
/// Visvalingam-Whyatt algorithm.
/// The Visvalingam-Whyatt algorithm simplifies geometry
/// by removing vertices while trying to minimize the area changed.
///
/// @version 1.7
///
class VWLineSimplifier {
  List<Coordinate> pts;
  double tolerance;

  VWLineSimplifier(this.pts, this.tolerance);

  List<Coordinate> simplify() {
    VWVertex vwLine = VWVertex.buildLine(pts);
    double minArea = tolerance;
    do {
      minArea = simplifyVertex(vwLine);
    } while (minArea < tolerance);
    List<Coordinate> simp = vwLine.getCoordinates();
    if (simp.length < 2) {
      List<Coordinate> list = [simp[0], Coordinate.fromCoordinate(simp[0])];
      return list;
    }
    return simp;
  }

  double simplifyVertex(VWVertex vwLine) {
    ///
    /// Scan vertices in line and remove the one with smallest effective area.
    ///
    // TODO: use an appropriate data structure to optimize finding the smallest area vertex
    VWVertex? curr = vwLine;
    double minArea = curr.getArea();
    VWVertex? minVertex;
    while (curr != null) {
      double area = curr.getArea();
      if (area < minArea) {
        minArea = area;
        minVertex = curr;
      }
      curr = curr.next;
    }
    if (minVertex != null && minArea < tolerance) {
      minVertex.remove();
    }
    if (!vwLine.isLive) return -1;
    return minArea;
  }
}

class VWVertex {
  static buildLine(List<Coordinate> pts) {
    VWVertex? first;
    VWVertex? prev;
    for (int i = 0; i < pts.length; i++) {
      var v = VWVertex(pts[i]);
      if (first == null) {
        first = v;
      }
      v.setPrev(prev!);

      prev.setNext(v);
      prev.updateArea();
      prev = v;
    }
    return first;
  }

  Coordinate pt;
  VWVertex? prev;
  VWVertex? next;
  double area = double.maxFinite;
  bool isLive = true;

  VWVertex(this.pt);

  void setPrev(VWVertex prev) {
    this.prev = prev;
  }

  void setNext(VWVertex next) {
    this.next = next;
  }

  void updateArea() {
    if (prev == null || next == null) {
      area = double.maxFinite;
      return;
    }
    area = (Triangle.areaStatic(prev!.pt, pt, next!.pt));
  }

  double getArea() {
    return area;
  }

  VWVertex? remove() {
    VWVertex? tmpPrev = prev;
    VWVertex? tmpNext = next;
    VWVertex? result;
    if (prev != null) {
      prev!.setNext(tmpNext!);
      prev!.updateArea();
      result = prev;
    }
    if (next != null) {
      next!.setPrev(tmpPrev!);
      next!.updateArea();
      if (result == null) result = next;
    }
    isLive = false;
    return result;
  }

  List<Coordinate> getCoordinates() {
    CoordinateList coords = CoordinateList();
    VWVertex? curr = this;
    do {
      coords.add(curr!.pt);
      curr = curr.next;
    } while (curr != null);
    return coords.toCoordinateArray();
  }
}
