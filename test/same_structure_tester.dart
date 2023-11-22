import 'package:dart_jts/dart_jts.dart';

class SameStructureTester {

  static bool isSameStructure(Geometry g1, Geometry g2)
  {
    if (!g1.isEquivalentClass(g2))
      return false;
    if (g1 is GeometryCollection)
      return isSameStructureCollection(g1, g2 as GeometryCollection);
    else if (g1 is Polygon)
      return isSameStructurePolygon(g1, g2 as Polygon);
    else if (g1 is LineString)
      return isSameStructureLineString(g1, g2 as LineString);
    else if (g1 is Point)
      return isSameStructurePoint(g1, g2 as Point);

    Assert.shouldNeverReachHere(
        "Unsupported Geometry class: " + g1.getGeometryType());
    return false;
  }

  static bool isSameStructureCollection(GeometryCollection g1, GeometryCollection g2)
  {
    if (g1.getNumGeometries() != g2.getNumGeometries())
      return false;
    for (int i = 0; i < g1.getNumGeometries(); i++) {
      if (! isSameStructure(g1.getGeometryN(i), g2.getGeometryN(i)))
        return false;
    }
    return true;
  }

  static bool isSameStructurePolygon(Polygon g1, Polygon g2)
  {
    if (g1.getNumInteriorRing() != g2.getNumInteriorRing())
      return false;
    // could check for both empty or nonempty here
    return true;
  }

  static bool isSameStructureLineString(LineString g1, LineString g2)
  {
    // could check for both empty or nonempty here
    return true;
  }

  static bool isSameStructurePoint(Point g1, Point g2)
  {
    // could check for both empty or nonempty here
    return true;
  }

}
