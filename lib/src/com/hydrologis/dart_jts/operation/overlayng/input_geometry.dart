part of dart_jts;

/**
 * Manages the input geometries for an overlay operation.
 * The second geometry is allowed to be null,
 * to support for instance precision reduction.
 *
 * @author Martin Davis
 *
 */
class InputGeometry {
  List<Geometry?> geom = List.empty(growable: true);
  PointOnGeometryLocator? ptLocatorA;
  PointOnGeometryLocator? ptLocatorB;
  List<bool> isCollapsed = List.filled(2, false);

  InputGeometry(Geometry geomA, Geometry? geomB) {
    geom.add(geomA);
    geom.add(geomB);
  }

  bool isSingle() {
    return geom.length < 2;
  }

  int getDimension(int index) {
    if (geom[index] == null) return -1;
    return geom[index]!.getDimension();
  }

  Geometry? getGeometry(int geomIndex) {
    return geom[geomIndex];
  }

  Envelope? getEnvelope(int geomIndex) {
    if (geom[geomIndex] == null) {
      return null;
    }
    return geom[geomIndex]!.getEnvelopeInternal();
  }

  bool isEmpty(int geomIndex) {
    return geom[geomIndex] != null && geom[geomIndex]!.isEmpty();
  }

  bool isArea(int geomIndex) {
    return geom[geomIndex] != null && geom[geomIndex]!.getDimension() == 2;
  }

  /**
   * Gets the index of an input which is an area,
   * if one exists.
   * Otherwise returns -1.
   * If both inputs are areas, returns the index of the first one (0).
   *
   * @return the index of an area input, or -1
   */
  int getAreaIndex() {
    if (getDimension(0) == 2) return 0;
    if (getDimension(1) == 2) return 1;
    return -1;
  }

  bool isLine(int geomIndex) {
    return getDimension(geomIndex) == 1;
  }

  bool isAllPoints() {
    return getDimension(0) == 0 && getDimension(1) == 0;
  }

  bool hasPoints() {
    return getDimension(0) == 0 || getDimension(1) == 0;
  }

  /**
   * Tests if an input geometry has edges.
   * This indicates that topology needs to be computed for them.
   *
   * @param geomIndex
   * @return true if the input geometry has edges
   */
  bool hasEdges(int geomIndex) {
    return geom[geomIndex] != null && geom[geomIndex]!.getDimension() > 0;
  }

  /**
   * Determines the location within an area geometry.
   * This allows disconnected edges to be fully
   * located.
   *
   * @param geomIndex the index of the geometry
   * @param pt the coordinate to locate
   * @return the location of the coordinate
   *
   * @see Location
   */
  int locatePointInArea(int geomIndex, Coordinate pt) {
    // Assert: only called if dimension(geomIndex) = 2

    if (isCollapsed[geomIndex]) return Location.EXTERIOR;

    //return ptLocator.locate(pt, geom[geomIndex]);

    //*
    // this check is required because IndexedPointInAreaLocator can't handle empty polygons
    if (getGeometry(geomIndex)!.isEmpty() || isCollapsed[geomIndex])
      return Location.EXTERIOR;

    PointOnGeometryLocator ptLocator = getLocator(geomIndex);
    return ptLocator.locate(pt);
    //*/
  }

  PointOnGeometryLocator getLocator(int geomIndex) {
    if (geomIndex == 0) {
      if (ptLocatorA == null) {
        ptLocatorA = IndexedPointInAreaLocator(getGeometry(geomIndex)!);
      }
      return ptLocatorA!;
    } else {
      if (ptLocatorB == null)
        ptLocatorB = IndexedPointInAreaLocator(getGeometry(geomIndex)!);
      return ptLocatorB!;
    }
  }

  void setCollapsed(int geomIndex, bool isGeomCollapsed) {
    isCollapsed[geomIndex] = isGeomCollapsed;
  }
}
