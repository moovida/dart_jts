part of dart_jts;

/**
 * Locates points on a linear geometry,
 * using a spatial index to provide good performance.
 *
 * @author mdavis
 *
 */
class IndexedPointOnLineLocator implements PointOnGeometryLocator {
  Geometry _inputGeom;

  IndexedPointOnLineLocator(this._inputGeom);

  @override
  int locate(Coordinate p) {
    // TODO: optimize this with a segment index
    PointLocator locator = PointLocator();
    return locator.locate(p, _inputGeom);
  }
}
