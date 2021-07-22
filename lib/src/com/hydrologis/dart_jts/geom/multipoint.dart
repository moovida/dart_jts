part of dart_jts;

/**
 * Models a collection of {@link Point}s.
 * <p>
 * Any collection of Points is a valid MultiPoint.
 *
 *@version 1.7
 */
class MultiPoint extends GeometryCollection implements Puntal {
  /**
   *  Constructs a <code>MultiPoint</code>.
   *
   *@param  points          the <code>Point</code>s for this <code>MultiPoint</code>
   *      , or <code>null</code> or an empty array to create the empty geometry.
   *      Elements may be empty <code>Point</code>s, but not <code>null</code>s.
   *@param  precisionModel  the specification of the grid of allowable points
   *      for this <code>MultiPoint</code>
   *@param  SRID            the ID of the Spatial Reference System used by this
   *      <code>MultiPoint</code>
   * @deprecated Use GeometryFactory instead
   */
  MultiPoint(List<Point> points, PrecisionModel precisionModel, int SRID)
      : super.withFactory(points,
            new GeometryFactory.withPrecisionModelSrid(precisionModel, SRID));

  /**
   *@param  points          the <code>Point</code>s for this <code>MultiPoint</code>
   *      , or <code>null</code> or an empty array to create the empty geometry.
   *      Elements may be empty <code>Point</code>s, but not <code>null</code>s.
   */
  MultiPoint.withFactory(List<Point>? points, GeometryFactory factory)
      : super.withFactory(points, factory);

  int getDimension() {
    return 0;
  }

  int getBoundaryDimension() {
    return Dimension.FALSE;
  }

  String getGeometryType() {
    return "MultiPoint";
  }

  /**
   * Gets the boundary of this geometry.
   * Zero-dimensional geometries have no boundary by definition,
   * so an empty GeometryCollection is returned.
   *
   * @return an empty GeometryCollection
   * @see Geometry#getBoundary
   */
  Geometry getBoundary() {
    return getFactory().createGeometryCollectionEmpty();
  }

  bool isValid() {
    return true;
  }

  bool equalsExactWithTol(Geometry other, double tolerance) {
    if (!isEquivalentClass(other)) {
      return false;
    }
    return super.equalsExactWithTol(other, tolerance);
  }

  /**
   *  Returns the <code>Coordinate</code> at the given position.
   *
   *@param  n  the index of the <code>Coordinate</code> to retrieve, beginning
   *      at 0
   *@return    the <code>n</code>th <code>Coordinate</code>
   */
  Coordinate? getCoordinateAt(int n) {
    return (geometries[n] as Point).getCoordinate();
  }

  MultiPoint copyInternal() {
    List<Point> points = []; //..length = this.geometries.length;
    for (int i = 0; i < this.geometries.length; i++) {
      points.add(this.geometries[i].copy() as Point);
    }
    return new MultiPoint.withFactory(points, geomFactory);
  }

  int getSortIndex() {
    return Geometry.SORTINDEX_MULTIPOINT;
  }
}
