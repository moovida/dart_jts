part of dart_jts;

/**
 * Models a collection of {@link LineString}s.
 * <p>
 * Any collection of LineStrings is a valid MultiLineString.
 *
 *@version 1.7
 */
class MultiLineString extends GeometryCollection implements Lineal {
  /**
   *  Constructs a <code>MultiLineString</code>.
   *
   *@param  lineStrings     the <code>LineString</code>s for this <code>MultiLineString</code>
   *      , or <code>null</code> or an empty array to create the empty geometry.
   *      Elements may be empty <code>LineString</code>s, but not <code>null</code>
   *      s.
   *@param  precisionModel  the specification of the grid of allowable points
   *      for this <code>MultiLineString</code>
   *@param  SRID            the ID of the Spatial Reference System used by this
   *      <code>MultiLineString</code>
   * @deprecated Use GeometryFactory instead
   */
  MultiLineString(
      List<LineString> lineStrings, PrecisionModel precisionModel, int SRID)
      : super.withFactory(lineStrings,
            new GeometryFactory.withPrecisionModelSrid(precisionModel, SRID));

  /**
   * @param lineStrings
   *            the <code>LineString</code>s for this <code>MultiLineString</code>,
   *            or <code>null</code> or an empty array to create the empty
   *            geometry. Elements may be empty <code>LineString</code>s,
   *            but not <code>null</code>s.
   */
  MultiLineString.withFactory(
      List<LineString>? lineStrings, GeometryFactory factory)
      : super.withFactory(lineStrings, factory);

  int getDimension() {
    return 1;
  }

  int getBoundaryDimension() {
    if (isClosed()) {
      return Dimension.FALSE;
    }
    return 0;
  }

  String getGeometryType() {
    return "MultiLineString";
  }

  bool isClosed() {
    if (isEmpty()) {
      return false;
    }
    for (int i = 0; i < geometries.length; i++) {
      if (!(geometries[i] as LineString).isClosed()) {
        return false;
      }
    }
    return true;
  }

  /**
   * Gets the boundary of this geometry.
   * The boundary of a lineal geometry is always a zero-dimensional geometry (which may be empty).
   *
   * @return the boundary geometry
   * @see Geometry#getBoundary
   */
  Geometry getBoundary() {
    return (new BoundaryOp(this)).getBoundary();
  }

  /**
   * Creates a {@link MultiLineString} in the reverse
   * order to this object.
   * Both the order of the component LineStrings
   * and the order of their coordinate sequences
   * are reversed.
   *
   * @return a {@link MultiLineString} in the reverse order
   */
  Geometry reverse() {
    // int nLines = geometries.length;
    List<LineString> revLines = []; //..length = nLines;
    geometries.reversed.forEach((geom) {
      revLines.add(geom.reverse() as LineString);
    });
    // for (int i = 0; i < geometries.length; i++) {
    //   revLines[nLines - 1 - i] = geometries[i].reverse() as LineString;
    // }
    return getFactory().createMultiLineString(revLines);
  }

  MultiLineString copyInternal() {
    List<LineString> lineStrings = []; //..length = this.geometries.length;
    for (int i = 0; i < this.geometries.length; i++) {
      // lineStrings[i] = this.geometries[i].copy() as LineString;
      lineStrings.add(this.geometries[i].copy() as LineString);
    }
    return new MultiLineString.withFactory(lineStrings, geomFactory);
  }

  bool equalsExact(Geometry other, double tolerance) {
    if (!isEquivalentClass(other)) {
      return false;
    }
    return super.equalsExactWithTol(other, tolerance);
  }

  int getSortIndex() {
    return Geometry.SORTINDEX_MULTILINESTRING;
  }
}
