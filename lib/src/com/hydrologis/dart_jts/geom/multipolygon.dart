part of dart_jts;

/**
 * Models a collection of {@link Polygon}s.
 * <p>
 * As per the OGC SFS specification,
 * the Polygons in a MultiPolygon may not overlap,
 * and may only touch at single points.
 * This allows the topological point-set semantics
 * to be well-defined.
 *
 *
 *@version 1.7
 */
class MultiPolygon extends GeometryCollection implements Polygonal {
  /**
   *  Constructs a <code>MultiPolygon</code>.
   *
   *@param  polygons        the <code>Polygon</code>s for this <code>MultiPolygon</code>
   *      , or <code>null</code> or an empty array to create the empty geometry.
   *      Elements may be empty <code>Polygon</code>s, but not <code>null</code>
   *      s. The polygons must conform to the assertions specified in the <A
   *      HREF="http://www.opengis.org/techno/specs.htm">OpenGIS Simple Features
   *      Specification for SQL</A> .
   *@param  precisionModel  the specification of the grid of allowable points
   *      for this <code>MultiPolygon</code>
   *@param  SRID            the ID of the Spatial Reference System used by this
   *      <code>MultiPolygon</code>
   * @deprecated Use GeometryFactory instead
   */
  MultiPolygon(List<Polygon> polygons, PrecisionModel precisionModel, int SRID)
      : this.withFactory(polygons,
            new GeometryFactory.withPrecisionModelSrid(precisionModel, SRID));

  /**
   * @param polygons
   *            the <code>Polygon</code>s for this <code>MultiPolygon</code>,
   *            or <code>null</code> or an empty array to create the empty
   *            geometry. Elements may be empty <code>Polygon</code>s, but
   *            not <code>null</code>s. The polygons must conform to the
   *            assertions specified in the <A
   *            HREF="http://www.opengis.org/techno/specs.htm">OpenGIS Simple
   *            Features Specification for SQL</A>.
   */
  MultiPolygon.withFactory(List<Polygon>? polygons, GeometryFactory factory)
      : super.withFactory(polygons, factory);

  int getDimension() {
    return 2;
  }

  int getBoundaryDimension() {
    return 1;
  }

  String getGeometryType() {
    return "MultiPolygon";
  }

  /*
   bool isSimple() {
    return true;
  }
*/

  /**
   * Computes the boundary of this geometry
   *
   * @return a lineal geometry (which may be empty)
   * @see Geometry#getBoundary
   */
  Geometry getBoundary() {
    if (isEmpty()) {
      return getFactory().createMultiLineStringEmpty();
    }
    List allRings = [];
    for (int i = 0; i < geometries.length; i++) {
      Polygon polygon = geometries[i] as Polygon;
      Geometry rings = polygon.getBoundary();
      for (int j = 0; j < rings.getNumGeometries(); j++) {
        allRings.add(rings.getGeometryN(j));
      }
    }
    return getFactory().createMultiLineString(allRings as List<LineString>);
  }

  bool equalsExactWithTol(Geometry other, double tolerance) {
    if (!isEquivalentClass(other)) {
      return false;
    }
    return super.equalsExactWithTol(other, tolerance);
  }

  /**
   * Creates a {@link MultiPolygon} with
   * every component reversed.
   * The order of the components in the collection are not reversed.
   *
   * @return a MultiPolygon in the reverse order
   */
  Geometry reverse() {
    int n = geometries.length;
    List<Polygon> revGeoms = []; //..length = n;
    for (int i = 0; i < n; i++) {
      revGeoms.add(geometries[i].reverse() as Polygon);
    }
    return getFactory().createMultiPolygon(revGeoms);
  }

  MultiPolygon copyInternal() {
    List<Polygon> polygons = []; //..length = this.geometries.length;
    for (int i = 0; i < this.geometries.length; i++) {
      polygons.add(this.geometries[i].copy() as Polygon);
    }
    return new MultiPolygon.withFactory(polygons, geomFactory);
  }

  int getSortIndex() {
    return Geometry.SORTINDEX_MULTIPOLYGON;
  }
}
