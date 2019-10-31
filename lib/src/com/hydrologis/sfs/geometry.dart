part of dart_sfs;

/// The abstract base class for all geometries.
abstract class Geometry {
  GeometryFactory _defaultGeometryFactory = GeometryFactory.defaultPrecision();

  GeometryComponentFilter _geometryChangedFilter = GeometryChangedFilter();

  int _srid;

  /// Creates a geometry from the WKT string [wkt].
  ///
  /// Throws a [WKTError] if [wkt] isn't a valid WKT geometry.
  factory Geometry.wkt(String wkt) => WKTReader().read(wkt);

  /// Creates a geometry from a GeoJSON string [json].
  ///
  /// Throws a  [FormatError] if [json] isn't valid.
  factory Geometry.geojson(String json) => parseGeoJson(json);

  Geometry();

  /// Returns true if this geometric object is the empty Geometry.
  @specification(name: "isEmpty()")
  bool get isEmpty;

  /// Returns true if this geometric object has z coordinate values.
  @specification(name: "is3D()")
  bool get is3D;

  /// Returns true if this geometric object has m coordinate values.
  @specification(name: "isMeasured()")
  bool get isMeasured;

  /// Returns the closure of the combinatorial boundary of this geometric object
  @specification(name: "boundary()")
  Geometry get boundary;

  ///  Returns true if this geometric object has no anomalous geometric points,
  ///  such as self intersection or self tangency.
  @specification(name: "isSimple()")
  bool get isSimple;

  /// Returns the name of the instantiable subtype of Geometry of which this
  /// geometric object is an instantiable member. The name of the subtype of
  /// Geometry is returned as a string.
  @specification(name: "geometryType()")
  String get geometryType;

  /// The inherent dimension of this geometric object, which must be less than
  /// or equal to the coordinate dimension. In non-homogeneous collections,
  /// this will return the largest topological dimension of the contained objects.
  @specification(name: "dimension()")
  int get dimension;

  /// Returns the Spatial Reference System ID for this geometric object.
  @specification(name: "srid()")
  int get SRID => _srid;

  /// the Spatial Reference System ID for this geometric object.
  set SRID(int value) => _srid = value;

  ///  Returns a vertex of this <code>Geometry</code>
  ///  (usually, but not necessarily, the first one).
  ///  The returned coordinate should not be assumed
  ///  to be an actual Coordinate object used in
  ///  the internal representation.
  ///
  ///@return    a [Coordinate] which is a vertex of this <code>Geometry</code> or null if this Geometry is empty
  Coordinate getCoordinate();

  ///  Returns an array containing the values of all the vertices for
  ///  this geometry.
  ///  If the geometry is a composite, the array will contain all the vertices
  ///  for the components, in the order in which the components occur in the geometry.
  ///  <p>
  ///  In general, the array cannot be assumed to be the actual internal
  ///  storage for the vertices.  Thus modifying the array
  ///  may not modify the geometry itself.
  ///
  ///@return    the vertices of this <code>Geometry</code>
  List<Coordinate> getCoordinates();

  /// Returns the number of {@link Geometry}s in a {@link GeometryCollection}
  /// (or 1, if the geometry is not a collection).
  ///
  /// @return the number of geometries contained in this geometry
  int getNumGeometries();

  /// Returns an element {@link Geometry} from a {@link GeometryCollection}
  /// (or <code>this</code>, if the geometry is not a collection).
  ///
  /// @param n the index of the geometry element
  /// @return the n'th geometry contained in this geometry
  Geometry getGeometryN(int n);

  Envelope get envelope {
    if (_cachedEnvelope == null) {
      _cachedEnvelope = _computeEnvelope();
    }
    return _cachedEnvelope;
  }

  Envelope _cachedEnvelope;

  Envelope _computeEnvelope();

  GeometryFactory getFactory() {
    return _defaultGeometryFactory;
  }

  ///  Returns the <code>PrecisionModel</code> used by the <code>Geometry</code>.
  ///
  ///@return    the specification of the grid of allowable points, for this
  ///      <code>Geometry</code> and all other <code>Geometry</code>s
  PrecisionModel getPrecisionModel() {
    return _defaultGeometryFactory.getPrecisionModel();
  }

  ///  Performs an operation with or on this Geometry and its
  ///  component Geometry's.  Only GeometryCollections and
  ///  Polygons have component Geometry's; for Polygons they are the LinearRings
  ///  of the shell and holes.
  ///
  ///@param  filter  the filter to apply to this <code>Geometry</code>.
  void applyGCF(GeometryComponentFilter filter);

  void applyCSF(CoordinateSequenceFilter filter);

  CoordinateSequence getCoordinateSequence() {
    return CoordinateArraySequence(getCoordinates());
  }

  /**
   * Notifies this geometry that its coordinates have been changed by an external
   * party (for example, via a {@link CoordinateFilter}).
   * When this method is called the geometry will flush
   * and/or update any derived information it has cached (such as its {@link Envelope} ).
   * The operation is applied to all component Geometries.
   */
  void geometryChanged() {
    applyGCF(_geometryChangedFilter);
  }

  /**
   * Notifies this Geometry that its Coordinates have been changed by an external
   * party. When #geometryChanged is called, this method will be called for
   * this Geometry and its component Geometries.
   *
   * @see #apply(GeometryComponentFilter)
   */
  void geometryChangedAction() {
    _cachedEnvelope = null;
  }

  ///  Returns the Well-known Text representation of this <code>Geometry</code>.
  ///  For a definition of the Well-known Text format, see the OpenGIS Simple
  ///  Features Specification.
  ///
  ///@return    the Well-known Text representation of this <code>Geometry</code>
  @specification(name: "toText()")
  String toText() {
    WKTWriter writer = WKTWriter();
    return writer.write(this);
  }

  _writeTaggedWKT(writer, {bool withZ: false, bool withM: false});

  /* ------------------------------------------------------------------- */
  //TODO: implement these
  /// Returns true if this geometric object is “spatially equal” to
  /// [other].
  @specification(name: "equals()")
  bool equals(Geometry other) {
    throw UnimplementedError();
  }

  /// Returns true if this geometric object is “spatially disjoint”
  /// from [other].
  @specification(name: "disjoint()")
  bool disjoint(Geometry other) {
    throw UnimplementedError();
  }

  /// Returns true if this geometric object “spatially intersects”
  /// [other].
  @specification(name: "intersects()")
  bool intersects(Geometry other) {
    throw UnimplementedError();
  }

  /// Returns true if this geometric object is “spatially touches”
  /// [other].
  @specification(name: "touches()")
  bool touches(Geometry other) {
    throw UnimplementedError();
  }

  /// Returns true if this geometric object is “spatially crosses”
  /// [other].
  @specification(name: "crosses()")
  bool crosses(Geometry other) {
    throw UnimplementedError();
  }

  /// Returns true if this geometric object is “spatially within”
  /// [other].
  @specification(name: "within()")
  bool within(Geometry other) {
    throw UnimplementedError();
  }

  /// Returns true if this geometric object is “spatially contains”
  /// [other].
  @specification(name: "contains()")
  bool contains(Geometry other) {
    throw UnimplementedError();
  }

  /// Returns true if this geometric object is “spatially overlaps”
  /// with [other].
  @specification(name: "overlaps()")
  bool overlaps(Geometry other) {
    throw UnimplementedError();
  }

  ///  Returns true if this geometric object is spatially related to
  ///  [other] by testing for intersections between the interior, boundary
  ///  and exterior of the two geometric objects as specified by the values
  ///  in the [pattern].
  ///
  ///  This returns false if all the tested intersections are empty except
  ///  exterior (this) intersect exterior (another).
  @specification(name: "relate()")
  bool relate(Geometry other, String pattern) {
    throw UnimplementedError();
  }
}

class GeometryFactory {
  CoordinateSequenceFactory _coordinateSequenceFactory;

  PrecisionModel _precisionModel;
  int _SRID;

  /// Gets the SRID value defined for this factory.
  ///
  /// @return the factory SRID value
  int getSRID() {
    return _SRID;
  }

  /// Constructs a GeometryFactory that generates Geometries having the given
  /// PrecisionModel, spatial-reference ID, and CoordinateSequence implementation.
  GeometryFactory(PrecisionModel precisionModel, int SRID, CoordinateSequenceFactory coordinateSequenceFactory) {
    this._precisionModel = precisionModel;
    this._coordinateSequenceFactory = coordinateSequenceFactory;
    this._SRID = SRID;
  }

  /// Constructs a GeometryFactory that generates Geometries having the given
  /// CoordinateSequence implementation, a double-precision floating PrecisionModel and a
  /// spatial-reference ID of 0.
  GeometryFactory.withCoordinateSequenceFactory(CoordinateSequenceFactory coordinateSequenceFactory) : this(PrecisionModel(), 0, coordinateSequenceFactory);

  /// Constructs a GeometryFactory that generates Geometries having the given
  /// {@link PrecisionModel} and the default CoordinateSequence
  /// implementation.
  ///
  /// @param precisionModel the PrecisionModel to use
  GeometryFactory.withPrecisionModel(PrecisionModel precisionModel) : this(precisionModel, 0, getDefaultCoordinateSequenceFactory());

  /// Constructs a GeometryFactory that generates Geometries having the given
  /// {@link PrecisionModel} and spatial-reference ID, and the default CoordinateSequence
  /// implementation.
  ///
  /// @param precisionModel the PrecisionModel to use
  /// @param SRID the SRID to use
  GeometryFactory.withPrecisionModelSrid(PrecisionModel precisionModel, int SRID) : this(precisionModel, SRID, getDefaultCoordinateSequenceFactory());

  /// Constructs a GeometryFactory that generates Geometries having a floating
  /// PrecisionModel and a spatial-reference ID of 0.
  GeometryFactory.defaultPrecision() : this.withPrecisionModelSrid(PrecisionModel(), 0);

  static CoordinateSequenceFactory getDefaultCoordinateSequenceFactory() {
    return CoordinateArraySequenceFactory();
  }

  /// Creates a {@link MultiPoint} using the given {@link Point}s.
  /// A null or empty array will create an empty MultiPoint.
  ///
  /// @param point an array of Points (without null elements) or an empty array, or <code>null</code>
  /// @return a MultiPoint object
  MultiPoint createMultiPoint([List<Point> points]) {
    if (points == null || points.isEmpty) return MultiPoint.empty();
    return MultiPoint(points);
  }

  MultiPoint createMultiPointFromCoords([List<Coordinate> coords]) {
    if (coords == null || coords.isEmpty) return MultiPoint.empty();
    return MultiPoint(coords.map((c) => Point(c.x, c.y)));
  }

  LineString createLineString(List<Coordinate> coordinates) {
    return LineString.fromCoordinates(coordinates);
  }

  LinearRing createLinearRing(List<Coordinate> coordinates) {
    return LinearRing.fromCoordinates(coordinates);
  }

  MultiLineString createMultiLineString(List<LineString> lines) {
    return MultiLineString(lines);
  }

  /// Creates a Point using the given Coordinate.
  /// A null Coordinate creates an empty Geometry.
  ///
  /// @param coordinate a Coordinate, or null
  /// @return the created Point
  Point createPoint(Coordinate coordinate) {
    if (coordinate == null) return Point.empty();
    return Point(coordinate.x, coordinate.y, z: coordinate.z, m: coordinate.getM());
  }

  Polygon createPolygon(LinearRing ring) {
    return Polygon(ring, null);
  }

  /// Creates a {@link Geometry} with the same extent as the given envelope.
  /// The Geometry returned is guaranteed to be valid.
  /// To provide this behaviour, the following cases occur:
  /// <p>
  /// If the <code>Envelope</code> is:
  /// <ul>
  /// <li>null : returns an empty {@link Point}
  /// <li>a point : returns a non-empty {@link Point}
  /// <li>a line : returns a two-point {@link LineString}
  /// <li>a rectangle : returns a {@link Polygon} whose points are (minx, miny),
  ///  (minx, maxy), (maxx, maxy), (maxx, miny), (minx, miny).
  /// </ul>
  ///
  ///@param  envelope the <code>Envelope</code> to convert
  ///@return an empty <code>Point</code> (for null <code>Envelope</code>s),
  ///	a <code>Point</code> (when min x = max x and min y = max y) or a
  ///      <code>Polygon</code> (in all other cases)
  Geometry toGeometry(Envelope envelope) {
    // null envelope - return empty point geometry
    if (envelope.isNull()) {
      return Point.empty();
    }

    // point?
    if (envelope.getMinX() == envelope.getMaxX() && envelope.getMinY() == envelope.getMaxY()) {
      return createPoint(Coordinate(envelope.getMinX(), envelope.getMinY()));
    }

    // vertical or horizontal line?
    if (envelope.getMinX() == envelope.getMaxX() || envelope.getMinY() == envelope.getMaxY()) {
      return createLineString([Coordinate(envelope.getMinX(), envelope.getMinY()), Coordinate(envelope.getMaxX(), envelope.getMaxY())]);
    }

    // create a CW ring for the polygon
    return createPolygon(createLinearRing([
      Coordinate(envelope.getMinX(), envelope.getMinY()),
      Coordinate(envelope.getMinX(), envelope.getMaxY()),
      Coordinate(envelope.getMaxX(), envelope.getMaxY()),
      Coordinate(envelope.getMaxX(), envelope.getMinY()),
      Coordinate(envelope.getMinX(), envelope.getMinY())
    ]));
  }

  CoordinateSequenceFactory getCoordinateSequenceFactory() {
    return _coordinateSequenceFactory;
  }

  /// Returns the PrecisionModel that Geometries created by this factory
  /// will be associated with.
  ///
  /// @return the PrecisionModel for this factory
  PrecisionModel getPrecisionModel() {
    return _precisionModel;
  }

  MultiPolygon createMultiPolygon(List<Polygon> polygons) {
    return MultiPolygon(polygons);
  }

  GeometryCollection createGeometryCollection(List<Geometry> geometries) {
    return GeometryCollection(geometries);
  }
}
