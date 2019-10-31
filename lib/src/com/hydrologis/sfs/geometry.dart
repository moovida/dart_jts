part of dart_sfs;

/// The abstract base class for all geometries.
abstract class Geometry {
  GeometryFactory _defaultGeometryFactory = new GeometryFactory();

  int _srid;

  /// Creates a geometry from the WKT string [wkt].
  ///
  /// Throws a [WKTError] if [wkt] isn't a valid WKT geometry.
  factory Geometry.wkt(String wkt) => parseWKT(wkt);

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

  ///  Performs an operation with or on this Geometry and its
  ///  component Geometry's.  Only GeometryCollections and
  ///  Polygons have component Geometry's; for Polygons they are the LinearRings
  ///  of the shell and holes.
  ///
  ///@param  filter  the filter to apply to this <code>Geometry</code>.
  void apply(GeometryComponentFilter filter);

  /// A WKT representation of the geometry
  @specification(name: "asText()")
  String get asText {
    var buffer = StringBuffer();
    var writer = _WKTWriter(buffer);
    _writeTaggedWKT(writer, withZ: is3D, withM: isMeasured);
    return buffer.toString();
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

  /// Creates a Point using the given Coordinate.
  /// A null Coordinate creates an empty Geometry.
  ///
  /// @param coordinate a Coordinate, or null
  /// @return the created Point
  Point createPoint(Coordinate coordinate) {
    if (coordinate == null) return Point.empty();
    return Point(coordinate.x, coordinate.y);
  }
}
