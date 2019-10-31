part of dart_sfs;

/// A Polygon is a planar Surface defined by 1 exterior boundary and 0 or
/// more interior boundaries. Each interior boundary defines a hole in the
/// Polygon.
class Polygon extends Surface {
  LineString _exterior;
  List<LineString> _interiors;

  /// Creates a new polygon.
  ///
  /// [exteriorRing] must not be null. If [exteriorRing] is empty,
  /// [interiorRings] must be null, the empty list, or consist of
  /// empty [LinearRing]s only.
  ///
  /// [interiorRings] can be null. If it is a non-empty list, then it
  /// must not contain null values.
  ///
  /// Note: According to the SFS, both the exterior and the interior must
  /// be valid [LinearRing]s, i.e. they must be simple, closed [LineString]s.
  /// Exterior and interior rings must not cross and a polygon must not have
  /// cuts, spikes, punctures (see the SFS for examples). These restrictions
  /// are currently not enforced when a [Polygon] is created.
  ///
  /// Throws [ArgumentError] if one of the preconditions is violated.
  Polygon(LineString exteriorRing, Iterable<LineString> interiorRings) {
    _require(exteriorRing != null);
    //TODO: valdiate that exteriorRing is indeed a ring
    if (interiorRings == null) interiorRings = [];
    _require(interiorRings.every((r) => r != null), "interior rings must not be null");
    if (exteriorRing.isEmpty) {
      _require(interiorRings.every((r) => r.isEmpty), "exterior ring is empty => all interior rings must be empty");
      _exterior = null;
      _interiors = null;
      return;
    }
    _exterior = exteriorRing;
    _interiors = List.from(interiorRings, growable: false);
    //TODO: check geometry of interiors. Must be valid rings, mutually not
    // intersection
    //TODO: check geometry of overall polygon. SFS mentions a set of edge
    // cases like spikes, cuts or punctures which are not allowed
  }

  /// Creates a new empty polygon.
  Polygon.empty()
      : _exterior = null,
        _interiors = null;

  /// Creates a new point from the WKT string [wkt].
  ///
  /// Throws a [WKTError] if [wkt] isn't a valid representation of
  /// a [Polygon].
  factory Polygon.wkt(String wkt) {
    var g = WKTReader().read(wkt);
    if (g is! Polygon) {
      throw ArgumentError("WKT string doesn't represent a Polygon");
    }
    return g;
  }

  /// Creates a triangle with the [exterior].
  ///
  /// [exterior] must be a non-null, closed linestring with exactly three
  /// distinct, non-colienar points.
  ///
  /// Throws [ArgumentError] if the preconditions aren't met.
  Polygon.triangle(LineString exterior) {
    _require(exterior != null);
    _require(exterior.map((p) => Coordinate(p.x, p.y)).toSet().length == 3, "a triangle must consist of three non-colinear nodes");
    _require(exterior.isClosed, "the exterior of a triangle must be closed");
    //TODO: check for colienarity of the three points

    _exterior = exterior;
  }

  @override
  Coordinate getCoordinate() {
    return _exterior.isEmpty ? null : _exterior.getCoordinate();
  }

  @override
  List<Coordinate> getCoordinates() {
    List<Coordinate> coords = [];
    coords.addAll(_exterior.getCoordinates());
    _interiors.forEach((inter) {
      coords.addAll(inter.getCoordinates());
    });
    return coords;
  }

  int getNumGeometries() {
    return 1;
  }

  Geometry getGeometryN(int n) {
    return this;
  }

  void applyGCF(GeometryComponentFilter filter) {
    filter.filter(this);
    _exterior.applyGCF(filter);
    _interiors.forEach((g) => g.applyGCF(filter));
  }

  void applyCSF(CoordinateSequenceFilter filter) {
    _exterior.applyCSF(filter);
    if (!filter.isDone()) {
      for (int i = 0; i < _interiors.length; i++) {
        _interiors[i].applyCSF(filter);
        if (filter.isDone()) {
          break;
        }
      }
    }
    if (filter.isGeometryChanged()) {
      geometryChanged();
    }
  }

  /// The exterior ring of this polygon.
  ///
  /// Replies an empty linestring if this polygon is empty.
  @specification(name: "exteriorRing()")
  LineString get exteriorRing => _exterior == null ? LineString.empty() : _exterior;

  /// The interior rings. Replies an empty iterable, if
  /// this polygon doesn't have interior rings.
  Iterable<LineString> get interiorRings => _interiors == null ? [] : _interiors;

  /// the number of interior rings
  @specification(name: "numInteriorRing()")
  int get numInteriorRing => _interiors == null ? 0 : _interiors.length;

  /// Replies the n-th interior ring.
  ///
  /// Throws a [RangeError] if [n] is out of range
  @specification(name: "interiorRingN()")
  LineString interiorRingN(int n) => interiorRings.elementAt(n);

  @override
  String get geometryType => "Polygon";

  /// Replies true if this polygon isn't empty and if both the
  /// exterior and each of the interior rings are 3D.
  @override
  bool get is3D {
    bool ret = _exterior == null ? false : _exterior.is3D;
    if (interiorRings.isEmpty) return ret;
    ret = ret && interiorRings.every((r) => r.is3D);
    return ret;
  }

  /// Replies true if this polygon isn't empty and if both the
  /// exterior and each of the interior rings are measured.
  @override
  bool get isMeasured {
    bool ret = _exterior == null ? false : _exterior.isMeasured;
    if (interiorRings.isEmpty) return ret;
    ret = ret && interiorRings.every((r) => r.isMeasured);
    return ret;
  }

  @override
  bool get isEmpty => _exterior == null;

  @override
  Geometry get boundary {
    if (isEmpty) return MultiLineString.empty();
    if (interiorRings.isEmpty) return exteriorRing;
    List<LineString> rings = [];
    rings.add(exteriorRing);
    rings.addAll(interiorRings);
    return MultiLineString(rings);
  }

  @override
  _writeTaggedWKT(writer, {bool withZ = false, bool withM = false}) {
    writer.write("POLYGON");
    writer.blank();
    if (!isEmpty) {
      writer.ordinateSpecification(withZ: is3D, withM: isMeasured);
    }
    _writeWKT(writer, withZ: withZ, withM: withM);
  }

  _writeWKT(writer, {bool withZ = false, bool withM = false}) {
    if (this.isEmpty) {
      writer.empty();
    } else {
      writer
        ..lparen()
        ..newline();
      writer
        ..incIdent()
        ..ident();
      _exterior._writeWKT(writer, withZ: withZ, withM: withM);
      if (_interiors.isNotEmpty) {
        writer
          ..comma()
          ..newline();
      }
      for (int i = 0; i < _interiors.length; i++) {
        if (i > 0) {
          writer
            ..comma()
            ..newline()
            ..ident();
        }
        _interiors[i]._writeWKT(writer, withZ: withZ, withM: withM);
      }
      writer..newline();
      writer
        ..decIdent()
        ..ident()
        ..rparen();
    }
  }

  @override
  Envelope _computeEnvelope() {
    if (this.isEmpty) return Envelope.empty();
    Envelope e = Envelope.empty();
    e.expandToIncludeEnvelope(_exterior.envelope);
    return e;
  }

  @override
  // TODO: implement area
  double get area => null;

  @override
  // TODO: implement centroid
  Point get centroid => null;

  @override
  // TODO: implement isSimple
  bool get isSimple => null;

  @override
  // TODO: implement pointOnSurface
  Point get pointOnSurface => null;
}

/// A Triangle is a polygon with 3 distinct, non-collinear vertices and no
/// interior boundary.
///
class Triangle extends Polygon {
  /// Creates a triangle with the [exterior].
  ///
  /// [exterior] must be a non-null, closed linestring with exactly three
  /// distinct, non-colienar points.
  ///
  /// Throws [ArgumentError] if the preconditions aren't met.
  Triangle(LineString exterior) : super.triangle(exterior);
}
