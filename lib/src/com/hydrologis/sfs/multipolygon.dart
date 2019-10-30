part of dart_sfs;

/// A MultiSurface is a 2-dimensional [GeometryCollection] whose elements
/// are [Surface]s, all using coordinates from the same coordinate reference
/// system.
abstract class MultiSurface extends GeometryCollection {
  MultiSurface(Iterable<Surface> surfaces) : super(surfaces);

  /// The mathematical centroid for this [Surface] as a [Point]. The result
  /// is not guaranteed to be on this [Surface].
  ///
  @specification(name: "centroid()")
  Point get centroid;

  /// A [Point] guaranteed to be on this [Surface].
  @specification(name: "pointOnSurface()")
  Point get pointOnSurface;

  /// The area of this [Surface], as measured in the spatial reference system
  /// of this [Surface].
  @specification(name: "area()")
  double get area;
}

final _EMPTY_MULTI_POLYGON = MultiPolygon(null);

/// A MultiPolygon is a MultiSurface whose elements are [Polygon]s.
class MultiPolygon extends MultiSurface {
  /// Creates a multipolygon.
  ///
  /// For polygons, a set of geometric invariants should hold, see
  /// the SFS:
  ///
  /// * the interiors of two polygons in this multi polygon may not intersect
  /// * the boundaries of two polygons may not "cross" and if they touch, then
  ///   only at a finite number of points
  /// * a multipolygon must not have cut lines, spikes or punctures
  ///
  /// Note, that none of these invariants is currently enforced when a
  /// polygon is created.
  MultiPolygon(Iterable<Polygon> polygons) : super(polygons);

  /// Creates an empty multipolygon.
  factory MultiPolygon.empty() => _EMPTY_MULTI_POLYGON;

  /// Creates a new multipolygon from the WKT string [wkt].
  ///
  /// Throws a [WKTError] if [wkt] isn't a valid representation of
  /// a [MultiPolygon].
  factory MultiPolygon.wkt(String wkt) {
    var g = parseWKT(wkt);
    if (g is! MultiPolygon) {
      throw WKTError("WKT string doesn't represent a MultiPolygon");
    }
    return g;
  }

  @override
  Geometry get boundary {
    if (this.isEmpty) return MultiLineString.empty();

    List<LineString> lines = [];
    _geometries.forEach((g) {
      MultiLineString mls = g.boundary as MultiLineString;
      lines.addAll(mls._geometries as List<LineString>);
    });
    return MultiLineString(lines);
  }

  @override
  int get dimension => 2;

  @override
  String get geometryType => "MultiPolygon";

  @override
  Point get centroid {
    throw UnimplementedError();
  }

  @override
  Point get pointOnSurface {
    throw UnimplementedError();
  }

  @override
  double get area {
    throw UnimplementedError();
  }

  @override
  _writeTaggedWKT(writer, {bool withZ: false, bool withM: false}) {
    writer.write("MULTIPOLYGON");
    writer.blank();
    if (!this.isEmpty) {
      writer.ordinateSpecification(withZ: withZ, withM: withM);
    }
    if (this.isEmpty) {
      writer.empty();
    } else {
      writer
        ..lparen()
        ..newline();
      writer
        ..incIdent()
        ..ident();
      for (int i = 0; i < length; i++) {
        if (i > 0) {
          writer
            ..comma()
            ..newline()
            ..ident();
        }
        (elementAt(i) as Polygon)._writeWKT(writer, withZ: withZ, withM: withM);
      }
      writer..newline();
      writer
        ..decIdent()
        ..ident()
        ..rparen();
    }
  }
}
