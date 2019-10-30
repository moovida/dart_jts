part of dart_sfs;

final _EMPTY_MULTIPOINT = MultiPoint(null);

/// A MultiPoint is a 0-dimensional GeometryCollection. The elements of a
/// MultiPoint are restricted to Points. The Points are not connected or
/// ordered in any semantically important way.
class MultiPoint extends GeometryCollection {
  /// Creates a new multipoint object [points].
  ///
  /// [points] must not include null values, otherwise throws an
  /// [ArgumentError].
  ///
  /// if [points] is null or empty, then an empty multipoint object
  /// is created.
  ///
  /// [points] don't have to be homogeneous with respect to the z- and
  /// m-coordinate. You can mix xy-, and xy{z,m}-points in a multipoint.
  /// However, [is3D] only returns true, iff all points have a z-coordinate.
  /// Similary, [isMeasured] only returns true, iff all points have an
  /// m-value.
  MultiPoint(List<Point> points) : super(points);

  /// Creates an empty multipoint object.
  factory MultiPoint.empty() => _EMPTY_MULTIPOINT;

  /// Creates a new multipoint from the WKT string [wkt].
  ///
  /// Throws a [WKTError] if [wkt] isn't a valid representation of
  /// a [MultiPoint].
  factory MultiPoint.wkt(String wkt) {
    var g = parseWKT(wkt);
    if (g is! MultiPoint) {
      throw WKTError("WKT string doesn't represent a MultiPoint");
    }
    return g;
  }

  @override
  int get dimension => 0;

  @override
  String get geometryType => "MultiPoint";

  @override
  bool get isValid => true;

  bool _isSimple;

  _computeIsSimple() {
    compare(Geometry p, Geometry q) {
      int c = (p as Point).x.compareTo((q as Point).x);
      return c != 0 ? c : (p as Point).y.compareTo((q as Point).y);
    }

    checkDuplicate(last, that) {
      if (last == null) return that; // that is the first element
      if (last == false) return false; // we already have a duplicate
      if (last.x == that.x && last.y == that.y) {
        // now we have a duplicate
        return false;
      }
      // no duplicate -> that becomes last in the next step
      return that;
    }

    if (this.isEmpty) {
      _isSimple = true;
      return;
    }
    _geometries.sort(compare);
    var ret = _geometries.fold(null, checkDuplicate);
    _isSimple = !(ret == false);
  }

  /// A MultiPoint is simple if no two points are identical.
  ///
  /// The value of this property is computed upon first access
  /// and then cached. Subsequent reads of the property
  /// efficiently reply the cached value.
  @override
  bool get isSimple {
    if (_isSimple == null) _computeIsSimple();
    return _isSimple;
  }

  /// The boundary of a [MultiPoint] is an empty [GeometryCollection]
  @override
  Geometry get boundary => GeometryCollection.empty();

  _writeTaggedWKT(writer, {bool withZ: false, bool withM: false}) {
    writer.write("MULTIPOINT");
    writer.blank();
    if (this.isNotEmpty) {
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
        if (i > 0) writer..comma();
        if (i % 10 == 0) {
          writer
            ..newline()
            ..ident();
        }
        writer.lparen();
        (elementAt(i) as Point)
            ._writeCoordinates(writer, withZ: withZ, withM: withM);
        writer.rparen();
      }
      writer..newline();
      writer
        ..decIdent()
        ..ident()
        ..rparen();
    }
  }
}
