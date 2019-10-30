part of simple_features;

final _EMPTY_MULTI_LINE_STRING = MultiLineString(null);

/// A MultiLineString is a MultiCurve whose elements are [LineString]s.
class MultiLineString extends GeometryCollection {
  /// Creates a multilinestring for a collection of [linestrings].
  ///
  /// If [linestrings] is null or empty, an empty [MultiLineString]
  /// is created.
  MultiLineString(Iterable<LineString> linestrings) : super(linestrings);

  /// Creates an empty multilinestring.
  factory MultiLineString.empty() => _EMPTY_MULTI_LINE_STRING;

  /// Creates a new multilinestring from the WKT string [wkt].
  ///
  /// Throws a [WKTError] if [wkt] isn't a valid representation of
  /// a [MultiLineString].
  factory MultiLineString.wkt(String wkt) {
    var g = parseWKT(wkt);
    if (g is! MultiLineString) {
      throw WKTError("WKT string doesn't represent a MultiLineString");
    }
    return g;
  }

  @override
  int get dimension => 1;

  @override
  String get geometryType => "MultiLineString";

  /// This multilinestring is closed if all child line strings are
  /// closed.
  bool get isClosed => _geometries.every((g) => (g as LineString).isClosed);

  /// Replies the spatial length of this multilinestring.
  @specification(name: "length()")
  num get spatialLength {
    throw UnimplementedError();
  }

  /// The boundary of a [MultiLineString] consists of the boundary
  /// points of the child geometries which occur an odd number of
  /// times in the boundaries.
  @override
  Geometry get boundary {
    var pointRefCounts = Map<DirectPosition2D, int>();
    countPosition(pos) {
      if (pointRefCounts.containsKey(pos)) {
        pointRefCounts[pos] = pointRefCounts[pos] + 1;
      } else {
        pointRefCounts[pos] = 1;
      }
    }

    // count the number of occurences for each boundary point
    forEach((child) {
      if (child.isEmpty) return;
      child.boundary.forEach((p) {
        countPosition(DirectPosition2D(p.x, p.y));
      });
    });

    // boundary points with odd occurences in the child boundaries
    // are considered to be boundary points of this MultiLineString
    // too
    var points = [];
    pointRefCounts.forEach((pos, count) {
      if (count % 2 == 0) return;
      points.add(Point(pos.x, pos.y));
    });
    return MultiPoint(points);
  }

  _writeTaggedWKT(writer, {bool withZ: false, bool withM: false}) {
    writer.write("MULTILINESTRING");
    writer.blank();
    if (this.isNotEmpty) {
      writer.ordinateSpecification(withZ: withZ, withM: withM);
    }
    if (isEmpty) {
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
        (elementAt(i) as LineString)
            ._writeWKT(writer, withZ: withZ, withM: withM);
      }
      writer..newline();
      writer
        ..decIdent()
        ..ident()
        ..rparen();
    }
  }
}
