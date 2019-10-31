part of dart_sfs;

/// the singleton empty line string
final _EMPTY_LINESTRING = LineString(null);

/// A LineString is a curve with linear interpolation between points.
///
class LineString extends Geometry
    with IterableMixin<Point>, _GeometryContainerMixin {
  List<Point> _points;

  _init(points) {
    if (points == null || points.isEmpty) {
      _points = null;
    } else {
      _require(
          points.length >= 2, "illegal number of points, got ${points.length}");
      _require(points.every((p) => p != null && !p.isEmpty),
          "points must not contain null values or empty points");
      for (int i = 1; i < points.length; i++) {
        _require(!points[i - 1].equals2D(points[i]),
            "idenical consequtive point ${points[i - 1]} at index ${i - 1}");
      }
      _points = List.from(points, growable: false);
    }
  }

  /// Creates a new linestring.
  ///
  /// Creates an empty linestring if [points] is null or empty.
  ///
  /// Throws an [ArgumentError] if
  /// * [points] contains only one point
  /// * or if it contains null values or empty points
  /// * or if it contains a run of 2 or more consequtive, identical
  ///   points
  LineString(List<Point> points) {
    _init(points);
  }

  /// Creates a new line.
  ///
  /// A line is a linestring with exactly two non-null, non-empty,
  /// non-equal [points].
  ///
  /// Throws an [ArgumentError] if one of these conditions isn't met.
  LineString.line(List<Point> points) {
    _require(points != null);
    _require(points.length == 2);
    _require(points.every((p) => p != null && !p.isEmpty));
    if (points.length == 2) {
      _require(!points.first.equals2D(points.last));
    }
    _points = List.from(points, growable: false);
  }

  /// Creates a new linear ring.
  ///
  /// A linear ring is a linestring that is both closed and simple.
  /// As a consequence, it must have at least four nodes.
  ///
  /// Throws an [ArgumentError] if the [points] passed in do not represent
  /// a linear ring.
  LineString.ring(List<Point> points) {
    _init(points);
    _require(this.isNotEmpty, "a ring can't be empty");
    _require(this.length >= 4, "a ring must have at least four nodes");
    _require(this.isClosed, "a ring must be closed");
    _require(this._checkIsSimple(), "a ring must be simple");
  }

  LineString.fromCoordinates(List<Coordinate> coordinates) {
    _require(coordinates != null);
    _require(coordinates.length == 2);
    _require(coordinates.every((p) => p != null));
    _points = coordinates.map((c) => Point(c.x, c.y)).toList();
  }

  /// Creates a new linestring from the WKT string [wkt].
  ///
  /// Throws a [WKTError] if [wkt] isn't a valid representation of
  /// a [LineString].
  factory LineString.wkt(String wkt) {
    var g = parseWKT(wkt);
    if (g is! LineString) {
      throw WKTError("WKT string doesn't represent a LineString");
    }
    return g;
  }

  /// Creates an empty linestring.
  factory LineString.empty() => _EMPTY_LINESTRING;

  Iterator<Point> get iterator =>
      _points == null ? <Point>[].iterator : _points.iterator;

  @override
  Coordinate getCoordinate() {
    return _points.isEmpty ? null : _points[0].toCoordinate();
  }

  @override
  List<Coordinate> getCoordinates() {
    return _points.isEmpty ? [] : _points.map((p) => p.toCoordinate()).toList();
  }

  int getNumGeometries() {
    return 1;
  }

  Geometry getGeometryN(int n) {
    return this;
  }

  Coordinate getCoordinateN(int n) {
    return elementAt(n).toCoordinate();
  }

  void apply(GeometryComponentFilter filter) {
    filter.filter(this);
  }

  @override
  Envelope _computeEnvelope() {
    if (isEmpty) return Envelope.empty();

    Envelope e = Envelope.empty();
    _points.forEach((p) {
      e.growTo(p.toCoordinate());
    });
    return e;
  }

  @override
  String get geometryType => "LineString";

  @override
  int get dimension => 1;

  /// Replies the number of points in this linestring.
  ///
  /// See also [length]
  @specification(name: "numPoints()")
  int numPoints() => length;

  /// Replies the n-th point this linestring.
  ///
  /// See also [elementAt]
  @specification(name: "pointN()")
  Point pointN(int n) => elementAt(n);

  /// Replies the (spatial) length of this line string.
  @specification(name: "length()")
  //TODO: implement
  num get spatialLength {
    throw UnimplementedError();
  }

  /// Replies the start point of this linestring.
  ///
  /// See also the Dart'ish property [first].
  ///
  /// Throws a [StateError] if this linestring is empty.
  @specification(name: "StartPoint()")
  Point get startPoint => first;

  /// Replies the end point of this linestring.
  ///
  /// See also the Dart'ish property [last].
  ///
  /// Throws a [StateError] if this linestring is empty.
  @specification(name: "EndPoint()")
  Point get endPoint => last;

  /// Replies true if this linestring isn't empty and its
  /// first and last points are equal (with respect to the xy-coordinates)
  bool get isClosed {
    if (this.isEmpty) return false;
    return first.equals2D(last);
  }

  /// Replies true if this linestring is closed and simple.
  ///
  //TODO: implement
  bool get isRing {
    throw UnimplementedError();
  }

  _writeWKT(writer, {bool withZ: false, bool withM: false}) {
    if (isEmpty) {
      writer.empty();
    } else {
      writer.lparen();
      for (int i = 0; i < length; i++) {
        if (i > 0) {
          writer.comma();
          writer.blank();
        }
        this.elementAt(i)._writeCoordinates(writer, withZ: withZ, withM: withM);
      }
      writer.rparen();
    }
  }

  _writeTaggedWKT(writer, {bool withZ: false, bool withM: false}) {
    writer.write("LINESTRING");
    writer.blank();
    if (!isEmpty) {
      writer.ordinateSpecification(withZ: withZ, withM: withM);
    }
    _writeWKT(writer, withZ: is3D, withM: isMeasured);
  }

  /// The boundary is empty if this linestring [isEmpty]
  /// or [isClosed]. Otherwise it consists of a
  /// [MultiPoint] with the two end points.
  @override
  Geometry get boundary {
    if (this.isEmpty) return GeometryCollection.empty();
    if (this.isClosed) return GeometryCollection.empty();
    return BoundaryOp(this).getBoundary();
  }

  bool _checkIsSimple() {
    //TODO: check spec - an empty linestring is always simple?
    if (isEmpty) return true;
    var pos = _points.map((p) => p.toCoordinate()).toList(growable: false);
    var segments = List.generate(
        length - 1, (i) => LineSegment(pos[i], pos[i + 1]),
        growable: false);

    var intersections = computeLineIntersections(segments);
    return intersections.every((intersection) {
      // only linesegment points can be intersections in order to be
      // simple
      if (!pos.contains(intersection.pos)) return false;
      // max. 2 segments can intersect in any point in order to be
      // simple
      if (intersection.intersecting.length > 2) return false;
      return true;
    });
  }

  /// This linesegment is simple if it doesn't have self intersections.
  ///
  /// An empty linestring is simple.
  //TODO: cache value for isSimple?
  @override
  bool get isSimple => _checkIsSimple();

  Geometry getStartPoint() {
    if (isEmpty) return null;
    return first;
  }

  getEndPoint() {
    if (isEmpty) return null;
    return last;
  }
}

/// A Line is a [LineString] with exactly 2 [Point]s.
class Line extends LineString {
  Line(List<Point> points) : super.line(points);

  // a line is always simple
  @override
  bool get isSimple => true;
}

/// A [LinearRing] is a [LineString] that is both closed and simple.
class LinearRing extends LineString {
  LinearRing(List<Point> points) : super.ring(points);

  // a ring is always simple, simplicity is enforced in the constructor
  @override
  bool get isSimple => true;
}
