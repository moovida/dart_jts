part of dart_sfs;

/// A PolyhedralSurface is a contiguous collection of polygons, which share
/// common boundary segments.
class PolyhedralSurface extends Geometry
    with IterableMixin<Polygon>, _GeometryContainerMixin {
  List<Polygon> _patches;

  /// Creates a polyhedral surface with the polygons [patches].
  ///
  /// If [patches] is null or empty, an empty polyhedral surface is
  /// created.
  PolyhedralSurface(List<Polygon> patches) {
    if (patches == null || patches.isEmpty) {
      _patches = null;
      return;
    }
    _require(patches.every((p) => p != null), "polygons must not be null");
    _patches = List.from(patches, growable: false);
  }

  /// Creates an empty polyhedral surface.
  PolyhedralSurface.empty() : this(null);

  /// Creates a new polyhedralsurface from the WKT string [wkt].
  ///
  /// Throws a [WKTError] if [wkt] isn't a valid representation of
  /// a [PolyhedralSurface].
  factory PolyhedralSurface.wkt(String wkt) {
    var g = parseWKT(wkt);
    if (g is! PolyhedralSurface) {
      throw WKTError("WKT string doesn't represent a PolyhedralSurface");
    }
    return g;
  }

  /// Returns the number of including polygons
  ///
  /// See also [length]
  @specification(name: "numPatches()")
  int get numPatches => length;

  /// From the SFS: "Returns a polygon in this surface, the order is arbitrary."
  ///
  /// See also [elementAt]
  @specification(name: "patchN()")
  Polygon patchN(int n) => elementAt(n);

  /// Returns the collection of polygons in this surface that bounds the given
  /// polygon “p” for any polygon “p” in the surface.
  @specification(name: "boundingPolygons()")
  MultiPolygon get boundingPolygons {
    throw UnimplementedError();
  }

  /// Returns true  if the polygon closes on itself, and thus has no boundary
  /// and encloses a solid
  @specification(name: "isClosed()")
  bool get isClosed {
    throw UnimplementedError();
  }

  Iterator<Polygon> get iterator =>
      _patches == null ? [].iterator : _patches.iterator;

  @override
  int get dimension => 2;

  @override
  String get geometryType => "PolyhedralSurface";

  String get _wktName => "POLYHEDRALSURFACE";

  _writeTaggedWKT(writer, {bool withZ = false, bool withM = false}) {
    writer.write(_wktName);
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
        if (i > 0) {
          writer
            ..comma()
            ..newline()
            ..ident();
        }
        elementAt(i)._writeWKT(writer, withZ: withZ, withM: withM);
      }
      writer..newline();
      writer
        ..decIdent()
        ..ident()
        ..rparen();
    }
  }

  @override
  _Envelope _computeEnvelope() {
    if (isEmpty) return _Envelope.empty();

    _Envelope e = _Envelope.empty();
    _patches.forEach((p) {
      e.growTo(p);
    });
    return e;
  }

  @override
  // TODO: implement boundary
  Geometry get boundary => null;

  @override
  Coordinate getCoordinate() {
    return _patches.isEmpty ? null : _patches[0].getCoordinate();
  }

  @override
  List<Coordinate> getCoordinates() {
    return _patches.isEmpty
        ? []
        : _patches.map((g) => g.getCoordinates()).expand((i) => i).toList();
  }

  int getNumGeometries() {
    return numPatches;
  }

  /// Replies the <em>n</em>-th geometry in this collection.
  @specification(name: "getGeometryN")
  Geometry getGeometryN(int n) => elementAt(n);

  @override
  // TODO: implement isSimple
  bool get isSimple => null;
}

/// A TIN (triangulated irregular network) is a [PolyhedralSurface] consisting
/// only of [Triangle] patches.
class Tin extends PolyhedralSurface {
  Tin(List<Polygon> patches) : super(patches);

  /// Creates an empty tin.
  Tin.empty() : super(null);

  /// Creates a new Tin from the WKT string [wkt].
  ///
  /// Throws a [WKTError] if [wkt] isn't a valid representation of
  /// a [Tin].
  factory Tin.wkt(String wkt) {
    var g = parseWKT(wkt);
    if (g is! Tin) {
      throw WKTError("WKT string doesn't represent a Tin");
    }
    return g;
  }

  @override
  String get geometryType => "Tin";

  String get _wktName => "TIN";
}
