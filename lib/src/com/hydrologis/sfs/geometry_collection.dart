part of dart_sfs;

/// the singleton empty geometry collection
final _EMPTY_GEOMETRY_COLLECTION = GeometryCollection(null);

/// A [GeometryCollection] is a geometric object that is a collection of
/// some number of geometric objects.
///
/// It implements the accesor methods [numGeometries] and [getGeometryN]
/// which are specified in the SFS. In addition, it provides the more
/// Dart'ish [length] property and an overloaded index operator. It also
/// implements the [Iterable] interface.
///
class GeometryCollection extends Geometry
    with IterableMixin<Geometry>, _GeometryContainerMixin {
  List<Geometry> _geometries;

  /// Creates a geometry collection given a collection of
  /// [geometries].
  GeometryCollection(Iterable<Geometry> geometries) {
    if (geometries == null || geometries.isEmpty) {
      _geometries = null;
    } else {
      _geometries = List<Geometry>.from(geometries, growable: false);
      _require(this.every((p) => p != null));
    }
  }

  /// Creates an empty geometry collection.
  factory GeometryCollection.empty() => _EMPTY_GEOMETRY_COLLECTION;

  /// Creates a new geometry collection from the WKT string [wkt].
  ///
  /// Throws a [WKTError] if [wkt] isn't a valid representation of
  /// a [GeometryCollection].
  factory GeometryCollection.wkt(String wkt) {
    var g = parseWKT(wkt);
    if (g is! GeometryCollection) {
      throw WKTError("WKT string doesn't represent a GeometryCollection");
    }
    return g;
  }

  @override
  Coordinate getCoordinate() {
    return _geometries.isEmpty ? null : _geometries[0].getCoordinate();
  }

  @override
  List<Coordinate> getCoordinates() {
    return _geometries.isEmpty
        ? []
        : _geometries.map((g) => g.getCoordinates()).expand((i) => i).toList();
  }

  int getNumGeometries() {
    return length;
  }

  /// Replies the <em>n</em>-th geometry in this collection.
  @specification(name: "getGeometryN")
  Geometry getGeometryN(int n) => elementAt(n);

  /// Replies the <em>n</em>-th geometry in this collection.
  ///
  /// This is the Dart'ish implemenation of `getGeometryN()` using
  /// operator overloading.
  @specification(name: "getGeometryN")
  operator [](int n) => elementAt(n);

  /// the iterator to access the geometry objects
  Iterator<Geometry> get iterator {
    if (_geometries == null) {
      return [].iterator;
    } else {
      return _geometries.iterator;
    }
  }

  /// A geometry collection is simple if all its child geometries are
  /// simple.
  @override
  bool get isSimple => every((g) => g.isSimple);

  @override
  _writeTaggedWKT(writer, {bool withZ: false, bool withM: false}) {
    writer.write("GEOMETRYCOLLECTION");
    writer.blank();
    if (!isEmpty) {
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
        elementAt(i)._writeTaggedWKT(writer, withZ: withZ, withM: withM);
      }
      writer..newline();
      writer
        ..decIdent()
        ..ident()
        ..rparen();
    }
  }

  @override
  String get geometryType => "GeometryCollection";

  @override
  int get dimension => fold(0, (prev, g) => math.max(prev, g.dimension));

  @override
  Geometry get boundary => throw UnsupportedError(
      "Operation does not support GeometryCollection arguments");
}

mixin _GeometryContainerMixin<E extends Geometry> on Iterable<E> {
  bool _is3D;

  _computeIs3D() {
    if (this.isEmpty) {
      _is3D = false;
    } else {
      _is3D = every((g) => g.is3D);
    }
  }

  /// A collection of geometries is considered 3D if *every* child geometry
  /// has a non-null z-component.
  ///
  /// The value of this property is computed upon first access and then
  /// cached. Subsequent reads of the property efficiently reply the cached
  /// value.
  @override
  bool get is3D {
    if (_is3D == null) _computeIs3D();
    return _is3D;
  }

  bool _isMeasured;

  _computeIsMeasured() {
    if (this.isEmpty) {
      _isMeasured = false;
    } else {
      _isMeasured = every((g) => g.isMeasured);
    }
  }

  /// A collection of geometries is considered *measured* if *every* child
  /// geometry has an m-component.
  ///
  /// The value of this property is computed upon first access and then
  /// cached. Subsequent reads of the property efficiently reply the cached
  /// value.
  @override
  bool get isMeasured {
    if (_isMeasured == null) _computeIsMeasured();
    return _isMeasured;
  }

  Envelope _computeEnvelope() {
    if (this.isEmpty) return Envelope.empty();
    Envelope e = Envelope.empty();
    forEach((p) => e.growTo(p));
    return e;
  }

  operator [](int n) => this.elementAt(n);
}
