part of dart_sfs;

/// A [Surface] is a 2-dimensional geometric object.
///
abstract class Surface extends Geometry {
  @override
  get dimension => 2;

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
