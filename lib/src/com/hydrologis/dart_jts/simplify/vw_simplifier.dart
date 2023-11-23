part of dart_jts;

///
/// Simplifies a {@link Geometry} using the Visvalingam-Whyatt area-based algorithm.
/// Ensures that any polygonal geometries returned are valid. Simple lines are not
/// guaranteed to remain simple after simplification. All geometry types are
/// handled. Empty and point geometries are returned unchanged. Empty geometry
/// components are deleted.
/// <p>
/// The simplification tolerance is specified as a distance.
/// This is converted to an area tolerance by squaring it.
/// <p>
/// Note that in general this algorithm does not preserve topology - e.g. polygons can be split,
/// collapse to lines or disappear holes can be created or disappear, and lines
/// can cross.
///
/// <h3>Known Bugs</h3>
/// <ul>
/// <li>Not yet optimized for performance
/// <li>Does not simplify the endpoint of rings
/// </ul>
/// <h3>To Do</h3>
/// <ul>
/// <li>Allow specifying desired number of vertices in the output
/// </ul>
///
/// @version 1.7
///
class VWSimplifier {
  Geometry inputGeom;
  double _distanceTolerance = 1.0;
  bool _isEnsureValidTopology = true;

  ///
  /// Creates a simplifier for a given geometry.
  ///
  /// @param inputGeom the geometry to simplify
  ///
  VWSimplifier(this.inputGeom);

  ///
  /// Simplifies a geometry using a given tolerance.
  ///
  /// @param geom geometry to simplify
  /// @param distanceTolerance the tolerance to use
  /// @return a simplified version of the geometry
  ///
  static Geometry simplify(Geometry geom, double distanceTolerance) {
    VWSimplifier simp = new VWSimplifier(geom);
    simp.setDistanceTolerance(distanceTolerance);
    return simp.getResultGeometry();
  }

  ///
  /// Sets the distance tolerance for the simplification.
  /// All vertices in the simplified geometry will be within this
  /// distance of the original geometry.
  /// The tolerance value must be non-negative.
  ///
  /// @param distanceTolerance the approximation tolerance to use
  ///
  void setDistanceTolerance(double distanceTolerance) {
    if (distanceTolerance < 0.0)
      throw ArgumentError("Tolerance must be non-negative");
    _distanceTolerance = distanceTolerance;
  }

  ///
  /// Controls whether simplified polygons will be "fixed" to have valid
  /// topology. The caller may choose to disable this because:
  /// <ul>
  /// <li>valid topology is not required
  /// <li>fixing topology is a relative expensive operation
  /// <li>in some pathological cases the topology fixing operation may either
  /// fail or run for too long
  /// </ul>
  ///
  /// The default is to fix polygon topology.
  ///
  /// @param isEnsureValidTopology
  ///
  void setEnsureValid(bool isEnsureValidTopology) {
    _isEnsureValidTopology = isEnsureValidTopology;
  }

  ///
  /// Gets the simplified geometry.
  ///
  /// @return the simplified geometry
  ///
  Geometry getResultGeometry() {
    // empty input produces an empty result
    if (inputGeom.isEmpty()) return inputGeom.copy();

    return (VWTransformer(_isEnsureValidTopology, _distanceTolerance))
        .transform(inputGeom);
  }
}

class VWTransformer extends GeometryTransformer {
  bool isEnsureValidTopology = true;
  double distanceTolerance;

  VWTransformer(this.isEnsureValidTopology, this.distanceTolerance);

  CoordinateSequence transformCoordinates(
      CoordinateSequence coords, Geometry parent) {
    List<Coordinate> inputPts = coords.toCoordinateArray();
    List<Coordinate> newPts;
    if (inputPts.length == 0) {
      newPts = List.empty(growable: true);
    } else {
      var simp = VWLineSimplifier(inputPts, distanceTolerance);
      newPts = simp.simplify();
    }
    return factory.getCoordinateSequenceFactory().create(newPts);
  }

  ///
  /// Simplifies a polygon, fixing it if required.
  ///
  Geometry transformPolygon(Polygon geom, Geometry parent) {
    if (geom.isEmpty()) {
      return geom;
    }
    Geometry rawGeom = super.transformPolygon(geom, parent);
    // don't try and correct if the parent is going to do this
    if (parent is MultiPolygon) {
      return rawGeom;
    }
    return createValidArea(rawGeom);
  }

  ///
  /// Simplifies a LinearRing. If the simplification results in a degenerate
  /// ring, remove the component.
  ///
  /// @return null if the simplification results in a degenerate ring
  ///
  Geometry transformLinearRing(LinearRing geom, Geometry parent) {
    bool removeDegenerateRings = parent is Polygon;
    Geometry simpResult = super.transformLinearRing(geom, parent);
    if (removeDegenerateRings && !(simpResult is LinearRing)) {
      return geom.geomFactory.createEmpty(geom.getDimension());
    }
    return simpResult;
  }

  ///
  /// Simplifies a MultiPolygon, fixing it if required.
  ///
  Geometry transformMultiPolygon(MultiPolygon geom, Geometry parent) {
    Geometry rawGeom = super.transformMultiPolygon(geom, parent);
    return createValidArea(rawGeom);
  }

  ///
  /// Creates a valid area geometry from one that possibly has bad topology
  /// (i.e. self-intersections). Since buffer can handle invalid topology, but
  /// always returns valid geometry, constructing a 0-width buffer "corrects"
  /// the topology. Note this only works for area geometries, since buffer
  /// always returns areas. This also may return empty geometries, if the input
  /// has no actual area.
  ///
  /// @param rawAreaGeom
  ///          an area geometry possibly containing self-intersections
  /// @return a valid area geometry
  ///
  Geometry createValidArea(Geometry rawAreaGeom) {
    // if geometry is invalid then make it valid
    if (isEnsureValidTopology && !rawAreaGeom.isValid())
      return rawAreaGeom.buffer(0.0);
    return rawAreaGeom;
  }
}
