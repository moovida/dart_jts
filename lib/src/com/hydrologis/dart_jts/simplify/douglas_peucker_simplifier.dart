part of dart_jts;

///
/// Simplifies a {@link Geometry} using the Douglas-Peucker algorithm.
/// Ensures that any polygonal geometries returned are valid.
/// Simple lines are not guaranteed to remain simple after simplification.
/// All geometry types are handled.
/// Empty and point geometries are returned unchanged.
/// Empty geometry components are deleted.
/// <p>
/// Note that in general D-P does not preserve topology -
/// e.g. polygons can be split, collapse to lines or disappear
/// holes can be created or disappear,
/// and lines can cross.
/// To simplify geometry while preserving topology use {@link TopologyPreservingSimplifier}.
/// (However, using D-P is significantly faster).
///<h2>KNOWN BUGS</h2>
///<ul>
///<li>In some cases the approach used to clean invalid simplified polygons
///can distort the output geometry severely.
///</ul>
///

class DouglasPeuckerSimplifier {
  Geometry inputGeom;
  double _distanceTolerance = 1;
  bool _isEnsureValidTopology = true;

  DouglasPeuckerSimplifier(this.inputGeom);

  ///
  /// Simplifies a geometry using a given tolerance.
  ///
  /// @param geom geometry to simplify
  /// @param distanceTolerance the tolerance to use
  /// @return a simplified version of the geometry
  ///

  static Geometry? simplify(Geometry geom, double tolerance) {
    DouglasPeuckerSimplifier tss = DouglasPeuckerSimplifier(geom);
    tss._distanceTolerance = tolerance;
    return tss.getResultGeometry();
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
  /// Controls whether simplified polygons will be "fixed"
  /// to have valid topology.
  /// The caller may choose to disable this because:
  /// <ul>
  /// <li>valid topology is not required
  /// <li>fixing topology is a relative expensive operation
  /// <li>in some pathological cases the topology fixing operation may either fail or run for too long
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
  Geometry? getResultGeometry() {
    // empty input produces an empty result
    if (inputGeom.isEmpty()) return inputGeom.copy();
    return (DPTransformer(_isEnsureValidTopology, _distanceTolerance))
        .transform(inputGeom);
  }
}

class DPTransformer extends GeometryTransformer {
  bool isEnsureValidTopology = true;
  double distanceTolerance;

  DPTransformer(this.isEnsureValidTopology, this.distanceTolerance);
  CoordinateSequence transformCoordinates(
      CoordinateSequence coords, Geometry parent) {
    List<Coordinate> inputPts = coords.toCoordinateArray();
    List<Coordinate> newPts = List.empty(growable: true);
    if (inputPts.length > 0) {
      DouglasPeuckerLineSimplifier simplifier =
          DouglasPeuckerLineSimplifier(inputPts);
      simplifier.setDistanceTolerance(distanceTolerance);
      newPts = simplifier.simplify();
    }
    return factory.getCoordinateSequenceFactory().create(newPts);
  }

  ///
  /// Simplifies a polygon, fixing it if required.
  ///
  Geometry transformPolygon(Polygon geom, Geometry parent) {
    // empty geometries are simply removed
    if (geom.isEmpty()) return geom;
    Geometry rawGeom = super.transformPolygon(geom, parent);
    // don't try and correct if the parent is going to do this
    if (parent is MultiPolygon) {
      return rawGeom;
    }
    return createValidArea(rawGeom);
  }

  ///
  /// Simplifies a LinearRing.  If the simplification results
  /// in a degenerate ring, remove the component.
  ///
  /// @return null if the simplification results in a degenerate ring
  ////
  //*
  Geometry transformLinearRing(LinearRing geom, Geometry parent) {
    bool removeDegenerateRings = parent is Polygon;
    Geometry? simpResult = super.transformLinearRing(geom, parent);
    if (removeDegenerateRings && !(simpResult is LinearRing)) {
      return geom.geomFactory.createEmpty(geom.getDimension());
    }
    return simpResult;
  }

  ///

  ///
  /// Simplifies a MultiPolygon, fixing it if required.
  ///
  Geometry transformMultiPolygon(MultiPolygon geom, Geometry parent) {
    Geometry rawGeom = super.transformMultiPolygon(geom, parent);
    return createValidArea(rawGeom);
  }

  ///
  /// Creates a valid area geometry from one that possibly has
  /// bad topology (i.e. self-intersections).
  /// Since buffer can handle invalid topology, but always returns
  /// valid geometry, constructing a 0-width buffer "corrects" the
  /// topology.
  /// Note this only works for area geometries, since buffer always returns
  /// areas.  This also may return empty geometries, if the input
  /// has no actual area.
  /// If the input is empty or is not polygonal,
  /// this ensures that POLYGON EMPTY is returned.
  ///
  /// @param rawAreaGeom an area geometry possibly containing self-intersections
  /// @return a valid area geometry
  ///
  Geometry createValidArea(Geometry rawAreaGeom) {
    bool isValidArea = rawAreaGeom.getDimension() == 2 && rawAreaGeom.isValid();
    // if geometry is invalid then make it valid
    if (isEnsureValidTopology && !isValidArea) {
      return rawAreaGeom.buffer(0.0);
    }
    return rawAreaGeom;
  }
}
