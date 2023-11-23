part of dart_jts;

class TopologyPreservingSimplifier {
  ///
  /// Simplifies a geometry and ensures that
  /// the result is a valid geometry having the
  /// same dimension and number of components as the input,
  /// and with the components having the same topological
  /// relationship.
  /// <p>
  /// If the input is a polygonal geometry
  /// ( {@link Polygon} or {@link MultiPolygon} ):
  /// <ul>
  /// <li>The result has the same number of shells and holes as the input,
  /// with the same topological structure
  /// <li>The result rings touch at <b>no more</b> than the number of touching points in the input
  /// (although they may touch at fewer points).
  /// The key implication of this statement is that if the
  /// input is topologically valid, so is the simplified output.
  /// </ul>
  /// For linear geometries, if the input does not contain
  /// any intersecting line segments, this property
  /// will be preserved in the output.
  /// <p>
  /// For all geometry types, the result will contain
  /// enough vertices to ensure validity.  For polygons
  /// and closed linear geometries, the result will have at
  /// least 4 vertices; for open linestrings the result
  /// will have at least 2 vertices.
  /// <p>
  /// All geometry types are handled.
  /// Empty and point geometries are returned unchanged.
  /// Empty geometry components are deleted.
  /// <p>
  /// The simplification uses a maximum-distance difference algorithm
  /// similar to the Douglas-Peucker algorithm.
  ///
  /// <h3>KNOWN BUGS</h3>
  /// <ul>
  /// <li>May create invalid topology if there are components which are
  /// small relative to the tolerance value.
  /// In particular, if a small hole is very near an edge, it is possible for the edge to be moved by
  /// a relatively large tolerance value and end up with the hole outside the result shell
  /// (or inside another hole).
  /// Similarly, it is possible for a small polygon component to end up inside
  /// a nearby larger polygon.
  /// A workaround is to test for this situation in post-processing and remove
  /// any invalid holes or polygons.
  /// </ul>
  ///
  /// @author Martin Davis
  /// @see DouglasPeuckerSimplifier
  ///
  ///
  Geometry inputGeom;
  TaggedLinesSimplifier lineSimplifier = TaggedLinesSimplifier();
  Map<Geometry, TaggedLineString> linestringMap = {};
  TopologyPreservingSimplifier(this.inputGeom);

  ///
  /// Sets the distance tolerance for the simplification.
  /// All vertices in the simplified geometry will be within this
  /// distance of the original geometry.
  /// The tolerance value must be non-negative.  A tolerance value
  /// of zero is effectively a no-op.
  ///
  /// @param distanceTolerance the approximation tolerance to use
  ///
  void setDistanceTolerance(double distanceTolerance) {
    if (distanceTolerance < 0.0)
      throw new ArgumentError("Tolerance must be non-negative");
    lineSimplifier.setDistanceTolerance(distanceTolerance);
  }

  static Geometry? simplify(Geometry geom, double distanceTolerance) {
    TopologyPreservingSimplifier tss = new TopologyPreservingSimplifier(geom);
    tss.setDistanceTolerance(distanceTolerance);
    return tss.getResultGeometry();
  }

  Geometry? getResultGeometry() {
    // empty input produces an empty result
    if (inputGeom.isEmpty()) return inputGeom.copy();

    linestringMap = new HashMap();
    inputGeom.applyGCF(LineStringMapBuilderFilter(this));
    lineSimplifier.simplify(linestringMap.values.toList());
    Geometry? result =
        LineStringTransformer(linestringMap).transform(inputGeom);
    return result;
  }
}

///
/// A filter to add linear geometries to the linestring map
/// with the appropriate minimum size constraint.
/// Closed {@link LineString}s (including {@link LinearRing}s
/// have a minimum output size constraint of 4,
/// to ensure the output is valid.
/// For all other linestrings, the minimum size is 2 points.
///
/// @author Martin Davis
///
///
class LineStringMapBuilderFilter extends GeometryComponentFilter {
  TopologyPreservingSimplifier tps;

  LineStringMapBuilderFilter(this.tps);

  ///
  /// Filters linear geometries.
  ///
  /// geom a geometry of any type
  ///
  void filter(Geometry geom) {
    if (geom is LineString) {
      LineString line = geom;
// skip empty geometries
      if (line.isEmpty()) return;

      int minSize = line.isClosed() ? 4 : 2;
      TaggedLineString taggedLine = new TaggedLineString(line, minSize);
      tps.linestringMap[line] = taggedLine;
    }
  }
}

class LineStringTransformer extends GeometryTransformer {
  Map<Geometry, TaggedLineString> linestringMap;

  LineStringTransformer(this.linestringMap);

  CoordinateSequence transformCoordinates(
      CoordinateSequence coords, Geometry parent) {
    // for linear components (including rings), simplify the linestring
    if (parent is LineString) {
      TaggedLineString? taggedLine = linestringMap[parent];
      if (taggedLine != null) {
        return createCoordinateSequence(taggedLine.getResultCoordinates());
      }
    }
    // for anything else (e.g. points) just copy the coordinates
    return super.transformCoordinates(coords, parent);
  }
}
