part of dart_jts;

////**
/// A transformer to reduce the precision of geometry in a
/// topologically valid way.
/// Repeated points are removed.
/// If geometry elements collapse below their valid length,
/// they may be removed
/// by specifying <code>isRemoveCollapsed</code as <code>true</code>.
///
/// @author mdavis
///
///
class PrecisionReducerTransformer extends GeometryTransformer {
  static Geometry reduce(
      Geometry geom, PrecisionModel targetPM, bool isRemoveCollapsed) {
    PrecisionReducerTransformer trans =
        new PrecisionReducerTransformer(targetPM, isRemoveCollapsed);
    return trans.transform(geom);
  }

  PrecisionModel targetPM;
  bool isRemoveCollapsed = false;

  PrecisionReducerTransformer(this.targetPM, this.isRemoveCollapsed);

  CoordinateSequence transformCoordinates(
      CoordinateSequence coordinates, Geometry parent) {
    List<Coordinate> coordsReduce = reduceCompress(coordinates);

    ///
    /// Check if the removal of repeated points collapsed the coordinate
    /// list to an invalid size for the type of the parent geometry. It is not
    /// necessary to check for Point collapses, since the coordinate list can
    /// never collapse to less than one point. If the size is invalid, return
    /// the full-size coordinate array first computed, or null if collapses are
    /// being removed. (This may create an invalid geometry - the client must
    /// handle this.)
    ///
    int minSize = 0;
    if (parent is LineString) minSize = 2;
    if (parent is LinearRing) minSize = LinearRing.MINIMUM_VALID_SIZE;

    if (coordsReduce.length < minSize) {
      if (isRemoveCollapsed) {
        return factory.getCoordinateSequenceFactory().create(List.empty());
      }
      coordsReduce = extend(coordsReduce, minSize);
    }
    return factory.getCoordinateSequenceFactory().create(coordsReduce);
  }

  List<Coordinate> extend(List<Coordinate> coords, int minLength) {
    if (coords.length >= minLength) return coords;
    List<Coordinate> exCoords = List.filled(minLength,Coordinate.empty2D());
    for (int i = 0; i < exCoords.length; i++) {
      int iSrc = i < coords.length ? i : coords.length - 1;
      exCoords[i] = coords[iSrc].copy();
    }
    return exCoords;
  }

  List<Coordinate> reduceCompress(CoordinateSequence coordinates) {
    CoordinateList noRepeatCoordList = new CoordinateList();
    // copy coordinates and reduce
    for (int i = 0; i < coordinates.size(); i++) {
      Coordinate coord = coordinates.getCoordinate(i).copy();
      targetPM.makeCoordinatePrecise(coord);
      noRepeatCoordList.addCoord(coord, false);
    }
    // remove repeated points, to simplify geometry as much as possible
    List<Coordinate> noRepeatCoords = noRepeatCoordList.toCoordinateArray(true);
    return noRepeatCoords;
  }

  Geometry transformPolygon(Polygon geom, Geometry parent) {
    return reduceArea(geom);
  }

  Geometry transformMultiPolygon(MultiPolygon geom, Geometry parent) {
    return reduceArea(geom);
  }

  Geometry reduceArea(Geometry geom) {
    Geometry? reduced = PrecisionReducer.reducePrecision(geom, targetPM);
    return reduced!;
  }
}
