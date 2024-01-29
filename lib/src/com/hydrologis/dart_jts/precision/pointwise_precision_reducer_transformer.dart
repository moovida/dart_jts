part of dart_jts;
/**
 * A transformer to reduce the precision of a geometry pointwise.
 *
 * @author mdavis
 *
 */
class PointwisePrecisionReducerTransformer extends GeometryTransformer {

  static Geometry reduce(Geometry geom, PrecisionModel targetPM) {
    PointwisePrecisionReducerTransformer trans = new PointwisePrecisionReducerTransformer(targetPM);
    return trans.transform(geom);
  }

  PrecisionModel targetPM;

  PointwisePrecisionReducerTransformer(this.targetPM);

  CoordinateSequence transformCoordinates(
      CoordinateSequence coordinates, Geometry parent) {
    List<Coordinate> coordsReduce = reducePointwise(coordinates);
    return factory.getCoordinateSequenceFactory().create(coordsReduce);
  }

  List<Coordinate> reducePointwise(CoordinateSequence coordinates) {
    List<Coordinate> coordReduce = List.empty(growable:true);
    // copy coordinates and reduce
    for (int i = 0; i < coordinates.size(); i++) {
      Coordinate coord = coordinates.getCoordinate(i).copy();
      targetPM.makeCoordinatePrecise(coord);
      coordReduce.add(coord);
    }
    return coordReduce;
  }

}