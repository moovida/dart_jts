import "package:test/test.dart";
import 'package:dart_sfs/dart_sfs.dart';
import 'dart:math' as math;

void main() {
  group("WKTWriterTest - ", () {
    PrecisionModel precisionModel = PrecisionModel.fixedPrecision(1);
    GeometryFactory geometryFactory = GeometryFactory.withPrecisionModelSrid(precisionModel, 0);
    WKTWriter writer = WKTWriter();
    WKTWriter writer3D = WKTWriter.withDimension(3);
    WKTWriter writer2DM = WKTWriter.withDimension(3);

    writer2DM.setOutputOrdinates(OrdinateSet_XYM);

    test("testProperties", () {
      expect(OrdinateSet_XY, writer.getOutputOrdinates());
      expect(OrdinateSet_XYZ, writer3D.getOutputOrdinates());
      expect(OrdinateSet_XYM, writer2DM.getOutputOrdinates());

      GeometryFactory gf = GeometryFactory.withCoordinateSequenceFactory(PackedCoordinateSequenceFactory.DOUBLE_FACTORY);
      WKTWriter writer3DM = WKTWriter.withDimension(4);
      expect(OrdinateSet_XYZM, writer3DM.getOutputOrdinates());

      writer3DM.setOutputOrdinates(OrdinateSet_XY);
      expect(OrdinateSet_XY, writer3DM.getOutputOrdinates());
      writer3DM.setOutputOrdinates(OrdinateSet_XYZ);
      expect(OrdinateSet_XYZ, writer3DM.getOutputOrdinates());
      writer3DM.setOutputOrdinates(OrdinateSet_XYM);
      expect(OrdinateSet_XYM, writer2DM.getOutputOrdinates());
      writer3DM.setOutputOrdinates(OrdinateSet_XYZM);
      expect(OrdinateSet_XYZM, writer3DM.getOutputOrdinates());
    });

    test("testWritePoint", () {
      Point point = geometryFactory.createPoint(Coordinate(10, 10));
      expect("POINT (10 10)", writer.write(point).toString());
    });
    test("testWriteLineString", () {
      List<Coordinate> coordinates = [Coordinate(10, 10, 0), Coordinate(20, 20, 0), Coordinate(30, 40, 0)];
      LineString lineString = geometryFactory.createLineString(coordinates);
      String written = writer.write(lineString);
      expect("LINESTRING (10 10, 20 20, 30 40)", written);
    });
    test("testWritePolygon", () {
      List<Coordinate> coordinates = [Coordinate(10, 10, 0), Coordinate(10, 20, 0), Coordinate(20, 20, 0), Coordinate(20, 15, 0), Coordinate(10, 10, 0)];
      LinearRing linearRing = geometryFactory.createLinearRing(coordinates);
      Polygon polygon = geometryFactory.createPolygon(linearRing);
      expect("POLYGON ((10 10, 10 20, 20 20, 20 15, 10 10))", writer.write(polygon));
    });
    test("testWriteMultiPoint", () {
      List<Point> points = [geometryFactory.createPoint(Coordinate(10, 10, 0)), geometryFactory.createPoint(Coordinate(20, 20, 0))];
      MultiPoint multiPoint = geometryFactory.createMultiPoint(points);
      expect("MULTIPOINT ((10 10), (20 20))", writer.write(multiPoint));
    });
    test("testWriteMultiLineString", () {
      List<Coordinate> coordinates1 = [Coordinate(10, 10, 0), Coordinate(20, 20, 0)];
      LineString lineString1 = geometryFactory.createLineString(coordinates1);
      List<Coordinate> coordinates2 = [Coordinate(15, 15, 0), Coordinate(30, 15, 0)];
      LineString lineString2 = geometryFactory.createLineString(coordinates2);
      List<LineString> lineStrings = [lineString1, lineString2];
      MultiLineString multiLineString = geometryFactory.createMultiLineString(lineStrings);
      expect("MULTILINESTRING ((10 10, 20 20), (15 15, 30 15))", writer.write(multiLineString));
    });
    test("testWriteMultiPolygon", () {
      List<Coordinate> coordinates1 = [Coordinate(10, 10, 0), Coordinate(10, 20, 0), Coordinate(20, 20, 0), Coordinate(20, 15, 0), Coordinate(10, 10, 0)];
      LinearRing linearRing1 = geometryFactory.createLinearRing(coordinates1);
      Polygon polygon1 = geometryFactory.createPolygon(linearRing1);
      List<Coordinate> coordinates2 = [Coordinate(60, 60, 0), Coordinate(70, 70, 0), Coordinate(80, 60, 0), Coordinate(60, 60, 0)];
      LinearRing linearRing2 = geometryFactory.createLinearRing(coordinates2);
      Polygon polygon2 = geometryFactory.createPolygon(linearRing2);
      List<Polygon> polygons = [polygon1, polygon2];
      MultiPolygon multiPolygon = geometryFactory.createMultiPolygon(polygons);
//    System.out.println("MULTIPOLYGON (((10 10, 10 20, 20 20, 20 15, 10 10)), ((60 60, 70 70, 80 60, 60 60)))");
//    System.out.println(writer.write(multiPolygon).toString());
      expect("MULTIPOLYGON (((10 10, 10 20, 20 20, 20 15, 10 10)), ((60 60, 70 70, 80 60, 60 60)))", writer.write(multiPolygon));
    });
    test("testWriteGeometryCollection", () {
      Point point1 = geometryFactory.createPoint(Coordinate(10, 10));
      Point point2 = geometryFactory.createPoint(Coordinate(30, 30));
      List<Coordinate> coordinates = [Coordinate(15, 15, 0), Coordinate(20, 20, 0)];
      LineString lineString1 = geometryFactory.createLineString(coordinates);
      List<Geometry> geometries = [point1, point2, lineString1];
      GeometryCollection geometryCollection = geometryFactory.createGeometryCollection(geometries);
      expect("GEOMETRYCOLLECTION (POINT (10 10), POINT (30 30), LINESTRING (15 15, 20 20))", writer.write(geometryCollection));
    });
    test("testWriteLargeNumbers1", () {
      PrecisionModel precisionModel = PrecisionModel.fixedPrecision(1E9);
      GeometryFactory geometryFactory = GeometryFactory.withPrecisionModelSrid(precisionModel, 0);
      Point point1 = geometryFactory.createPoint(Coordinate(123456789012345680, 10E9));
      expect("POINT (123456789012345680 10000000000)", point1.toText());
    });
    test("testWrite3D", () {
      GeometryFactory geometryFactory = GeometryFactory.defaultPrecision();
      Point point = geometryFactory.createPoint(Coordinate(1, 1, 1));
      String wkt = writer3D.write(point);
      expect("POINT Z(1 1 1)", wkt);
      wkt = writer2DM.write(point);
      expect("POINT (1 1)", wkt);
    });
    test("testWrite3D_withNaN", () {
      GeometryFactory geometryFactory = GeometryFactory.defaultPrecision();
      List<Coordinate> coordinates = [Coordinate(1, 1), Coordinate(2, 2, 2)];
      LineString line = geometryFactory.createLineString(coordinates);
      String wkt = writer3D.write(line);
      expect("LINESTRING Z(1 1 NaN, 2 2 2)", wkt);
      wkt = writer2DM.write(line);
      expect("LINESTRING (1 1, 2 2)", wkt);
    });
  });

  WKTReader reader2D = getWKTReaderFromOrdinateSetAndScale(OrdinateSet_XY, 1.0);
  reader2D.setIsOldJtsCoordinateSyntaxAllowed(false);
  WKTReader reader2DOld = getWKTReaderFromOrdinateSetAndScale(OrdinateSet_XY, 1.0);
  reader2DOld.setIsOldJtsCoordinateSyntaxAllowed(true);
  WKTReader reader3D = getWKTReaderFromOrdinateSetAndScale(OrdinateSet_XYZ, 1.0);
  WKTReader reader2DM = getWKTReaderFromOrdinateSetAndScale(OrdinateSet_XYM, 1.0);
  WKTReader reader3DM = getWKTReaderFromOrdinateSetAndScale(OrdinateSet_XYZM, 1.0);
  group("WKTReaderTest - ", () {
    test("", () {
      // arrange
      CoordinateSequence seq = createSequence(OrdinateSet_XYZ, [10, 10]);
      seq.setOrdinate(0, CoordinateSequence.Z, double.nan);

      // act
      Point pt1 = reader2DOld.read("POINT (10 10 NaN)");
      Point pt2 = reader2DOld.read("POINT (10 10 nan)");
      Point pt3 = reader2DOld.read("POINT (10 10 NAN)");

      // assert
      expect(checkEqual(seq, pt1.getCoordinateSequence()), true);
      expect(checkEqual(seq, pt2.getCoordinateSequence()), true);
      expect(checkEqual(seq, pt3.getCoordinateSequence()), true);
    });
    test("", () {

    });
    test("", () {

    });
    test("", () {

    });
    test("", () {

    });
    test("", () {

    });
    test("", () {

    });
  });
}

CoordinateSequence createSequence(List<Ordinate> ordinateFlags, List<double> xy) {
// get the number of dimension to verify size of provided ordinate values array
  int dimension = requiredDimension(ordinateFlags);

// inject additional values
  List<double> ordinateValues = injectZM(ordinateFlags, xy);

  if ((ordinateValues.length % dimension) != 0) throw new ArgumentError("ordinateFlags and number of provided ordinate values don't match");

// get the required size of the sequence
  int size = ordinateValues.length ~/ dimension;

// create a sequence capable of storing all ordinate values.
  CoordinateSequence res = getCSFactory(ordinateFlags).createSizeDim(size, requiredDimension(ordinateFlags));

// fill in values
  int k = 0;
  for (int i = 0; i < ordinateValues.length; i += dimension) {
    for (int j = 0; j < dimension; j++) res.setOrdinate(k, j, ordinateValues[i + j]);
    k++;
  }

  return res;
}

int requiredDimension(List<Ordinate> ordinateFlags) {
  return ordinateFlags.length;
}

List<double> injectZM(List<Ordinate> ordinateFlags, List<double> xy) {
  int size = xy.length ~/ 2;
  int dimension = requiredDimension(ordinateFlags);
  List<double> res = List(size * dimension);
  int k = 0;
  for (int i = 0; i < xy.length; i += 2) {
    res[k++] = xy[i];
    res[k++] = xy[i + 1];
    if (ordinateFlags.contains(Ordinate.Z)) res[k++] = 10;
    if (ordinateFlags.contains(Ordinate.M)) res[k++] = 11;
  }
  return res;
}

/// Gets a {@link CoordinateSequenceFactory} that can create sequences
/// for ordinates defined in the provided bit-pattern.
/// @param ordinateFlags a bit-pattern of ordinates
/// @return a {@code CoordinateSequenceFactory}
CoordinateSequenceFactory getCSFactory(List<Ordinate> ordinateFlags) {
  if (ordinateFlags.contains(Ordinate.M)) {
    return PackedCoordinateSequenceFactory.DOUBLE_FACTORY;
  }

  return CoordinateArraySequenceFactory();
}

/// Gets a {@link WKTReader} to read geometries from WKT with expected ordinates.
///
/// @param ordinateFlags a set of expected ordinates
/// @param precisionModel a precision model
///
/// @return a {@code WKTReader}
WKTReader getWKTReader(List<Ordinate> ordinateFlags, PrecisionModel precisionModel) {
  WKTReader result;

  if (!ordinateFlags.contains(Ordinate.X)) ordinateFlags.add(Ordinate.X);
  if (!ordinateFlags.contains(Ordinate.Y)) ordinateFlags.add(Ordinate.Y);

  if (ordinateFlags.length == 2) {
    result = WKTReader.withFactory(GeometryFactory(precisionModel, 0, CoordinateArraySequenceFactory()));
    result.setIsOldJtsCoordinateSyntaxAllowed(false);
  } else if (ordinateFlags.contains(Ordinate.Z)) {
    result = WKTReader.withFactory(GeometryFactory(precisionModel, 0, CoordinateArraySequenceFactory()));
  } else if (ordinateFlags.contains(Ordinate.M)) {
    result = WKTReader.withFactory(GeometryFactory(precisionModel, 0, PackedCoordinateSequenceFactory.DOUBLE_FACTORY));
    result.setIsOldJtsCoordinateSyntaxAllowed(false);
  } else {
    result = WKTReader.withFactory(GeometryFactory(precisionModel, 0, PackedCoordinateSequenceFactory.DOUBLE_FACTORY));
  }

  return result;
}

/// Gets a {@link WKTReader} to read geometries from WKT with expected ordinates.
///
/// @param ordinateFlags a set of expected ordinates
/// @return a {@code WKTReader}
WKTReader getWKTReaderFromOrdinateSet(List<Ordinate> ordinateFlags) {
  return getWKTReader(ordinateFlags, PrecisionModel());
}

/// Gets a {@link WKTReader} to read geometries from WKT with expected ordinates.
///
/// @param ordinateFlags a set of expected ordinates
/// @param scale         a scale value to create a {@link PrecisionModel}
///
/// @return a {@code WKTReader}
WKTReader getWKTReaderFromOrdinateSetAndScale(List<Ordinate> ordinateFlags, double scale) {
  return getWKTReader(ordinateFlags, PrecisionModel.fixedPrecision(scale));
}

/// Checks two {@link CoordinateSequence}s for equality. The following items are checked:
/// <ul>
///   <li>size</li><li>dimension</li><li>ordinate values</li>
/// </ul>
/// @param seq1 a sequence
/// @param seq2 another sequence
/// @return {@code true} if both sequences are equal
bool checkEqual(CoordinateSequence seq1, CoordinateSequence seq2) {
  return checkEqualWithTolerance(seq1, seq2, 0.0);
}

/// Checks two {@link CoordinateSequence}s for equality. The following items are checked:
/// <ul>
///   <li>size</li><li>dimension</li><li>ordinate values with {@code tolerance}</li>
/// </ul>
/// @param seq1 a sequence
/// @param seq2 another sequence
/// @return {@code true} if both sequences are equal
bool checkEqualWithTolerance(CoordinateSequence seq1, CoordinateSequence seq2, double tolerance) {
  if (seq1.getDimension() != seq2.getDimension()) return false;
  return checkEqualWithDimTol(seq1, seq2, seq1.getDimension(), tolerance);
}

/// Checks two {@link CoordinateSequence}s for equality. The following items are checked:
/// <ul>
///   <li>size</li><li>dimension up to {@code dimension}</li><li>ordinate values</li>
/// </ul>
/// @param seq1 a sequence
/// @param seq2 another sequence
/// @return {@code true} if both sequences are equal
bool checkEqualWithDim(CoordinateSequence seq1, CoordinateSequence seq2, int dimension) {
  return checkEqualWithDimTol(seq1, seq2, dimension, 0.0);
}

/// Checks two {@link CoordinateSequence}s for equality. The following items are checked:
/// <ul>
///   <li>size</li><li>dimension up to {@code dimension}</li><li>ordinate values with {@code tolerance}</li>
/// </ul>
/// @param seq1 a sequence
/// @param seq2 another sequence
/// @return {@code true} if both sequences are equal
bool checkEqualWithDimTol(CoordinateSequence seq1, CoordinateSequence seq2, int dimension, double tolerance) {
  if (seq1 != null && seq2 == null) return false;
  if (seq1 == null && seq2 != null) return false;

  if (seq1.size() != seq2.size()) return false;

  if (seq1.getDimension() < dimension) throw ArgumentError("dimension too high for seq1");
  if (seq2.getDimension() < dimension) throw ArgumentError("dimension too high for seq2");

  for (int i = 0; i < seq1.size(); i++) {
    for (int j = 0; j < dimension; j++) {
      double val1 = seq1.getOrdinate(i, j);
      double val2 = seq2.getOrdinate(i, j);
      if (val1.isNaN) {
        if (!val2.isNaN) return false;
      } else if ((val1 - val2).abs() > tolerance) return false;
    }
  }

  return true;
}
