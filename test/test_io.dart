import "package:test/test.dart";
import 'package:dart_sfs/dart_sfs.dart';
import 'dart:math' as math;

void main() {
//  group("WKTReaderTest - ", () {
//    test("testRead", () {
//      var seq = CoordinateArraySequenceFactory().create([Coordinate(10, 10)]);
//
//      WKTReader reader = WKTReader();
//
//      // act
//      Point pt1 = reader.read("POINT (10 10)");
//
//      // assert
//      expect(checkEqual(seq, CoordinateArraySequenceFactory().create(pt1.getCoordinates())), true);
//    });
//  });

  PrecisionModel precisionModel = PrecisionModel.fixedPrecision(1);
  GeometryFactory geometryFactory = GeometryFactory.withPrecisionModelSrid(precisionModel, 0);
  WKTWriter writer = WKTWriter();
  WKTWriter writer3D = WKTWriter.withDimension(3);
  WKTWriter writer2DM = WKTWriter.withDimension(3);

  group("WKTWriterTest - ", () {
    writer2DM.setOutputOrdinates(XYM);

    test("testProperties", () {
      expect(XY, writer.getOutputOrdinates());
      expect(XYZ, writer3D.getOutputOrdinates());
      expect(XYM, writer2DM.getOutputOrdinates());

      GeometryFactory gf = GeometryFactory.withCoordinateSequenceFactory(PackedCoordinateSequenceFactory.DOUBLE_FACTORY);
      WKTWriter writer3DM = WKTWriter.withDimension(4);
      expect(XYZM, writer3DM.getOutputOrdinates());

      writer3DM.setOutputOrdinates(XY);
      expect(XY, writer3DM.getOutputOrdinates());
      writer3DM.setOutputOrdinates(XYZ);
      expect(XYZ, writer3DM.getOutputOrdinates());
      writer3DM.setOutputOrdinates(XYM);
      expect(XYM, writer2DM.getOutputOrdinates());
      writer3DM.setOutputOrdinates(XYZM);
      expect(XYZM, writer3DM.getOutputOrdinates());
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

//  group("StringUtils - ", () {
//    test("", () {
//    });
//  });
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
