import "package:test/test.dart";
import 'package:dart_jts/dart_jts.dart';
import 'dart:math' as math;
import 'dart:typed_data';
import 'testing_utilities.dart';

void main() {
  group("WKBTest - ", () {
    WKBTest test = WKBTest();


  });


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
      List<Coordinate> coordinates = [Coordinate.fromXYZ(10, 10, 0), Coordinate.fromXYZ(20, 20, 0), Coordinate.fromXYZ(30, 40, 0)];
      LineString lineString = geometryFactory.createLineString(coordinates);
      String written = writer.write(lineString);
      expect("LINESTRING (10 10, 20 20, 30 40)", written);
    });
    test("testWritePolygon", () {
      List<Coordinate> coordinates = [
        Coordinate.fromXYZ(10, 10, 0),
        Coordinate.fromXYZ(10, 20, 0),
        Coordinate.fromXYZ(20, 20, 0),
        Coordinate.fromXYZ(20, 15, 0),
        Coordinate.fromXYZ(10, 10, 0)
      ];
      LinearRing linearRing = geometryFactory.createLinearRing(coordinates);
      Polygon polygon = geometryFactory.createPolygonFromRing(linearRing);
      expect("POLYGON ((10 10, 10 20, 20 20, 20 15, 10 10))", writer.write(polygon));
    });
    test("testWriteMultiPoint", () {
      List<Point> points = [geometryFactory.createPoint(Coordinate.fromXYZ(10, 10, 0)), geometryFactory.createPoint(Coordinate.fromXYZ(20, 20, 0))];
      MultiPoint multiPoint = geometryFactory.createMultiPoint(points);
      expect("MULTIPOINT ((10 10), (20 20))", writer.write(multiPoint));
    });
    test("testWriteMultiLineString", () {
      List<Coordinate> coordinates1 = [Coordinate.fromXYZ(10, 10, 0), Coordinate.fromXYZ(20, 20, 0)];
      LineString lineString1 = geometryFactory.createLineString(coordinates1);
      List<Coordinate> coordinates2 = [Coordinate.fromXYZ(15, 15, 0), Coordinate.fromXYZ(30, 15, 0)];
      LineString lineString2 = geometryFactory.createLineString(coordinates2);
      List<LineString> lineStrings = [lineString1, lineString2];
      MultiLineString multiLineString = geometryFactory.createMultiLineString(lineStrings);
      expect("MULTILINESTRING ((10 10, 20 20), (15 15, 30 15))", writer.write(multiLineString));
    });
    test("testWriteMultiPolygon", () {
      List<Coordinate> coordinates1 = [
        Coordinate.fromXYZ(10, 10, 0),
        Coordinate.fromXYZ(10, 20, 0),
        Coordinate.fromXYZ(20, 20, 0),
        Coordinate.fromXYZ(20, 15, 0),
        Coordinate.fromXYZ(10, 10, 0)
      ];
      LinearRing linearRing1 = geometryFactory.createLinearRing(coordinates1);
      Polygon polygon1 = geometryFactory.createPolygonFromRing(linearRing1);
      List<Coordinate> coordinates2 = [
        Coordinate.fromXYZ(60, 60, 0),
        Coordinate.fromXYZ(70, 70, 0),
        Coordinate.fromXYZ(80, 60, 0),
        Coordinate.fromXYZ(60, 60, 0)
      ];
      LinearRing linearRing2 = geometryFactory.createLinearRing(coordinates2);
      Polygon polygon2 = geometryFactory.createPolygonFromRing(linearRing2);
      List<Polygon> polygons = [polygon1, polygon2];
      MultiPolygon multiPolygon = geometryFactory.createMultiPolygon(polygons);
//    System.out.println("MULTIPOLYGON (((10 10, 10 20, 20 20, 20 15, 10 10)), ((60 60, 70 70, 80 60, 60 60)))");
//    System.out.println(writer.write(multiPolygon).toString());
      expect("MULTIPOLYGON (((10 10, 10 20, 20 20, 20 15, 10 10)), ((60 60, 70 70, 80 60, 60 60)))", writer.write(multiPolygon));
    });
    test("testWriteGeometryCollection", () {
      Point point1 = geometryFactory.createPoint(Coordinate(10, 10));
      Point point2 = geometryFactory.createPoint(Coordinate(30, 30));
      List<Coordinate> coordinates = [Coordinate.fromXYZ(15, 15, 0), Coordinate.fromXYZ(20, 20, 0)];
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
      Point point = geometryFactory.createPoint(Coordinate.fromXYZ(1, 1, 1));
      String wkt = writer3D.write(point);
      expect("POINT Z(1 1 1)", wkt);
      wkt = writer2DM.write(point);
      expect("POINT (1 1)", wkt);
    });
    test("testWrite3D_withNaN", () {
      GeometryFactory geometryFactory = GeometryFactory.defaultPrecision();
      List<Coordinate> coordinates = [Coordinate(1, 1), Coordinate.fromXYZ(2, 2, 2)];
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
    test("testReadNaN", () {
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
    test("testReadPoint", () {
      // arrange
      List<double> coordinates = [10, 10];
      CoordinateSequence seqPt2D = createSequence(OrdinateSet_XY, coordinates);
      CoordinateSequence seqPt2DE = createSequence(OrdinateSet_XY, <double>[]);
      CoordinateSequence seqPt3D = createSequence(OrdinateSet_XYZ, coordinates);
      CoordinateSequence seqPt2DM = createSequence(OrdinateSet_XYM, coordinates);
      CoordinateSequence seqPt3DM = createSequence(OrdinateSet_XYZM, coordinates);

      // act
      Point pt2D = reader2D.read("POINT (10 10)");
      Point pt2DE = reader2D.read("POINT EMPTY");
      Point pt3D = reader3D.read("POINT Z(10 10 10)");
      Point pt2DM = reader2DM.read("POINT M(10 10 11)");
      Point pt3DM = reader3DM.read("POINT ZM(10 10 10 11)");

      // assert
      var pt2DCS = pt2D.getCoordinateSequence();
      expect(checkEqual(seqPt2D, pt2DCS), true);
      var pt2DECS = pt2DE.getCoordinateSequence();
      expect(checkEqual(seqPt2DE, pt2DECS), true);
      expect(checkEqual(seqPt3D, pt3D.getCoordinateSequence()), true);
      expect(checkEqual(seqPt2DM, pt2DM.getCoordinateSequence()), true);
      expect(checkEqual(seqPt3DM, pt3DM.getCoordinateSequence()), true);
    });
    test("testReadLineString", () {
      // arrange
      List<double> coordinates = [10, 10, 20, 20, 30, 40];
      CoordinateSequence seqLs2D = createSequence(OrdinateSet_XY, coordinates);
      CoordinateSequence seqLs2DE = createSequence(OrdinateSet_XY, <double>[]);
      CoordinateSequence seqLs3D = createSequence(OrdinateSet_XYZ, coordinates);
      CoordinateSequence seqLs2DM = createSequence(OrdinateSet_XYM, coordinates);
      CoordinateSequence seqLs3DM = createSequence(OrdinateSet_XYZM, coordinates);

      // act
      LineString ls2D = reader2D.read("LINESTRING (10 10, 20 20, 30 40)");
      LineString ls2DE = reader2D.read("LINESTRING EMPTY");
      LineString ls3D = reader3D.read("LINESTRING Z(10 10 10, 20 20 10, 30 40 10)");
      LineString ls2DM = reader2DM.read("LINESTRING M(10 10 11, 20 20 11, 30 40 11)");
      LineString ls3DM = reader3DM.read("LINESTRING ZM(10 10 10 11, 20 20 10 11, 30 40 10 11)");

      // assert
      var ls2DCS = ls2D.getCoordinateSequence();
      expect(checkEqual(seqLs2D, ls2DCS), true);
      expect(checkEqual(seqLs2DE, ls2DE.getCoordinateSequence()), true);
      expect(checkEqual(seqLs3D, ls3D.getCoordinateSequence()), true);
      expect(checkEqual(seqLs2DM, ls2DM.getCoordinateSequence()), true);
      var ls3DMCS = ls3DM.getCoordinateSequence();
      expect(checkEqual(seqLs3DM, ls3DMCS), true);
    });
    test("testReadLinearRing", () {
      List<double> coordinates = [10, 10, 20, 20, 30, 40, 10, 10];
      CoordinateSequence seqLs2D = createSequence(OrdinateSet_XY, coordinates);
      CoordinateSequence seqLs2DE = createSequence(OrdinateSet_XY, <double>[]);
      CoordinateSequence seqLs3D = createSequence(OrdinateSet_XYZ, coordinates);
      CoordinateSequence seqLs2DM = createSequence(OrdinateSet_XYM, coordinates);
      CoordinateSequence seqLs3DM = createSequence(OrdinateSet_XYZM, coordinates);

      // act
      LineString ls2D = reader2D.read("LINEARRING (10 10, 20 20, 30 40, 10 10)");
      LineString ls2DE = reader2D.read("LINEARRING EMPTY");
      LineString ls3D = reader3D.read("LINEARRING Z(10 10 10, 20 20 10, 30 40 10, 10 10 10)");
      LineString ls2DM = reader2DM.read("LINEARRING M(10 10 11, 20 20 11, 30 40 11, 10 10 11)");
      LineString ls3DM = reader3DM.read("LINEARRING ZM(10 10 10 11, 20 20 10 11, 30 40 10 11, 10 10 10 11)");

      // assert
      expect(checkEqual(seqLs2D, ls2D.getCoordinateSequence()), true);
      expect(checkEqual(seqLs2DE, ls2DE.getCoordinateSequence()), true);
      expect(checkEqual(seqLs3D, ls3D.getCoordinateSequence()), true);
      expect(checkEqual(seqLs2DM, ls2DM.getCoordinateSequence()), true);
      expect(checkEqual(seqLs3DM, ls3DM.getCoordinateSequence()), true);

      expect(() => reader2D.read("LINEARRING (10 10, 20 20, 30 40, 10 99)"), throwsArgumentError);
    });
    test("testReadPolygon", () {
      List<double> shell = [10, 10, 10, 20, 20, 20, 20, 15, 10, 10];
      List<double> ring1 = [11, 11, 12, 11, 12, 12, 12, 11, 11, 11];
      List<double> ring2 = [11, 19, 11, 18, 12, 18, 12, 19, 11, 19];

      List<CoordinateSequence> csPoly2D = [createSequence(OrdinateSet_XY, shell), createSequence(OrdinateSet_XY, ring1), createSequence(OrdinateSet_XY, ring2)];
      CoordinateSequence csPoly2DE = createSequence(OrdinateSet_XY, <double>[]);
      List<CoordinateSequence> csPoly3D = [
        createSequence(OrdinateSet_XYZ, shell),
        createSequence(OrdinateSet_XYZ, ring1),
        createSequence(OrdinateSet_XYZ, ring2)
      ];
      List<CoordinateSequence> csPoly2DM = [
        createSequence(OrdinateSet_XYM, shell),
        createSequence(OrdinateSet_XYM, ring1),
        createSequence(OrdinateSet_XYM, ring2)
      ];
      List<CoordinateSequence> csPoly3DM = [
        createSequence(OrdinateSet_XYZM, shell),
        createSequence(OrdinateSet_XYZM, ring1),
        createSequence(OrdinateSet_XYZM, ring2)
      ];

      WKTReader rdr = reader2D;
      List<Polygon> poly2D = [
        rdr.read("POLYGON ((10 10, 10 20, 20 20, 20 15, 10 10))"),
        rdr.read("POLYGON ((10 10, 10 20, 20 20, 20 15, 10 10), (11 11, 12 11, 12 12, 12 11, 11 11))"),
        rdr.read("POLYGON ((10 10, 10 20, 20 20, 20 15, 10 10), (11 11, 12 11, 12 12, 12 11, 11 11), (11 19, 11 18, 12 18, 12 19, 11 19))")
      ];
      Polygon poly2DE = rdr.read("POLYGON EMPTY");
      rdr = reader3D;
      List<Polygon> poly3D = [
        rdr.read("POLYGON Z((10 10 10, 10 20 10, 20 20 10, 20 15 10, 10 10 10))"),
        rdr.read("POLYGON Z((10 10 10, 10 20 10, 20 20 10, 20 15 10, 10 10 10), (11 11 10, 12 11 10, 12 12 10, 12 11 10, 11 11 10))"),
        rdr.read(
            "POLYGON Z((10 10 10, 10 20 10, 20 20 10, 20 15 10, 10 10 10), (11 11 10, 12 11 10, 12 12 10, 12 11 10, 11 11 10), (11 19 10, 11 18 10, 12 18 10, 12 19 10, 11 19 10))")
      ];
      rdr = reader2DM;
      List<Polygon> poly2DM = [
        rdr.read("POLYGON M((10 10 11, 10 20 11, 20 20 11, 20 15 11, 10 10 11))"),
        rdr.read("POLYGON M((10 10 11, 10 20 11, 20 20 11, 20 15 11, 10 10 11), (11 11 11, 12 11 11, 12 12 11, 12 11 11, 11 11 11))"),
        rdr.read(
            "POLYGON M((10 10 11, 10 20 11, 20 20 11, 20 15 11, 10 10 11), (11 11 11, 12 11 11, 12 12 11, 12 11 11, 11 11 11), (11 19 11, 11 18 11, 12 18 11, 12 19 11, 11 19 11))")
      ];
      rdr = reader3DM;
      List<Polygon> poly3DM = [
        rdr.read("POLYGON ZM((10 10 10 11, 10 20 10 11, 20 20 10 11, 20 15 10 11, 10 10 10 11))"),
        rdr.read(
            "POLYGON ZM((10 10 10 11, 10 20 10 11, 20 20 10 11, 20 15 10 11, 10 10 10 11), (11 11 10 11, 12 11 10 11, 12 12 10 11, 12 11 10 11, 11 11 10 11))"),
        rdr.read(
            "POLYGON ZM((10 10 10 11, 10 20 10 11, 20 20 10 11, 20 15 10 11, 10 10 10 11), (11 11 10 11, 12 11 10 11, 12 12 10 11, 12 11 10 11, 11 11 10 11), (11 19 10 11, 11 18 10 11, 12 18 10 11, 12 19 10 11, 11 19 10 11))")
      ];
      // assert
      expect(checkEqual(csPoly2D[0], poly2D[2].getExteriorRing().getCoordinateSequence()), true);
      expect(checkEqual(csPoly2D[1], poly2D[2].getInteriorRingN(0).getCoordinateSequence()), true);
      expect(checkEqual(csPoly2D[2], poly2D[2].getInteriorRingN(1).getCoordinateSequence()), true);
      expect(checkEqualWithDim(csPoly2DE, poly2DE.getExteriorRing().getCoordinateSequence(), 2), true);
      expect(checkEqual(csPoly3D[0], poly3D[2].getExteriorRing().getCoordinateSequence()), true);
      expect(checkEqual(csPoly3D[1], poly3D[2].getInteriorRingN(0).getCoordinateSequence()), true);
      expect(checkEqual(csPoly3D[2], poly3D[2].getInteriorRingN(1).getCoordinateSequence()), true);
      expect(checkEqual(csPoly2DM[0], poly2DM[2].getExteriorRing().getCoordinateSequence()), true);
      expect(checkEqual(csPoly2DM[1], poly2DM[2].getInteriorRingN(0).getCoordinateSequence()), true);
      expect(checkEqual(csPoly2DM[2], poly2DM[2].getInteriorRingN(1).getCoordinateSequence()), true);
      expect(checkEqual(csPoly3DM[0], poly3DM[2].getExteriorRing().getCoordinateSequence()), true);
      expect(checkEqual(csPoly3DM[1], poly3DM[2].getInteriorRingN(0).getCoordinateSequence()), true);
      expect(checkEqual(csPoly3DM[2], poly3DM[2].getInteriorRingN(1).getCoordinateSequence()), true);
    });
    test("testReadMultiPoint", () {
      // arrange
      List<List<double>> coordinates = [
        [10, 10],
        [20, 20]
      ];
      List<CoordinateSequence> csMP2D = [createSequence(OrdinateSet_XY, coordinates[0]), createSequence(OrdinateSet_XY, coordinates[1])];
      List<CoordinateSequence> csMP3D = [createSequence(OrdinateSet_XYZ, coordinates[0]), createSequence(OrdinateSet_XYZ, coordinates[1])];
      List<CoordinateSequence> csMP2DM = [createSequence(OrdinateSet_XYM, coordinates[0]), createSequence(OrdinateSet_XYM, coordinates[1])];
      List<CoordinateSequence> csMP3DM = [createSequence(OrdinateSet_XYZM, coordinates[0]), createSequence(OrdinateSet_XYZM, coordinates[1])];

      // act
      WKTReader rdr = reader2D;
      MultiPoint mP2D = rdr.read("MULTIPOINT ((10 10), (20 20))");
      MultiPoint mP2DE = rdr.read("MULTIPOINT EMPTY");
      rdr = reader3D;
      MultiPoint mP3D = rdr.read("MULTIPOINT Z((10 10 10), (20 20 10))");
      rdr = reader2DM;
      MultiPoint mP2DM = rdr.read("MULTIPOINT M((10 10 11), (20 20 11))");
      rdr = reader3DM;
      MultiPoint mP3DM = rdr.read("MULTIPOINT ZM((10 10 10 11), (20 20 10 11))");

      // assert
      expect(checkEqual(csMP2D[0], (mP2D.getGeometryN(0) as Point).getCoordinateSequence()), true);
      expect(checkEqual(csMP2D[1], (mP2D.getGeometryN(1) as Point).getCoordinateSequence()), true);
      expect(mP2DE.isEmpty(), true);
      expect(mP2DE.getNumGeometries() == 0, true);
      expect(checkEqual(csMP3D[0], (mP3D.getGeometryN(0) as Point).getCoordinateSequence()), true);
      expect(checkEqual(csMP3D[1], (mP3D.getGeometryN(1) as Point).getCoordinateSequence()), true);
      expect(checkEqual(csMP2DM[0], (mP2DM.getGeometryN(0) as Point).getCoordinateSequence()), true);
      expect(checkEqual(csMP2DM[1], (mP2DM.getGeometryN(1) as Point).getCoordinateSequence()), true);
      var mp3DMCS = (mP3DM.getGeometryN(0) as Point).getCoordinateSequence();
      expect(checkEqual(csMP3DM[0], mp3DMCS), true);
      expect(checkEqual(csMP3DM[1], (mP3DM.getGeometryN(1) as Point).getCoordinateSequence()), true);
    });
    test("testReadMultiLineString", () {
      // arrange
      List<List<double>> coordinates = [
        [10, 10, 20, 20],
        [15, 15, 30, 15]
      ];
      List<CoordinateSequence> csMls2D = [createSequence(OrdinateSet_XY, coordinates[0]), createSequence(OrdinateSet_XY, coordinates[1])];
      List<CoordinateSequence> csMls3D = [createSequence(OrdinateSet_XYZ, coordinates[0]), createSequence(OrdinateSet_XYZ, coordinates[1])];
      List<CoordinateSequence> csMls2DM = [createSequence(OrdinateSet_XYM, coordinates[0]), createSequence(OrdinateSet_XYM, coordinates[1])];
      List<CoordinateSequence> csMls3DM = [createSequence(OrdinateSet_XYZM, coordinates[0]), createSequence(OrdinateSet_XYZM, coordinates[1])];

      // act
      WKTReader rdr = reader2D;
      MultiLineString mLs2D = rdr.read("MULTILINESTRING ((10 10, 20 20), (15 15, 30 15))");
      MultiLineString mLs2DE = rdr.read("MULTILINESTRING EMPTY");
      rdr = reader3D;
      MultiLineString mLs3D = rdr.read("MULTILINESTRING Z((10 10 10, 20 20 10), (15 15 10, 30 15 10))");
      rdr = reader2DM;
      MultiLineString mLs2DM = rdr.read("MULTILINESTRING M((10 10 11, 20 20 11), (15 15 11, 30 15 11))");
      rdr = reader3DM;
      MultiLineString mLs3DM = rdr.read("MULTILINESTRING ZM((10 10 10 11, 20 20 10 11), (15 15 10 11, 30 15 10 11))");

      // assert
      expect(checkEqual(csMls2D[0], (mLs2D.getGeometryN(0) as LineString).getCoordinateSequence()), true);
      expect(checkEqual(csMls2D[1], (mLs2D.getGeometryN(1) as LineString).getCoordinateSequence()), true);
      expect(mLs2DE.isEmpty(), true);
      expect(mLs2DE.getNumGeometries() == 0, true);
      expect(checkEqual(csMls3D[0], (mLs3D.getGeometryN(0) as LineString).getCoordinateSequence()), true);
      expect(checkEqual(csMls3D[1], (mLs3D.getGeometryN(1) as LineString).getCoordinateSequence()), true);
      expect(checkEqual(csMls2DM[0], (mLs2DM.getGeometryN(0) as LineString).getCoordinateSequence()), true);
      expect(checkEqual(csMls2DM[1], (mLs2DM.getGeometryN(1) as LineString).getCoordinateSequence()), true);
      expect(checkEqual(csMls3DM[0], (mLs3DM.getGeometryN(0) as LineString).getCoordinateSequence()), true);
      expect(checkEqual(csMls3DM[1], (mLs3DM.getGeometryN(1) as LineString).getCoordinateSequence()), true);
    });
    test("testReadMultiPolygon", () {
      List<double> shell1 = [10, 10, 10, 20, 20, 20, 20, 15, 10, 10];
      List<double> ring1 = [11, 11, 12, 11, 12, 12, 12, 11, 11, 11];
      List<double> shell2 = [60, 60, 70, 70, 80, 60, 60, 60];

      List<CoordinateSequence> csPoly2D = [
        createSequence(OrdinateSet_XY, shell1),
        createSequence(OrdinateSet_XY, ring1),
        createSequence(OrdinateSet_XY, shell2)
      ];
      List<CoordinateSequence> csPoly3D = [
        createSequence(OrdinateSet_XYZ, shell1),
        createSequence(OrdinateSet_XYZ, ring1),
        createSequence(OrdinateSet_XYZ, shell2)
      ];
      List<CoordinateSequence> csPoly2DM = [
        createSequence(OrdinateSet_XYM, shell1),
        createSequence(OrdinateSet_XYM, ring1),
        createSequence(OrdinateSet_XYM, shell2)
      ];
      List<CoordinateSequence> csPoly3DM = [
        createSequence(OrdinateSet_XYZM, shell1),
        createSequence(OrdinateSet_XYZM, ring1),
        createSequence(OrdinateSet_XYZM, shell2)
      ];

      WKTReader rdr = reader2D;
      List<MultiPolygon> poly2D = [
        rdr.read("MULTIPOLYGON (((10 10, 10 20, 20 20, 20 15, 10 10)))"),
        rdr.read("MULTIPOLYGON (((10 10, 10 20, 20 20, 20 15, 10 10), (11 11, 12 11, 12 12, 12 11, 11 11)))"),
        rdr.read("MULTIPOLYGON (((10 10, 10 20, 20 20, 20 15, 10 10), (11 11, 12 11, 12 12, 12 11, 11 11)), ((60 60, 70 70, 80 60, 60 60)))")
      ];
      MultiPolygon poly2DE = rdr.read("MULTIPOLYGON EMPTY");
      rdr = reader3D;
      List<MultiPolygon> poly3D = [
        rdr.read("MULTIPOLYGON Z(((10 10 10, 10 20 10, 20 20 10, 20 15 10, 10 10 10)))"),
        rdr.read("MULTIPOLYGON Z(((10 10 10, 10 20 10, 20 20 10, 20 15 10, 10 10 10), (11 11 10, 12 11 10, 12 12 10, 12 11 10, 11 11 10)))"),
        rdr.read(
            "MULTIPOLYGON Z(((10 10 10, 10 20 10, 20 20 10, 20 15 10, 10 10 10), (11 11 10, 12 11 10, 12 12 10, 12 11 10, 11 11 10)), ((60 60 10, 70 70 10, 80 60 10, 60 60 10)))")
      ];
      List<MultiPolygon> poly2DM = [
        rdr.read("MULTIPOLYGON M(((10 10 11, 10 20 11, 20 20 11, 20 15 11, 10 10 11)))"),
        rdr.read("MULTIPOLYGON M(((10 10 11, 10 20 11, 20 20 11, 20 15 11, 10 10 11), (11 11 11, 12 11 11, 12 12 11, 12 11 11, 11 11 11)))"),
        rdr.read(
            "MULTIPOLYGON M(((10 10 11, 10 20 11, 20 20 11, 20 15 11, 10 10 11), (11 11 11, 12 11 11, 12 12 11, 12 11 11, 11 11 11)), ((60 60 11, 70 70 11, 80 60 11, 60 60 11)))")
      ];
      rdr = reader3DM;
      List<MultiPolygon> poly3DM = [
        rdr.read("MULTIPOLYGON ZM(((10 10 10 11, 10 20 10 11, 20 20 10 11, 20 15 10 11, 10 10 10 11)))"),
        rdr.read(
            "MULTIPOLYGON ZM(((10 10 10 11, 10 20 10 11, 20 20 10 11, 20 15 10 11, 10 10 10 11), (11 11 10 11, 12 11 10 11, 12 12 10 11, 12 11 10 11, 11 11 10 11)))"),
        rdr.read(
            "MULTIPOLYGON ZM(((10 10 10 11, 10 20 10 11, 20 20 10 11, 20 15 10 11, 10 10 10 11), (11 11 10 11, 12 11 10 11, 12 12 10 11, 12 11 10 11, 11 11 10 11)), ((60 60 10 11, 70 70 10 11, 80 60 10 11, 60 60 10 11)))")
      ];

      // assert
      expect(checkEqual(csPoly2D[0], (poly2D[2].getGeometryN(0) as Polygon).getExteriorRing().getCoordinateSequence()), true);
      expect(checkEqual(csPoly2D[1], (poly2D[2].getGeometryN(0) as Polygon).getInteriorRingN(0).getCoordinateSequence()), true);
      expect(checkEqual(csPoly2D[2], (poly2D[2].getGeometryN(1) as Polygon).getExteriorRing().getCoordinateSequence()), true);
      expect(poly2DE.isEmpty(), true);
      expect(poly2DE.getNumGeometries() == 0, true);

      expect(checkEqual(csPoly3D[0], (poly3D[2].getGeometryN(0) as Polygon).getExteriorRing().getCoordinateSequence()), true);
      expect(checkEqual(csPoly3D[1], (poly3D[2].getGeometryN(0) as Polygon).getInteriorRingN(0).getCoordinateSequence()), true);
      expect(checkEqual(csPoly3D[2], (poly3D[2].getGeometryN(1) as Polygon).getExteriorRing().getCoordinateSequence()), true);

      var poly2DMExtRing = (poly2DM[2].getGeometryN(0) as Polygon).getExteriorRing();
      var poly2DMExtRingCS = poly2DMExtRing.getCoordinateSequence();
      expect(checkEqual(csPoly2DM[0], poly2DMExtRingCS), true);
      expect(checkEqual(csPoly2DM[1], (poly2DM[2].getGeometryN(0) as Polygon).getInteriorRingN(0).getCoordinateSequence()), true);
      expect(checkEqual(csPoly2DM[2], (poly2DM[2].getGeometryN(1) as Polygon).getExteriorRing().getCoordinateSequence()), true);

      expect(checkEqual(csPoly3DM[0], (poly3DM[2].getGeometryN(0) as Polygon).getExteriorRing().getCoordinateSequence()), true);
      expect(checkEqual(csPoly3DM[1], (poly3DM[2].getGeometryN(0) as Polygon).getInteriorRingN(0).getCoordinateSequence()), true);
      expect(checkEqual(csPoly3DM[2], (poly3DM[2].getGeometryN(1) as Polygon).getExteriorRing().getCoordinateSequence()), true);
    });
    test("", () {});
    test("", () {});
    test("", () {});
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
      if (val1 != null && val2 == null) {
        return false;
      } else if (val1 == null && val2 != null) {
        return false;
      } else if (val1.isNaN) {
        if (!val2.isNaN) return false;
      } else if ((val1 - val2).abs() > tolerance) return false;
    }
  }

  return true;
}

class WKBTest {
  static GeometryFactory geomFactory = new GeometryFactory.defaultPrecision();
  WKTReader rdr = new WKTReader.withFactory(geomFactory);

  void testFirst() {
    runWKBTestWKT("MULTIPOINT ((0 0), (1 4), (100 200))");
  }

  void testPointPCS() {
    runWKBTestPackedCoordinate("POINT (1 2)");
  }

  void testPoint() {
    runWKBTestWKT("POINT (1 2)");
  }

  void testLineString() {
    runWKBTestWKT("LINESTRING (1 2, 10 20, 100 200)");
  }

  void testPolygon() {
    runWKBTestWKT("POLYGON ((0 0, 100 0, 100 100, 0 100, 0 0))");
  }

  void testPolygonWithHole() {
    runWKBTestWKT("POLYGON ((0 0, 100 0, 100 100, 0 100, 0 0), (1 1, 1 10, 10 10, 10 1, 1 1) )");
  }

  void testMultiPoint() {
    runWKBTestWKT("MULTIPOINT ((0 0), (1 4), (100 200))");
  }

  void testMultiLineString() {
    runWKBTestWKT("MULTILINESTRING ((0 0, 1 10), (10 10, 20 30), (123 123, 456 789))");
  }

  void testMultiPolygon() {
    runWKBTestWKT("MULTIPOLYGON ( ((0 0, 100 0, 100 100, 0 100, 0 0), (1 1, 1 10, 10 10, 10 1, 1 1) ), ((200 200, 200 250, 250 250, 250 200, 200 200)) )");
  }

  void testGeometryCollection() {
    runWKBTestWKT("GEOMETRYCOLLECTION ( POINT ( 1 1), LINESTRING (0 0, 10 10), POLYGON ((0 0, 100 0, 100 100, 0 100, 0 0)) )");
  }

  void testNestedGeometryCollection() {
    runWKBTestWKT(
        "GEOMETRYCOLLECTION ( POINT (20 20), GEOMETRYCOLLECTION ( POINT ( 1 1), LINESTRING (0 0, 10 10), POLYGON ((0 0, 100 0, 100 100, 0 100, 0 0)) ) )");
  }

  void testLineStringEmpty() {
    runWKBTestWKT("LINESTRING EMPTY");
  }

// void testBigPolygon()
//
//{
//GeometricShapeFactory shapeFactory = new GeometricShapeFactory(geomFactory);
//shapeFactory.setBase(new Coordinate(0,0));
//shapeFactory.setSize(1000);
//shapeFactory.setNumPoints(1000);
//Geometry geom = shapeFactory.createRectangle();
//runWKBTest(geom, 2, false);
//}

  void testPolygonEmpty() {
    runWKBTestWKT("POLYGON EMPTY");
  }

  void testMultiPointEmpty() {
    runWKBTestWKT("MULTIPOINT EMPTY");
  }

  void testMultiLineStringEmpty() {
    runWKBTestWKT("MULTILINESTRING EMPTY");
  }

  void testMultiPolygonEmpty() {
    runWKBTestWKT("MULTIPOLYGON EMPTY");
  }

  void testGeometryCollectionEmpty() {
    runWKBTestWKT("GEOMETRYCOLLECTION EMPTY");
  }

  void runWKBTestWKT(String wkt) {
    runWKBTestCoordinateArray(wkt);
    runWKBTestPackedCoordinate(wkt);
  }

  void runWKBTestPackedCoordinate(String wkt) {
    GeometryFactory geomFactory =
        new GeometryFactory.withCoordinateSequenceFactory(new PackedCoordinateSequenceFactory.withType(PackedCoordinateSequenceFactory.DOUBLE));
    WKTReader rdr = new WKTReader.withFactory(geomFactory);
    Geometry g = rdr.read(wkt);

// Since we are using a PCS of dim=2, only check 2-dimensional storage
    runWKBTest(g, 2, true);
    runWKBTest(g, 2, false);
  }

  void runWKBTestCoordinateArray(String wkt) {
    GeometryFactory geomFactory = new GeometryFactory.defaultPrecision();
    WKTReader rdr = new WKTReader.withFactory(geomFactory);
    Geometry g = rdr.read(wkt);

// CoordinateArrays support dimension 3, so test both dimensions
    runWKBTest(g, 2, true);
    runWKBTest(g, 2, false);
    runWKBTest(g, 3, true);
    runWKBTest(g, 3, false);
  }

  void runWKBTest(Geometry g, int dimension, bool toHex) {
    setZ(g);
    runWKBTestWithBO(g, dimension, Endian.little, toHex);
    runWKBTestWithBO(g, dimension, Endian.big, toHex);
  }

  void runWKBTestWithBO(Geometry g, int dimension, Endian byteOrder, bool toHex) {
    runGeometry(g, dimension, byteOrder, toHex, 100);
    runGeometry(g, dimension, byteOrder, toHex, 0);
    runGeometry(g, dimension, byteOrder, toHex, 101010);
    runGeometry(g, dimension, byteOrder, toHex, -1);
  }

  void setZ(Geometry g) {
    g.applyCF((coord) {
      coord.setZ((coord.x + coord.y) / 2);
    });
  }

//static Comparator comp2D = new Coordinate.DimensionalComparator();
//static Comparator comp3D = new Coordinate.DimensionalComparator(3);

  static Comparator<CoordinateSequence> comp2 = CoordinateSequenceComparatorBuilder.withLimit(2);
  static Comparator<CoordinateSequence> comp3 = CoordinateSequenceComparatorBuilder.withLimit(3);

  /**
   * Use single WKB reader, to ensure it can be used for multiple input geometries
   */
  WKBReader wkbReader = new WKBReader.withFactory(geomFactory);

  void runGeometry(Geometry g, int dimension, Endian byteOrder, bool toHex, int srid) {
    bool includeSRID = false;
    if (srid >= 0) {
      includeSRID = true;
      g.setSRID(srid);
    }

    WKBWriter wkbWriter = new WKBWriter.withDimOrderSrid(dimension, byteOrder, includeSRID);
    List<int> wkb = wkbWriter.write(g);
    String wkbHex = null;
    if (toHex) wkbHex = WKBWriter.toHex(wkb);

    if (toHex) wkb = WKBReader.hexToBytes(wkbHex);
    Geometry g2 = wkbReader.read(wkb);

    Comparator<CoordinateSequence> comp = (dimension == 2) ? comp2 : comp3;
    bool isEqual = (g.compareToWithComparator(g2, comp) == 0);
    assertTrue(isEqual);

    if (includeSRID) {
      bool isSRIDEqual = g.getSRID() == g2.getSRID();
      assertTrue(isSRIDEqual);
    }
  }
}
