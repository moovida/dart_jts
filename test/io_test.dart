import "package:test/test.dart";
import 'package:dart_jts/dart_jts.dart';
import 'dart:math' as math;
import 'dart:typed_data';
import 'testing_utilities.dart';

void main() {
  group("WKBReaderTest - ", () {
    WKBReaderTest testClass = WKBReaderTest();

    test("testShortPolygons", () => testClass.testShortPolygons());
    test("testSinglePointLineString",
        () => testClass.testSinglePointLineString());
    test("testSpatialiteMultiGeometry",
        () => testClass.testSpatialiteMultiGeometry());
    test("test2dSpatialiteWKB", () => testClass.test2dSpatialiteWKB());
    test("testSpatialiteWKB_Z", () => testClass.testSpatialiteWKB_Z());
    test("testSpatialiteWKB_M", () => testClass.testSpatialiteWKB_M());
    test("testSpatialiteWKB_ZM", () => testClass.testSpatialiteWKB_ZM());
//    test("XXtestIllFormedWKB", () => testClass.XXtestIllFormedWKB());
  });

  group("WKBWriterTest - ", () {
    test("testSRID", () {
      GeometryFactory gf = new GeometryFactory.defaultPrecision();
      Point p1 = gf.createPoint(new Coordinate(1, 2));
      p1.setSRID(1234);

      //first write out without srid set
      WKBWriter w = new WKBWriter();
      List<int> wkb = w.write(p1);

      //check the 3rd bit of the second byte, should be unset
      int b = (wkb[1] & 0x20);
      assertEquals(0, b);

      //read geometry back in
      WKBReader r = new WKBReader.withFactory(gf);
      Geometry p2 = r.read(Uint8List.fromList(wkb));

      assertTrue(p1.equalsExactGeom(p2));
      assertEquals(0, p2.getSRID());

      //not write out with srid set
      w = new WKBWriter.withDimSrid(2, true);
      wkb = w.write(p1);

      //check the 3rd bit of the second byte, should be set
      b = (wkb[1] & 0x20);
      assertEquals(0x20, b);

      int srid = ((wkb[5] & 0xff) << 24) |
          ((wkb[6] & 0xff) << 16) |
          ((wkb[7] & 0xff) << 8) |
          ((wkb[8] & 0xff));

      assertEquals(1234, srid);

      r = new WKBReader.withFactory(gf);
      p2 = r.read(Uint8List.fromList(wkb));

      //read the geometry back in
      assertTrue(p1.equalsExactGeom(p2));
      assertEquals(1234, p2.getSRID());
    });
  });

  group("WKBTest - ", () {
    WKBTest testClass = WKBTest();
    test("testFirst", () => testClass.testFirst());
    test("testPointPCS", () => testClass.testPointPCS());
    test("testPoint", () => testClass.testPoint());
    test("testLineString", () => testClass.testLineString());
    test("testPolygon", () => testClass.testPolygon());
    test("testPolygonWithHole", () => testClass.testPolygonWithHole());
    test("testMultiPoint", () => testClass.testMultiPoint());
    test("testMultiLineString", () => testClass.testMultiLineString());
    test("testMultiPolygon", () => testClass.testMultiPolygon());
    test("testGeometryCollection", () => testClass.testGeometryCollection());
    test("testNestedGeometryCollection",
        () => testClass.testNestedGeometryCollection());
    test("testLineStringEmpty", () => testClass.testLineStringEmpty());
    test("testPolygonEmpty", () => testClass.testPolygonEmpty());
    test("testMultiPointEmpty", () => testClass.testMultiPointEmpty());
    test(
        "testMultiLineStringEmpty", () => testClass.testMultiLineStringEmpty());
    test("testMultiPolygonEmpty", () => testClass.testMultiPolygonEmpty());
    test("testGeometryCollectionEmpty",
        () => testClass.testGeometryCollectionEmpty());

    test("test Spatialite WKB Geometries reading and writing", () {
      var geomStrings = [
        "POLYGON ((71 70, 40 70, 40 40, 5 40, 5 15, 15 15, 15 4, 50 4, 71 70))",
        "POLYGON ((10 42, 11.9 42, 11.9 40, 10 40, 10 42))",
        "POLYGON ((11.1 43.2, 11.3 41.3, 13.9 41, 13.8 43.2, 11.1 43.2))",
        "LINESTRING (11.3 44.3, 8.3 41.4, 11.4 38.1, 14.9 41.3)",
        "POINT (12.7 44.2)",
        "POINT (15.1 43.3)",
        "POINT (15 40.4)",
        "POINT (13.2 38.4)",
        "MULTIPOLYGON (((6.9 45.9, 8.4 45.9, 8.4 44.3, 6.9 44.3, 6.9 45.9)), ((9.1 46.3, 10.8 46.3, 10.8 44.6, 9.1 44.6, 9.1 46.3)))",
        "MULTILINESTRING ((7.4 42.6, 7.4 39, 8.6 38.5), (8 40.3, 9.5 38.6, 8.4 37.5))",
        "MULTIPOINT ((6.8 42.5), (6.8 41.4), (6.6 40.2))",
      ];

      for (var geomString in geomStrings) {
//        print(geomString);
        var geom = WKTReader().read(geomString)!;
        var bytes = WKBWriter().write(geom, doSpatialite: true);
        var read = WKBReader().read(bytes, doSpatialite: true);
        expect(read.equalsExactGeom(read), true);
      }
    });
  });

  group("WKTWriterTest - ", () {
    PrecisionModel precisionModel = PrecisionModel.fixedPrecision(1);
    GeometryFactory geometryFactory =
        GeometryFactory.withPrecisionModelSrid(precisionModel, 0);
    WKTWriter writer = WKTWriter();
    WKTWriter writer3D = WKTWriter.withDimension(3);
    WKTWriter writer2DM = WKTWriter.withDimension(3);

    writer2DM.setOutputOrdinates(OrdinateSet_XYM);

    test("testProperties", () {
      expect(OrdinateSet_XY, writer.getOutputOrdinates());
      expect(OrdinateSet_XYZ, writer3D.getOutputOrdinates());
      expect(OrdinateSet_XYM, writer2DM.getOutputOrdinates());

      GeometryFactory gf = GeometryFactory.withCoordinateSequenceFactory(
          PackedCoordinateSequenceFactory.DOUBLE_FACTORY);
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
      List<Coordinate> coordinates = [
        Coordinate.fromXYZ(10, 10, 0),
        Coordinate.fromXYZ(20, 20, 0),
        Coordinate.fromXYZ(30, 40, 0)
      ];
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
      expect("POLYGON ((10 10, 10 20, 20 20, 20 15, 10 10))",
          writer.write(polygon));
    });
    test("testWriteMultiPoint", () {
      List<Point> points = [
        geometryFactory.createPoint(Coordinate.fromXYZ(10, 10, 0)),
        geometryFactory.createPoint(Coordinate.fromXYZ(20, 20, 0))
      ];
      MultiPoint multiPoint = geometryFactory.createMultiPoint(points);
      expect("MULTIPOINT ((10 10), (20 20))", writer.write(multiPoint));
    });
    test("testWriteMultiLineString", () {
      List<Coordinate> coordinates1 = [
        Coordinate.fromXYZ(10, 10, 0),
        Coordinate.fromXYZ(20, 20, 0)
      ];
      LineString lineString1 = geometryFactory.createLineString(coordinates1);
      List<Coordinate> coordinates2 = [
        Coordinate.fromXYZ(15, 15, 0),
        Coordinate.fromXYZ(30, 15, 0)
      ];
      LineString lineString2 = geometryFactory.createLineString(coordinates2);
      List<LineString> lineStrings = [lineString1, lineString2];
      MultiLineString multiLineString =
          geometryFactory.createMultiLineString(lineStrings);
      expect("MULTILINESTRING ((10 10, 20 20), (15 15, 30 15))",
          writer.write(multiLineString));
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
      expect(
          "MULTIPOLYGON (((10 10, 10 20, 20 20, 20 15, 10 10)), ((60 60, 70 70, 80 60, 60 60)))",
          writer.write(multiPolygon));
    });
    test("testWriteGeometryCollection", () {
      Point point1 = geometryFactory.createPoint(Coordinate(10, 10));
      Point point2 = geometryFactory.createPoint(Coordinate(30, 30));
      List<Coordinate> coordinates = [
        Coordinate.fromXYZ(15, 15, 0),
        Coordinate.fromXYZ(20, 20, 0)
      ];
      LineString lineString1 = geometryFactory.createLineString(coordinates);
      List<Geometry> geometries = [point1, point2, lineString1];
      GeometryCollection geometryCollection =
          geometryFactory.createGeometryCollection(geometries);
      expect(
          "GEOMETRYCOLLECTION (POINT (10 10), POINT (30 30), LINESTRING (15 15, 20 20))",
          writer.write(geometryCollection));
    });
    test("testWriteLargeNumbers1", () {
      PrecisionModel precisionModel = PrecisionModel.fixedPrecision(1E9);
      GeometryFactory geometryFactory =
          GeometryFactory.withPrecisionModelSrid(precisionModel, 0);
      Point point1 =
          geometryFactory.createPoint(Coordinate(123456789012345680, 10E9));
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
      List<Coordinate> coordinates = [
        Coordinate(1, 1),
        Coordinate.fromXYZ(2, 2, 2)
      ];
      LineString line = geometryFactory.createLineString(coordinates);
      String wkt = writer3D.write(line);
      expect("LINESTRING Z(1 1 NaN, 2 2 2)", wkt);
      wkt = writer2DM.write(line);
      expect("LINESTRING (1 1, 2 2)", wkt);
    });
  });

  WKTReader reader2D = getWKTReaderFromOrdinateSetAndScale(OrdinateSet_XY, 1.0);
  reader2D.setIsOldJtsCoordinateSyntaxAllowed(false);
  WKTReader reader2DOld =
      getWKTReaderFromOrdinateSetAndScale(OrdinateSet_XY, 1.0);
  reader2DOld.setIsOldJtsCoordinateSyntaxAllowed(true);
  WKTReader reader3D =
      getWKTReaderFromOrdinateSetAndScale(OrdinateSet_XYZ, 1.0);
  WKTReader reader2DM =
      getWKTReaderFromOrdinateSetAndScale(OrdinateSet_XYM, 1.0);
  WKTReader reader3DM =
      getWKTReaderFromOrdinateSetAndScale(OrdinateSet_XYZM, 1.0);
  group("WKTReaderTest - ", () {
    test("testReadNaN", () {
      // arrange
      CoordinateSequence seq = createSequence(OrdinateSet_XYZ, [10, 10]);
      seq.setOrdinate(0, CoordinateSequence.Z, double.nan);

      // act
      Point pt1 = reader2DOld.read("POINT (10 10 NaN)")! as Point;
      Point pt2 = reader2DOld.read("POINT (10 10 nan)")! as Point;
      Point pt3 = reader2DOld.read("POINT (10 10 NAN)")! as Point;

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
      CoordinateSequence seqPt2DM =
          createSequence(OrdinateSet_XYM, coordinates);
      CoordinateSequence seqPt3DM =
          createSequence(OrdinateSet_XYZM, coordinates);

      // act
      Point pt2D = reader2D.read("POINT (10 10)")! as Point;
      Point pt2DE = reader2D.read("POINT EMPTY")! as Point;
      Point pt3D = reader3D.read("POINT Z(10 10 10)")! as Point;
      Point pt2DM = reader2DM.read("POINT M(10 10 11)")! as Point;
      Point pt3DM = reader3DM.read("POINT ZM(10 10 10 11)")! as Point;

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
      CoordinateSequence seqLs2DM =
          createSequence(OrdinateSet_XYM, coordinates);
      CoordinateSequence seqLs3DM =
          createSequence(OrdinateSet_XYZM, coordinates);

      // act
      LineString ls2D =
          reader2D.read("LINESTRING (10 10, 20 20, 30 40)")! as LineString;
      LineString ls2DE = reader2D.read("LINESTRING EMPTY")! as LineString;
      LineString ls3D = reader3D
          .read("LINESTRING Z(10 10 10, 20 20 10, 30 40 10)")! as LineString;
      LineString ls2DM = reader2DM
          .read("LINESTRING M(10 10 11, 20 20 11, 30 40 11)")! as LineString;
      LineString ls3DM = reader3DM
              .read("LINESTRING ZM(10 10 10 11, 20 20 10 11, 30 40 10 11)")!
          as LineString;

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
      CoordinateSequence seqLs2DM =
          createSequence(OrdinateSet_XYM, coordinates);
      CoordinateSequence seqLs3DM =
          createSequence(OrdinateSet_XYZM, coordinates);

      // act
      LineString ls2D = reader2D
          .read("LINEARRING (10 10, 20 20, 30 40, 10 10)")! as LineString;
      LineString ls2DE = reader2D.read("LINEARRING EMPTY")! as LineString;
      LineString ls3D =
          reader3D.read("LINEARRING Z(10 10 10, 20 20 10, 30 40 10, 10 10 10)")!
              as LineString;
      LineString ls2DM = reader2DM
              .read("LINEARRING M(10 10 11, 20 20 11, 30 40 11, 10 10 11)")!
          as LineString;
      LineString ls3DM = reader3DM.read(
              "LINEARRING ZM(10 10 10 11, 20 20 10 11, 30 40 10 11, 10 10 10 11)")!
          as LineString;

      // assert
      expect(checkEqual(seqLs2D, ls2D.getCoordinateSequence()), true);
      expect(checkEqual(seqLs2DE, ls2DE.getCoordinateSequence()), true);
      expect(checkEqual(seqLs3D, ls3D.getCoordinateSequence()), true);
      expect(checkEqual(seqLs2DM, ls2DM.getCoordinateSequence()), true);
      expect(checkEqual(seqLs3DM, ls3DM.getCoordinateSequence()), true);

      expect(() => reader2D.read("LINEARRING (10 10, 20 20, 30 40, 10 99)"),
          throwsArgumentError);
    });
    test("testReadPolygon", () {
      List<double> shell = [10, 10, 10, 20, 20, 20, 20, 15, 10, 10];
      List<double> ring1 = [11, 11, 12, 11, 12, 12, 12, 11, 11, 11];
      List<double> ring2 = [11, 19, 11, 18, 12, 18, 12, 19, 11, 19];

      List<CoordinateSequence> csPoly2D = [
        createSequence(OrdinateSet_XY, shell),
        createSequence(OrdinateSet_XY, ring1),
        createSequence(OrdinateSet_XY, ring2)
      ];
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
        rdr.read("POLYGON ((10 10, 10 20, 20 20, 20 15, 10 10))")! as Polygon,
        rdr.read(
                "POLYGON ((10 10, 10 20, 20 20, 20 15, 10 10), (11 11, 12 11, 12 12, 12 11, 11 11))")!
            as Polygon,
        rdr.read(
                "POLYGON ((10 10, 10 20, 20 20, 20 15, 10 10), (11 11, 12 11, 12 12, 12 11, 11 11), (11 19, 11 18, 12 18, 12 19, 11 19))")!
            as Polygon
      ];
      Polygon poly2DE = rdr.read("POLYGON EMPTY")! as Polygon;
      rdr = reader3D;
      List<Polygon> poly3D = [
        rdr.read(
                "POLYGON Z((10 10 10, 10 20 10, 20 20 10, 20 15 10, 10 10 10))")!
            as Polygon,
        rdr.read(
                "POLYGON Z((10 10 10, 10 20 10, 20 20 10, 20 15 10, 10 10 10), (11 11 10, 12 11 10, 12 12 10, 12 11 10, 11 11 10))")!
            as Polygon,
        rdr.read(
                "POLYGON Z((10 10 10, 10 20 10, 20 20 10, 20 15 10, 10 10 10), (11 11 10, 12 11 10, 12 12 10, 12 11 10, 11 11 10), (11 19 10, 11 18 10, 12 18 10, 12 19 10, 11 19 10))")!
            as Polygon
      ];
      rdr = reader2DM;
      List<Polygon> poly2DM = [
        rdr.read(
                "POLYGON M((10 10 11, 10 20 11, 20 20 11, 20 15 11, 10 10 11))")!
            as Polygon,
        rdr.read(
                "POLYGON M((10 10 11, 10 20 11, 20 20 11, 20 15 11, 10 10 11), (11 11 11, 12 11 11, 12 12 11, 12 11 11, 11 11 11))")!
            as Polygon,
        rdr.read(
                "POLYGON M((10 10 11, 10 20 11, 20 20 11, 20 15 11, 10 10 11), (11 11 11, 12 11 11, 12 12 11, 12 11 11, 11 11 11), (11 19 11, 11 18 11, 12 18 11, 12 19 11, 11 19 11))")!
            as Polygon
      ];
      rdr = reader3DM;
      List<Polygon> poly3DM = [
        rdr.read(
                "POLYGON ZM((10 10 10 11, 10 20 10 11, 20 20 10 11, 20 15 10 11, 10 10 10 11))")!
            as Polygon,
        rdr.read(
                "POLYGON ZM((10 10 10 11, 10 20 10 11, 20 20 10 11, 20 15 10 11, 10 10 10 11), (11 11 10 11, 12 11 10 11, 12 12 10 11, 12 11 10 11, 11 11 10 11))")!
            as Polygon,
        rdr.read(
                "POLYGON ZM((10 10 10 11, 10 20 10 11, 20 20 10 11, 20 15 10 11, 10 10 10 11), (11 11 10 11, 12 11 10 11, 12 12 10 11, 12 11 10 11, 11 11 10 11), (11 19 10 11, 11 18 10 11, 12 18 10 11, 12 19 10 11, 11 19 10 11))")!
            as Polygon
      ];
      // assert
      expect(
          checkEqual(
              csPoly2D[0], poly2D[2].getExteriorRing().getCoordinateSequence()),
          true);
      expect(
          checkEqual(csPoly2D[1],
              poly2D[2].getInteriorRingN(0).getCoordinateSequence()),
          true);
      expect(
          checkEqual(csPoly2D[2],
              poly2D[2].getInteriorRingN(1).getCoordinateSequence()),
          true);
      expect(
          checkEqualWithDim(
              csPoly2DE, poly2DE.getExteriorRing().getCoordinateSequence(), 2),
          true);
      expect(
          checkEqual(
              csPoly3D[0], poly3D[2].getExteriorRing().getCoordinateSequence()),
          true);
      expect(
          checkEqual(csPoly3D[1],
              poly3D[2].getInteriorRingN(0).getCoordinateSequence()),
          true);
      expect(
          checkEqual(csPoly3D[2],
              poly3D[2].getInteriorRingN(1).getCoordinateSequence()),
          true);
      expect(
          checkEqual(csPoly2DM[0],
              poly2DM[2].getExteriorRing().getCoordinateSequence()),
          true);
      expect(
          checkEqual(csPoly2DM[1],
              poly2DM[2].getInteriorRingN(0).getCoordinateSequence()),
          true);
      expect(
          checkEqual(csPoly2DM[2],
              poly2DM[2].getInteriorRingN(1).getCoordinateSequence()),
          true);
      expect(
          checkEqual(csPoly3DM[0],
              poly3DM[2].getExteriorRing().getCoordinateSequence()),
          true);
      expect(
          checkEqual(csPoly3DM[1],
              poly3DM[2].getInteriorRingN(0).getCoordinateSequence()),
          true);
      expect(
          checkEqual(csPoly3DM[2],
              poly3DM[2].getInteriorRingN(1).getCoordinateSequence()),
          true);
    });
    test("testReadMultiPoint", () {
      // arrange
      List<List<double>> coordinates = [
        [10, 10],
        [20, 20]
      ];
      List<CoordinateSequence> csMP2D = [
        createSequence(OrdinateSet_XY, coordinates[0]),
        createSequence(OrdinateSet_XY, coordinates[1])
      ];
      List<CoordinateSequence> csMP3D = [
        createSequence(OrdinateSet_XYZ, coordinates[0]),
        createSequence(OrdinateSet_XYZ, coordinates[1])
      ];
      List<CoordinateSequence> csMP2DM = [
        createSequence(OrdinateSet_XYM, coordinates[0]),
        createSequence(OrdinateSet_XYM, coordinates[1])
      ];
      List<CoordinateSequence> csMP3DM = [
        createSequence(OrdinateSet_XYZM, coordinates[0]),
        createSequence(OrdinateSet_XYZM, coordinates[1])
      ];

      // act
      WKTReader rdr = reader2D;
      MultiPoint mP2D =
          rdr.read("MULTIPOINT ((10 10), (20 20))")! as MultiPoint;
      MultiPoint mP2DE = rdr.read("MULTIPOINT EMPTY")! as MultiPoint;
      rdr = reader3D;
      MultiPoint mP3D =
          rdr.read("MULTIPOINT Z((10 10 10), (20 20 10))")! as MultiPoint;
      rdr = reader2DM;
      MultiPoint mP2DM =
          rdr.read("MULTIPOINT M((10 10 11), (20 20 11))")! as MultiPoint;
      rdr = reader3DM;
      MultiPoint mP3DM = rdr
          .read("MULTIPOINT ZM((10 10 10 11), (20 20 10 11))")! as MultiPoint;

      // assert
      expect(
          checkEqual(csMP2D[0],
              (mP2D.getGeometryN(0) as Point).getCoordinateSequence()),
          true);
      expect(
          checkEqual(csMP2D[1],
              (mP2D.getGeometryN(1) as Point).getCoordinateSequence()),
          true);
      expect(mP2DE.isEmpty(), true);
      expect(mP2DE.getNumGeometries() == 0, true);
      expect(
          checkEqual(csMP3D[0],
              (mP3D.getGeometryN(0) as Point).getCoordinateSequence()),
          true);
      expect(
          checkEqual(csMP3D[1],
              (mP3D.getGeometryN(1) as Point).getCoordinateSequence()),
          true);
      expect(
          checkEqual(csMP2DM[0],
              (mP2DM.getGeometryN(0) as Point).getCoordinateSequence()),
          true);
      expect(
          checkEqual(csMP2DM[1],
              (mP2DM.getGeometryN(1) as Point).getCoordinateSequence()),
          true);
      var mp3DMCS = (mP3DM.getGeometryN(0) as Point).getCoordinateSequence();
      expect(checkEqual(csMP3DM[0], mp3DMCS), true);
      expect(
          checkEqual(csMP3DM[1],
              (mP3DM.getGeometryN(1) as Point).getCoordinateSequence()),
          true);
    });
    test("testReadMultiLineString", () {
      // arrange
      List<List<double>> coordinates = [
        [10, 10, 20, 20],
        [15, 15, 30, 15]
      ];
      List<CoordinateSequence> csMls2D = [
        createSequence(OrdinateSet_XY, coordinates[0]),
        createSequence(OrdinateSet_XY, coordinates[1])
      ];
      List<CoordinateSequence> csMls3D = [
        createSequence(OrdinateSet_XYZ, coordinates[0]),
        createSequence(OrdinateSet_XYZ, coordinates[1])
      ];
      List<CoordinateSequence> csMls2DM = [
        createSequence(OrdinateSet_XYM, coordinates[0]),
        createSequence(OrdinateSet_XYM, coordinates[1])
      ];
      List<CoordinateSequence> csMls3DM = [
        createSequence(OrdinateSet_XYZM, coordinates[0]),
        createSequence(OrdinateSet_XYZM, coordinates[1])
      ];

      // act
      WKTReader rdr = reader2D;
      MultiLineString mLs2D =
          rdr.read("MULTILINESTRING ((10 10, 20 20), (15 15, 30 15))")!
              as MultiLineString;
      MultiLineString mLs2DE =
          rdr.read("MULTILINESTRING EMPTY")! as MultiLineString;
      rdr = reader3D;
      MultiLineString mLs3D = rdr.read(
              "MULTILINESTRING Z((10 10 10, 20 20 10), (15 15 10, 30 15 10))")!
          as MultiLineString;
      rdr = reader2DM;
      MultiLineString mLs2DM = rdr.read(
              "MULTILINESTRING M((10 10 11, 20 20 11), (15 15 11, 30 15 11))")!
          as MultiLineString;
      rdr = reader3DM;
      MultiLineString mLs3DM = rdr.read(
              "MULTILINESTRING ZM((10 10 10 11, 20 20 10 11), (15 15 10 11, 30 15 10 11))")!
          as MultiLineString;

      // assert
      expect(
          checkEqual(csMls2D[0],
              (mLs2D.getGeometryN(0) as LineString).getCoordinateSequence()),
          true);
      expect(
          checkEqual(csMls2D[1],
              (mLs2D.getGeometryN(1) as LineString).getCoordinateSequence()),
          true);
      expect(mLs2DE.isEmpty(), true);
      expect(mLs2DE.getNumGeometries() == 0, true);
      expect(
          checkEqual(csMls3D[0],
              (mLs3D.getGeometryN(0) as LineString).getCoordinateSequence()),
          true);
      expect(
          checkEqual(csMls3D[1],
              (mLs3D.getGeometryN(1) as LineString).getCoordinateSequence()),
          true);
      expect(
          checkEqual(csMls2DM[0],
              (mLs2DM.getGeometryN(0) as LineString).getCoordinateSequence()),
          true);
      expect(
          checkEqual(csMls2DM[1],
              (mLs2DM.getGeometryN(1) as LineString).getCoordinateSequence()),
          true);
      expect(
          checkEqual(csMls3DM[0],
              (mLs3DM.getGeometryN(0) as LineString).getCoordinateSequence()),
          true);
      expect(
          checkEqual(csMls3DM[1],
              (mLs3DM.getGeometryN(1) as LineString).getCoordinateSequence()),
          true);
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
        rdr.read("MULTIPOLYGON (((10 10, 10 20, 20 20, 20 15, 10 10)))")!
            as MultiPolygon,
        rdr.read(
                "MULTIPOLYGON (((10 10, 10 20, 20 20, 20 15, 10 10), (11 11, 12 11, 12 12, 12 11, 11 11)))")!
            as MultiPolygon,
        rdr.read(
                "MULTIPOLYGON (((10 10, 10 20, 20 20, 20 15, 10 10), (11 11, 12 11, 12 12, 12 11, 11 11)), ((60 60, 70 70, 80 60, 60 60)))")!
            as MultiPolygon
      ];
      MultiPolygon poly2DE = rdr.read("MULTIPOLYGON EMPTY")! as MultiPolygon;
      rdr = reader3D;
      List<MultiPolygon> poly3D = [
        rdr.read(
                "MULTIPOLYGON Z(((10 10 10, 10 20 10, 20 20 10, 20 15 10, 10 10 10)))")!
            as MultiPolygon,
        rdr.read(
                "MULTIPOLYGON Z(((10 10 10, 10 20 10, 20 20 10, 20 15 10, 10 10 10), (11 11 10, 12 11 10, 12 12 10, 12 11 10, 11 11 10)))")!
            as MultiPolygon,
        rdr.read(
                "MULTIPOLYGON Z(((10 10 10, 10 20 10, 20 20 10, 20 15 10, 10 10 10), (11 11 10, 12 11 10, 12 12 10, 12 11 10, 11 11 10)), ((60 60 10, 70 70 10, 80 60 10, 60 60 10)))")!
            as MultiPolygon
      ];
      List<MultiPolygon> poly2DM = [
        rdr.read(
                "MULTIPOLYGON M(((10 10 11, 10 20 11, 20 20 11, 20 15 11, 10 10 11)))")!
            as MultiPolygon,
        rdr.read(
                "MULTIPOLYGON M(((10 10 11, 10 20 11, 20 20 11, 20 15 11, 10 10 11), (11 11 11, 12 11 11, 12 12 11, 12 11 11, 11 11 11)))")!
            as MultiPolygon,
        rdr.read(
                "MULTIPOLYGON M(((10 10 11, 10 20 11, 20 20 11, 20 15 11, 10 10 11), (11 11 11, 12 11 11, 12 12 11, 12 11 11, 11 11 11)), ((60 60 11, 70 70 11, 80 60 11, 60 60 11)))")!
            as MultiPolygon
      ];
      rdr = reader3DM;
      List<MultiPolygon> poly3DM = [
        rdr.read(
                "MULTIPOLYGON ZM(((10 10 10 11, 10 20 10 11, 20 20 10 11, 20 15 10 11, 10 10 10 11)))")!
            as MultiPolygon,
        rdr.read(
                "MULTIPOLYGON ZM(((10 10 10 11, 10 20 10 11, 20 20 10 11, 20 15 10 11, 10 10 10 11), (11 11 10 11, 12 11 10 11, 12 12 10 11, 12 11 10 11, 11 11 10 11)))")!
            as MultiPolygon,
        rdr.read(
                "MULTIPOLYGON ZM(((10 10 10 11, 10 20 10 11, 20 20 10 11, 20 15 10 11, 10 10 10 11), (11 11 10 11, 12 11 10 11, 12 12 10 11, 12 11 10 11, 11 11 10 11)), ((60 60 10 11, 70 70 10 11, 80 60 10 11, 60 60 10 11)))")!
            as MultiPolygon
      ];

      // assert
      expect(
          checkEqual(
              csPoly2D[0],
              (poly2D[2].getGeometryN(0) as Polygon)
                  .getExteriorRing()
                  .getCoordinateSequence()),
          true);
      expect(
          checkEqual(
              csPoly2D[1],
              (poly2D[2].getGeometryN(0) as Polygon)
                  .getInteriorRingN(0)
                  .getCoordinateSequence()),
          true);
      expect(
          checkEqual(
              csPoly2D[2],
              (poly2D[2].getGeometryN(1) as Polygon)
                  .getExteriorRing()
                  .getCoordinateSequence()),
          true);
      expect(poly2DE.isEmpty(), true);
      expect(poly2DE.getNumGeometries() == 0, true);

      expect(
          checkEqual(
              csPoly3D[0],
              (poly3D[2].getGeometryN(0) as Polygon)
                  .getExteriorRing()
                  .getCoordinateSequence()),
          true);
      expect(
          checkEqual(
              csPoly3D[1],
              (poly3D[2].getGeometryN(0) as Polygon)
                  .getInteriorRingN(0)
                  .getCoordinateSequence()),
          true);
      expect(
          checkEqual(
              csPoly3D[2],
              (poly3D[2].getGeometryN(1) as Polygon)
                  .getExteriorRing()
                  .getCoordinateSequence()),
          true);

      var poly2DMExtRing =
          (poly2DM[2].getGeometryN(0) as Polygon).getExteriorRing();
      var poly2DMExtRingCS = poly2DMExtRing.getCoordinateSequence();
      expect(checkEqual(csPoly2DM[0], poly2DMExtRingCS), true);
      expect(
          checkEqual(
              csPoly2DM[1],
              (poly2DM[2].getGeometryN(0) as Polygon)
                  .getInteriorRingN(0)
                  .getCoordinateSequence()),
          true);
      expect(
          checkEqual(
              csPoly2DM[2],
              (poly2DM[2].getGeometryN(1) as Polygon)
                  .getExteriorRing()
                  .getCoordinateSequence()),
          true);

      expect(
          checkEqual(
              csPoly3DM[0],
              (poly3DM[2].getGeometryN(0) as Polygon)
                  .getExteriorRing()
                  .getCoordinateSequence()),
          true);
      expect(
          checkEqual(
              csPoly3DM[1],
              (poly3DM[2].getGeometryN(0) as Polygon)
                  .getInteriorRingN(0)
                  .getCoordinateSequence()),
          true);
      expect(
          checkEqual(
              csPoly3DM[2],
              (poly3DM[2].getGeometryN(1) as Polygon)
                  .getExteriorRing()
                  .getCoordinateSequence()),
          true);
    });
    test("", () {});
    test("", () {});
    test("", () {});
  });
}

CoordinateSequence createSequence(
    List<Ordinate> ordinateFlags, List<double> xy) {
// get the number of dimension to verify size of provided ordinate values array
  int dimension = requiredDimension(ordinateFlags);

// inject additional values
  List<double> ordinateValues = injectZM(ordinateFlags, xy);

  if ((ordinateValues.length % dimension) != 0)
    throw new ArgumentError(
        "ordinateFlags and number of provided ordinate values don't match");

// get the required size of the sequence
  int size = ordinateValues.length ~/ dimension;

// create a sequence capable of storing all ordinate values.
  CoordinateSequence res = getCSFactory(ordinateFlags)
      .createSizeDim(size, requiredDimension(ordinateFlags));

// fill in values
  int k = 0;
  for (int i = 0; i < ordinateValues.length; i += dimension) {
    for (int j = 0; j < dimension; j++)
      res.setOrdinate(k, j, ordinateValues[i + j]);
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
  List<double> res = List.filled(size * dimension, 0.0);
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
WKTReader getWKTReader(
    List<Ordinate> ordinateFlags, PrecisionModel precisionModel) {
  WKTReader result;

  if (!ordinateFlags.contains(Ordinate.X)) ordinateFlags.add(Ordinate.X);
  if (!ordinateFlags.contains(Ordinate.Y)) ordinateFlags.add(Ordinate.Y);

  if (ordinateFlags.length == 2) {
    result = WKTReader.withFactory(
        GeometryFactory(precisionModel, 0, CoordinateArraySequenceFactory()));
    result.setIsOldJtsCoordinateSyntaxAllowed(false);
  } else if (ordinateFlags.contains(Ordinate.Z)) {
    result = WKTReader.withFactory(
        GeometryFactory(precisionModel, 0, CoordinateArraySequenceFactory()));
  } else if (ordinateFlags.contains(Ordinate.M)) {
    result = WKTReader.withFactory(GeometryFactory(
        precisionModel, 0, PackedCoordinateSequenceFactory.DOUBLE_FACTORY));
    result.setIsOldJtsCoordinateSyntaxAllowed(false);
  } else {
    result = WKTReader.withFactory(GeometryFactory(
        precisionModel, 0, PackedCoordinateSequenceFactory.DOUBLE_FACTORY));
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
WKTReader getWKTReaderFromOrdinateSetAndScale(
    List<Ordinate> ordinateFlags, double scale) {
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
bool checkEqualWithTolerance(
    CoordinateSequence seq1, CoordinateSequence seq2, double tolerance) {
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
bool checkEqualWithDim(
    CoordinateSequence seq1, CoordinateSequence seq2, int dimension) {
  return checkEqualWithDimTol(seq1, seq2, dimension, 0.0);
}

/// Checks two {@link CoordinateSequence}s for equality. The following items are checked:
/// <ul>
///   <li>size</li><li>dimension up to {@code dimension}</li><li>ordinate values with {@code tolerance}</li>
/// </ul>
/// @param seq1 a sequence
/// @param seq2 another sequence
/// @return {@code true} if both sequences are equal
bool checkEqualWithDimTol(CoordinateSequence seq1, CoordinateSequence seq2,
    int dimension, double tolerance) {
  if (seq1 != null && seq2 == null) return false;
  if (seq1 == null && seq2 != null) return false;

  if (seq1.size() != seq2.size()) return false;

  if (seq1.getDimension() < dimension)
    throw ArgumentError("dimension too high for seq1");
  if (seq2.getDimension() < dimension)
    throw ArgumentError("dimension too high for seq2");

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
    runWKBTestWKT(
        "POLYGON ((0 0, 100 0, 100 100, 0 100, 0 0), (1 1, 1 10, 10 10, 10 1, 1 1) )");
  }

  void testMultiPoint() {
    runWKBTestWKT("MULTIPOINT ((0 0), (1 4), (100 200))");
  }

  void testMultiLineString() {
    runWKBTestWKT(
        "MULTILINESTRING ((0 0, 1 10), (10 10, 20 30), (123 123, 456 789))");
  }

  void testMultiPolygon() {
    runWKBTestWKT(
        "MULTIPOLYGON ( ((0 0, 100 0, 100 100, 0 100, 0 0), (1 1, 1 10, 10 10, 10 1, 1 1) ), ((200 200, 200 250, 250 250, 250 200, 200 200)) )");
  }

  void testGeometryCollection() {
    runWKBTestWKT(
        "GEOMETRYCOLLECTION ( POINT ( 1 1), LINESTRING (0 0, 10 10), POLYGON ((0 0, 100 0, 100 100, 0 100, 0 0)) )");
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
        new GeometryFactory.withCoordinateSequenceFactory(
            new PackedCoordinateSequenceFactory.withType(
                PackedCoordinateSequenceFactory.DOUBLE));
    WKTReader rdr = new WKTReader.withFactory(geomFactory);
    Geometry g = rdr.read(wkt)!;

// Since we are using a PCS of dim=2, only check 2-dimensional storage
    runWKBTest(g, 2, true);
    runWKBTest(g, 2, false);
  }

  void runWKBTestCoordinateArray(String wkt) {
    GeometryFactory geomFactory = new GeometryFactory.defaultPrecision();
    WKTReader rdr = new WKTReader.withFactory(geomFactory);
    Geometry g = rdr.read(wkt)!;

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

  void runWKBTestWithBO(
      Geometry g, int dimension, Endian byteOrder, bool toHex) {
    runGeometry(g, dimension, byteOrder, toHex, 100);
    runGeometry(g, dimension, byteOrder, toHex, 0);
    runGeometry(g, dimension, byteOrder, toHex, 101010);
    runGeometry(g, dimension, byteOrder, toHex, -1);
  }

  void setZ(Geometry g) {
    g.applyCF(Cf());
  }

//static Comparator comp2D = new Coordinate.DimensionalComparator();
//static Comparator comp3D = new Coordinate.DimensionalComparator(3);

  static Comparator<CoordinateSequence> comp2 =
      CoordinateSequenceComparatorBuilder.withLimit(2);
  static Comparator<CoordinateSequence> comp3 =
      CoordinateSequenceComparatorBuilder.withLimit(3);

  /**
   * Use single WKB reader, to ensure it can be used for multiple input geometries
   */
  WKBReader wkbReader = new WKBReader.withFactory(geomFactory);

  void runGeometry(
      Geometry g, int dimension, Endian byteOrder, bool toHex, int srid) {
    bool includeSRID = false;
    if (srid >= 0) {
      includeSRID = true;
      g.setSRID(srid);
    }

    WKBWriter wkbWriter =
        new WKBWriter.withDimOrderSrid(dimension, byteOrder, includeSRID);
    List<int> wkb = wkbWriter.write(g);
    Uint8List wkbU8 = Uint8List.fromList(wkb);

    String? wkbHex = null;
    if (toHex) wkbHex = WKBWriter.toHex(wkbU8);

    if (toHex) wkb = WKBReader.hexToBytes(wkbHex!);
    Geometry g2 = wkbReader.read(wkbU8);

    Comparator<CoordinateSequence> comp = (dimension == 2) ? comp2 : comp3;
    bool isEqual = (g.compareToWithComparator(g2, comp) == 0);
    assertTrue(isEqual);

    if (includeSRID) {
      bool isSRIDEqual = g.getSRID() == g2.getSRID();
      assertTrue(isSRIDEqual);
    }
  }
}

class Cf implements CoordinateFilter {
  @override
  void filter(Coordinate? coord) {
    coord!.setZ((coord.x + coord.y) / 2);
  }
}

/**
 * Tests for reading WKB.
 *
 * @author Martin Davis
 *
 */
class WKBReaderTest {
  static GeometryFactory geomFactory = new GeometryFactory.defaultPrecision();
  WKTReader rdr = new WKTReader.withFactory(geomFactory);
  WKTReader rdrM = new WKTReader.withFactory(
      new GeometryFactory.withCoordinateSequenceFactory(
          PackedCoordinateSequenceFactory.DOUBLE_FACTORY));

  void testShortPolygons() {
    // one point
    checkWKBGeometry(
        "0000000003000000010000000140590000000000004069000000000000",
        "POLYGON ((100 200, 100 200, 100 200, 100 200))");
    // two point
    checkWKBGeometry(
        "000000000300000001000000024059000000000000406900000000000040590000000000004069000000000000",
        "POLYGON ((100 200, 100 200, 100 200, 100 200))");
  }

  void testSinglePointLineString() {
    checkWKBGeometry("00000000020000000140590000000000004069000000000000",
        "LINESTRING (100 200, 100 200)");
  }

  /**
   * After removing the 39 bytes of MBR info at the front, and the
   * end-of-geometry byte, * Spatialite native BLOB is very similar
   * to WKB, except instead of a endian marker at the start of each
   * geometry in a multi-geometry, it has a start marker of 0x69.
   * Endianness is determined by the endian value of the multigeometry.
   *
   * @
   */
  void testSpatialiteMultiGeometry() {
//multipolygon
    checkWKBGeometry(
        "01060000000200000069030000000100000004000000000000000000444000000000000044400000000000003440000000000080464000000000008046400000000000003E4000000000000044400000000000004440690300000001000000040000000000000000003E40000000000000344000000000000034400000000000002E40000000000000344000000000000039400000000000003E400000000000003440",
        "MULTIPOLYGON (((40 40, 20 45, 45 30, 40 40)), ((30 20, 20 15, 20 25, 30 20)))'");

//multipoint
    checkWKBGeometry(
        "0104000000020000006901000000000000000000F03F000000000000F03F690100000000000000000000400000000000000040",
        "MULTIPOINT(1 1,2 2)'");

//multiline
    checkWKBGeometry(
        "010500000002000000690200000003000000000000000000244000000000000024400000000000003440000000000000344000000000000024400000000000004440690200000004000000000000000000444000000000000044400000000000003E400000000000003E40000000000000444000000000000034400000000000003E400000000000002440",
        "MULTILINESTRING ((10 10, 20 20, 10 40), (40 40, 30 30, 40 20, 30 10))");

//geometrycollection
    checkWKBGeometry(
        "010700000002000000690100000000000000000010400000000000001840690200000002000000000000000000104000000000000018400000000000001C400000000000002440",
        "GEOMETRYCOLLECTION(POINT(4 6),LINESTRING(4 6,7 10))");
  }

  void test2dSpatialiteWKB() {
// Point
    checkWKBGeometry(
        "0101000020E6100000000000000000F03F0000000000000040", "POINT(1 2)");
// LineString
    checkWKBGeometry(
        "0102000020E610000002000000000000000000F03F000000000000004000000000000008400000000000001040",
        "LINESTRING(1 2, 3 4)");
// Polygon
    checkWKBGeometry(
        "0103000020E61000000200000005000000000000000000000000000000000000000000000000000000000000000000244000000000000024400000000000002440000000000000244000000000000000000000000000000000000000000000000005000000000000000000F03F000000000000F03F000000000000F03F0000000000002240000000000000224000000000000022400000000000002240000000000000F03F000000000000F03F000000000000F03F",
        "POLYGON((0 0,0 10,10 10,10 0,0 0),(1 1,1 9,9 9,9 1,1 1))");
// MultiPoint
    checkWKBGeometry(
        "0104000020E61000000200000001010000000000000000000000000000000000F03F010100000000000000000000400000000000000840",
        "MULTIPOINT(0 1,2 3)");
// MultiLineString
    checkWKBGeometry(
        "0105000020E6100000020000000102000000020000000000000000000000000000000000F03F000000000000004000000000000008400102000000020000000000000000001040000000000000144000000000000018400000000000001C40",
        "MULTILINESTRING((0 1,2 3),(4 5,6 7))");
// MultiPolygon
    checkWKBGeometry(
        "0106000020E61000000200000001030000000200000005000000000000000000000000000000000000000000000000000000000000000000244000000000000024400000000000002440000000000000244000000000000000000000000000000000000000000000000005000000000000000000F03F000000000000F03F000000000000F03F0000000000002240000000000000224000000000000022400000000000002240000000000000F03F000000000000F03F000000000000F03F0103000000010000000500000000000000000022C0000000000000000000000000000022C00000000000002440000000000000F0BF0000000000002440000000000000F0BF000000000000000000000000000022C00000000000000000",
        "MULTIPOLYGON(((0 0,0 10,10 10,10 0,0 0),(1 1,1 9,9 9,9 1,1 1)),((-9 0,-9 10,-1 10,-1 0,-9 0)))");
// GeometryCollection
    checkWKBGeometry(
        "0107000020E61000000900000001010000000000000000000000000000000000F03F01010000000000000000000000000000000000F03F01010000000000000000000040000000000000084001020000000200000000000000000000400000000000000840000000000000104000000000000014400102000000020000000000000000000000000000000000F03F000000000000004000000000000008400102000000020000000000000000001040000000000000144000000000000018400000000000001C4001030000000200000005000000000000000000000000000000000000000000000000000000000000000000244000000000000024400000000000002440000000000000244000000000000000000000000000000000000000000000000005000000000000000000F03F000000000000F03F000000000000F03F0000000000002240000000000000224000000000000022400000000000002240000000000000F03F000000000000F03F000000000000F03F01030000000200000005000000000000000000000000000000000000000000000000000000000000000000244000000000000024400000000000002440000000000000244000000000000000000000000000000000000000000000000005000000000000000000F03F000000000000F03F000000000000F03F0000000000002240000000000000224000000000000022400000000000002240000000000000F03F000000000000F03F000000000000F03F0103000000010000000500000000000000000022C0000000000000000000000000000022C00000000000002440000000000000F0BF0000000000002440000000000000F0BF000000000000000000000000000022C00000000000000000",
        "GEOMETRYCOLLECTION(POINT(0 1),POINT(0 1),POINT(2 3),LINESTRING(2 3,4 5),LINESTRING(0 1,2 3),LINESTRING(4 5,6 7),POLYGON((0 0,0 10,10 10,10 0,0 0),(1 1,1 9,9 9,9 1,1 1)),POLYGON((0 0,0 10,10 10,10 0,0 0),(1 1,1 9,9 9,9 1,1 1)),POLYGON((-9 0,-9 10,-1 10,-1 0,-9 0)))");
  }

  void testSpatialiteWKB_Z() {
// PointZ
    checkWKBGeometry(
        "01010000A0E6100000000000000000F03F00000000000000400000000000000840",
        "POINT Z(1 2 3)");
// LineStringZ
    checkWKBGeometry(
        "01020000A0E610000002000000000000000000F03F00000000000000400000000000000840000000000000104000000000000014400000000000001840",
        "LINESTRING Z(1 2 3, 4 5 6)");
// PolygonZ
    checkWKBGeometry(
        "01030000A0E6100000020000000500000000000000000000000000000000000000000000000000594000000000000000000000000000002440000000000000594000000000000024400000000000002440000000000000594000000000000024400000000000000000000000000000594000000000000000000000000000000000000000000000594005000000000000000000F03F000000000000F03F0000000000005940000000000000F03F000000000000224000000000000059400000000000002240000000000000224000000000000059400000000000002240000000000000F03F0000000000005940000000000000F03F000000000000F03F0000000000005940",
        "POLYGON Z((0 0 100,0 10 100,10 10 100,10 0 100,0 0 100),(1 1 100,1 9 100,9 9 100,9 1 100,1 1 100))");
// MultiPointZ
    checkWKBGeometry(
        "01040000A0E61000000200000001010000800000000000000000000000000000F03F00000000000000400101000080000000000000084000000000000010400000000000001440",
        "MULTIPOINTS Z(0 1 2, 3 4 5)");
// MultiLineStringZ
    checkWKBGeometry(
        "01050000A0E6100000020000000102000080020000000000000000000000000000000000F03F000000000000004000000000000008400000000000001040000000000000144001020000800200000000000000000018400000000000001C400000000000002040000000000000224000000000000024400000000000002640",
        "MULTILINESTRING Z((0 1 2,3 4 5),(6 7 8,9 10 11))");
// MultiPolygonZ
    checkWKBGeometry(
        "01060000A0E6100000020000000103000080020000000500000000000000000000000000000000000000000000000000594000000000000000000000000000002440000000000000594000000000000024400000000000002440000000000000594000000000000024400000000000000000000000000000594000000000000000000000000000000000000000000000594005000000000000000000F03F000000000000F03F0000000000005940000000000000F03F000000000000224000000000000059400000000000002240000000000000224000000000000059400000000000002240000000000000F03F0000000000005940000000000000F03F000000000000F03F00000000000059400103000080010000000500000000000000000022C00000000000000000000000000000494000000000000022C000000000000024400000000000004940000000000000F0BF00000000000024400000000000004940000000000000F0BF0000000000000000000000000000494000000000000022C000000000000000000000000000004940",
        "MULTIPOLYGON Z(((0 0 100,0 10 100,10 10 100,10 0 100,0 0 100),(1 1 100,1 9 100,9 9 100,9 1 100,1 1 100)),((-9 0 50,-9 10 50,-1 10 50,-1 0 50,-9 0 50)))");
// GeometryCollectionZ
  }

  void testSpatialiteWKB_M() {
// PointM
    checkWKBGeometry(
        "0101000060E6100000000000000000F03F00000000000000400000000000000840",
        "POINT M(1 2 3)");
// LineStringM
    checkWKBGeometry(
        "0102000060E610000002000000000000000000F03F00000000000000400000000000000840000000000000104000000000000014400000000000001840",
        "LINESTRING M(1 2 3,4 5 6)");
// PolygonM
    checkWKBGeometry(
        "0103000060E6100000020000000500000000000000000000000000000000000000000000000000594000000000000000000000000000002440000000000000594000000000000024400000000000002440000000000000594000000000000024400000000000000000000000000000594000000000000000000000000000000000000000000000594005000000000000000000F03F000000000000F03F0000000000005940000000000000F03F000000000000224000000000000059400000000000002240000000000000224000000000000059400000000000002240000000000000F03F0000000000005940000000000000F03F000000000000F03F0000000000005940",
        "POLYGON M((0 0 100,0 10 100,10 10 100,10 0 100,0 0 100),(1 1 100,1 9 100,9 9 100,9 1 100,1 1 100))");
// MultiPointM
    checkWKBGeometry(
        "01040000A0E61000000200000001010000800000000000000000000000000000F03F00000000000000400101000080000000000000084000000000000010400000000000001440",
        "MULTIPOINT M(0 1 2,3 4 5)");
// MultiLineStringM
    checkWKBGeometry(
        "0105000060E6100000020000000102000040020000000000000000000000000000000000F03F000000000000004000000000000008400000000000001040000000000000144001020000400200000000000000000018400000000000001C400000000000002040000000000000224000000000000024400000000000002640",
        "MULTILINESTRING M((0 1 2,3 4 5),(6 7 8,9 10 11))");
// MultiPolygonM
    checkWKBGeometry(
        "0106000060E6100000020000000103000040020000000500000000000000000000000000000000000000000000000000594000000000000000000000000000002440000000000000594000000000000024400000000000002440000000000000594000000000000024400000000000000000000000000000594000000000000000000000000000000000000000000000594005000000000000000000F03F000000000000F03F0000000000005940000000000000F03F000000000000224000000000000059400000000000002240000000000000224000000000000059400000000000002240000000000000F03F0000000000005940000000000000F03F000000000000F03F00000000000059400103000040010000000500000000000000000022C00000000000000000000000000000494000000000000022C000000000000024400000000000004940000000000000F0BF00000000000024400000000000004940000000000000F0BF0000000000000000000000000000494000000000000022C000000000000000000000000000004940",
        "MULTIPOLYGON M(((0 0 100,0 10 100,10 10 100,10 0 100,0 0 100),(1 1 100,1 9 100,9 9 100,9 1 100,1 1 100)),((-9 0 50,-9 10 50,-1 10 50,-1 0 50,-9 0 50)))");
// GeometryCollectionM
//checkWKBGeometry("0107000020E61000000900000001010000000000000000000000000000000000F03F01010000000000000000000000000000000000F03F01010000000000000000000040000000000000084001020000000200000000000000000000400000000000000840000000000000104000000000000014400102000000020000000000000000000000000000000000F03F000000000000004000000000000008400102000000020000000000000000001040000000000000144000000000000018400000000000001C4001030000000200000005000000000000000000000000000000000000000000000000000000000000000000244000000000000024400000000000002440000000000000244000000000000000000000000000000000000000000000000005000000000000000000F03F000000000000F03F000000000000F03F0000000000002240000000000000224000000000000022400000000000002240000000000000F03F000000000000F03F000000000000F03F01030000000200000005000000000000000000000000000000000000000000000000000000000000000000244000000000000024400000000000002440000000000000244000000000000000000000000000000000000000000000000005000000000000000000F03F000000000000F03F000000000000F03F0000000000002240000000000000224000000000000022400000000000002240000000000000F03F000000000000F03F000000000000F03F0103000000010000000500000000000000000022C0000000000000000000000000000022C00000000000002440000000000000F0BF0000000000002440000000000000F0BF000000000000000000000000000022C00000000000000000",
//    "MULTIPOLYGONM(((0 0 100,0 10 100,10 10 100,10 0 100,0 0 100),(1 1 100,1 9 100,9 9 100,9 1 100,1 1 100)),((-9 0 50,-9 10 50,-1 10 50,-1 0 50,-9 0 50)))");
  }

  void testSpatialiteWKB_ZM() {
// PointZM
    checkWKBGeometry(
        "01010000E0E6100000000000000000F03F000000000000004000000000000008400000000000006940",
        "POINT ZM (1 2 3 200)");
// LineStringZM
    checkWKBGeometry(
        "01020000E0E610000002000000000000000000F03F0000000000000040000000000000084000000000000069400000000000001040000000000000144000000000000018400000000000006940",
        "LINESTRING ZM (1 2 3 200,4 5 6 200)");
// PolygonZM
    checkWKBGeometry(
        "01030000E0E610000002000000050000000000000000000000000000000000000000000000000059400000000000006940000000000000000000000000000024400000000000005940000000000000694000000000000024400000000000002440000000000000594000000000000069400000000000002440000000000000000000000000000059400000000000006940000000000000000000000000000000000000000000005940000000000000694005000000000000000000F03F000000000000F03F00000000000059400000000000006940000000000000F03F00000000000022400000000000005940000000000000694000000000000022400000000000002240000000000000594000000000000069400000000000002240000000000000F03F00000000000059400000000000006940000000000000F03F000000000000F03F00000000000059400000000000006940",
        "POLYGON ZM ((0 0 100 200,0 10 100 200,10 10 100 200,10 0 100 200,0 0 100 200),(1 1 100 200,1 9 100 200,9 9 100 200,9 1 100 200,1 1 100 200))");
// MultiPointZM
    checkWKBGeometry(
        "01040000E0E61000000200000001010000C00000000000000000000000000000F03F0000000000000040000000000000694001010000C00000000000000840000000000000104000000000000014400000000000006940",
        "MULTIPOINT ZM (0 1 2 200,3 4 5 200)");
// MultiLineStringZM
    checkWKBGeometry(
        "01050000E0E61000000200000001020000C0020000000000000000000000000000000000F03F00000000000000400000000000006940000000000000084000000000000010400000000000001440000000000000694001020000C00200000000000000000018400000000000001C40000000000000204000000000000069400000000000002240000000000000244000000000000026400000000000006940",
        "MULTILINESTRING ZM ((0 1 2 200,3 4 5 200),(6 7 8 200,9 10 11 200))");
// MultiPolygonZM
    checkWKBGeometry(
        "01060000E0E61000000200000001030000C002000000050000000000000000000000000000000000000000000000000059400000000000006940000000000000000000000000000024400000000000005940000000000000694000000000000024400000000000002440000000000000594000000000000069400000000000002440000000000000000000000000000059400000000000006940000000000000000000000000000000000000000000005940000000000000694005000000000000000000F03F000000000000F03F00000000000059400000000000006940000000000000F03F00000000000022400000000000005940000000000000694000000000000022400000000000002240000000000000594000000000000069400000000000002240000000000000F03F00000000000059400000000000006940000000000000F03F000000000000F03F0000000000005940000000000000694001030000C0010000000500000000000000000022C000000000000000000000000000004940000000000000694000000000000022C0000000000000244000000000000049400000000000006940000000000000F0BF000000000000244000000000000049400000000000006940000000000000F0BF00000000000000000000000000004940000000000000694000000000000022C0000000000000000000000000000049400000000000006940",
        "MULTIPOLYGON ZM (((0 0 100 200,0 10 100 200,10 10 100 200,10 0 100 200,0 0 100 200),(1 1 100 200,1 9 100 200,9 9 100 200,9 1 100 200,1 1 100 200)),((-9 0 50 200,-9 10 50 200,-1 10 50 200,-1 0 50 200,-9 0 50 200)))");
// GeometryCollectionZM
//checkWKBGeometry("0107000020E61000000900000001010000000000000000000000000000000000F03F01010000000000000000000000000000000000F03F01010000000000000000000040000000000000084001020000000200000000000000000000400000000000000840000000000000104000000000000014400102000000020000000000000000000000000000000000F03F000000000000004000000000000008400102000000020000000000000000001040000000000000144000000000000018400000000000001C4001030000000200000005000000000000000000000000000000000000000000000000000000000000000000244000000000000024400000000000002440000000000000244000000000000000000000000000000000000000000000000005000000000000000000F03F000000000000F03F000000000000F03F0000000000002240000000000000224000000000000022400000000000002240000000000000F03F000000000000F03F000000000000F03F01030000000200000005000000000000000000000000000000000000000000000000000000000000000000244000000000000024400000000000002440000000000000244000000000000000000000000000000000000000000000000005000000000000000000F03F000000000000F03F000000000000F03F0000000000002240000000000000224000000000000022400000000000002240000000000000F03F000000000000F03F000000000000F03F0103000000010000000500000000000000000022C0000000000000000000000000000022C00000000000002440000000000000F0BF0000000000002440000000000000F0BF000000000000000000000000000022C00000000000000000",
//    "MULTIPOLYGONM(((0 0 100,0 10 100,10 10 100,10 0 100,0 0 100),(1 1 100,1 9 100,9 9 100,9 1 100,1 1 100)),((-9 0 50,-9 10 50,-1 10 50,-1 0 50,-9 0 50)))");
  }

//  /**
//   * Not yet implemented satisfactorily.
//   *
//   * @
//   */
//  void XXtestIllFormedWKB() {
//// WKB is missing LinearRing entry
//    checkWKBGeometry("00000000030000000140590000000000004069000000000000", "POLYGON ((100 200, 100 200, 100 200, 100 200)");
//  }

  static Comparator<CoordinateSequence> comp2 =
      CoordinateSequenceComparatorBuilder.withLimit(2);

  void checkWKBGeometry(String wkbHex, String expectedWKT) {
    WKBReader wkbReader = new WKBReader.withFactory(geomFactory);
    List<int> wkb = WKBReader.hexToBytes(wkbHex);
    Geometry g2 = wkbReader.read(wkb);

    WKTReader useRdr = rdr;
    if (expectedWKT.contains("ZM"))
      useRdr = rdrM;
    else if (expectedWKT.contains("M(") || expectedWKT.contains("M ("))
      useRdr = rdrM;

    Geometry expected = useRdr.read(expectedWKT)!;

    bool isEqual = (expected.compareToWithComparator(g2, comp2) == 0);
    if (!isEqual) {
      print(g2);
      print(expected);
    }
    assertTrue(isEqual);
  }
}
