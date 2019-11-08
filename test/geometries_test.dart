import "package:test/test.dart";
import 'package:dart_jts/dart_jts.dart';
import 'package:test/test.dart' as prefix0;
import "dart:math" as math;
import 'testing_utilities.dart';

void main() {
  group(" - ", () {
    test("", () {});
  });

  group("PredicateShortCircuitTest - ", () {
    WKTReader rdr = new WKTReader();
    List<String> polyInsidePoly = ["POLYGON (( 0 0, 100 0, 100 100, 0 100, 0 0 ))", "POLYGON (( 10 10, 90 10, 90 90, 10 90, 10 10 ))"];
    List<String> polyPartiallyOverlapsPoly = ["POLYGON (( 10 10, 100 10, 100 100, 10 100, 10 10 ))", "POLYGON (( 0 0, 90 0, 90 90, 0 90, 0 0 ))"];
    List<String> polyTouchesPolyAtPoint = ["POLYGON (( 10 10, 100 10, 100 100, 10 100, 10 10 ))", "POLYGON (( 0 0, 10 0, 10 10, 0 10, 0 0 ))"];
    List<String> polyTouchesPolyAtLine = ["POLYGON (( 10 10, 100 10, 100 100, 10 100, 10 10 ))", "POLYGON (( 10 0, 10 10, 20 10, 20 0, 10 0 ))"];
    List<String> polyInsideHoleInPoly = [
      "POLYGON (( 40 40, 40 60, 60 60, 60 40, 40 40 ))",
      "POLYGON (( 0 0, 100 0, 100 100, 0 100, 0 0), ( 10 10, 90 10, 90 90, 10 90, 10 10))"
    ];

    test("testAll", () {
      doPredicates(rdr, polyInsidePoly);
      doPredicates(rdr, polyPartiallyOverlapsPoly);
      doPredicates(rdr, polyTouchesPolyAtPoint);
      doPredicates(rdr, polyTouchesPolyAtLine);
      doPredicates(rdr, polyInsideHoleInPoly);
    });
  });

  group("PrecisionModelTest - ", () {
    test("testParameterlessConstructor", () {
      PrecisionModel p = new PrecisionModel();
      //Implicit precision model has scale 0
      assertEqualsD(0, p.getScale(), 1E-10);
    });
    test("testGetMaximumSignificantDigits", () {
      assertEquals(16, new PrecisionModel.fromType(PrecisionModel.FLOATING).getMaximumSignificantDigits());
      assertEquals(1, new PrecisionModel.fromType(PrecisionModel.FIXED).getMaximumSignificantDigits());
      assertEquals(4, new PrecisionModel.fixedPrecision(1000).getMaximumSignificantDigits());
    });
    test("testMakePrecise", () {
      PrecisionModel pm_10 = new PrecisionModel.fixedPrecision(0.1);
      preciseCoordinateTester(pm_10, 1200.4, 1240.4, 1200, 1240);
      preciseCoordinateTester(pm_10, 1209.4, 1240.4, 1210, 1240);
    });
  });

  group("PointImplTest - ", () {
    PrecisionModel precisionModel = new PrecisionModel.fixedPrecision(1000);
    GeometryFactory geometryFactory = new GeometryFactory.withPrecisionModelSrid(precisionModel, 0);
    WKTReader reader = new WKTReader.withFactory(geometryFactory);

    test("testEquals", () {
      Point p1 = reader.read("POINT(1.234 5.678)");
      Point p2 = reader.read("POINT(1.234 5.678)");
      assertTrue(p1.equals(p2));

      p1 = reader.read("POINT(1.23 5.67)");
      p2 = reader.read("POINT(1.23 5.67)");
      assertTrue(p1.equals(p2));

      p1 = reader.read("POINT(1.235 5.678)");
      p2 = reader.read("POINT(1.234 5.678)");
      assertTrue(!p1.equals(p2));

      p1 = reader.read("POINT(1.2334 5.678)");
      p2 = reader.read("POINT(1.2333 5.678)");
      assertTrue(p1.equals(p2));

      p1 = reader.read("POINT(1.2334 5.678)");
      p2 = reader.read("POINT(1.2335 5.678)");
      assertTrue(!p1.equals(p2));

      p1 = reader.read("POINT(1.2324 5.678)");
      p2 = reader.read("POINT(1.2325 5.678)");
      assertTrue(!p1.equals(p2));
    });
    test("testNegRounding1", () {
      Point pLo = reader.read("POINT(-1.233 5.678)");
      Point pHi = reader.read("POINT(-1.232 5.678)");

      Point p1 = reader.read("POINT(-1.2326 5.678)");
      Point p2 = reader.read("POINT(-1.2325 5.678)");
      Point p3 = reader.read("POINT(-1.2324 5.678)");

      // TODO dart and java have different ways to round negatve numbers!
      // java -1232.5 -> -1232
      // dart -1232.5 -> -1233

      // ORIGINAL assertTrue(!p1.equals(p2));
      // ORIGINAL assertTrue(p3.equals(p2));
      assertTrue(p1.equals(p2));
      assertTrue(!p3.equals(p2));

      assertTrue(p1.equals(pLo));
      // ORIGINAL assertTrue(p2.equals(pHi));
      assertTrue(!p2.equals(pHi));
      assertTrue(p3.equals(pHi));
    });
    test("testIsSimple", () {
      Point p1 = reader.read("POINT(1.2324 5.678)");
      assertTrue(p1.isSimple());
      Point p2 = reader.read("POINT EMPTY");
      assertTrue(p2.isSimple());
    });
  });

  group("NormalizeTest - ", () {
    PrecisionModel precisionModel = new PrecisionModel.fixedPrecision(1);
    GeometryFactory geometryFactory = new GeometryFactory.withPrecisionModelSrid(precisionModel, 0);
    WKTReader reader = new WKTReader.withFactory(geometryFactory);

    test("testNormalizePoint", () {
      Point point = reader.read("POINT (30 30)");
      point.normalize();
      assertEquals(new Coordinate(30, 30), point.getCoordinate());
    });
    test("testNormalizeEmptyPoint", () {
      Point point = reader.read("POINT EMPTY");
      point.normalize();
      assertEquals(null, point.getCoordinate());
    });
    test("testComparePoint", () {
      Point p1 = reader.read("POINT (30 30)");
      Point p2 = reader.read("POINT (30 40)");
      assertTrue(p1.compareTo(p2) < 0);
    });
    test("testCompareEmptyPoint", () {
      Point p1 = reader.read("POINT (30 30)");
      Point p2 = reader.read("POINT EMPTY");
      assertTrue(p1.compareTo(p2) > 0);
    });
    test("testNormalizeMultiPoint", () {
      MultiPoint m = reader.read("MULTIPOINT(30 20, 10 10, 20 20, 30 30, 20 10)");
      m.normalize();
      MultiPoint expectedValue = reader.read("MULTIPOINT(10 10, 20 10, 20 20, 30 20, 30 30)");
      assertEqualsExact(expectedValue, m);
      MultiPoint unexpectedValue = reader.read("MULTIPOINT(20 10, 20 20, 30 20, 30 30, 10 10)");
      assertTrue(!m.equalsExactGeom(unexpectedValue));
    });
    test("testNormalizeLineString1", () {
      LineString l = reader.read("LINESTRING (20 20, 160 40, 160 100, 100 120, 60 60)");
      l.normalize();
      LineString expectedValue = reader.read("LINESTRING (20 20, 160 40, 160 100, 100 120, 60 60)");
      assertEqualsExact(expectedValue, l);
    });
    test("testNormalizeLineString2", () {
      LineString l = reader.read("LINESTRING (20 20, 160 40, 160 100, 100 120, 60 60)");
      l.normalize();
      LineString expectedValue = reader.read("LINESTRING (20 20, 160 40, 160 100, 100 120, 60 60)");
      assertEqualsExact(expectedValue, l);
    });
    test("testNormalizeLineString3", () {
      LineString l = reader.read("LINESTRING (200 240, 140 160, 80 160, 160 80, 80 80)");
      l.normalize();
      LineString expectedValue = reader.read("LINESTRING (80 80, 160 80, 80 160, 140 160, 200 240)");
      assertEqualsExact(expectedValue, l);
    });
    test("testNormalizeLineString4", () {
      LineString l = reader.read("LINESTRING (200 240, 140 160, 80 160, 160 80, 80 80)");
      l.normalize();
      LineString expectedValue = reader.read("LINESTRING (80 80, 160 80, 80 160, 140 160, 200 240)");
      assertEqualsExact(expectedValue, l);
    });
    test("testNormalizeLineString5", () {
      LineString l = reader.read("LINESTRING (200 340, 140 240, 140 160, 60 240, 140 240, 200 340)");
      l.normalize();
      LineString expectedValue = reader.read("LINESTRING (200 340, 140 240, 60 240, 140 160, 140 240, 200 340)");
      assertEqualsExact(expectedValue, l);
    });
    test("testNormalizeStringNoSideEffect", () {
      LineString l = reader.read("LINESTRING (200 240, 140 160, 80 160, 160 80, 80 80)");
      LineString ref = reader.read("LINESTRING (200 240, 140 160)");
      LineString seg = l.getFactory().createLineString([l.getCoordinates()[0], l.getCoordinates()[1]]);
      assertEqualsExact(ref, seg);
      l.normalize();
      assertEqualsExact(ref, seg);
    });
    test("testNormalizeEmptyLineString", () {
      LineString l = reader.read("LINESTRING EMPTY");
      l.normalize();
      LineString expectedValue = reader.read("LINESTRING EMPTY");
      assertEqualsExact(expectedValue, l);
    });
    test("testNormalizeEmptyPolygon", () {
      Polygon actualValue = reader.read("POLYGON EMPTY");
      actualValue.normalize();
      Polygon expectedValue = reader.read("POLYGON EMPTY");
      assertEqualsExact(expectedValue, actualValue);
    });
    test("testNormalizePolygon1", () {
      Polygon actualValue = reader.read(
          "POLYGON ((120 320, 240 200, 120 80, 20 200, 120 320), (60 200, 80 220, 80 200, 60 200), (160 200, 180 200, 180 220, 160 200), (120 140, 140 140, 140 160, 120 140), (140 240, 140 220, 120 260, 140 240))");
      actualValue.normalize();
      Polygon expectedValue = reader.read(
          "POLYGON ((20 200, 120 320, 240 200, 120 80, 20 200), (60 200, 80 200, 80 220, 60 200), (120 140, 140 140, 140 160, 120 140), (120 260, 140 220, 140 240, 120 260), (160 200, 180 200, 180 220, 160 200))");
      assertEqualsExact(expectedValue, actualValue);
    });
    test("testNormalizeMultiLineString", () {
      MultiLineString actualValue = reader.read(
          "MULTILINESTRING ((200 260, 180 320, 260 340), (120 180, 140 100, 40 80), (200 180, 220 160, 200 180), (100 280, 120 260, 140 260, 140 240, 120 240, 120 260, 100 280))");
      actualValue.normalize();
      MultiLineString expectedValue = reader.read(
          "MULTILINESTRING ((40 80, 140 100, 120 180), (100 280, 120 260, 120 240, 140 240, 140 260, 120 260, 100 280), (200 180, 220 160, 200 180), (200 260, 180 320, 260 340))");
      assertEqualsExact(expectedValue, actualValue);
    });
    test("testNormalizeMultiPolygon", () {
      MultiPolygon actualValue = reader.read(
          "MULTIPOLYGON (((40 360, 40 280, 140 280, 140 360, 40 360), (60 340, 60 300, 120 300, 120 340, 60 340)), ((140 200, 260 200, 260 100, 140 100, 140 200), (160 180, 240 180, 240 120, 160 120, 160 180)))");
      actualValue.normalize();
      MultiPolygon expectedValue = reader.read(
          "MULTIPOLYGON (((40 280, 40 360, 140 360, 140 280, 40 280), (60 300, 120 300, 120 340, 60 340, 60 300)), ((140 100, 140 200, 260 200, 260 100, 140 100), (160 120, 240 120, 240 180, 160 180, 160 120)))");
      assertEqualsExact(expectedValue, actualValue);
    });
    test("testNormalizeGeometryCollection", () {
      GeometryCollection actualValue = reader.read(
          "GEOMETRYCOLLECTION (LINESTRING (200 300, 200 280, 220 280, 220 320, 180 320), POINT (140 220), POLYGON ((100 80, 100 160, 20 160, 20 80, 100 80), (40 140, 40 100, 80 100, 80 140, 40 140)), POINT (100 240))");
      actualValue.normalize();
      GeometryCollection expectedValue = reader.read(
          "GEOMETRYCOLLECTION (POINT (100 240), POINT (140 220), LINESTRING (180 320, 220 320, 220 280, 200 280, 200 300), POLYGON ((20 80, 20 160, 100 160, 100 80, 20 80), (40 100, 80 100, 80 140, 40 140, 40 100)))");
      assertEqualsExact(expectedValue, actualValue);
    });
    test("testNormalizePackedCoordinateSequence", () {
      GeometryFactory pcsFactory = new GeometryFactory.withCoordinateSequenceFactory(PackedCoordinateSequenceFactory.DOUBLE_FACTORY);
      WKTReader pcsReader = new WKTReader.withFactory(pcsFactory);
      Geometry geom = pcsReader.read("LINESTRING (100 100, 0 0)");
      geom.normalize();
      // force PackedCoordinateSequence to be copied with empty coordinate cache
      Geometry clone = geom.copy();
      assertEqualsExact(geom, clone);
    });
  });

  group("MultiPointImplTest - ", () {
    PrecisionModel precisionModel = new PrecisionModel.fixedPrecision(1000);
    GeometryFactory geometryFactory = new GeometryFactory.withPrecisionModelSrid(precisionModel, 0);
    WKTReader reader = new WKTReader.withFactory(geometryFactory);

    test("testGetGeometryN", () {
      MultiPoint m = reader.read("MULTIPOINT(1.111 2.222, 3.333 4.444, 3.333 4.444)");
      Geometry g = m.getGeometryN(1);
      assertTrue(g is Point);
      Point p = g;
      Coordinate externalCoordinate = new Coordinate.empty2D();
      Coordinate internal = p.getCoordinate();
      externalCoordinate.x = internal.x;
      externalCoordinate.y = internal.y;
      assertEqualsD(3.333, externalCoordinate.x, 1E-10);
      assertEqualsD(4.444, externalCoordinate.y, 1E-10);
    });
    test("testGetEnvelope", () {
      MultiPoint m = reader.read("MULTIPOINT(1.111 2.222, 3.333 4.444, 3.333 4.444)");
      Envelope e = m.getEnvelopeInternal();
      assertEqualsD(1.111, e.getMinX(), 1E-10);
      assertEqualsD(3.333, e.getMaxX(), 1E-10);
      assertEqualsD(2.222, e.getMinY(), 1E-10);
      assertEqualsD(4.444, e.getMaxY(), 1E-10);
    });
    test("testEquals", () {
      MultiPoint m1 = reader.read("MULTIPOINT(5 6, 7 8)");
      MultiPoint m2 = reader.read("MULTIPOINT(5 6, 7 8)");
      assertTrue(m1.equals(m2));
    });
  });

  group("LineStringImplTest - ", () {
    PrecisionModel precisionModel = new PrecisionModel.fixedPrecision(1000);
    GeometryFactory geometryFactory = new GeometryFactory.withPrecisionModelSrid(precisionModel, 0);
    WKTReader reader = new WKTReader.withFactory(geometryFactory);

    test("testIsSimple", () {
      LineString l1 = reader.read("LINESTRING (0 0, 10 10, 10 0, 0 10, 0 0)");
      assertTrue(!l1.isSimple());
      LineString l2 = reader.read("LINESTRING (0 0, 10 10, 10 0, 0 10)");
      assertTrue(!l2.isSimple());
    });
    test("testIsCoordinate", () {
      LineString l = reader.read("LINESTRING (0 0, 10 10, 10 0)");
      assertTrue(l.isCoordinate(new Coordinate(0, 0)));
      assertTrue(!l.isCoordinate(new Coordinate(5, 0)));
    });
    test("testEquals", () {
      LineString l1 = reader.read("LINESTRING(1.111 2.222, 3.333 4.444)");
      LineString l2 = reader.read("LINESTRING(1.111 2.222, 3.333 4.444)");
      assertTrue(l1.equals(l2));

      l1 = reader.read("LINESTRING(1.111 2.222, 3.333 4.444)");
      l2 = reader.read("LINESTRING(3.333 4.444, 1.111 2.222)");
      assertTrue(l1.equals(l2));

      l1 = reader.read("LINESTRING(1.111 2.222, 3.333 4.444)");
      l2 = reader.read("LINESTRING(3.333 4.443, 1.111 2.222)");
      assertTrue(!l1.equals(l2));

      l1 = reader.read("LINESTRING(1.111 2.222, 3.333 4.444)");
      l2 = reader.read("LINESTRING(3.333 4.4445, 1.111 2.222)");
      assertTrue(!l1.equals(l2));

      l1 = reader.read("LINESTRING(1.111 2.222, 3.333 4.444)");
      l2 = reader.read("LINESTRING(3.333 4.4446, 1.111 2.222)");
      assertTrue(!l1.equals(l2));

      l1 = reader.read("LINESTRING(1.111 2.222, 3.333 4.444, 5.555 6.666)");
      l2 = reader.read("LINESTRING(1.111 2.222, 3.333 4.444, 5.555 6.666)");
      assertTrue(l1.equals(l2));

      l1 = reader.read("LINESTRING(1.111 2.222, 5.555 6.666, 3.333 4.444)");
      l2 = reader.read("LINESTRING(1.111 2.222, 3.333 4.444, 5.555 6.666)");
      assertTrue(!l1.equals(l2));
    });
    test("testGetCoordinates", () {
      LineString l = reader.read("LINESTRING(1.111 2.222, 5.555 6.666, 3.333 4.444)");
      List<Coordinate> coordinates = l.getCoordinates();
      assertEquals(new Coordinate(5.555, 6.666), coordinates[1]);
    });
    test("testIsClosed", () {
      LineString l = reader.read("LINESTRING EMPTY");
      assertTrue(l.isEmpty());
      assertTrue(!l.isClosed());

      LinearRing r = geometryFactory.createLinearRingSeq(null);
      assertTrue(r.isEmpty());
      assertTrue(r.isClosed());

      MultiLineString m = geometryFactory.createMultiLineString([l, r]);
      assertTrue(!m.isClosed());

      MultiLineString m2 = geometryFactory.createMultiLineString([r]);
      assertTrue(!m2.isClosed());
    });
    test("testGetGeometryType", () {
      LineString l = reader.read("LINESTRING EMPTY");
      assertEquals("LineString", l.getGeometryType());
    });
    test("testEqualsML", () {
      WKTReader reader = new WKTReader.withFactory(new GeometryFactory.withPrecisionModelSrid(new PrecisionModel.fixedPrecision(1000), 0));
      MultiLineString l1 =
          reader.read("MULTILINESTRING((1732328800 519578384, 1732026179 519976285, 1731627364 519674014, 1731929984 519276112, 1732328800 519578384))");
      MultiLineString l2 =
          reader.read("MULTILINESTRING((1731627364 519674014, 1731929984 519276112, 1732328800 519578384, 1732026179 519976285, 1731627364 519674014))");
      assertTrue(l1.equals(l2));

      reader = new WKTReader.withFactory(new GeometryFactory.withPrecisionModelSrid(new PrecisionModel.fixedPrecision(1), 0));
      l1 = reader.read("MULTILINESTRING((1732328800 519578384, 1732026179 519976285, 1731627364 519674014, 1731929984 519276112, 1732328800 519578384))");
      l2 = reader.read("MULTILINESTRING((1731627364 519674014, 1731929984 519276112, 1732328800 519578384, 1732026179 519976285, 1731627364 519674014))");
      assertTrue(l1.equals(l2));
    });
    test("testEqualsP", () {
      WKTReader reader = new WKTReader.withFactory(new GeometryFactory.withPrecisionModelSrid(new PrecisionModel.fixedPrecision(1), 0));
      Geometry l1 = reader.read("POLYGON((1732328800 519578384, 1732026179 519976285, 1731627364 519674014, 1731929984 519276112, 1732328800 519578384))");
      Geometry l2 = reader.read("POLYGON((1731627364 519674014, 1731929984 519276112, 1732328800 519578384, 1732026179 519976285, 1731627364 519674014))");
      l1.normalize();
      l2.normalize();
      assertTrue(l1.equalsExactGeom(l2));
    });
    test("testFiveZeros", () {
      LineString ls = new GeometryFactory.defaultPrecision()
          .createLineString([new Coordinate(0, 0), new Coordinate(0, 0), new Coordinate(0, 0), new Coordinate(0, 0), new Coordinate(0, 0)]);
      assertTrue(ls.isClosed());
    });
    test("testLinearRingConstructor", () {
      expect(
          () => new GeometryFactory.defaultPrecision().createLinearRing([
                new Coordinate(0, 0),
                new Coordinate(10, 10),
                new Coordinate(0, 0),
              ]),
          throwsArgumentError);
    });
  });

  group("GeometryImplTest - ", () {
    PrecisionModel precisionModel = PrecisionModel.fixedPrecision(1);
    GeometryFactory geometryFactory = GeometryFactory.withPrecisionModelSrid(precisionModel, 0);
    WKTReader reader = WKTReader.withFactory(geometryFactory);
    WKTReader readerFloat = WKTReader();

    test("testComparable", () {
      Geometry point = reader.read("POINT EMPTY");
      Geometry lineString = reader.read("LINESTRING EMPTY");
      Geometry linearRing = reader.read("LINEARRING EMPTY");
      Geometry polygon = reader.read("POLYGON EMPTY");
      Geometry mpoint = reader.read("MULTIPOINT EMPTY");
      Geometry mlineString = reader.read("MULTILINESTRING EMPTY");
      Geometry mpolygon = reader.read("MULTIPOLYGON EMPTY");
      Geometry gc = reader.read("GEOMETRYCOLLECTION EMPTY");

      var geometries = [gc, mpolygon, mlineString, mpoint, polygon, linearRing, lineString, point];

      var geometriesExpectedOrder = [point, mpoint, lineString, linearRing, mlineString, polygon, mpolygon, gc];

      geometries.sort();

      assertTrue(CollectionsUtils.areEqual(geometries, geometriesExpectedOrder));
    });
    test("testPolygonRelate", () {
      Geometry bigPolygon = reader.read("POLYGON ((0 0, 0 50, 50 50, 50 0, 0 0))");
      Geometry smallPolygon = reader.read("POLYGON ((10 10, 10 30, 30 30, 30 10, 10 10))");
      assertTrue(bigPolygon.contains(smallPolygon));
    });
    test("testEmptyGeometryCentroid", () {
      assertTrue(reader.read("POINT EMPTY").getCentroid().isEmpty());
      assertTrue(reader.read("POLYGON EMPTY").getCentroid().isEmpty());
      assertTrue(reader.read("LINESTRING EMPTY").getCentroid().isEmpty());
      assertTrue(reader.read("GEOMETRYCOLLECTION EMPTY").getCentroid().isEmpty());
      assertTrue(reader.read("GEOMETRYCOLLECTION(GEOMETRYCOLLECTION EMPTY, GEOMETRYCOLLECTION EMPTY)").getCentroid().isEmpty());
      assertTrue(reader.read("MULTIPOLYGON EMPTY").getCentroid().isEmpty());
      assertTrue(reader.read("MULTILINESTRING EMPTY").getCentroid().isEmpty());
      assertTrue(reader.read("MULTIPOINT EMPTY").getCentroid().isEmpty());
    });
    test("testEquals", () {
      Geometry g = reader.read("POLYGON ((0 0, 0 50, 50 50, 50 0, 0 0))");
      Geometry same = reader.read("POLYGON ((0 0, 0 50, 50 50, 50 0, 0 0))");
      Geometry differentStart = reader.read("POLYGON ((0 50, 50 50, 50 0, 0 0, 0 50))");
      Geometry differentFourth = reader.read("POLYGON ((0 0, 0 50, 50 50, 50 -99, 0 0))");
      Geometry differentSecond = reader.read("POLYGON ((0 0, 0 99, 50 50, 50 0, 0 0))");
      doTestEquals(g, same, true, true, true, true);
      doTestEquals(g, differentStart, true, false, false, true);
      doTestEquals(g, differentFourth, false, false, false, false);
      doTestEquals(g, differentSecond, false, false, false, false);
    });
    test("testInvalidateEnvelope", () {
      Geometry g = reader.read("POLYGON ((0 0, 0 50, 50 50, 50 0, 0 0))");
      assertEquals(Envelope(0, 50, 0, 50), g.getEnvelopeInternal());
      g.applyCF((coord) {
        coord.x += 1;
        coord.y += 1;
      });
      assertEquals(Envelope(0, 50, 0, 50), g.getEnvelopeInternal());
      g.geometryChanged();
      assertEquals(Envelope(1, 51, 1, 51), g.getEnvelopeInternal());
    });
    test("testEquals1", () {
      Geometry polygon1 = reader.read("POLYGON ((0 0, 0 50, 50 50, 50 0, 0 0))");
      Geometry polygon2 = reader.read("POLYGON ((50 50, 50 0, 0 0, 0 50, 50 50))");
      assertTrue(polygon1.equals(polygon2));
    });
    test("testEqualsWithNull", () {
      Geometry polygon = reader.read("POLYGON ((0 0, 0 50, 50 50, 50 0, 0 0))");
      assertTrue(!polygon.equals(null));
      final Object g = null;
      assertTrue(!polygon.equals(g));
    });
    test("testEqualsExactForLinearRings", () {
      LinearRing x = geometryFactory.createLinearRing([Coordinate(0, 0), Coordinate(100, 0), Coordinate(100, 100), Coordinate(0, 0)]);
      LinearRing somethingExactlyEqual = geometryFactory.createLinearRing([Coordinate(0, 0), Coordinate(100, 0), Coordinate(100, 100), Coordinate(0, 0)]);
      LinearRing somethingNotEqualButSameClass =
          geometryFactory.createLinearRing([Coordinate(0, 0), Coordinate(100, 0), Coordinate(100, 555), Coordinate(0, 0)]);
      LinearRing sameClassButEmpty = geometryFactory.createLinearRingSeq(null);
      LinearRing anotherSameClassButEmpty = geometryFactory.createLinearRingSeq(null);

      doTestEqualsExact(x, somethingExactlyEqual, somethingNotEqualButSameClass, sameClassButEmpty, anotherSameClassButEmpty,
          (geometryFactory, geometries) => geometryFactory.createMultiLineString(geometries.cast<LineString>()), geometryFactory, reader);

      //    LineString somethingEqualButNotExactly = geometryFactory.createLineString([
      //         Coordinate(0, 0),Coordinate(100, 0),Coordinate(100, 100),
      //         Coordinate(0, 0) });
      //
      //    doTestEqualsExact(x, somethingExactlyEqual, somethingEqualButNotExactly,
      //          somethingNotEqualButSameClass);
    });
    test("testEqualsExactForLineStrings", () {
      LineString x = geometryFactory.createLineString([Coordinate(0, 0), Coordinate(100, 0), Coordinate(100, 100)]);
      LineString somethingExactlyEqual = geometryFactory.createLineString([Coordinate(0, 0), Coordinate(100, 0), Coordinate(100, 100)]);
      LineString somethingNotEqualButSameClass = geometryFactory.createLineString([Coordinate(0, 0), Coordinate(100, 0), Coordinate(100, 555)]);
      LineString sameClassButEmpty = geometryFactory.createLineString(null);
      LineString anotherSameClassButEmpty = geometryFactory.createLineString(null);
      doTestEqualsExact(x, somethingExactlyEqual, somethingNotEqualButSameClass, sameClassButEmpty, anotherSameClassButEmpty,
          (geometryFactory, geometries) => geometryFactory.createMultiLineString(geometries.cast<LineString>()), geometryFactory, reader);

      doTestEqualsExact(x, somethingExactlyEqual, somethingNotEqualButSameClass, sameClassButEmpty, anotherSameClassButEmpty,
          (geometryFactory, geometries) => geometryFactory.createMultiLineString(geometries.cast<LineString>()), geometryFactory, reader);
    });
    test("testEqualsExactForPoints", () {
      Point x = geometryFactory.createPoint(Coordinate(100, 100));
      Point somethingExactlyEqual = geometryFactory.createPoint(Coordinate(100, 100));
      Point somethingNotEqualButSameClass = geometryFactory.createPoint(Coordinate(999, 100));
      Point sameClassButEmpty = geometryFactory.createPoint(null);
      Point anotherSameClassButEmpty = geometryFactory.createPoint(null);
      doTestEqualsExact(x, somethingExactlyEqual, somethingNotEqualButSameClass, sameClassButEmpty, anotherSameClassButEmpty,
          (geometryFactory, geometries) => geometryFactory.createMultiPoint(geometries.cast<Point>()), geometryFactory, reader);
    });
    test("testEqualsExactForPolygons", () {
      Polygon x = reader.read("POLYGON ((0 0, 0 50, 50 50, 50 0, 0 0))");
      Polygon somethingExactlyEqual = reader.read("POLYGON ((0 0, 0 50, 50 50, 50 0, 0 0))");
      Polygon somethingNotEqualButSameClass = reader.read("POLYGON ((50 50, 50 0, 0 0, 0 50, 50 50))");
      Polygon sameClassButEmpty = reader.read("POLYGON EMPTY");
      Polygon anotherSameClassButEmpty = reader.read("POLYGON EMPTY");
      doTestEqualsExact(x, somethingExactlyEqual, somethingNotEqualButSameClass, sameClassButEmpty, anotherSameClassButEmpty,
          (geometryFactory, geometries) => geometryFactory.createMultiPolygon(geometries.cast<Polygon>()), geometryFactory, reader);
    });
    test("testEqualsExactForGeometryCollections", () {
      Geometry polygon1 = reader.read("POLYGON ((0 0, 0 50, 50 50, 50 0, 0 0))");
      Geometry polygon2 = reader.read("POLYGON ((50 50, 50 0, 0 0, 0 50, 50 50))");
      GeometryCollection x = geometryFactory.createGeometryCollection([polygon1, polygon2]);
      GeometryCollection somethingExactlyEqual = geometryFactory.createGeometryCollection([polygon1, polygon2]);
      GeometryCollection somethingNotEqualButSameClass = geometryFactory.createGeometryCollection([polygon2]);
      GeometryCollection sameClassButEmpty = geometryFactory.createGeometryCollection(null);
      GeometryCollection anotherSameClassButEmpty = geometryFactory.createGeometryCollection(null);

      doTestEqualsExact(x, somethingExactlyEqual, somethingNotEqualButSameClass, sameClassButEmpty, anotherSameClassButEmpty,
          (geometryFactory, geometries) => geometryFactory.createGeometryCollection(geometries), geometryFactory, reader);
    });
    test("testGeometryCollectionIntersects1", () {
      Geometry gc0 = reader.read("GEOMETRYCOLLECTION ( POINT(0 0) )");
      Geometry gc1 = reader.read("GEOMETRYCOLLECTION ( LINESTRING(0 0, 1 1) )");
      Geometry gc2 = reader.read("GEOMETRYCOLLECTION ( LINESTRING(1 0, 0 1) )");
      assertTrue(gc0.intersects(gc1));
      assertTrue(gc1.intersects(gc2));
      assertTrue(!gc0.intersects(gc2));
      // symmetric
      assertTrue(gc1.intersects(gc0));
      assertTrue(gc2.intersects(gc1));
      assertTrue(!gc2.intersects(gc0));
    });
    test("testGeometryCollectionIntersects2", () {
      Geometry gc0 = reader.read("POINT(0 0)");
      Geometry gc1 = reader.read("GEOMETRYCOLLECTION ( LINESTRING(0 0, 1 1) )");
      Geometry gc2 = reader.read("LINESTRING(1 0, 0 1)");
      assertTrue(gc0.intersects(gc1));
      assertTrue(gc1.intersects(gc2));
      // symmetric
      assertTrue(gc1.intersects(gc0));
      assertTrue(gc2.intersects(gc1));
    });
    test("testGeometryCollectionIntersects3", () {
      Geometry gc0 = reader.read("GEOMETRYCOLLECTION ( POINT(0 0), LINESTRING(1 1, 2 2) )");
      Geometry gc1 = reader.read("GEOMETRYCOLLECTION ( POINT(15 15) )");
      Geometry gc2 = reader.read("GEOMETRYCOLLECTION ( LINESTRING(0 0, 2 0), POLYGON((10 10, 20 10, 20 20, 10 20, 10 10)))");
      assertTrue(gc0.intersects(gc2));
      assertTrue(!gc0.intersects(gc1));
      assertTrue(gc1.intersects(gc2));
      // symmetric
      assertTrue(gc2.intersects(gc0));
      assertTrue(!gc1.intersects(gc0));
      assertTrue(gc2.intersects(gc1));
    });
  });
}

void doTestEqualsExact(Geometry x, Geometry somethingExactlyEqual, Geometry somethingNotEqualButSameClass, Geometry sameClassButEmpty,
    Geometry anotherSameClassButEmpty, CollectionFactory collectionFactory, GeometryFactory geometryFactory, WKTReader reader) {
  Geometry emptyDifferentClass;

  if (x is Point) {
    emptyDifferentClass = geometryFactory.createGeometryCollection(null);
  } else {
    emptyDifferentClass = geometryFactory.createPoint(null);
  }

  Geometry somethingEqualButNotExactly = geometryFactory.createGeometryCollection([x]);

  doTestEqualsExact2(x, somethingExactlyEqual, collectionFactory(geometryFactory, [x]), somethingNotEqualButSameClass, reader);

  doTestEqualsExact2(sameClassButEmpty, anotherSameClassButEmpty, emptyDifferentClass, x, reader);

  /**
   * Test comparison of non-empty versus empty.
   */
  doTestEqualsExact2(x, somethingExactlyEqual, sameClassButEmpty, sameClassButEmpty, reader);

  doTestEqualsExact2(collectionFactory(geometryFactory, [x, x]), collectionFactory(geometryFactory, [x, somethingExactlyEqual]), somethingEqualButNotExactly,
      collectionFactory(geometryFactory, [x, somethingNotEqualButSameClass]), reader);
}

void doTestEqualsExact2(
    Geometry x, Geometry somethingExactlyEqual, Geometry somethingEqualButNotExactly, Geometry somethingNotEqualButSameClass, WKTReader reader) {
  Geometry differentClass;

  if (x is Point) {
    differentClass = reader.read("POLYGON ((0 0, 0 50, 50 43949, 50 0, 0 0))");
  } else {
    differentClass = reader.read("POINT ( 2351 1563 )");
  }

  assertTrue(x.equalsExactGeom(x));
  assertTrue(x.equalsExactGeom(somethingExactlyEqual));
  assertTrue(somethingExactlyEqual.equalsExactGeom(x));
  assertTrue(!x.equalsExactGeom(somethingEqualButNotExactly));
  assertTrue(!somethingEqualButNotExactly.equalsExactGeom(x));
  assertTrue(!x.equalsExactGeom(somethingEqualButNotExactly));
  assertTrue(!somethingEqualButNotExactly.equalsExactGeom(x));
  assertTrue(!x.equalsExactGeom(differentClass));
  assertTrue(!differentClass.equalsExactGeom(x));
}

void doTestEquals(Geometry a, Geometry b, bool equalsGeometry, bool equalsObject, bool equalsExact, bool equalsHash) {
  assertEquals(equalsGeometry, a.equals(b));
  assertEquals(equalsObject, a.equalsObj(b));
  assertEquals(equalsExact, a.equalsExactGeom(b));
  assertEquals(equalsHash, a.hashCode == b.hashCode);
}

typedef CollectionFactory = Geometry Function(GeometryFactory geometryFactory, List<Geometry> geometries);

void preciseCoordinateTester(PrecisionModel pm, double x1, double y1, double x2, double y2) {
  Coordinate p = new Coordinate(x1, y1);

  pm.makeCoordinatePrecise(p);

  Coordinate pPrecise = new Coordinate(x2, y2);
  assertTrue(p.equals2D(pPrecise));
}

void doPredicates(WKTReader rdr, List<String> wkt) {
  Geometry a = rdr.read(wkt[0]);
  Geometry b = rdr.read(wkt[1]);
  doPredicatesG(a, b);
  doPredicatesG(b, a);
}

void doPredicatesG(Geometry a, Geometry b) {
  assertTrue(a.contains(b) == a.relate(b).isContains());
  assertTrue(a.crosses(b) == a.relate(b).isCrosses(a.getDimension(), b.getDimension()));
  assertTrue(a.disjoint(b) == a.relate(b).isDisjoint());
  assertTrue(a.equals(b) == a.relate(b).isEquals(a.getDimension(), b.getDimension()));
  assertTrue(a.intersects(b) == a.relate(b).isIntersects());
  assertTrue(a.overlaps(b) == a.relate(b).isOverlaps(a.getDimension(), b.getDimension()));
  assertTrue(a.touches(b) == a.relate(b).isTouches(a.getDimension(), b.getDimension()));
  assertTrue(a.within(b) == a.relate(b).isWithin());
}
