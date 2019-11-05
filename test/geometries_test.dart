import "package:test/test.dart";
import 'package:dart_sfs/dart_sfs.dart';
import "dart:math" as math;
import 'testing_utilities.dart';

void main() {
  group("envelope - ", () {
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
//    test("testEmptyGeometryCentroid", () {
//      assertTrue(reader.read("POINT EMPTY").getCentroid().isEmpty());
//      assertTrue(reader.read("POLYGON EMPTY").getCentroid().isEmpty());
//      assertTrue(reader.read("LINESTRING EMPTY").getCentroid().isEmpty());
//      assertTrue(reader.read("GEOMETRYCOLLECTION EMPTY").getCentroid().isEmpty());
//      assertTrue(reader.read("GEOMETRYCOLLECTION(GEOMETRYCOLLECTION EMPTY, GEOMETRYCOLLECTION EMPTY)").getCentroid().isEmpty());
//      assertTrue(reader.read("MULTIPOLYGON EMPTY").getCentroid().isEmpty());
//      assertTrue(reader.read("MULTILINESTRING EMPTY").getCentroid().isEmpty());
//      assertTrue(reader.read("MULTIPOINT EMPTY").getCentroid().isEmpty());
//    });
//    test("testEquals", () {
//      Geometry g = reader.read("POLYGON ((0 0, 0 50, 50 50, 50 0, 0 0))");
//      Geometry same = reader.read("POLYGON ((0 0, 0 50, 50 50, 50 0, 0 0))");
//      Geometry differentStart = reader.read("POLYGON ((0 50, 50 50, 50 0, 0 0, 0 50))");
//      Geometry differentFourth = reader.read("POLYGON ((0 0, 0 50, 50 50, 50 -99, 0 0))");
//      Geometry differentSecond = reader.read("POLYGON ((0 0, 0 99, 50 50, 50 0, 0 0))");
//      doTestEquals(g, same, true, true, true, true);
//      doTestEquals(g, differentStart, true, false, false, true);
//      doTestEquals(g, differentFourth, false, false, false, false);
//      doTestEquals(g, differentSecond, false, false, false, false);
//    });
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
//    test("testEquals1", () {
//      Geometry polygon1 = reader.read("POLYGON ((0 0, 0 50, 50 50, 50 0, 0 0))");
//      Geometry polygon2 = reader.read("POLYGON ((50 50, 50 0, 0 0, 0 50, 50 50))");
//      assertTrue(polygon1.equals(polygon2));
//    });
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
//    test("testGeometryCollectionIntersects1", () {
//      Geometry gc0 = reader.read("GEOMETRYCOLLECTION ( POINT(0 0) )");
//      Geometry gc1 = reader.read("GEOMETRYCOLLECTION ( LINESTRING(0 0, 1 1) )");
//      Geometry gc2 = reader.read("GEOMETRYCOLLECTION ( LINESTRING(1 0, 0 1) )");
//      assertTrue(gc0.intersects(gc1));
//      assertTrue(gc1.intersects(gc2));
//      assertTrue(!gc0.intersects(gc2));
//      // symmetric
//      assertTrue(gc1.intersects(gc0));
//      assertTrue(gc2.intersects(gc1));
//      assertTrue(!gc2.intersects(gc0));
//    });
//    test("testGeometryCollectionIntersects2", () {
//      Geometry gc0 = reader.read("POINT(0 0)");
//      Geometry gc1 = reader.read("GEOMETRYCOLLECTION ( LINESTRING(0 0, 1 1) )");
//      Geometry gc2 = reader.read("LINESTRING(1 0, 0 1)");
//      assertTrue(gc0.intersects(gc1));
//      assertTrue(gc1.intersects(gc2));
//      // symmetric
//      assertTrue(gc1.intersects(gc0));
//      assertTrue(gc2.intersects(gc1));
//    });
//    test("testGeometryCollectionIntersects3", () {
//      Geometry gc0 = reader.read("GEOMETRYCOLLECTION ( POINT(0 0), LINESTRING(1 1, 2 2) )");
//      Geometry gc1 = reader.read("GEOMETRYCOLLECTION ( POINT(15 15) )");
//      Geometry gc2 = reader.read("GEOMETRYCOLLECTION ( LINESTRING(0 0, 2 0), POLYGON((10 10, 20 10, 20 20, 10 20, 10 10)))");
//      assertTrue(gc0.intersects(gc2));
//      assertTrue(!gc0.intersects(gc1));
//      assertTrue(gc1.intersects(gc2));
//      // symmetric
//      assertTrue(gc2.intersects(gc0));
//      assertTrue(!gc1.intersects(gc0));
//      assertTrue(gc2.intersects(gc1));
//    });
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
