import "package:test/test.dart";
import 'package:dart_jts/dart_jts.dart';
import "dart:math" as math;
import 'testing_utilities.dart';

double TOLERANCE = 1E-5;

PrecisionModel precisionModel = PrecisionModel.fixedPrecision(1);

GeometryFactory geometryFactory = GeometryFactory.withPrecisionModelSrid(precisionModel, 0);

WKTReader reader = WKTReader.withFactory(geometryFactory);

void main() {
  group("GeometryFactoryTest - ", () {
    PrecisionModel precisionModel = new PrecisionModel();
    GeometryFactory geometryFactory = new GeometryFactory.withPrecisionModelSrid(precisionModel, 0);
    WKTReader reader = new WKTReader.withFactory(geometryFactory);

    test("testCreateGeometry", () {
      checkCreateGeometryExact("POINT EMPTY");
      checkCreateGeometryExact("POINT ( 10 20 )");
      checkCreateGeometryExact("LINESTRING EMPTY");
      checkCreateGeometryExact("LINESTRING(0 0, 10 10)");
      checkCreateGeometryExact("MULTILINESTRING ((50 100, 100 200), (100 100, 150 200))");
      checkCreateGeometryExact("POLYGON ((100 200, 200 200, 200 100, 100 100, 100 200))");
      checkCreateGeometryExact("MULTIPOLYGON (((100 200, 200 200, 200 100, 100 100, 100 200)), ((300 200, 400 200, 400 100, 300 100, 300 200)))");
      checkCreateGeometryExact("GEOMETRYCOLLECTION (POLYGON ((100 200, 200 200, 200 100, 100 100, 100 200)), LINESTRING (250 100, 350 200), POINT (350 150))");
    });
    test("testCreateEmpty", () {
      checkEmpty(geometryFactory.createEmpty(0), geometryFactory.createPointEmpty().runtimeType.toString());
      checkEmpty(geometryFactory.createEmpty(1), geometryFactory.createLineStringEmpty().runtimeType.toString());
      checkEmpty(geometryFactory.createEmpty(2), geometryFactory.createPolygonEmpty().runtimeType.toString());

      checkEmpty(geometryFactory.createPointEmpty(), geometryFactory.createPointEmpty().runtimeType.toString());
      checkEmpty(geometryFactory.createLineStringEmpty(), geometryFactory.createLineStringEmpty().runtimeType.toString());
      checkEmpty(geometryFactory.createPolygonEmpty(), geometryFactory.createPolygonEmpty().runtimeType.toString());

      checkEmpty(geometryFactory.createMultiPointEmpty(), geometryFactory.createMultiPointEmpty().runtimeType.toString());
      checkEmpty(geometryFactory.createMultiLineStringEmpty(), geometryFactory.createMultiLineStringEmpty().runtimeType.toString());
      checkEmpty(geometryFactory.createMultiPolygonEmpty(), geometryFactory.createMultiPolygonEmpty().runtimeType.toString());
      checkEmpty(geometryFactory.createGeometryCollectionEmpty(), geometryFactory.createGeometryCollectionEmpty().runtimeType.toString());
    });
    test("testDeepCopy", () {
      Point g = read("POINT ( 10 10) ");
      Geometry g2 = geometryFactory.createGeometry(g);
      g.getCoordinateSequence().setOrdinate(0, 0, 99);
      assertTrue(!g.equalsExactGeom(g2));
    });
    test("testMultiPointCS", () {
      GeometryFactory gf = GeometryFactory.withCoordinateSequenceFactory(PackedCoordinateSequenceFactory());
      CoordinateSequence mpSeq = gf.getCoordinateSequenceFactory().createSizeDim(1, 4);
      mpSeq.setOrdinate(0, 0, 50);
      mpSeq.setOrdinate(0, 1, -2);
      mpSeq.setOrdinate(0, 2, 10);
      mpSeq.setOrdinate(0, 3, 20);

      MultiPoint mp = gf.createMultiPointSeq(mpSeq);
      CoordinateSequence pSeq = (mp.getGeometryN(0) as Point).getCoordinateSequence();
      assertEquals(4, pSeq.getDimension());
      for (int i = 0; i < 4; i++) assertEquals(mpSeq.getOrdinate(0, i), pSeq.getOrdinate(0, i));
    });
    test("testCopyGeometryWithNonDefaultDimension", () {
      GeometryFactory gf = GeometryFactory.withCoordinateSequenceFactory(CoordinateArraySequenceFactory());
      CoordinateSequence mpSeq = gf.getCoordinateSequenceFactory().createSizeDim(1, 2);
      mpSeq.setOrdinate(0, 0, 50);
      mpSeq.setOrdinate(0, 1, -2);

      Point g = gf.createPointSeq(mpSeq);
      CoordinateSequence pSeq = (g.getGeometryN(0) as Point).getCoordinateSequence();
      assertEquals(2, pSeq.getDimension());

      Point g2 = geometryFactory.createGeometry(g);
      assertEquals(2, g2.getCoordinateSequence().getDimension());
    });
  });
  group("GeometryCopyTest - ", () {
    test("testCopy", () {
      checkCopy(read(WKT_POINT));
      checkCopy(read(WKT_LINESTRING));
      checkCopy(read(WKT_LINEARRING));
      checkCopy(read(WKT_POLY));
      checkCopy(read(WKT_MULTIPOINT));
      checkCopy(read(WKT_MULTILINESTRING));
      checkCopy(read(WKT_MULTIPOLYGON));
      checkCopy(read(WKT_GC));
    });
  });
  group("GeometryCollectionIteratorTest - ", () {
    test("testGeometryCollection", () {
      GeometryCollection g = read("GEOMETRYCOLLECTION (GEOMETRYCOLLECTION (POINT (10 10)))");
      GeometryCollectionIterator i = GeometryCollectionIterator(g);
      assertTrue(i.hasNext());
      assertTrue(i.next() is GeometryCollection);
      assertTrue(i.hasNext());
      assertTrue(i.next() is GeometryCollection);
      assertTrue(i.hasNext());
      assertTrue(i.next() is Point);
      assertTrue(!i.hasNext());
    });
    test("testAtomic", () {
      Polygon g = read("POLYGON ((1 9, 9 9, 9 1, 1 1, 1 9))");
      GeometryCollectionIterator i = GeometryCollectionIterator(g);
      assertTrue(i.hasNext());
      assertTrue(i.next() is Polygon);
      assertTrue(!i.hasNext());
    });
  });
  group("GeometryCollectionImplTest - ", () {
    PrecisionModel precisionModel1 = PrecisionModel.fixedPrecision(1000);
    GeometryFactory geometryFactory1 = GeometryFactory.withPrecisionModelSrid(precisionModel1, 0);
    WKTReader reader1 = WKTReader.withFactory(geometryFactory1);
    test("testGetDimension", () {
      GeometryCollection g = reader1.read("GEOMETRYCOLLECTION (POINT (10 10), POINT (30 30), LINESTRING (15 15, 20 20))");
      assertEquals(1, g.getDimension());
    });
    test("testGetCoordinates", () {
      GeometryCollection g = reader1.read("GEOMETRYCOLLECTION (POINT (10 10), POINT (30 30), LINESTRING (15 15, 20 20))");
      var coordinates = g.getCoordinates();
      assertEquals(4, g.getNumPoints());
      assertEquals(4, coordinates.length);
      assertEquals(new Coordinate(10, 10), coordinates[0]);
      assertEquals(new Coordinate(20, 20), coordinates[3]);
    });
    test("testGeometryCollectionIterator", () {
      GeometryCollection g = reader1.read("GEOMETRYCOLLECTION (GEOMETRYCOLLECTION (POINT (10 10)))");
      GeometryCollectionIterator i = GeometryCollectionIterator(g);
      assertTrue(i.hasNext());
      assertTrue(i.next() is GeometryCollection);
      assertTrue(i.next() is GeometryCollection);
      assertTrue(i.next() is Point);
    });
    test("testGetLength", () {
      GeometryCollection g = WKTReader().read("MULTIPOLYGON(" +
          "((0 0, 10 0, 10 10, 0 10, 0 0), (3 3, 3 7, 7 7, 7 3, 3 3))," +
          "((100 100, 110 100, 110 110, 100 110, 100 100), (103 103, 103 107, 107 107, 107 103, 103 103)))");
      assertEqualsD(112, g.getLength(), 1E-15);
    });
  });
  group("CoordinateTest - ", () {
    test("testConstructor3D", () {
      Coordinate c = Coordinate.fromXYZ(350.2, 4566.8, 5266.3);
      assertEquals(c.x, 350.2);
      assertEquals(c.y, 4566.8);
      assertEquals(c.getZ(), 5266.3);
    });
    test("testConstructor2D", () {
      Coordinate c = Coordinate(350.2, 4566.8);
      assertEquals(c.x, 350.2);
      assertEquals(c.y, 4566.8);
      assertEquals(c.getZ(), Coordinate.NULL_ORDINATE);
    });
    test("testDefaultConstructor", () {
      Coordinate c = Coordinate.empty2D();
      assertEquals(c.x, 0.0);
      assertEquals(c.y, 0.0);
      assertEquals(c.getZ(), Coordinate.NULL_ORDINATE);
    });
    test("testCopyConstructor3D", () {
      Coordinate orig = Coordinate.fromXYZ(350.2, 4566.8, 5266.3);
      Coordinate c = Coordinate.fromCoordinate(orig);
      assertEquals(c.x, 350.2);
      assertEquals(c.y, 4566.8);
      assertEquals(c.getZ(), 5266.3);
    });
    test("testSetCoordinate", () {
      Coordinate orig = Coordinate.fromXYZ(350.2, 4566.8, 5266.3);
      Coordinate c = Coordinate.empty2D();
      c.setCoordinate(orig);
      assertEquals(c.x, 350.2);
      assertEquals(c.y, 4566.8);
      assertEquals(c.getZ(), 5266.3);
    });
    test("testGetOrdinate", () {
      Coordinate c = Coordinate.fromXYZ(350.2, 4566.8, 5266.3);
      assertEquals(c.getOrdinate(Coordinate.X), 350.2);
      assertEquals(c.getOrdinate(Coordinate.Y), 4566.8);
      assertEquals(c.getOrdinate(Coordinate.Z), 5266.3);
    });
    test("testSetOrdinate", () {
      Coordinate c = Coordinate.empty2D();
      c.setOrdinate(Coordinate.X, 111);
      c.setOrdinate(Coordinate.Y, 222);
      c.setOrdinate(Coordinate.Z, 333);
      assertEquals(c.getOrdinate(Coordinate.X), 111.0);
      assertEquals(c.getOrdinate(Coordinate.Y), 222.0);
      assertEquals(c.getOrdinate(Coordinate.Z), 333.0);
    });
    test("testEquals", () {
      Coordinate c1 = Coordinate.fromXYZ(1, 2, 3);
      String s = "Not a coordinate";
      assertTrue(!c1.equals(s));

      Coordinate c2 = Coordinate.fromXYZ(1, 2, 3);
      assertTrue(c1.equals2D(c2));

      Coordinate c3 = Coordinate.fromXYZ(1, 22, 3);
      assertTrue(!c1.equals2D(c3));
    });
    test("testEquals2D", () {
      Coordinate c1 = Coordinate.fromXYZ(1, 2, 3);
      Coordinate c2 = Coordinate.fromXYZ(1, 2, 3);
      assertTrue(c1.equals2D(c2));

      Coordinate c3 = Coordinate.fromXYZ(1, 22, 3);
      assertTrue(!c1.equals2D(c3));
    });
    test("testEquals3D", () {
      Coordinate c1 = Coordinate.fromXYZ(1, 2, 3);
      Coordinate c2 = Coordinate.fromXYZ(1, 2, 3);
      assertTrue(c1.equals3D(c2));

      Coordinate c3 = Coordinate.fromXYZ(1, 22, 3);
      assertTrue(!c1.equals3D(c3));
    });
    test("testEquals2DWithinTolerance", () {
      Coordinate c = Coordinate.fromXYZ(100.0, 200.0, 50.0);
      Coordinate aBitOff = Coordinate.fromXYZ(100.1, 200.1, 50.0);
      assertTrue(c.equals2DWithTolerance(aBitOff, 0.2));
    });
    test("testEqualsInZ", () {
      Coordinate c = Coordinate.fromXYZ(100.0, 200.0, 50.0);
      Coordinate withSameZ = Coordinate.fromXYZ(100.1, 200.1, 50.1);
      assertTrue(c.equalInZ(withSameZ, 0.2));
    });
    test("testCompareTo", () {
      Coordinate lowest = Coordinate.fromXYZ(10.0, 100.0, 50.0);
      Coordinate highest = Coordinate.fromXYZ(20.0, 100.0, 50.0);
      Coordinate equalToHighest = Coordinate.fromXYZ(20.0, 100.0, 50.0);
      Coordinate higherStill = Coordinate.fromXYZ(20.0, 200.0, 50.0);

      assertEquals(-1, lowest.compareTo(highest));
      assertEquals(1, highest.compareTo(lowest));
      assertEquals(-1, highest.compareTo(higherStill));
      assertEquals(0, highest.compareTo(equalToHighest));
    });
    test("testToString", () {
      String expectedResult = "(100.0, 200.0, 50.0)";
      String actualResult = Coordinate.fromXYZ(100.0, 200.0, 50.0).toString();
      assertEquals(expectedResult, actualResult);
    });
    test("testClone", () {
      Coordinate c = Coordinate.fromXYZ(100.0, 200.0, 50.0);
      Coordinate clone = c.clone();
      assertTrue(c.equals3D(clone));
    });
    test("testDistance", () {
      Coordinate coord1 = Coordinate.fromXYZ(0.0, 0.0, 0.0);
      Coordinate coord2 = Coordinate.fromXYZ(100.0, 200.0, 50.0);
      double distance = coord1.distance(coord2);
      assertEqualsD(distance, 223.60679774997897, 0.00001);
    });
    test("testDistance3D", () {
      Coordinate coord1 = Coordinate.fromXYZ(0.0, 0.0, 0.0);
      Coordinate coord2 = Coordinate.fromXYZ(100.0, 200.0, 50.0);
      double distance = coord1.distance3D(coord2);
      assertEqualsD(distance, 229.128784747792, 0.000001);
    });
    test("testCoordinateXY", () {
      Coordinate xy = CoordinateXY();
      checkZUnsupported(xy);
      checkMUnsupported(xy);

      xy = CoordinateXY.fromXY(1.0, 1.0); // 2D
      Coordinate coord = Coordinate.fromCoordinate(xy); // copy
      assertEquals(xy, coord);
      assertTrue(!xy.equalInZ(coord, 0.000001));

      coord = Coordinate.fromXYZ(1.0, 1.0, 1.0); // 2.5d
      xy = CoordinateXY.fromCoordinate(coord); // copy
      assertEquals(xy, coord);
      assertTrue(!xy.equalInZ(coord, 0.000001));
    });
    test("testCoordinateXYM", () {
      Coordinate xym = CoordinateXYM.empty();
      checkZUnsupported(xym);

      xym.setM(1.0);
      assertEquals(1.0, xym.getM());

      Coordinate coord = Coordinate.fromCoordinate(xym); // copy
      assertEquals(xym, coord);
      assertTrue(!xym.equalInZ(coord, 0.000001));

      coord = Coordinate.fromXYZ(1.0, 1.0, 1.0); // 2.5d
      xym = CoordinateXYM.fromCoordinate(coord); // copy
      assertEquals(xym, coord);
      assertTrue(!xym.equalInZ(coord, 0.000001));
    });
    test("testCoordinateXYZM", () {
      Coordinate xyzm = CoordinateXYZM.empty();
      xyzm.setZ(1.0);
      assertEquals(1.0, xyzm.getZ());
      xyzm.setM(1.0);
      assertEquals(1.0, xyzm.getM());

      Coordinate coord = Coordinate.fromCoordinate(xyzm); // copy
      assertEquals(xyzm, coord);
      assertTrue(xyzm.equalInZ(coord, 0.000001));
      assertTrue(coord.getM().isNaN);

      coord = Coordinate.fromXYZ(1.0, 1.0, 1.0); // 2.5d
      xyzm = CoordinateXYZM.fromCoordinate(coord); // copy
      assertEquals(xyzm, coord);
      assertTrue(xyzm.equalInZ(coord, 0.000001));
    });
    test("testCoordinateHash", () {
      doTestCoordinateHash(true, Coordinate(1, 2), Coordinate(1, 2));
      doTestCoordinateHash(false, Coordinate(1, 2), Coordinate(3, 4));
      doTestCoordinateHash(false, Coordinate(1, 2), Coordinate(1, 4));
      doTestCoordinateHash(false, Coordinate(1, 2), Coordinate(3, 2));
      doTestCoordinateHash(false, Coordinate(1, 2), Coordinate(2, 1));
    });
  });
  group("AreaLengthTest - ", () {
    test("testLength", () {
      checkLength("MULTIPOINT (220 140, 180 280)", 0.0);
      checkLength("LINESTRING (220 140, 180 280)", 145.6021977);
      checkLength("LINESTRING (0 0, 100 100)", 141.4213562373095);
      checkLength("POLYGON ((20 20, 40 20, 40 40, 20 40, 20 20))", 80.0);
      checkLength("POLYGON ((20 20, 40 20, 40 40, 20 40, 20 20), (25 35, 35 35, 35 25, 25 25, 25 35))", 120.0);
    });
    test("testArea", () {
      checkArea("MULTIPOINT (220 140, 180 280)", 0.0);
      checkArea("LINESTRING (220 140, 180 280)", 0.0);
      checkArea("POLYGON ((20 20, 40 20, 40 40, 20 40, 20 20))", 400.0);
      checkArea("POLYGON ((20 20, 40 20, 40 40, 20 40, 20 20), (25 35, 35 35, 35 25, 25 25, 25 35))", 300.0);
    });
  });
  group("BidirectionalComparatorTest - ", () {
    test("testLineString1", () {
      expect(
          0 ==
              compareBiDir(
                  "LINESTRING ( 1388155.775 794886.703, 1388170.712 794887.346, 1388185.425 794892.987, 1388195.167 794898.409, 1388210.091 794899.06, 1388235.117 794900.145, 1388250.276 794895.796, 1388270.174 794896.648, 1388280.138 794897.079, 1388295.063 794897.731, 1388310.348 794893.382, 1388330.479 794889.255, 1388345.617 794884.895, 1388360.778 794880.538, 1388366.184 794870.766, 1388366.62 794860.776, 1388362.086 794850.563, 1388357.761 794835.234, 1388343.474 794819.588, 1388339.151 794804.386, 1388320.114 794783.54, 1388310.597 794773.107, 1388301.155 794757.682, 1388286.452 794751.914, 1388282.129 794736.7, 1388273.037 794716.275, 1388278.444 794706.504, 1388293.603 794702.155, 1388303.994 794692.585, 1388319.278 794688.247, 1388339.4 794684.108, 1388369.486 794680.401, 1388394.513 794681.487, 1388409.429 794682.126, 1388433.884 794693.192, 1388454.204 794698.202 )",
                  "LINESTRING ( 1388454.204 794698.202, 1388433.884 794693.192, 1388409.429 794682.126, 1388394.513 794681.487, 1388369.486 794680.401, 1388339.4 794684.108, 1388319.278 794688.247, 1388303.994 794692.585, 1388293.603 794702.155, 1388278.444 794706.504, 1388273.037 794716.275, 1388282.129 794736.7, 1388286.452 794751.914, 1388301.155 794757.682, 1388310.597 794773.107, 1388320.114 794783.54, 1388339.151 794804.386, 1388343.474 794819.588, 1388357.761 794835.234, 1388362.086 794850.563, 1388366.62 794860.776, 1388366.184 794870.766, 1388360.778 794880.538, 1388345.617 794884.895, 1388330.479 794889.255, 1388310.348 794893.382, 1388295.063 794897.731, 1388280.138 794897.079, 1388270.174 794896.648, 1388250.276 794895.796, 1388235.117 794900.145, 1388210.091 794899.06, 1388195.167 794898.409, 1388185.425 794892.987, 1388170.712 794887.346, 1388155.775 794886.703 )"),
          true);
    });
    test("testLineString2", () {
      expect(0 == compareBiDir("LINESTRING (1389103.293 794193.755, 1389064.931 794188.991)", "LINESTRING (1389064.931 794188.991, 1389103.293 794193.755)"),
          true);
    });
  });
  group("CoordinateArraysTest - ", () {
    var COORDS_1 = [Coordinate(1, 1), Coordinate(2, 2), Coordinate(3, 3)];
    var COORDS_EMPTY = <Coordinate>[];

    test("testPtNotInList1", () {
      expect(
          CoordinateArrays.ptNotInList([Coordinate(1, 1), Coordinate(2, 2), Coordinate(3, 3)], [Coordinate(1, 1), Coordinate(1, 2), Coordinate(1, 3)])
              .equals2D(Coordinate(2, 2)),
          true);
    });
    test("testPtNotInList2", () {
      expect(
          CoordinateArrays.ptNotInList([Coordinate(1, 1), Coordinate(2, 2), Coordinate(3, 3)], [Coordinate(1, 1), Coordinate(2, 2), Coordinate(3, 3)]) == null,
          true);
    });
    test("testEnvelope1", () {
      expect(CoordinateArrays.envelope(COORDS_1), Envelope(1, 3, 1, 3));
    });
    test("testEnvelopeEmpty", () {
      expect(CoordinateArrays.envelope(COORDS_EMPTY), Envelope.empty());
    });
    test("testIntersection_envelope1", () {
      expect(CoordinateArrays.equals(CoordinateArrays.intersection(COORDS_1, Envelope(1, 2, 1, 2)), [Coordinate(1, 1), Coordinate(2, 2)]), true);
    });
    test("testIntersection_envelopeDisjoint", () {
      expect(CoordinateArrays.equals(CoordinateArrays.intersection(COORDS_1, Envelope(10, 20, 10, 20)), COORDS_EMPTY), true);
    });
    test("testIntersection_empty_envelope", () {
      expect(CoordinateArrays.equals(CoordinateArrays.intersection(COORDS_EMPTY, Envelope(1, 2, 1, 2)), COORDS_EMPTY), true);
    });
    test("testIntersection_coords_emptyEnvelope", () {
      expect(CoordinateArrays.equals(CoordinateArrays.intersection(COORDS_1, Envelope.empty()), COORDS_EMPTY), true);
    });
  });

  group("CoordinateListTest - ", () {
    test("testForward and reverse", () {
      checkValue(coordList([0, 0, 1, 1, 2, 2]), [0, 0, 1, 1, 2, 2]);
      checkValue(coordList([0, 0, 1, 1, 2, 2]).reversed.toList(), [2, 2, 1, 1, 0, 0]);
    });
  });

  group("CoordinateSequencesTest - ", () {
    test("testCopyToLargerDim", () {
      PackedCoordinateSequenceFactory csFactory = PackedCoordinateSequenceFactory();
      CoordinateSequence cs2D = createTestSequence(csFactory, 10, 2);
      CoordinateSequence cs3D = csFactory.createSizeDim(10, 3);
      CoordinateSequences.copy(cs2D, 0, cs3D, 0, cs3D.size());
      expect(CoordinateSequences.isEqual(cs2D, cs3D), true);
    });
    test("testCopyToSmallerDim", () {
      PackedCoordinateSequenceFactory csFactory = PackedCoordinateSequenceFactory();
      CoordinateSequence cs3D = createTestSequence(csFactory, 10, 3);
      CoordinateSequence cs2D = csFactory.createSizeDim(10, 2);
      CoordinateSequences.copy(cs3D, 0, cs2D, 0, cs2D.size());
      expect(CoordinateSequences.isEqual(cs2D, cs3D), true);
    });
    test("testScrollRing", () {
      doTestScrollRing(CoordinateArraySequenceFactory(), 2);
      doTestScrollRing(CoordinateArraySequenceFactory(), 3);
      doTestScrollRing(PackedCoordinateSequenceFactory.DOUBLE_FACTORY, 2);
      doTestScrollRing(PackedCoordinateSequenceFactory.DOUBLE_FACTORY, 4);
    });
    test("testScroll", () {
      doTestScroll(CoordinateArraySequenceFactory(), 2);
      doTestScroll(CoordinateArraySequenceFactory(), 3);
      doTestScroll(PackedCoordinateSequenceFactory.DOUBLE_FACTORY, 2);
      doTestScroll(PackedCoordinateSequenceFactory.DOUBLE_FACTORY, 4);
    });
    test("testIndexOf", () {
      doTestIndexOf(CoordinateArraySequenceFactory(), 2);
      doTestIndexOf(PackedCoordinateSequenceFactory.DOUBLE_FACTORY, 5);
    });
    test("testMinCoordinateIndex", () {
      doTestMinCoordinateIndex(CoordinateArraySequenceFactory(), 2);
      doTestMinCoordinateIndex(PackedCoordinateSequenceFactory.DOUBLE_FACTORY, 5);
    });
    test("testIsRing", () {
      doTestIsRing(CoordinateArraySequenceFactory(), 2);
      doTestIsRing(PackedCoordinateSequenceFactory.DOUBLE_FACTORY, 5);
    });
    test("testCopy", () {
      doTestCopy(CoordinateArraySequenceFactory(), 2);
      doTestCopy(PackedCoordinateSequenceFactory.DOUBLE_FACTORY, 5);
    });
    test("testReverse", () {
      doTestReverse(CoordinateArraySequenceFactory(), 2);
      doTestReverse(PackedCoordinateSequenceFactory.DOUBLE_FACTORY, 5);
    });
  });
}

int compareBiDir(String wkt0, String wkt1) {
  LineString g0 = reader.read(wkt0) as LineString;
  LineString g1 = reader.read(wkt1) as LineString;
  var pts0 = g0.getCoordinates();
  var pts1 = g1.getCoordinates();
  return CoordinateArrays.bidirectionalComparator(pts0, pts1);
}

void checkLength(String wkt, double expectedValue) {
  Geometry g = reader.read(wkt);
  double len = g.getLength();
//		//System.out.println(len);
  expect(NumberUtils.equalsWithTolerance(expectedValue, len, TOLERANCE), true);
}

void checkCopy(final Geometry g) {
  int SRID = 123;
  g.setSRID(SRID);

  Object DATA = 999;
  g.setUserData(DATA);

  Geometry copy = g.copy();

  assertEquals(g.getSRID(), copy.getSRID());
  assertEquals(g.getUserData(), copy.getUserData());

  //TODO: use a test which checks all ordinates of CoordinateSequences
  assertTrue(g.equalsExactGeom(copy));
}

void checkArea(String wkt, double expectedValue) {
  Geometry g = reader.read(wkt);
  expect(NumberUtils.equalsWithTolerance(expectedValue, g.getArea(), TOLERANCE), true);
}

void checkValue(List<Coordinate> coordArray, List<double> ords) {
  expect(coordArray.length * 2, ords.length);

  for (int i = 0; i < coordArray.length; i += 2) {
    Coordinate pt = coordArray[i];
    expect(pt.x, ords[2 * i]);
    expect(pt.y, ords[2 * i + 1]);
  }
}

List<Coordinate> coordList(List<double> ords) {
  List<Coordinate> cl = [];
  for (int i = 0; i < ords.length; i += 2) {
    cl.add(Coordinate(ords[i], ords[i + 1]));
//    CollectionsUtils.addIfNotEqualToLast(cl, Coordinate(ords[i], ords[i + 1]));
  }
  return cl;
}

final List<List<double>> ordinateValues = [
  [75.76, 77.43],
  [41.35, 90.75],
  [73.74, 41.67],
  [20.87, 86.49],
  [17.49, 93.59],
  [67.75, 80.63],
  [63.01, 52.57],
  [32.9, 44.44],
  [79.36, 29.8],
  [38.17, 88.0],
  [19.31, 49.71],
  [57.03, 19.28],
  [63.76, 77.35],
  [45.26, 85.15],
  [51.71, 50.38],
  [92.16, 19.85],
  [64.18, 27.7],
  [64.74, 65.1],
  [80.07, 13.55],
  [55.54, 94.07]
];

CoordinateSequence createSequenceFromOrdinates(CoordinateSequenceFactory csFactory, int dim) {
  CoordinateSequence sequence = csFactory.createSizeDim(ordinateValues.length, dim);
  for (int i = 0; i < ordinateValues.length; i++) {
    sequence.setOrdinate(i, 0, ordinateValues[i][0]);
    sequence.setOrdinate(i, 1, ordinateValues[i][1]);
  }
  return fillNonPlanarDimensions(sequence);
}

CoordinateSequence createTestSequence(CoordinateSequenceFactory csFactory, int size, int dim) {
  CoordinateSequence cs = csFactory.createSizeDim(size, dim);
  // initialize with a data signature where coords look like [1, 10, 100, ...]
  for (int i = 0; i < size; i++) {
    for (int d = 0; d < dim; d++) {
      cs.setOrdinate(i, d, (i * math.pow(10, d)).toDouble());
    }
  }
  return cs;
}

/// @deprecated only use to update in conjunction with {@link this.ttestCreateRandomOrdinates}
CoordinateSequence createRandomTestSequence(CoordinateSequenceFactory csFactory, int size, int dim, math.Random rnd, Envelope range, PrecisionModel pm) {
  CoordinateSequence cs = csFactory.createSizeDim(size, dim);
  for (int i = 0; i < size; i++) {
    cs.setOrdinate(i, 0, pm.makePrecise(range.getWidth() * rnd.nextDouble() + range.getMinX()));
    cs.setOrdinate(i, 1, pm.makePrecise(range.getHeight() * rnd.nextDouble() + range.getMinY()));
  }

  return fillNonPlanarDimensions(cs);
}

void doTestReverse(CoordinateSequenceFactory factory, int dimension) {
  // arrange
  CoordinateSequence sequence = createSequenceFromOrdinates(factory, dimension);
  CoordinateSequence reversed = sequence.copy();

  // act
  CoordinateSequences.reverse(reversed);

  // assert
  for (int i = 0; i < sequence.size(); i++) checkCoordinateAt(sequence, i, reversed, sequence.size() - i - 1, dimension);
}

void doTestCopy(CoordinateSequenceFactory factory, int dimension) {
  // arrange
  CoordinateSequence sequence = createSequenceFromOrdinates(factory, dimension);
  if (sequence.size() <= 7) {
    print("sequence has a size of ${sequence.size()}. Execution of this test needs a sequence " + "with more than 6 coordinates.");
    return;
  }

  CoordinateSequence fullCopy = factory.createSizeDim(sequence.size(), dimension);
  CoordinateSequence partialCopy = factory.createSizeDim(sequence.size() - 5, dimension);

  // act
  CoordinateSequences.copy(sequence, 0, fullCopy, 0, sequence.size());
  CoordinateSequences.copy(sequence, 2, partialCopy, 0, partialCopy.size());

  // assert
  for (int i = 0; i < fullCopy.size(); i++) checkCoordinateAt(sequence, i, fullCopy, i, dimension);
  for (int i = 0; i < partialCopy.size(); i++) checkCoordinateAt(sequence, 2 + i, partialCopy, i, dimension);

  // ToDo test if dimensions don't match
}

void doTestIsRing(CoordinateSequenceFactory factory, int dimension) {
  // arrange
  CoordinateSequence ring = createCircle(factory, dimension, Coordinate.empty2D(), 5);
  CoordinateSequence noRing = createCircularString(factory, dimension, Coordinate.empty2D(), 5, 0.1, 22);
  CoordinateSequence empty = createAlmostRing(factory, dimension, 0);
  CoordinateSequence incomplete1 = createAlmostRing(factory, dimension, 1);
  CoordinateSequence incomplete2 = createAlmostRing(factory, dimension, 2);
  CoordinateSequence incomplete3 = createAlmostRing(factory, dimension, 3);
  CoordinateSequence incomplete4a = createAlmostRing(factory, dimension, 4);
  CoordinateSequence incomplete4b = CoordinateSequences.ensureValidRing(factory, incomplete4a);

  // act
  bool isRingRing = CoordinateSequences.isRing(ring);
  bool isRingNoRing = CoordinateSequences.isRing(noRing);
  bool isRingEmpty = CoordinateSequences.isRing(empty);
  bool isRingIncomplete1 = CoordinateSequences.isRing(incomplete1);
  bool isRingIncomplete2 = CoordinateSequences.isRing(incomplete2);
  bool isRingIncomplete3 = CoordinateSequences.isRing(incomplete3);
  bool isRingIncomplete4a = CoordinateSequences.isRing(incomplete4a);
  bool isRingIncomplete4b = CoordinateSequences.isRing(incomplete4b);

  // assert
  expect(isRingRing, true);
  expect(!isRingNoRing, true);
  expect(isRingEmpty, true);
  expect(!isRingIncomplete1, true);
  expect(!isRingIncomplete2, true);
  expect(!isRingIncomplete3, true);
  expect(!isRingIncomplete4a, true);
  expect(isRingIncomplete4b, true);
}

void doTestIndexOf(CoordinateSequenceFactory factory, int dimension) {
  // arrange
  CoordinateSequence sequence = createSequenceFromOrdinates(factory, dimension);

  // act & assert
  List<Coordinate> coordinates = sequence.toCoordinateArray();
  for (int i = 0; i < sequence.size(); i++) expect(i, CoordinateSequences.indexOf(coordinates[i], sequence));
}

void doTestMinCoordinateIndex(CoordinateSequenceFactory factory, int dimension) {
  CoordinateSequence sequence = createSequenceFromOrdinates(factory, dimension);
  if (sequence.size() <= 6) {
    print("sequence has a size of ${sequence.size()}. Execution of this test needs a sequence " + "with more than 5 coordinates.");
    return;
  }

  int minIndex = (sequence.size() / 2).round();
  sequence.setOrdinate(minIndex, 0, 5);
  sequence.setOrdinate(minIndex, 1, 5);

  expect(minIndex, CoordinateSequences.minCoordinateIndex(sequence));
  expect(minIndex, CoordinateSequences.minCoordinateIndexWithRange(sequence, 2, sequence.size() - 2));
}

void doTestScroll(CoordinateSequenceFactory factory, int dimension) {
  // arrange
  CoordinateSequence sequence = createCircularString(factory, dimension, Coordinate(20, 20), 7.0, 0.1, 22);
  CoordinateSequence scrolled = sequence.copy();

  // act
  CoordinateSequences.scrollWithIndex(scrolled, 12);

  // assert
  int io = 12;
  for (int iis = 0; iis < scrolled.size() - 1; iis++) {
    checkCoordinateAt(sequence, io, scrolled, iis, dimension);
    io++;
    io %= scrolled.size();
  }
}

void doTestScrollRing(CoordinateSequenceFactory factory, int dimension) {
  // arrange
  //System.out.println("Testing '" + factory.getClass().getSimpleName() + "' with dim=" +dimension );
  CoordinateSequence sequence = createCircle(factory, dimension, Coordinate(10, 10), 9.0);
  CoordinateSequence scrolled = sequence.copy();

  // act
  CoordinateSequences.scrollWithIndex(scrolled, 12);

  // assert
  int io = 12;
  for (int iis = 0; iis < scrolled.size() - 1; iis++) {
    checkCoordinateAt(sequence, io, scrolled, iis, dimension);
    io++;
    io %= scrolled.size() - 1;
  }
  checkCoordinateAt(scrolled, 0, scrolled, scrolled.size() - 1, dimension);
}

void checkCoordinateAt(CoordinateSequence seq1, int pos1, CoordinateSequence seq2, int pos2, int dim) {
  expect(seq1.getOrdinate(pos1, 0), seq2.getOrdinate(pos2, 0));
  expect(seq1.getOrdinate(pos1, 1), seq2.getOrdinate(pos2, 1));

  // check additional ordinates
  for (int j = 2; j < dim; j++) {
    expect(seq1.getOrdinate(pos1, j), seq2.getOrdinate(pos2, j));
  }
}

CoordinateSequence createAlmostRing(CoordinateSequenceFactory factory, int dimension, int num) {
  if (num > 4) num = 4;

  CoordinateSequence sequence = factory.createSizeDim(num, dimension);
  if (num == 0) return fillNonPlanarDimensions(sequence);

  sequence.setOrdinate(0, 0, 10);
  sequence.setOrdinate(0, 0, 10);
  if (num == 1) return fillNonPlanarDimensions(sequence);

  sequence.setOrdinate(0, 0, 20);
  sequence.setOrdinate(0, 0, 10);
  if (num == 2) return fillNonPlanarDimensions(sequence);

  sequence.setOrdinate(0, 0, 20);
  sequence.setOrdinate(0, 0, 20);
  if (num == 3) return fillNonPlanarDimensions(sequence);

  sequence.setOrdinate(0, 0, 10.0000000000001);
  sequence.setOrdinate(0, 0, 9.9999999999999);
  return fillNonPlanarDimensions(sequence);
}

CoordinateSequence fillNonPlanarDimensions(CoordinateSequence seq) {
  if (seq.getDimension() < 3) return seq;

  for (int i = 0; i < seq.size(); i++) {
    for (int j = 2; j < seq.getDimension(); j++) {
      seq.setOrdinate(i, j, (i * math.pow(10, j - 1)).toDouble());
    }
  }

  return seq;
}

CoordinateSequence createCircle(CoordinateSequenceFactory factory, int dimension, Coordinate center, double radius) {
  // Get a complete circular string
  CoordinateSequence res = createCircularString(factory, dimension, center, radius, 0, 49);

  // ensure it is closed
  for (int i = 0; i < dimension; i++) res.setOrdinate(48, i, res.getOrdinate(0, i));

  return res;
}

CoordinateSequence createCircularString(CoordinateSequenceFactory factory, int dimension, Coordinate center, double radius, double startAngle, int numPoints) {
  final int numSegmentsCircle = 48;
  final double angleCircle = 2 * math.pi;
  final double angleStep = angleCircle / numSegmentsCircle;

  CoordinateSequence sequence = factory.createSizeDim(numPoints, dimension);
  PrecisionModel pm = PrecisionModel.fixedPrecision(100);
  double angle = startAngle;
  for (int i = 0; i < numPoints; i++) {
    double dx = math.cos(angle) * radius;
    sequence.setOrdinate(i, 0, pm.makePrecise(center.x + dx));
    double dy = math.sin(angle) * radius;
    sequence.setOrdinate(i, 1, pm.makePrecise(center.y + dy));

    // set other ordinate values to predictable values
    for (int j = 2; j < dimension; j++) {
      sequence.setOrdinate(i, j, (math.pow(10, j - 1) * i).toDouble());
    }

    angle += angleStep;
    angle %= angleCircle;
  }

  return sequence;
}

void doTestCoordinateHash(bool equal, Coordinate a, Coordinate b) {
  assertEquals(equal, a.equals(b));
  assertEquals(equal, a.hashCode == b.hashCode);
}

/**
 * Confirm the z field is not supported by getZ and setZ.
 */
void checkZUnsupported(Coordinate coord) {
  try {
    coord.setZ(0.0);
    fail(coord.runtimeType.toString() + " does not support Z");
  } catch (expected) {}
  assertTrue(coord.z.isNaN);
  coord.z = 0.0; // field still public
  assertTrueMsg("z field not used", coord.getZ().isNaN); // but not used
}

/**
 * Confirm the z field is not supported by getZ and setZ.
 */
void checkMUnsupported(Coordinate coord) {
  try {
    coord.setM(0.0);
    fail(coord.runtimeType.toString() + " does not support M");
  } catch (expected) {}
  assertTrue(coord.getM().isNaN);
}

void checkEmpty(Geometry geom, String clz) {
  assertTrue(geom.isEmpty());
  assertTrue(geom.runtimeType.toString() == clz);
}

/**
 * CoordinateArraySequences default their dimension to 3 unless explicitly told otherwise.
 * This test ensures that GeometryFactory.createGeometry() recreates the input dimension properly.
 *
 * @throws ParseException
 */
void testCopyGeometryWithNonDefaultDimension() {
  GeometryFactory gf = new GeometryFactory.withCoordinateSequenceFactory(CoordinateArraySequenceFactory());
  CoordinateSequence mpSeq = gf.getCoordinateSequenceFactory().createSizeDim(1, 2);
  mpSeq.setOrdinate(0, 0, 50);
  mpSeq.setOrdinate(0, 1, -2);

  Point g = gf.createPointSeq(mpSeq);
  CoordinateSequence pSeq = (g.getGeometryN(0) as Point).getCoordinateSequence();
  assertEquals(2, pSeq.getDimension());

  Point g2 = geometryFactory.createGeometry(g) as Point;
  assertEquals(2, g2.getCoordinateSequence().getDimension());
}

void checkCreateGeometryExact(String wkt) {
  Geometry g = read(wkt);
  Geometry g2 = geometryFactory.createGeometry(g);
  assertTrue(g.equalsExactGeom(g2));
}
