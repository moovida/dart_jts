import "package:test/test.dart";
import 'package:dart_sfs/dart_sfs.dart';

PrecisionModel precisionModel = new PrecisionModel.fixedPrecision(1);

GeometryFactory geometryFactory = new GeometryFactory.withPrecisionModelSrid(precisionModel, 0);

WKTReader reader = new WKTReader.withFactory(geometryFactory);

void main() {
  group("envelope - ", () {
    test("testEverything", () {
      Envelope e1 = Envelope.empty();
      expect(e1.isNull(), true);
      expect(0, e1.getWidth());
      expect(0, e1.getHeight());
      e1.expandToInclude(100, 101);
      e1.expandToInclude(200, 202);
      e1.expandToInclude(150, 151);
      expect(200, e1.getMaxX());
      expect(202, e1.getMaxY());
      expect(100, e1.getMinX());
      expect(101, e1.getMinY());
      expect(e1.contains(120, 120), true);
      expect(e1.contains(120, 101), true);
      expect(!e1.contains(120, 100), true);
      expect(101, e1.getHeight());
      expect(100, e1.getWidth());
      expect(!e1.isNull(), true);

      Envelope e2 = Envelope(499, 500, 500, 501);
      expect(!e1.containsEnvelope(e2), true);
      expect(!e1.intersectsEnvelope(e2), true);
      e1.expandToIncludeEnvelope(e2);
      expect(e1.containsEnvelope(e2), true);
      expect(e1.intersectsEnvelope(e2), true);
      expect(500, e1.getMaxX());
      expect(501, e1.getMaxY());
      expect(100, e1.getMinX());
      expect(101, e1.getMinY());

      Envelope e3 = Envelope(300, 700, 300, 700);
      expect(!e1.containsEnvelope(e3), true);
      expect(e1.intersectsEnvelope(e3), true);

      Envelope e4 = Envelope(300, 301, 300, 301);
      expect(e1.containsEnvelope(e4), true);
      expect(e1.intersectsEnvelope(e4), true);
    });
    test("testIntersects", () {
      checkIntersectsPermuted(1, 1, 2, 2, 2, 2, 3, 3, true);
      checkIntersectsPermuted(1, 1, 2, 2, 3, 3, 4, 4, false);
    });
    test("testIntersectsEmpty", () {
      expect(!Envelope(-5, 5, -5, 5).intersectsEnvelope(Envelope.empty()), true);
      expect(!Envelope.empty().intersectsEnvelope(Envelope(-5, 5, -5, 5)), true);
      expect(!Envelope.empty().intersectsEnvelope(Envelope(100, 101, 100, 101)), true);
      expect(!Envelope(100, 101, 100, 101).intersectsEnvelope(Envelope.empty()), true);
    });
    test("testDisjointEmpty", () {
      expect(Envelope(-5, 5, -5, 5).disjoint(Envelope.empty()), true);
      expect(Envelope.empty().disjoint(Envelope(-5, 5, -5, 5)), true);
      expect(Envelope.empty().disjoint(Envelope(100, 101, 100, 101)), true);
      expect(Envelope(100, 101, 100, 101).disjoint(Envelope.empty()), true);
    });
    test("testContainsEmpty", () {
      expect(!Envelope(-5, 5, -5, 5).containsEnvelope(Envelope.empty()), true);
      expect(!Envelope.empty().containsEnvelope(Envelope(-5, 5, -5, 5)), true);
      expect(!Envelope.empty().containsEnvelope(Envelope(100, 101, 100, 101)), true);
      expect(!Envelope(100, 101, 100, 101).containsEnvelope(Envelope.empty()), true);
    });
    test("testExpandToIncludeEmpty", () {
      expect(Envelope(-5, 5, -5, 5), expandToInclude(Envelope(-5, 5, -5, 5), Envelope.empty()));
      expect(Envelope(-5, 5, -5, 5), expandToInclude(Envelope.empty(), Envelope(-5, 5, -5, 5)));
      expect(Envelope(100, 101, 100, 101), expandToInclude(Envelope.empty(), Envelope(100, 101, 100, 101)));
      expect(Envelope(100, 101, 100, 101), expandToInclude(Envelope(100, 101, 100, 101), Envelope.empty()));
    });
    test("testEmpty", () {
      expect(0, Envelope.empty().getHeight());
      expect(0, Envelope.empty().getWidth());
      expect(Envelope.empty(), Envelope.empty());
      Envelope e = Envelope(100, 101, 100, 101);
      e.initFromEnvelope(Envelope.empty());
      expect(Envelope.empty(), e);
    });
    test("testAsGeometry", () {
      expect(geometryFactory.createPoint(null).getEnvelope().isEmpty(), true);

      Geometry g = geometryFactory.createPoint(Coordinate(5, 6)).getEnvelope();
      expect(!g.isEmpty(), true);
      expect(g is Point, true);

      Point p = g as Point;
      expect(NumberUtils.equalsWithTolerance(5, p.getX(), 1E-1), true);
      expect(NumberUtils.equalsWithTolerance(6, p.getY(), 1E-1), true);

      LineString l = reader.read("LINESTRING(10 10, 20 20, 30 40)") as LineString;
      Geometry g2 = l.getEnvelope();
      expect(!g2.isEmpty(), true);
      expect(g2 is Polygon, true);

      Polygon poly = g2 as Polygon;
      poly.normalize();
      expect(5, poly.getExteriorRing().getNumPoints());
      expectCoords(Coordinate(10, 10), poly.getExteriorRing().getCoordinateN(0));
      expectCoords(Coordinate(10, 40), poly.getExteriorRing().getCoordinateN(1));
      expectCoords(Coordinate(30, 40), poly.getExteriorRing().getCoordinateN(2));
      expectCoords(Coordinate(30, 10), poly.getExteriorRing().getCoordinateN(3));
      expectCoords(Coordinate(10, 10), poly.getExteriorRing().getCoordinateN(4));
    });

    test("testSetToNull", () {
      Envelope e1 = Envelope.empty();
      expect(e1.isNull(), true);
      e1.expandToInclude(5, 5);
      expect(!e1.isNull(), true);
      e1.setToNull();
      expect(e1.isNull(), true);
    });
    test("testEquals", () {
      Envelope e1 = Envelope(1, 2, 3, 4);
      Envelope e2 = Envelope(1, 2, 3, 4);
      expect(e1, e2);
      expect(e1.hashCode, e2.hashCode);

      Envelope e3 = Envelope(1, 2, 3, 5);
      expect(e1 != e3, true);
      expect(e1.hashCode != e3.hashCode, true);
      e1.setToNull();
      expect(e1 != e2, true);
      expect(e1.hashCode != e2.hashCode, true);
      e2.setToNull();
      expect(e1, e2);
      expect(e1.hashCode, e2.hashCode);
    });
    test("testEquals2", () {
      expect(Envelope.empty() == Envelope.empty(), true);
      expect(Envelope(1, 2, 1, 2) == Envelope(1, 2, 1, 2), true);
      expect(Envelope(1, 2, 1.5, 2) != Envelope(1, 2, 1, 2), true);
    });
    test("testCopyConstructor", () {
      Envelope e1 = Envelope(1, 2, 3, 4);
      Envelope e2 = Envelope.fromEnvelope(e1);
      expect(1, e2.getMinX());
      expect(2, e2.getMaxX());
      expect(3, e2.getMinY());
      expect(4, e2.getMaxY());
    });
    test("testCopy", () {
      Envelope e1 = Envelope(1, 2, 3, 4);
      Envelope e2 = e1.copy();
      expect(1, e2.getMinX());
      expect(2, e2.getMaxX());
      expect(3, e2.getMinY());
      expect(4, e2.getMaxY());

      Envelope eNull = Envelope.empty();
      Envelope eNullCopy = eNull.copy();
      expect(eNullCopy.isNull(), true);
    });
    test("testGeometryFactoryCreateEnvelope", () {
      checkExpectedEnvelopeGeometry("POINT (0 0)");
      checkExpectedEnvelopeGeometry("POINT (100 13)");
      checkExpectedEnvelopeGeometry("LINESTRING (0 0, 0 10)");
      checkExpectedEnvelopeGeometry("LINESTRING (0 0, 10 0)");

      String poly10 = "POLYGON ((0 10, 10 10, 10 0, 0 0, 0 10))";
      checkExpectedEnvelopeGeometry(poly10);

      checkExpectedEnvelopeGeometry2("LINESTRING (0 0, 10 10)", poly10);
      checkExpectedEnvelopeGeometry2("POLYGON ((5 10, 10 6, 5 0, 0 6, 5 10))", poly10);
    });
    test("testMetrics", () {
      Envelope env = Envelope(0, 4, 0, 3);
      expect(env.getWidth(), 4.0);
      expect(env.getHeight(), 3.0);
      expect(env.getDiameter(), 5.0);
    });
    test("testEmptyMetrics", () {
      Envelope env = Envelope.empty();
      expect(env.getWidth(), 0.0);
      expect(env.getHeight(), 0.0);
      expect(env.getDiameter(), 0.0);
    });
    test("testCompareTo", () {
      checkCompareTo(0, Envelope.empty(), Envelope.empty());
      checkCompareTo(0, Envelope(1, 2, 1, 2), Envelope(1, 2, 1, 2));
      checkCompareTo(1, Envelope(2, 3, 1, 2), Envelope(1, 2, 1, 2));
      checkCompareTo(-1, Envelope(1, 2, 1, 2), Envelope(2, 3, 1, 2));
      checkCompareTo(1, Envelope(1, 2, 1, 3), Envelope(1, 2, 1, 2));
      checkCompareTo(1, Envelope(2, 3, 1, 3), Envelope(1, 3, 1, 2));
    });
    test("empty", () {});
  });
}

void expectCoords(Coordinate c1, Coordinate c2) {
  expect(c1.equals3D(c2), true);
}

void checkIntersectsPermuted(double a1x, double a1y, double a2x, double a2y, double b1x, double b1y, double b2x, double b2y, bool expected) {
  checkIntersects(a1x, a1y, a2x, a2y, b1x, b1y, b2x, b2y, expected);
  checkIntersects(a1x, a2y, a2x, a1y, b1x, b1y, b2x, b2y, expected);
  checkIntersects(a1x, a1y, a2x, a2y, b1x, b2y, b2x, b1y, expected);
  checkIntersects(a1x, a2y, a2x, a1y, b1x, b2y, b2x, b1y, expected);
}

void checkIntersects(double a1x, double a1y, double a2x, double a2y, double b1x, double b1y, double b2x, double b2y, bool expected) {
  Envelope a = Envelope(a1x, a2x, a1y, a2y);
  Envelope b = Envelope(b1x, b2x, b1y, b2y);
  expect(expected, a.intersectsEnvelope(b));
  expect(expected, !a.disjoint(b));

  Coordinate a1 = Coordinate(a1x, a1y);
  Coordinate a2 = Coordinate(a2x, a2y);
  Coordinate b1 = Coordinate(b1x, b1y);
  Coordinate b2 = Coordinate(b2x, b2y);
  expect(expected, Envelope.intersectsEnvelopeCoords(a1, a2, b1, b2));

  expect(expected, a.intersectsEnvelopeCoordinates(b1, b2));
}

void checkExpectedEnvelopeGeometry(String wktInput) {
  checkExpectedEnvelopeGeometry2(wktInput, wktInput);
}

void checkExpectedEnvelopeGeometry2(String wktInput, String wktEnvGeomExpected) {
//  Geometry input = WKTReader(wktInput).parseGeometryCollection();
//  Geometry envGeomExpected = WKTReader(wktEnvGeomExpected).parseGeometryCollection();
//
//  Envelope env = input.envelope;
//  Geometry envGeomActual = input.getFactory().toGeometry(env);
//  bool isEqual = envGeomActual.equalsNorm(envGeomExpected);
//  expect(isEqual, true);
}

void checkCompareTo(int expected, Envelope env1, Envelope env2) {
  expect(expected == env1.compareTo(env2), true);
  expect(-expected == env2.compareTo(env1), true);
}

Envelope expandToInclude(Envelope a, Envelope b) {
  a.expandToIncludeEnvelope(b);
  return a;
}
