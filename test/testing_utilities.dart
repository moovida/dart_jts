import 'package:dart_jts/dart_jts.dart';
import 'package:test/test.dart';

const String WKT_POINT = "POINT ( 10 10)";

const String WKT_LINESTRING = "LINESTRING (10 10, 20 20, 30 40)";

const String WKT_LINEARRING = "LINEARRING (10 10, 20 20, 30 40, 10 10)";

const String WKT_POLY = "POLYGON ((50 50, 50 150, 150 150, 150 50, 50 50))";

const String WKT_MULTIPOINT = "MULTIPOINT ((10 10), (20 20))";

const String WKT_MULTILINESTRING = "MULTILINESTRING ((10 10, 20 20), (15 15, 30 15))";

const String WKT_MULTIPOLYGON = "MULTIPOLYGON (((10 10, 10 20, 20 20, 20 15, 10 10)), ((60 60, 70 70, 80 60, 60 60)))";

const String WKT_GC = "GEOMETRYCOLLECTION (POLYGON ((100 200, 200 200, 200 100, 100 100, 100 200)), LINESTRING (150 250, 250 250))";

assertEquals(actual, matcher) {
  if (actual is double && matcher is double) {
    if (actual.isNaN && matcher.isNaN) return;
  }
  if (actual is Coordinate && matcher is Coordinate) {
    if (actual.equals(matcher)) return;
  }
  expect(actual, matcher);
}

assertEqualsD(double n1, double n2, double tolerance) {
  expect(NumberUtils.equalsWithTolerance(n1, n2, tolerance), true);
}

assertTrue(actual) {
  expect(actual, true);
}

assertTrueMsg(String msg, actual) {
  expect(actual, true, reason: msg);
}

GeometryFactory geomFactory = GeometryFactory.withCoordinateSequenceFactory(CoordinateArraySequenceFactory());
WKTReader readerWKT = WKTReader.withFactory(geomFactory);

Geometry read(String wkt) {
  //return read(readerWKT, wkt);
  return WKTReader.withFactory(geomFactory).read(wkt);
}
