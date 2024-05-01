import 'package:dart_jts/dart_jts.dart';
import 'package:test/expect.dart';
import 'package:test/scaffolding.dart';

import 'testing_utilities.dart';

void main() {
  var rdr = WKTReader();
  void checkEqualMsg(String msg, Geometry expected, Geometry actual) {
    Geometry actualNorm = actual.norm();
    Geometry expectedNorm = expected.norm();
    bool equal = actualNorm.equalsExactGeom(expectedNorm);
    if (!equal) {
      print(
          'FAIL - $msg: Expected = ${expectedNorm.toString()} -- Actual = ${actualNorm.toString()}');
    }
    expect(equal, true);
  }

  void checkEqual(Geometry expected, Geometry actual) {
    checkEqualMsg('', expected, actual);
  }

  void assertEqualsExactAndHasSameFactory(Geometry expected, Geometry actual) {
    checkEqual(expected, actual);
    assertTrue(expected.getFactory() == actual.getFactory());
  }

  void checkReducePointwise(String wkt, String wktExpected) {
    Geometry? g = rdr.read(wkt);
    Geometry? gExpected = rdr.read(wktExpected);
    PrecisionModel pm = PrecisionModel.fixedPrecision(1);
    Geometry? gReduce = GeometryPrecisionReducer.reducePointwise(g!, pm);
    assertEqualsExactAndHasSameFactory(gExpected!, gReduce!);
  }
  test('testLineWithCollapse', () {
    checkReducePointwise(
        "LINESTRING (0 0,  0.1 0,  1 0)",
        "LINESTRING (0 0,  0   0,  1 0)");
  });
  test('testLineDuplicatePointsPreserved',(){
    checkReducePointwise(
        "LINESTRING (0 0,  0.1 0,  0.1 0,  1 0, 1 0)",
        "LINESTRING (0 0,  0   0,  0   0,  1 0, 1 0)");
  });
  test('testLineFullCollapse',(){
    checkReducePointwise(
        "LINESTRING (0 0,  0.1 0)",
        "LINESTRING (0 0,  0   0)");
  });
  test('testPolygonFullCollapse',(){
    checkReducePointwise(
        "POLYGON ((0.1 0.3, 0.3 0.3, 0.3 0.1, 0.1 0.1, 0.1 0.3))",
        "POLYGON ((0 0, 0 0, 0 0, 0 0, 0 0))");
  });
  test('testPolygonWithCollapsedLine',(){
    checkReducePointwise(
        "POLYGON ((10 10, 100 100, 200 10.1, 300 10, 10 10))",
        "POLYGON ((10 10, 100 100, 200 10,   300 10, 10 10))");
  });
  test('testPolygonWithCollapsedPoint',(){
    checkReducePointwise(
        "POLYGON ((10 10, 100 100, 200 10.1, 300 100, 400 10, 10 10))",
        "POLYGON ((10 10, 100 100, 200 10,   300 100, 400 10, 10 10))");
  });
}
