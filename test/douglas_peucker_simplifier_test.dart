import 'package:test/test.dart';
import 'package:dart_jts/dart_jts.dart';

import 'geometry_operation_validator.dart';

void main() {
  var rdr = WKTReader();
  List<Geometry?> getSimplifierResult(String wkt, double tolerance) {
    List<Geometry?> ioGeom = List.empty(growable: true);
    ioGeom.add(rdr.read(wkt));
    ioGeom.add(DouglasPeuckerSimplifier.simplify(ioGeom[0]!, tolerance));
    return ioGeom;
  }

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

  void checkDP(String wkt, double tolerance, String wktExpected) {
    Geometry? geom = rdr.read(wkt);
    Geometry? result = DouglasPeuckerSimplifier.simplify(geom!, tolerance);
    Geometry? expected = rdr.read(wktExpected);
    checkEqual(expected!, result!);
  }

  test('testEmptyPolygon', () {
//    print('testEmptyPolygon');
    String geomStr = 'POLYGON(EMPTY)';
    GeometryOperationValidator(getSimplifierResult(geomStr, 1.0))
        .setExpectedResult(geomStr)
        .test();
  });
  test('testInvalidPolygon', () {
//    print('testInvalidPolygon');
    WKTReader rdr = new WKTReader();
    Geometry? geom = rdr.read(
        'POLYGON ((21.32686 47.78723, 21.31486 47.81023,21.32786 47.81123, '
        '21.33986 47.80223, 21.32586 47.82723,21.32786 47.82323, 21.33886 47.82623, '
        '21.34186 47.82123,21.40686 47.81723, 21.32686 47.78723))');
    IsValidOp op = new IsValidOp(geom!);
    if (!op.isValid()) {
//      print(op.getValidationError()!.getMessage());
      Geometry geomFixed = geom.buffer(0.0);
      Geometry? expected = rdr.read('POLYGON ((21.32686 47.78723, 21.31486 47.81023, ' +
          '21.32786 47.81123, 21.33986 47.80223, 21.328068201892744 47.823286782334385, ' +
          '21.33886 47.82623, 21.34186 47.82123, 21.40686 47.81723, 21.32686 47.78723))');
      checkEqual(expected!,geomFixed);
    }
  });

  test('testPoint', () {
//    print('testPoint');
    String geomStr = 'POINT (10 10)';
    GeometryOperationValidator(getSimplifierResult(geomStr, 1.0))
        .setExpectedResult(geomStr)
        .test();
  });
  test('testPolygonNoReduction', () {
//    print('testPolygonNoReduction');
    String geomStr =
        'POLYGON ((20 220, 40 220, 60 220, 80 220, 100 220, 120 220, 140 220, 140 180, 100 180, 60 180,     20 180, 20 220))';
    GeometryOperationValidator(getSimplifierResult(geomStr, 1.0)).test();
  });
  test('testPolygonReductionWithSplit', () {
//    print('testPolygonReductionWithSplit');
    String geomStr =
        'POLYGON ((40 240, 160 241, 280 240, 280 160, 160 240, 40 140, 40 240))';
    GeometryOperationValidator(getSimplifierResult(geomStr, 1.0)).test();
  });
  test('testPolygonReduction', () {
//    print('testPolygonReduction');
    GeometryOperationValidator(getSimplifierResult(
            'POLYGON ((120 120, 121 121, 122 122, 220 120, 180 199, 160 200, 140 199, 120 120))',
            10.0))
        .test();
  });
  test('testPolygonWithTouchingHole', () {
//    print('testPolygonWithTouchingHole');
    GeometryOperationValidator(getSimplifierResult(
            'POLYGON ((80 200, 240 200, 240 60, 80 60, 80 200), (120 120, 220 120, 180 199, 160 200, 140 199, 120 120))',
            10.0))
        .setExpectedResult(
            'POLYGON ((80 200, 240 200, 240 60, 80 60, 80 200), (120 120, 220 120, 180 199, 160 200, 140 199, 120 120))')
        .test();
  });
  test('testFlattishPolygon', () {
//    print('testFlattishPolygon');
    GeometryOperationValidator(getSimplifierResult(
            'POLYGON ((0 0, 50 0, 53 0, 55 0, 100 0, 70 1,  60 1, 50 1, 40 1, 0 0))',
            10.0))
        .setExpectedResult("POLYGON EMPTY")
        .test();
  });
  test('testTinySquare', () {
//    print('testTinySquare');
    GeometryOperationValidator(getSimplifierResult(
            'POLYGON ((0 5, 5 5, 5 0, 0 0, 0 1, 0 5))', 10.0))
        .test();
  });
  test('testTinyHole', () {
//    print('testTinySquare');
    GeometryOperationValidator(getSimplifierResult(
            'POLYGON ((10 10, 10 310, 370 310, 370 10, 10 10), (160 190, 180 190, 180 170, 160 190))',
            30.0))
        .testEmpty(false);
  });
  test('testTinyLineString', () {
//    print('testTinyLineString');
    GeometryOperationValidator(
            getSimplifierResult('LINESTRING (0 5, 1 5, 2 5, 5 5)', 10.0))
        .test();
  });
  test('testMultiPoint', () {
//    print('testMultiPoint');
    String geomStr =
        'MULTIPOINT(80 200, 240 200, 240 60, 80 60, 80 200, 140 199, 120 120)';
    GeometryOperationValidator(getSimplifierResult(geomStr, 10.0))
        .setExpectedResult(geomStr)
        .test();
  });
  test('testMultiLineString', () {
//    print('testMultiLineString');
    GeometryOperationValidator(getSimplifierResult(
            'MULTILINESTRING( (0 0, 50 0, 70 0, 80 0, 100 0), (0 0, 50 1, 60 1, 100 0) )',
            10.0))
        .test();
  });
  test('testMultiLineStringWithEmpty', () {
//    print('testMultiLineStringWithEmpty');
    GeometryOperationValidator(getSimplifierResult(
            'MULTILINESTRING( EMPTY, (0 0, 50 0, 70 0, 80 0, 100 0), (0 0, 50 1, 60 1, 100 0) )',
            10.0))
        .test();
  });
  test('testMultiPolygonWithEmpty', () {
//    print('testMultiPolygonWithEmpty');
    GeometryOperationValidator(getSimplifierResult(
            'MULTIPOLYGON (EMPTY, ((-36 91.5, 4.5 91.5, 4.5 57.5, -36 57.5, -36 91.5)), ((25.5 57.5, 61.5 57.5, 61.5 23.5, 25.5 23.5, 25.5 57.5)))',
            10.0))
        .test();
  });
  test('testGeometryCollection', () {
//    print('testGeometryCollection');
    GeometryOperationValidator(getSimplifierResult(
            'GEOMETRYCOLLECTION (' +
                'MULTIPOINT (80 200, 240 200, 240 60, 80 60, 80 200, 140 199, 120 120),' +
                'POLYGON ((80 200, 240 200, 240 60, 80 60, 80 200)),' +
                'LINESTRING (80 200, 240 200, 240 60, 80 60, 80 200, 140 199, 120 120)' +
                ')',
            10.0))
        .test();
  });
  test('testInvalidPolygonFixed', () {
//    print('testInvalidPolygonFixed');
    checkDP(
        'POLYGON ((21.32686 47.78723, 21.32386 47.79023, 21.32186 47.80223, 21.31486 47.81023, 21.32786 47.81123, 21.33986 47.80223, 21.33886 47.81123, 21.32686 47.82023, 21.32586 47.82723, 21.32786 47.82323, 21.33886 47.82623, 21.34186 47.82123, 21.36386 47.82223, 21.40686 47.81723, 21.32686 47.78723))',
        0.0036,
        'POLYGON ((21.32686 47.78723, 21.31486 47.81023, 21.32786 47.81123, 21.33986 47.80223, 21.328068201892744 47.823286782334385, 21.33886 47.82623, 21.34186 47.82123, 21.40686 47.81723, 21.32686 47.78723))');
  });
  test('testPolygonCollapseRemoved', () {
//    print('testPolygonCollapseRemoved');
    checkDP(
        'MULTIPOLYGON (((-76.02716827 36.55671692, -75.99866486 36.55665207, -75.91191864 36.54253006, -75.92480469 36.47397614, -75.97727966 36.4780159, -75.97628784 36.51792526, -76.02716827 36.55671692)), ((-75.90198517 36.55619812, -75.8781662 36.55587387, -75.77315521 36.22925568, -75.78317261 36.22519302, -75.90198517 36.55619812)))',
        0.05,
        'POLYGON ((-76.02716827 36.55671692, -75.91191864 36.54253006, -75.92480469 36.47397614, -76.02716827 36.55671692))');
  });
}
