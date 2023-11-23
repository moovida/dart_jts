import 'package:dart_jts/dart_jts.dart';
import 'package:test/test.dart';

import 'geometry_operation_validator.dart';

void main() {
  var rdr = WKTReader();
  List<Geometry?> getVWSimplifierResult(String wkt, double tolerance) {
    List<Geometry?> ioGeom = List.empty(growable: true);
    ioGeom.add(rdr.read(wkt));
    ioGeom.add(VWSimplifier.simplify(ioGeom[0]!, tolerance));
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

  void checkVWS(String wkt, double tolerance, String wktExpected) {
    Geometry? geom = rdr.read(wkt);
    Geometry? actual = VWSimplifier.simplify(geom!, tolerance);
    Geometry? expected = rdr.read(wktExpected);
    checkEqual(expected!, actual);
  }

  test('testEmptyPolygon', ()  {
    GeometryOperationValidator(getVWSimplifierResult(
            'POLYGON ((20 220, 40 220, 60 220, 80 220, 100 220, 120 220, 140 220, 140 180, 100 180, 60 180, 20 180, 20 220))',
            1.0))
        .test();
  });
  test('testPolygonNoReduction', ()  {
    GeometryOperationValidator(getVWSimplifierResult(
            'POLYGON ((20 220, 40 220, 60 220, 80 220, 100 220, 120 220, 140 220, 140 180, 100 180, 60 180, 20 180, 20 220))',
            1.0))
        .test();
  });
  test('testPolygonSpikeInShell', ()  {
    GeometryOperationValidator(getVWSimplifierResult(
            'POLYGON ((1721355.3 693015.146, 1721318.687 693046.251, 1721306.747 693063.038, 1721367.025 692978.29, 1721355.3 693015.146))',
            10.0))
        .setExpectedResult(
            'POLYGON ((1721355.3 693015.146, 1721318.687 693046.251, 1721367.025 692978.29, 1721355.3 693015.146))')
        .test();
  });
  test('testPolygonSpikeInHole', ()  {
    GeometryOperationValidator(getVWSimplifierResult(
            'POLYGON ((1721270 693090, 1721400 693090, 1721400 692960, 1721270 692960, 1721270 693090), (1721355.3 693015.146, 1721318.687 693046.251, 1721306.747 693063.038, 1721367.025 692978.29, 1721355.3 693015.146))',
            10.0))
        .setExpectedResult(
            'POLYGON ((1721270 693090, 1721400 693090, 1721400 692960, 1721270 692960, 1721270 693090), (1721355.3 693015.146, 1721318.687 693046.251, 1721367.025 692978.29, 1721355.3 693015.146))')
        .test();
  });
}
