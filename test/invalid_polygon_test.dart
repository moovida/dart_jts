import 'package:test/test.dart';
import 'package:dart_jts/dart_jts.dart';
void main() {
  void checkEqualMsg(String msg, Geometry expected, Geometry actual) {
    // Geometry actualNorm = actual.norm();
    // Geometry expectedNorm = expected.norm();
    bool equal = actual.equalsExactGeom(expected);
    if (!equal) {
      print(
          'FAIL - $msg: Expected = ${expected.toString()} -- Actual = ${actual.toString()}');
    }
    expect(equal, true);
  }

  test('testInvalidPolygon', () {
    WKTReader rdr = WKTReader();
    Geometry? geom = rdr.read(
        'POLYGON ((21.32686 47.78723, 21.31486 47.81023,21.32786 47.81123, '
            '21.33986 47.80223, 21.32586 47.82723,21.32786 47.82323, 21.33886 47.82623, '
            '21.34186 47.82123,21.40686 47.81723, 21.32686 47.78723))');
    IsValidOp op = IsValidOp(geom!);
    if (!op.isValid()) {
      print(op.getValidationError()!.getMessage());
      Geometry geomFixed = geom.buffer(0.0);
      var op1 = IsValidOp(geomFixed);
      if(!op1.isValid()){
        print(op1.getValidationError()!.getMessage());
      }
      Geometry? expected = rdr.read('POLYGON ((21.32686 47.78723, 21.31486 47.81023, ' +
          '21.32786 47.81123, 21.33986 47.80223, 21.328068201892744 47.823286782334385, ' +
          '21.33886 47.82623, 21.34186 47.82123, 21.40686 47.81723, 21.32686 47.78723))');
      checkEqualMsg('Geometries are not equal',expected!,geomFixed);
    }
  });

}