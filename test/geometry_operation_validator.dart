import 'package:dart_jts/dart_jts.dart';
import 'package:test/test.dart';
import 'same_structure_tester.dart';

class GeometryOperationValidator {
  final WKTReader rdr = new WKTReader();
  List<Geometry?> ioGeometry;
  bool expectedSameStructure = false;
  String wktExpected = '';

  GeometryOperationValidator(this.ioGeometry);

  GeometryOperationValidator setExpectedResult(String wktExpected) {
    this.wktExpected = wktExpected;
    return this;
  }

  GeometryOperationValidator setExpectedSameStructure() {
    this.expectedSameStructure = true;
    return this;
  }

  bool isAllTestsPassed() {
    try {
      test();
    } catch (e) {
      return false;
    }
    return true;
  }

  ///
  /// Tests if the result is valid.
  /// Throws an exception if result is not valid.
  /// This allows chaining multiple tests together.
  ///
  /// @throws Exception if the result is not valid.
  ///
  void test() {
    testSameStructure();
    try {
      testValid();
    }catch(e){
      throw e;
    }
    testExpectedResult();
  }

  GeometryOperationValidator testSameStructure() {
    if (!expectedSameStructure) {
      return this;
    }
    expect(SameStructureTester.isSameStructure(ioGeometry[0]!, ioGeometry[1]!),
        true);
    return this;
  }

  GeometryOperationValidator testValid() {
    IsValidOp op = IsValidOp(ioGeometry[1]!);
    bool isValid = op.isValid();
    // if(!isValid){
    //   print(op.validErr!.getMessage());
    // }
    assert(isValid,'Simplified geometry is invalid');
    return this;
  }

  GeometryOperationValidator testEmpty(bool isEmpty) {
    String failureCondition = isEmpty ? "not empty" : "empty";
    Assert.isTrue(ioGeometry[1]!.isEmpty() == isEmpty,
        'simplified geometry is $failureCondition');
    return this;
  }

  void testExpectedResult() {
    if (wktExpected.isEmpty) {
      return;
    }
    Geometry? expectedGeom = rdr.read(wktExpected);
    bool isEqual = expectedGeom!.equalsExactGeom(ioGeometry[1]!);
    if (!isEqual) {
      print('Result not expected: ${ioGeometry[1]}');
    }
    Assert.isTrue(isEqual, 'Expected result not found');
  }
}
