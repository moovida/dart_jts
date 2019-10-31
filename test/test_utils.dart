import "package:test/test.dart";
import 'package:dart_sfs/dart_sfs.dart';

void main() {
  group("numeric utils - ", () {
    test("intgrid", () {
      var intMatrix = MatrixUtils.createIntMatrix(3, 5);

      expect(intMatrix.length, 3);
      expect(intMatrix[0].length, 5);
    });
  });
}
