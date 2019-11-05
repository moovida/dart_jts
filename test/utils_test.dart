import "package:test/test.dart";
import 'package:dart_jts/dart_jts.dart';

void main() {
  group("numeric utils - ", () {
    test("intgrid", () {
      var intMatrix = MatrixUtils.createIntMatrix(3, 5);

      expect(intMatrix.length, 3);
      expect(intMatrix[0].length, 5);
    });
  });

  group("collection utils - ", () {
    test("shift list", () {
      var list = ["a", "b", "c", 1, 2, 3];
      var expectedList1 = [1, 2, 3, "a", "b", "c"];
      var shiftedList = CollectionsUtils.shiftToFirst(list, 3);

      expect(expectedList1, shiftedList);
    });
  });

  group("StringUtils - ", () {
    test("replaceCharAt", () {
      String string = "abcdef";
      String expectedString = "abadef";
      var newString = StringUtils.replaceCharAt(string, 2, 'a');

      expect(expectedString, newString);
    });
  });
}
