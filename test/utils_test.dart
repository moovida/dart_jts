import 'dart:typed_data';

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
  group("ByteUtils - ", () {
    test('Float64', () {
      var byteOrder = Endian.little;
      List<int> buf = List.from([0, 0, 0, 0, 0, 0, 0, 0]);
      double d1 = 1.246e370;
      double d2 = -0.0000745687436849437;

      Byteutils.putFloat64(d1, buf, byteOrder);
      var float64 = Byteutils.getFloat64(buf, byteOrder);
      expect(float64, d1);

      Byteutils.putFloat64(d2, buf, byteOrder);
      float64 = Byteutils.getFloat64(buf, byteOrder);
      expect(float64, d2);
    });
    test('Int32', () {
      var byteOrder = Endian.little;
      List<int> buf = List.from([0, 0, 0, 0, 0, 0, 0, 0]);
      int i1 = 192233720;
      Byteutils.putInt32(i1, buf, byteOrder);
      var int32 = Byteutils.getInt32(buf, byteOrder);
      expect(int32, i1);
    });
  });

  group("ByteStream - ", () {
    test('ByteOrderDataInStream', () {
      List<int> buf = List.from([71, 80, 0, 1, 0, 0, 0, 0]);
      var din = ByteOrderDataInStream(buf);
      var b1 = din.readByte();
      var b2 = din.readByte();
      var b3 = din.readByte();
      var b4 = din.readByte();

      expect(b1, 71);
      expect(b2, 80);
      expect(b3, 0);
      expect(b4, 1);

      din.setOrder(Endian.little);

      var int32 = din.readInt();
      expect(int32, 0);
    });
  });
}
