import "package:test/test.dart";
import 'package:dart_jts/dart_jts.dart';

void main() {
  int A = Dimension.A;
  int L = Dimension.L;
  int P = Dimension.P;

  group("intersectionmatrix - ", () {
    test("test 1", () {
      IntersectionMatrix i = IntersectionMatrix();
      i.setDimensionSimbols("012*TF012");
      expect("012*TF012", i.toString());
      IntersectionMatrix c = IntersectionMatrix.fromMatrix(i);
      expect("012*TF012", c.toString());
    });
    test("test transpose", () {
      IntersectionMatrix x = IntersectionMatrix.fromDimensionSymbols("012*TF012");
      IntersectionMatrix i = IntersectionMatrix.fromMatrix(x);
      IntersectionMatrix j = i.transpose();
      expect(i, j);
      expect("0*01T12F2", i.toString());
      expect("012*TF012", x.toString());
    });

    test("testIsDisjoint", () {
      expect((IntersectionMatrix.fromDimensionSymbols("FF*FF****")).isDisjoint(), true);
      expect((IntersectionMatrix.fromDimensionSymbols("FF1FF2T*0")).isDisjoint(), true);
      expect(!(IntersectionMatrix.fromDimensionSymbols("*F*FF****")).isDisjoint(), true);
    });

    test("testIsTouches", () {
      expect((IntersectionMatrix.fromDimensionSymbols("FT*******")).isTouches(P, A), true);
      expect((IntersectionMatrix.fromDimensionSymbols("FT*******")).isTouches(A, P), true);
      expect(!(IntersectionMatrix.fromDimensionSymbols("FT*******")).isTouches(P, P), true);
    });
    test("testIsIntersects", () {
      expect(!(IntersectionMatrix.fromDimensionSymbols("FF*FF****")).isIntersects(), true);
      expect(!(IntersectionMatrix.fromDimensionSymbols("FF1FF2T*0")).isIntersects(), true);
      expect((IntersectionMatrix.fromDimensionSymbols("*F*FF****")).isIntersects(), true);
    });
    test("testIsIntersects", () {
      expect(!(IntersectionMatrix.fromDimensionSymbols("FF*FF****")).isIntersects(), true);
      expect(!(IntersectionMatrix.fromDimensionSymbols("FF1FF2T*0")).isIntersects(), true);
      expect((IntersectionMatrix.fromDimensionSymbols("*F*FF****")).isIntersects(), true);
    });
    test("testIsCrosses", () {
      expect((IntersectionMatrix.fromDimensionSymbols("TFTFFFFFF")).isCrosses(P, L), true);
      expect(!(IntersectionMatrix.fromDimensionSymbols("TFTFFFFFF")).isCrosses(L, P), true);
      expect(!(IntersectionMatrix.fromDimensionSymbols("TFFFFFTFF")).isCrosses(P, L), true);
      expect((IntersectionMatrix.fromDimensionSymbols("TFFFFFTFF")).isCrosses(L, P), true);
      expect((IntersectionMatrix.fromDimensionSymbols("0FFFFFFFF")).isCrosses(L, L), true);
      expect(!(IntersectionMatrix.fromDimensionSymbols("1FFFFFFFF")).isCrosses(L, L), true);
    });
    test("testIsWithin", () {
      expect((IntersectionMatrix.fromDimensionSymbols("T0F00F000")).isWithin(), true);
      expect(!(IntersectionMatrix.fromDimensionSymbols("T00000FF0")).isWithin(), true);
    });
    test("testIsContains", () {
      expect(!(IntersectionMatrix.fromDimensionSymbols("T0F00F000")).isContains(), true);
      expect((IntersectionMatrix.fromDimensionSymbols("T00000FF0")).isContains(), true);
    });
    test("testIsOverlaps", () {
      expect((IntersectionMatrix.fromDimensionSymbols("2*2***2**")).isOverlaps(P, P), true);
      expect((IntersectionMatrix.fromDimensionSymbols("2*2***2**")).isOverlaps(A, A), true);
      expect(!(IntersectionMatrix.fromDimensionSymbols("2*2***2**")).isOverlaps(P, A), true);
      expect(!(IntersectionMatrix.fromDimensionSymbols("2*2***2**")).isOverlaps(L, L), true);
      expect((IntersectionMatrix.fromDimensionSymbols("1*2***2**")).isOverlaps(L, L), true);

      expect(!(IntersectionMatrix.fromDimensionSymbols("0FFFFFFF2")).isOverlaps(P, P), true);
      expect(!(IntersectionMatrix.fromDimensionSymbols("1FFF0FFF2")).isOverlaps(L, L), true);
      expect(!(IntersectionMatrix.fromDimensionSymbols("2FFF1FFF2")).isOverlaps(A, A), true);
    });
    test("testIsEquals", () {
      expect((IntersectionMatrix.fromDimensionSymbols("0FFFFFFF2")).isEquals(P, P), true);
      expect((IntersectionMatrix.fromDimensionSymbols("1FFF0FFF2")).isEquals(L, L), true);
      expect((IntersectionMatrix.fromDimensionSymbols("2FFF1FFF2")).isEquals(A, A), true);

      expect(!(IntersectionMatrix.fromDimensionSymbols("0F0FFFFF2")).isEquals(P, P), true);
      expect((IntersectionMatrix.fromDimensionSymbols("1FFF1FFF2")).isEquals(L, L), true);
      expect(!(IntersectionMatrix.fromDimensionSymbols("2FFF1*FF2")).isEquals(A, A), true);

      expect(!(IntersectionMatrix.fromDimensionSymbols("0FFFFFFF2")).isEquals(P, L), true);
      expect(!(IntersectionMatrix.fromDimensionSymbols("1FFF0FFF2")).isEquals(L, A), true);
      expect(!(IntersectionMatrix.fromDimensionSymbols("2FFF1FFF2")).isEquals(A, P), true);
    });
  });
}
