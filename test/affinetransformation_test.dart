import "package:test/test.dart";
import 'package:dart_jts/dart_jts.dart';
import "dart:math" as math;
import 'testing_utilities.dart';

const double TOLERANCE = 1e-10;

void main() {
  group("Affinetransformation Tests", () {
    test("testRotate1", () {
      AffineTransformation t =
          AffineTransformation.rotationInstance(math.pi / 2);
      checkTransformation(10, 0, t, 0, 10);
      checkTransformation(0, 10, t, -10, 0);
      checkTransformation(-10, -10, t, 10, -10);
    });
    test("testRotate2", () {
      AffineTransformation t =
          AffineTransformation.rotationInstanceSinCos(1, 0);
      checkTransformation(10, 0, t, 0, 10);
      checkTransformation(0, 10, t, -10, 0);
      checkTransformation(-10, -10, t, 10, -10);
    });
    test("testRotateAroundPoint1", () {
      AffineTransformation t =
          AffineTransformation.rotationInstanceTXY(math.pi / 2, 1, 1);
      checkTransformation(1, 1, t, 1, 1);
      checkTransformation(10, 0, t, 2, 10);
      checkTransformation(0, 10, t, -8, 0);
      checkTransformation(-10, -10, t, 12, -10);
    });
    test("testRotateAroundPoint2", () {
      AffineTransformation t =
          AffineTransformation.rotationInstanceSinCosXY(1, 0, 1, 1);
      checkTransformation(1, 1, t, 1, 1);
      checkTransformation(10, 0, t, 2, 10);
      checkTransformation(0, 10, t, -8, 0);
      checkTransformation(-10, -10, t, 12, -10);
    });
    test("testReflectXY1", () {
      AffineTransformation t = AffineTransformation.reflectionInstanceXY(1, 1);
      checkTransformation(10, 0, t, 0, 10);
      checkTransformation(0, 10, t, 10, 0);
      checkTransformation(-10, -10, t, -10, -10);
      checkTransformation(-3, -4, t, -4, -3);
    });
    test("testReflectXY2", () {
      AffineTransformation t = AffineTransformation.reflectionInstanceXY(1, -1);
      checkTransformation(10, 0, t, 0, -10);
      checkTransformation(0, 10, t, -10, 0);
      checkTransformation(-10, -10, t, 10, 10);
      checkTransformation(-3, -4, t, 4, 3);
    });
    test("testReflectXYXY1", () {
      AffineTransformation t =
          AffineTransformation.reflectionInstance(0, 5, 5, 0);
      checkTransformation(5, 0, t, 5, 0);
      checkTransformation(0, 0, t, 5, 5);
      checkTransformation(-10, -10, t, 15, 15);
    });
    test("testScale1", () {
      AffineTransformation t = AffineTransformation.scaleInstance(2, 3);
      checkTransformation(10, 0, t, 20, 0);
      checkTransformation(0, 10, t, 0, 30);
      checkTransformation(-10, -10, t, -20, -30);
    });
    test("testShear1", () {
      AffineTransformation t = AffineTransformation.shearInstance(2, 3);
      checkTransformation(10, 0, t, 10, 30);
    });
    test("testTranslate1", () {
      AffineTransformation t = AffineTransformation.translationInstance(2, 3);
      checkTransformation(1, 0, t, 3, 3);
      checkTransformation(0, 0, t, 2, 3);
      checkTransformation(-10, -5, t, -8, -2);
    });
    test("testTranslateRotate1", () {
      AffineTransformation t =
          AffineTransformation.translationInstance(3, 3).rotate(math.pi / 2);
      checkTransformation(10, 0, t, -3, 13);
      checkTransformation(-10, -10, t, 7, -7);
    });
    test("testCompose1", () {
      AffineTransformation t0 = AffineTransformation.translationInstance(10, 0);
      t0.rotate(math.pi / 2);
      t0.translate(0, -10);

      AffineTransformation t1 = AffineTransformation.translationInstance(0, 0);
      t1.rotate(math.pi / 2);

      checkTransformation2Trans(t0, t1);
    });
    test("testCompose2", () {
      AffineTransformation t0 =
          AffineTransformation.reflectionInstance(0, 0, 1, 0);
      t0.reflect(0, 0, 0, -1);

      AffineTransformation t1 = AffineTransformation.rotationInstance(math.pi);

      checkTransformation2Trans(t0, t1);
    });
    test("testComposeRotation1", () {
      AffineTransformation t0 =
          AffineTransformation.rotationInstanceTXY(1, 10, 10);

      AffineTransformation t1 =
          AffineTransformation.translationInstance(-10, -10);
      t1.rotate(1);
      t1.translate(10, 10);

      checkTransformation2Trans(t0, t1);
    });
    test("testLineString", () {
      checkTransformationGeomStr("LINESTRING (1 2, 10 20, 100 200)");
    });
    test("testPolygon", () {
      checkTransformationGeomStr("POLYGON ((0 0, 100 0, 100 100, 0 100, 0 0))");
    });
    test("testPolygonWthHole", () {
      checkTransformationGeomStr(
          "POLYGON ((0 0, 100 0, 100 100, 0 100, 0 0), (1 1, 1 10, 10 10, 10 1, 1 1) )");
    });
    test("testMultiPoint", () {
      checkTransformationGeomStr("MULTIPOINT (0 0, 1 4, 100 200)");
    });
    test("testMultiLineString", () {
      checkTransformationGeomStr(
          "MULTILINESTRING ((0 0, 1 10), (10 10, 20 30), (123 123, 456 789))");
    });
    test("testMultiPolygon", () {
      checkTransformationGeomStr(
          "MULTIPOLYGON ( ((0 0, 100 0, 100 100, 0 100, 0 0), (1 1, 1 10, 10 10, 10 1, 1 1) ), ((200 200, 200 250, 250 250, 250 200, 200 200)) )");
    });
    test("testGeometryCollection", () {
      checkTransformationGeomStr(
          "GEOMETRYCOLLECTION ( POINT ( 1 1), LINESTRING (0 0, 10 10), POLYGON ((0 0, 100 0, 100 100, 0 100, 0 0)) )");
    });
    test("testNestedGeometryCollection", () {
      checkTransformationGeomStr(
          "GEOMETRYCOLLECTION ( POINT (20 20), GEOMETRYCOLLECTION ( POINT ( 1 1), LINESTRING (0 0, 10 10), POLYGON ((0 0, 100 0, 100 100, 0 100, 0 0)) ) )");
    });
  });

  group("AffinetransformationBuiler Tests", () {
    test("testRotate1", () {
      AffineTransformation t =
          AffineTransformation.rotationInstance(math.pi / 2);
      checkTransformation(10, 0, t, 0, 10);
      checkTransformation(0, 10, t, -10, 0);
      checkTransformation(-10, -10, t, 10, -10);
    });

    test("testRotate1", () {
      // run3(0, 0,    1, 0,    0, 1,
      //     0, 0,    0, 1,    -1, 0);
      run3(0, 0, 1, 0, 0, 1, 0, 0, 0, 1, -1, 0);
    });

    test("testRotate2", () {
      // run3(0, 0,    1, 0,    0, 1,
      //     0, 0,    1, 1,    -1, 1);
      run3(0, 0, 1, 0, 0, 1, 0, 0, 1, 1, -1, 1);
    });

    test("testScale1", () {
      // run3(0, 0,    1, 0,    0, 1,
      //     0, 0,    2, 0,    0, 2);
      run3(0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 2);
    });

    test("testTranslate1", () {
      // run3(0, 0,    1, 0,    0, 1,
      //     5, 6,    6, 6,    5, 7);
      run3(0, 0, 1, 0, 0, 1, 5, 6, 6, 6, 5, 7);
    });

    test("testLinear1", () {
      // run3(0, 0,    1, 0,    0, 1,
      //     0, 0,    0, 0,    5, 7);
      run3(0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 5, 7);
    });

    test("testSingular2", () {
      // points on a line mapping to collinear points - not uniquely specified
      // runSingular(0, 0,    1,   1,     2, 2,
      //             0, 0,    10, 10,    30, 30);
      runSingular(0, 0, 1, 1, 2, 2, 0, 0, 10, 10, 30, 30);
    });

    test("testSingular3", () {
      // points on a line mapping to collinear points - not uniquely specified
      // runSingular(0, 0,    1,   1,     2, 2,
      //             0, 0,    10, 10,    20, 20);
      runSingular(0, 0, 1, 1, 2, 2, 0, 0, 10, 10, 20, 20);
    });

    test("testSingular1", () {
      // points on a line mapping to non-collinear points - no solution
      // runSingular(0, 0,    1, 1,    2, 2,
      //             0, 0,    1, 2,    1, 3);
      runSingular(0, 0, 1, 1, 2, 2, 0, 0, 1, 2, 1, 3);
    });

    test("testSingleControl1", () {
      // run(0, 0,
      //     5, 6);
      run(0, 0, 5, 6);
    });

    test("testDualControl_Translation", () {
      // run2(0, 0,    1, 1,
      // 		5, 5,    6, 6);
      run2(0, 0, 1, 1, 5, 5, 6, 6);
    });

    test("testDualControl_General", () {
      // run2(0, 0,    1, 1,
      // 		5, 5,    6, 9);
      run2(0, 0, 1, 1, 5, 5, 6, 9);
    });
  });
}

/**
   * Checks that a transformation produces the expected result
   * @param x the input pt x
   * @param y the input pt y
   * @param trans the transformation
   * @param xp the expected output x
   * @param yp the expected output y
   */
void checkTransformation(
    double x, double y, AffineTransformation trans, double xp, double yp) {
  Coordinate p = new Coordinate(x, y);
  Coordinate p2 = new Coordinate.empty2D();
  trans.transform(p, p2);
  assertEqualsD(xp, p2.x, .00005);
  assertEqualsD(yp, p2.y, .00005);

  // if the transformation is invertible, test the inverse

  try {
    AffineTransformation invTrans = trans.getInverse();
    Coordinate pInv = new Coordinate.empty2D();
    invTrans.transform(p2, pInv);
    assertEqualsD(x, pInv.x, .00005);
    assertEqualsD(y, pInv.y, .00005);

    double det = trans.getDeterminant();
    double detInv = invTrans.getDeterminant();
    assertEqualsD(det, 1.0 / detInv, .00005);
  } on Exception catch (e) {}
}

WKTReader rdr = new WKTReader();

void checkTransformationGeomStr(String geomStr) {
  Geometry geom = rdr.read(geomStr);
  AffineTransformation trans =
      AffineTransformation.rotationInstance(math.pi / 2);
  AffineTransformation inv = trans.getInverse();
  Geometry transGeom = geom.copy();
  transGeom.applyCSF(trans);
  // System.out.println(transGeom);
  transGeom.applyCSF(inv);
  // check if transformed geometry is equal to original
  bool isEqual = geom.equalsExactWithTol(transGeom, 0.0005);
  assertTrue(isEqual);
}

void checkTransformation2Trans(
    AffineTransformation trans0, AffineTransformation trans1) {
  List<double> m0 = trans0.getMatrixEntries();
  List<double> m1 = trans1.getMatrixEntries();
  for (int i = 0; i < m0.length; i++) {
    assertEqualsD(m0[i], m1[i], 0.000005);
  }
}

void run3(
    double p0x,
    double p0y,
    double p1x,
    double p1y,
    double p2x,
    double p2y,
    double pp0x,
    double pp0y,
    double pp1x,
    double pp1y,
    double pp2x,
    double pp2y) {
  Coordinate p0 = new Coordinate(p0x, p0y);
  Coordinate p1 = new Coordinate(p1x, p1y);
  Coordinate p2 = new Coordinate(p2x, p2y);

  Coordinate pp0 = new Coordinate(pp0x, pp0y);
  Coordinate pp1 = new Coordinate(pp1x, pp1y);
  Coordinate pp2 = new Coordinate(pp2x, pp2y);

  AffineTransformationBuilder atb =
      new AffineTransformationBuilder(p0, p1, p2, pp0, pp1, pp2);
  AffineTransformation trans = atb.getTransformation();

  Coordinate dest = new Coordinate.empty2D();
  assertEqualPoint(pp0, trans.transform(p0, dest));
  assertEqualPoint(pp1, trans.transform(p1, dest));
  assertEqualPoint(pp2, trans.transform(p2, dest));
}

void run2(double p0x, double p0y, double p1x, double p1y, double pp0x,
    double pp0y, double pp1x, double pp1y) {
  Coordinate p0 = new Coordinate(p0x, p0y);
  Coordinate p1 = new Coordinate(p1x, p1y);

  Coordinate pp0 = new Coordinate(pp0x, pp0y);
  Coordinate pp1 = new Coordinate(pp1x, pp1y);

  AffineTransformation trans =
      AffineTransformationFactory.createFromControlVectors2(p0, p1, pp0, pp1);

  Coordinate dest = new Coordinate.empty2D();
  assertEqualPoint(pp0, trans.transform(p0, dest));
  assertEqualPoint(pp1, trans.transform(p1, dest));
}

void run(double p0x, double p0y, double pp0x, double pp0y) {
  Coordinate p0 = new Coordinate(p0x, p0y);

  Coordinate pp0 = new Coordinate(pp0x, pp0y);

  AffineTransformation trans =
      AffineTransformationFactory.createFromControlVectors(p0, pp0);

  Coordinate dest = new Coordinate.empty2D();
  assertEqualPoint(pp0, trans.transform(p0, dest));
}

void runSingular(
    double p0x,
    double p0y,
    double p1x,
    double p1y,
    double p2x,
    double p2y,
    double pp0x,
    double pp0y,
    double pp1x,
    double pp1y,
    double pp2x,
    double pp2y) {
  Coordinate p0 = new Coordinate(p0x, p0y);
  Coordinate p1 = new Coordinate(p1x, p1y);
  Coordinate p2 = new Coordinate(p2x, p2y);

  Coordinate pp0 = new Coordinate(pp0x, pp0y);
  Coordinate pp1 = new Coordinate(pp1x, pp1y);
  Coordinate pp2 = new Coordinate(pp2x, pp2y);

  AffineTransformationBuilder atb =
      new AffineTransformationBuilder(p0, p1, p2, pp0, pp1, pp2);
  AffineTransformation trans = atb.getTransformation();
  assertEquals(trans, null);
}

Coordinate ctl0 = new Coordinate(-10, -10);
Coordinate ctl1 = new Coordinate(10, 20);
Coordinate ctl2 = new Coordinate(10, -20);

void testTransform1() {
  AffineTransformation trans = new AffineTransformation();
  trans.rotate(1);
  trans.translate(10, 10);
  trans.scale(2, 2);
  runTransform(trans, ctl0, ctl1, ctl2);
}

void testTransform2() {
  AffineTransformation trans = new AffineTransformation();
  trans.rotate(3);
  trans.translate(10, 10);
  trans.scale(2, 10);
  trans.shear(5, 2);
  trans.reflect(5, 8, 10, 2);
  runTransform(trans, ctl0, ctl1, ctl2);
}

void runTransform(
    AffineTransformation trans, Coordinate p0, Coordinate p1, Coordinate p2) {
  Coordinate pp0 = trans.transform(p0, new Coordinate.empty2D());
  Coordinate pp1 = trans.transform(p1, new Coordinate.empty2D());
  Coordinate pp2 = trans.transform(p2, new Coordinate.empty2D());

  AffineTransformationBuilder atb =
      new AffineTransformationBuilder(p0, p1, p2, pp0, pp1, pp2);
  AffineTransformation atbTrans = atb.getTransformation();

  Coordinate dest = new Coordinate.empty2D();
  assertEqualPoint(pp0, atbTrans.transform(p0, dest));
  assertEqualPoint(pp1, atbTrans.transform(p1, dest));
  assertEqualPoint(pp2, atbTrans.transform(p2, dest));
}

void assertEqualPoint(Coordinate p, Coordinate q) {
  assertEqualsD(p.x, q.x, 0.00005);
  assertEqualsD(p.y, q.y, 0.00005);
}
