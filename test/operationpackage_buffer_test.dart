import "package:test/test.dart";
import 'package:dart_jts/dart_jts.dart';
import "dart:math" as math;
import 'testing_utilities.dart';

double TOLERANCE = 1E-5;

PrecisionModel precisionModel = PrecisionModel.fixedPrecision(1);

GeometryFactory geometryFactory = GeometryFactory.withPrecisionModelSrid(precisionModel, 0);

WKTReader reader = WKTReader.withFactory(geometryFactory);

void main() {
  group("LineSegmentTest - ", () {
    test("testProjectionFactor", () {
      // zero-length line
      LineSegment seg = new LineSegment(10, 0, 10, 0);
      assertTrue(seg.projectionFactor(new Coordinate(11, 0)).isNaN);

      LineSegment seg2 = new LineSegment(10, 0, 20, 0);
      assertTrue(seg2.projectionFactor(new Coordinate(11, 0)) == 0.1);
    });
  });
}
