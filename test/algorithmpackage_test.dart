import 'dart:io';

import "package:test/test.dart";
import 'package:dart_jts/dart_jts.dart';
import "dart:math" as math;
import 'testing_utilities.dart';

const double TOLERANCE = 1e-10;

void main() {
  group("Centroid Tests - ", () {
    test("testCentroidMultiPolygon", () {
      // Verify that the computed centroid of a MultiPolygon is equivalent to the
      // area-weighted average of its components.
      Geometry g = WKTReader().read(
          "MULTIPOLYGON ((( -92.661322 36.58994900000003, -92.66132199999993 36.58994900000005, -92.66132199999993 36.589949000000004, -92.661322 36.589949, -92.661322 36.58994900000003)), (( -92.65560500000008 36.58708800000005, -92.65560499999992 36.58708800000005, -92.65560499998745 36.587087999992576, -92.655605 36.587088, -92.65560500000008 36.58708800000005 )), (( -92.65512450000065 36.586800000000466, -92.65512449999994 36.58680000000004, -92.65512449998666 36.5867999999905, -92.65512450000065 36.586800000000466 )))")!;

      assertTrue(areaWeightedCentroid(g)
          .equals2DWithTolerance(g.getCentroid().getCoordinate()!, TOLERANCE));
    });
  });

  group("Convex Hull Tests - ", () {
    PrecisionModel precisionModel = PrecisionModel.fixedPrecision(1000);
    GeometryFactory geometryFactory =
        GeometryFactory.withPrecisionModelSrid(precisionModel, 0);

    test("testManyIdenticalPoints", () {
      List<Coordinate> pts = [];
      for (int i = 0; i < 99; i++) pts.add(Coordinate(0, 0));
      pts.add(Coordinate(1, 1));
      ConvexHull ch = ConvexHull.fromPoints(pts, geometryFactory);
      Geometry actualGeometry = ch.getConvexHull();
      Geometry? expectedGeometry = WKTReader().read("LINESTRING (0 0, 1 1)");
      assertTrue(expectedGeometry?.equalsExactGeom(actualGeometry));
    });
    test("testAllIdenticalPoints", () {
      List<Coordinate> pts = [];
      for (int i = 0; i < 100; i++) pts.add(Coordinate(0, 0));
      ConvexHull ch = ConvexHull.fromPoints(pts, geometryFactory);
      Geometry actualGeometry = ch.getConvexHull();
      Geometry? expectedGeometry = WKTReader().read("POINT (0 0)");
      assertTrue(expectedGeometry?.equalsExactGeom(actualGeometry));
    });
    test("test1", () {
      WKTReader reader = WKTReader.withFactory(
          GeometryFactory.withPrecisionModelSrid(
              PrecisionModel.fixedPrecision(1), 0));
      LineString lineString =
          reader.read("LINESTRING (30 220, 240 220, 240 220)") as LineString;
      LineString convexHull =
          reader.read("LINESTRING (30 220, 240 220)") as LineString;
      assertTrue(convexHull.equalsExactGeom(lineString.convexHull()));
    });
    test("test2", () {
      WKTReader reader = WKTReader.withFactory(
          GeometryFactory.withPrecisionModelSrid(
              PrecisionModel.fixedPrecision(1), 0));
      Geometry? geometry = reader.read(
          "MULTIPOINT (130 240, 130 240, 130 240, 570 240, 570 240, 570 240, 650 240)");
      LineString convexHull =
          reader.read("LINESTRING (130 240, 650 240)") as LineString;
      assertTrue(convexHull.equalsExactGeom(geometry!.convexHull()));
    });
    test("test3", () {
      WKTReader reader = WKTReader.withFactory(
          GeometryFactory.withPrecisionModelSrid(
              PrecisionModel.fixedPrecision(1), 0));
      Geometry? geometry = reader.read("MULTIPOINT (0 0, 0 0, 10 0)");
      LineString convexHull =
          reader.read("LINESTRING (0 0, 10 0)") as LineString;
      assertTrue(convexHull.equalsExactGeom(geometry!.convexHull()));
    });
    test("test4", () {
      WKTReader reader = WKTReader.withFactory(
          GeometryFactory.withPrecisionModelSrid(
              PrecisionModel.fixedPrecision(1), 0));
      Geometry? geometry = reader.read("MULTIPOINT (0 0, 10 0, 10 0)");
      LineString convexHull =
          reader.read("LINESTRING (0 0, 10 0)") as LineString;
      assertTrue(convexHull.equalsExactGeom(geometry!.convexHull()));
    });
    test("test5", () {
      WKTReader reader = WKTReader.withFactory(
          GeometryFactory.withPrecisionModelSrid(
              PrecisionModel.fixedPrecision(1), 0));
      Geometry? geometry = reader.read("MULTIPOINT (0 0, 5 0, 10 0)");
      LineString convexHull =
          reader.read("LINESTRING (0 0, 10 0)") as LineString;
      assertTrue(convexHull.equalsExactGeom(geometry!.convexHull()));
    });
    test("test6", () {
      WKTReader reader = WKTReader.withFactory(
          GeometryFactory.withPrecisionModelSrid(
              PrecisionModel.fixedPrecision(1), 0));
      Geometry? geometry = reader.read("MULTIPOINT (0 0, 5 1, 10 0)");
      Geometry? convexHull = reader.read("POLYGON ((0 0, 5 1, 10 0, 0 0))");
      assertTrue(convexHull!.equalsExactGeom(geometry!.convexHull()));
    });
    test("test7", () {
      WKTReader reader = WKTReader.withFactory(
          GeometryFactory.withPrecisionModelSrid(
              PrecisionModel.fixedPrecision(1), 0));
      Geometry? geometry =
          reader.read("MULTIPOINT (0 0, 0 0, 5 0, 5 0, 10 0, 10 0)");
      LineString convexHull =
          reader.read("LINESTRING (0 0, 10 0)") as LineString;
      assertTrue(convexHull.equalsExactGeom(geometry!.convexHull()));
    });

    // test convexhull from multipoint found in wkt file data/convexhull.wkt
    test("testConvexHullFromMultipoint", () {
      // read file data/convexhull.wkt
      String filename = "./test//data/convexhull.wkt";
      var file = new File(filename);
      var wkt = file.readAsStringSync();
      WKTReader reader = WKTReader();
      Geometry? geometry = reader.read(wkt);
      Geometry? convexHull = geometry!.convexHull();

      var expectedWkt =
          "POLYGON ((12.19209 44.41011016666666, 11.429490500000002 44.52205366666667, 11.429474833333332 44.52207583333334, 12.212261666666667 44.43709083333333, 12.212392833333334 44.43656716666666, 12.212399333333334 44.43644266666667, 12.202527513493703 44.41216171584291, 12.192167333333336 44.41011516666667, 12.19209 44.41011016666666))";
      Geometry? expectedGeometry = reader.read(expectedWkt);
      assertTrue(expectedGeometry!.equalsExactGeom(convexHull));
    });
  });
}

/** Compute the centroid of a geometry as an area-weighted average of the centroids
 * of its components.
 *
 * @param g a polygonal geometry
 * @return Coordinate of the geometry's centroid
 */
Coordinate areaWeightedCentroid(Geometry g) {
  double totalArea = g.getArea();
  double cx = 0;
  double cy = 0;

  for (int i = 0; i < g.getNumGeometries(); i++) {
    Geometry component = g.getGeometryN(i);
    double areaFraction = component.getArea() / totalArea;

    Coordinate componentCentroid = component.getCentroid().getCoordinate()!;

    cx += areaFraction * componentCentroid.x;
    cy += areaFraction * componentCentroid.y;
  }

  return new Coordinate(cx, cy);
}
