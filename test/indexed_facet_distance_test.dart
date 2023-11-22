import 'package:test/test.dart';
import 'package:dart_jts/dart_jts.dart';

void main() {
  final double TOLERANCE = 1E-10;
  var casFactory = CoordinateArraySequenceFactory();
  var geomFactory = GeometryFactory.withCoordinateSequenceFactory(casFactory);
  var rdr = WKTReader();
  List<Coordinate>? nearestPoints(Geometry g1, Geometry g2) {
    return IndexedFacetDistance.nearestPoints(g1, g2);
  }

  double distance(Geometry? g1, Geometry? g2) {
    return IndexedFacetDistance.distance(g1!, g2!);
  }

  bool isWithinDistance(Geometry g1, Geometry g2, double distance) {
    return IndexedFacetDistance.isWithinDistance(g1, g2, distance);
  }

  void checkDistanceNearestPoints(
      String wkt0, String wkt1, double distance, Coordinate p0, Coordinate p1) {
    Geometry? g0 = rdr.read(wkt0);
    Geometry? g1 = rdr.read(wkt1);
    List<Coordinate>? pointList = nearestPoints(g0!, g1!);
    expect(distance, closeTo(pointList![0].distance(pointList[1]), TOLERANCE));
    expect(p0.x, closeTo(pointList[0].x, TOLERANCE));
    expect(p0.y, closeTo(pointList[0].y, TOLERANCE));
    expect(p1.x, closeTo(pointList[1].x, TOLERANCE));
    expect(p1.y, closeTo(pointList[1].y, TOLERANCE));
  }

  test('testDisjointCollinearSegments', ()  {
    Geometry? g1 = rdr.read("LINESTRING (0.0 0.0, 9.9 1.4)");
    Geometry? g2 = rdr.read("LINESTRING (11.88 1.68, 21.78 3.08)");
    double dist = distance(g1, g2);
    expect(1.9996999774966246, closeTo(dist, 0.0001));
    expect(!isWithinDistance(g1!, g2!, 1), true);
    expect(isWithinDistance(g1, g2, 3), true);
  });
  test('testPolygonsDisjoint', () {
    Geometry? g1 = rdr.read(
        "POLYGON ((40 320, 200 380, 320 80, 40 40, 40 320),  (180 280, 80 280, 100 100, 220 140, 180 280))");
    Geometry? g2 =
        rdr.read("POLYGON ((160 240, 120 240, 120 160, 160 140, 160 240))");
    expect(18.97366596, closeTo(distance(g1, g2), 1E-5));
    expect(!isWithinDistance(g1!, g2!, 0), true);
    expect(!isWithinDistance(g1, g2, 10), true);
    expect(isWithinDistance(g1, g2, 20), true);
  });
  test('testPolygonsOverlapping', ()  {
    Geometry? g1 = rdr.read(
        "POLYGON ((40 320, 200 380, 320 80, 40 40, 40 320),  (180 280, 80 280, 100 100, 220 140, 180 280))");
    Geometry? g3 =
        rdr.read("POLYGON ((160 240, 120 240, 120 160, 180 100, 160 240))");
    expect(0.0, closeTo(distance(g1, g3), 1E-9));
    expect(isWithinDistance(g1!, g3!, 0.0), true);
  });
  test('testLinesIdentical', ()  {
    LineString? l1 = rdr.read("LINESTRING(10 10, 20 20, 30 40)") as LineString;
    expect(0.0, closeTo(distance(l1, l1), 1E-5));
    expect(isWithinDistance(l1, l1, 0), true);
  });
  test('testEmpty', ()  {
    Geometry? g1 = rdr.read("POINT (0 0)");
    Geometry? g2 = rdr.read("POLYGON EMPTY");
    expect(0.0 == g1!.distance(g2!), true);
  });
  test('testClosestPoints1', () {
    checkDistanceNearestPoints(
        "POLYGON ((200 180, 60 140, 60 260, 200 180))",
        "POINT (140 280)",
        57.05597791103589,
        Coordinate(111.6923076923077, 230.46153846153845),
        Coordinate(140, 280));
  });
  test('testClosestPoints2', ()  {
    checkDistanceNearestPoints(
        "POLYGON ((200 180, 60 140, 60 260, 200 180))",
        "MULTIPOINT ((140 280), (140 320))",
        57.05597791103589,
        Coordinate(111.6923076923077, 230.46153846153845),
        Coordinate(140, 280));
  });
  test('testClosestPoints3', ()  {
    checkDistanceNearestPoints(
        "LINESTRING (100 100, 200 100, 200 200, 100 200, 100 100)",
        "POINT (10 10)",
        127.27922061357856,
        Coordinate(100, 100),
        Coordinate(10, 10));
  });
  test('testClosestPoints4', ()  {
    checkDistanceNearestPoints(
        "LINESTRING (100 100, 200 200)",
        "LINESTRING (100 200, 200 100)",
        0.0,
        Coordinate(150, 150),
        Coordinate(150, 150));
  });
  test('testClosestPoints5', ()  {
    checkDistanceNearestPoints(
        "LINESTRING (100 100, 200 200)",
        "LINESTRING (150 121, 200 0)",
        20.506096654409877,
        Coordinate(135.5, 135.5),
        Coordinate(150, 121));
  });
  test('testClosestPoints6', ()  {
    checkDistanceNearestPoints(
        "POLYGON ((76 185, 125 283, 331 276, 324 122, 177 70, 184 155, 69 123, 76 185), (267 237, 148 248, 135 185, 223 189, 251 151, 286 183, 267 237))",
        "LINESTRING (153 204, 185 224, 209 207, 238 222, 254 186)",
        13.788860460124573,
        Coordinate(139.4956500724988, 206.78661188980183),
        Coordinate(153, 204));
  });
  test('testClosestPoints7', ()  {
    checkDistanceNearestPoints(
        "POLYGON ((76 185, 125 283, 331 276, 324 122, 177 70, 184 155, 69 123, 76 185), (267 237, 148 248, 135 185, 223 189, 251 151, 286 183, 267 237))",
        "LINESTRING (120 215, 185 224, 209 207, 238 222, 254 186)",
        0.0,
        Coordinate(120, 215),
        Coordinate(120, 215));
  });
}
