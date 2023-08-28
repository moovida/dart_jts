import 'package:dart_jts/dart_jts.dart';
import 'package:test/test.dart';

/// This class is an adapted copy of the original https://github.com/wingkwong/geodesy
/// project that is no longer maintained. Changes have been done to work with
/// coordinates instead of latlng.
void main() {
  var geodesy = Geodesy();

  test('distanceBetweenTwoGeoPoints', () async {
    final l1 = Coordinate.fromYX(50.06638889, 5.71472222);
    final l2 = Coordinate.fromYX(58.64388889, 3.07000000);
    final distance = geodesy.distanceBetweenTwoGeoPoints(l1, l2);
    expect(distance, 969938.8877266833);
  });

  test('destinationPointByDistanceAndBearing', () async {
    final l3 = Coordinate.fromYX(51.4778, -0.0015);
    final destinationPoint =
        geodesy.destinationPointByDistanceAndBearing(l3, 7794.0, 300.7);
    expect(destinationPoint,
        Coordinate.fromYX(51.51350573766377, -0.09823692862482858));
  });

  test('bearingBetweenTwoGeoPoints', () async {
    final l4 = Coordinate.fromYX(52.205, 0.119);
    final l5 = Coordinate.fromYX(48.857, 2.351);
    final bearing = geodesy.bearingBetweenTwoGeoPoints(l4, l5);
    expect(bearing.toStringAsFixed(6), 156.16658258153166.toStringAsFixed(6));
  });

  test('bearingBetweenTwoGeoPoints', () async {
    final l4 = Coordinate.fromYX(52.205, 0.119);
    final l5 = Coordinate.fromYX(48.857, 2.351);
    final finalBearing = geodesy.finalBearingBetweenTwoGeoPoints(l4, l5);
    expect(finalBearing, 157.89044019049243);
  });

  test('midPointBetweenTwoGeoPoints', () async {
    final l4 = Coordinate.fromYX(52.205, 0.119);
    final l5 = Coordinate.fromYX(48.857, 2.351);
    final midpoint = geodesy.midPointBetweenTwoGeoPoints(l4, l5);
    expect(
        midpoint.equals2DWithTolerance(
            Coordinate.fromYX(50.53632687827432, 1.2746141006782636),
            0.0000001),
        true);
  });

  test('midPointBetweenTwoGeoPoints2', () async {
    final l4 = Coordinate.fromYX(52.205, 0.119);
    final l5 = Coordinate.fromYX(48.857, 2.351);
    final midpoint = geodesy.midPointBetweenTwoGeoPoints(l4, l5);
    expect(
        midpoint.equals2DWithTolerance(
            Coordinate.fromYX(50.53632687827432, 1.2746141006782636), 0.000001),
        true);
  });

  test('isGeoPointInBoudingBox', () async {
    final l3 = Coordinate.fromYX(51.4778, -0.0015);
    final l4 = Coordinate.fromYX(52.205, 0.119);
    final l5 = Coordinate.fromYX(48.857, 2.351);
    final inBoudingBox = geodesy.isGeoPointInBoudingBox(l3, l5, l4);
    expect(inBoudingBox, false);
  });

  test('intersectionByPaths', () async {
    final l4 = Coordinate.fromYX(52.205, 0.119);
    final l5 = Coordinate.fromYX(48.857, 2.351);
    final b1 = 108.547;
    final b2 = 32.435;
    final intersectionByPaths = geodesy.intersectionByPaths(l4, l5, b1, b2);
    expect(
        intersectionByPaths!.equals2DWithTolerance(
            Coordinate.fromYX(51.15151144654275, 4.698604211862175), 0.0000001),
        true);
  });

  test('crossTrackDistanceTo', () async {
    final l4 = Coordinate.fromYX(52.205, 0.119);
    final l5 = Coordinate.fromYX(48.857, 2.351);
    final l6 = Coordinate.fromYX(50.587, 1.231);
    final distanceToGreatCircle = geodesy.crossTrackDistanceTo(l4, l5, l6);
    expect(distanceToGreatCircle, 1241.7274005073216);
  });

  test('crossTrackDistanceTo', () async {
    var poly = <Coordinate>[
      Coordinate.fromYX(1.0, 1.0),
      Coordinate.fromYX(1.0, 2.0),
      Coordinate.fromYX(2.0, 2.0),
      Coordinate.fromYX(2.0, 1.0)
    ];
    var l7 = Coordinate.fromYX(1.5, 1.5);
    var isGeoPointInPolygon = geodesy.isGeoPointInPolygon(l7, poly);
    expect(isGeoPointInPolygon, true);
  });

  test('pointsInRange', () async {
    final point = Coordinate.fromYX(51.0, 0);
    final distance = 10000;
    final pointNotInRange = geodesy.destinationPointByDistanceAndBearing(
        point, distance + 10, 420.0);
    final pointInRange = geodesy.destinationPointByDistanceAndBearing(
        point, distance - 10, 420.0);
    final pointsToCheck = <Coordinate>[pointInRange, pointNotInRange];
    final geofencedPoints =
        geodesy.pointsInRange(point, pointsToCheck, distance);
    expect((geofencedPoints.contains(pointInRange)), true);
    expect((geofencedPoints.contains(pointNotInRange)), false);
  });

  test('area and length', () async {
    const world_area = 511207893395811.06;
    const world_perim = 40075016.69;
    const wkt = "POLYGON((-180 -90,-180 90,180 90,180 -90,-180 -90))";

    Geometry polygon = WKTReader().read(wkt)!;
    var area = geodesy.area(polygon);
    expect(area, world_area);

    var length = geodesy.length(polygon);
    expect(length, world_perim);
  });
}
