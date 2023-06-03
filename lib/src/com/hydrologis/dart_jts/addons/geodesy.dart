part of dart_jts;

/// The main geodesy class
///
/// This class is an adapted copy of the original https://github.com/wingkwong/geodesy
/// project that is no longer maintained. Changes have been done to work with
/// coordinates instead of latlng.
class Geodesy {
  final num _RADIUS = 6371e3; // meters
  final num _PI = math.pi;

  /// calculate a destination point given the distance and bearing
  Coordinate destinationPointByDistanceAndBearing(
      Coordinate l, num distance, num bearing,
      [num? radius]) {
    radius = radius ?? _RADIUS;

    num angularDistanceRadius = distance / radius;
    num bearingRadians = degToRadian(bearing as double);

    num latRadians = degToRadian(l.y);
    num lngRadians = degToRadian(l.x);

    num sinLatRadians = math.sin(latRadians);
    num cosLatRadians = math.cos(latRadians);
    num sinAngularDistanceRadius = math.sin(angularDistanceRadius);
    num cosAngularDistanceRadius = math.cos(angularDistanceRadius);
    num sinBearingRadians = math.sin(bearingRadians);
    num cosBearingRadians = math.cos(bearingRadians);

    var sinLatRadians2 = sinLatRadians * cosAngularDistanceRadius +
        cosLatRadians * sinAngularDistanceRadius * cosBearingRadians;
    num latRadians2 = math.asin(sinLatRadians2);
    var y = sinBearingRadians * sinAngularDistanceRadius * cosLatRadians;
    var x = cosAngularDistanceRadius - sinLatRadians * sinLatRadians2;
    num lngRadians2 = lngRadians + math.atan2(y, x);
    return Coordinate((radianToDeg(lngRadians2 as double) + 540) % 360 - 180,
        radianToDeg(latRadians2 as double));
  }

  /// calcuate the midpoint bewteen teo geo points
  Coordinate midPointBetweenTwoGeoPoints(Coordinate l1, Coordinate l2) {
    num l1LatRadians = degToRadian(l1.y);
    num l1LngRadians = degToRadian(l1.x);
    num l2LatRadians = degToRadian(l2.y);
    num lngRadiansDiff = degToRadian(l2.x - l1.x);

    num vectorX = math.cos(l2LatRadians) * math.cos(lngRadiansDiff);
    num vectorY = math.cos(l2LatRadians) * math.sin(lngRadiansDiff);

    num x = math.sqrt((math.cos(l1LatRadians) + vectorX) *
            (math.cos(l1LatRadians) + vectorX) +
        vectorY * vectorY);
    num y = math.sin(l1LatRadians) + math.sin(l2LatRadians);
    num latRadians = math.atan2(y, x);
    num lngRadians =
        l1LngRadians + math.atan2(vectorY, math.cos(l1LatRadians) + vectorX);

    return Coordinate((radianToDeg(lngRadians as double) + 540) % 360 - 180,
        radianToDeg(latRadians as double));
  }

  /// calculate the geo point of intersection of two given paths
  Coordinate? intersectionByPaths(
      Coordinate l1, Coordinate l2, num b1, num b2) {
    num l1LatRadians = degToRadian(l1.y);
    num l1LngRadians = degToRadian(l1.x);
    num l2LatRadians = degToRadian(l2.y);
    num l2LngRadians = degToRadian(l2.x);
    num b1Radians = degToRadian(b1 as double);
    num b2Radians = degToRadian(b2 as double);
    var latRadiansDiff = l2LatRadians - l1LatRadians;
    var lngRadiansDiff = l2LngRadians - l1LngRadians;

    num angularDistance = 2 *
        math.asin(math.sqrt(
            math.sin(latRadiansDiff / 2) * math.sin(latRadiansDiff / 2) +
                math.cos(l1LatRadians) *
                    math.cos(l2LatRadians) *
                    math.sin(lngRadiansDiff / 2) *
                    math.sin(lngRadiansDiff / 2)));
    if (angularDistance == 0) return null;

    num initBearingX = math.acos((math.sin(l2LatRadians) -
            math.sin(l1LatRadians) * math.cos(angularDistance)) /
        (math.sin(angularDistance) * math.cos(l1LatRadians)));
    if (initBearingX.isNaN) initBearingX = 0;
    num initBearingY = math.acos((math.sin(l1LatRadians) -
            math.sin(l2LatRadians) * math.cos(angularDistance)) /
        (math.sin(angularDistance) * math.cos(l2LatRadians)));

    var finalBearingX = math.sin(l2LngRadians - l1LngRadians) > 0
        ? initBearingX
        : 2 * _PI - initBearingX;
    var finalBearingY = math.sin(l2LngRadians - l1LngRadians) > 0
        ? 2 * _PI - initBearingY
        : initBearingY;

    var angle1 = b1Radians - finalBearingX;
    var angle2 = finalBearingY - b2Radians;

    if (math.sin(angle1) == 0 && math.sin(angle2) == 0) return null;
    if (math.sin(angle1) * math.sin(angle2) < 0) return null;

    num angle3 = math.acos(-math.cos(angle1) * math.cos(angle2) +
        math.sin(angle1) * math.sin(angle2) * math.cos(angularDistance));
    num dst13 = math.atan2(
        math.sin(angularDistance) * math.sin(angle1) * math.sin(angle2),
        math.cos(angle2) + math.cos(angle1) * math.cos(angle3));
    num lat3 = math.asin(math.sin(l1LatRadians) * math.cos(dst13) +
        math.cos(l1LatRadians) * math.sin(dst13) * math.cos(b1Radians));
    num lngRadiansDiff13 = math.atan2(
        math.sin(b1Radians) * math.sin(dst13) * math.cos(l1LatRadians),
        math.cos(dst13) - math.sin(l1LatRadians) * math.sin(lat3));
    var l3LngRadians = l1LngRadians + lngRadiansDiff13;

    return Coordinate((radianToDeg(l3LngRadians as double) + 540) % 360 - 180,
        radianToDeg(lat3 as double));
  }

  /// calculate the distance in meters between two geo points
  num distanceBetweenTwoGeoPoints(Coordinate l1, Coordinate l2, [num? radius]) {
    radius = radius ?? _RADIUS;
    var R = radius;
    num l1LatRadians = degToRadian(l1.y);
    num l1LngRadians = degToRadian(l1.x);
    num l2LatRadians = degToRadian(l2.y);
    num l2LngRadians = degToRadian(l2.x);
    var latRadiansDiff = l2LatRadians - l1LatRadians;
    var lngRadiansDiff = l2LngRadians - l1LngRadians;

    num a = math.sin(latRadiansDiff / 2) * math.sin(latRadiansDiff / 2) +
        math.cos(l1LatRadians) *
            math.cos(l2LatRadians) *
            math.sin(lngRadiansDiff / 2) *
            math.sin(lngRadiansDiff / 2);
    num c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a));
    var distance = R * c;

    return distance;
  }

  /// calculate the bearing from point l1 to point l2
  num bearingBetweenTwoGeoPoints(Coordinate l1, Coordinate l2) {
    num l1LatRadians = degToRadian(l1.y);
    num l2LatRadians = degToRadian(l2.y);
    num lngRadiansDiff = degToRadian(l2.x - l1.x);
    num y = math.sin(lngRadiansDiff) * math.cos(l2LatRadians);
    num x = math.cos(l1LatRadians) * math.sin(l2LatRadians) -
        math.sin(l1LatRadians) *
            math.cos(l2LatRadians) *
            math.cos(lngRadiansDiff);
    num radians = math.atan2(y, x);

    return (radianToDeg(radians as double) + 360) % 360;
  }

  /// calculate the final bearing from point l1 to point l2
  num finalBearingBetweenTwoGeoPoints(Coordinate l1, Coordinate l2) {
    return (bearingBetweenTwoGeoPoints(l2, l1) + 180) % 360;
  }

  /// calculate signed distance from a geo point
  /// to greate circle with start and end points
  num crossTrackDistanceTo(Coordinate l1, Coordinate start, Coordinate end,
      [num? radius]) {
    var R = radius ?? _RADIUS;

    num distStartL1 = distanceBetweenTwoGeoPoints(start, l1, R) / R;
    num radiansStartL1 =
        degToRadian(bearingBetweenTwoGeoPoints(start, l1) as double);
    num radiansEndL1 =
        degToRadian(bearingBetweenTwoGeoPoints(start, end) as double);

    num x = math
        .asin(math.sin(distStartL1) * math.sin(radiansStartL1 - radiansEndL1));

    return x * R;
  }

  /// check if a given geo point is in the bouding box
  bool isGeoPointInBoudingBox(
      Coordinate l, Coordinate topLeft, Coordinate bottomRight) {
    return topLeft.y <= l.y &&
            l.y <= bottomRight.y &&
            topLeft.x <= l.x &&
            l.x <= bottomRight.x
        ? true
        : false;
  }

  /// check if a given geo point is in the a polygon
  /// using even-odd rule algorithm
  bool isGeoPointInPolygon(Coordinate l, List<Coordinate> polygon) {
    var isInPolygon = false;

    for (var i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
      if ((((polygon[i].y <= l.y) && (l.y < polygon[j].y)) ||
              ((polygon[j].y <= l.y) && (l.y < polygon[i].y))) &&
          (l.x <
              (polygon[j].x - polygon[i].x) *
                      (l.y - polygon[i].y) /
                      (polygon[j].y - polygon[i].y) +
                  polygon[i].x)) isInPolygon = !isInPolygon;
    }
    return isInPolygon;
  }

  /// Get a list of [LatLng] points within a distance from
  /// a given point
  ///
  /// Distance is in meters
  List<Coordinate> pointsInRange(
      Coordinate point, List<Coordinate> pointsToCheck, num distance) {
    final geoFencedPoints = <Coordinate>[];
    for (final p in pointsToCheck) {
      final distanceFromCenter = distanceBetweenTwoGeoPoints(point, p);
      if (distanceFromCenter <= distance) {
        geoFencedPoints.add(p);
      }
    }
    return geoFencedPoints;
  }

  /// Converts degree to radian
  double degToRadian(final double deg) => deg * (_PI / 180.0);

  /// Radian to degree
  double radianToDeg(final double rad) => rad * (180.0 / _PI);
}
