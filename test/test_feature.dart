import "package:test/test.dart";
import 'package:dart_sfs/dart_sfs.dart';

main() {
  group("Feature -", () {
    group("constructors", () {
      test("standard - accepts a geomtry and properties", () {
        var geometry = Point(1, 2);
        var properties = {"k1": "v1", "k2": 122};
        var feature = Feature(geometry, properties);
        expect(feature.geometry.getCoordinate().x, 1);
        expect(feature.properties["k1"], "v1");
        expect(feature.properties.length, 2);
      });

      test("accepts a geometry only", () {
        var geometry = Point(1, 2);
        var feature = Feature(geometry);
        expect(feature.geometry.getCoordinate().x, 1);
        expect(feature.properties.isEmpty, true);
      });

      test("accepts a geomtry and null properties", () {
        var geometry = Point(1, 2);
        var feature = Feature(geometry, null);
        expect(feature.geometry.getCoordinate().x, 1);
        expect(feature.properties.isEmpty, true);
      });

      test("accepts a geomtry and empty properties", () {
        var geometry = Point(1, 2);
        var feature = Feature(geometry, {});
        expect(feature.geometry.getCoordinate().x, 1);
        expect(feature.properties.isEmpty, true);
      });
    });

    group("json parsing", () {
      test("simple point featue", () {
        var json = """
        {"type": "Feature",
            "geometry":  {"type": "Point", "coordinates": [1,2]},
            "properties": {
               "k1": "string value",
               "k2": 123
            }
        }
        """;
        var feature = parseGeoJson(json);
        expect(feature.geometry.y, 2);
        expect(feature.properties.length, 2);
        expect(feature.properties["k2"], 123);
      });
    });
  });

  group("FeatureCollection -", () {
    group("constructors", () {
      test("standard without parameters", () {
        var fc = FeatureCollection();
        expect(fc.isEmpty, true);
      });
      test("standard with null features", () {
        var fc = FeatureCollection(null);
        expect(fc.isEmpty, true);
      });
      test("standard with empty features", () {
        var fc = FeatureCollection([]);
        expect(fc.isEmpty, true);
      });

      test("with two features", () {
        var f1 = Feature(Point(1, 2));
        var f2 = Feature(Point(3, 4), {"k1": "v1"});
        var fc = FeatureCollection([f1, f2]);
        expect(fc.isEmpty, false);
        expect(fc.length, 2);
        expect(fc.first.geometry.getCoordinate().x, 1);
        expect(fc.last.geometry.getCoordinate().y, 4);
      });

      test("accepts a geometry only", () {
        var geometry = Point(1, 2);
        var feature = Feature(geometry);
        expect(feature.geometry.getCoordinate().x, 1);
        expect(feature.properties.isEmpty, true);
      });

      test("accepts a geomtry and null properties", () {
        var geometry = Point(1, 2);
        var feature = Feature(geometry, null);
        expect(feature.geometry.getCoordinate().x, 1);
        expect(feature.properties.isEmpty, true);
      });

      test("accepts a geomtry and empty properties", () {
        var geometry = Point(1, 2);
        var feature = Feature(geometry, {});
        expect(feature.geometry.getCoordinate().x, 1);
        expect(feature.properties.isEmpty, true);
      });
    });

    group("json parsing", () {
      test("a feature collection with two features", () {
        var json = """
            {"type": "FeatureCollection",
            "features": [
            {"type": "Feature",
            "geometry":  {"type": "Point", "coordinates": [1,2]},
            "properties": {
            "k1": "string value",
            "k2": 123
            }
            },
            {"type": "Feature",
            "geometry":  {"type": "MultiPoint", "coordinates": [[1,2], [3,4], [5,6]]}
            }
            ]
            }
            """;
        var fc = parseGeoJson(json);
        expect(fc is FeatureCollection, true);
        expect(fc.length, 2);
      });
    });
  });
}
