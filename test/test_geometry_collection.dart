import "package:test/test.dart";
import 'package:dart_sfs/dart_sfs.dart';

main() {
  group("constructors -", () {
    test("create an empty collection", () {
      var g = GeometryCollection.empty();
      expect(g.isEmpty, true);
    });

    test("create a geometry collectioin with a point and a line string", () {
      var p = Point(1, 2);
      var ls = LineString([Point(1, 2), Point(3, 4)]);
      var g = GeometryCollection([p, ls]);
      expect(g.isEmpty, false);
      expect(g.length, 2);
    });

    test("from WKT - empty collection", () {
      var g = GeometryCollection.wkt("geometrycollection empty");
      expect(g.isEmpty, true);
    });

    test("from WKT - non empty collection", () {
      var wkt = """geometrycollection (
         point (1 2),
         linestring (1 2, 3 4, 5 6)
      )
      """;
      var g = GeometryCollection.wkt(wkt);
      expect(g.isEmpty, false);
      expect(g.length, 2);
      expect(g.first.runtimeType, Point);
      expect(g[1].runtimeType, LineString);
    });
  });
  /* ----------------------------------------------------------------- */
  group("geojson -", () {
    test("- deserialize a GeometryCollection", () {
      var gjson = """
      {"type": "GeometryCollection", "geometries": [
        {"type":"Point", "coordinates":[1,2]},
        {"type":"MultiPoint", "coordinates": [[1,2],[3,4]]},
        {"type":"LineString", "coordinates": [[1,2],[3,4],[5,6],[7,8]]}
        ]
      }
      """;
      var o = parseGeoJson(gjson);
      expect(o is GeometryCollection, true);
      o = (o as GeometryCollection);
      expect(o.length, 3);
      expect(o[0] is Point, true);
      expect(o[1] is MultiPoint, true);
      expect(o[2] is LineString, true);
    });
  });
}
