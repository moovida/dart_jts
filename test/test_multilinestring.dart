import 'package:dart_sfs/dart_sfs.dart';
import "package:test/test.dart";

main() {
//  group("constructor -", () {
//    test("empty() constructor", () {
//      var mls = MultiLineString.empty();
//      expect(mls.isEmpty, true);
//    });
//    test("null list leads to empty multilinestring", () {
//      var mls = MultiLineString(null);
//      expect(mls.isEmpty, true);
//    });
//    test("empty list of linestrings leads to empty multilinestring", () {
//      var mls = MultiLineString([]);
//      expect(mls.isEmpty, true);
//    });
//
//    test("create a multilinestring with three children", () {
//      List<LineString> children = [
//        parseWKT("linestring empty"),
//        parseWKT("linestring (1 2, 3 4, 5 6, 7 8)"),
//        parseWKT("linestring z (10 11 12, 20 21 22)"),
//      ];
//      var mls = MultiLineString(children);
//      expect(mls.isEmpty, false);
//      expect(mls.length, 3);
//    });
//    test("don't allow nulls in child list", () {
//      List<LineString> children = [
//        parseWKT("linestring empty"),
//        parseWKT("linestring (1 2, 3 4, 5 6, 7 8)"),
//        parseWKT("linestring z (10 11 12, 20 21 22)"),
//        null
//      ];
//      var mls;
//      expect(() => mls = MultiLineString(children),
//          throwsA(isInstanceOf<ArgumentError>()));
//    });
//  });
//
//  /* ---------------------------------------------------------------- */
//  group("boundary -", () {
//    test("of an empty multilinestring is empty", () {
//      var mls = MultiLineString.empty();
//      expect(mls.boundary.isEmpty, true);
//    });
//    test("one non closed child -> the boundary consists of the endpoints", () {
//      List<LineString> children = [parseWKT("linestring (1 2, 3 4, 5 6, 7 8)")];
//      var mls = MultiLineString(children);
//      expect(mls.boundary.getCoordinates().length, 2);
//      expect(mls.boundary.getCoordinates().where((p) => p.x == 1).length, 1);
//      expect(mls.boundary.getCoordinates().where((p) => p.x == 7).length, 1);
//    });
//    test("one closed child -> the boundary is empty", () {
//      List<LineString> children = [
//        parseWKT("linestring (1 2, 3 4, 5 6, 7 8, 1 2)")
//      ];
//      var mls = MultiLineString(children);
//      expect(mls.boundary.isEmpty, true);
//    });
//
//    test(
//        "two open children, one shared endpoint -> 2 of 4 endpoints in the boundary",
//        () {
//      List<LineString> children = [
//        parseWKT("linestring (1 2, 3 4, 5 6)"),
//        parseWKT("linestring (5 6, 7 8, 9 10)")
//      ];
//      var mls = MultiLineString(children);
//      expect(mls.boundary.getCoordinates().length, 2);
//      expect(mls.boundary.getCoordinates().where((p) => p.x == 1).length, 1);
//      expect(mls.boundary.getCoordinates().where((p) => p.x == 9).length, 1);
//    });
//  });
//
//  group("asText -", () {
//    test("of an empty multilinestring", () {
//      var mls = MultiLineString.empty();
//      mls = parseWKT(mls.asText);
//      expect(mls is MultiLineString, true);
//      expect(mls.isEmpty, true);
//    });
//
//    test("multi line string with three line strings", () {
//      List<LineString> children = [
//        parseWKT("linestring empty"),
//        parseWKT("linestring (1 2, 3 4, 5 6, 7 8)"),
//        parseWKT("linestring (10 11, 20 21)")
//      ];
//      var mls = MultiLineString(children);
//      mls = parseWKT(mls.asText);
//      expect(mls is MultiLineString, true);
//      expect(mls.length, 3);
//      expect(mls.first.isEmpty, true);
//    });
//  });
//
//  /* --------------------------------------------------------------------- */
//  group("geojson -", () {
//    test("- deserialize a multilinestring", () {
//      var gjson = """
//      {"type": "MultiLineString", "coordinates": [
//        [[1,2], [3,4], [5,6]],
//        [[11,12], [13,14], [15,16]],
//        [[21,22], [22,24], [25,26]]
//      ]}
//      """;
//      var o = parseGeoJson(gjson);
//      expect(o is MultiLineString, true);
//      expect(o.length, 3);
//      for (int i = 0; i < o.length; i++) {
//        expect(o[i] is LineString, true);
//      }
//      expectPoints(ls, points) {
//        for (int i = 0; i < points.length; i++) {
//          expect(ls[i].x, points[i][0]);
//          expect(ls[i].y, points[i][1]);
//        }
//      }
//
//      expectPoints(o[0], [
//        [1, 2],
//        [3, 4],
//        [5, 6]
//      ]);
//      expectPoints(o[1], [
//        [11, 12],
//        [13, 14],
//        [15, 16]
//      ]);
//      expectPoints(o[2], [
//        [21, 22],
//        [22, 24],
//        [25, 26]
//      ]);
//    });
//  });
}
