import "package:test/test.dart";
import 'package:dart_sfs/dart_sfs.dart';

main() {
  group("orientation -", () {
    test("counter-clock-wise", () {
      var p = Coordinate(0, 0);
      var q = Coordinate(10, 0);
      var r = Coordinate(5, 2);
      var o = orientation(p, q, r);
      expect(o, 1);
    });
    test("clock-wise", () {
      var p = Coordinate(0, 0);
      var q = Coordinate(10, 0);
      var r = Coordinate(5, -2);
      var o = orientation(p, q, r);
      expect(o, -1);
    });

    test("collinear on half-line", () {
      var p = Coordinate(0, 0);
      var q = Coordinate(10, 0);
      var r = Coordinate(20, 0);
      var o = orientation(p, q, r);
      expect(o, 0);
    });

    test("collinear on segment", () {
      var p = Coordinate(0, 0);
      var q = Coordinate(10, 0);
      var r = Coordinate(5, 0);
      var o = orientation(p, q, r);
      expect(o, 0);
    });

    test("collinear on end point", () {
      var p = Coordinate(0, 0);
      var q = Coordinate(10, 0);
      var r = Coordinate(10, 0);
      var o = orientation(p, q, r);
      expect(o, 0);
    });

    test("slight distortion from collinearity", () {
      var p = Coordinate(0, 0);
      var q = Coordinate(10, 0);
      var r = Coordinate(20, 1E10 - 10);
      var o = orientation(p, q, r);
      expect(o, 1);
    });
  });

  group("segment relations -", () {
    test("two intersecting segments", () {
      var s1 = LineSegment.from([0, 0], [10, 0]);
      var s2 = LineSegment.from([5, 5], [5, -5]);
      expect(s1.intersects(s2), true);
    });

    test("two non-intersecting segments", () {
      var s1 = LineSegment.from([0, 0], [10, 0]);
      var s2 = LineSegment.from([11, 5], [11, -5]);
      expect(s1.intersects(s2), false);
    });

    test("two intersecting segments", () {
      var s1 = LineSegment.from([0, 0], [10, 0]);
      var s2 = LineSegment.from([9, 10], [11, -10]);
      expect(s1.intersects(s2), true);
    });

    test("two non-intersecting segments", () {
      var s1 = LineSegment.from([0, 0], [10, 0]);
      var s2 = LineSegment.from([9, 10], [12, -10]);
      expect(s1.intersects(s2), false);
    });

    test("two colinear segments, but no overlap", () {
      var s1 = LineSegment.from([0, 0], [1, 0]);
      var s2 = LineSegment.from([2, 0], [3, 0]);
      expect(s1.isColinear(s2), true);
      expect(s2.isColinear(s1), true);
      expect(s1.overlaps(s2), false);
      expect(s1.connectsTo(s2), false);
      expect(s1.intersects(s2), false);
    });

    test("two colinear segments, with overlap", () {
      var s1 = LineSegment.from([0, 0], [2, 0]);
      var s2 = LineSegment.from([1, 0], [4, 0]);
      expect(s1.isColinear(s2), true);
      expect(s2.isColinear(s1), true);
      expect(s1.overlaps(s2), true);
      expect(s1.connectsTo(s2), false);
      expect(s1.intersects(s2), false);
    });

    test("two colinear segments, not overlapping, but connecting", () {
      var s1 = LineSegment.from([0, 0], [2, 0]);
      var s2 = LineSegment.from([2, 0], [4, 0]);
      expect(s1.isColinear(s2), true);
      expect(s2.isColinear(s1), true);
      expect(s1.overlaps(s2), false);
      expect(s1.connectsTo(s2), true);
      expect(s1.intersects(s2), true);
    });

    test("two identical line segments", () {
      var s1 = LineSegment.from([0, 0], [2, 0]);
      var s2 = LineSegment.from([0, 0], [2, 0]);
      expect(s1.isColinear(s2), true);
      expect(s2.isColinear(s1), true);
      expect(s1.overlaps(s2), true);
      expect(s1.connectsTo(s2), false);
      expect(s1.intersects(s2), false);
    });
  });

  group("segment intersection -", () {
    test("two intersecting segments", () {
      var s1 = LineSegment.from([-2, 0], [2, 0]);
      var s2 = LineSegment.from([0, 2], [0, -2]);
      var intersection = s1.intersectionWith(s2);
      expect(intersection, equals(Coordinate(0, 0)));
    });

    test("two non-intersecting segments", () {
      var s1 = LineSegment.from([-2, 0], [2, 0]);
      var s2 = LineSegment.from([1, 1], [5, 5]);
      var intersection = s1.intersectionWith(s2);
      expect(intersection, isNull);
    });

    test("two colinear, overlapping segments", () {
      var s1 = LineSegment.from([-2, 0], [2, 0]);
      var s2 = LineSegment.from([0, 0], [4, 0]);
      var intersection = s1.intersectionWith(s2);
      expect(intersection, isNull);
    });

    test("two connected segments, colinaear, connected at end", () {
      var s1 = LineSegment.from([-2, 0], [2, 0]);
      var s2 = LineSegment.from([2, 0], [4, 0]);
      var intersection = s1.intersectionWith(s2);
      expect(intersection, Coordinate(2, 0));
    });

    test("two connected segments, non-colienar, connected at start", () {
      var s1 = LineSegment.from([-2, -2], [0, 0]);
      var s2 = LineSegment.from([0, 0], [4, 7]);
      var intersection = s1.intersectionWith(s2);
      expect(intersection, Coordinate(0, 0));
    });
  });

  group("_EventQueue -", () {
    group("constructor -", () {
      test("create an event queue", () {
        var queue = EventQueue();
        expect(queue.isEmpty, true);
      });
    });

    group("addEvent -", () {
      test("add a non-exitent event", () {
        var queue = EventQueue();
        queue.addEvent(Coordinate(0, 1));
        expect(queue.length, 1);
        expect(queue.first.pos.x, 0);
        expect(queue.first.pos.y, 1);
        expect(queue.first.segments.isEmpty, true);
      });

      test("add an already existing event", () {
        var queue = EventQueue();
        queue.addEvent(Coordinate(0, 1));
        queue.addEvent(Coordinate(0, 1));
        expect(queue.length, 1);
        expect(queue.first.pos.x, 0);
        expect(queue.first.pos.y, 1);
        expect(queue.first.segments.isEmpty, true);
      });
    });

    group("addLineSegmentEvents -", () {
      test("for a segment", () {
        var queue = EventQueue();
        var ls = LineSegment(Coordinate(0, 0), Coordinate(10, 10));
        queue.addLineSegmentEvents(ls);
        expect(queue.length, 2);
        expect(queue.first.pos.x, 10);
        expect(queue.first.pos.y, 10);
        expect(queue.first.segments.length, 1);
        expect(queue.first.segments.first, ls);
      });
    });

    group("unshift -", () {
      test("from a queue with one event", () {
        var queue = EventQueue();
        queue.addEvent(Coordinate(0, 1));

        var event = queue.unshift();
        expect(event.pos.x, 0);
        expect(event.segments.isEmpty, true);
        expect(queue.isEmpty, true);
      });

      test("from a queue with one event", () {
        var queue = EventQueue();
        var ls = LineSegment(Coordinate(0, 0), Coordinate(10, 10));
        queue.addLineSegmentEvents(ls);
        var event = queue.unshift();
        expect(event.pos.x, 10);
        expect(event.segments.length, 1);
        expect(queue.length, 1);
        expect(queue.first.pos.x, 0);
        expect(queue.first.segments.isEmpty, true);
      });
    });
  });

  group("flatten -", () {
    Iterable flatten(e) =>
        e is Iterable ? e.expand(flatten) : List.filled(1, e);

    test("flatten", () {
      expect(flatten([1, 2]).toList(), equals([1, 2]));
      expect(
          flatten([
            1,
            [2, 3],
            4
          ]).toList(),
          equals([1, 2, 3, 4]));
      expect(
          flatten([
            1,
            [
              2,
              [3, 4, 5],
              [6],
              []
            ],
            7
          ]).toList(),
          equals([1, 2, 3, 4, 5, 6, 7]));
    });
  });

  group("line intersection -", () {
    test("simple case - two lines intersecting in one point", () {
      var s1 = LineSegment.from([-5, 5], [5, -5]);
      var s2 = LineSegment.from([2, 2], [-2, -2]);
      var intersections = computeLineIntersections([s1, s2]);
      expect(intersections.length, 1);
      expect(intersections.first.pos, equals(Coordinate(0, 0)));
      expect(intersections.first.intersecting, unorderedEquals([s1, s2]));
    });

    test("simple case - two lines being connected in one point", () {
      var s1 = LineSegment.from([0, 0], [5, 0]);
      var s2 = LineSegment.from([5, 0], [7, 7]);
      var intersections = computeLineIntersections([s1, s2]);
      expect(intersections.length, 1);
      expect(intersections.first.pos, equals(Coordinate(5, 0)));
      expect(intersections.first.intersecting, unorderedEquals([s1, s2]));
    });

    test("two overlapping lines", () {
      var s1 = LineSegment.from([5, 5], [3, 3]);
      var s2 = LineSegment.from([4, 4], [2, 2]);
      var intersections = computeLineIntersections([s1, s2]);
      expect(intersections.length, 2);
      expect(intersections.elementAt(0).pos, equals(Coordinate(4, 4)));
      expect(
          intersections.elementAt(0).intersecting, unorderedEquals([s1, s2]));
      expect(intersections.elementAt(1).pos, equals(Coordinate(3, 3)));
      expect(
          intersections.elementAt(1).intersecting, unorderedEquals([s1, s2]));
    });
  });
}
