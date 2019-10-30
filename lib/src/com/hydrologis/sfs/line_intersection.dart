part of dart_sfs;

// tolerance interval [-_EPSILON .. _EPSILON]  ~ 0
const _EPSILON = 10E-10;

/// Computes the counterclockwise orientation of point [r] with respect
/// to the line given by the start and end points [p] and [q].
///
/// Returns -1 if, if [r] is left to [p]-[q], 0 if is colinear with [p]-[q]
/// and 1 if it is to the right of [p]-[q].
int _orientation(Coordinate p, Coordinate q, Coordinate r) {
  var v = (p.x - r.x) * (q.y - r.y) - (q.x - r.x) * (p.y - r.y);
  if (v < -_EPSILON) return -1;
  if (v > _EPSILON) return 1;
  return 0;
}

/// determinante
num _det(p, q) => p.x * q.y - p.y * q.x;

/// Represents closed line segment given by two 2D points.
///
/// [start] is always the end point with the larger y-coordinate, or,
/// if y is equal, the lower x-coordinate. This is the order used
/// on events in the line intersection algorithm.
///
/// LineSegments are never collapsed, they always have exactly two distinct
/// end points.
class LineSegment {
  Coordinate _start;
  Coordinate _end;

  /// Creates a line segment with start position [start] and end position
  /// [end]. [start] and [end] don't have to be ordered, but the position
  /// with the higher y-coordinate, or, if y is equal, the lower x-coordinate
  /// becomes the [start] event. [start] and [end] must not be equal.
  ///
  LineSegment(Coordinate start, Coordinate end) {
    _require(start != end, "line segments must not be collapsed: $start, $end");
    int c = _comparePositionsInEventOrder(start, end);
    if (c <= 0) {
      _start = start;
      _end = end;
    } else {
      _start = end;
      _end = start;
    }
  }

  /// Creates a line segment from a [start] and [end] event, where
  /// [start] and [end] are list of [num]s with two elements.
  ///
  /// # Example
  ///     var s = new LineSegment.from([0,0], [2,2]);
  factory LineSegment.from(start, end) {
    _require(start is List);
    _require(end is List);
    _require(start.length == 2);
    _require(end.length == 2);
    return LineSegment(
        Coordinate(start[0], start[1]), Coordinate(end[0], end[1]));
  }

  /// returns true, if [p] is either the start or the end point
  /// of this segment
  bool isEndPoint(Coordinate p) => p == _start || p == _end;

  bool _fastExcludeIntersect(LineSegment other) {
    var minx = _start.x;
    var miny = _start.y < _end.y ? _start.y : _end.y;
    var maxx = _end.x;
    var maxy = _start.y >= _end.y ? _start.y : _end.y;
    if (other._start.x < minx && other._end.x < minx) return true;
    if (other._start.x > maxx && other._end.x > maxx) return true;
    if (other._start.y < miny && other._end.y < miny) return true;
    if (other._start.y > maxy && other._end.y > maxy) return true;
    return false;
  }

  /// This [LineSegment] intersects with an [other] iff the intersection
  /// consists of exactly one point.
  bool intersects(LineSegment other) {
    if (this == other) return false;
    if (_fastExcludeIntersect(other)) return false;
    if (this.connectsTo(other) && !this.overlaps(other)) return true;
    var o1 = _orientation(_start, _end, other._start);
    var o2 = _orientation(_start, _end, other._end);
    if (o1 == o2) return false;
    o1 = _orientation(other._start, other._end, _start);
    o2 = _orientation(other._start, other._end, _end);
    return o1 != o2;
  }

  /// This [LineSegment] connects to [other] iff the intersection consists
  /// of exactly one end point.
  bool connectsTo(LineSegment other) {
    return this != other &&
        (other.isEndPoint(_start) || other.isEndPoint(_end));
  }

  /// This [LineSegment] is colinear with [other] iff this and [other]
  /// lie on a straight line.
  ///
  bool isColinear(LineSegment other) {
    var o1 = _orientation(_start, _end, other._start);
    var o2 = _orientation(_start, _end, other._end);
    return o1 == 0 && o2 == 0;
  }

  /// This [LineSegment] and [other] overlap iff their intersection
  /// consists of more than one point.
  bool overlaps(LineSegment other) {
    if (this == other) return true;
    if (!isColinear(other)) return false;
    if (other.start < start && start < other.end) return true;
    if (other.start < end && end < other.end) return true;
    return false;
  }

  bool operator ==(other) => start == other.start && end == other.end;

  int get hashCode {
    var result = 17;
    result = 37 * result + start.hashCode;
    result = 37 * result + end.hashCode;
    return result;
  }

  bool get isHorizontal => _start.y == _end.y;

  /// Returns the intersection point of this segment with the segment
  /// [other] or null, if the segments don't intersection or if they
  /// overlap
  Coordinate intersectionWith(LineSegment other) {
    // formula and nomenclatura from
    // from http://bloggingmath.wordpress.com/2009/05/29/line-segment-intersection/
    var p = start;
    var r = end - start;
    var q = other.start;
    var s = other.end - other.start;
    var rs = _det(r, s);
    if (rs == 0) {
      if (isEndPoint(other.start)) return other.start;
      if (isEndPoint(other.end)) return other.end;
      // no proper intersection, because this and other are
      // colinear
      return null;
    }
    // fast check whether this segments are connected
    // (a common szenario in our data)
    if (isEndPoint(other.start)) return other.start;
    if (isEndPoint(other.end)) return other.end;

    // properly compute the intersection point, if any.
    var t = _det(q - p, s) / rs;
    var u = _det(q - p, r) / rs;
    if (t < 0 || t > 1 || u < 0 || u > 1) {
      // no intersection
      return null;
    }
    var intersection = p + r.scale(t);
    return intersection;
  }

  double get slope =>
      (_end.y - _start.y).toDouble() / (_end.x - _start.x).toDouble();

  /// the start point of this segment
  Coordinate get start => _start;

  ///  the end point of this segment
  Coordinate get end => _end;

  /// caches the orientation value
  num _ccwOrientation;

  /// Returns a number in the range [-1,1]. A horizontal line
  /// has the orientation 1, a vertical line has the orientation 0.
  ///
  /// If -1 < value < 0, then the segment has a positive slope, the smaller
  /// the value, the flatter the line segment.
  ///
  /// If 0 < value < 1, then the segment has a negative slope, the
  /// higher the value, the latter the line segment.
  num get counterclockwiseOrientation {
    if (_ccwOrientation != null) return _ccwOrientation;
    if (isHorizontal) return 1;
    var dx = end.x - start.x;
    var dy = end.y - start.y;
    var c = math.sqrt(dx * dx + dy * dy);
    _ccwOrientation = dx / c;
    return _ccwOrientation;
  }

  String toString() => "{LineSegment: start=$start, end=$end}";
}

/// Represent an event to be processed in the Bentley-Ottman-Algorithm
class _Event implements Comparable<_Event> {
  final Coordinate pos;
  List<LineSegment> _segments;

  _Event(this.pos, [this._segments = null]);

  int compareTo(_Event other) => _comparePositionsInEventOrder(pos, other.pos);

  add(LineSegment segment) {
    if (_segments == null) _segments = List<LineSegment>();
    if (!_segments.contains(segment)) _segments.add(segment);
  }

  List<LineSegment> get segments => _segments == null ? [] : _segments;

  bool operator ==(other) => compareTo(other) == 0;

  int get hashCode => pos.hashCode;

  toString() => "{event: pos=$pos, segments=$_segments}";
}

_comparePositionsInEventOrder(p1, p2) {
  _require(p1 != null);
  _require(p2 != null);
  int c = p1.y.compareTo(p2.y);
  // The "higher" the event, the "earlier" it is processed
  // p1 < p2 if p1.y > p2.y
  // p1 > p2 if p1.y < p2.y
  if (c != 0) return -c;

  // The "more left" the event, the "earlier" it is processed
  return p1.x.compareTo(p2.x);
}

class _EventQueue {
  /// both key and value are events. Having events as keys
  /// ensure that they are ordered in
  /// comparePositionsInEventOrder. Having events as values ensures
  /// that we can maintain a list of segments starting at the
  /// events position
  var _events = SplayTreeMap<_Event, _Event>();

  _EventQueue();

  /// Adds a new event at position [pos], unless there is already
  /// an event at this position in the queue.
  addEvent(Coordinate pos) {
    _require(pos != null);
    var evt = _Event(pos);
    _events.putIfAbsent(evt, () => evt);
  }

  /// Adds two events at the start and end position of [segment].
  addLineSegmentEvents(LineSegment segment) {
    _require(segment != null);
    addEvent(segment.start);
    addEvent(segment.end);
    _events[_Event(segment.start)].add(segment);
  }

  /// Returns true if this queue is empty
  bool get isEmpty => _events.isEmpty;

  /// Removes the first event from the queue and returns it.
  /// Throws an [StateError] if the queue is empty.
  _Event unshift() {
    if (isEmpty) throw StateError("queue is empty");
    var key = _events.firstKey();
    var ret = _events[key];
    _events.remove(key);
    return ret;
  }

  /// Returns the first event in this queue without removing it
  /// from this queue.
  ///
  /// Throws a [StateError] if the queue is empty.
  _Event get first {
    if (this.isEmpty) throw StateError("queue is empty");
    var evt = _events.firstKey();
    return _events[evt];
  }

  int get length => _events.length;
}

class _SweepLineCompareFunction {
  Coordinate event;
  int sign = 1;

  /// the compare function used for compare operations just *before*
  /// the sweep line goes through the event position
  _SweepLineCompareFunction.atMinusEpsilon(this.event) : sign = -1;

  /// the compare function used for compare operations just *after*
  /// the sweep line went through the event position
  _SweepLineCompareFunction.atPlusEpsilon(this.event) : sign = 1;

  int call(LineSegment value, LineSegment other) {
    // if other is left or right of the reference point
    // 'event' (which a priori is known to be on the segment
    // 'value' too), the ordering is clear
    var o = _orientation(other.start, other.end, event);
    if (o != 0) return o;

    // otherwise 'value' and 'other' intersect in the
    // point event and are ordered according to the counter clockwise orientation
    // value. If both values are identical, then 'value'
    // and 'other' overlap (i.e. all their endpoints are colineaer and the
    // intersection isn't empty). They are considered to be in an
    // equivalence class.
    return sign *
        value.counterclockwiseOrientation
            .compareTo(other.counterclockwiseOrientation);
  }
}

/// A [LineIntersection] represents te interesection of two or more
/// [LineSegment]s.
class LineIntersection {
  /// the position where the line segments intersect
  final Coordinate pos;
  List<LineSegment> _intersecting;

  /// the list of intersecting line segments
  List<LineSegment> get intersecting => UnmodifiableListView(_intersecting);

  LineIntersection._(this.pos, Iterable<LineSegment> intersecting) {
    _require(intersecting != null);
    this._intersecting = List.from(intersecting, growable: false);
  }

  toString() => "{intersection: pos=$pos, segments=${_intersecting.join(",")}";
}

/// Computes the line intersections for the collection of line segments in
/// [segments].
///
/// Returns an empty [Iterable], if no interesections have been detected in
/// [segments].
Iterable<LineIntersection> computeLineIntersections(
    Iterable<LineSegment> segments) {
  /*
   * Implementation and terminology closely follows
   *    M. van de Berg
   *    Computational Geometry
   *    Springer
   *    1998
   *    Chapter 2, Line Segment Intersection
   */

  final AvlTree<LineSegment> sweepLine =
      AvlTree<LineSegment>(withEquivalenceClasses: true);
  final eventQueue = _EventQueue();
  final intersections = List<LineIntersection>();

  /// initializes the event queue with the start events of the segments
  initEventQueue() {
    segments.forEach((s) {
      eventQueue.addLineSegmentEvents(s);
    });
  }

  Iterable flatten(e) =>
      e is Iterable ? e.expand(flatten) : new List.filled(1, e);

  /**
   * [sl] the left neighbours to look for intersections
   * [sr] the right neighbours to look for intersections
   *
   * This methods finds all intersections between line segments
   * in sl and line segments in sr. If an identified intersection
   * is identified as "future" event in the line sweeping algorithm,
   * it is added to the event queue.
   */
  findAndRememberEvent(sl, sr, pos) {
    if (sl == null || sr == null) return;
    sl.forEach((left) => sr.forEach((right) {
          var intersection = left.intersectionWith(right);
          if (intersection == null) return;
          if (_comparePositionsInEventOrder(intersection, pos) > 0) {
            eventQueue.addEvent(intersection);
          }
        }));
  }

  handleEvent(_Event event) {
    // the list of events starting at this point
    var U = List.from(event.segments, growable: false);

    // the list of segments intersecting with the sweepline at point
    // event.pos
    var t = sweepLine.inorder
        .where((s) => _orientation(s.first.start, s.first.end, event.pos) == 0)
        .expand((s) => s) // 'unwrap' the lists
        .toList(growable: false);

    // the segments ending in event.pos
    var L = t.where((s) => s.end == event.pos).toList(growable: false);

    // the segments properly containing event.pos (event.pos is neither
    // the start nor the end point)
    var C = t.where((s) => s.end != event.pos).toList(growable: false);

    if (U.length + L.length + C.length > 1) {
      // create a new line intersection at position event.pos
      // and the intersecting lines
      Set<LineSegment> intersecting = flatten([U, L, C]).toSet();
      intersections.add(LineIntersection._(event.pos, intersecting));
    }

    Function compare = _SweepLineCompareFunction.atMinusEpsilon(event.pos);

    L.forEach((s) => sweepLine.remove(s, compare: compare));
    C.forEach((s) => sweepLine.remove(s, compare: compare));

    compare = _SweepLineCompareFunction.atPlusEpsilon(event.pos);
    U.forEach((s) => sweepLine.add(s, compare: compare));
    C.forEach((s) => sweepLine.add(s, compare: compare));

    compare =
        (LineSegment other) => _orientation(other.start, other.end, event.pos);

    if (U.length + C.length == 0) {
      var sl = sweepLine.leftNeighbour(compare);
      var sr = sweepLine.rightNeighbour(compare);
      findAndRememberEvent(sl, sr, event.pos);
    } else {
      var UC = flatten([U, C]).toList(growable: false);
      var sl = sweepLine.leftNeighbour(compare);
      var sr = sweepLine.rightNeighbour(compare);
      compare = _SweepLineCompareFunction.atPlusEpsilon(event.pos);
      UC.sort(compare);
      findAndRememberEvent(sl, [UC.first], event.pos);
      findAndRememberEvent(sr, [UC.last], event.pos);
    }
  }

  initEventQueue();
  while (!eventQueue.isEmpty) {
    var event = eventQueue.unshift();
    handleEvent(event);
  }
  return intersections;
}
