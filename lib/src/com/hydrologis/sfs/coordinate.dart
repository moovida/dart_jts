part of dart_sfs;

class Coordinate implements Comparable {
  final num x;
  final num y;

  const Coordinate(this.x, this.y);

  int compareTo(other) {
    _require(other is Coordinate);
    int ret = x.compareTo(other.x);
    return ret != 0 ? ret : y.compareTo(other.y);
  }

  String toString() => "($x, $y)";

  int get hashCode {
    const prime = 31;
    int result = 1;
    result = prime * result + x.hashCode;
    result = prime * result + y.hashCode;
    return result;
  }

  bool operator ==(other) => compareTo(other) == 0;

  bool operator <(other) => compareTo(other) == -1;

  bool operator <=(other) => compareTo(other) <= 0;

  bool operator >(other) => compareTo(other) == 1;

  bool operator >=(other) => compareTo(other) >= 0;

  /// Returns the (vector-) difference of this and [other].
  Coordinate operator -(Coordinate other) {
    return Coordinate(x - other.x, y - other.y);
  }

  /// Returns the (vector-) sum of this and [other].
  Coordinate operator +(Coordinate other) {
    return Coordinate(x + other.x, y + other.y);
  }

  /// Replies this position scaled by [factor]
  Coordinate scale(num factor) {
    return Coordinate(x * factor, y * factor);
  }
}
