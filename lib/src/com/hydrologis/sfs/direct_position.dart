part of dart_sfs;

class DirectPosition2D implements Comparable {
  final num  x;
  final num y;
  const DirectPosition2D(this.x, this.y);

  int compareTo(other) {
    _require(other is DirectPosition2D);
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
  DirectPosition2D operator -(DirectPosition2D other) {
    return new DirectPosition2D(x - other.x, y - other.y);
  }

  /// Returns the (vector-) sum of this and [other].
  DirectPosition2D operator +(DirectPosition2D other) {
    return new DirectPosition2D(x + other.x, y + other.y);
  }

  /// Replies this position scaled by [factor]
  DirectPosition2D scale(num factor) {
    return new DirectPosition2D(x * factor, y * factor);
  }
}