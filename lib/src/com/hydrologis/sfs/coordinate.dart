part of dart_sfs;

class Coordinate implements Comparable {
  /// The value used to indicate a null or missing ordinate value.
  /// In particular, used for the value of ordinates for dimensions
  /// greater than the defined dimension of a coordinate.
  static const double NULL_ORDINATE = double.nan;

  /// Standard ordinate index value for, where X is 0 */
  static const int X = 0;

  /// Standard ordinate index value for, where Y is 1 */
  static const int Y = 1;

  /// Standard ordinate index value for, where Z is 2.
  ///
  /// <p>This constant assumes XYZM coordinate sequence definition, please check this assumption
  /// using {@link #getDimension()} and {@link #getMeasures()} before use.
  static const int Z = 2;

  /// Standard ordinate index value for, where M is 3.
  ///
  /// <p>This constant assumes XYZM coordinate sequence definition, please check this assumption
  /// using {@link #getDimension()} and {@link #getMeasures()} before use.
  static const int M = 3;

  num _x;
  num _y;
  num _z;

  Coordinate(this._x, this._y);

  Coordinate.empty() {
    _x = 0.0;
    _y = 0.0;
    _z = NULL_ORDINATE;
  }

  Coordinate.fromCoordinate(Coordinate c) {
    _x = c._x;
    _y = c._y;
    _z = c._z;
  }

  num get x => _x;

  set x(newX) => _x = newX;

  num get y => _y;

  set y(newY) => _y = newY;

  num get z => _z;

  set z(newZ) => _z = newZ;


  void setCoordinate(Coordinate other)
  {
    x = other.x;
    y = other.y;
    z = other.z;
  }


  /// Gets the ordinate value for the given index.
  ///
  /// The base implementation supports values for the index are
  /// {@link X}, {@link Y}, and {@link Z}.
  ///
  /// @param ordinateIndex the ordinate index
  ///\ @return the value of the ordinate
  /// @throws IllegalArgumentException if the index is not valid
  num getOrdinate(int ordinateIndex)
  {
    switch (ordinateIndex) {
      case X: return _x;
      case Y: return _y;
      case Z: return _z; // sure to delegate to subclass rather than offer direct field access
    }
    throw ArgumentError("Invalid ordinate index: $ordinateIndex");
  }

  /// Sets the ordinate for the given index
  /// to a given value.
  ///
  /// The base implementation supported values for the index are
  /// {@link X}, {@link Y}, and {@link Z}.
  ///
  /// @param ordinateIndex the ordinate index
  /// @param value the value to set
  /// @throws IllegalArgumentException if the index is not valid
  void setOrdinate(int ordinateIndex, double value)
  {
    switch (ordinateIndex) {
      case X:
        x = value;
        break;
      case Y:
        y = value;
        break;
      case Z:
        _z = value;
        break;
      default:
        throw ArgumentError("Invalid ordinate index: $ordinateIndex");
    }
  }

  int compareTo(other) {
    _require(other is Coordinate);
    int ret = _x.compareTo(other._x);
    return ret != 0 ? ret : _y.compareTo(other.y_);
  }

  String toString() => "($_x, $_y, $_z)";

  int get hashCode {
    const prime = 31;
    int result = 1;
    result = prime * result + _x.hashCode;
    result = prime * result + _y.hashCode;
    return result;
  }

  bool operator ==(other) => compareTo(other) == 0;

  bool operator <(other) => compareTo(other) == -1;

  bool operator <=(other) => compareTo(other) <= 0;

  bool operator >(other) => compareTo(other) == 1;

  bool operator >=(other) => compareTo(other) >= 0;

  /// Returns the (vector-) difference of this and [other].
  Coordinate operator -(Coordinate other) {
    return Coordinate(_x - other._x, _y - other._y);
  }

  /// Returns the (vector-) sum of this and [other].
  Coordinate operator +(Coordinate other) {
    return Coordinate(_x + other._x, _y + other._y);
  }

  /// Replies this position scaled by [factor]
  Coordinate scale(num factor) {
    return Coordinate(_x * factor, _y * factor);
  }
}
