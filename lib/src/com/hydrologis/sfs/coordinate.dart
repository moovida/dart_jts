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

  double _x;
  double _y;
  double _z;

  Coordinate(this._x, this._y, [this._z]);

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

  void setCoordinate(Coordinate other) {
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
  num getOrdinate(int ordinateIndex) {
    switch (ordinateIndex) {
      case X:
        return _x;
      case Y:
        return _y;
      case Z:
        return _z; // sure to delegate to subclass rather than offer direct field access
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
  void setOrdinate(int ordinateIndex, double value) {
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
    return ret != 0 ? ret : _y.compareTo(other._y);
  }

  ///  Returns whether the planar projections of the two <code>Coordinate</code>s
  ///  are equal.
  ///
  ///@param  other  a <code>Coordinate</code> with which to do the 2D comparison.
  ///@return        <code>true</code> if the x- and y-coordinates are equal; the
  ///      z-coordinates do not have to be equal.
  bool equals2D(Coordinate other) {
    if (x != other.x) {
      return false;
    }
    if (y != other.y) {
      return false;
    }
    return true;
  }

  Coordinate copy() {
    return Coordinate.fromCoordinate(this);
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

  double getM() {
    return double.nan;
  }

  bool hasM() => false;

  void setM(double m) {}
}

//class CoordinateList<E> extends ListBase<E> {
//  final List<E> l = [];
//  CoordinateList();
//
//  void set length(int newLength) { l.length = newLength; }
//  int get length => l.length;
//  E operator [](int index) => l[index];
//  void operator []=(int index, E value) { l[index] = value; }
//
//
//
//  Coordinate getCoordinate(int i) { return (Coordinate) [i]; }
//
//
//  /**
//   * Adds a section of an array of coordinates to the list.
//   * @param coord The coordinates
//   * @param allowRepeated if set to false, repeated coordinates are collapsed
//   * @param start the index to start from
//   * @param end the index to add up to but not including
//   * @return true (as by general collection contract)
//   */
//   bool addSection(List<Coordinate> coord, bool allowRepeated, int start, int end)
//  {
//  int inc = 1;
//  if (start > end) inc = -1;
//
//  for (int i = start; i != end; i += inc) {
//  add(coord[i], allowRepeated);
//  }
//  return true;
//  }
//
//  /// Adds an array of coordinates to the list.
//  /// @param coord The coordinates
//  /// @param allowRepeated if set to false, repeated coordinates are collapsed
//  /// @param direction if false, the array is added in reverse order
//  /// @return true (as by general collection contract)
//  bool addListAllowDirect(List<Coordinate> coord, bool allowRepeated, bool direction)
//  {
//  if (direction) {
//  for (int i = 0; i < coord.length; i++) {
//  add(coord[i], allowRepeated);
//  }
//  }
//  else {
//  for (int i = coord.length - 1; i >= 0; i--) {
//  add(coord[i], allowRepeated);
//  }
//  }
//  return true;
//  }
//
//
//  /// Adds an array of coordinates to the list.
//  /// @param coord The coordinates
//  /// @param allowRepeated if set to false, repeated coordinates are collapsed
//  /// @return true (as by general collection contract)
//   bool addAllow(List<Coordinate> coord, bool allowRepeated)
//  {
//    addListAllowDirect(coord, allowRepeated, true);
//  return true;
//  }
//
//  /**
//   * Adds a coordinate to the list.
//   * @param obj The coordinate to add
//   * @param allowRepeated if set to false, repeated coordinates are collapsed
//   * @return true (as by general collection contract)
//   */
//   bool addObjectAllow(Object obj, bool allowRepeated)
//  {
//    add((Coordinate) obj, allowRepeated);
//    return true;
//  }
//
//  /**
//   * Adds a coordinate to the end of the list.
//   *
//   * @param coord The coordinates
//   * @param allowRepeated if set to false, repeated coordinates are collapsed
//   */
//   void addCoordinate(Coordinate coord, bool allowRepeated)
//  {
//    // don't add duplicate coordinates
//    if (! allowRepeated) {
//      if (size() >= 1) {
//        Coordinate last = (Coordinate) get(size() - 1);
//        if (last.equals2D(coord)) return;
//      }
//    }
//    super.add(coord);
//  }
//
//  /**
//   * Inserts the specified coordinate at the specified position in this list.
//   *
//   * @param i the position at which to insert
//   * @param coord the coordinate to insert
//   * @param allowRepeated if set to false, repeated coordinates are collapsed
//   */
//   void addAtPosition(int i, Coordinate coord, bool allowRepeated)
//  {
//    // don't add duplicate coordinates
//    if (! allowRepeated) {
//      int size = size();
//      if (size > 0) {
//        if (i > 0) {
//          Coordinate prev = (Coordinate) get(i - 1);
//          if (prev.equals2D(coord)) return;
//        }
//        if (i < size) {
//          Coordinate next = (Coordinate) get(i);
//          if (next.equals2D(coord)) return;
//        }
//      }
//    }
//    super.add(i, coord);
//  }
//
//  /** Add an array of coordinates
//   * @param coll The coordinates
//   * @param allowRepeated if set to false, repeated coordinates are collapsed
//   * @return true (as by general collection contract)
//   */
//   bool addAllList(List<E extends Coordinate> coll, bool allowRepeated)
//  {
//  bool isChanged = false;
//  for (Iterator<? extends Coordinate> i = coll.iterator(); i.hasNext(); ) {
//  add(i.next(), allowRepeated);
//  isChanged = true;
//  }
//  return isChanged;
//  }
//
//  /**
//   * Ensure this coordList is a ring, by adding the start point if necessary
//   */
//   void closeRing()
//  {
//    if (length > 0) {
//      Coordinate duplicate = [0].copy();
//      addCoordinate(duplicate, false);
//    }
//  }
//
//  /**
//   * Creates an array containing the coordinates in this list,
//   * oriented in the given direction (forward or reverse).
//   *
//   * @param direction the direction value: true for forward, false for reverse
//   * @return an oriented array of coordinates
//   */
//   List<Coordinate> toCoordinateArray(bool isForward)
//  {
//    if (isForward) {
//      return (Coordinate[]) toArray(coordArrayType);
//    }
//    // construct reversed array
//    int size = size();
//    Coordinate[] pts = new Coordinate[size];
//    for (int i = 0; i < size; i++) {
//      pts[i] = get(size - i - 1);
//    }
//    return pts;
//  }
//
//  /**
//   * Returns a deep copy of this <tt>CoordinateList</tt> instance.
//   *
//   * @return a clone of this <tt>CoordinateList</tt> instance
//   */
//   Object clone() {
//    CoordinateList clone = (CoordinateList) super.clone();
//    for (int i = 0; i < this.size(); i++) {
//      clone.add(i, (Coordinate) this.get(i).clone());
//  }
//    return clone;
//  }
//}
