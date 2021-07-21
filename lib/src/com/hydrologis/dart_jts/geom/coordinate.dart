part of dart_jts;

/// A lightweight class used to store coordinates on the 2-dimensional Cartesian plane.
/// <p>
/// It is distinct from {@link Point}, which is a subclass of {@link Geometry}.
/// Unlike objects of type {@link Point} (which contain additional
/// information such as an envelope, a precision model, and spatial reference
/// system information), a <code>Coordinate</code> only contains ordinate values
/// and accessor methods. </p>
/// <p>
/// <code>Coordinate</code>s are two-dimensional points, with an additional Z-ordinate.
/// If an Z-ordinate value is not specified or not defined,
/// constructed coordinates have a Z-ordinate of <code>NaN</code>
/// (which is also the value of <code>NULL_ORDINATE</code>).
/// The standard comparison functions ignore the Z-ordinate.
/// Apart from the basic accessor functions, JTS supports
/// only specific operations involving the Z-ordinate.</p>
/// <p>
/// Implementations may optionally support Z-ordinate and M-measure values
/// as appropriate for a {@link CoordinateSequence}.
/// Use of {@link #getZ()} and {@link #getM()}
/// accessors, or {@link #getOrdinate(int)} are recommended.</p>
///
/// @version 1.16
class Coordinate implements Comparable<Coordinate> {
  /// The value used to indicate a null or missing ordinate value.
  /// In particular, used for the value of ordinates for dimensions
  /// greater than the defined dimension of a coordinate.
  static final double NULL_ORDINATE = double.nan;

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

  /// The x-ordinate.
  late double x;

  /// The y-ordinate.
  late double y;

  /// The z-ordinate.
  /// <p>
  /// Direct access to this field is discouraged; use {@link #getZ()}.
  late double z;

  ///  Constructs a <code>Coordinate</code> at (x,y,z).
  ///
  ///@param  x  the x-ordinate
  ///@param  y  the y-ordinate
  ///@param  z  the z-ordinate
  Coordinate.fromXYZ(double x, double y, double z) {
    this.x = x;
    this.y = y;
    this.z = z;
  }

  ///  Constructs a <code>Coordinate</code> at (0,0,NaN).
  Coordinate.empty2D() : this(0.0, 0.0);

  ///  Constructs a <code>Coordinate</code> having the same (x,y,z) values as
  ///  <code>other</code>.
  ///
  ///@param  c  the <code>Coordinate</code> to copy.
  Coordinate.fromCoordinate(Coordinate c) : this.fromXYZ(c.x, c.y, c.getZ());

  ///  Constructs a <code>Coordinate</code> at (x,y,NaN).
  ///
  ///@param  x  the x-value
  ///@param  y  the y-value
  Coordinate(double x, double y) : this.fromXYZ(x, y, NULL_ORDINATE);

  ///  Sets this <code>Coordinate</code>s (x,y,z) values to that of <code>other</code>.
  ///
  ///@param  other  the <code>Coordinate</code> to copy
  void setCoordinate(Coordinate other) {
    x = other.x;
    y = other.y;
    z = other.getZ();
  }

  ///  Retrieves the value of the X ordinate.
  ///
  ///  @return the value of the X ordinate
  double getX() {
    return x;
  }

  /// Sets the X ordinate value.
  ///
  /// @param x the value to set as X
  void setX(double x) {
    this.x = x;
  }

  ///  Retrieves the value of the Y ordinate.
  ///
  ///  @return the value of the Y ordinate
  double getY() {
    return y;
  }

  /// Sets the Y ordinate value.
  ///
  /// @param y the value to set as Y
  void setY(double y) {
    this.y = y;
  }

  ///  Retrieves the value of the Z ordinate, if present.
  ///  If no Z value is present returns <tt>NaN</tt>.
  ///
  ///  @return the value of the Z ordinate, or <tt>NaN</tt>
  double getZ() {
    return z;
  }

  /// Sets the Z ordinate value.
  ///
  /// @param z the value to set as Z
  void setZ(double z) {
    this.z = z;
  }

  ///  Retrieves the value of the measure, if present.
  ///  If no measure value is present returns <tt>NaN</tt>.
  ///
  ///  @return the value of the measure, or <tt>NaN</tt>
  double getM() {
    return double.nan;
  }

  /// Sets the measure value, if supported.
  ///
  /// @param m the value to set as M
  void setM(double m) {
    throw ArgumentError("Invalid ordinate index: $M");
  }

  /// Gets the ordinate value for the given index.
  ///
  /// The base implementation supports values for the index are
  /// {@link X}, {@link Y}, and {@link Z}.
  ///
  /// @param ordinateIndex the ordinate index
  /// @return the value of the ordinate
  /// @throws IllegalArgumentException if the index is not valid
  double getOrdinate(int ordinateIndex) {
    switch (ordinateIndex) {
      case X:
        return x;
      case Y:
        return y;
      case Z:
        return getZ(); // sure to delegate to subclass rather than offer direct field access
    }
    throw new ArgumentError("Invalid ordinate index: $ordinateIndex");
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
        setZ(
            value); // delegate to subclass rather than offer direct field access
        break;
      default:
        throw new ArgumentError("Invalid ordinate index: $ordinateIndex");
    }
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

  /// Tests if another Coordinate has the same values for the X and Y ordinates,
  /// within a specified tolerance value.
  /// The Z ordinate is ignored.
  ///
  ///@param c a <code>Coordinate</code> with which to do the 2D comparison.
  ///@param tolerance the tolerance value to use
  ///@return true if <code>other</code> is a <code>Coordinate</code>
  ///      with the same values for X and Y.
  bool equals2DWithTolerance(Coordinate c, double tolerance) {
    if (!NumberUtils.equalsWithTolerance(this.x, c.x, tolerance)) {
      return false;
    }
    if (!NumberUtils.equalsWithTolerance(this.y, c.y, tolerance)) {
      return false;
    }
    return true;
  }

  /// Tests if another coordinate has the same values for the X, Y and Z ordinates.
  ///
  ///@param other a <code>Coordinate</code> with which to do the 3D comparison.
  ///@return true if <code>other</code> is a <code>Coordinate</code>
  ///      with the same values for X, Y and Z.
  bool equals3D(Coordinate other) {
    return (x == other.x) &&
        (y == other.y) &&
        ((getZ() == other.getZ()) || (getZ().isNaN && other.getZ().isNaN));
  }

  /// Tests if another coordinate has the same value for Z, within a tolerance.
  ///
  /// @param c a coordinate
  /// @param tolerance the tolerance value
  /// @return true if the Z ordinates are within the given tolerance
  bool equalInZ(Coordinate c, double tolerance) {
    return NumberUtils.equalsWithTolerance(this.getZ(), c.getZ(), tolerance);
  }

  ///  Returns <code>true</code> if <code>other</code> has the same values for
  ///  the x and y ordinates.
  ///  Since Coordinates are 2.5D, this routine ignores the z value when making the comparison.
  ///
  ///@param  other  a <code>Coordinate</code> with which to do the comparison.
  ///@return        <code>true</code> if <code>other</code> is a <code>Coordinate</code>
  ///      with the same values for the x and y ordinates.
  bool equals(Object other) {
    if (!(other is Coordinate)) {
      return false;
    }
    return equals2D(other);
  }

  ///  Compares this {@link Coordinate} with the specified {@link Coordinate} for order.
  ///  This method ignores the z value when making the comparison.
  ///  Returns:
  ///  <UL>
  ///    <LI> -1 : this.x &lt; other.x || ((this.x == other.x) &amp;&amp; (this.y &lt; other.y))
  ///    <LI> 0 : this.x == other.x &amp;&amp; this.y = other.y
  ///    <LI> 1 : this.x &gt; other.x || ((this.x == other.x) &amp;&amp; (this.y &gt; other.y))
  ///
  ///  </UL>
  ///  Note: This method assumes that ordinate values
  /// are valid numbers.  NaN values are not handled correctly.
  ///
  ///@param  o  the <code>Coordinate</code> with which this <code>Coordinate</code>
  ///      is being compared
  ///@return    -1, zero, or 1 as this <code>Coordinate</code>
  ///      is less than, equal to, or greater than the specified <code>Coordinate</code>
  int compareTo(Coordinate other) {
    if (x < other.x) return -1;
    if (x > other.x) return 1;
    if (y < other.y) return -1;
    if (y > other.y) return 1;
    return 0;
  }

  ///  Returns a <code>String</code> of the form <I>(x,y,z)</I> .
  ///
  ///@return    a <code>String</code> of the form <I>(x,y,z)</I>
  String toString() {
    return "($x, $y, ${getZ()})";
  }

  Object clone() {
    Coordinate coord = Coordinate.fromCoordinate(this);
    return coord; // return the clone
  }

  /// Creates a copy of this Coordinate.
  ///
  /// @return a copy of this coordinate.
  Coordinate copy() {
    return Coordinate.fromCoordinate(this);
  }

  /// Computes the 2-dimensional Euclidean distance to another location.
  /// The Z-ordinate is ignored.
  ///
  /// @param c a point
  /// @return the 2-dimensional Euclidean distance between the locations
  double distance(Coordinate c) {
    double dx = x - c.x;
    double dy = y - c.y;
    return math.sqrt(dx * dx + dy * dy);
  }

  /// Computes the 3-dimensional Euclidean distance to another location.
  ///
  /// @param c a coordinate
  /// @return the 3-dimensional Euclidean distance between the locations
  double distance3D(Coordinate c) {
    double dx = x - c.x;
    double dy = y - c.y;
    double dz = getZ() - c.getZ();
    return math.sqrt(dx * dx + dy * dy + dz * dz);
  }

  int get hashCode {
    //Algorithm from Effective Java by Joshua Bloch [Jon Aquino]
    const prime = 31;
    int result = 1;
    result = prime * result + x.hashCode;
    result = prime * result + y.hashCode;
    return result;
  }

  Coordinate operator -(Coordinate other) {
    return Coordinate(x - other.x, y - other.y);
  }

  /// Returns the (vector-) sum of this and [other].
  Coordinate operator +(Coordinate other) {
    return Coordinate(x + other.x, y + other.y);
  }

  bool operator ==(o) {
    if (o is Coordinate) {
      if (!equals3D(o)) {
        return false;
      }
      if (((getM() == o.getM()) || (getM().isNaN && o.getM().isNaN))) {
        return true;
      }
    }
    return false;
  }
}
