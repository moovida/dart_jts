part of dart_jts;

/// Indicates an invalid or inconsistent topological situation encountered during processing
///
/// @version 1.7
class TopologyException implements Exception {
  static String msgWithCoord(String msg, Coordinate? pt) {
    if (pt != null) {
      return msg + " [ $pt ]";
    }
    return msg;
  }

  Coordinate? pt;

  String? msg;

  TopologyException(this.msg);

  TopologyException.withCoord(String msg, Coordinate pt) {
    this.msg = msgWithCoord(msg, pt);
    this.pt = Coordinate.fromCoordinate(pt);
  }

  Coordinate? getCoordinate() {
    return pt;
  }
}

/// Compares two {@link CoordinateSequence}s.
/// For sequences of the same dimension, the ordering is lexicographic.
/// Otherwise, lower dimensions are sorted before higher.
/// The dimensions compared can be limited; if this is done
/// ordinate dimensions above the limit will not be compared.
/// <p>
/// If different behaviour is required for comparing size, dimension, or
/// coordinate values, any or all methods can be overridden.
///
class CoordinateSequenceComparatorBuilder {
  static Comparator<CoordinateSequence> regular() {
    return (cs1, cs2) {
      return compare(cs1, cs2, NumberUtils.MAX_INT);
    };
  }

  static Comparator<CoordinateSequence> withLimit(int dimensionLimit) {
    return (cs1, cs2) {
      return compare(cs1, cs2, dimensionLimit);
    };
  }

  /// Compare two <code>double</code>s, allowing for NaN values.
  /// NaN is treated as being less than any valid number.
  ///
  /// @param a a <code>double</code>
  /// @param b a <code>double</code>
  /// @return -1, 0, or 1 depending on whether a is less than, equal to or greater than b
  static int compareStatic(double a, double b) {
    if (a < b) return -1;
    if (a > b) return 1;

    if (a.isNaN) {
      if (b.isNaN) return 0;
      return -1;
    }

    if (b.isNaN) return 1;
    return 0;
  }

  /// Compares two {@link CoordinateSequence}s for relative order.
  ///
  /// @param o1 a {@link CoordinateSequence}
  /// @param o2 a {@link CoordinateSequence}
  /// @return -1, 0, or 1 depending on whether o1 is less than, equal to, or greater than o2
  static int compare(Object o1, Object o2, int dimensionLimit) {
    CoordinateSequence s1 = o1 as CoordinateSequence;
    CoordinateSequence s2 = o2 as CoordinateSequence;

    int size1 = s1.size();
    int size2 = s2.size();

    int dim1 = s1.getDimension();
    int dim2 = s2.getDimension();

    int minDim = dim1;
    if (dim2 < minDim) minDim = dim2;
    bool dimLimited = false;
    if (dimensionLimit <= minDim) {
      minDim = dimensionLimit;
      dimLimited = true;
    }

    // lower dimension is less than higher
    if (!dimLimited) {
      if (dim1 < dim2) return -1;
      if (dim1 > dim2) return 1;
    }

    // lexicographic ordering of point sequences
    int i = 0;
    while (i < size1 && i < size2) {
      int ptComp = compareCoordinate(s1, s2, i, minDim);
      if (ptComp != 0) return ptComp;
      i++;
    }
    if (i < size1) return 1;
    if (i < size2) return -1;

    return 0;
  }

  /// Compares the same coordinate of two {@link CoordinateSequence}s
  /// along the given number of dimensions.
  ///
  /// @param s1 a {@link CoordinateSequence}
  /// @param s2 a {@link CoordinateSequence}
  /// @param i the index of the coordinate to test
  /// @param dimension the number of dimensions to test
  /// @return -1, 0, or 1 depending on whether s1[i] is less than, equal to, or greater than s2[i]
  static int compareCoordinate(
      CoordinateSequence s1, CoordinateSequence s2, int i, int dimension) {
    for (int d = 0; d < dimension; d++) {
      double ord1 = s1.getOrdinate(i, d);
      double ord2 = s2.getOrdinate(i, d);
      int comp = compareStatic(ord1, ord2);
      if (comp != 0) return comp;
    }
    return 0;
  }
}

///  An interface for classes which process the coordinates in a {@link CoordinateSequence}.
///  A filter can either record information about each coordinate,
///  or change the value of the coordinate.
///  Filters can be
///  used to implement operations such as coordinate transformations, centroid and
///  envelope computation, and many other functions.
///  {@link Geometry} classes support the concept of applying a
///  <code>CoordinateSequenceFilter</code> to each
///  {@link CoordinateSequence}s they contain.
///  <p>
///  For maximum efficiency, the execution of filters can be short-circuited by using the {@link #isDone} method.
///  <p>
///  <code>CoordinateSequenceFilter</code> is
///  an example of the Gang-of-Four Visitor pattern.
///  <p>
/// <b>Note</b>: In general, it is preferable to treat Geometrys as immutable.
/// Mutation should be performed by creating a new Geometry object (see {@link GeometryEditor}
/// and {@link GeometryTransformer} for convenient ways to do this).
/// An exception to this rule is when a new Geometry has been created via {@link Geometry#copy()}.
/// In this case mutating the Geometry will not cause aliasing issues,
/// and a filter is a convenient way to implement coordinate transformation.
///
/// @see Geometry#apply(CoordinateFilter)
/// @see GeometryTransformer
/// @see GeometryEditor
///
///@see Geometry#apply(CoordinateSequenceFilter)
///@author Martin Davis
///@version 1.7
abstract class CoordinateSequenceFilter {
  /// Performs an operation on a coordinate in a {@link CoordinateSequence}.
  ///
  ///@param seq  the <code>CoordinateSequence</code> to which the filter is applied
  ///@param i the index of the coordinate to apply the filter to
  void filter(CoordinateSequence seq, int i);

  /// Reports whether the application of this filter can be terminated.
  /// Once this method returns <tt>true</tt>, it must
  /// continue to return <tt>true</tt> on every subsequent call.
  ///
  /// @return true if the application of this filter can be terminated.
  bool isDone();

  /// Reports whether the execution of this filter
  /// has modified the coordinates of the geometry.
  /// If so, {@link Geometry#geometryChanged} will be executed
  /// after this filter has finished being executed.
  /// <p>
  /// Most filters can simply return a constant value reflecting
  /// whether they are able to change the coordinates.
  ///
  /// @return true if this filter has changed the coordinates of the geometry
  bool isGeometryChanged();
}

///  An interface for classes which use the values of the coordinates in a {@link Geometry}.
/// Coordinate filters can be used to implement centroid and
/// envelope computation, and many other functions.
/// <p>
/// <code>CoordinateFilter</code> is
/// an example of the Gang-of-Four Visitor pattern.
/// <p>
/// <b>Note</b>: it is not recommended to use these filters to mutate the coordinates.
/// There is no guarantee that the coordinate is the actual object stored in the source geometry.
/// In particular, modified values may not be preserved if the source Geometry uses a non-default {@link CoordinateSequence}.
/// If in-place mutation is required, use {@link CoordinateSequenceFilter}.
///
/// @see Geometry#apply(CoordinateFilter)
/// @see CoordinateSequenceFilter
///
///@version 1.7
//abstract class CoordinateFilter {
//  /// Performs an operation with the provided <code>coord</code>.
//  /// Note that there is no guarantee that the input coordinate
//  /// is the actual object stored in the source geometry,
//  /// so changes to the coordinate object may not be persistent.
//  ///
//  ///@param  coord  a <code>Coordinate</code> to which the filter is applied.
//  void filter(Coordinate coord);
//}

//typedef CoordinateFilter<Coordinate> = void Function(Coordinate coordinate);
abstract class CoordinateFilter {
  void filter(Coordinate? coordinate);
}

///  <code>GeometryCollection</code> classes support the concept of
///  applying a <code>GeometryFilter</code> to the <code>Geometry</code>.
///  The filter is applied to every element <code>Geometry</code>.
///  A <code>GeometryFilter</code> can either record information about the <code>Geometry</code>
///  or change the <code>Geometry</code> in some way.
///  <code>GeometryFilter</code>
///  is an example of the Gang-of-Four Visitor pattern.
///
///@version 1.7
abstract class GeometryFilter {
  ///  Performs an operation with or on <code>geom</code>.
  ///
  ///@param  geom  a <code>Geometry</code> to which the filter is applied.
  void filter(Geometry geom);
}

/// Identifies {@link Geometry} subclasses which
/// are 0-dimensional and with components which are {@link Point}s.
///
/// @author Martin Davis
///
abstract class Puntal {}

/// Identifies {@link Geometry} subclasses which
/// are 1-dimensional and have components which are {@link LineString}s.
///
/// @author Martin Davis
///
abstract class Lineal {}

/// Identifies {@link Geometry} subclasses which
/// are 2-dimensional
/// and have components which have {@link Lineal} boundaries.
///
/// @author Martin Davis
///
abstract class Polygonal {}

///  Iterates over all {@link Geometry}s in a {@link Geometry},
///  (which may be either a collection or an atomic geometry).
///  The iteration sequence follows a pre-order, depth-first traversal of the
///  structure of the <code>GeometryCollection</code>
///  (which may be nested). The original <code>Geometry</code> object is
///  returned as well (as the first object), as are all sub-collections and atomic elements.
///  It is  simple to ignore the intermediate <code>GeometryCollection</code> objects if they are not
///  needed.
///
///@version 1.7
class GeometryCollectionIterator implements Iterator {
  ///  The <code>Geometry</code> being iterated over.
  late Geometry parent;

  ///  Indicates whether or not the first element
  ///  (the root <code>GeometryCollection</code>) has been returned.
  late bool atStart;

  ///  The number of <code>Geometry</code>s in the the <code>GeometryCollection</code>.
  late int max;

  ///  The index of the <code>Geometry</code> that will be returned when <code>next</code>
  ///  is called.
  late int index;

  ///  The iterator over a nested <code>Geometry</code>, or <code>null</code>
  ///  if this <code>GeometryCollectionIterator</code> is not currently iterating
  ///  over a nested <code>GeometryCollection</code>.
  GeometryCollectionIterator? subcollectionIterator;

  ///  Constructs an iterator over the given <code>Geometry</code>.
  ///
  ///@param  parent  the geometry over which to iterate; also, the first
  ///      element returned by the iterator.
  GeometryCollectionIterator(Geometry parent) {
    this.parent = parent;
    atStart = true;
    index = 0;
    max = parent.getNumGeometries();
  }

  @override
  get current => next();

  @override
  bool moveNext() {
    return hasNext();
  }

  /// Tests whether any geometry elements remain to be returned.
  ///
  /// @return true if more geometry elements remain
  bool hasNext() {
    if (atStart) {
      return true;
    }
    if (subcollectionIterator != null) {
      if (subcollectionIterator!.hasNext()) {
        return true;
      }
      subcollectionIterator = null;
    }
    if (index >= max) {
      return false;
    }
    return true;
  }

  /// Gets the next geometry in the iteration sequence.
  ///
  /// @return the next geometry in the iteration
  Object next() {
    // the parent GeometryCollection is the first object returned
    if (atStart) {
      atStart = false;
      if (isAtomic(parent)) index++;
      return parent;
    }
    if (subcollectionIterator != null) {
      if (subcollectionIterator!.hasNext()) {
        return subcollectionIterator!.next();
      } else {
        subcollectionIterator = null;
      }
    }
    if (index >= max) {
      throw RangeError("index >= max");
    }
    Geometry obj = parent.getGeometryN(index++);
    if (obj is GeometryCollection) {
      subcollectionIterator = GeometryCollectionIterator(obj);
      // there will always be at least one element in the sub-collection
      return subcollectionIterator!.next();
    }
    return obj;
  }

  static bool isAtomic(Geometry geom) {
    return !(geom is GeometryCollection);
  }

  /// Removal is not supported.
  ///
  /// @throws  UnsupportedOperationException  This method is not implemented.
  void remove() {
    throw UnsupportedError(this.runtimeType.toString());
  }
}

/**
 * Coordinate subclass supporting XY ordinates.
 * <p>
 * This data object is suitable for use with coordinate sequences with <tt>dimension</tt> = 2.
 * <p>
 * The {@link Coordinate#z} field is visible, but intended to be ignored.
 *
 * @since 1.16
 */
class CoordinateXY extends Coordinate {
  /** Standard ordinate index value for X */
  static const int X = 0;

  /** Standard ordinate index value for Y */
  static const int Y = 1;

  /** CoordinateXY does not support Z values. */
  static const int Z = -1;

  /** CoordinateXY does not support M measures. */
  static const int M = -1;

  /** Default constructor */
  CoordinateXY() : super.empty2D();

  /**
   * Constructs a CoordinateXY instance with the given ordinates.
   *
   * @param x the X ordinate
   * @param y the Y ordinate
   */
  CoordinateXY.fromXY(double x, double y)
      : super.fromXYZ(x, y, Coordinate.NULL_ORDINATE);

  /**
   * Constructs a CoordinateXY instance with the x and y ordinates of the given Coordinate.
   *
   * @param coord the Coordinate providing the ordinates
   */
  CoordinateXY.fromCoordinate(Coordinate coord) : super(coord.x, coord.y);

  /**
   * Constructs a CoordinateXY instance with the x and y ordinates of the given CoordinateXY.
   *
   * @param coord the CoordinateXY providing the ordinates
   */
  CoordinateXY.fromCoordinateXY(CoordinateXY coord) : super(coord.x, coord.y);

  /**
   * Creates a copy of this CoordinateXY.
   *
   * @return a copy of this CoordinateXY
   */
  CoordinateXY copy() {
    return new CoordinateXY.fromCoordinateXY(this);
  }

  /** The z-ordinate is not supported */
  double getZ() {
    return Coordinate.NULL_ORDINATE;
  }

  /** The z-ordinate is not supported */
  void setZ(double z) {
    throw new ArgumentError(
        "CoordinateXY dimension 2 does not support z-ordinate");
  }

  void setCoordinate(Coordinate other) {
    x = other.x;
    y = other.y;
    z = other.getZ();
  }

  double getOrdinate(int ordinateIndex) {
    switch (ordinateIndex) {
      case X:
        return x;
      case Y:
        return y;
    }
    throw new ArgumentError("Invalid ordinate index: $ordinateIndex");
  }

  void setOrdinate(int ordinateIndex, double value) {
    switch (ordinateIndex) {
      case X:
        x = value;
        break;
      case Y:
        y = value;
        break;
      default:
        throw new ArgumentError("Invalid ordinate index: $ordinateIndex");
    }
  }

  String toString() {
    return "(" + x.toString() + ", " + y.toString() + ")";
  }
}

/**
 * Coordinate subclass supporting XYM ordinates.
 * <p>
 * This data object is suitable for use with coordinate sequences with <tt>dimension</tt> = 3 and <tt>measures</tt> = 1.
 * <p>
 * The {@link Coordinate#z} field is visible, but intended to be ignored.
 *
 * @since 1.16
 */
class CoordinateXYM extends Coordinate {
  /** Standard ordinate index value for X */
  static const int X = 0;

  /** Standard ordinate index value for Y */
  static const int Y = 1;

  /** CoordinateXYM does not support Z values. */
  static const int Z = -1;

  /**
   * Standard ordinate index value for M in XYM sequences.
   *
   * <p>This constant assumes XYM coordinate sequence definition.  Check this assumption using
   * {@link #getDimension()} and {@link #getMeasures()} before use.
   */
  static const int M = 2;

  /** Default constructor */
  CoordinateXYM.empty() : super.empty2D() {
    this.m = 0.0;
  }

  /**
   * Constructs a CoordinateXYM instance with the given ordinates and measure.
   *
   * @param x the X ordinate
   * @param y the Y ordinate
   * @param m the M measure value
   */
  CoordinateXYM(double x, double y, double m)
      : super.fromXYZ(x, y, Coordinate.NULL_ORDINATE) {
    this.m = m;
  }

  /**
   * Constructs a CoordinateXYM instance with the x and y ordinates of the given Coordinate.
   *
   * @param coord the coordinate providing the ordinates
   */
  CoordinateXYM.fromCoordinate(Coordinate coord) : super(coord.x, coord.y) {
    m = getM();
  }

  /**
   * Constructs a CoordinateXY instance with the x and y ordinates of the given CoordinateXYM.
   *
   * @param coord the coordinate providing the ordinates
   */
  CoordinateXYM.fromCoordinateXYM(CoordinateXYM coord)
      : super(coord.x, coord.y) {
    m = coord.m;
  }

  /**
   * Creates a copy of this CoordinateXYM.
   *
   * @return a copy of this CoordinateXYM
   */
  CoordinateXYM copy() {
    return new CoordinateXYM.fromCoordinateXYM(this);
  }

  /** The m-measure. */
  double m = 0;

  /** The m-measure, if available. */
  double getM() {
    return m;
  }

  void setM(double m) {
    this.m = m;
  }

  /** The z-ordinate is not supported */
  double getZ() {
    return Coordinate.NULL_ORDINATE;
  }

  /** The z-ordinate is not supported */
  void setZ(double z) {
    throw new ArgumentError(
        "CoordinateXY dimension 2 does not support z-ordinate");
  }

  void setCoordinate(Coordinate other) {
    x = other.x;
    y = other.y;
    z = other.getZ();
    m = other.getM();
  }

  double getOrdinate(int ordinateIndex) {
    switch (ordinateIndex) {
      case X:
        return x;
      case Y:
        return y;
      case M:
        return m;
    }
    throw new ArgumentError("Invalid ordinate index: $ordinateIndex");
  }

  void setOrdinate(int ordinateIndex, double value) {
    switch (ordinateIndex) {
      case X:
        x = value;
        break;
      case Y:
        y = value;
        break;
      case M:
        m = value;
        break;
      default:
        throw new ArgumentError("Invalid ordinate index: $ordinateIndex");
    }
  }

  String toString() {
    return "(" +
        x.toString() +
        ", " +
        y.toString() +
        " m=" +
        getM().toString() +
        ")";
  }
}

/**
 * Coordinate subclass supporting XYZM ordinates.
 * <p>
 * This data object is suitable for use with coordinate sequences with <tt>dimension</tt> = 4 and <tt>measures</tt> = 1.
 *
 * @since 1.16
 */
class CoordinateXYZM extends Coordinate {
  /** Default constructor */
  CoordinateXYZM.empty() : super.empty2D() {
    this.m = 0.0;
  }

  /**
   * Constructs a CoordinateXYZM instance with the given ordinates and measure.
   *
   * @param x the X ordinate
   * @param y the Y ordinate
   * @param z the Z ordinate
   * @param m the M measure value
   */
  CoordinateXYZM(double x, double y, double z, double m)
      : super.fromXYZ(x, y, z) {
    this.m = m;
  }

  /**
   * Constructs a CoordinateXYZM instance with the ordinates of the given Coordinate.
   *
   * @param coord the coordinate providing the ordinates
   */
  CoordinateXYZM.fromCoordinate(Coordinate coord)
      : super.fromCoordinate(coord) {
    m = getM();
  }

  /**
   * Constructs a CoordinateXYZM instance with the ordinates of the given CoordinateXYZM.
   *
   * @param coord the coordinate providing the ordinates
   */
  CoordinateXYZM.fromCoordinateXYZM(CoordinateXYZM coord)
      : super.fromCoordinate(coord) {
    m = coord.m;
  }

  /**
   * Creates a copy of this CoordinateXYZM.
   *
   * @return a copy of this CoordinateXYZM
   */
  CoordinateXYZM copy() {
    return new CoordinateXYZM.fromCoordinateXYZM(this);
  }

  /** The m-measure. */
  double m = 0;

  /** The m-measure, if available. */
  double getM() {
    return m;
  }

  void setM(double m) {
    this.m = m;
  }

  double getOrdinate(int ordinateIndex) {
    switch (ordinateIndex) {
      case Coordinate.X:
        return x;
      case Coordinate.Y:
        return y;
      case Coordinate.Z:
        return getZ(); // sure to delegate to subclass rather than offer direct field access
      case Coordinate.M:
        return getM(); // sure to delegate to subclass rather than offer direct field access
    }
    throw new ArgumentError("Invalid ordinate index: $ordinateIndex");
  }

  void setCoordinate(Coordinate other) {
    x = other.x;
    y = other.y;
    z = other.getZ();
    m = other.getM();
  }

  void setOrdinate(int ordinateIndex, double value) {
    switch (ordinateIndex) {
      case Coordinate.X:
        x = value;
        break;
      case Coordinate.Y:
        y = value;
        break;
      case Coordinate.Z:
        z = value;
        break;
      case Coordinate.M:
        m = value;
        break;
      default:
        throw new ArgumentError("Invalid ordinate index: $ordinateIndex");
    }
  }

  String toString() {
    return "(" +
        x.toString() +
        ", " +
        y.toString() +
        ", " +
        getZ().toString() +
        " m=" +
        getM().toString() +
        ")";
  }
}

/**
 * A list of {@link Coordinate}s, which may
 * be set to prevent repeated coordinates from occurring in the list.
 *
 *
 * @version 1.7
 */
class CoordinateList {
  late List<Coordinate> _backingList;

  /**
   * Constructs a new list without any coordinates
   */
  CoordinateList() {
    _backingList = [];
  }

  /**
   * Constructs a new list from an array of Coordinates, allowing repeated points.
   * (I.e. this constructor produces a {@link CoordinateList} with exactly the same set of points
   * as the input array.)
   *
   * @param coord the initial coordinates
   */
  CoordinateList.fromList(List<Coordinate> coord) {
    _backingList = List.from(coord);
    addList(coord, true);
  }

  /**
   * Constructs a new list from an array of Coordinates,
   * allowing caller to specify if repeated points are to be removed.
   *
   * @param coord the array of coordinates to load into the list
   * @param allowRepeated if <code>false</code>, repeated points are removed
   */
  CoordinateList.fromListAllowRepeat(
      List<Coordinate> coord, bool allowRepeated) {
    _backingList = List.from(coord);
    addList(coord, allowRepeated);
  }

  bool add(Coordinate coord) {
    _backingList.add(coord);
    return true;
  }

  Coordinate getCoordinate(int i) {
    return _backingList[i];
  }

  /**
   * Adds a section of an array of coordinates to the list.
   * @param coord The coordinates
   * @param allowRepeated if set to false, repeated coordinates are collapsed
   * @param start the index to start from
   * @param end the index to add up to but not including
   * @return true (as by general collection contract)
   */
  bool add4(List<Coordinate> coord, bool allowRepeated, int start, int end) {
    int inc = 1;
    if (start > end) inc = -1;

    for (int i = start; i != end; i += inc) {
      addCoord(coord[i], allowRepeated);
    }
    return true;
  }

  /**
   * Adds an array of coordinates to the list.
   * @param coord The coordinates
   * @param allowRepeated if set to false, repeated coordinates are collapsed
   * @param direction if false, the array is added in reverse order
   * @return true (as by general collection contract)
   */
  bool add3(List<Coordinate> coord, bool allowRepeated, bool direction) {
    if (direction) {
      for (int i = 0; i < coord.length; i++) {
        addCoord(coord[i], allowRepeated);
      }
    } else {
      for (int i = coord.length - 1; i >= 0; i--) {
        addCoord(coord[i], allowRepeated);
      }
    }
    return true;
  }

  /**
   * Adds an array of coordinates to the list.
   * @param coord The coordinates
   * @param allowRepeated if set to false, repeated coordinates are collapsed
   * @return true (as by general collection contract)
   */
  bool addList(List<Coordinate> coord, bool allowRepeated) {
    add3(coord, allowRepeated, true);
    return true;
  }

//  /**
//   * Adds a coordinate to the list.
//   * @param obj The coordinate to add
//   * @param allowRepeated if set to false, repeated coordinates are collapsed
//   * @return true (as by general collection contract)
//   */
//   bool add(Object obj, bool allowRepeated)
//  {
//    add((Coordinate) obj, allowRepeated);
//    return true;
//  }

  /**
   * Adds a coordinate to the end of the list.
   *
   * @param coord The coordinates
   * @param allowRepeated if set to false, repeated coordinates are collapsed
   */
  void addCoord(Coordinate coord, bool allowRepeated) {
    // don't add duplicate coordinates
    if (!allowRepeated) {
      if (_backingList.length >= 1) {
        Coordinate last = _backingList.last;
        if (last.equals2D(coord)) return;
      }
    }
    _backingList.add(coord);
  }

  /**
   * Inserts the specified coordinate at the specified position in this list.
   *
   * @param i the position at which to insert
   * @param coord the coordinate to insert
   * @param allowRepeated if set to false, repeated coordinates are collapsed
   */
  void add33(int i, Coordinate coord, bool allowRepeated) {
    // don't add duplicate coordinates
    if (!allowRepeated) {
      int size = _backingList.length;
      if (size > 0) {
        if (i > 0) {
          Coordinate prev = _backingList[i - 1];
          if (prev.equals2D(coord)) return;
        }
        if (i < size) {
          Coordinate next = _backingList[i];
          if (next.equals2D(coord)) return;
        }
      }
    }
    _backingList.insert(i, coord);
  }

  /** Add an array of coordinates
   * @param coll The coordinates
   * @param allowRepeated if set to false, repeated coordinates are collapsed
   * @return true (as by general collection contract)
   */
  bool addAll(List<Coordinate> coll, bool allowRepeated) {
    bool isChanged = false;
    for (var c in coll) {
      addCoord(c, allowRepeated);
      isChanged = true;
    }

    return isChanged;
  }

  /**
   * Ensure this coordList is a ring, by adding the start point if necessary
   */
  void closeRing() {
    if (_backingList.length > 0) {
      Coordinate duplicate = _backingList[0].copy();
      addCoord(duplicate, false);
    }
  }

  /** Returns the Coordinates in this collection.
   *
   * @return the coordinates
   */
  List<Coordinate> toCoordinateArray() {
    return _backingList;
  }

  /**
   * Creates an array containing the coordinates in this list,
   * oriented in the given direction (forward or reverse).
   *
   * @param direction the direction value: true for forward, false for reverse
   * @return an oriented array of coordinates
   */
  List<Coordinate> toCoordinateArrayWithCheck(bool isForward) {
    if (isForward) {
      return _backingList;
    }
    // construct reversed array
    List<Coordinate> pts = [];
    _backingList.reversed.forEach((element) {
      pts.add(element);
    });
    return pts;
  }

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
}

/**
 * Represents a planar triangle, and provides methods for calculating various
 * properties of triangles.
 *
 * @version 1.7
 */
class Triangle {
  /**
   * Tests whether a triangle is acute. A triangle is acute iff all interior
   * angles are acute. This is a strict test - right triangles will return
   * <tt>false</tt> A triangle which is not acute is either right or obtuse.
   * <p>
   * Note: this implementation is not robust for angles very close to 90
   * degrees.
   *
   * @param a
   *          a vertex of the triangle
   * @param b
   *          a vertex of the triangle
   * @param c
   *          a vertex of the triangle
   * @return true if the triangle is acute
   */
  static bool isAcuteStatic(Coordinate a, Coordinate b, Coordinate c) {
    if (!Angle.isAcute(a, b, c)) return false;
    if (!Angle.isAcute(b, c, a)) return false;
    if (!Angle.isAcute(c, a, b)) return false;
    return true;
  }

  /**
   * Computes the line which is the perpendicular bisector of the line segment
   * a-b.
   *
   * @param a
   *          a point
   * @param b
   *          another point
   * @return the perpendicular bisector, as an HCoordinate
   */
  static HCoordinate perpendicularBisector(Coordinate a, Coordinate b) {
    // returns the perpendicular bisector of the line segment ab
    double dx = b.x - a.x;
    double dy = b.y - a.y;
    HCoordinate l1 = new HCoordinate.xyw(a.x + dx / 2.0, a.y + dy / 2.0, 1.0);
    HCoordinate l2 =
        new HCoordinate.xyw(a.x - dy + dx / 2.0, a.y + dx + dy / 2.0, 1.0);
    return new HCoordinate.from2Hc(l1, l2);
  }

  /**
   * Computes the circumcentre of a triangle. The circumcentre is the centre of
   * the circumcircle, the smallest circle which encloses the triangle. It is
   * also the common intersection point of the perpendicular bisectors of the
   * sides of the triangle, and is the only point which has equal distance to
   * all three vertices of the triangle.
   *
   * @param a
   *          a vertex of the triangle
   * @param b
   *          a vertex of the triangle
   * @param c
   *          a vertex of the triangle
   * @return the circumcentre of the triangle
   */
  /*
   * // original non-robust algorithm  static Coordinate
   * circumcentre(Coordinate a, Coordinate b, Coordinate c) { // compute the
   * perpendicular bisector of chord ab HCoordinate cab =
   * perpendicularBisector(a, b); // compute the perpendicular bisector of chord
   * bc HCoordinate cbc = perpendicularBisector(b, c); // compute the
   * intersection of the bisectors (circle radii) HCoordinate hcc = new
   * HCoordinate(cab, cbc); Coordinate cc = null; try { cc = new
   * Coordinate(hcc.getX(), hcc.getY()); } catch (NotRepresentableException ex)
   * { // MD - not sure what we can do to prevent this (robustness problem) //
   * Idea - can we condition which edges we choose? throw new
   * IllegalStateException(ex.getMessage()); }
   * 
   * //System.out.println("Acc = " + a.distance(cc) + ", Bcc = " +
   * b.distance(cc) + ", Ccc = " + c.distance(cc) );
   * 
   * return cc; }
   */

  /**
   * Computes the circumcentre of a triangle. The circumcentre is the centre of
   * the circumcircle, the smallest circle which encloses the triangle. It is
   * also the common intersection point of the perpendicular bisectors of the
   * sides of the triangle, and is the only point which has equal distance to
   * all three vertices of the triangle.
   * <p>
   * The circumcentre does not necessarily lie within the triangle. For example,
   * the circumcentre of an obtuse isosceles triangle lies outside the triangle.
   * <p>
   * This method uses an algorithm due to J.R.Shewchuk which uses normalization
   * to the origin to improve the accuracy of computation. (See <i>Lecture Notes
   * on Geometric Robustness</i>, Jonathan Richard Shewchuk, 1999).
   *
   * @param a
   *          a vertex of the triangle
   * @param b
   *          a vertex of the triangle
   * @param c
   *          a vertex of the triangle
   * @return the circumcentre of the triangle
   */
  static Coordinate circumcentreStatic(
      Coordinate a, Coordinate b, Coordinate c) {
    double cx = c.x;
    double cy = c.y;
    double ax = a.x - cx;
    double ay = a.y - cy;
    double bx = b.x - cx;
    double by = b.y - cy;

    double denom = 2 * det(ax, ay, bx, by);
    double numx = det(ay, ax * ax + ay * ay, by, bx * bx + by * by);
    double numy = det(ax, ax * ax + ay * ay, bx, bx * bx + by * by);

    double ccx = cx - numx / denom;
    double ccy = cy + numy / denom;

    return new Coordinate(ccx, ccy);
  }

  /**
   * Computes the circumcentre of a triangle. The circumcentre is the centre of
   * the circumcircle, the smallest circle which encloses the triangle. It is
   * also the common intersection point of the perpendicular bisectors of the
   * sides of the triangle, and is the only point which has equal distance to
   * all three vertices of the triangle.
   * <p>
   * The circumcentre does not necessarily lie within the triangle. For example,
   * the circumcentre of an obtuse isosceles triangle lies outside the triangle.
   * <p>
   * This method uses {@link DD} extended-precision arithmetic to
   * provide more accurate results than {@link #circumcentre(Coordinate, Coordinate, Coordinate)}
   *
   * @param a
   *          a vertex of the triangle
   * @param b
   *          a vertex of the triangle
   * @param c
   *          a vertex of the triangle
   * @return the circumcentre of the triangle
   */
  static Coordinate circumcentreStaticDD(
      Coordinate a, Coordinate b, Coordinate c) {
    DD ax = DD.valueOf(a.x).subtract(c.x);
    DD ay = DD.valueOf(a.y).subtract(c.y);
    DD bx = DD.valueOf(b.x).subtract(c.x);
    DD by = DD.valueOf(b.y).subtract(c.y);

    DD denom = DD.determinantDD(ax, ay, bx, by).multiply(2);
    DD asqr = ax.sqr().addDD(ay.sqr());
    DD bsqr = bx.sqr().addDD(by.sqr());
    DD numx = DD.determinantDD(ay, asqr, by, bsqr);
    DD numy = DD.determinantDD(ax, asqr, bx, bsqr);

    double ccx = DD.valueOf(c.x).subtractDD(numx.divideDD(denom)).doubleValue();
    double ccy = DD.valueOf(c.y).addDD(numy.divideDD(denom)).doubleValue();

    return new Coordinate(ccx, ccy);
  }

  /**
   * Computes the determinant of a 2x2 matrix. Uses standard double-precision
   * arithmetic, so is susceptible to round-off error.
   *
   * @param m00
   *          the [0,0] entry of the matrix
   * @param m01
   *          the [0,1] entry of the matrix
   * @param m10
   *          the [1,0] entry of the matrix
   * @param m11
   *          the [1,1] entry of the matrix
   * @return the determinant
   */
  static double det(double m00, double m01, double m10, double m11) {
    return m00 * m11 - m01 * m10;
  }

  /**
   * Computes the incentre of a triangle. The <i>inCentre</i> of a triangle is
   * the point which is equidistant from the sides of the triangle. It is also
   * the point at which the bisectors of the triangle's angles meet. It is the
   * centre of the triangle's <i>incircle</i>, which is the unique circle that
   * is tangent to each of the triangle's three sides.
   * <p>
   * The incentre always lies within the triangle.
   *
   * @param a
   *          a vertex of the triangle
   * @param b
   *          a vertex of the triangle
   * @param c
   *          a vertex of the triangle
   * @return the point which is the incentre of the triangle
   */
  static Coordinate inCentreStatic(Coordinate a, Coordinate b, Coordinate c) {
    // the lengths of the sides, labelled by their opposite vertex
    double len0 = b.distance(c);
    double len1 = a.distance(c);
    double len2 = a.distance(b);
    double circum = len0 + len1 + len2;

    double inCentreX = (len0 * a.x + len1 * b.x + len2 * c.x) / circum;
    double inCentreY = (len0 * a.y + len1 * b.y + len2 * c.y) / circum;
    return new Coordinate(inCentreX, inCentreY);
  }

  /**
   * Computes the centroid (centre of mass) of a triangle. This is also the
   * point at which the triangle's three medians intersect (a triangle median is
   * the segment from a vertex of the triangle to the midpoint of the opposite
   * side). The centroid divides each median in a ratio of 2:1.
   * <p>
   * The centroid always lies within the triangle.
   *
   *
   * @param a
   *          a vertex of the triangle
   * @param b
   *          a vertex of the triangle
   * @param c
   *          a vertex of the triangle
   * @return the centroid of the triangle
   */
  static Coordinate centroidStatic(Coordinate a, Coordinate b, Coordinate c) {
    double x = (a.x + b.x + c.x) / 3;
    double y = (a.y + b.y + c.y) / 3;
    return new Coordinate(x, y);
  }

  /**
   * Computes the length of the longest side of a triangle
   *
   * @param a
   *          a vertex of the triangle
   * @param b
   *          a vertex of the triangle
   * @param c
   *          a vertex of the triangle
   * @return the length of the longest side of the triangle
   */
  static double longestSideLengthStatic(
      Coordinate a, Coordinate b, Coordinate c) {
    double lenAB = a.distance(b);
    double lenBC = b.distance(c);
    double lenCA = c.distance(a);
    double maxLen = lenAB;
    if (lenBC > maxLen) maxLen = lenBC;
    if (lenCA > maxLen) maxLen = lenCA;
    return maxLen;
  }

  /**
   * Computes the point at which the bisector of the angle ABC cuts the segment
   * AC.
   *
   * @param a
   *          a vertex of the triangle
   * @param b
   *          a vertex of the triangle
   * @param c
   *          a vertex of the triangle
   * @return the angle bisector cut point
   */
  static Coordinate angleBisector(Coordinate a, Coordinate b, Coordinate c) {
    /**
     * Uses the fact that the lengths of the parts of the split segment are
     * proportional to the lengths of the adjacent triangle sides
     */
    double len0 = b.distance(a);
    double len2 = b.distance(c);
    double frac = len0 / (len0 + len2);
    double dx = c.x - a.x;
    double dy = c.y - a.y;

    Coordinate splitPt = new Coordinate(a.x + frac * dx, a.y + frac * dy);
    return splitPt;
  }

  /**
   * Computes the 2D area of a triangle. The area value is always non-negative.
   *
   * @param a
   *          a vertex of the triangle
   * @param b
   *          a vertex of the triangle
   * @param c
   *          a vertex of the triangle
   * @return the area of the triangle
   *
   * @see #signedArea(Coordinate, Coordinate, Coordinate)
   */
  static double areaStatic(Coordinate a, Coordinate b, Coordinate c) {
    return (((c.x - a.x) * (b.y - a.y) - (b.x - a.x) * (c.y - a.y)) / 2).abs();
  }

  /**
   * Computes the signed 2D area of a triangle. The area value is positive if
   * the triangle is oriented CW, and negative if it is oriented CCW.
   * <p>
   * The signed area value can be used to determine point orientation, but the
   * implementation in this method is susceptible to round-off errors. Use
   * {@link Orientation#index(Coordinate, Coordinate, Coordinate)}
   * for robust orientation calculation.
   *
   * @param a
   *          a vertex of the triangle
   * @param b
   *          a vertex of the triangle
   * @param c
   *          a vertex of the triangle
   * @return the signed 2D area of the triangle
   *
   * @see Orientation#index(Coordinate, Coordinate, Coordinate)
   */
  static double signedAreaStatic(Coordinate a, Coordinate b, Coordinate c) {
    /**
     * Uses the formula 1/2 * | u x v | where u,v are the side vectors of the
     * triangle x is the vector cross-product For 2D vectors, this formula
     * simplifies to the expression below
     */
    return ((c.x - a.x) * (b.y - a.y) - (b.x - a.x) * (c.y - a.y)) / 2;
  }

  /**
   * Computes the 3D area of a triangle. The value computed is always
   * non-negative.
   *
   * @param a
   *          a vertex of the triangle
   * @param b
   *          a vertex of the triangle
   * @param c
   *          a vertex of the triangle
   * @return the 3D area of the triangle
   */
  static double area3DStatic(Coordinate a, Coordinate b, Coordinate c) {
    /**
     * Uses the formula 1/2 * | u x v | where u,v are the side vectors of the
     * triangle x is the vector cross-product
     */
    // side vectors u and v
    double ux = b.x - a.x;
    double uy = b.y - a.y;
    double uz = b.getZ() - a.getZ();

    double vx = c.x - a.x;
    double vy = c.y - a.y;
    double vz = c.getZ() - a.getZ();

    // cross-product = u x v
    double crossx = uy * vz - uz * vy;
    double crossy = uz * vx - ux * vz;
    double crossz = ux * vy - uy * vx;

    // tri area = 1/2 * | u x v |
    double absSq = crossx * crossx + crossy * crossy + crossz * crossz;
    double area3D = math.sqrt(absSq) / 2;

    return area3D;
  }

  /**
   * Computes the Z-value (elevation) of an XY point on a three-dimensional
   * plane defined by a triangle whose vertices have Z-values. The defining
   * triangle must not be degenerate (in other words, the triangle must enclose
   * a non-zero area), and must not be parallel to the Z-axis.
   * <p>
   * This method can be used to interpolate the Z-value of a point inside a
   * triangle (for example, of a TIN facet with elevations on the vertices).
   *
   * @param p
   *          the point to compute the Z-value of
   * @param v0
   *          a vertex of a triangle, with a Z ordinate
   * @param v1
   *          a vertex of a triangle, with a Z ordinate
   * @param v2
   *          a vertex of a triangle, with a Z ordinate
   * @return the computed Z-value (elevation) of the point
   */
  static double interpolateZStatic(
      Coordinate p, Coordinate v0, Coordinate v1, Coordinate v2) {
    double x0 = v0.x;
    double y0 = v0.y;
    double a = v1.x - x0;
    double b = v2.x - x0;
    double c = v1.y - y0;
    double d = v2.y - y0;
    double det = a * d - b * c;
    double dx = p.x - x0;
    double dy = p.y - y0;
    double t = (d * dx - b * dy) / det;
    double u = (-c * dx + a * dy) / det;
    double z =
        v0.getZ() + t * (v1.getZ() - v0.getZ()) + u * (v2.getZ() - v0.getZ());
    return z;
  }

  /**
   * The coordinates of the vertices of the triangle
   */
  late Coordinate p0, p1, p2;

  /**
   * Creates a new triangle with the given vertices.
   *
   * @param p0
   *          a vertex
   * @param p1
   *          a vertex
   * @param p2
   *          a vertex
   */
  Triangle(Coordinate p0, Coordinate p1, Coordinate p2) {
    this.p0 = p0;
    this.p1 = p1;
    this.p2 = p2;
  }

  /**
   * Computes the incentre of this triangle. The <i>incentre</i> of a triangle
   * is the point which is equidistant from the sides of the triangle. It is
   * also the point at which the bisectors of the triangle's angles meet. It is
   * the centre of the triangle's <i>incircle</i>, which is the unique circle
   * that is tangent to each of the triangle's three sides.
   *
   * @return the point which is the inCentre of this triangle
   */
  Coordinate inCentre() {
    return inCentreStatic(p0, p1, p2);
  }

  /**
   * Tests whether this triangle is acute. A triangle is acute iff all interior
   * angles are acute. This is a strict test - right triangles will return
   * <tt>false</tt> A triangle which is not acute is either right or obtuse.
   * <p>
   * Note: this implementation is not robust for angles very close to 90
   * degrees.
   *
   * @return true if this triangle is acute
   */
  bool isAcute() {
    return isAcuteStatic(this.p0, this.p1, this.p2);
  }

  /**
   * Computes the circumcentre of this triangle. The circumcentre is the centre
   * of the circumcircle, the smallest circle which encloses the triangle. It is
   * also the common intersection point of the perpendicular bisectors of the
   * sides of the triangle, and is the only point which has equal distance to
   * all three vertices of the triangle.
   * <p>
   * The circumcentre does not necessarily lie within the triangle.
   * <p>
   * This method uses an algorithm due to J.R.Shewchuk which uses normalization
   * to the origin to improve the accuracy of computation. (See <i>Lecture Notes
   * on Geometric Robustness</i>, Jonathan Richard Shewchuk, 1999).
   *
   * @return the circumcentre of this triangle
   */
  Coordinate circumcentre() {
    return circumcentreStatic(this.p0, this.p1, this.p2);
  }

  /**
   * Computes the centroid (centre of mass) of this triangle. This is also the
   * point at which the triangle's three medians intersect (a triangle median is
   * the segment from a vertex of the triangle to the midpoint of the opposite
   * side). The centroid divides each median in a ratio of 2:1.
   * <p>
   * The centroid always lies within the triangle.
   *
   * @return the centroid of this triangle
   */
  Coordinate centroid() {
    return centroidStatic(this.p0, this.p1, this.p2);
  }

  /**
   * Computes the length of the longest side of this triangle
   *
   * @return the length of the longest side of this triangle
   */
  double longestSideLength() {
    return longestSideLengthStatic(this.p0, this.p1, this.p2);
  }

  /**
   * Computes the 2D area of this triangle. The area value is always
   * non-negative.
   *
   * @return the area of this triangle
   *
   * @see #signedArea()
   */
  double area() {
    return areaStatic(this.p0, this.p1, this.p2);
  }

  /**
   * Computes the signed 2D area of this triangle. The area value is positive if
   * the triangle is oriented CW, and negative if it is oriented CCW.
   * <p>
   * The signed area value can be used to determine point orientation, but the
   * implementation in this method is susceptible to round-off errors. Use
   * {@link Orientation#index(Coordinate, Coordinate, Coordinate)}
   * for robust orientation calculation.
   *
   * @return the signed 2D area of this triangle
   *
   * @see Orientation#index(Coordinate, Coordinate, Coordinate)
   */
  double signedArea() {
    return signedAreaStatic(this.p0, this.p1, this.p2);
  }

  /**
   * Computes the 3D area of this triangle. The value computed is always
   * non-negative.
   *
   * @return the 3D area of this triangle
   */
  double area3D() {
    return area3DStatic(this.p0, this.p1, this.p2);
  }

  /**
   * Computes the Z-value (elevation) of an XY point on a three-dimensional
   * plane defined by this triangle (whose vertices must have Z-values). This
   * triangle must not be degenerate (in other words, the triangle must enclose
   * a non-zero area), and must not be parallel to the Z-axis.
   * <p>
   * This method can be used to interpolate the Z-value of a point inside this
   * triangle (for example, of a TIN facet with elevations on the vertices).
   *
   * @param p
   *          the point to compute the Z-value of
   * @return the computed Z-value (elevation) of the point
   */
  double interpolateZ(Coordinate p) {
    if (p == null) throw ArgumentError("Supplied point is null.");
    return interpolateZStatic(p, this.p0, this.p1, this.p2);
  }
}
