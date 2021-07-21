part of dart_jts;

/**
 * Represents a single point.
 *
 * A <code>Point</code> is topologically valid if and only if:
 * <ul>
 * <li>the coordinate which defines it (if any) is a valid coordinate
 * (i.e. does not have an <code>NaN</code> X or Y ordinate)
 * </ul>
 *
 *@version 1.7
 */
class Point extends Geometry implements Puntal {
  /**
   *  The <code>Coordinate</code> wrapped by this <code>Point</code>.
   */
  CoordinateSequence? coordinates;

  /**
   *  Constructs a <code>Point</code> with the given coordinate.
   *
   *@param  coordinate      the coordinate on which to base this <code>Point</code>
   *      , or <code>null</code> to create the empty geometry.
   *@param  precisionModel  the specification of the grid of allowable points
   *      for this <code>Point</code>
   *@param  SRID            the ID of the Spatial Reference System used by this
   *      <code>Point</code>
   * @deprecated Use GeometryFactory instead
   */
  Point(Coordinate? coordinate, PrecisionModel precisionModel, int SRID)
      : super(
            new GeometryFactory.withPrecisionModelSrid(precisionModel, SRID)) {
    init(getFactory()
        .getCoordinateSequenceFactory()
        .create(coordinate != null ? [coordinate] : []));
  }

  /**
   *@param  coordinates      contains the single coordinate on which to base this <code>Point</code>
   *      , or <code>null</code> to create the empty geometry.
   */
  Point.fromSequence(CoordinateSequence? coordinates, GeometryFactory factory)
      : super(factory) {
    init(coordinates);
  }

  void init(CoordinateSequence? coordinates) {
    if (coordinates == null) {
      coordinates = getFactory().getCoordinateSequenceFactory().create([]);
    }
    Assert.isTrue(coordinates.size() <= 1);
    this.coordinates = coordinates;
  }

  List<Coordinate> getCoordinates() {
    return isEmpty() ? [] : [getCoordinate()!];
  }

  int getNumPoints() {
    return isEmpty() ? 0 : 1;
  }

  bool isEmpty() {
    return coordinates!.size() == 0;
  }

  bool isSimple() {
    return true;
  }

  int getDimension() {
    return 0;
  }

  int getBoundaryDimension() {
    return Dimension.FALSE;
  }

  double getX() {
    if (getCoordinate() == null) {
      throw new StateError("getX called on empty Point");
    }
    return getCoordinate()!.x;
  }

  double getY() {
    if (getCoordinate() == null) {
      throw new StateError("getY called on empty Point");
    }
    return getCoordinate()!.y;
  }

  Coordinate? getCoordinate() {
    return coordinates!.size() != 0 ? coordinates!.getCoordinate(0) : null;
  }

  String getGeometryType() {
    return "Point";
  }

  /**
   * Gets the boundary of this geometry.
   * Zero-dimensional geometries have no boundary by definition,
   * so an empty GeometryCollection is returned.
   *
   * @return an empty GeometryCollection
   * @see Geometry#getBoundary
   */
  Geometry getBoundary() {
    return getFactory().createGeometryCollectionEmpty();
  }

  Envelope computeEnvelopeInternal() {
    if (isEmpty()) {
      return new Envelope.empty();
    }
    Envelope env = new Envelope.empty();
    env.expandToInclude(coordinates!.getX(0), coordinates!.getY(0));
    return env;
  }

  bool equalsExactWithTol(Geometry other, double tolerance) {
    if (!isEquivalentClass(other)) {
      return false;
    }
    if (isEmpty() && other.isEmpty()) {
      return true;
    }
    if (isEmpty() != other.isEmpty()) {
      return false;
    }
    return equal(
        (other as Point).getCoordinate()!, this.getCoordinate()!, tolerance);
  }

  void applyCF(CoordinateFilter cf) {
    if (isEmpty()) {
      return;
    }
    cf.filter(getCoordinate());
  }

  void applyCSF(CoordinateSequenceFilter filter) {
    if (isEmpty()) return;
    filter.filter(coordinates!, 0);
    if (filter.isGeometryChanged()) geometryChanged();
  }

  void applyGF(GeometryFilter filter) {
    filter.filter(this);
  }

  void applyGCF(GeometryComponentFilter filter) {
    filter.filter(this);
  }

  /**
   * Creates and returns a full copy of this {@link Point} object.
   * (including all coordinates contained by it).
   *
   * @return a clone of this instance
   * @deprecated
   */
  Object clone() {
    return copy();
  }

  Point copyInternal() {
    return new Point.fromSequence(coordinates!.copy(), geomFactory);
  }

  Geometry reverse() {
    return copy();
  }

  void normalize() {
    // a Point is always in normalized form
  }

  int compareToSameClass(Object other) {
    Point point = other as Point;
    return getCoordinate()!.compareTo(point.getCoordinate()!);
  }

  int compareToSameClassWithComparator(
      Object other, Comparator<CoordinateSequence> comp) {
    Point point = other as Point;
    return comp(this.coordinates!, point.coordinates!);
  }

  int getSortIndex() {
    return Geometry.SORTINDEX_POINT;
  }

  CoordinateSequence getCoordinateSequence() {
    return coordinates!;
  }
}
