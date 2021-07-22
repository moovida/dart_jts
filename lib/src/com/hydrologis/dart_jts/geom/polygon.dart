part of dart_jts;

/**
 * Represents a polygon with linear edges, which may include holes.
 * The outer boundary (shell)
 * and inner boundaries (holes) of the polygon are represented by {@link LinearRing}s.
 * The boundary rings of the polygon may have any orientation.
 * Polygons are closed, simple geometries by definition.
 * <p>
 * The polygon model conforms to the assertions specified in the
 * <A HREF="http://www.opengis.org/techno/specs.htm">OpenGIS Simple Features
 * Specification for SQL</A>.
 * <p>
 * A <code>Polygon</code> is topologically valid if and only if:
 * <ul>
 * <li>the coordinates which define it are valid coordinates
 * <li>the linear rings for the shell and holes are valid
 * (i.e. are closed and do not self-intersect)
 * <li>holes touch the shell or another hole at at most one point
 * (which implies that the rings of the shell and holes must not cross)
 * <li>the interior of the polygon is connected,
 * or equivalently no sequence of touching holes
 * makes the interior of the polygon disconnected
 * (i.e. effectively split the polygon into two pieces).
 * </ul>
 *
 *@version 1.7
 */
class Polygon extends Geometry implements Polygonal {
  /**
   *  The exterior boundary,
   * or <code>null</code> if this <code>Polygon</code>
   *  is empty.
   */
  LinearRing? shell = null;

  /**
   * The interior boundaries, if any.
   * This instance var is never null.
   * If there are no holes, the array is of zero length.
   */
  List<LinearRing>? holes = null;

  /**
   *  Constructs a <code>Polygon</code> with the given exterior boundary.
   *
   *@param  shell           the outer boundary of the new <code>Polygon</code>,
   *      or <code>null</code> or an empty <code>LinearRing</code> if the empty
   *      geometry is to be created.
   *@param  precisionModel  the specification of the grid of allowable points
   *      for this <code>Polygon</code>
   *@param  SRID            the ID of the Spatial Reference System used by this
   *      <code>Polygon</code>
   * @deprecated Use GeometryFactory instead
   */
  Polygon(LinearRing shell, PrecisionModel precisionModel, int SRID)
      : this.withFactory(shell, <LinearRing>[],
            new GeometryFactory.withPrecisionModelSrid(precisionModel, SRID));

  /**
   *  Constructs a <code>Polygon</code> with the given exterior boundary and
   *  interior boundaries.
   *
   *@param  shell           the outer boundary of the new <code>Polygon</code>,
   *      or <code>null</code> or an empty <code>LinearRing</code> if the empty
   *      geometry is to be created.
   *@param  holes           the inner boundaries of the new <code>Polygon</code>
   *      , or <code>null</code> or empty <code>LinearRing</code>s if the empty
   *      geometry is to be created.
   *@param  precisionModel  the specification of the grid of allowable points
   *      for this <code>Polygon</code>
   *@param  SRID            the ID of the Spatial Reference System used by this
   *      <code>Polygon</code>
   * @deprecated Use GeometryFactory instead
   */
  Polygon.withPrecisionModelSrid(LinearRing shell, List<LinearRing> holes,
      PrecisionModel precisionModel, int SRID)
      : this.withFactory(shell, holes,
            new GeometryFactory.withPrecisionModelSrid(precisionModel, SRID));

  /**
   *  Constructs a <code>Polygon</code> with the given exterior boundary and
   *  interior boundaries.
   *
   *@param  shell           the outer boundary of the new <code>Polygon</code>,
   *      or <code>null</code> or an empty <code>LinearRing</code> if the empty
   *      geometry is to be created.
   *@param  holes           the inner boundaries of the new <code>Polygon</code>
   *      , or <code>null</code> or empty <code>LinearRing</code>s if the empty
   *      geometry is to be created.
   */
  Polygon.withFactory(
      LinearRing? shell, List<LinearRing>? holes, GeometryFactory factory)
      : super(factory) {
    if (shell == null) {
      shell = getFactory().createLinearRingEmpty();
    }
    if (holes == null) {
      holes = <LinearRing>[];
    }
    if (Geometry.hasNullElements(holes)) {
      throw new ArgumentError("holes must not contain null elements");
    }
    if (shell.isEmpty() && Geometry.hasNonEmptyElements(holes)) {
      throw new ArgumentError("shell is empty but holes are not");
    }
    this.shell = shell;
    this.holes = holes;
  }

  Coordinate? getCoordinate() {
    return shell!.getCoordinate();
  }

  List<Coordinate> getCoordinates() {
    if (isEmpty()) {
      return <Coordinate>[];
    }
    List<Coordinate> coordinates = []; //..length = getNumPoints();
    // int k = -1;
    List<Coordinate> shellCoordinates = shell!.getCoordinates();
    for (int x = 0; x < shellCoordinates.length; x++) {
      // k++;
      // coordinates[k] = shellCoordinates[x];
      coordinates.add(shellCoordinates[x]);
    }
    for (int i = 0; i < holes!.length; i++) {
      List<Coordinate> childCoordinates = holes![i].getCoordinates();
      for (int j = 0; j < childCoordinates.length; j++) {
        // k++;
        // coordinates[k] = childCoordinates[j];
        coordinates.add(childCoordinates[j]);
      }
    }
    return coordinates;
  }

  int getNumPoints() {
    int numPoints = shell!.getNumPoints();
    for (int i = 0; i < holes!.length; i++) {
      numPoints += holes![i].getNumPoints();
    }
    return numPoints;
  }

  int getDimension() {
    return 2;
  }

  int getBoundaryDimension() {
    return 1;
  }

  bool isEmpty() {
    return shell!.isEmpty();
  }

  bool isRectangle() {
    if (getNumInteriorRing() != 0) return false;
    if (shell == null) return false;
    if (shell!.getNumPoints() != 5) return false;

    CoordinateSequence seq = shell!.getCoordinateSequence();

    // check vertices have correct values
    Envelope env = getEnvelopeInternal();
    for (int i = 0; i < 5; i++) {
      double x = seq.getX(i);
      if (!(x == env.getMinX() || x == env.getMaxX())) return false;
      double y = seq.getY(i);
      if (!(y == env.getMinY() || y == env.getMaxY())) return false;
    }

    // check vertices are in right order
    double prevX = seq.getX(0);
    double prevY = seq.getY(0);
    for (int i = 1; i <= 4; i++) {
      double x = seq.getX(i);
      double y = seq.getY(i);
      bool xChanged = x != prevX;
      bool yChanged = y != prevY;
      if (xChanged == yChanged) return false;
      prevX = x;
      prevY = y;
    }
    return true;
  }

  LinearRing getExteriorRing() {
    return shell!;
  }

  int getNumInteriorRing() {
    return holes!.length;
  }

  LinearRing getInteriorRingN(int n) {
    return holes![n];
  }

  String getGeometryType() {
    return "Polygon";
  }

  /**
   *  Returns the area of this <code>Polygon</code>
   *
   *@return the area of the polygon
   */
  double getArea() {
    double area = 0.0;
    area += Area.ofRingSeq(shell!.getCoordinateSequence());
    for (int i = 0; i < holes!.length; i++) {
      area -= Area.ofRingSeq(holes![i].getCoordinateSequence());
    }
    return area;
  }

  /**
   *  Returns the perimeter of this <code>Polygon</code>
   *
   *@return the perimeter of the polygon
   */
  double getLength() {
    double len = 0.0;
    len += shell!.getLength();
    for (int i = 0; i < holes!.length; i++) {
      len += holes![i].getLength();
    }
    return len;
  }

  /**
   * Computes the boundary of this geometry
   *
   * @return a lineal geometry (which may be empty)
   * @see Geometry#getBoundary
   */
  Geometry getBoundary() {
    if (isEmpty()) {
      return getFactory().createMultiLineStringEmpty();
    }
    List<LinearRing> rings = []; //..length = holes!.length + 1;
    rings.add(shell!);
    // rings[0] = shell!;
    for (int i = 0; i < holes!.length; i++) {
      rings.add(holes![i]);
      // rings[i + 1] = holes![i];
    }
    // create LineString or MultiLineString as appropriate
    if (rings.length <= 1)
      return getFactory().createLinearRingSeq(rings[0].getCoordinateSequence());
    return getFactory().createMultiLineString(rings);
  }

  Envelope computeEnvelopeInternal() {
    return shell!.getEnvelopeInternal();
  }

  bool equalsExactWithTol(Geometry other, double tolerance) {
    if (!isEquivalentClass(other)) {
      return false;
    }
    Polygon otherPolygon = other as Polygon;
    Geometry thisShell = shell!;
    Geometry otherPolygonShell = otherPolygon.shell!;
    if (!thisShell.equalsExactWithTol(otherPolygonShell, tolerance)) {
      return false;
    }
    if (holes!.length != otherPolygon.holes!.length) {
      return false;
    }
    for (int i = 0; i < holes!.length; i++) {
      if (!holes![i].equalsExactWithTol(otherPolygon.holes![i], tolerance)) {
        return false;
      }
    }
    return true;
  }

  void applyCF(CoordinateFilter filter) {
    shell!.applyCF(filter);
    for (int i = 0; i < holes!.length; i++) {
      holes![i].applyCF(filter);
    }
  }

  void applyCSF(CoordinateSequenceFilter filter) {
    shell!.applyCSF(filter);
    if (!filter.isDone()) {
      for (int i = 0; i < holes!.length; i++) {
        holes![i].applyCSF(filter);
        if (filter.isDone()) break;
      }
    }
    if (filter.isGeometryChanged()) geometryChanged();
  }

  void applyGF(GeometryFilter filter) {
    filter.filter(this);
  }

  void applyGCF(GeometryComponentFilter filter) {
    filter.filter(this);
    shell!.applyGCF(filter);
    for (int i = 0; i < holes!.length; i++) {
      holes![i].applyGCF(filter);
    }
  }

  /**
   * Creates and returns a full copy of this {@link Polygon} object.
   * (including all coordinates contained by it).
   *
   * @return a clone of this instance
   * @deprecated
   */
  Object clone() {
    return copy();
  }

  Polygon copyInternal() {
    LinearRing shellCopy = shell!.copy() as LinearRing;
    List<LinearRing> holeCopies = []; //..length = this.holes!.length;
    for (int i = 0; i < holes!.length; i++) {
      holeCopies.add(holes![i].copy() as LinearRing);
      // holeCopies[i] = holes![i].copy() as LinearRing;
    }
    return new Polygon.withFactory(shellCopy, holeCopies, geomFactory);
  }

  Geometry convexHull() {
    return getExteriorRing().convexHull();
  }

  void normalize() {
    shell = normalized(shell!, true);
    for (int i = 0; i < holes!.length; i++) {
      holes![i] = normalized(holes![i], false);
    }
    holes!.sort();
  }

  int compareToSameClass(Object o) {
    LinearRing thisShell = shell!;
    LinearRing otherShell = (o as Polygon).shell!;
    return thisShell.compareToSameClass(otherShell);
  }

  int compareToSameClassWithComparator(
      Object o, Comparator<CoordinateSequence> comp) {
    Polygon poly = o as Polygon;

    LinearRing thisShell = shell!;
    LinearRing otherShell = poly.shell!;
    int shellComp =
        thisShell.compareToSameClassWithComparator(otherShell, comp);
    if (shellComp != 0) return shellComp;

    int nHole1 = getNumInteriorRing();
    int nHole2 = poly.getNumInteriorRing();
    int i = 0;
    while (i < nHole1 && i < nHole2) {
      LinearRing thisHole = getInteriorRingN(i);
      LinearRing otherHole = poly.getInteriorRingN(i);
      int holeComp = thisHole.compareToSameClassWithComparator(otherHole, comp);
      if (holeComp != 0) return holeComp;
      i++;
    }
    if (i < nHole1) return 1;
    if (i < nHole2) return -1;
    return 0;
  }

  int getSortIndex() {
    return Geometry.SORTINDEX_POLYGON;
  }

  LinearRing normalized(LinearRing ring, bool clockwise) {
    LinearRing res = ring.copy() as LinearRing;
    normalizeRing(res, clockwise);
    return res;
  }

  void normalizeRing(LinearRing ring, bool clockwise) {
    if (ring.isEmpty()) {
      return;
    }

    CoordinateSequence seq = ring.getCoordinateSequence();
    int minCoordinateIndex =
        CoordinateSequences.minCoordinateIndexWithRange(seq, 0, seq.size() - 2);
    CoordinateSequences.scrollWithIndexAndRingcheck(
        seq, minCoordinateIndex, true);
    if (Orientation.isCCWFromSeq(seq) == clockwise)
      CoordinateSequences.reverse(seq);
  }

  Geometry reverse() {
    Polygon poly = copy() as Polygon;
    poly.shell = shell!.copy().reverse() as LinearRing;
    poly.holes = []; //..length = holes!.length;
    for (int i = 0; i < holes!.length; i++) {
      poly.holes!.add(holes![i].copy().reverse() as LinearRing);
      // poly.holes![i] = holes![i].copy().reverse() as LinearRing;
    }
    return poly; // return the clone
  }
}
