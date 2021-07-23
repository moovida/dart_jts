part of dart_jts;

/**
 * Models a collection of {@link Geometry}s of
 * arbitrary type and dimension.
 *
 *
 *@version 1.7
 */
class GeometryCollection extends Geometry {
//  With contributions from Markus Schaber [schabios@logi-track.com] 2004-03-26
  /**
   *  Internal representation of this <code>GeometryCollection</code>.
   */
  late List<Geometry> geometries;

  /** @deprecated Use GeometryFactory instead */
  GeometryCollection(
      List<Geometry> geometries, PrecisionModel precisionModel, int SRID)
      : this.withFactory(geometries,
            new GeometryFactory.withPrecisionModelSrid(precisionModel, SRID));

  /**
   * @param geometries
   *            the <code>Geometry</code>s for this <code>GeometryCollection</code>,
   *            or <code>null</code> or an empty array to create the empty
   *            geometry. Elements may be empty <code>Geometry</code>s,
   *            but not <code>null</code>s.
   */
  GeometryCollection.withFactory(
      List<Geometry>? geometries, GeometryFactory factory)
      : super(factory) {
    if (geometries == null) {
      geometries = [];
    }
    if (Geometry.hasNullElements(geometries)) {
      throw new ArgumentError("geometries must not contain null elements");
    }
    this.geometries = geometries;
  }

  Coordinate? getCoordinate() {
    if (isEmpty()) return null;
    return geometries[0].getCoordinate();
  }

  /**
   * Collects all coordinates of all subgeometries into an Array.
   *
   * Note that while changes to the coordinate objects themselves
   * may modify the Geometries in place, the returned Array as such
   * is only a temporary container which is not synchronized back.
   *
   * @return the collected coordinates
   *    */
  List<Coordinate> getCoordinates() {
    List<Coordinate> coordinates = []; //..length = getNumPoints();
    // int k = -1;
    for (int i = 0; i < geometries.length; i++) {
      List<Coordinate> childCoordinates = geometries[i].getCoordinates();
      for (int j = 0; j < childCoordinates.length; j++) {
        // k++;
        // coordinates[k] = childCoordinates[j];
        coordinates.add(childCoordinates[j]);
      }
    }
    return coordinates;
  }

  bool isEmpty() {
    for (int i = 0; i < geometries.length; i++) {
      if (!geometries[i].isEmpty()) {
        return false;
      }
    }
    return true;
  }

  int getDimension() {
    int dimension = Dimension.FALSE;
    for (int i = 0; i < geometries.length; i++) {
      dimension = math.max(dimension, geometries[i].getDimension());
    }
    return dimension;
  }

  int getBoundaryDimension() {
    int dimension = Dimension.FALSE;
    for (int i = 0; i < geometries.length; i++) {
      dimension = math.max(
          dimension, (geometries[i] as Geometry).getBoundaryDimension());
    }
    return dimension;
  }

  int getNumGeometries() {
    return geometries.length;
  }

  Geometry getGeometryN(int n) {
    return geometries[n];
  }

  int getNumPoints() {
    int numPoints = 0;
    for (int i = 0; i < geometries.length; i++) {
      numPoints += (geometries[i] as Geometry).getNumPoints();
    }
    return numPoints;
  }

  String getGeometryType() {
    return "GeometryCollection";
  }

  Geometry getBoundary() {
    Geometry.checkNotGeometryCollection(this);
    Assert.shouldNeverReachHere();
    throw StateError("Should never reach here");
  }

  /**
   *  Returns the area of this <code>GeometryCollection</code>
   *
   * @return the area of the polygon
   */
  double getArea() {
    double area = 0.0;
    for (int i = 0; i < geometries.length; i++) {
      area += geometries[i].getArea();
    }
    return area;
  }

  double getLength() {
    double sum = 0.0;
    for (int i = 0; i < geometries.length; i++) {
      sum += (geometries[i]).getLength();
    }
    return sum;
  }

  bool equalsExactWithTol(Geometry other, double tolerance) {
    if (!isEquivalentClass(other)) {
      return false;
    }
    GeometryCollection otherCollection = other as GeometryCollection;
    if (geometries.length != otherCollection.geometries.length) {
      return false;
    }
    for (int i = 0; i < geometries.length; i++) {
      if (!(geometries[i] as Geometry)
          .equalsExactWithTol(otherCollection.geometries[i], tolerance)) {
        return false;
      }
    }
    return true;
  }

  void applyCF(CoordinateFilter filter) {
    for (int i = 0; i < geometries.length; i++) {
      geometries[i].applyCF(filter);
    }
  }

  void applyCSF(CoordinateSequenceFilter filter) {
    if (geometries.length == 0) return;
    for (int i = 0; i < geometries.length; i++) {
      geometries[i].applyCSF(filter);
      if (filter.isDone()) {
        break;
      }
    }
    if (filter.isGeometryChanged()) geometryChanged();
  }

  void applyGF(GeometryFilter filter) {
    filter.filter(this);
    for (int i = 0; i < geometries.length; i++) {
      geometries[i].applyGF(filter);
    }
  }

  void applyGCF(GeometryComponentFilter filter) {
    filter.filter(this);
    for (int i = 0; i < geometries.length; i++) {
      geometries[i].applyGCF(filter);
    }
  }

  /**
   * Creates and returns a full copy of this {@link GeometryCollection} object.
   * (including all coordinates contained by it).
   *
   * @return a clone of this instance
   * @deprecated
   */
  Object clone() {
    return copy();
  }

  GeometryCollection copyInternal() {
    List<Geometry> geometries = []; //..length = (this.geometries.length);
    for (int i = 0; i < this.geometries.length; i++) {
      // geometries[i] = this.geometries[i].copy();
      geometries.add(this.geometries[i].copy());
    }
    return new GeometryCollection.withFactory(geometries, geomFactory);
  }

  void normalize() {
    for (int i = 0; i < geometries.length; i++) {
      geometries[i].normalize();
    }
    geometries.sort();
  }

  Envelope computeEnvelopeInternal() {
    Envelope envelope = new Envelope.empty();
    for (int i = 0; i < geometries.length; i++) {
      envelope.expandToIncludeEnvelope(geometries[i].getEnvelopeInternal());
    }
    return envelope;
  }

  int compareToSameClass(Object o) {
    Set theseElements = SplayTreeSet.from(geometries);
    Set otherElements = SplayTreeSet.from((o as GeometryCollection).geometries);
    return compare(theseElements.toList(), otherElements.toList());
  }

  int compareToSameClassWithComparator(
      Object o, Comparator<CoordinateSequence> comp) {
    GeometryCollection gc = o as GeometryCollection;

    int n1 = getNumGeometries();
    int n2 = gc.getNumGeometries();
    int i = 0;
    while (i < n1 && i < n2) {
      Geometry thisGeom = getGeometryN(i);
      Geometry otherGeom = gc.getGeometryN(i);
      int holeComp = thisGeom.compareToSameClassWithComparator(otherGeom, comp);
      if (holeComp != 0) return holeComp;
      i++;
    }
    if (i < n1) return 1;
    if (i < n2) return -1;
    return 0;
  }

  int getSortIndex() {
    return Geometry.SORTINDEX_GEOMETRYCOLLECTION;
  }

  /**
   * Creates a {@link GeometryCollection} with
   * every component reversed.
   * The order of the components in the collection are not reversed.
   *
   * @return a {@link GeometryCollection} in the reverse order
   */
  Geometry reverse() {
    int n = geometries.length;
    List<Geometry> revGeoms = []; //..length = n;
    for (int i = 0; i < n; i++) {
      revGeoms.add(geometries[i].reverse());
      // revGeoms[i] = geometries[i].reverse();
    }
    return getFactory().createGeometryCollection(revGeoms);
  }
}
