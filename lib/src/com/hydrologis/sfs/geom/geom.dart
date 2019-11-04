part of dart_sfs;

/**
 * Indicates an invalid or inconsistent topological situation encountered during processing
 *
 * @version 1.7
 */
class TopologyException implements Exception {
  static String msgWithCoord(String msg, Coordinate pt) {
    if (pt != null) {
      return msg + " [ $pt ]";
    }
    return msg;
  }

  Coordinate pt;

  String msg;

  TopologyException(this.msg);

  TopologyException.withCoord(String msg, Coordinate pt) {
    this.msg = msgWithCoord(msg, pt);
    this.pt = Coordinate.fromCoordinate(pt);
  }

  Coordinate getCoordinate() {
    return pt;
  }
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
  Geometry parent;

  ///  Indicates whether or not the first element
  ///  (the root <code>GeometryCollection</code>) has been returned.
  bool atStart;

  ///  The number of <code>Geometry</code>s in the the <code>GeometryCollection</code>.
  int max;

  ///  The index of the <code>Geometry</code> that will be returned when <code>next</code>
  ///  is called.
  int index;

  ///  The iterator over a nested <code>Geometry</code>, or <code>null</code>
  ///  if this <code>GeometryCollectionIterator</code> is not currently iterating
  ///  over a nested <code>GeometryCollection</code>.
  GeometryCollectionIterator subcollectionIterator;

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
      if (subcollectionIterator.hasNext()) {
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
      if (subcollectionIterator.hasNext()) {
        return subcollectionIterator.next();
      } else {
        subcollectionIterator = null;
      }
    }
    if (index >= max) {
      throw RangeError("index >= max");
    }
    Geometry obj = parent.getGeometryN(index++);
    if (obj is GeometryCollection) {
      subcollectionIterator = new GeometryCollectionIterator(obj);
      // there will always be at least one element in the sub-collection
      return subcollectionIterator.next();
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
    throw new UnsupportedError(this.runtimeType.toString());
  }
}
