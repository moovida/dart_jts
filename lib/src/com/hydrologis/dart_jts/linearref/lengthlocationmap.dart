import '../geom/coordinate.dart';
import '../geom/geometry.dart';
import 'linearlocation.dart';
import 'lineariterator.dart';

/**
 * Computes the {@link LinearLocation} for a given length
 * along a linear {@link Geometry}.
 * Negative lengths are measured in reverse from end of the linear geometry.
 * Out-of-range values are clamped.
 */
class LengthLocationMap {
  // TODO: cache computed cumulative length for each vertex
  // TODO: support user-defined measures
  // TODO: support measure index for fast mapping to a location

  /**
   * Computes the {@link LinearLocation} for a
   * given length along a linear {@link Geometry}.
   *
   * @param linearGeom the linear geometry to use
   * @param length the length index of the location
   * @return the {@link LinearLocation} for the length
   */
  static LinearLocation getLocationStatic(Geometry linearGeom, double length) {
    LengthLocationMap locater = new LengthLocationMap(linearGeom);
    return locater.getLocation(length);
  }

  /**
   * Computes the {@link LinearLocation} for a
   * given length along a linear {@link Geometry},
   * with control over how the location
   * is resolved at component endpoints.
   *
   * @param linearGeom the linear geometry to use
   * @param length the length index of the location
   * @param resolveLower if true lengths are resolved to the lowest possible index
   * @return the {@link LinearLocation} for the length
   */
  static LinearLocation getLocationStaticWithResolver(Geometry linearGeom, double length, bool resolveLower) {
    LengthLocationMap locater = new LengthLocationMap(linearGeom);
    return locater.getLocationWithResolver(length, resolveLower);
  }

  /**
   * Computes the length for a given {@link LinearLocation}
   * on a linear {@link Geometry}.
   *
   * @param linearGeom the linear geometry to use
   * @param loc the {@link LinearLocation} index of the location
   * @return the length for the {@link LinearLocation}
   */
  static double getLengthStatic(Geometry linearGeom, LinearLocation loc) {
    LengthLocationMap locater = new LengthLocationMap(linearGeom);
    return locater.getLength(loc);
  }

  Geometry linearGeom;

  LengthLocationMap(this.linearGeom);

  /**
   * Compute the {@link LinearLocation} corresponding to a length.
   * Negative lengths are measured in reverse from end of the linear geometry.
   * Out-of-range values are clamped.
   * Ambiguous indexes are resolved to the lowest possible location value.
   *
   * @param length the length index
   * @return the corresponding LinearLocation
   */
  LinearLocation getLocation(double length) {
    return getLocationWithResolver(length, true);
  }

  /**
   * Compute the {@link LinearLocation} corresponding to a length.
   * Negative lengths are measured in reverse from end of the linear geometry.
   * Out-of-range values are clamped.
   * Ambiguous indexes are resolved to the lowest or highest possible location value,
   * depending on the value of <tt>resolveLower</tt>
   *
   * @param length the length index
   * @return the corresponding LinearLocation
   */
  LinearLocation getLocationWithResolver(double length, bool resolveLower) {
    double forwardLength = length;

    // negative values are measured from end of geometry
    if (length < 0.0) {
      double lineLen = linearGeom.getLength();
      forwardLength = lineLen + length;
    }
    LinearLocation loc = getLocationForward(forwardLength);
    if (resolveLower) {
      return loc;
    }
    return resolveHigher(loc);
  }

  LinearLocation getLocationForward(double length) {
    if (length <= 0.0) return new LinearLocation();

    double totalLength = 0.0;

    LinearIterator it = new LinearIterator(linearGeom);
    while (it.hasNext()) {
      /**
       * Special handling is required for the situation when the 
       * length references exactly to a component endpoint.
       * In this case, the endpoint location of the current component 
       * is returned,
       * rather than the startpoint location of the next component.
       * This produces consistent behaviour with the project method.
       */
      if (it.isEndOfLine()) {
        if (totalLength == length) {
          int compIndex = it.getComponentIndex();
          int segIndex = it.getVertexIndex();
          return new LinearLocation.fromComponentSegmentIndexFraction(compIndex, segIndex, 0.0);
        }
      } else {
        Coordinate p0 = it.getSegmentStart();
        Coordinate p1 = it.getSegmentEnd()!;
        double segLen = p1.distance(p0);
        // length falls in this segment
        if (totalLength + segLen > length) {
          double frac = (length - totalLength) / segLen;
          int compIndex = it.getComponentIndex();
          int segIndex = it.getVertexIndex();
          return new LinearLocation.fromComponentSegmentIndexFraction(compIndex, segIndex, frac);
        }
        totalLength += segLen;
      }

      it.next();
    }
    // length is longer than line - return end location
    return LinearLocation.getEndLocation(linearGeom);
  }

  LinearLocation resolveHigher(LinearLocation loc) {
    if (!loc.isEndpoint(linearGeom)) return loc;
    int compIndex = loc.getComponentIndex();
    // if last component can't resolve any higher
    if (compIndex >= linearGeom.getNumGeometries() - 1) return loc;

    do {
      compIndex++;
    } while (compIndex < linearGeom.getNumGeometries() - 1 && linearGeom.getGeometryN(compIndex).getLength() == 0);
    // resolve to next higher location
    return new LinearLocation.fromComponentSegmentIndexFraction(compIndex, 0, 0.0);
  }

  double getLength(LinearLocation loc) {
    double totalLength = 0.0;

    LinearIterator it = new LinearIterator(linearGeom);
    while (it.hasNext()) {
      if (!it.isEndOfLine()) {
        Coordinate p0 = it.getSegmentStart();
        Coordinate p1 = it.getSegmentEnd()!;
        double segLen = p1.distance(p0);
        // length falls in this segment
        if (loc.getComponentIndex() == it.getComponentIndex() && loc.getSegmentIndex() == it.getVertexIndex()) {
          return totalLength + segLen * loc.getSegmentFraction();
        }
        totalLength += segLen;
      }
      it.next();
    }
    return totalLength;
  }
}
