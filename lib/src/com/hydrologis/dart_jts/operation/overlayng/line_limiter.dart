part of dart_jts;

/**
 * Limits the segments in a list of segments
 * to those which intersect an envelope.
 * This creates zero or more sections of the input segment sequences,
 * containing only line segments which intersect the limit envelope.
 * Segments are not clipped, since that can move
 * line segments enough to alter topology,
 * and it happens in the overlay in any case.
 * This can substantially reduce the number of vertices which need to be
 * processed during overlay.
 * <p>
 * This optimization is only applicable to Line geometries,
 * since it does not maintain the closed topology of rings.
 * Polygonal geometries are optimized using the {@link RingClipper}.
 *
 * @author Martin Davis
 *
 * @see RingClipper
 */
class LineLimiter {
  Envelope limitEnv;
  CoordinateList? ptList;
  Coordinate? lastOutside = null;
  List<List<Coordinate>>? sections = null;

  /**
   * Creates a new limiter for a given envelope.
   *
   * @param env the envelope to limit to
   */
  LineLimiter(this.limitEnv);

  /**
   * Limits a list of segments.
   *
   * @param pts the segment sequence to limit
   * @return the sections which intersect the limit envelope
   */
  List<List<Coordinate>> limit(List<Coordinate> pts) {
    lastOutside = null;
    ptList = null;
    sections = List.empty(growable: true);

    for (int i = 0; i < pts.length; i++) {
      Coordinate p = pts[i];
      if (limitEnv.intersects(p.x, p.y))
        addPoint(p);
      else {
        addOutside(p);
      }
    }
    // finish last section, if any
    _finishSection();
    return sections!;
  }

  void addPoint(Coordinate? p) {
    if (p == null) return;
    _startSection();
    ptList!.addCoord(p, false);
  }

  void addOutside(Coordinate p) {
    bool segIntersects = _isLastSegmentIntersecting(p);
    if (!segIntersects) {
      _finishSection();
    } else {
      addPoint(lastOutside);
      addPoint(p);
    }
    lastOutside = p;
  }

  bool _isLastSegmentIntersecting(Coordinate p) {
    if (lastOutside == null) {
      // last point must have been inside
      if (_isSectionOpen()) return true;
      return false;
    }
    return limitEnv.intersectsEnvelopeCoordinates(lastOutside!, p);
  }

  bool _isSectionOpen() {
    return ptList != null;
  }

  void _startSection() {
    if (ptList == null) {
      ptList = CoordinateList();
    }
    if (lastOutside != null) {
      ptList!.addCoord(lastOutside!, false);
    }
    lastOutside = null;
  }

  void _finishSection() {
    if (ptList == null) return;
    // finish off this section
    if (lastOutside != null) {
      ptList!.addCoord(lastOutside!, false);
      lastOutside = null;
    }

    List<Coordinate> section = ptList!.toCoordinateArray(true);
    sections!.add(section);
    ptList = null;
  }
}
