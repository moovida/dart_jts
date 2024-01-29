part of dart_jts;

/**
 * Records topological information about an
 * edge representing a piece of linework (lineString or polygon ring)
 * from a single source geometry.
 * This information is carried through the noding process
 * (which may result in many noded edges sharing the same information object).
 * It is then used to populate the topology info fields
 * in {@link EdgeNG}s (possibly via merging).
 * That information is used to construct the topology graph {@link OverlayLabel}s.
 *
 * @author mdavis
 *
 */
class EdgeSourceInfo {
  int index;
  int dim = -999;
  bool isHole = false;
  int depthDelta = 0;

  EdgeSourceInfo(this.index, this.depthDelta, this.isHole) {
    this.dim = Dimension.A;
  }

  factory EdgeSourceInfo.fromIndex(int dex) {
    return EdgeSourceInfo(dex, 0, false);
  }

  int getIndex() {
    return index;
  }

  int getDimension() {
    return dim;
  }

  int getDepthDelta() {
    return depthDelta;
  }

  String toString() {
    return EdgeNG.infoString(index, dim, isHole, depthDelta);
  }
}
