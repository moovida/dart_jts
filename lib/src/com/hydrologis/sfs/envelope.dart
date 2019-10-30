part of dart_sfs;

class _Envelope {
  num minx;
  num miny;
  num maxx;
  num maxy;

  _Envelope(this.minx, this.miny, this.maxx, this.maxy);

  _Envelope.empty() : this(null, null, null, null);

  _Envelope.collapsed(num minx, num miny) : this(minx, miny, minx, miny);

  bool get isEmpty => minx == null;

  Geometry toGeometry() {
    if (isEmpty) {
      return Point.empty();
    }
    if (minx == maxx && miny == maxy) return Point(minx, maxy);

    //TODO: implement other kind of envelopes
    throw UnimplementedError();
  }

  growTo(other) {
    if (other is _Envelope) {
    } else if (other is Geometry) {
      other = other._envelope;
    } else if (other is Coordinate) {
      if (minx == null || other.x < minx) minx = other.x;
      if (miny == null || other.y < miny) miny = other.y;
      if (maxx == null || other.x > maxx) maxx = other.x;
      if (maxy == null || other.y > maxy) maxy = other.y;
      return;
    } else {
      _require(false, "unexpected type of argument, got ${other}");
    }
    if (minx == null || other.minx < minx) minx = other.minx;
    if (miny == null || other.miny < miny) miny = other.miny;
    if (maxx == null || other.maxx > maxx) maxx = other.maxx;
    if (maxy == null || other.maxy > maxy) maxy = other.maxy;
  }
}
