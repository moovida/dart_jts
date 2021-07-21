part of dart_jts;

class Envelope {
  ///  the minimum x-coordinate
  double _minx = 0;

  ///  the maximum x-coordinate
  double _maxx = 0;

  ///  the minimum y-coordinate
  double _miny = 0;

  ///  the maximum y-coordinate
  double _maxy = 0;

  /// Test the point q to see whether it intersects the Envelope defined by p1-p2
  /// @param p1 one extremal point of the envelope
  /// @param p2 another extremal point of the envelope
  /// @param q the point to test for intersection
  /// @return <code>true</code> if q intersects the envelope p1-p2
  static bool intersectsPoint(Coordinate p1, Coordinate p2, Coordinate q) {
    //OptimizeIt shows that Math#min and Math#max here are a bottleneck.
    //Replace with direct comparisons. [Jon Aquino]
    if (((q.x >= (p1.x < p2.x ? p1.x : p2.x)) &&
            (q.x <= (p1.x > p2.x ? p1.x : p2.x))) &&
        ((q.y >= (p1.y < p2.y ? p1.y : p2.y)) &&
            (q.y <= (p1.y > p2.y ? p1.y : p2.y)))) {
      return true;
    }
    return false;
  }

  /// Tests whether the envelope defined by p1-p2
  /// and the envelope defined by q1-q2
  /// intersect.
  ///
  /// @param p1 one extremal point of the envelope P
  /// @param p2 another extremal point of the envelope P
  /// @param q1 one extremal point of the envelope Q
  /// @param q2 another extremal point of the envelope Q
  /// @return <code>true</code> if Q intersects P
  static bool intersectsEnvelopeCoords(
      Coordinate p1, Coordinate p2, Coordinate q1, Coordinate q2) {
    double minq = math.min(q1.x, q2.x);
    double maxq = math.max(q1.x, q2.x);
    double minp = math.min(p1.x, p2.x);
    double maxp = math.max(p1.x, p2.x);

    if (minp > maxq) return false;
    if (maxp < minq) return false;

    minq = math.min(q1.y, q2.y);
    maxq = math.max(q1.y, q2.y);
    minp = math.min(p1.y, p2.y);
    maxp = math.max(p1.y, p2.y);

    if (minp > maxq) return false;
    if (maxp < minq) return false;
    return true;
  }

  ///  Creates a null <code>Envelope</code>.
  Envelope.empty() {
    initEmpty();
  }

  ///  Creates an <code>Envelope</code> for a region defined by maximum and minimum values.
  ///
  ///@param  x1  the first x-value
  ///@param  x2  the second x-value
  ///@param  y1  the first y-value
  ///@param  y2  the second y-value
  Envelope(double x1, double x2, double y1, double y2) {
    init(x1, x2, y1, y2);
  }

  ///  Creates an <code>Envelope</code> for a region defined by two Coordinates.
  ///
  ///@param  p1  the first Coordinate
  ///@param  p2  the second Coordinate
  Envelope.fromCoordinates(Coordinate p1, Coordinate p2) {
    init(p1.x, p2.x, p1.y, p2.y);
  }

  ///  Creates an <code>Envelope</code> for a region defined by a single Coordinate.
  ///
  ///@param  p  the Coordinate
  Envelope.fromCoordinate(Coordinate p) {
    init(p.x, p.x, p.y, p.y);
  }

  ///  Create an <code>Envelope</code> from an existing Envelope.
  ///
  ///@param  env  the Envelope to initialize from
  Envelope.fromEnvelope(Envelope env) {
    initFromEnvelope(env);
  }

  ///  Initialize to a null <code>Envelope</code>.
  void initEmpty() {
    setToNull();
  }

  ///  Initialize an <code>Envelope</code> for a region defined by maximum and minimum values.
  ///
  ///@param  x1  the first x-value
  ///@param  x2  the second x-value
  ///@param  y1  the first y-value
  ///@param  y2  the second y-value
  void init(double x1, double x2, double y1, double y2) {
    if (x1 < x2) {
      _minx = x1;
      _maxx = x2;
    } else {
      _minx = x2;
      _maxx = x1;
    }
    if (y1 < y2) {
      _miny = y1;
      _maxy = y2;
    } else {
      _miny = y2;
      _maxy = y1;
    }
  }

  /// Creates a copy of this envelope object.
  ///
  /// @return a copy of this envelope
  Envelope copy() {
    return Envelope.fromEnvelope(this);
  }

  ///  Initialize an <code>Envelope</code> to a region defined by two Coordinates.
  ///
  ///@param  p1  the first Coordinate
  ///@param  p2  the second Coordinate
  void initFromCoordinates(Coordinate p1, Coordinate p2) {
    init(p1.x, p2.x, p1.y, p2.y);
  }

  ///  Initialize an <code>Envelope</code> to a region defined by a single Coordinate.
  ///
  ///@param  p  the coordinate
  void initFromCoordinate(Coordinate p) {
    init(p.x, p.x, p.y, p.y);
  }

  ///  Initialize an <code>Envelope</code> from an existing Envelope.
  ///
  ///@param  env  the Envelope to initialize from
  void initFromEnvelope(Envelope env) {
    this._minx = env._minx;
    this._maxx = env._maxx;
    this._miny = env._miny;
    this._maxy = env._maxy;
  }

  ///  Makes this <code>Envelope</code> a "null" envelope, that is, the envelope
  ///  of the empty geometry.
  void setToNull() {
    _minx = 0;
    _maxx = -1;
    _miny = 0;
    _maxy = -1;
  }

  ///  Returns <code>true</code> if this <code>Envelope</code> is a "null"
  ///  envelope.
  ///
  ///@return    <code>true</code> if this <code>Envelope</code> is uninitialized
  ///      or is the envelope of the empty geometry.
  bool isNull() {
    return _maxx < _minx;
  }

  ///  Returns the difference between the maximum and minimum x values.
  ///
  ///@return    max x - min x, or 0 if this is a null <code>Envelope</code>
  double getWidth() {
    if (isNull()) {
      return 0;
    }
    return _maxx - _minx;
  }

  ///  Returns the difference between the maximum and minimum y values.
  ///
  ///@return    max y - min y, or 0 if this is a null <code>Envelope</code>
  double getHeight() {
    if (isNull()) {
      return 0;
    }
    return _maxy - _miny;
  }

  /// Gets the length of the diameter (diagonal) of the envelope.
  ///
  /// @return the diameter length
  double getDiameter() {
    if (isNull()) {
      return 0;
    }
    double w = getWidth();
    double h = getHeight();
    return math.sqrt(w * w + h * h);
  }

  ///  Returns the <code>Envelope</code>s minimum x-value. min x &gt; max x
  ///  indicates that this is a null <code>Envelope</code>.
  ///
  ///@return    the minimum x-coordinate
  double getMinX() {
    return _minx;
  }

  ///  Returns the <code>Envelope</code>s maximum x-value. min x &gt; max x
  ///  indicates that this is a null <code>Envelope</code>.
  ///
  ///@return    the maximum x-coordinate
  double getMaxX() {
    return _maxx;
  }

  ///  Returns the <code>Envelope</code>s minimum y-value. min y &gt; max y
  ///  indicates that this is a null <code>Envelope</code>.
  ///
  ///@return    the minimum y-coordinate
  double getMinY() {
    return _miny;
  }

  ///  Returns the <code>Envelope</code>s maximum y-value. min y &gt; max y
  ///  indicates that this is a null <code>Envelope</code>.
  ///
  ///@return    the maximum y-coordinate
  double getMaxY() {
    return _maxy;
  }

  /// Gets the area of this envelope.
  ///
  /// @return the area of the envelope
  /// @return 0.0 if the envelope is null
  double getArea() {
    return getWidth() * getHeight();
  }

  /// Gets the minimum extent of this envelope across both dimensions.
  ///
  /// @return the minimum extent of this envelope
  double minExtent() {
    if (isNull()) return 0.0;
    double w = getWidth();
    double h = getHeight();
    if (w < h) return w;
    return h;
  }

  /// Gets the maximum extent of this envelope across both dimensions.
  ///
  /// @return the maximum extent of this envelope
  double maxExtent() {
    if (isNull()) return 0.0;
    double w = getWidth();
    double h = getHeight();
    if (w > h) return w;
    return h;
  }

  ///  Enlarges this <code>Envelope</code> so that it contains
  ///  the given {@link Coordinate}.
  ///  Has no effect if the point is already on or within the envelope.
  ///
  ///@param  p  the Coordinate to expand to include
  void expandToIncludeCoordinate(Coordinate p) {
    expandToInclude(p.x, p.y);
  }

  /// Expands this envelope by a given distance in all directions.
  /// Both positive and negative distances are supported.
  ///
  /// @param distance the distance to expand the envelope
  void expandByDistance(double distance) {
    expandBy(distance, distance);
  }

  /// Expands this envelope by a given distance in all directions.
  /// Both positive and negative distances are supported.
  ///
  /// @param deltaX the distance to expand the envelope along the the X axis
  /// @param deltaY the distance to expand the envelope along the the Y axis
  void expandBy(double deltaX, double deltaY) {
    if (isNull()) return;

    _minx -= deltaX;
    _maxx += deltaX;
    _miny -= deltaY;
    _maxy += deltaY;

    // check for envelope disappearing
    if (_minx > _maxx || _miny > _maxy) {
      setToNull();
    }
  }

  ///  Enlarges this <code>Envelope</code> so that it contains
  ///  the given point.
  ///  Has no effect if the point is already on or within the envelope.
  ///
  ///@param  x  the value to lower the minimum x to or to raise the maximum x to
  ///@param  y  the value to lower the minimum y to or to raise the maximum y to
  void expandToInclude(double x, double y) {
    if (isNull()) {
      _minx = x;
      _maxx = x;
      _miny = y;
      _maxy = y;
    } else {
      if (x < _minx) {
        _minx = x;
      }
      if (x > _maxx) {
        _maxx = x;
      }
      if (y < _miny) {
        _miny = y;
      }
      if (y > _maxy) {
        _maxy = y;
      }
    }
  }

  ///  Enlarges this <code>Envelope</code> so that it contains
  ///  the <code>other</code> Envelope.
  ///  Has no effect if <code>other</code> is wholly on or
  ///  within the envelope.
  ///
  ///@param  other  the <code>Envelope</code> to expand to include
  void expandToIncludeEnvelope(Envelope other) {
    if (other.isNull()) {
      return;
    }
    if (isNull()) {
      _minx = other.getMinX();
      _maxx = other.getMaxX();
      _miny = other.getMinY();
      _maxy = other.getMaxY();
    } else {
      if (other._minx < _minx) {
        _minx = other._minx;
      }
      if (other._maxx > _maxx) {
        _maxx = other._maxx;
      }
      if (other._miny < _miny) {
        _miny = other._miny;
      }
      if (other._maxy > _maxy) {
        _maxy = other._maxy;
      }
    }
  }

  /// Translates this envelope by given amounts in the X and Y direction.
  ///
  /// @param transX the amount to translate along the X axis
  /// @param transY the amount to translate along the Y axis
  void translate(double transX, double transY) {
    if (isNull()) {
      return;
    }
    init(getMinX() + transX, getMaxX() + transX, getMinY() + transY,
        getMaxY() + transY);
  }

  /// Computes the coordinate of the centre of this envelope (as long as it is non-null
  ///
  /// @return the centre coordinate of this envelope
  /// <code>null</code> if the envelope is null
  Coordinate? centre() {
    if (isNull()) return null;
    return Coordinate(
        (getMinX() + getMaxX()) / 2.0, (getMinY() + getMaxY()) / 2.0);
  }

  /// Computes the intersection of two {@link Envelope}s.
  ///
  /// @param env the envelope to intersect with
  /// @return a new Envelope representing the intersection of the envelopes (this will be
  /// the null envelope if either argument is null, or they do not intersect
  Envelope intersection(Envelope env) {
    if (isNull() || env.isNull() || !intersectsEnvelope(env))
      return Envelope.empty();

    double intMinX = _minx > env._minx ? _minx : env._minx;
    double intMinY = _miny > env._miny ? _miny : env._miny;
    double intMaxX = _maxx < env._maxx ? _maxx : env._maxx;
    double intMaxY = _maxy < env._maxy ? _maxy : env._maxy;
    return Envelope(intMinX, intMaxX, intMinY, intMaxY);
  }

  /// Tests if the region defined by <code>other</code>
  /// intersects the region of this <code>Envelope</code>.
  ///
  ///@param  other  the <code>Envelope</code> which this <code>Envelope</code> is
  ///          being checked for intersecting
  ///@return        <code>true</code> if the <code>Envelope</code>s intersect
  bool intersectsEnvelope(Envelope other) {
    if (isNull() || other.isNull()) {
      return false;
    }
    return !(other._minx > _maxx ||
        other._maxx < _minx ||
        other._miny > _maxy ||
        other._maxy < _miny);
  }

  /// Tests if the extent defined by two extremal points
  /// intersects the extent of this <code>Envelope</code>.
  ///
  ///@param a a point
  ///@param b another point
  ///@return   <code>true</code> if the extents intersect
  bool intersectsEnvelopeCoordinates(Coordinate a, Coordinate b) {
    if (isNull()) {
      return false;
    }

    double envminx = (a.x < b.x) ? a.x : b.x;
    if (envminx > _maxx) return false;

    double envmaxx = (a.x > b.x) ? a.x : b.x;
    if (envmaxx < _minx) return false;

    double envminy = (a.y < b.y) ? a.y : b.y;
    if (envminy > _maxy) return false;

    double envmaxy = (a.y > b.y) ? a.y : b.y;
    if (envmaxy < _miny) return false;

    return true;
  }

  /// Tests if the region defined by <code>other</code>
  /// is disjoint from the region of this <code>Envelope</code>.
  ///
  ///@param  other  the <code>Envelope</code> being checked for disjointness
  ///@return        <code>true</code> if the <code>Envelope</code>s are disjoint
  ///
  ///@see #intersects(Envelope)
  bool disjoint(Envelope other) {
    if (isNull() || other.isNull()) {
      return true;
    }
    return other._minx > _maxx ||
        other._maxx < _minx ||
        other._miny > _maxy ||
        other._maxy < _miny;
  }

  /// Tests if the point <code>p</code>
  /// intersects (lies inside) the region of this <code>Envelope</code>.
  ///
  ///@param  p  the <code>Coordinate</code> to be tested
  ///@return <code>true</code> if the point intersects this <code>Envelope</code>
  bool intersectsCoordinate(Coordinate p) {
    return intersects(p.x, p.y);
  }

  ///  Check if the point <code>(x, y)</code>
  ///  intersects (lies inside) the region of this <code>Envelope</code>.
  ///
  ///@param  x  the x-ordinate of the point
  ///@param  y  the y-ordinate of the point
  ///@return        <code>true</code> if the point overlaps this <code>Envelope</code>
  bool intersects(double x, double y) {
    if (isNull()) return false;
    return !(x > _maxx || x < _minx || y > _maxy || y < _miny);
  }

  /// Tests if the <code>Envelope other</code>
  /// lies wholely inside this <code>Envelope</code> (inclusive of the boundary).
  /// <p>
  /// Note that this is <b>not</b> the same definition as the SFS <tt>contains</tt>,
  /// which would exclude the envelope boundary.
  ///
  ///@param  other the <code>Envelope</code> to check
  ///@return true if <code>other</code> is contained in this <code>Envelope</code>
  ///
  ///@see #covers(Envelope)
  bool containsEnvelope(Envelope other) {
    return coversEnvelope(other);
  }

  /// Tests if the given point lies in or on the envelope.
  /// <p>
  /// Note that this is <b>not</b> the same definition as the SFS <tt>contains</tt>,
  /// which would exclude the envelope boundary.
  ///
  ///@param  p  the point which this <code>Envelope</code> is
  ///      being checked for containing
  ///@return    <code>true</code> if the point lies in the interior or
  ///      on the boundary of this <code>Envelope</code>.
  ///
  ///@see #covers(Coordinate)
  bool containsCoordinate(Coordinate p) {
    return coversCoordinate(p);
  }

  /// Tests if the given point lies in or on the envelope.
  /// <p>
  /// Note that this is <b>not</b> the same definition as the SFS <tt>contains</tt>,
  /// which would exclude the envelope boundary.
  ///
  ///@param  x  the x-coordinate of the point which this <code>Envelope</code> is
  ///      being checked for containing
  ///@param  y  the y-coordinate of the point which this <code>Envelope</code> is
  ///      being checked for containing
  ///@return    <code>true</code> if <code>(x, y)</code> lies in the interior or
  ///      on the boundary of this <code>Envelope</code>.
  ///
  ///@see #covers(double, double)
  bool contains(double x, double y) {
    return covers(x, y);
  }

  /// Tests if the given point lies in or on the envelope.
  ///
  ///@param  x  the x-coordinate of the point which this <code>Envelope</code> is
  ///      being checked for containing
  ///@param  y  the y-coordinate of the point which this <code>Envelope</code> is
  ///      being checked for containing
  ///@return    <code>true</code> if <code>(x, y)</code> lies in the interior or
  ///      on the boundary of this <code>Envelope</code>.
  bool covers(double x, double y) {
    if (isNull()) return false;
    return x >= _minx && x <= _maxx && y >= _miny && y <= _maxy;
  }

  /// Tests if the given point lies in or on the envelope.
  ///
  ///@param  p  the point which this <code>Envelope</code> is
  ///      being checked for containing
  ///@return    <code>true</code> if the point lies in the interior or
  ///      on the boundary of this <code>Envelope</code>.
  bool coversCoordinate(Coordinate p) {
    return covers(p.x, p.y);
  }

  /// Tests if the <code>Envelope other</code>
  /// lies wholely inside this <code>Envelope</code> (inclusive of the boundary).
  ///
  ///@param  other the <code>Envelope</code> to check
  ///@return true if this <code>Envelope</code> covers the <code>other</code>
  bool coversEnvelope(Envelope other) {
    if (isNull() || other.isNull()) {
      return false;
    }
    return other.getMinX() >= _minx &&
        other.getMaxX() <= _maxx &&
        other.getMinY() >= _miny &&
        other.getMaxY() <= _maxy;
  }

  /// Computes the distance between this and another
  /// <code>Envelope</code>.
  /// The distance between overlapping Envelopes is 0.  Otherwise, the
  /// distance is the Euclidean distance between the closest points.
  double distance(Envelope env) {
    if (intersectsEnvelope(env)) return 0;

    double dx = 0.0;
    if (_maxx < env._minx) {
      dx = env._minx - _maxx;
    } else if (_minx > env._maxx) dx = _minx - env._maxx;

    double dy = 0.0;
    if (_maxy < env._miny) {
      dy = env._miny - _maxy;
    } else if (_miny > env._maxy) dy = _miny - env._maxy;

    // if either is zero, the envelopes overlap either vertically or horizontally
    if (dx == 0.0) return dy;
    if (dy == 0.0) return dx;
    return math.sqrt(dx * dx + dy * dy);
  }

  bool operator ==(other) {
    if (!(other is Envelope)) {
      return false;
    }
    Envelope otherEnvelope = other;
    if (isNull()) {
      return otherEnvelope.isNull();
    }
    return _maxx == otherEnvelope.getMaxX() &&
        _maxy == otherEnvelope.getMaxY() &&
        _minx == otherEnvelope.getMinX() &&
        _miny == otherEnvelope.getMinY();
  }

  int get hashCode {
    //Algorithm from Effective Java by Joshua Bloch [Jon Aquino]
    int result = 17;
    result = 37 * result + _minx.hashCode;
    result = 37 * result + _maxx.hashCode;
    result = 37 * result + _miny.hashCode;
    result = 37 * result + _maxy.hashCode;
    return result;
  }

  String toString() {
    return "Env[$_minx : $_maxx , $_miny : $_maxy ]";
  }

  /// Compares two envelopes using lexicographic ordering.
  /// The ordering comparison is based on the usual numerical
  /// comparison between the sequence of ordinates.
  /// Null envelopes are less than all non-null envelopes.
  ///
  /// @param o an Envelope object
  int compareTo(Object o) {
    Envelope env = o as Envelope;
    // compare nulls if present
    if (isNull()) {
      if (env.isNull()) return 0;
      return -1;
    } else {
      if (env.isNull()) return 1;
    }
    // compare based on numerical ordering of ordinates
    if (_minx < env._minx) return -1;
    if (_minx > env._minx) return 1;
    if (_miny < env._miny) return -1;
    if (_miny > env._miny) return 1;
    if (_maxx < env._maxx) return -1;
    if (_maxx > env._maxx) return 1;
    if (_maxy < env._maxy) return -1;
    if (_maxy > env._maxy) return 1;
    return 0;
  }
}
