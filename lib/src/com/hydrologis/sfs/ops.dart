part of dart_sfs;
/*
 * Copyright (c) 2016 Vivid Solutions.
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * and Eclipse Distribution License v. 1.0 which accompanies this distribution.
 * The Eclipse Public License is available at http://www.eclipse.org/legal/epl-v10.html
 * and the Eclipse Distribution License is available at
 *
 * http://www.eclipse.org/org/documents/edl-v10.php.
 */

/// A {@link BoundaryNodeRule} specifies that points are in the
/// boundary of a lineal geometry iff
/// the point lies on the boundary of an odd number
/// of components.
/// Under this rule {@link LinearRing}s and closed
/// {@link LineString}s have an empty boundary.
/// <p>
/// This is the rule specified by the <i>OGC SFS</i>,
/// and is the default rule used in JTS.
///
/// @author Martin Davis
/// @version 1.7
class Mod2BoundaryNodeRule {
  isInBoundary(int boundaryCount) {
    // the "Mod-2 Rule"
    return boundaryCount % 2 == 1;
  }
}

/// Computes the boundary of a {@link Geometry}.
/// This operation will always return a {@link Geometry} of the appropriate
/// dimension for the boundary (even if the input geometry is empty).
/// The boundary of zero-dimensional geometries (Points) is
/// always the empty {@link GeometryCollection}.
///
/// @author Martin Davis
/// @author Andrea Antonello - dart port
/// @version 1.7

class BoundaryOp {
  static Geometry getBoundaryFromGeometry(Geometry g) {
    BoundaryOp bop = BoundaryOp(g);
    return bop.getBoundary();
  }

  Geometry geom;
  GeometryFactory geomFact = GeometryFactory();
  Mod2BoundaryNodeRule bnRule = Mod2BoundaryNodeRule();

  BoundaryOp(Geometry geom) {
    this.geom = geom;
  }

  Geometry getBoundary() {
    if (geom is LineString) return boundaryLineString(geom);
    if (geom is MultiLineString) {
      return boundaryMultiLineString(geom);
    }
    return geom.boundary;
  }

  MultiPoint getEmptyMultiPoint() {
    return geomFact.createMultiPoint();
  }

  Geometry boundaryMultiLineString(MultiLineString mLine) {
    if (geom.isEmpty) {
      return getEmptyMultiPoint();
    }

    List<Coordinate> bdyPts = computeBoundaryCoordinates(mLine);

    // return Point or MultiPoint
    if (bdyPts.length == 1) {
      return geomFact.createPoint(bdyPts[0]);
    }
    // this handles 0 points case as well
    return geomFact.createMultiPointFromCoords(bdyPts);
  }

/*
// MD - superseded
  private Coordinate[] computeBoundaryFromGeometryGraph(MultiLineString mLine)
  {
    GeometryGraph g = new GeometryGraph(0, mLine, bnRule);
    Coordinate[] bdyPts = g.getBoundaryPoints();
    return bdyPts;
  }
*/

  Map<Coordinate, Counter> endpointMap;

  List<Coordinate> computeBoundaryCoordinates(MultiLineString mLine) {
    List<Coordinate> bdyPts = [];
    endpointMap = SplayTreeMap();
    for (int i = 0; i < mLine.getNumGeometries(); i++) {
      LineString line = mLine.getGeometryN(i);
      if (line.isEmpty) {
        continue;
      }
      addEndpoint(line.getCoordinateN(0));
      addEndpoint(line.getCoordinateN(line.length - 1));
    }

    endpointMap.forEach((coord, counter) {
      int valence = counter.count;
      if (bnRule.isInBoundary(valence)) {
        bdyPts.add(coord);
      }
    });

    return bdyPts;
  }

  void addEndpoint(Coordinate pt) {
    Counter counter = endpointMap[pt];
    if (counter == null) {
      counter = Counter();
      endpointMap[pt] = counter;
    }
    counter.count++;
  }

  Geometry boundaryLineString(LineString line) {
    if (geom.isEmpty) {
      return getEmptyMultiPoint();
    }

    if (line.isClosed) {
      // check whether endpoints of valence 2 are on the boundary or not
      bool closedEndpointOnBoundary = bnRule.isInBoundary(2);
      if (closedEndpointOnBoundary) {
        return line.getStartPoint();
      } else {
        return geomFact.createMultiPoint();
      }
    }
    return geomFact
        .createMultiPoint([line.getStartPoint(), line.getEndPoint()]);
  }
}

/// Stores an integer count, for use as a Map entry.
///
/// @author Martin Davis
/// @version 1.7
class Counter {
  /// The value of the count
  int count;
}
