part of dart_jts;
/*
 * Copyright (c) 2016 Vivid Solutions.
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse  License v1.0
 * and Eclipse Distribution License v. 1.0 which accompanies this distribution.
 * The Eclipse  License is available at http://www.eclipse.org/legal/epl-v10.html
 * and the Eclipse Distribution License is available at
 *
 * http://www.eclipse.org/org/documents/edl-v10.php.
 */

/**
 * An iterator over the components and coordinates of a linear geometry
 * ({@link LineString}s and {@link MultiLineString}s.
 *
 * The standard usage pattern for a {@link LinearIterator} is:
 *
 * <pre>
 * for (LinearIterator it = new LinearIterator(...); it.hasNext(); it.next()) {
 *   ...
 *   int ci = it.getComponentIndex();   // for example
 *   int vi = it.getVertexIndex();      // for example
 *   ...
 * }
 * </pre>
 *
 * @version 1.7
 */
class LinearIterator {
  static int segmentEndVertexIndex(LinearLocation loc) {
    if (loc.getSegmentFraction() > 0.0) return loc.getSegmentIndex() + 1;
    return loc.getSegmentIndex();
  }

  Geometry linearGeom;
  int numLines = 0;

  /**
   * Invariant: currentLine <> null if the iterator is pointing at a valid coordinate
   */
  LineString? currentLine;
  int componentIndex = 0;
  int vertexIndex = 0;

  /**
   * Creates an iterator initialized to the start of a linear {@link Geometry}
   *
   * @param linear the linear geometry to iterate over
   * @throws IllegalArgumentException if linearGeom is not lineal
   */
  LinearIterator(Geometry linear) : this.withIndexes(linear, 0, 0);

  /**
   * Creates an iterator starting at
   * a {@link LinearLocation} on a linear {@link Geometry}
   *
   * @param linear the linear geometry to iterate over
   * @param start the location to start at
   * @throws IllegalArgumentException if linearGeom is not lineal
   */
  LinearIterator.withStart(Geometry linear, LinearLocation start)
      : this.withIndexes(
            linear, start.getComponentIndex(), segmentEndVertexIndex(start));

  /**
   * Creates an iterator starting at
   * a specified component and vertex in a linear {@link Geometry}
   *
   * @param linearGeom the linear geometry to iterate over
   * @param componentIndex the component to start at
   * @param vertexIndex the vertex to start at
   * @throws IllegalArgumentException if linearGeom is not lineal
   */
  LinearIterator.withIndexes(
      this.linearGeom, int componentIndex, int vertexIndex) {
    if (!(linearGeom is Lineal))
      throw new ArgumentError("Lineal geometry is required");
    numLines = linearGeom.getNumGeometries();
    this.componentIndex = componentIndex;
    this.vertexIndex = vertexIndex;
    loadCurrentLine();
  }

  void loadCurrentLine() {
    if (componentIndex >= numLines) {
      currentLine = null;
      return;
    }
    currentLine = linearGeom.getGeometryN(componentIndex) as LineString;
  }

  /**
   * Tests whether there are any vertices left to iterator over.
   * Specifically, hasNext() return <tt>true</tt> if the
   * current state of the iterator represents a valid location
   * on the linear geometry. 
   * 
   * @return <code>true</code> if there are more vertices to scan
   */
  bool hasNext() {
    if (componentIndex >= numLines) return false;
    if (componentIndex == numLines - 1 &&
        vertexIndex >= currentLine!.getNumPoints()) return false;
    return true;
  }

  /**
   * Moves the iterator ahead to the next vertex and (possibly) linear component.
   */
  void next() {
    if (!hasNext()) return;

    vertexIndex++;
    if (vertexIndex >= currentLine!.getNumPoints()) {
      componentIndex++;
      loadCurrentLine();
      vertexIndex = 0;
    }
  }

  /**
   * Checks whether the iterator cursor is pointing to the
   * endpoint of a component {@link LineString}.
   *
   * @return <code>true</code> if the iterator is at an endpoint
   */
  bool isEndOfLine() {
    if (componentIndex >= numLines) return false;
    //LineString currentLine = (LineString) linear.getGeometryN(componentIndex);
    if (vertexIndex < currentLine!.getNumPoints() - 1) return false;
    return true;
  }

  /**
   * The component index of the vertex the iterator is currently at.
   * @return the current component index
   */
  int getComponentIndex() {
    return componentIndex;
  }

  /**
   * The vertex index of the vertex the iterator is currently at.
   * @return the current vertex index
   */
  int getVertexIndex() {
    return vertexIndex;
  }

  /**
   * Gets the {@link LineString} component the iterator is current at.
   * @return a linestring
   */
  LineString getLine() {
    return currentLine!;
  }

  /**
   * Gets the first {@link Coordinate} of the current segment.
   * (the coordinate of the current vertex).
   * @return a {@link Coordinate}
   */
  Coordinate getSegmentStart() {
    return currentLine!.getCoordinateN(vertexIndex);
  }

  /**
   * Gets the second {@link Coordinate} of the current segment.
   * (the coordinate of the next vertex).
   * If the iterator is at the end of a line, <code>null</code> is returned.
   *
   * @return a {@link Coordinate} or <code>null</code>
   */
  Coordinate? getSegmentEnd() {
    if (vertexIndex < getLine().getNumPoints() - 1)
      return currentLine!.getCoordinateN(vertexIndex + 1);
    return null;
  }
}
