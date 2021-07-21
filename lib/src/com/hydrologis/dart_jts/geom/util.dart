part of dart_jts;

/**
 * A visitor to {@link Geometry} components, which
 * allows short-circuiting when a defined condition holds.
 *
 * @version 1.7
 */
abstract class ShortCircuitedGeometryVisitor {
  bool _isDone = false;

  ShortCircuitedGeometryVisitor() {}

  void applyTo(Geometry geom) {
    for (int i = 0; i < geom.getNumGeometries() && !_isDone; i++) {
      Geometry element = geom.getGeometryN(i);
      if (!(element is GeometryCollection)) {
        visit(element);
        if (isDone()) {
          _isDone = true;
          return;
        }
      } else
        applyTo(element);
    }
  }

  void visit(Geometry element);

  /**
   * Reports whether visiting components can be terminated.
   * Once this method returns <tt>true</tt>, it must
   * continue to return <tt>true</tt> on every subsequent call.
   *
   * @return true if visiting can be terminated.
   */
  bool isDone();
}

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

///  Constants representing the different topological locations
///  which can occur in a {@link Geometry}.
///  The constants are also used as the row and column indices
///  of DE-9IM {@link IntersectionMatrix}es.
///
///@version 1.7
class Location {
  /// The location value for the interior of a geometry.
  /// Also, DE-9IM row index of the interior of the first geometry and column index of
  ///  the interior of the second geometry.
  static const int INTERIOR = 0;

  /// The location value for the boundary of a geometry.
  /// Also, DE-9IM row index of the boundary of the first geometry and column index of
  ///  the boundary of the second geometry.
  static const int BOUNDARY = 1;

  /// The location value for the exterior of a geometry.
  /// Also, DE-9IM row index of the exterior of the first geometry and column index of
  ///  the exterior of the second geometry.
  static const int EXTERIOR = 2;

  ///  Used for uninitialized location values.
  static const int NONE = -1;

  ///  Converts the location value to a location symbol, for example, <code>EXTERIOR =&gt; 'e'</code>
  ///  .
  ///
  ///@param  locationValue  either EXTERIOR, BOUNDARY, INTERIOR or NONE
  ///@return                either 'e', 'b', 'i' or '-'
  static String toLocationSymbol(int locationValue) {
    switch (locationValue) {
      case EXTERIOR:
        return 'e';
      case BOUNDARY:
        return 'b';
      case INTERIOR:
        return 'i';
      case NONE:
        return '-';
    }
    throw ArgumentError("Unknown location value: $locationValue");
  }
}

/// Provides constants representing the dimensions of a point, a curve and a surface.
/// Also provides constants representing the dimensions of the empty geometry and
/// non-empty geometries, and the wildcard constant {@link #DONTCARE} meaning "any dimension".
/// These constants are used as the entries in {@link IntersectionMatrix}s.
///
/// @version 1.7
class Dimension {
  ///  Dimension value of a point (0).
  static const int P = 0;

  ///  Dimension value of a curve (1).
  static const int L = 1;

  ///  Dimension value of a surface (2).
  static const int A = 2;

  ///  Dimension value of the empty geometry (-1).
  static const int FALSE = -1;

  ///  Dimension value of non-empty geometries (= {P, L, A}).
  static const int TRUE = -2;

  ///  Dimension value for any dimension (= {FALSE, TRUE}).
  static const int DONTCARE = -3;

  /// Symbol for the FALSE pattern matrix entry
  static const String SYM_FALSE = 'F';

  /// Symbol for the TRUE pattern matrix entry
  static const String SYM_TRUE = 'T';

  /// Symbol for the DONTCARE pattern matrix entry
  static const String SYM_DONTCARE = '*';

  /// Symbol for the P (dimension 0) pattern matrix entry
  static const String SYM_P = '0';

  /// Symbol for the L (dimension 1) pattern matrix entry
  static const String SYM_L = '1';

  /// Symbol for the A (dimension 2) pattern matrix entry
  static const String SYM_A = '2';

  ///  Converts the dimension value to a dimension symbol, for example, <code>TRUE =&gt; 'T'</code>
  ///  .
  ///
  ///@param  dimensionValue  a number that can be stored in the <code>IntersectionMatrix</code>
  ///      . Possible values are <code>{TRUE, FALSE, DONTCARE, 0, 1, 2}</code>.
  ///@return                 a character for use in the string representation of
  ///      an <code>IntersectionMatrix</code>. Possible values are <code>{T, F, * , 0, 1, 2}</code>
  ///      .
  static String toDimensionSymbol(int dimensionValue) {
    switch (dimensionValue) {
      case FALSE:
        return SYM_FALSE;
      case TRUE:
        return SYM_TRUE;
      case DONTCARE:
        return SYM_DONTCARE;
      case P:
        return SYM_P;
      case L:
        return SYM_L;
      case A:
        return SYM_A;
    }
    throw ArgumentError("Unknown dimension value: $dimensionValue");
  }

  ///  Converts the dimension symbol to a dimension value, for example, <code>'*' =&gt; DONTCARE</code>
  ///  .
  ///
  ///@param  dimensionSymbol  a character for use in the string representation of
  ///      an <code>IntersectionMatrix</code>. Possible values are <code>{T, F, * , 0, 1, 2}</code>
  ///      .
  ///@return a number that can be stored in the <code>IntersectionMatrix</code>
  ///      . Possible values are <code>{TRUE, FALSE, DONTCARE, 0, 1, 2}</code>.
  static int toDimensionValue(String dimensionSymbol) {
    switch (dimensionSymbol.toUpperCase()) {
      case SYM_FALSE:
        return FALSE;
      case SYM_TRUE:
        return TRUE;
      case SYM_DONTCARE:
        return DONTCARE;
      case SYM_P:
        return P;
      case SYM_L:
        return L;
      case SYM_A:
        return A;
    }
    throw ArgumentError("Unknown dimension symbol: $dimensionSymbol");
  }
}

/// Models a <b>Dimensionally Extended Nine-Intersection Model (DE-9IM)</b> matrix.
/// DE-9IM matrices (such as "212FF1FF2")
/// specify the topological relationship between two {@link Geometry}s.
/// This class can also represent matrix patterns (such as "T*T******")
/// which are used for matching instances of DE-9IM matrices.
///
///  Methods are provided to:
///  <UL>
///    <LI> set and query the elements of the matrix in a convenient fashion
///    <LI> convert to and from the standard string representation (specified in
///    SFS Section 2.1.13.2).
///    <LI> test to see if a matrix matches a given pattern string.
///  </UL>
///  <P>
///
///  For a description of the DE-9IM and the spatial predicates derived from it,
///  see the <i><A
///  HREF="http://www.opengis.org/techno/specs.htm">OGC 99-049 OpenGIS Simple Features
///  Specification for SQL</A></i>, as well as
///  <i>OGC 06-103r4 OpenGIS
///  Implementation Standard for Geographic information -
///  Simple feature access - Part 1: Common architecture</i>
///  (which provides some further details on certain predicate specifications).
/// <p>
/// The entries of the matrix are defined by the constants in the {@link Dimension} class.
/// The indices of the matrix represent the topological locations
/// that occur in a geometry (Interior, Boundary, Exterior).
/// These are provided as constants in the {@link Location} class.
///
///
///@version 1.7
class IntersectionMatrix {
  ///  Internal representation of this <code>IntersectionMatrix</code>.
  late List<List<int>> matrix;

  ///  Creates an <code>IntersectionMatrix</code> with <code>FALSE</code>
  ///  dimension values.
  IntersectionMatrix() {
    _init();
  }

  void _init() {
    matrix = MatrixUtils.createIntMatrix(3, 3);
    setAll(Dimension.FALSE);
  }

  ///  Creates an <code>IntersectionMatrix</code> with the given dimension
  ///  symbols.
  ///
  ///@param  elements  a String of nine dimension symbols in row major order
  IntersectionMatrix.fromDimensionSymbols(String elements) {
    _init();
    setDimensionSimbols(elements);
  }

  ///  Creates an <code>IntersectionMatrix</code> with the same elements as
  ///  <code>other</code>.
  ///
  ///@param  other  an <code>IntersectionMatrix</code> to copy
  IntersectionMatrix.fromMatrix(IntersectionMatrix other) {
    _init();
    matrix[Location.INTERIOR][Location.INTERIOR] =
        other.matrix[Location.INTERIOR][Location.INTERIOR];
    matrix[Location.INTERIOR][Location.BOUNDARY] =
        other.matrix[Location.INTERIOR][Location.BOUNDARY];
    matrix[Location.INTERIOR][Location.EXTERIOR] =
        other.matrix[Location.INTERIOR][Location.EXTERIOR];
    matrix[Location.BOUNDARY][Location.INTERIOR] =
        other.matrix[Location.BOUNDARY][Location.INTERIOR];
    matrix[Location.BOUNDARY][Location.BOUNDARY] =
        other.matrix[Location.BOUNDARY][Location.BOUNDARY];
    matrix[Location.BOUNDARY][Location.EXTERIOR] =
        other.matrix[Location.BOUNDARY][Location.EXTERIOR];
    matrix[Location.EXTERIOR][Location.INTERIOR] =
        other.matrix[Location.EXTERIOR][Location.INTERIOR];
    matrix[Location.EXTERIOR][Location.BOUNDARY] =
        other.matrix[Location.EXTERIOR][Location.BOUNDARY];
    matrix[Location.EXTERIOR][Location.EXTERIOR] =
        other.matrix[Location.EXTERIOR][Location.EXTERIOR];
  }

  /// Adds one matrix to another.
  /// Addition is defined by taking the maximum dimension value of each position
  /// in the summand matrices.
  ///
  /// @param im the matrix to add
  void add(IntersectionMatrix im) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        setAtLeast(i, j, im.get(i, j));
      }
    }
  }

  ///  Tests if the dimension value matches <tt>TRUE</tt>
  ///  (i.e.  has value 0, 1, 2 or TRUE).
  ///
  ///@param  actualDimensionValue     a number that can be stored in the <code>IntersectionMatrix</code>
  ///      . Possible values are <code>{TRUE, FALSE, DONTCARE, 0, 1, 2}</code>.
  ///@return true if the dimension value matches TRUE
  static bool isTrue(int actualDimensionValue) {
    if (actualDimensionValue >= 0 || actualDimensionValue == Dimension.TRUE) {
      return true;
    }
    return false;
  }

  ///  Tests if the dimension value satisfies the dimension symbol.
  ///
  ///@param  actualDimensionValue     a number that can be stored in the <code>IntersectionMatrix</code>
  ///      . Possible values are <code>{TRUE, FALSE, DONTCARE, 0, 1, 2}</code>.
  ///@param  requiredDimensionSymbol  a character used in the string
  ///      representation of an <code>IntersectionMatrix</code>. Possible values
  ///      are <code>{T, F, * , 0, 1, 2}</code>.
  ///@return                          true if the dimension symbol matches
  ///      the dimension value
  static bool matchesDimValue(
      int actualDimensionValue, String requiredDimensionSymbol) {
    if (requiredDimensionSymbol == Dimension.SYM_DONTCARE) {
      return true;
    }
    if (requiredDimensionSymbol == Dimension.SYM_TRUE &&
        (actualDimensionValue >= 0 || actualDimensionValue == Dimension.TRUE)) {
      return true;
    }
    if (requiredDimensionSymbol == Dimension.SYM_FALSE &&
        actualDimensionValue == Dimension.FALSE) {
      return true;
    }
    if (requiredDimensionSymbol == Dimension.SYM_P &&
        actualDimensionValue == Dimension.P) {
      return true;
    }
    if (requiredDimensionSymbol == Dimension.SYM_L &&
        actualDimensionValue == Dimension.L) {
      return true;
    }
    if (requiredDimensionSymbol == Dimension.SYM_A &&
        actualDimensionValue == Dimension.A) {
      return true;
    }
    return false;
  }

  ///  Tests if each of the actual dimension symbols in a matrix string satisfies the
  ///  corresponding required dimension symbol in a pattern string.
  ///
  ///@param  actualDimensionSymbols    nine dimension symbols to validate.
  ///      Possible values are <code>{T, F, * , 0, 1, 2}</code>.
  ///@param  requiredDimensionSymbols  nine dimension symbols to validate
  ///      against. Possible values are <code>{T, F, * , 0, 1, 2}</code>.
  ///@return                           true if each of the required dimension
  ///      symbols encompass the corresponding actual dimension symbol
  static bool matchesDimSymbols(
      String actualDimensionSymbols, String requiredDimensionSymbols) {
    IntersectionMatrix m =
        IntersectionMatrix.fromDimensionSymbols(actualDimensionSymbols);
    return m.matches(requiredDimensionSymbols);
  }

  ///  Changes the value of one of this <code>IntersectionMatrix</code>s
  ///  elements.
  ///
  ///@param  row             the row of this <code>IntersectionMatrix</code>,
  ///      indicating the interior, boundary or exterior of the first <code>Geometry</code>
  ///@param  column          the column of this <code>IntersectionMatrix</code>,
  ///      indicating the interior, boundary or exterior of the second <code>Geometry</code>
  ///@param  dimensionValue  the new value of the element
  void set(int row, int column, int dimensionValue) {
    matrix[row][column] = dimensionValue;
  }

  ///  Changes the elements of this <code>IntersectionMatrix</code> to the
  ///  dimension symbols in <code>dimensionSymbols</code>.
  ///
  ///@param  dimensionSymbols  nine dimension symbols to which to set this <code>IntersectionMatrix</code>
  ///      s elements. Possible values are <code>{T, F, * , 0, 1, 2}</code>
  void setDimensionSimbols(String dimensionSymbols) {
    for (int i = 0; i < dimensionSymbols.length; i++) {
      int row = (i / 3).floor();
      int col = i % 3;
      matrix[row][col] = Dimension.toDimensionValue(dimensionSymbols[i]);
    }
  }

  ///  Changes the specified element to <code>minimumDimensionValue</code> if the
  ///  element is less.
  ///
  ///@param  row                    the row of this <code>IntersectionMatrix</code>
  ///      , indicating the interior, boundary or exterior of the first <code>Geometry</code>
  ///@param  column                 the column of this <code>IntersectionMatrix</code>
  ///      , indicating the interior, boundary or exterior of the second <code>Geometry</code>
  ///@param  minimumDimensionValue  the dimension value with which to compare the
  ///      element. The order of dimension values from least to greatest is
  ///      <code>{DONTCARE, TRUE, FALSE, 0, 1, 2}</code>.
  void setAtLeast(int row, int column, int minimumDimensionValue) {
    if (matrix[row][column] < minimumDimensionValue) {
      matrix[row][column] = minimumDimensionValue;
    }
  }

  ///  If row &gt;= 0 and column &gt;= 0, changes the specified element to <code>minimumDimensionValue</code>
  ///  if the element is less. Does nothing if row &lt;0 or column &lt; 0.
  ///
  ///@param  row                    the row of this <code>IntersectionMatrix</code>
  ///      , indicating the interior, boundary or exterior of the first <code>Geometry</code>
  ///@param  column                 the column of this <code>IntersectionMatrix</code>
  ///      , indicating the interior, boundary or exterior of the second <code>Geometry</code>
  ///@param  minimumDimensionValue  the dimension value with which to compare the
  ///      element. The order of dimension values from least to greatest is
  ///      <code>{DONTCARE, TRUE, FALSE, 0, 1, 2}</code>.
  void setAtLeastIfValid(int row, int column, int minimumDimensionValue) {
    if (row >= 0 && column >= 0) {
      setAtLeast(row, column, minimumDimensionValue);
    }
  }

  ///  For each element in this <code>IntersectionMatrix</code>, changes the
  ///  element to the corresponding minimum dimension symbol if the element is
  ///  less.
  ///
  ///@param  minimumDimensionSymbols  nine dimension symbols with which to
  ///      compare the elements of this <code>IntersectionMatrix</code>. The
  ///      order of dimension values from least to greatest is <code>{DONTCARE, TRUE, FALSE, 0, 1, 2}</code>
  ///      .
  void setAtLeastDimensionSymbols(String minimumDimensionSymbols) {
    for (int i = 0; i < minimumDimensionSymbols.length; i++) {
      int row = (i / 3).floor();
      int col = i % 3;
      setAtLeast(
          row, col, Dimension.toDimensionValue(minimumDimensionSymbols[i]));
    }
  }

  ///  Changes the elements of this <code>IntersectionMatrix</code> to <code>dimensionValue</code>
  ///  .
  ///
  ///@param  dimensionValue  the dimension value to which to set this <code>IntersectionMatrix</code>
  ///      s elements. Possible values <code>{TRUE, FALSE, DONTCARE, 0, 1, 2}</code>
  ///      .
  void setAll(int dimensionValue) {
    for (int ai = 0; ai < 3; ai++) {
      for (int bi = 0; bi < 3; bi++) {
        matrix[ai][bi] = dimensionValue;
      }
    }
  }

  ///  Returns the value of one of this matrix
  ///  entries.
  ///  The value of the provided index is one of the
  ///  values from the {@link Location} class.
  ///  The value returned is a constant
  ///  from the {@link Dimension} class.
  ///
  ///@param  row     the row of this <code>IntersectionMatrix</code>, indicating
  ///      the interior, boundary or exterior of the first <code>Geometry</code>
  ///@param  column  the column of this <code>IntersectionMatrix</code>,
  ///      indicating the interior, boundary or exterior of the second <code>Geometry</code>
  ///@return         the dimension value at the given matrix position.
  int get(int row, int column) {
    return matrix[row][column];
  }

  ///  Returns <code>true</code> if this <code>IntersectionMatrix</code> is
  ///  FF*FF****.
  ///
  ///@return    <code>true</code> if the two <code>Geometry</code>s related by
  ///      this <code>IntersectionMatrix</code> are disjoint
  bool isDisjoint() {
    return matrix[Location.INTERIOR][Location.INTERIOR] == Dimension.FALSE &&
        matrix[Location.INTERIOR][Location.BOUNDARY] == Dimension.FALSE &&
        matrix[Location.BOUNDARY][Location.INTERIOR] == Dimension.FALSE &&
        matrix[Location.BOUNDARY][Location.BOUNDARY] == Dimension.FALSE;
  }

  ///  Returns <code>true</code> if <code>isDisjoint</code> returns false.
  ///
  ///@return    <code>true</code> if the two <code>Geometry</code>s related by
  ///      this <code>IntersectionMatrix</code> intersect
  bool isIntersects() {
    return !isDisjoint();
  }

  ///  Returns <code>true</code> if this <code>IntersectionMatrix</code> is
  ///  FT*******, F**T***** or F***T****.
  ///
  ///@param  dimensionOfGeometryA  the dimension of the first <code>Geometry</code>
  ///@param  dimensionOfGeometryB  the dimension of the second <code>Geometry</code>
  ///@return                       <code>true</code> if the two <code>Geometry</code>
  ///      s related by this <code>IntersectionMatrix</code> touch; Returns false
  ///      if both <code>Geometry</code>s are points.
  bool isTouches(int dimensionOfGeometryA, int dimensionOfGeometryB) {
    if (dimensionOfGeometryA > dimensionOfGeometryB) {
      //no need to get transpose because pattern matrix is symmetrical
      return isTouches(dimensionOfGeometryB, dimensionOfGeometryA);
    }
    if ((dimensionOfGeometryA == Dimension.A &&
            dimensionOfGeometryB == Dimension.A) ||
        (dimensionOfGeometryA == Dimension.L &&
            dimensionOfGeometryB == Dimension.L) ||
        (dimensionOfGeometryA == Dimension.L &&
            dimensionOfGeometryB == Dimension.A) ||
        (dimensionOfGeometryA == Dimension.P &&
            dimensionOfGeometryB == Dimension.A) ||
        (dimensionOfGeometryA == Dimension.P &&
            dimensionOfGeometryB == Dimension.L)) {
      return matrix[Location.INTERIOR][Location.INTERIOR] == Dimension.FALSE &&
          (isTrue(matrix[Location.INTERIOR][Location.BOUNDARY]) ||
              isTrue(matrix[Location.BOUNDARY][Location.INTERIOR]) ||
              isTrue(matrix[Location.BOUNDARY][Location.BOUNDARY]));
    }
    return false;
  }

  /// Tests whether this geometry crosses the
  /// specified geometry.
  /// <p>
  /// The <code>crosses</code> predicate has the following equivalent definitions:
  /// <ul>
  /// <li>The geometries have some but not all interior points in common.
  /// <li>The DE-9IM Intersection Matrix for the two geometries is
  ///   <ul>
  ///    <li>T*T****** (for P/L, P/A, and L/A situations)
  ///    <li>T*****T** (for L/P, L/A, and A/L situations)
  ///    <li>0******** (for L/L situations)
  ///   </ul>
  /// </ul>
  /// For any other combination of dimensions this predicate returns <code>false</code>.
  /// <p>
  /// The SFS defined this predicate only for P/L, P/A, L/L, and L/A situations.
  /// JTS extends the definition to apply to L/P, A/P and A/L situations as well.
  /// This makes the relation symmetric.
  ///
  ///@param  dimensionOfGeometryA  the dimension of the first <code>Geometry</code>
  ///@param  dimensionOfGeometryB  the dimension of the second <code>Geometry</code>
  ///@return                       <code>true</code> if the two <code>Geometry</code>s
  ///      related by this <code>IntersectionMatrix</code> cross.
  bool isCrosses(int dimensionOfGeometryA, int dimensionOfGeometryB) {
    if ((dimensionOfGeometryA == Dimension.P &&
            dimensionOfGeometryB == Dimension.L) ||
        (dimensionOfGeometryA == Dimension.P &&
            dimensionOfGeometryB == Dimension.A) ||
        (dimensionOfGeometryA == Dimension.L &&
            dimensionOfGeometryB == Dimension.A)) {
      return isTrue(matrix[Location.INTERIOR][Location.INTERIOR]) &&
          isTrue(matrix[Location.INTERIOR][Location.EXTERIOR]);
    }
    if ((dimensionOfGeometryA == Dimension.L &&
            dimensionOfGeometryB == Dimension.P) ||
        (dimensionOfGeometryA == Dimension.A &&
            dimensionOfGeometryB == Dimension.P) ||
        (dimensionOfGeometryA == Dimension.A &&
            dimensionOfGeometryB == Dimension.L)) {
      return isTrue(matrix[Location.INTERIOR][Location.INTERIOR]) &&
          isTrue(matrix[Location.EXTERIOR][Location.INTERIOR]);
    }
    if (dimensionOfGeometryA == Dimension.L &&
        dimensionOfGeometryB == Dimension.L) {
      return matrix[Location.INTERIOR][Location.INTERIOR] == 0;
    }
    return false;
  }

  ///  Tests whether this <code>IntersectionMatrix</code> is
  ///  T*F**F***.
  ///
  ///@return    <code>true</code> if the first <code>Geometry</code> is within
  ///      the second
  bool isWithin() {
    return isTrue(matrix[Location.INTERIOR][Location.INTERIOR]) &&
        matrix[Location.INTERIOR][Location.EXTERIOR] == Dimension.FALSE &&
        matrix[Location.BOUNDARY][Location.EXTERIOR] == Dimension.FALSE;
  }

  ///  Tests whether this <code>IntersectionMatrix</code> is
  ///  T*****FF*.
  ///
  ///@return    <code>true</code> if the first <code>Geometry</code> contains the
  ///      second
  bool isContains() {
    return isTrue(matrix[Location.INTERIOR][Location.INTERIOR]) &&
        matrix[Location.EXTERIOR][Location.INTERIOR] == Dimension.FALSE &&
        matrix[Location.EXTERIOR][Location.BOUNDARY] == Dimension.FALSE;
  }

  ///  Returns <code>true</code> if this <code>IntersectionMatrix</code> is
  ///    <code>T*****FF*</code>
  /// or <code>*T****FF*</code>
  /// or <code>***T**FF*</code>
  /// or <code>****T*FF*</code>
  ///
  ///@return    <code>true</code> if the first <code>Geometry</code> covers the
  ///      second
  bool isCovers() {
    bool hasPointInCommon =
        isTrue(matrix[Location.INTERIOR][Location.INTERIOR]) ||
            isTrue(matrix[Location.INTERIOR][Location.BOUNDARY]) ||
            isTrue(matrix[Location.BOUNDARY][Location.INTERIOR]) ||
            isTrue(matrix[Location.BOUNDARY][Location.BOUNDARY]);

    return hasPointInCommon &&
        matrix[Location.EXTERIOR][Location.INTERIOR] == Dimension.FALSE &&
        matrix[Location.EXTERIOR][Location.BOUNDARY] == Dimension.FALSE;
  }

  ///  Returns <code>true</code> if this <code>IntersectionMatrix</code> is
  ///    <code>T*F**F***</code>
  /// or <code>*TF**F***</code>
  /// or <code>**FT*F***</code>
  /// or <code>**F*TF***</code>
  ///
  ///@return    <code>true</code> if the first <code>Geometry</code>
  /// is covered by the second
  bool isCoveredBy() {
    bool hasPointInCommon =
        isTrue(matrix[Location.INTERIOR][Location.INTERIOR]) ||
            isTrue(matrix[Location.INTERIOR][Location.BOUNDARY]) ||
            isTrue(matrix[Location.BOUNDARY][Location.INTERIOR]) ||
            isTrue(matrix[Location.BOUNDARY][Location.BOUNDARY]);

    return hasPointInCommon &&
        matrix[Location.INTERIOR][Location.EXTERIOR] == Dimension.FALSE &&
        matrix[Location.BOUNDARY][Location.EXTERIOR] == Dimension.FALSE;
  }

  ///  Tests whether the argument dimensions are equal and
  ///  this <code>IntersectionMatrix</code> matches
  ///  the pattern <tt>T*F**FFF*</tt>.
  ///  <p>
  ///  <b>Note:</b> This pattern differs from the one stated in
  ///  <i>Simple feature access - Part 1: Common architecture</i>.
  ///  That document states the pattern as <tt>TFFFTFFFT</tt>.  This would
  ///  specify that
  ///  two identical <tt>POINT</tt>s are not equal, which is not desirable behaviour.
  ///  The pattern used here has been corrected to compute equality in this situation.
  ///
  ///@param  dimensionOfGeometryA  the dimension of the first <code>Geometry</code>
  ///@param  dimensionOfGeometryB  the dimension of the second <code>Geometry</code>
  ///@return                       <code>true</code> if the two <code>Geometry</code>s
  ///      related by this <code>IntersectionMatrix</code> are equal; the
  ///      <code>Geometry</code>s must have the same dimension to be equal
  bool isEquals(int dimensionOfGeometryA, int dimensionOfGeometryB) {
    if (dimensionOfGeometryA != dimensionOfGeometryB) {
      return false;
    }
    return isTrue(matrix[Location.INTERIOR][Location.INTERIOR]) &&
        matrix[Location.INTERIOR][Location.EXTERIOR] == Dimension.FALSE &&
        matrix[Location.BOUNDARY][Location.EXTERIOR] == Dimension.FALSE &&
        matrix[Location.EXTERIOR][Location.INTERIOR] == Dimension.FALSE &&
        matrix[Location.EXTERIOR][Location.BOUNDARY] == Dimension.FALSE;
  }

  ///  Returns <code>true</code> if this <code>IntersectionMatrix</code> is
  ///  <UL>
  ///    <LI> T*T***T** (for two points or two surfaces)
  ///    <LI> 1*T***T** (for two curves)
  ///  </UL>.
  ///
  ///@param  dimensionOfGeometryA  the dimension of the first <code>Geometry</code>
  ///@param  dimensionOfGeometryB  the dimension of the second <code>Geometry</code>
  ///@return                       <code>true</code> if the two <code>Geometry</code>s
  ///      related by this <code>IntersectionMatrix</code> overlap. For this
  ///      function to return <code>true</code>, the <code>Geometry</code>s must
  ///      be two points, two curves or two surfaces.
  bool isOverlaps(int dimensionOfGeometryA, int dimensionOfGeometryB) {
    if ((dimensionOfGeometryA == Dimension.P &&
            dimensionOfGeometryB == Dimension.P) ||
        (dimensionOfGeometryA == Dimension.A &&
            dimensionOfGeometryB == Dimension.A)) {
      return isTrue(matrix[Location.INTERIOR][Location.INTERIOR]) &&
          isTrue(matrix[Location.INTERIOR][Location.EXTERIOR]) &&
          isTrue(matrix[Location.EXTERIOR][Location.INTERIOR]);
    }
    if (dimensionOfGeometryA == Dimension.L &&
        dimensionOfGeometryB == Dimension.L) {
      return matrix[Location.INTERIOR][Location.INTERIOR] == 1 &&
          isTrue(matrix[Location.INTERIOR][Location.EXTERIOR]) &&
          isTrue(matrix[Location.EXTERIOR][Location.INTERIOR]);
    }
    return false;
  }

  ///  Returns whether the elements of this <code>IntersectionMatrix</code>
  ///  satisfies the required dimension symbols.
  ///
  ///@param  requiredDimensionSymbols  nine dimension symbols with which to
  ///      compare the elements of this <code>IntersectionMatrix</code>. Possible
  ///      values are <code>{T, F, * , 0, 1, 2}</code>.
  ///@return                           <code>true</code> if this <code>IntersectionMatrix</code>
  ///      matches the required dimension symbols
  bool matches(String requiredDimensionSymbols) {
    if (requiredDimensionSymbols.length != 9) {
      throw ArgumentError("Should be length 9: $requiredDimensionSymbols");
    }
    for (int ai = 0; ai < 3; ai++) {
      for (int bi = 0; bi < 3; bi++) {
        if (!matchesDimValue(
            matrix[ai][bi], requiredDimensionSymbols[3 * ai + bi])) {
          return false;
        }
      }
    }
    return true;
  }

  ///  Transposes this IntersectionMatrix.
  ///
  ///@return    this <code>IntersectionMatrix</code> as a convenience
  IntersectionMatrix transpose() {
    int temp = matrix[1][0];
    matrix[1][0] = matrix[0][1];
    matrix[0][1] = temp;
    temp = matrix[2][0];
    matrix[2][0] = matrix[0][2];
    matrix[0][2] = temp;
    temp = matrix[2][1];
    matrix[2][1] = matrix[1][2];
    matrix[1][2] = temp;
    return this;
  }

  ///  Returns a nine-character <code>String</code> representation of this <code>IntersectionMatrix</code>
  ///  .
  ///
  ///@return    the nine dimension symbols of this <code>IntersectionMatrix</code>
  ///      in row-major order.
  String toString() {
    String base = "123456789";
    for (int ai = 0; ai < 3; ai++) {
      for (int bi = 0; bi < 3; bi++) {
        base = StringUtils.replaceCharAt(
            base, 3 * ai + bi, Dimension.toDimensionSymbol(matrix[ai][bi]));
      }
    }
    return base;
  }
}

///// Coordinate subclass supporting XY ordinate.
///// <p>
///// This data object is suitable for use with coordinate sequences dimension 3, measures 1.
///// The {@link Coordinate#Z} field is visible, but intended to be ignored.
/////
///// @since 1.16
//class CoordinateXY extends Coordinate {
//  CoordinateXY.empty() : super.empty();
//
//  CoordinateXY(double x, double y) : super(x, y);
//
//  CoordinateXY.fromCoordinate(Coordinate coord) : super.fromCoordinate(coord);
//
//  CoordinateXY.fromCoordinateXY(CoordinateXY coord) : super(coord.x, coord.y);
//
//  CoordinateXY copy() {
//    return CoordinateXY.fromCoordinateXY(this);
//  }
//
//  num get x => _x;
//
//  set x(newX) => _x = newX;
//
//  num get y => _y;
//
//  set y(newY) => _y = newY;
//
//  /// The z-ordinate is not supported
//  num get z => Coordinate.NULL_ORDINATE;
//
//  /// The z-ordinate is not supported
//  set z(newZ) => ArgumentError("CoordinateXY dimension 2 does not support z-ordinate");
//
//  double getOrdinate(int ordinateIndex) {
//    switch (ordinateIndex) {
//      case Coordinate.X:
//        return _x;
//      case Coordinate.Y:
//        return _y;
//    }
//    throw ArgumentError("Invalid ordinate index: $ordinateIndex");
//  }
//
//  void setOrdinate(int ordinateIndex, double value) {
//    switch (ordinateIndex) {
//      case Coordinate.X:
//        _x = value;
//        break;
//      case Coordinate.Y:
//        _y = value;
//        break;
//      default:
//        throw ArgumentError("Invalid ordinate index: $ordinateIndex");
//    }
//  }
//
//  String toString() {
//    return "($_x, $_y)";
//  }
//}

/// The internal representation of a list of coordinates inside a Geometry.
/// <p>
/// This allows Geometries to store their
/// points using something other than the JTS {@link Coordinate} class.
/// For example, a storage-efficient implementation
/// might store coordinate sequences as an array of x's
/// and an array of y's.
/// Or a custom coordinate class might support extra attributes like M-values.
/// <p>
/// Implementing a custom coordinate storage structure
/// requires implementing the {@link CoordinateSequence} and
/// {@link CoordinateSequenceFactory} interfaces.
/// To use the custom CoordinateSequence, create a
/// new {@link GeometryFactory} parameterized by the CoordinateSequenceFactory
/// The {@link GeometryFactory} can then be used to create new {@link Geometry}s.
/// The new Geometries
/// will use the custom CoordinateSequence implementation.
/// <p>
/// For an example, see the code for
/// {@link ExtendedCoordinateExample}.
///
/// @see CoordinateArraySequenceFactory
/// @see PackedCoordinateSequenceFactory
/// @see ExtendedCoordinateExample
///
/// @version 1.7
abstract class CoordinateSequence {
  /// Standard ordinate index value for, where X is 0 */
  static const int X = 0;

  /// Standard ordinate index value for, where Y is 1 */
  static const int Y = 1;

  /// Standard ordinate index value for, where Z is 2.
  ///
  /// <p>This constant assumes XYZM coordinate sequence definition, please check this assumption
  /// using {@link #getDimension()} and {@link #getMeasures()} before use.
  static const int Z = 2;

  /// Standard ordinate index value for, where M is 3.
  ///
  /// <p>This constant assumes XYZM coordinate sequence definition, please check this assumption
  /// using {@link #getDimension()} and {@link #getMeasures()} before use.
  static const int M = 3;

  /// Returns the dimension (number of ordinates in each coordinate) for this sequence.
  ///
  /// <p>This total includes any measures, indicated by non-zero {@link #getMeasures()}.
  ///
  /// @return the dimension of the sequence.
  int getDimension();

  /// Returns the number of measures included in {@link #getDimension()} for each coordinate for this
  /// sequence.
  ///
  /// For a measured coordinate sequence a non-zero value is returned.
  /// <ul>
  /// <li>For XY sequence measures is zero</li>
  /// <li>For XYM sequence measure is one<li>
  /// <li>For XYZ sequence measure is zero</li>
  /// <li>For XYZM sequence measure is one</li>
  /// <li>Values greater than one are supported</li>
  /// </ul>
  ///
  /// @return the number of measures included in dimension
  int getMeasures() {
    return 0;
  }

  /// Checks {@link #getDimension()} and {@link #getMeasures()} to determine if {@link #getZ(int)}
  /// is supported.
  ///
  /// @return true if {@link #getZ(int)} is supported.
  bool hasZ() {
    return (getDimension() - getMeasures()) > 2;
  }

  /// Checks {@link #getMeasures()} to determine if {@link #getM(int)}
  /// is supported.
  ///
  /// @return true if {@link #getM(int)} is supported.
  bool hasM() {
    return getDimension() > 2 && getMeasures() > 0;
  }

  /// Creates a coordinate for use in this sequence.
  /// <p>
  /// The coordinate is created supporting the same number of {@link #getDimension()} and {@link #getMeasures()}
  /// as this sequence and is suitable for use with {@link #getCoordinate(int, Coordinate)}.
  /// </p>
  /// @return coordinate for use with this sequence
  Coordinate createCoordinate() {
    return Coordinates.createWithMeasure(getDimension(), getMeasures());
  }

  /// Returns (possibly a copy of) the i'th coordinate in this sequence.
  /// Whether or not the Coordinate returned is the actual underlying
  /// Coordinate or merely a copy depends on the implementation.
  /// <p>
  /// Note that in the future the semantics of this method may change
  /// to guarantee that the Coordinate returned is always a copy.
  /// Callers should not to assume that they can modify a CoordinateSequence by
  /// modifying the object returned by this method.
  ///
  /// @param i the index of the coordinate to retrieve
  /// @return the i'th coordinate in the sequence
  Coordinate getCoordinate(int i);

  /// Returns a copy of the i'th coordinate in this sequence.
  /// This method optimizes the situation where the caller is
  /// going to make a copy anyway - if the implementation
  /// has already created a new Coordinate object, no further copy is needed.
  ///
  /// @param i the index of the coordinate to retrieve
  /// @return a copy of the i'th coordinate in the sequence
  Coordinate getCoordinateCopy(int i);

  /// Copies the i'th coordinate in the sequence to the supplied
  /// {@link Coordinate}.  Only the first two dimensions are copied.
  ///
  /// @param index the index of the coordinate to copy
  /// @param coord a {@link Coordinate} to receive the value
  void getCoordinateInto(int index, Coordinate coord);

  /// Returns ordinate X (0) of the specified coordinate.
  ///
  /// @param index
  /// @return the value of the X ordinate in the index'th coordinate
  double getX(int index);

  /// Returns ordinate Y (1) of the specified coordinate.
  ///
  /// @param index
  /// @return the value of the Y ordinate in the index'th coordinate
  double getY(int index);

  /// Returns ordinate Z of the specified coordinate if available.
  ///
  /// @param index
  /// @return the value of the Z ordinate in the index'th coordinate, or Double.NaN if not defined.
  double getZ(int index) {
    if (hasZ()) {
      return getOrdinate(index, 2);
    } else {
      return double.nan;
    }
  }

  /// Returns ordinate M of the specified coordinate if available.
  ///
  /// @param index
  /// @return the value of the M ordinate in the index'th coordinate, or Double.NaN if not defined.
  double getM(int index) {
    if (hasM()) {
      final int mIndex = getDimension() - getMeasures();
      return getOrdinate(index, mIndex);
    } else {
      return double.nan;
    }
  }

  /// Returns the ordinate of a coordinate in this sequence.
  /// Ordinate indices 0 and 1 are assumed to be X and Y.
  /// <p>
  /// Ordinates indices greater than 1 have user-defined semantics
  /// (for instance, they may contain other dimensions or measure
  /// values as described by {@link #getDimension()} and {@link #getMeasures()}).
  ///
  /// @param index  the coordinate index in the sequence
  /// @param ordinateIndex the ordinate index in the coordinate (in range [0, dimension-1])
  double getOrdinate(int index, int ordinateIndex);

  /// Returns the number of coordinates in this sequence.
  /// @return the size of the sequence
  int size();

  /// Sets the value for a given ordinate of a coordinate in this sequence.
  ///
  /// @param index  the coordinate index in the sequence
  /// @param ordinateIndex the ordinate index in the coordinate (in range [0, dimension-1])
  /// @param value  the new ordinate value
  void setOrdinate(int index, int ordinateIndex, double value);

  /// Returns (possibly copies of) the Coordinates in this collection.
  /// Whether or not the Coordinates returned are the actual underlying
  /// Coordinates or merely copies depends on the implementation. Note that
  /// if this implementation does not store its data as an array of Coordinates,
  /// this method will incur a performance penalty because the array needs to
  /// be built from scratch.
  ///
  /// @return a array of coordinates containing the point values in this sequence
  List<Coordinate> toCoordinateArray();

  /// Expands the given {@link Envelope} to include the coordinates in the sequence.
  /// Allows implementing classes to optimize access to coordinate values.
  ///
  /// @param env the envelope to expand
  /// @return a ref to the expanded envelope
  Envelope expandEnvelope(Envelope env);

  /// Returns a deep copy of this collection.
  /// Called by Geometry#clone.
  ///
  /// @return a copy of the coordinate sequence containing copies of all points
  /// @deprecated Recommend {@link #copy()}
  Object clone();

  /// Returns a deep copy of this collection.
  ///
  /// @return a copy of the coordinate sequence containing copies of all points
  CoordinateSequence copy();
}

/// A factory to create concrete instances of {@link CoordinateSequence}s.
/// Used to configure {@link GeometryFactory}s
/// to provide specific kinds of CoordinateSequences.
///
/// @version 1.7
abstract class CoordinateSequenceFactory {
  /// Returns a {@link CoordinateSequence} based on the given array.
  /// Whether the array is copied or simply referenced
  /// is implementation-dependent.
  /// This method must handle null arguments by creating an empty sequence.
  ///
  /// @param coordinates the coordinates
  CoordinateSequence create(List<Coordinate> coordinates);

  /// Creates a {@link CoordinateSequence} which is a copy
  /// of the given {@link CoordinateSequence}.
  /// This method must handle null arguments by creating an empty sequence.
  ///
  /// @param coordSeq the coordinate sequence to copy
  CoordinateSequence createFromSequence(CoordinateSequence coordSeq);

  /// Creates a {@link CoordinateSequence} of the specified size and dimension.
  /// For this to be useful, the {@link CoordinateSequence} implementation must
  /// be mutable.
  /// <p>
  /// If the requested dimension is larger than the CoordinateSequence implementation
  /// can provide, then a sequence of maximum possible dimension should be created.
  /// An error should not be thrown.
  ///
  /// @param size the number of coordinates in the sequence
  /// @param dimension the dimension of the coordinates in the sequence (if user-specifiable,
  /// otherwise ignored)
  CoordinateSequence createSizeDim(int size, int dimension);

  /// Creates a {@link CoordinateSequence} of the specified size and dimension with measure support.
  /// For this to be useful, the {@link CoordinateSequence} implementation must
  /// be mutable.
  /// <p>
  /// If the requested dimension or measures are larger than the CoordinateSequence implementation
  /// can provide, then a sequence of maximum possible dimension should be created.
  /// An error should not be thrown.
  ///
  /// @param size the number of coordinates in the sequence
  /// @param dimension the dimension of the coordinates in the sequence (if user-specifiable,
  /// otherwise ignored)
  /// @param measures the number of measures of the coordinates in the sequence (if user-specifiable,
  /// otherwise ignored)
  CoordinateSequence createSizeDimMeas(int size, int dimension, int measures) {
    return createSizeDim(size, dimension);
  }
}

/// The types of Precision Model which JTS supports.
class Type {
  static Map<String, Type> nameToTypeMap = {};
  String name;

  Type(this.name) {
    nameToTypeMap[name] = this;
  }

  String toString() {
    return name;
  }

  ///see http://www.javaworld.com/javaworld/javatips/jw-javatip122.html
  Object readResolve() {
    return nameToTypeMap[name] as Object;
  }
}

/// Specifies the precision model of the {@link Coordinate}s in a {@link Geometry}.
/// In other words, specifies the grid of allowable
///  points for all <code>Geometry</code>s.
/// <p>
/// The {@link #makePrecise(Coordinate)} method allows rounding a coordinate to
/// a "precise" value; that is, one whose
///  precision is known exactly.
///<p>
/// Coordinates are assumed to be precise in geometries.
/// That is, the coordinates are assumed to be rounded to the
/// precision model given for the geometry.
/// JTS input routines automatically round coordinates to the precision model
/// before creating Geometries.
/// All internal operations
/// assume that coordinates are rounded to the precision model.
/// Constructive methods (such as bool operations) always round computed
/// coordinates to the appropriate precision model.
/// <p>
/// Currently three types of precision model are supported:
/// <ul>
/// <li>FLOATING - represents full double precision floating point.
/// This is the default precision model used in JTS
/// <li>FLOATING_SINGLE - represents single precision floating point.
/// <li>FIXED - represents a model with a fixed number of decimal places.
///  A Fixed Precision Model is specified by a scale factor.
///  The scale factor specifies the size of the grid which numbers are rounded to.
///  Input coordinates are mapped to fixed coordinates according to the following
///  equations:
///    <UL>
///      <LI> jtsPt.x = round( (inputPt.x * scale ) / scale
///      <LI> jtsPt.y = round( (inputPt.y * scale ) / scale
///    </UL>
/// </ul>
/// For example, to specify 3 decimal places of precision, use a scale factor
/// of 1000. To specify -3 decimal places of precision (i.e. rounding to
/// the nearest 1000), use a scale factor of 0.001.
/// <p>
/// Coordinates are represented internally as Java double-precision values.
/// Since Java uses the IEEE-394 floating point standard, this
/// provides 53 bits of precision. (Thus the maximum precisely representable
/// <i>integer</i> is 9,007,199,254,740,992 - or almost 16 decimal digits of precision).
/// <p>
/// JTS binary methods currently do not handle inputs which have different precision models.
/// The precision model of any constructed geometric value is undefined.
///
///@version 1.7
class PrecisionModel implements Comparable {
  /// Determines which of two {@link PrecisionModel}s is the most precise
  /// (allows the greatest number of significant digits).
  ///
  /// @param pm1 a PrecisionModel
  /// @param pm2 a PrecisionModel
  /// @return the PrecisionModel which is most precise
  static PrecisionModel mostPrecise(PrecisionModel pm1, PrecisionModel pm2) {
    if (pm1.compareTo(pm2) >= 0) return pm1;
    return pm2;
  }

  /// Fixed Precision indicates that coordinates have a fixed number of decimal places.
  /// The number of decimal places is determined by the log10 of the scale factor.
  static final Type FIXED = Type("FIXED");

  /// Floating precision corresponds to the standard Java
  /// double-precision floating-point representation, which is
  /// based on the IEEE-754 standard
  static final Type FLOATING = Type("FLOATING");

  ///  The maximum precise value representable in a double. Since IEE754
  ///  double-precision numbers allow 53 bits of mantissa, the value is equal to
  ///  2^53 - 1.  This provides <i>almost</i> 16 decimal digits of precision.
  static final double maximumPreciseValue = 9007199254740992.0;

  /// The type of PrecisionModel this represents.
  late Type modelType;

  /// The scale factor which determines the number of decimal places in fixed precision.
  double scale = 0.0;

  /// Creates a <code>PrecisionModel</code> with a default precision
  /// of FLOATING.
  PrecisionModel() {
    // default is floating precision
    modelType = FLOATING;
  }

  /// Creates a <code>PrecisionModel</code> that specifies
  /// an explicit precision model type.
  /// If the model type is FIXED the scale factor will default to 1.
  ///
  /// @param modelType the type of the precision model
  PrecisionModel.fromType(Type modelType) {
    this.modelType = modelType;
    if (modelType == FIXED) {
      setScale(1.0);
    }
  }

  ///  Creates a <code>PrecisionModel</code> that specifies Fixed precision.
  ///  Fixed-precision coordinates are represented as precise internal coordinates,
  ///  which are rounded to the grid defined by the scale factor.
  ///
  ///@param  scale    amount by which to multiply a coordinate after subtracting
  ///      the offset, to obtain a precise coordinate
  PrecisionModel.fixedPrecision(double scale) {
    modelType = FIXED;
    setScale(scale);
  }

  ///  Copy constructor to create a new <code>PrecisionModel</code>
  ///  from an existing one.
  PrecisionModel.fromPrecisionModel(PrecisionModel pm) {
    modelType = pm.modelType;
    scale = pm.scale;
  }

  /// Tests whether the precision model supports floating point
  /// @return <code>true</code> if the precision model supports floating point
  bool isFloating() {
    return modelType == FLOATING;
  }

  /// Returns the maximum number of significant digits provided by this
  /// precision model.
  /// Intended for use by routines which need to print out
  /// decimal representations of precise values (such as {@link WKTWriter}).
  /// <p>
  /// This method would be more correctly called
  /// <tt>getMinimumDecimalPlaces</tt>,
  /// since it actually computes the number of decimal places
  /// that is required to correctly display the full
  /// precision of an ordinate value.
  /// <p>
  /// Since it is difficult to compute the required number of
  /// decimal places for scale factors which are not powers of 10,
  /// the algorithm uses a very rough approximation in this case.
  /// This has the side effect that for scale factors which are
  /// powers of 10 the value returned is 1 greater than the true value.
  ///
  ///
  /// @return the maximum number of decimal places provided by this precision model
  int getMaximumSignificantDigits() {
    int maxSigDigits = 16;
    if (modelType == FLOATING) {
      maxSigDigits = 16;
    } else if (modelType == FIXED) {
      double n = math.log(getScale()) / math.log(10);
      maxSigDigits = 1 + n.ceil();
    }
    return maxSigDigits;
  }

  /// Returns the scale factor used to specify a fixed precision model.
  /// The number of decimal places of precision is
  /// equal to the base-10 logarithm of the scale factor.
  /// Non-integral and negative scale factors are supported.
  /// Negative scale factors indicate that the places
  /// of precision is to the left of the decimal point.
  ///
  ///@return the scale factor for the fixed precision model
  double getScale() {
    return scale;
  }

  /// Gets the type of this precision model
  /// @return the type of this precision model
  /// @see Type
  Type getType() {
    return modelType;
  }

  ///  Sets the multiplying factor used to obtain a precise coordinate.
  /// This method is  because PrecisionModel is an immutable (value) type.
  void setScale(double scale) {
    this.scale = scale.abs();
  }

  /// Rounds a numeric value to the PrecisionModel grid.
  /// Asymmetric Arithmetic Rounding is used, to provide
  /// uniform rounding behaviour no matter where the number is
  /// on the number line.
  /// <p>
  /// This method has no effect on NaN values.
  /// <p>
  /// <b>Note:</b> Java's <code>Math#rint</code> uses the "Banker's Rounding" algorithm,
  /// which is not suitable for precision operations elsewhere in JTS.
  double makePrecise(double val) {
    // don't change NaN values
    if (val.isNaN) return val;

    if (modelType == FIXED) {
      return (val * scale).round() / scale;
    }
    // modelType == FLOATING - no rounding necessary
    return val;
  }

  /// Rounds a Coordinate to the PrecisionModel grid.
  void makeCoordinatePrecise(Coordinate coord) {
    // optimization for full precision
    if (modelType == FLOATING) return;

    coord.x = makePrecise(coord.x);
    coord.y = makePrecise(coord.y);
    //MD says it's OK that we're not makePrecise'ing the z [Jon Aquino]
  }

  String toString() {
    String description = "UNKNOWN";
    if (modelType == FLOATING) {
      description = "Floating";
    } else if (modelType == FIXED) {
      description = "Fixed (Scale=${getScale()})";
    }
    return description;
  }

  bool equals(Object other) {
    if (!(other is PrecisionModel)) {
      return false;
    }
    PrecisionModel otherPrecisionModel = other;
    return modelType == otherPrecisionModel.modelType &&
        scale == otherPrecisionModel.scale;
  }

  ///  Compares this {@link PrecisionModel} object with the specified object for order.
  /// A PrecisionModel is greater than another if it provides greater precision.
  /// The comparison is based on the value returned by the
  /// {@link #getMaximumSignificantDigits} method.
  /// This comparison is not strictly accurate when comparing floating precision models
  /// to fixed models; however, it is correct when both models are either floating or fixed.
  ///
  ///@param  o  the <code>PrecisionModel</code> with which this <code>PrecisionModel</code>
  ///      is being compared
  ///@return    a negative integer, zero, or a positive integer as this <code>PrecisionModel</code>
  ///      is less than, equal to, or greater than the specified <code>PrecisionModel</code>
  int compareTo(dynamic o) {
    PrecisionModel other = o;

    int sigDigits = getMaximumSignificantDigits();
    int otherSigDigits = other.getMaximumSignificantDigits();
    return sigDigits.compareTo(otherSigDigits);
//    if (sigDigits > otherSigDigits)
//      return 1;
//    else if
//    if (modelType == FLOATING && other.modelType == FLOATING) return 0;
//    if (modelType == FLOATING && other.modelType != FLOATING) return 1;
//    if (modelType != FLOATING && other.modelType == FLOATING) return -1;
//    if (modelType == FIXED && other.modelType == FIXED) {
//      if (scale > other.scale)
//        return 1;
//      else if (scale < other.scale)
//        return -1;
//      else
//        return 0;
//    }
//    Assert.shouldNeverReachHere("Unknown Precision Model type encountered");
//    return 0;
  }
}

///  <code>Geometry</code> classes support the concept of applying
///  a <code>GeometryComponentFilter</code>
///  filter to the <code>Geometry</code>.
///  The filter is applied to every component of the <code>Geometry</code>
///  which is itself a <code>Geometry</code>
///  and which does not itself contain any components.
/// (For instance, all the {@link LinearRing}s in {@link Polygon}s are visited,
/// but in a {@link MultiPolygon} the {@link Polygon}s themselves are not visited.)
/// Thus the only classes of Geometry which must be
/// handled as arguments to {@link #filter}
/// are {@link LineString}s, {@link LinearRing}s and {@link Point}s.
///  <p>
///  A <code>GeometryComponentFilter</code> filter can either
///  record information about the <code>Geometry</code>
///  or change the <code>Geometry</code> in some way.
///  <code>GeometryComponentFilter</code>
///  is an example of the Gang-of-Four Visitor pattern.
///
///@version 1.7
abstract class GeometryComponentFilter {
  ///  Performs an operation with or on <code>geom</code>.
  ///
  ///@param  geom  a <code>Geometry</code> to which the filter is applied.
  void filter(Geometry geom);
}

class GeometryChangedFilter implements GeometryComponentFilter {
  void filter(Geometry geom) {
    geom.geometryChangedAction();
  }
}

/**
 * Extracts all the 1-dimensional ({@link LineString}) components from a {@link Geometry}.
 * For polygonal geometries, this will extract all the component {@link LinearRing}s.
 * If desired, <code>LinearRing</code>s can be forced to be returned as <code>LineString</code>s.
 *
 * @version 1.7
 */
class LinearComponentExtracter implements GeometryComponentFilter {
  /**
   * Extracts the linear components from a single {@link Geometry}
   * and adds them to the provided {@link Collection}.
   *
   * @param geoms the collection of geometries from which to extract linear components
   * @param lines the collection to add the extracted linear components to
   * @return the collection of linear components (LineStrings or LinearRings)
   */
  static List getLinesLL(List geoms, List<LineString> lines) {
    for (Iterator i = geoms.iterator; i.moveNext();) {
      Geometry g = i.current as Geometry;
      getLinesGL(g, lines);
    }
    return lines;
  }

  /**
   * Extracts the linear components from a single {@link Geometry}
   * and adds them to the provided {@link List}.
   *
   * @param geoms the List of geometries from which to extract linear components
   * @param lines the collection to add the extracted linear components to
   * @param forceToLineString true if LinearRings should be converted to LineStrings
   * @return the collection of linear components (LineStrings or LinearRings)
   */
  static List getLinesLLF(
      List geoms, List<LineString> lines, bool forceToLineString) {
    for (Iterator i = geoms.iterator; i.moveNext();) {
      Geometry g = i.current as Geometry;
      getLinesGLF(g, lines, forceToLineString);
    }
    return lines;
  }

  /**
   * Extracts the linear components from a single {@link Geometry}
   * and adds them to the provided {@link List}.
   *
   * @param geom the geometry from which to extract linear components
   * @param lines the List to add the extracted linear components to
   * @return the List of linear components (LineStrings or LinearRings)
   */
  static List getLinesGL(Geometry geom, List<LineString> lines) {
    if (geom is LineString) {
      lines.add(geom);
    } else {
      geom.applyGCF(new LinearComponentExtracter(lines));
    }
    return lines;
  }

  /**
   * Extracts the linear components from a single {@link Geometry}
   * and adds them to the provided {@link List}.
   *
   * @param geom the geometry from which to extract linear components
   * @param lines the List to add the extracted linear components to
   * @param forceToLineString true if LinearRings should be converted to LineStrings
   * @return the List of linear components (LineStrings or LinearRings)
   */
  static List getLinesGLF(
      Geometry geom, List<LineString> lines, bool forceToLineString) {
    geom.applyGCF(
        new LinearComponentExtracter.withForced(lines, forceToLineString));
    return lines;
  }

  /**
   * Extracts the linear components from a single geometry.
   * If more than one geometry is to be processed, it is more
   * efficient to create a single {@link LinearComponentExtracter} instance
   * and pass it to multiple geometries.
   *
   * @param geom the geometry from which to extract linear components
   * @return the list of linear components
   */
  static List<Geometry> getLines(Geometry geom) {
    return getLinesGF(geom, false);
  }

  /**
   * Extracts the linear components from a single geometry.
   * If more than one geometry is to be processed, it is more
   * efficient to create a single {@link LinearComponentExtracter} instance
   * and pass it to multiple geometries.
   *
   * @param geom the geometry from which to extract linear components
   * @param forceToLineString true if LinearRings should be converted to LineStrings
   * @return the list of linear components
   */
  static List<Geometry> getLinesGF(Geometry geom, bool forceToLineString) {
    List<LineString> lines = [];
    geom.applyGCF(
        new LinearComponentExtracter.withForced(lines, forceToLineString));
    return lines;
  }

  /**
   * Extracts the linear components from a single {@link Geometry}
   * and returns them as either a {@link LineString} or {@link MultiLineString}.
   *
   * @param geom the geometry from which to extract
   * @return a linear geometry
   */
  static Geometry getGeometry(Geometry geom) {
    return geom.getFactory().buildGeometry(getLines(geom));
  }

  /**
   * Extracts the linear components from a single {@link Geometry}
   * and returns them as either a {@link LineString} or {@link MultiLineString}.
   *
   * @param geom the geometry from which to extract
   * @param forceToLineString true if LinearRings should be converted to LineStrings
   * @return a linear geometry
   */
  static Geometry getGeometryWithForce(Geometry geom, bool forceToLineString) {
    return geom.getFactory().buildGeometry(getLinesGF(geom, forceToLineString));
  }

  List<LineString> lines;
  bool isForcedToLineString = false;

  /**
   * Constructs a LineExtracterFilter with a list in which to store LineStrings found.
   */
  LinearComponentExtracter(this.lines);

  /**
   * Constructs a LineExtracterFilter with a list in which to store LineStrings found.
   */
  LinearComponentExtracter.withForced(this.lines, bool isForcedToLineString) {
    this.isForcedToLineString = isForcedToLineString;
  }

  /**
   * Indicates that LinearRing components should be
   * converted to pure LineStrings.
   *
   * @param isForcedToLineString true if LinearRings should be converted to LineStrings
   */
  void setForceToLineString(bool isForcedToLineString) {
    this.isForcedToLineString = isForcedToLineString;
  }

  void filter(Geometry geom) {
    if (isForcedToLineString && geom is LinearRing) {
      LineString line =
          geom.getFactory().createLineStringSeq(geom.getCoordinateSequence());
      lines.add(line);
      return;
    }
    // if not being forced, and this is a linear component
    if (geom is LineString) lines.add(geom);

    // else this is not a linear component, so skip it
  }
}

/// Extracts all the 0-dimensional ({@link Point}) components from a {@link Geometry}.
///
/// @version 1.7
/// @see GeometryExtracter
class PointExtracter implements GeometryFilter {
  /// Extracts the {@link Point} elements from a single {@link Geometry}
  /// and adds them to the provided {@link List}.
  ///
  /// @param geom the geometry from which to extract
  /// @param list the list to add the extracted elements to
  static List<Point> getPointsWithList(Geometry geom, List<Point> list) {
    if (geom is Point) {
      list.add(geom);
    } else if (geom is GeometryCollection) {
      geom.applyGF(PointExtracter(list));
    }
    // skip non-Polygonal elemental geometries

    return list;
  }

  /// Extracts the {@link Point} elements from a single {@link Geometry}
  /// and returns them in a {@link List}.
  ///
  /// @param geom the geometry from which to extract
  static List getPoints(Geometry geom) {
    if (geom is Point) {
      return [geom];
    }
    return getPointsWithList(geom, <Point>[]);
  }

  List pts;

  /// Constructs a PointExtracterFilter with a list in which to store Points found.
  PointExtracter(this.pts);

  void filter(Geometry geom) {
    if (geom is Point) pts.add(geom);
  }
}

/// Extracts all the {@link Polygon} elements from a {@link Geometry}.
///
/// @version 1.7
/// @see GeometryExtracter
class PolygonExtracter implements GeometryFilter {
  /// Extracts the {@link Polygon} elements from a single {@link Geometry}
  /// and adds them to the provided {@link List}.
  ///
  /// @param geom the geometry from which to extract
  /// @param list the list to add the extracted elements to
  static List<Polygon> getPolygonsWithList(Geometry geom, List<Polygon> list) {
    if (geom is Polygon) {
      list.add(geom);
    } else if (geom is GeometryCollection) {
      geom.applyGF(PolygonExtracter(list));
    }
    // skip non-Polygonal elemental geometries

    return list;
  }

  /// Extracts the {@link Polygon} elements from a single {@link Geometry}
  /// and returns them in a {@link List}.
  ///
  /// @param geom the geometry from which to extract
  static List getPolygons(Geometry geom) {
    return getPolygonsWithList(geom, <Polygon>[]);
  }

  List comps;

  /// Constructs a PolygonExtracterFilter with a list in which to store Polygons found.
  PolygonExtracter(this.comps);

  void filter(Geometry geom) {
    if (geom is Polygon) comps.add(geom);
  }
}

/**
 * A class which supports creating new {@link Geometry}s
 * which are modifications of existing ones,
 * maintaining the same type structure.
 * Geometry objects are intended to be treated as immutable.
 * This class "modifies" Geometrys
 * by traversing them, applying a user-defined
 * {@link GeometryEditorOperation}, {@link CoordinateSequenceOperation} or {@link CoordinateOperation}
 * and creating new Geometrys with the same structure but
 * (possibly) modified components.
 * <p>
 * Examples of the kinds of modifications which can be made are:
 * <ul>
 * <li>the values of the coordinates may be changed.
 *     The editor does not check whether changing coordinate values makes the result Geometry invalid
 * <li>the coordinate lists may be changed
 *     (e.g. by adding, deleting or modifying coordinates).
 *     The modified coordinate lists must be consistent with their original parent component
 *     (e.g. a <tt>LinearRing</tt> must always have at least 4 coordinates, and the first and last
 *     coordinate must be equal)
 * <li>components of the original geometry may be deleted
 *    (e.g. holes may be removed from a Polygon, or LineStrings removed from a MultiLineString).
 *     Deletions will be propagated up the component tree appropriately.
 * </ul>
 * All changes must be consistent with the original Geometry's structure
 * (e.g. a <tt>Polygon</tt> cannot be collapsed into a <tt>LineString</tt>).
 * If changing the structure is required, use a {@link GeometryTransformer}.
 * <p>
 * This class supports creating an edited Geometry
 * using a different <code>GeometryFactory</code> via the {@link #GeometryEditor(GeometryFactory)}
 * constructor.
 * Examples of situations where this is required is if the geometry is
 * transformed to a new SRID and/or a new PrecisionModel.
 * <p>
 * <b>Usage Notes</b>
 * <ul>
 * <li>The resulting Geometry is not checked for validity.
 * If validity needs to be enforced, the new Geometry's
 * {@link Geometry#isValid} method should be called.
 * <li>By default the UserData of the input geometry is not copied to the result.
 * </ul>
 *
 * @see GeometryTransformer
 * @see Geometry#isValid
 *
 * @version 1.7
 */
class GeometryEditor {
  /**
   * The factory used to create the modified Geometry.
   * If <tt>null</tt> the GeometryFactory of the input is used.
   */
  GeometryFactory? _geomFactory = null;
  bool isUserDataCopied = false;

  /**
   * Creates a new GeometryEditor object which will create
   * edited {@link Geometry}s with the same {@link GeometryFactory} as the input Geometry.
   */
  GeometryEditor.empty() {}

  /**
   * Creates a new GeometryEditor object which will create
   * edited {@link Geometry}s with the given {@link GeometryFactory}.
   *
   * @param factory the GeometryFactory to create  edited Geometrys with
   */
  GeometryEditor(GeometryFactory factory) {
    this._geomFactory = factory;
  }

  /**
   * Sets whether the User Data is copied to the edit result.
   * Only the object reference is copied.
   *
   * @param isUserDataCopied true if the input user data should be copied.
   */
  void setCopyUserData(bool isUserDataCopied) {
    this.isUserDataCopied = isUserDataCopied;
  }

  /**
   * Edit the input {@link Geometry} with the given edit operation.
   * Clients can create subclasses of {@link GeometryEditorOperation} or
   * {@link CoordinateOperation} to perform required modifications.
   *
   * @param geometry the Geometry to edit
   * @param operation the edit operation to carry out
   * @return a new {@link Geometry} which is the result of the editing (which may be empty)
   */
  Geometry? edit(Geometry? geometry, GeometryEditorOperation operation) {
    // nothing to do
    if (geometry == null) return null;

    Geometry? result = editInternal(geometry, operation);
    if (isUserDataCopied) {
      result?.setUserData(geometry.getUserData());
    }
    return result;
  }

  Geometry? editInternal(Geometry geometry, GeometryEditorOperation operation) {
    // if client did not supply a GeometryFactory, use the one from the input Geometry
    if (_geomFactory == null) _geomFactory = geometry.getFactory();

    if (geometry is GeometryCollection) {
      return editGeometryCollection(geometry, operation);
    }

    if (geometry is Polygon) {
      return editPolygon(geometry, operation);
    }

    if (geometry is Point) {
      return operation.edit(geometry, _geomFactory);
    }

    if (geometry is LineString) {
      return operation.edit(geometry, _geomFactory);
    }

    var msg = "Unsupported Geometry class: ${geometry.runtimeType.toString()}";
    Assert.shouldNeverReachHere(msg);
    throw StateError(msg);
  }

  Polygon editPolygon(Polygon? polygon, GeometryEditorOperation operation) {
    Polygon? newPolygon = operation.edit(polygon, _geomFactory) as Polygon;
    // create one if needed
    if (newPolygon == null) newPolygon = _geomFactory!.createPolygonEmpty();
    if (newPolygon.isEmpty()) {
      //RemoveSelectedPlugIn relies on this behaviour. [Jon Aquino]
      return newPolygon;
    }

    LinearRing? shell =
        edit(newPolygon.getExteriorRing(), operation) as LinearRing;
    if (shell == null || shell.isEmpty()) {
      //RemoveSelectedPlugIn relies on this behaviour. [Jon Aquino]
      return _geomFactory!.createPolygonEmpty();
    }

    List<LinearRing> holes = [];
    for (int i = 0; i < newPolygon.getNumInteriorRing(); i++) {
      LinearRing? hole =
          edit(newPolygon.getInteriorRingN(i), operation) as LinearRing;
      if (hole == null || hole.isEmpty()) {
        continue;
      }
      holes.add(hole);
    }

    return _geomFactory!.createPolygon(shell, holes);
  }

  GeometryCollection editGeometryCollection(
      GeometryCollection collection, GeometryEditorOperation operation) {
    // first edit the entire collection
    // MD - not sure why this is done - could just check original collection?
    GeometryCollection collectionForType =
        operation.edit(collection, _geomFactory) as GeometryCollection;

    // edit the component geometries
    List geometries = [];
    for (int i = 0; i < collectionForType.getNumGeometries(); i++) {
      Geometry? geometry = edit(collectionForType.getGeometryN(i), operation);
      if (geometry == null || geometry.isEmpty()) {
        continue;
      }
      geometries.add(geometry);
    }

    if (collectionForType is MultiPoint) {
      return _geomFactory!.createMultiPoint(geometries.cast<Point>());
    }
    if (collectionForType is MultiLineString) {
      return _geomFactory!.createMultiLineString(geometries.cast<LineString>());
    }
    if (collectionForType is MultiPolygon) {
      return _geomFactory!.createMultiPolygon(geometries.cast<Polygon>());
    }
    return _geomFactory!.createGeometryCollection(geometries.cast<Geometry>());
  }
}

/**
 * A interface which specifies an edit operation for Geometries.
 *
 * @version 1.7
 */
abstract class GeometryEditorOperation {
  /**
   * Edits a Geometry by returning a new Geometry with a modification.
   * The returned geometry may be:
   * <ul>
   * <li>the input geometry itself.
   * The returned Geometry might be the same as the Geometry passed in.
   * <li><code>null</code> if the geometry is to be deleted.
   * </ul>
   *
   * @param geometry the Geometry to modify
   * @param factory the factory with which to construct the modified Geometry
   * (may be different to the factory of the input geometry)
   * @return a new Geometry which is a modification of the input Geometry
   * @return null if the Geometry is to be deleted completely
   */
  Geometry? edit(Geometry? geometry, GeometryFactory? gfactory);
}

/**
 * A GeometryEditorOperation which does not modify
 * the input geometry.
 * This can be used for simple changes of
 * GeometryFactory (including PrecisionModel and SRID).
 *
 * @author mbdavis
 *
 */
class NoOpGeometryOperation implements GeometryEditorOperation {
  Geometry? edit(Geometry? geometry, GeometryFactory? gfactory) {
    return geometry;
  }
}

/**
 * A {@link GeometryEditorOperation} which edits the coordinate list of a {@link Geometry}.
 * Operates on Geometry subclasses which contains a single coordinate list.
 */
abstract class CoordinateOperation implements GeometryEditorOperation {
  Geometry? edit(Geometry? geometry, GeometryFactory? gfactory) {
    if (geometry is LinearRing) {
      return gfactory
          ?.createLinearRing(editCoords(geometry.getCoordinates(), geometry));
    }

    if (geometry is LineString) {
      return gfactory
          ?.createLineString(editCoords(geometry.getCoordinates(), geometry));
    }

    if (geometry is Point) {
      List<Coordinate> newCoordinates =
          editCoords(geometry.getCoordinates(), geometry);

      return gfactory
          ?.createPoint((newCoordinates.length > 0) ? newCoordinates[0] : null);
    }

    return geometry;
  }

  /**
   * Edits the array of {@link Coordinate}s from a {@link Geometry}.
   * <p>
   * If it is desired to preserve the immutability of Geometrys,
   * if the coordinates are changed a new array should be created
   * and returned.
   *
   * @param coordinates the coordinate array to operate on
   * @param geometry the geometry containing the coordinate list
   * @return an edited coordinate array (which may be the same as the input)
   */
  List<Coordinate> editCoords(List<Coordinate> coordinates, Geometry geometry);
}

/**
 * A {@link GeometryEditorOperation} which edits the {@link CoordinateSequence}
 * of a {@link Geometry}.
 * Operates on Geometry subclasses which contains a single coordinate list.
 */
abstract class CoordinateSequenceOperation implements GeometryEditorOperation {
  Geometry? edit(Geometry? geometry, GeometryFactory? gfactory) {
    if (geometry is LinearRing) {
      return gfactory?.createLinearRingSeq(
          editSeq((geometry).getCoordinateSequence(), geometry));
    }

    if (geometry is LineString) {
      return gfactory?.createLineStringSeq(
          editSeq((geometry).getCoordinateSequence(), geometry));
    }

    if (geometry is Point) {
      return gfactory?.createPointSeq(
          editSeq((geometry).getCoordinateSequence(), geometry));
    }

    return geometry;
  }

  /**
   * Edits a {@link CoordinateSequence} from a {@link Geometry}.
   *
   * @param coordSeq the coordinate array to operate on
   * @param geometry the geometry containing the coordinate list
   * @return an edited coordinate sequence (which may be the same as the input)
   */
  CoordinateSequence editSeq(CoordinateSequence coordSeq, Geometry geometry);
}
