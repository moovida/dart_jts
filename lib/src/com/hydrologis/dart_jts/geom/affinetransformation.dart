part of dart_jts;

/*
 * Copyright (c) 2016 Vivid Solutions.
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse  License 2.0
 * and Eclipse Distribution License v. 1.0 which accompanies this distribution.
 * The Eclipse  License is available at http://www.eclipse.org/legal/epl-v20.html
 * and the Eclipse Distribution License is available at
 *
 * http://www.eclipse.org/org/documents/edl-v10.php.
 */

/**
 * Represents an affine transformation on the 2D Cartesian plane. 
 * It can be used to transform a {@link Coordinate} or {@link Geometry}.
 * An affine transformation is a mapping of the 2D plane into itself
 * via a series of transformations of the following basic types:
 * <ul>
 * <li>reflection (through a line)
 * <li>rotation (around the origin)
 * <li>scaling (relative to the origin)
 * <li>shearing (in both the X and Y directions)
 * <li>translation 
 * </ul>
 * In general, affine transformations preserve straightness and parallel lines,
 * but do not preserve distance or shape.
 * <p>
 * An affine transformation can be represented by a 3x3 
 * matrix in the following form:
 * <blockquote><pre>
 * T = | m00 m01 m02 |
 *     | m10 m11 m12 |
 *     |  0   0   1  |
 * </pre></blockquote>
 * A coordinate P = (x, y) can be transformed to a new coordinate P' = (x', y')
 * by representing it as a 3x1 matrix and using matrix multiplication to compute:
 * <blockquote><pre>
 * | x' |  = T x | x |
 * | y' |        | y |
 * | 1  |        | 1 |
 * </pre></blockquote>
 * <h3>Transformation Composition</h3>
 * Affine transformations can be composed using the {@link #compose} method.
 * Composition is computed via multiplication of the 
 * transformation matrices, and is defined as:
 * <blockquote><pre>
 * A.compose(B) = T<sub>B</sub> x T<sub>A</sub>
 * </pre></blockquote>
 * This produces a transformation whose effect is that of A followed by B.
 * The methods {@link #reflect}, {@link #rotate}, 
 * {@link #scale}, {@link #shear}, and {@link #translate} 
 * have the effect of composing a transformation of that type with
 * the transformation they are invoked on.  
 * <p>
 * The composition of transformations is in general <i>not</i> commutative.
 * 
 * <h3>Transformation Inversion</h3>
 * Affine transformations may be invertible or non-invertible.  
 * If a transformation is invertible, then there exists 
 * an inverse transformation which when composed produces 
 * the identity transformation.  
 * The {@link #getInverse} method
 * computes the inverse of a transformation, if one exists.
 * 
 * @author Martin Davis
 *
 */
class AffineTransformation implements CoordinateSequenceFilter {
  /**
   * Creates a transformation for a reflection about the 
   * line (x0,y0) - (x1,y1).
   * 
   * @param x0 the x-ordinate of a point on the reflection line
   * @param y0 the y-ordinate of a point on the reflection line
   * @param x1 the x-ordinate of a another point on the reflection line
   * @param y1 the y-ordinate of a another point on the reflection line
   * @return a transformation for the reflection
   */
  static AffineTransformation reflectionInstance(
      double x0, double y0, double x1, double y1) {
    AffineTransformation trans = new AffineTransformation();
    trans.setToReflection(x0, y0, x1, y1);
    return trans;
  }

  /**
   * Creates a transformation for a reflection about the 
   * line (0,0) - (x,y).
   * 
   * @param x the x-ordinate of a point on the reflection line
   * @param y the y-ordinate of a point on the reflection line
   * @return a transformation for the reflection
   */
  static AffineTransformation reflectionInstanceXY(double x, double y) {
    AffineTransformation trans = new AffineTransformation();
    trans.setToReflectionXY(x, y);
    return trans;
  }

  /**
   * Creates a transformation for a rotation
   * about the origin 
   * by an angle <i>theta</i>.
   * Positive angles correspond to a rotation 
   * in the counter-clockwise direction.
   * 
   * @param theta the rotation angle, in radians
   * @return a transformation for the rotation
   */
  static AffineTransformation rotationInstance(double theta) {
    return rotationInstanceSinCos(math.sin(theta), math.cos(theta));
  }

  /**
   * Creates a transformation for a rotation 
   * by an angle <i>theta</i>,
   * specified by the sine and cosine of the angle.
   * This allows providing exact values for sin(theta) and cos(theta)
   * for the common case of rotations of multiples of quarter-circles. 
   * 
   * @param sinTheta the sine of the rotation angle
   * @param cosTheta the cosine of the rotation angle
   * @return a transformation for the rotation
   */
  static AffineTransformation rotationInstanceSinCos(
      double sinTheta, double cosTheta) {
    AffineTransformation trans = new AffineTransformation();
    trans.setToRotationSinCos(sinTheta, cosTheta);
    return trans;
  }

  /**
   * Creates a transformation for a rotation
   * about the point (x,y) by an angle <i>theta</i>.
   * Positive angles correspond to a rotation 
   * in the counter-clockwise direction.
   * 
   * @param theta the rotation angle, in radians
   * @param x the x-ordinate of the rotation point
   * @param y the y-ordinate of the rotation point
   * @return a transformation for the rotation
   */
  static AffineTransformation rotationInstanceTXY(
      double theta, double x, double y) {
    return rotationInstanceSinCosXY(math.sin(theta), math.cos(theta), x, y);
  }

  /**
   * Creates a transformation for a rotation 
   * about the point (x,y) by an angle <i>theta</i>,
   * specified by the sine and cosine of the angle.
   * This allows providing exact values for sin(theta) and cos(theta)
   * for the common case of rotations of multiples of quarter-circles. 
   * 
   * @param sinTheta the sine of the rotation angle
   * @param cosTheta the cosine of the rotation angle
   * @param x the x-ordinate of the rotation point
   * @param y the y-ordinate of the rotation point
   * @return a transformation for the rotation
   */
  static AffineTransformation rotationInstanceSinCosXY(
      double sinTheta, double cosTheta, double x, double y) {
    AffineTransformation trans = new AffineTransformation();
    trans.setToRotationSinCosXY(sinTheta, cosTheta, x, y);
    return trans;
  }

  /**
   * Creates a transformation for a scaling relative to the origin.
   * 
   * @param xScale the value to scale by in the x direction
   * @param yScale the value to scale by in the y direction
   * @return a transformation for the scaling
   */
  static AffineTransformation scaleInstance(double xScale, double yScale) {
    AffineTransformation trans = new AffineTransformation();
    trans.setToScale(xScale, yScale);
    return trans;
  }

  /**
   * Creates a transformation for a scaling relative to the point (x,y).
   * 
   * @param xScale the value to scale by in the x direction
   * @param yScale the value to scale by in the y direction
   * @param x the x-ordinate of the point to scale around
   * @param y the y-ordinate of the point to scale around
   * @return a transformation for the scaling
   */
  static AffineTransformation scaleInstanceScaleXY(
      double xScale, double yScale, double x, double y) {
    AffineTransformation trans = new AffineTransformation();
    trans.translate(-x, -y);
    trans.scale(xScale, yScale);
    trans.translate(x, y);
    return trans;
  }

  /**
   * Creates a transformation for a shear.
   * 
   * @param xShear the value to shear by in the x direction
   * @param yShear the value to shear by in the y direction
   * @return a transformation for the shear
   */
  static AffineTransformation shearInstance(double xShear, double yShear) {
    AffineTransformation trans = new AffineTransformation();
    trans.setToShear(xShear, yShear);
    return trans;
  }

  /**
   * Creates a transformation for a translation.
   * 
   * @param x the value to translate by in the x direction
   * @param y the value to translate by in the y direction
   * @return a transformation for the translation
   */
  static AffineTransformation translationInstance(double x, double y) {
    AffineTransformation trans = new AffineTransformation();
    trans.setToTranslation(x, y);
    return trans;
  }

  // affine matrix entries
  // (bottom row is always [ 0 0 1 ])
  late double m00;
  late double m01;
  late double m02;
  late double m10;
  late double m11;
  late double m12;

  /**
   * Constructs a new identity transformation
   */
  AffineTransformation() {
    setToIdentity();
  }

  /**
   * Constructs a new transformation whose 
   * matrix has the specified values.
   * 
   * @param matrix an array containing the 6 values { m00, m01, m02, m10, m11, m12 }
   * @throws NullPointerException if matrix is null
   * @throws ArrayIndexOutOfBoundsException if matrix is too small 
   */
  AffineTransformation.fromMatrix(List<double> matrix) {
    m00 = matrix[0];
    m01 = matrix[1];
    m02 = matrix[2];
    m10 = matrix[3];
    m11 = matrix[4];
    m12 = matrix[5];
  }

  /**
   * Constructs a new transformation whose 
   * matrix has the specified values.
   * 
   * @param m00 the entry for the [0, 0] element in the transformation matrix 
   * @param m01 the entry for the [0, 1] element in the transformation matrix
   * @param m02 the entry for the [0, 2] element in the transformation matrix
   * @param m10 the entry for the [1, 0] element in the transformation matrix
   * @param m11 the entry for the [1, 1] element in the transformation matrix
   * @param m12 the entry for the [1, 2] element in the transformation matrix
   */
  AffineTransformation.fromMatrixValues(
      double m00, double m01, double m02, double m10, double m11, double m12) {
    setTransformation(m00, m01, m02, m10, m11, m12);
  }

  /**
   * Constructs a transformation which is
   * a copy of the given one.
   * 
   * @param trans the transformation to copy
   */
  AffineTransformation.fromTransformation(AffineTransformation trans) {
    setTransformationFromTransformation(trans);
  }

  /**
   * Constructs a transformation
   * which maps the given source
   * points into the given destination points.
   * 
   * @param src0 source point 0
   * @param src1 source point 1
   * @param src2 source point 2
   * @param dest0 the mapped point for source point 0
   * @param dest1 the mapped point for source point 1
   * @param dest2 the mapped point for source point 2
   * 
   */
  AffineTransformation.fromCoordinates(Coordinate src0, Coordinate src1,
      Coordinate src2, Coordinate dest0, Coordinate dest1, Coordinate dest2) {}

  /**
   * Sets this transformation to be the identity transformation.
   * The identity transformation has the matrix:
   * <blockquote><pre>
   * | 1 0 0 |
   * | 0 1 0 |
   * | 0 0 1 |
   * </pre></blockquote>
   * @return this transformation, with an updated matrix
   */
  AffineTransformation setToIdentity() {
    m00 = 1.0;
    m01 = 0.0;
    m02 = 0.0;
    m10 = 0.0;
    m11 = 1.0;
    m12 = 0.0;
    return this;
  }

  /**
   * Sets this transformation's matrix to have the given values.
   * 
   * @param m00 the entry for the [0, 0] element in the transformation matrix 
   * @param m01 the entry for the [0, 1] element in the transformation matrix
   * @param m02 the entry for the [0, 2] element in the transformation matrix
   * @param m10 the entry for the [1, 0] element in the transformation matrix
   * @param m11 the entry for the [1, 1] element in the transformation matrix
   * @param m12 the entry for the [1, 2] element in the transformation matrix
   * @return this transformation, with an updated matrix
   */
  AffineTransformation setTransformation(
      double m00, double m01, double m02, double m10, double m11, double m12) {
    this.m00 = m00;
    this.m01 = m01;
    this.m02 = m02;
    this.m10 = m10;
    this.m11 = m11;
    this.m12 = m12;
    return this;
  }

  /**
   * Sets this transformation to be a copy of the given one
   * 
   * @param trans a transformation to copy
   * @return this transformation, with an updated matrix
   */
  AffineTransformation setTransformationFromTransformation(
      AffineTransformation trans) {
    m00 = trans.m00;
    m01 = trans.m01;
    m02 = trans.m02;
    m10 = trans.m10;
    m11 = trans.m11;
    m12 = trans.m12;
    return this;
  }

  /**
   * Gets an array containing the entries
   * of the transformation matrix.
   * Only the 6 non-trivial entries are returned,
   * in the sequence:
   * <pre>
   * m00, m01, m02, m10, m11, m12
   * </pre>
   * 
   * @return an array of length 6
   */
  List<double> getMatrixEntries() {
    return [m00, m01, m02, m10, m11, m12];
  }

  /**
  * Computes the determinant of the transformation matrix. 
  * The determinant is computed as:
  * <blockquote><pre>
  * | m00 m01 m02 |
  * | m10 m11 m12 | = m00 * m11 - m01 * m10
  * |  0   0   1  |
  * </pre></blockquote>
  * If the determinant is zero, 
  * the transform is singular (not invertible), 
  * and operations which attempt to compute
  * an inverse will throw a <tt>NoninvertibleTransformException</tt>. 

  * @return the determinant of the transformation
  * @see #getInverse()
  */
  double getDeterminant() {
    return m00 * m11 - m01 * m10;
  }

  /**
   * Computes the inverse of this transformation, if one
   * exists.
   * The inverse is the transformation which when 
   * composed with this one produces the identity 
   * transformation.
   * A transformation has an inverse if and only if it
   * is not singular (i.e. its
   * determinant is non-zero).
   * Geometrically, an transformation is non-invertible
   * if it maps the plane to a line or a point.
   * If no inverse exists this method
   * will throw a <tt>NoninvertibleTransformationException</tt>.
   * <p>
   * The matrix of the inverse is equal to the 
   * inverse of the matrix for the transformation.
   * It is computed as follows:
   * <blockquote><pre>  
   *                 1    
   * inverse(A)  =  ---   x  adjoint(A) 
   *                det 
   *
   *
   *             =   1       |  m11  -m01   m01*m12-m02*m11  |
   *                ---   x  | -m10   m00  -m00*m12+m10*m02  |
   *                det      |  0     0     m00*m11-m10*m01  |
   *
   *
   *
   *             = |  m11/det  -m01/det   m01*m12-m02*m11/det |
   *               | -m10/det   m00/det  -m00*m12+m10*m02/det |
   *               |   0           0          1               |
   *
   * </pre></blockquote>  
   *  
   * @return a new inverse transformation
   * @throws NoninvertibleTransformationException
   * @see #getDeterminant()
   */
  AffineTransformation getInverse() {
    double det = getDeterminant();
    if (det == 0) throw StateError("Transformation is non-invertible");

    double im00 = m11 / det;
    double im10 = -m10 / det;
    double im01 = -m01 / det;
    double im11 = m00 / det;
    double im02 = (m01 * m12 - m02 * m11) / det;
    double im12 = (-m00 * m12 + m10 * m02) / det;

    return AffineTransformation.fromMatrixValues(
        im00, im01, im02, im10, im11, im12);
  }

  /**
   * Explicitly computes the math for a reflection.  May not work.
   * @param x0 the X ordinate of one point on the reflection line
   * @param y0 the Y ordinate of one point on the reflection line
   * @param x1 the X ordinate of another point on the reflection line
   * @param y1 the Y ordinate of another point on the reflection line
   * @return this transformation, with an updated matrix
   */
  AffineTransformation setToReflectionBasic(
      double x0, double y0, double x1, double y1) {
    if (x0 == x1 && y0 == y1) {
      throw ArgumentError("Reflection line points must be distinct");
    }
    double dx = x1 - x0;
    double dy = y1 - y0;
    double d = math.sqrt(dx * dx + dy * dy);
    double sin = dy / d;
    double cos = dx / d;
    double cs2 = 2 * sin * cos;
    double c2s2 = cos * cos - sin * sin;
    m00 = c2s2;
    m01 = cs2;
    m02 = 0.0;
    m10 = cs2;
    m11 = -c2s2;
    m12 = 0.0;
    return this;
  }

  /**
   * Sets this transformation to be a reflection 
   * about the line defined by a line <tt>(x0,y0) - (x1,y1)</tt>.
   * 
   * @param x0 the X ordinate of one point on the reflection line
   * @param y0 the Y ordinate of one point on the reflection line
   * @param x1 the X ordinate of another point on the reflection line
   * @param y1 the Y ordinate of another point on the reflection line
   * @return this transformation, with an updated matrix
   */
  AffineTransformation setToReflection(
      double x0, double y0, double x1, double y1) {
    if (x0 == x1 && y0 == y1) {
      throw ArgumentError("Reflection line points must be distinct");
    }
    // translate line vector to origin
    setToTranslation(-x0, -y0);

    // rotate vector to positive x axis direction
    double dx = x1 - x0;
    double dy = y1 - y0;
    double d = math.sqrt(dx * dx + dy * dy);
    double sin = dy / d;
    double cos = dx / d;
    rotateSinCos(-sin, cos);
    // reflect about the x axis
    scale(1, -1);
    // rotate back
    rotateSinCos(sin, cos);
    // translate back
    translate(x0, y0);
    return this;
  }

  /**
   * Sets this transformation to be a reflection 
   * about the line defined by vector (x,y).
   * The transformation for a reflection
   * is computed by:
   * <blockquote><pre>
   * d = sqrt(x<sup>2</sup> + y<sup>2</sup>)  
   * sin = y / d;
   * cos = x / d;
   * 
   * T<sub>ref</sub> = T<sub>rot(sin, cos)</sub> x T<sub>scale(1, -1)</sub> x T<sub>rot(-sin, cos)</sub>
   * </pre></blockquote> 
   * 
   * @param x the x-component of the reflection line vector
   * @param y the y-component of the reflection line vector
   * @return this transformation, with an updated matrix
   */
  AffineTransformation setToReflectionXY(double x, double y) {
    if (x == 0.0 && y == 0.0) {
      throw ArgumentError("Reflection vector must be non-zero");
    }

    /**
     * Handle special case - x = y.
     * This case is specified explicitly to avoid roundoff error.
     */
    if (x == y) {
      m00 = 0.0;
      m01 = 1.0;
      m02 = 0.0;
      m10 = 1.0;
      m11 = 0.0;
      m12 = 0.0;
      return this;
    }

    // rotate vector to positive x axis direction
    double d = math.sqrt(x * x + y * y);
    double sin = y / d;
    double cos = x / d;
    rotateSinCos(-sin, cos);
    // reflect about the x-axis
    scale(1, -1);
    // rotate back
    rotateSinCos(sin, cos);
    return this;
  }

  /**
   * Sets this transformation to be a rotation around the origin.
   * A positive rotation angle corresponds 
   * to a counter-clockwise rotation.
   * The transformation matrix for a rotation
   * by an angle <tt>theta</tt>
   * has the value:
   * <blockquote><pre>  
   * |  cos(theta)  -sin(theta)   0 |
   * |  sin(theta)   cos(theta)   0 |
   * |           0            0   1 |
   * </pre></blockquote> 
   * 
   * @param theta the rotation angle, in radians
   * @return this transformation, with an updated matrix
   */
  AffineTransformation setToRotation(double theta) {
    setToRotationSinCos(math.sin(theta), math.cos(theta));
    return this;
  }

  /**
   * Sets this transformation to be a rotation around the origin
   * by specifying the sin and cos of the rotation angle directly.
   * The transformation matrix for the rotation
   * has the value:
   * <blockquote><pre>  
   * |  cosTheta  -sinTheta   0 |
   * |  sinTheta   cosTheta   0 |
   * |         0          0   1 |
   * </pre></blockquote> 
   * 
   * @param sinTheta the sine of the rotation angle
   * @param cosTheta the cosine of the rotation angle
   * @return this transformation, with an updated matrix
   */
  AffineTransformation setToRotationSinCos(double sinTheta, double cosTheta) {
    m00 = cosTheta;
    m01 = -sinTheta;
    m02 = 0.0;
    m10 = sinTheta;
    m11 = cosTheta;
    m12 = 0.0;
    return this;
  }

  /**
   * Sets this transformation to be a rotation
   * around a given point (x,y).
   * A positive rotation angle corresponds 
   * to a counter-clockwise rotation.
   * The transformation matrix for a rotation
   * by an angle <tt>theta</tt>
   * has the value:
   * <blockquote><pre>  
   * |  cosTheta  -sinTheta   x-x*cos+y*sin |
   * |  sinTheta   cosTheta   y-x*sin-y*cos |
   * |           0            0   1 |
   * </pre></blockquote> 
   * 
   * @param theta the rotation angle, in radians
   * @param x the x-ordinate of the rotation point
   * @param y the y-ordinate of the rotation point
   * @return this transformation, with an updated matrix
   */
  AffineTransformation setToRotationTXY(double theta, double x, double y) {
    setToRotationSinCosXY(math.sin(theta), math.cos(theta), x, y);
    return this;
  }

  /**
   * Sets this transformation to be a rotation
   * around a given point (x,y)
   * by specifying the sin and cos of the rotation angle directly.
   * The transformation matrix for the rotation
   * has the value:
   * <blockquote><pre>  
   * |  cosTheta  -sinTheta   x-x*cos+y*sin |
   * |  sinTheta   cosTheta   y-x*sin-y*cos |
   * |         0          0         1       |
   * </pre></blockquote> 
   * 
   * @param sinTheta the sine of the rotation angle
   * @param cosTheta the cosine of the rotation angle
   * @param x the x-ordinate of the rotation point
   * @param y the y-ordinate of the rotation point
   * @return this transformation, with an updated matrix
   */
  AffineTransformation setToRotationSinCosXY(
      double sinTheta, double cosTheta, double x, double y) {
    m00 = cosTheta;
    m01 = -sinTheta;
    m02 = x - x * cosTheta + y * sinTheta;
    m10 = sinTheta;
    m11 = cosTheta;
    m12 = y - x * sinTheta - y * cosTheta;
    return this;
  }

  /**
   * Sets this transformation to be a scaling.
   * The transformation matrix for a scale
   * has the value:
   * <blockquote><pre>  
   * |  xScale      0  dx |
   * |  1      yScale  dy |
   * |  0           0   1 |
   * </pre></blockquote> 
   * 
   * @param xScale the amount to scale x-ordinates by
   * @param yScale the amount to scale y-ordinates by
   * @return this transformation, with an updated matrix
   */
  AffineTransformation setToScale(double xScale, double yScale) {
    m00 = xScale;
    m01 = 0.0;
    m02 = 0.0;
    m10 = 0.0;
    m11 = yScale;
    m12 = 0.0;
    return this;
  }

  /**
   * Sets this transformation to be a shear.
   * The transformation matrix for a shear 
   * has the value:
   * <blockquote><pre>  
   * |  1      xShear  0 |
   * |  yShear      1  0 |
   * |  0           0  1 |
   * </pre></blockquote> 
   * Note that a shear of (1, 1) is <i>not</i> 
   * equal to shear(1, 0) composed with shear(0, 1).
   * Instead, shear(1, 1) corresponds to a mapping onto the 
   * line x = y.
   * 
   * @param xShear the x component to shear by
   * @param yShear the y component to shear by
   * @return this transformation, with an updated matrix
   */
  AffineTransformation setToShear(double xShear, double yShear) {
    m00 = 1.0;
    m01 = xShear;
    m02 = 0.0;
    m10 = yShear;
    m11 = 1.0;
    m12 = 0.0;
    return this;
  }

  /**
   * Sets this transformation to be a translation.
   * For a translation by the vector (x, y)
   * the transformation matrix has the value:
   * <blockquote><pre>  
   * |  1  0  dx |
   * |  1  0  dy |
   * |  0  0   1 |
   * </pre></blockquote> 
   * @param dx the x component to translate by
   * @param dy the y component to translate by
   * @return this transformation, with an updated matrix
   */
  AffineTransformation setToTranslation(double dx, double dy) {
    m00 = 1.0;
    m01 = 0.0;
    m02 = dx;
    m10 = 0.0;
    m11 = 1.0;
    m12 = dy;
    return this;
  }

  /**
   * Updates the value of this transformation
   * to that of a reflection transformation composed 
   * with the current value.
   * 
   * @param x0 the x-ordinate of a point on the line to reflect around
   * @param y0 the y-ordinate of a point on the line to reflect around
   * @param x1 the x-ordinate of a point on the line to reflect around
   * @param y1 the y-ordinate of a point on the line to reflect around
   * @return this transformation, with an updated matrix
   */
  AffineTransformation reflect(double x0, double y0, double x1, double y1) {
    compose(reflectionInstance(x0, y0, x1, y1));
    return this;
  }

  /**
   * Updates the value of this transformation
   * to that of a reflection transformation composed 
   * with the current value.
   * 
   * @param x the x-ordinate of the line to reflect around
   * @param y the y-ordinate of the line to reflect around
   * @return this transformation, with an updated matrix
   */
  AffineTransformation reflectXY(double x, double y) {
    compose(reflectionInstanceXY(x, y));
    return this;
  }

  /**
   * Updates the value of this transformation
   * to that of a rotation transformation composed 
   * with the current value.
   * Positive angles correspond to a rotation 
   * in the counter-clockwise direction.
   * 
   * @param theta the angle to rotate by, in radians
   * @return this transformation, with an updated matrix
   */
  AffineTransformation rotate(double theta) {
    compose(rotationInstance(theta));
    return this;
  }

  /**
   * Updates the value of this transformation
   * to that of a rotation around the origin composed 
   * with the current value,
   * with the sin and cos of the rotation angle specified directly.
   * 
   * @param sinTheta the sine of the angle to rotate by
   * @param cosTheta the cosine of the angle to rotate by
   * @return this transformation, with an updated matrix
   */
  AffineTransformation rotateSinCos(double sinTheta, double cosTheta) {
    compose(rotationInstanceSinCos(sinTheta, cosTheta));
    return this;
  }

  /**
   * Updates the value of this transformation
   * to that of a rotation around a given point composed 
   * with the current value.
   * Positive angles correspond to a rotation 
   * in the counter-clockwise direction.
   * 
   * @param theta the angle to rotate by, in radians
   * @param x the x-ordinate of the rotation point
   * @param y the y-ordinate of the rotation point
   * @return this transformation, with an updated matrix
   */
  AffineTransformation rotateTXY(double theta, double x, double y) {
    compose(rotationInstanceTXY(theta, x, y));
    return this;
  }

  /**
   * Updates the value of this transformation
   * to that of a rotation around a given point composed 
   * with the current value,
   * with the sin and cos of the rotation angle specified directly.
   * 
   * @param sinTheta the sine of the angle to rotate by
   * @param cosTheta the cosine of the angle to rotate by
   * @param x the x-ordinate of the rotation point
   * @param y the y-ordinate of the rotation point
   * @return this transformation, with an updated matrix
   */
  AffineTransformation rotateSinCosXY(
      double sinTheta, double cosTheta, double x, double y) {
    compose(rotationInstanceSinCosXY(sinTheta, cosTheta, x, y));
    return this;
  }

  /**
   * Updates the value of this transformation
   * to that of a scale transformation composed 
   * with the current value.
   * 
   * @param xScale the value to scale by in the x direction
   * @param yScale the value to scale by in the y direction
   * @return this transformation, with an updated matrix
   */
  AffineTransformation scale(double xScale, double yScale) {
    compose(scaleInstance(xScale, yScale));
    return this;
  }

  /**
   * Updates the value of this transformation
   * to that of a shear transformation composed 
   * with the current value.
   * 
   * @param xShear the value to shear by in the x direction
   * @param yShear the value to shear by in the y direction
   * @return this transformation, with an updated matrix
   */
  AffineTransformation shear(double xShear, double yShear) {
    compose(shearInstance(xShear, yShear));
    return this;
  }

  /**
   * Updates the value of this transformation
   * to that of a translation transformation composed 
   * with the current value.
   * 
   * @param x the value to translate by in the x direction
   * @param y the value to translate by in the y direction
   * @return this transformation, with an updated matrix
   */
  AffineTransformation translate(double x, double y) {
    compose(translationInstance(x, y));
    return this;
  }

  /**
   * Updates this transformation to be
   * the composition of this transformation with the given {@link AffineTransformation}. 
   * This produces a transformation whose effect 
   * is equal to applying this transformation 
   * followed by the argument transformation.
   * mathematically,
   * <blockquote><pre>
   * A.compose(B) = T<sub>B</sub> x T<sub>A</sub>
   * </pre></blockquote>
   * 
   * @param trans an affine transformation
   * @return this transformation, with an updated matrix
   */
  AffineTransformation compose(AffineTransformation trans) {
    double mp00 = trans.m00 * m00 + trans.m01 * m10;
    double mp01 = trans.m00 * m01 + trans.m01 * m11;
    double mp02 = trans.m00 * m02 + trans.m01 * m12 + trans.m02;
    double mp10 = trans.m10 * m00 + trans.m11 * m10;
    double mp11 = trans.m10 * m01 + trans.m11 * m11;
    double mp12 = trans.m10 * m02 + trans.m11 * m12 + trans.m12;
    m00 = mp00;
    m01 = mp01;
    m02 = mp02;
    m10 = mp10;
    m11 = mp11;
    m12 = mp12;
    return this;
  }

  /**
   * Updates this transformation to be the composition 
   * of a given {@link AffineTransformation} with this transformation.
   * This produces a transformation whose effect 
   * is equal to applying the argument transformation 
   * followed by this transformation.
   * mathematically,
   * <blockquote><pre>
   * A.composeBefore(B) = T<sub>A</sub> x T<sub>B</sub>
   * </pre></blockquote>
   * 
   * @param trans an affine transformation
   * @return this transformation, with an updated matrix
   */
  AffineTransformation composeBefore(AffineTransformation trans) {
    double mp00 = m00 * trans.m00 + m01 * trans.m10;
    double mp01 = m00 * trans.m01 + m01 * trans.m11;
    double mp02 = m00 * trans.m02 + m01 * trans.m12 + m02;
    double mp10 = m10 * trans.m00 + m11 * trans.m10;
    double mp11 = m10 * trans.m01 + m11 * trans.m11;
    double mp12 = m10 * trans.m02 + m11 * trans.m12 + m12;
    m00 = mp00;
    m01 = mp01;
    m02 = mp02;
    m10 = mp10;
    m11 = mp11;
    m12 = mp12;
    return this;
  }

  /**
   * Applies this transformation to the <tt>src</tt> coordinate
   * and places the results in the <tt>dest</tt> coordinate
   * (which may be the same as the source).
   * 
   * @param src the coordinate to transform
   * @param dest the coordinate to accept the results 
   * @return the <tt>dest</tt> coordinate
   */
  Coordinate transform(Coordinate src, Coordinate dest) {
    double xp = m00 * src.x + m01 * src.y + m02;
    double yp = m10 * src.x + m11 * src.y + m12;
    dest.x = xp;
    dest.y = yp;
    return dest;
  }

  /**
   * Creates a new {@link Geometry} which is the result
   * of this transformation applied to the input Geometry.
   * 
   *@param g  a <code>Geometry</code>
   *@return a transformed Geometry
   */
  Geometry transformGeom(Geometry g) {
    Geometry g2 = g.copy();
    g2.applyCSF(this);
    return g2;
  }

  /**
   * Applies this transformation to the i'th coordinate
   * in the given CoordinateSequence.
   * 
   *@param seq  a <code>CoordinateSequence</code>
   *@param i the index of the coordinate to transform
   */
  void transformCS(CoordinateSequence seq, int i) {
    double xp = m00 * seq.getOrdinate(i, 0) + m01 * seq.getOrdinate(i, 1) + m02;
    double yp = m10 * seq.getOrdinate(i, 0) + m11 * seq.getOrdinate(i, 1) + m12;
    seq.setOrdinate(i, 0, xp);
    seq.setOrdinate(i, 1, yp);
  }

  /**
   * Transforms the i'th coordinate in the input sequence
   * 
   *@param seq  a <code>CoordinateSequence</code>
   *@param i the index of the coordinate to transform
   */
  void filter(CoordinateSequence seq, int i) {
    transformCS(seq, i);
  }

  bool isGeometryChanged() {
    return true;
  }

  /**
   * Reports that this filter should continue to be executed until 
   * all coordinates have been transformed.
   * 
   * @return false
   */
  bool isDone() {
    return false;
  }

  /**
  * Tests if this transformation is the identity transformation.
  *
  * @return true if this is the identity transformation
  */
  bool isIdentity() {
    return (m00 == 1 &&
        m01 == 0 &&
        m02 == 0 &&
        m10 == 0 &&
        m11 == 1 &&
        m12 == 0);
  }

  /**
  * Tests if an object is an
  * <tt>AffineTransformation</tt>
  * and has the same matrix as 
  * this transformation.
  * 
  * @param obj an object to test
  * @return true if the given object is equal to this object
  */
  bool equals(Object? obj) {
    if (obj == null) return false;
    if (!(obj is AffineTransformation)) return false;

    AffineTransformation trans = obj;
    return m00 == trans.m00 &&
        m01 == trans.m01 &&
        m02 == trans.m02 &&
        m10 == trans.m10 &&
        m11 == trans.m11 &&
        m12 == trans.m12;
  }

  /**
   * Gets a text representation of this transformation.
   * The string is of the form:
   * <pre>
   * AffineTransformation[[m00, m01, m02], [m10, m11, m12]]
   * </pre>
   * 
   * @return a string representing this transformation
   * 
   */
  String toString() {
    return "AffineTransformation[[$m00, $m01, $m02], [$m10, $m11, $m12]]";
  }
}

/**
 * Builds an {@link AffineTransformation} defined by a set of control vectors. 
 * A control vector consists of a source point and a destination point, 
 * which is the image of the source point under the desired transformation.
 * <p>
 * A transformation is well-defined 
 * by a set of three control vectors 
 * if and only if the source points are not collinear. 
 * (In particular, the degenerate situation
 * where two or more source points are identical will not produce a well-defined transformation).
 * A well-defined transformation exists and is unique.
 * If the control vectors are not well-defined, the system of equations
 * defining the transformation matrix entries is not solvable,
 * and no transformation can be determined.
 * <p>
 * No such restriction applies to the destination points.
 * However, if the destination points are collinear or non-unique,
 * a non-invertible transformations will be generated.
 * <p>
 * This technique of recovering a transformation
 * from its effect on known points is used in the Bilinear Interpolated Triangulation
 * algorithm for warping planar surfaces.
 *
 * @author Martin Davis
 */
class AffineTransformationBuilder {
  Coordinate src0;
  Coordinate src1;
  Coordinate src2;
  Coordinate dest0;
  Coordinate dest1;
  Coordinate dest2;

  // the matrix entries for the transformation
  late double m00, m01, m02, m10, m11, m12;

  /**
   * Constructs a new builder for
   * the transformation defined by the given 
   * set of control point mappings.
   * 
   * @param src0 a control point
   * @param src1 a control point
   * @param src2 a control point
   * @param dest0 the image of control point 0 under the required transformation
   * @param dest1 the image of control point 1 under the required transformation
   * @param dest2 the image of control point 2 under the required transformation
   */
  AffineTransformationBuilder(
      this.src0, this.src1, this.src2, this.dest0, this.dest1, this.dest2);

  /**
   * Computes the {@link AffineTransformation}
   * determined by the control point mappings,
   * or <code>null</code> if the control vectors do not determine a well-defined transformation.
   * 
   * @return an affine transformation, or null if the control vectors do not determine a well-defined transformation
   */
  AffineTransformation? getTransformation() {
    // compute full 3-point transformation
    bool isSolvable = compute();
    if (isSolvable)
      return AffineTransformation.fromMatrixValues(
          m00, m01, m02, m10, m11, m12);
    return null;
  }

  /**
   * Computes the transformation matrix by 
   * solving the two systems of linear equations
   * defined by the control point mappings,
   * if this is possible.
   * 
   * @return true if the transformation matrix is solvable
   */
  bool compute() {
    List<double> bx = [dest0.x, dest1.x, dest2.x];
    List<double?>? row0 = solve(bx);
    if (row0 == null) return false;
    m00 = row0[0]!;
    m01 = row0[1]!;
    m02 = row0[2]!;

    List<double> by = [dest0.y, dest1.y, dest2.y];
    List<double?>? row1 = solve(by);
    if (row1 == null) return false;
    m10 = row1[0]!;
    m11 = row1[1]!;
    m12 = row1[2]!;
    return true;
  }

  /**
   * Solves the transformation matrix system of linear equations
   * for the given right-hand side vector.
   * 
   * @param b the vector for the right-hand side of the system
   * @return the solution vector, or <code>null</code> if no solution could be determined
   */
  List<double?>? solve(List<double> b) {
    List<List<double>> a = [
      [src0.x, src0.y, 1],
      [src1.x, src1.y, 1],
      [src2.x, src2.y, 1]
    ];
    return Matrix.solve(a, b);
  }
}

/**
 * Supports creating {@link AffineTransformation}s defined by various kinds of
 * inputs and transformation mapping rules.
 * 
 * @author Martin Davis
 * 
 */
class AffineTransformationFactory {
  /**
	 * Creates a transformation from a set of three control vectors. A control
	 * vector consists of a source point and a destination point, which is the
	 * image of the source point under the desired transformation. Three control
	 * vectors allows defining a fully general affine transformation.
	 * 
	 * @param src0
	 * @param src1
	 * @param src2
	 * @param dest0
	 * @param dest1
	 * @param dest2
	 * @return the computed transformation
	 */
  static AffineTransformation? createFromControlVectors3(
      Coordinate src0,
      Coordinate src1,
      Coordinate src2,
      Coordinate dest0,
      Coordinate dest1,
      Coordinate dest2) {
    AffineTransformationBuilder builder =
        new AffineTransformationBuilder(src0, src1, src2, dest0, dest1, dest2);
    return builder.getTransformation();
  }

  /**
	 * Creates an AffineTransformation defined by a pair of control vectors. A
	 * control vector consists of a source point and a destination point, which is
	 * the image of the source point under the desired transformation. The
	 * computed transformation is a combination of one or more of a uniform scale,
	 * a rotation, and a translation (i.e. there is no shear component and no
	 * reflection)
	 * 
	 * @param src0
	 * @param src1
	 * @param dest0
	 * @param dest1
	 * @return the computed transformation, or null if the control vectors do not determine a well-defined transformation
	 */
  static AffineTransformation? createFromControlVectors2(
      Coordinate src0, Coordinate src1, Coordinate dest0, Coordinate dest1) {
    Coordinate rotPt = new Coordinate(dest1.x - dest0.x, dest1.y - dest0.y);

    double ang = Angle.angleBetweenOriented(src1, src0, rotPt);

    double srcDist = src1.distance(src0);
    double destDist = dest1.distance(dest0);

    if (srcDist == 0.0) return null;

    double scale = destDist / srcDist;

    AffineTransformation trans =
        AffineTransformation.translationInstance(-src0.x, -src0.y);
    trans.rotate(ang);
    trans.scale(scale, scale);
    trans.translate(dest0.x, dest0.y);
    return trans;
  }

  /**
	 * Creates an AffineTransformation defined by a single control vector. A
	 * control vector consists of a source point and a destination point, which is
	 * the image of the source point under the desired transformation. This
	 * produces a translation.
	 * 
	 * @param src0
	 *          the start point of the control vector
	 * @param dest0
	 *          the end point of the control vector
	 * @return the computed transformation
	 */
  static AffineTransformation createFromControlVectors(
      Coordinate? src0, Coordinate? dest0) {
    double dx = dest0!.x - src0!.x;
    double dy = dest0.y - src0.y;
    return AffineTransformation.translationInstance(dx, dy);
  }

  /**
	 * Creates an AffineTransformation defined by a set of control vectors.
	 * Between one and three vectors must be supplied.
	 * 
	 * @param src
	 *          the source points of the vectors
	 * @param dest
	 *          the destination points of the vectors
	 * @return the computed transformation
	 * @throws IllegalArgumentException
	 *           if the control vector arrays are too short, long or of different
	 *           lengths
	 */
  static AffineTransformation? createFromControlVectorsCL(
      List<Coordinate?> src, List<Coordinate?> dest) {
    if (src.length != dest.length)
      throw ArgumentError("Src and Dest arrays are not the same length");
    if (src.length <= 0) throw ArgumentError("Too few control points");
    if (src.length > 3) throw ArgumentError("Too many control points");

    if (src.length == 1) return createFromControlVectors(src[0], dest[0]);
    if (src.length == 2)
      return createFromControlVectors2(src[0]!, src[1]!, dest[0]!, dest[1]!);

    return createFromControlVectors3(
        src[0]!, src[1]!, src[2]!, dest[0]!, dest[1]!, dest[2]!);
  }

  /**
	 * Creates an AffineTransformation defined by a mapping between two baselines.
	 * The computed transformation consists of:
	 * <ul>
	 * <li>a translation 
	 * from the start point of the source baseline to the start point of the destination baseline,
	 * <li>a rotation through the angle between the baselines about the destination start point,
	 * <li>and a scaling equal to the ratio of the baseline lengths.
	 * </ul>
	 * If the source baseline has zero length, an identity transformation is returned.
	 * 
	 * @param src0 the start point of the source baseline
	 * @param src1 the end point of the source baseline
	 * @param dest0 the start point of the destination baseline
	 * @param dest1 the end point of the destination baseline
	 * @return the computed transformation
	 */
  static AffineTransformation createFromBaseLines(
      Coordinate src0, Coordinate src1, Coordinate dest0, Coordinate dest1) {
    Coordinate rotPt =
        new Coordinate(src0.x + dest1.x - dest0.x, src0.y + dest1.y - dest0.y);

    double ang = Angle.angleBetweenOriented(src1, src0, rotPt);

    double srcDist = src1.distance(src0);
    double destDist = dest1.distance(dest0);

    // return identity if transformation would be degenerate
    if (srcDist == 0.0) return new AffineTransformation();

    double scale = destDist / srcDist;

    AffineTransformation trans =
        AffineTransformation.translationInstance(-src0.x, -src0.y);
    trans.rotate(ang);
    trans.scale(scale, scale);
    trans.translate(dest0.x, dest0.y);
    return trans;
  }
}

class AffineTransformationFunctions {
  ///Transforms a geometry using one to three control vectors
  static Geometry transformByVectors(Geometry g, Geometry control) {
    int nControl = control.getNumGeometries();
    List<Coordinate?> src = []..length = nControl;
    List<Coordinate?> dest = []..length = nControl;
    for (int i = 0; i < nControl; i++) {
      Geometry contComp = control.getGeometryN(i);
      List<Coordinate> pts = contComp.getCoordinates();
      src[i] = pts[0];
      dest[i] = pts[1];
    }
    AffineTransformation? trans =
        AffineTransformationFactory.createFromControlVectorsCL(src, dest);
    return trans!.transformGeom(g);
  }

  /// Transforms a geometry by mapping envelope baseline to target vector
  static Geometry transformByBaseline(Geometry g, Geometry destBaseline) {
    Envelope env = g.getEnvelopeInternal();
    Coordinate src0 = new Coordinate(env.getMinX(), env.getMinY());
    Coordinate src1 = new Coordinate(env.getMaxX(), env.getMinY());

    var destPts = destBaseline.getCoordinates();
    Coordinate dest0 = destPts[0];
    Coordinate dest1 = destPts[1];
    AffineTransformation trans =
        AffineTransformationFactory.createFromBaseLines(
            src0, src1, dest0, dest1);
    return trans.transformGeom(g);
  }

  static Coordinate? envelopeCentre(Geometry g) {
    return g.getEnvelopeInternal().centre();
  }

  static Coordinate envelopeLowerLeft(Geometry g) {
    Envelope env = g.getEnvelopeInternal();
    return new Coordinate(env.getMinX(), env.getMinY());
  }

  static Geometry transformToViewport(Geometry g, Geometry gViewport) {
    Envelope viewEnv = gViewport.getEnvelopeInternal();
    Envelope env = g.getEnvelopeInternal();
    AffineTransformation trans = viewportTrans(env, viewEnv);
    return trans.transformGeom(g);
  }

  static AffineTransformation viewportTrans(Envelope srcEnv, Envelope viewEnv) {
    // works even if W or H are zero, thanks to Java infinity value.
    double scaleW = viewEnv.getWidth() / srcEnv.getWidth();
    double scaleH = viewEnv.getHeight() / srcEnv.getHeight();
    // choose minimum scale to ensure source fits viewport
    double scale = math.min(scaleW, scaleH);

    Coordinate centre = srcEnv.centre()!;
    Coordinate viewCentre = viewEnv.centre()!;

    // isotropic scaling
    AffineTransformation trans = AffineTransformation.scaleInstanceScaleXY(
        scale, scale, centre.x, centre.y);
    // translate using envelope centres
    trans.translate(viewCentre.x - centre.x, viewCentre.y - centre.y);
    return trans;
  }

  static Geometry scale(Geometry g, double scale) {
    Coordinate centre = envelopeCentre(g)!;
    AffineTransformation trans = AffineTransformation.scaleInstanceScaleXY(
        scale, scale, centre.x, centre.y);
    return trans.transformGeom(g);
  }

  static Geometry reflectInX(Geometry g) {
    Coordinate centre = envelopeCentre(g)!;
    AffineTransformation trans =
        AffineTransformation.scaleInstanceScaleXY(1, -1, centre.x, centre.y);
    return trans.transformGeom(g);
  }

  static Geometry reflectInY(Geometry g) {
    Coordinate centre = envelopeCentre(g)!;
    AffineTransformation trans =
        AffineTransformation.scaleInstanceScaleXY(-1, 1, centre.x, centre.y);
    return trans.transformGeom(g);
  }

  /// Rotate a geometry by an multiple of Pi radians
  static Geometry rotateByPiMultiple(
      Geometry g,
      // Angle (multiple of Pi)
      double multipleOfPi) {
    Coordinate centre = envelopeCentre(g)!;
    AffineTransformation trans = AffineTransformation.rotationInstanceTXY(
        multipleOfPi * math.pi, centre.x, centre.y);
    return trans.transformGeom(g);
  }

  /// Rotate a geometry around a point by an multiple of Pi radians
  static Geometry rotateByPiMultipleAroundPoint(
      Geometry g,
      Geometry pt,

      /// Angle (multiple of Pi)
      double multipleOfPi) {
    Coordinate loc;
    if (pt == null) {
      loc = new Coordinate(0, 0);
    } else {
      loc = pt.getCoordinates()[0];
    }
    AffineTransformation trans = AffineTransformation.rotationInstanceTXY(
        multipleOfPi * math.pi, loc.x, loc.y);
    return trans.transformGeom(g);
  }

  /// Rotate a geometry by an angle in radians
  static Geometry rotate(
      Geometry g,
      // Angle (radians)
      double angle) {
    Coordinate centre = envelopeCentre(g)!;
    AffineTransformation trans =
        AffineTransformation.rotationInstanceTXY(angle, centre.x, centre.y);
    return trans.transformGeom(g);
  }

  /// Rotate a geometry around a point by an angle in radians
  static Geometry rotateAroundPoint(
      Geometry g,
      Geometry pt,

      /// Angle (radians)
      double angle) {
    Coordinate loc;
    if (pt == null) {
      loc = new Coordinate(0, 0);
    } else {
      loc = pt.getCoordinates()[0];
    }
    AffineTransformation trans =
        AffineTransformation.rotationInstanceTXY(angle, loc.x, loc.y);
    return trans.transformGeom(g);
  }

  static Geometry translateCentreToOrigin(Geometry g) {
    Coordinate centre = envelopeCentre(g)!;
    AffineTransformation trans =
        AffineTransformation.translationInstance(-centre.x, -centre.y);
    return trans.transformGeom(g);
  }

  static Geometry translateToOrigin(Geometry g) {
    Coordinate lowerLeft = envelopeLowerLeft(g);
    AffineTransformation trans =
        AffineTransformation.translationInstance(-lowerLeft.x, -lowerLeft.y);
    return trans.transformGeom(g);
  }
}
