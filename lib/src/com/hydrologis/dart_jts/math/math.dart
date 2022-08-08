import 'dart:math' as math;
import '../util.dart';

/**
 * Implements extended-precision floating-point numbers
 * which maintain 106 bits (approximately 30 decimal digits) of precision.
 * <p>
 * A DoubleDouble uses a representation containing two double-precision values.
 * A number x is represented as a pair of doubles, x.hi and x.lo,
 * such that the number represented by x is x.hi + x.lo, where
 * <pre>
 *    |x.lo| &lt;= 0.5*ulp(x.hi)
 * </pre>
 * and ulp(y) means "unit in the last place of y".
 * The basic arithmetic operations are implemented using
 * convenient properties of IEEE-754 floating-point arithmetic.
 * <p>
 * The range of values which can be represented is the same as in IEEE-754.
 * The precision of the representable numbers
 * is twice as great as IEEE-754 double precision.
 * <p>
 * The correctness of the arithmetic algorithms relies on operations
 * being performed with standard IEEE-754 double precision and rounding.
 * This is the Java standard arithmetic model, but for performance reasons
 * Java implementations are not
 * constrained to using this standard by default.
 * Some processors (notably the Intel Pentium architecture) perform
 * floating point operations in (non-IEEE-754-standard) extended-precision.
 * A JVM implementation may choose to use the non-standard extended-precision
 * as its default arithmetic mode.
 * To prevent this from happening, this code uses the
 * Java <tt>strictfp</tt> modifier,
 * which forces all operations to take place in the standard IEEE-754 rounding model.
 * <p>
 * The API provides both a set of value-oriented operations
 * and a set of mutating operations.
 * Value-oriented operations treat DoubleDouble values as
 * immutable; operations on them return new objects carrying the result
 * of the operation.  This provides a simple and safe semantics for
 * writing DoubleDouble expressions.  However, there is a performance
 * penalty for the object allocations required.
 * The mutable interface updates object values in-place.
 * It provides optimum memory performance, but requires
 * care to ensure that aliasing errors are not created
 * and constant values are not changed.
 * <p>
 * For example, the following code example constructs three DD instances:
 * two to hold the input values and one to hold the result of the addition.
 * <pre>
 *     DD a = new DD(2.0);
 *     DD b = new DD(3.0);
 *     DD c = a.add(b);
 * </pre>
 * In contrast, the following approach uses only one object:
 * <pre>
 *     DD a = new DD(2.0);
 *     a.selfAdd(3.0);
 * </pre>
 * <p>
 * This implementation uses algorithms originally designed variously by
 * Knuth, Kahan, Dekker, and Linnainmaa.
 * Douglas Priest developed the first C implementation of these techniques.
 * Other more recent C++ implementation are due to Keith M. Briggs and David Bailey et al.
 *
 * <h3>References</h3>
 * <ul>
 * <li>Priest, D., <i>Algorithms for Arbitrary Precision Floating Point Arithmetic</i>,
 * in P. Kornerup and D. Matula, Eds., Proc. 10th Symposium on Computer Arithmetic,
 * IEEE Computer Society Press, Los Alamitos, Calif., 1991.
 * <li>Yozo Hida, Xiaoye S. Li and David H. Bailey,
 * <i>Quad-Double Arithmetic: Algorithms, Implementation, and Application</i>,
 * manuscript, Oct 2000; Lawrence Berkeley National Laboratory Report BNL-46996.
 * <li>David Bailey, <i>High Precision Software Directory</i>;
 * <tt>http://crd.lbl.gov/~dhbailey/mpdist/index.html</tt>
 * </ul>
 *
 *
 * @author Martin Davis
 *
 */
class DD implements Comparable {
  /**
   * The value nearest to the constant Pi.
   */
  static final DD PI = new DD.withHiLo(3.141592653589793116e+00, 1.224646799147353207e-16);

  /**
   * The value nearest to the constant 2 * Pi.
   */
  static final DD TWO_PI = new DD.withHiLo(6.283185307179586232e+00, 2.449293598294706414e-16);

  /**
   * The value nearest to the constant Pi / 2.
   */
  static final DD PI_2 = new DD.withHiLo(1.570796326794896558e+00, 6.123233995736766036e-17);

  /**
   * The value nearest to the constant e (the natural logarithm base).
   */
  static final DD E = new DD.withHiLo(2.718281828459045091e+00, 1.445646891729250158e-16);

  /**
   * A value representing the result of an operation which does not return a valid number.
   */
  static final DD NaN = new DD.withHiLo(double.nan, double.nan);

  /**
   * The smallest representable relative difference between two {link @ DoubleDouble} values
   */
  static final double EPS = 1.23259516440783e-32;

  /* = 2^-106 */

  static DD createNaN() {
    return new DD.withHiLo(double.nan, double.nan);
  }

  /**
   * Converts the string argument to a DoubleDouble number.
   *
   * @param str a string containing a representation of a numeric value
   * @return the extended precision version of the value
   * @throws NumberFormatException if <tt>s</tt> is not a valid representation of a number
   */
  static DD valueOfStr(String str) {
    return parse(str);
  }

  /**
   * Converts the <tt>double</tt> argument to a DoubleDouble number.
   *
   * @param x a numeric value
   * @return the extended precision version of the value
   */
  static DD valueOf(double x) {
    return new DD(x);
  }

  /**
   * The value to split a double-precision value on during multiplication
   */
  static final double SPLIT = 134217729.0; // 2^27+1, for IEEE double

  /**
   * The high-order component of the double-double precision value.
   */
  double hi = 0.0;

  /**
   * The low-order component of the double-double precision value.
   */
  double lo = 0.0;

  /**
   * Creates a new DoubleDouble with value 0.0.
   */
  DD.empty() {
    init(0.0);
  }

  /**
   * Creates a new DoubleDouble with value x.
   *
   * @param x the value to initialize
   */
  DD(double x) {
    init(x);
  }

  /**
   * Creates a new DoubleDouble with value (hi, lo).
   *
   * @param hi the high-order component
   * @param lo the high-order component
   */
  DD.withHiLo(double hi, double lo) {
    iniWithHiLo(hi, lo);
  }

  /**
   * Creates a new DoubleDouble with value equal to the argument.
   *
   * @param dd the value to initialize
   */
  DD.fromDD(DD dd) {
    initFromDD(dd);
  }

  /**
   * Creates a new DoubleDouble with value equal to the argument.
   *
   * @param str the value to initialize by
   * @throws NumberFormatException if <tt>str</tt> is not a valid representation of a number
   */
  DD.withString(String str) : this.fromDD(parse(str));

  /**
   * Creates a new DoubleDouble with the value of the argument.
   *
   * @param dd the DoubleDouble value to copy
   * @return a copy of the input value
   */
  static DD copy(DD dd) {
    return DD.fromDD(dd);
  }

  ///**
// * Creates and returns a copy of this value.
// *
// * @return a copy of this value
// */
// Object clone()
//{
//try {
//return super.clone();
//}
//catch (CloneNotSupportedException ex) {
//// should never reach here
//return null;
//}
//}

  void init(double x) {
    this.hi = x;
    this.lo = 0.0;
  }

  void iniWithHiLo(double hi, double lo) {
    this.hi = hi;
    this.lo = lo;
  }

  void initFromDD(DD dd) {
    hi = dd.hi;
    lo = dd.lo;
  }

/*
  double getHighComponent() { return hi; }
  
  double getLowComponent() { return lo; }
  */

// Testing only - should not be
/*
   void RENORM()
  {
    double s = hi + lo;
    double err = lo - (s - hi);
    hi = s;
    lo = err;
  }
  */

  /**
   * Set the value for the DD object. This method supports the mutating
   * operations concept described in the class documentation (see above).
   * @param value a DD instance supplying an extended-precision value.
   * @return a self-reference to the DD instance.
   */
  DD setValueDD(DD value) {
    initFromDD(value);
    return this;
  }

  /**
   * Set the value for the DD object. This method supports the mutating
   * operations concept described in the class documentation (see above).
   * @param value a floating point value to be stored in the instance.
   * @return a self-reference to the DD instance.
   */
  DD setValue(double value) {
    init(value);
    return this;
  }

  /**
   * Returns a new DoubleDouble whose value is <tt>(this + y)</tt>.
   *
   * @param y the addend
   * @return <tt>(this + y)</tt>
   */
  DD addDD(DD y) {
    return copy(this).selfAddDD(y);
  }

  /**
   * Returns a new DoubleDouble whose value is <tt>(this + y)</tt>.
   *
   * @param y the addend
   * @return <tt>(this + y)</tt>
   */
  DD add(double y) {
    return copy(this).selfAdd(y);
  }

  /**
   * Adds the argument to the value of <tt>this</tt>.
   * To prevent altering constants,
   * this method <b>must only</b> be used on values known to
   * be newly created.
   *
   * @param y the addend
   * @return this object, increased by y
   */
  DD selfAddDD(DD y) {
    return selfAddHiLo(y.hi, y.lo);
  }

  /**
   * Adds the argument to the value of <tt>this</tt>.
   * To prevent altering constants,
   * this method <b>must only</b> be used on values known to
   * be newly created.
   *
   * @param y the addend
   * @return this object, increased by y
   */
  DD selfAdd(double y) {
    double H, h, S, s, e, f;
    S = hi + y;
    e = S - hi;
    s = S - e;
    s = (y - e) + (hi - s);
    f = s + lo;
    H = S + f;
    h = f + (S - H);
    hi = H + h;
    lo = h + (H - hi);
    return this;
// return selfAdd(y, 0.0);
  }

  DD selfAddHiLo(double yhi, double ylo) {
    double H, h, T, t, S, s, e, f;
    S = hi + yhi;
    T = lo + ylo;
    e = S - hi;
    f = T - lo;
    s = S - e;
    t = T - f;
    s = (yhi - e) + (hi - s);
    t = (ylo - f) + (lo - t);
    e = s + T;
    H = S + e;
    h = e + (S - H);
    e = t + h;

    double zhi = H + e;
    double zlo = e + (H - zhi);
    hi = zhi;
    lo = zlo;
    return this;
  }

  /**
   * Computes a new DoubleDouble object whose value is <tt>(this - y)</tt>.
   *
   * @param y the subtrahend
   * @return <tt>(this - y)</tt>
   */
  DD subtractDD(DD y) {
    return addDD(y.negate());
  }

  /**
   * Computes a new DoubleDouble object whose value is <tt>(this - y)</tt>.
   *
   * @param y the subtrahend
   * @return <tt>(this - y)</tt>
   */
  DD subtract(double y) {
    return add(-y);
  }

  /**
   * Subtracts the argument from the value of <tt>this</tt>.
   * To prevent altering constants,
   * this method <b>must only</b> be used on values known to
   * be newly created.
   *
   * @param y the addend
   * @return this object, decreased by y
   */
  DD selfSubtractDD(DD y) {
    if (isNaN()) return this;
    return selfAddHiLo(-y.hi, -y.lo);
  }

  /**
   * Subtracts the argument from the value of <tt>this</tt>.
   * To prevent altering constants,
   * this method <b>must only</b> be used on values known to
   * be newly created.
   *
   * @param y the addend
   * @return this object, decreased by y
   */
  DD selfSubtract(double y) {
    if (isNaN()) return this;
    return selfAddHiLo(-y, 0.0);
  }

  /**
   * Returns a new DoubleDouble whose value is <tt>-this</tt>.
   *
   * @return <tt>-this</tt>
   */
  DD negate() {
    if (isNaN()) return this;
    return DD.withHiLo(-hi, -lo);
  }

  /**
   * Returns a new DoubleDouble whose value is <tt>(this * y)</tt>.
   *
   * @param y the multiplicand
   * @return <tt>(this * y)</tt>
   */
  DD multiplyDD(DD y) {
    if (y.isNaN()) return createNaN();
    return copy(this).selfMultiplyDD(y);
  }

  /**
   * Returns a new DoubleDouble whose value is <tt>(this * y)</tt>.
   *
   * @param y the multiplicand
   * @return <tt>(this * y)</tt>
   */
  DD multiply(double y) {
    if (y.isNaN) return createNaN();
    return copy(this).selfMultiplyHiLo(y, 0.0);
  }

  /**
   * Multiplies this object by the argument, returning <tt>this</tt>.
   * To prevent altering constants,
   * this method <b>must only</b> be used on values known to
   * be newly created.
   *
   * @param y the value to multiply by
   * @return this object, multiplied by y
   */
  DD selfMultiplyDD(DD y) {
    return selfMultiplyHiLo(y.hi, y.lo);
  }

  /**
   * Multiplies this object by the argument, returning <tt>this</tt>.
   * To prevent altering constants,
   * this method <b>must only</b> be used on values known to
   * be newly created.
   *
   * @param y the value to multiply by
   * @return this object, multiplied by y
   */
  DD selfMultiply(double y) {
    return selfMultiplyHiLo(y, 0.0);
  }

  DD selfMultiplyHiLo(double yhi, double ylo) {
    double hx, tx, hy, ty, C, c;
    C = SPLIT * hi;
    hx = C - hi;
    c = SPLIT * yhi;
    hx = C - hx;
    tx = hi - hx;
    hy = c - yhi;
    C = hi * yhi;
    hy = c - hy;
    ty = yhi - hy;
    c = ((((hx * hy - C) + hx * ty) + tx * hy) + tx * ty) + (hi * ylo + lo * yhi);
    double zhi = C + c;
    hx = C - zhi;
    double zlo = c + hx;
    hi = zhi;
    lo = zlo;
    return this;
  }

  /**
   * Computes a new DoubleDouble whose value is <tt>(this / y)</tt>.
   *
   * @param y the divisor
   * @return a new object with the value <tt>(this / y)</tt>
   */
  DD divideDD(DD y) {
    double hc, tc, hy, ty, C, c, U, u;
    C = hi / y.hi;
    c = SPLIT * C;
    hc = c - C;
    u = SPLIT * y.hi;
    hc = c - hc;
    tc = C - hc;
    hy = u - y.hi;
    U = C * y.hi;
    hy = u - hy;
    ty = y.hi - hy;
    u = (((hc * hy - U) + hc * ty) + tc * hy) + tc * ty;
    c = ((((hi - U) - u) + lo) - C * y.lo) / y.hi;
    u = C + c;

    double zhi = u;
    double zlo = (C - u) + c;
    return new DD.withHiLo(zhi, zlo);
  }

  /**
   * Computes a new DoubleDouble whose value is <tt>(this / y)</tt>.
   *
   * @param y the divisor
   * @return a new object with the value <tt>(this / y)</tt>
   */
  DD divide(double y) {
    if (y.isNaN) return createNaN();
    return copy(this).selfDivideHiLo(y, 0.0);
  }

  /**
   * Divides this object by the argument, returning <tt>this</tt>.
   * To prevent altering constants,
   * this method <b>must only</b> be used on values known to
   * be newly created.
   *
   * @param y the value to divide by
   * @return this object, divided by y
   */
  DD selfDivideDD(DD y) {
    return selfDivideHiLo(y.hi, y.lo);
  }

  /**
   * Divides this object by the argument, returning <tt>this</tt>.
   * To prevent altering constants,
   * this method <b>must only</b> be used on values known to
   * be newly created.
   *
   * @param y the value to divide by
   * @return this object, divided by y
   */
  DD selfDivide(double y) {
    return selfDivideHiLo(y, 0.0);
  }

  DD selfDivideHiLo(double yhi, double ylo) {
    double hc, tc, hy, ty, C, c, U, u;
    C = hi / yhi;
    c = SPLIT * C;
    hc = c - C;
    u = SPLIT * yhi;
    hc = c - hc;
    tc = C - hc;
    hy = u - yhi;
    U = C * yhi;
    hy = u - hy;
    ty = yhi - hy;
    u = (((hc * hy - U) + hc * ty) + tc * hy) + tc * ty;
    c = ((((hi - U) - u) + lo) - C * ylo) / yhi;
    u = C + c;

    hi = u;
    lo = (C - u) + c;
    return this;
  }

  /**
   * Returns a DoubleDouble whose value is  <tt>1 / this</tt>.
   *
   * @return the reciprocal of this value
   */
  DD reciprocal() {
    double hc, tc, hy, ty, C, c, U, u;
    C = 1.0 / hi;
    c = SPLIT * C;
    hc = c - C;
    u = SPLIT * hi;
    hc = c - hc;
    tc = C - hc;
    hy = u - hi;
    U = C * hi;
    hy = u - hy;
    ty = hi - hy;
    u = (((hc * hy - U) + hc * ty) + tc * hy) + tc * ty;
    c = ((((1.0 - U) - u)) - C * lo) / hi;

    double zhi = C + c;
    double zlo = (C - zhi) + c;
    return new DD.withHiLo(zhi, zlo);
  }

  /**
   * Returns the largest (closest to positive infinity)
   * value that is not greater than the argument
   * and is equal to a mathematical integer.
   * Special cases:
   * <ul>
   * <li>If this value is NaN, returns NaN.
   * </ul>
   *
   * @return the largest (closest to positive infinity)
   * value that is not greater than the argument
   * and is equal to a mathematical integer.
   */
  DD floor() {
    if (isNaN()) return NaN;
    double fhi = hi.floorToDouble();
    double flo = 0.0;
// Hi is already integral.  Floor the low word
    if (fhi == hi) {
      flo = lo.floorToDouble();
    }
// do we need to renormalize here?
    return new DD.withHiLo(fhi, flo);
  }

  /**
   * Returns the smallest (closest to negative infinity) value
   * that is not less than the argument and is equal to a mathematical integer.
   * Special cases:
   * <ul>
   * <li>If this value is NaN, returns NaN.
   * </ul>
   *
   * @return the smallest (closest to negative infinity) value
   * that is not less than the argument and is equal to a mathematical integer.
   */
  DD ceil() {
    if (isNaN()) return NaN;
    double fhi = hi.ceilToDouble();
    double flo = 0.0;
// Hi is already integral.  Ceil the low word
    if (fhi == hi) {
      flo = lo.ceilToDouble();
// do we need to renormalize here?
    }
    return new DD.withHiLo(fhi, flo);
  }

  /**
   * Returns an integer indicating the sign of this value.
   * <ul>
   * <li>if this value is &gt; 0, returns 1
   * <li>if this value is &lt; 0, returns -1
   * <li>if this value is = 0, returns 0
   * <li>if this value is NaN, returns 0
   * </ul>
   *
   * @return an integer indicating the sign of this value
   */
  int signum() {
    if (hi > 0) return 1;
    if (hi < 0) return -1;
    if (lo > 0) return 1;
    if (lo < 0) return -1;
    return 0;
  }

  /**
   * Rounds this value to the nearest integer.
   * The value is rounded to an integer by adding 1/2 and taking the floor of the result.
   * Special cases:
   * <ul>
   * <li>If this value is NaN, returns NaN.
   * </ul>
   *
   * @return this value rounded to the nearest integer
   */
  DD rint() {
    if (isNaN()) return this;
// may not be 100% correct
    DD plus5 = this.add(0.5);
    return plus5.floor();
  }

  /**
   * Returns the integer which is largest in absolute value and not further
   * from zero than this value.
   * Special cases:
   * <ul>
   * <li>If this value is NaN, returns NaN.
   * </ul>
   *
   * @return the integer which is largest in absolute value and not further from zero than this value
   */
  DD trunc() {
    if (isNaN()) return NaN;
    if (isPositive())
      return floor();
    else
      return ceil();
  }

  /**
   * Returns the absolute value of this value.
   * Special cases:
   * <ul>
   * <li>If this value is NaN, it is returned.
   * </ul>
   *
   * @return the absolute value of this value
   */
  DD abs() {
    if (isNaN()) return NaN;
    if (isNegative()) return negate();
    return new DD.fromDD(this);
  }

  /**
   * Computes the square of this value.
   *
   * @return the square of this value.
   */
  DD sqr() {
    return this.multiplyDD(this);
  }

  /**
   * Squares this object.
   * To prevent altering constants,
   * this method <b>must only</b> be used on values known to
   * be newly created.
   *
   * @return the square of this value.
   */
  DD selfSqr() {
    return this.selfMultiplyDD(this);
  }

  /**
   * Computes the square of this value.
   *
   * @return the square of this value.
   */
  static DD sqrDouble(double x) {
    return valueOf(x).selfMultiply(x);
  }

  /**
   * Computes the positive square root of this value.
   * If the number is NaN or negative, NaN is returned.
   *
   * @return the positive square root of this number.
   * If the argument is NaN or less than zero, the result is NaN.
   */
  DD sqrt() {
/* Strategy:  Use Karp's trick:  if x is an approximation
    to sqrt(a), then

       sqrt(a) = a*x + [a - (a*x)^2] * x / 2   (approx)

    The approximation is accurate to twice the accuracy of x.
    Also, the multiplication (a*x) and [-]*x can be done with
    only half the precision.
 */

    if (isZero()) return valueOf(0.0);

    if (isNegative()) {
      return NaN;
    }

    double x = 1.0 / math.sqrt(hi);
    double ax = hi * x;

    DD axdd = valueOf(ax);
    DD diffSq = this.subtractDD(axdd.sqr());
    double d2 = diffSq.hi * (x * 0.5);

    return axdd.add(d2);
  }

  static DD sqrtDouble(double x) {
    return valueOf(x).sqrt();
  }

  /**
   * Computes the value of this number raised to an integral power.
   * Follows semantics of Java Math.pow as closely as possible.
   *
   * @param exp the integer exponent
   * @return x raised to the integral power exp
   */
  DD pow(int exp) {
    if (exp == 0.0) return valueOf(1.0);

    DD r = new DD.fromDD(this);
    DD s = valueOf(1.0);
    int n = exp.abs();

    if (n > 1) {
/* Use binary exponentiation */
      while (n > 0) {
        if (n % 2 == 1) {
          s.selfMultiplyDD(r);
        }
        n = n ~/ 2;
        if (n > 0) r = r.sqr();
      }
    } else {
      s = r;
    }

/* Compute the reciprocal if n is negative. */
    if (exp < 0) return s.reciprocal();
    return s;
  }

  /**
   * Computes the determinant of the 2x2 matrix with the given entries.
   *
   * @param x1 a double value
   * @param y1 a double value
   * @param x2 a double value
   * @param y2 a double value
   * @return the determinant of the values
   */
  static DD determinant(double x1, double y1, double x2, double y2) {
    return determinantDD(valueOf(x1), valueOf(y1), valueOf(x2), valueOf(y2));
  }

  /**
   * Computes the determinant of the 2x2 matrix with the given entries.
   *
   * @param x1 a matrix entry
   * @param y1 a matrix entry
   * @param x2 a matrix entry
   * @param y2 a matrix entry
   * @return the determinant of the matrix of values
   */
  static DD determinantDD(DD x1, DD y1, DD x2, DD y2) {
    DD det = x1.multiplyDD(y2).selfSubtractDD(y1.multiplyDD(x2));
    return det;
  }

/*------------------------------------------------------------
   *   Ordering Functions
   *------------------------------------------------------------
   */

  /**
   * Computes the minimum of this and another DD number.
   *
   * @param x a DD number
   * @return the minimum of the two numbers
   */
  DD min(DD x) {
    if (this.le(x)) {
      return this;
    } else {
      return x;
    }
  }

  /**
   * Computes the maximum of this and another DD number.
   *
   * @param x a DD number
   * @return the maximum of the two numbers
   */
  DD max(DD x) {
    if (this.ge(x)) {
      return this;
    } else {
      return x;
    }
  }

/*------------------------------------------------------------
   *   Conversion Functions
   *------------------------------------------------------------
   */

  /**
   * Converts this value to the nearest double-precision number.
   *
   * @return the nearest double-precision number to this value
   */
  double doubleValue() {
    return hi + lo;
  }

  /**
   * Converts this value to the nearest integer.
   *
   * @return the nearest integer to this value
   */
  int intValue() {
    return hi.round();
  }

/*------------------------------------------------------------
   *   Predicates
   *------------------------------------------------------------
   */

  /**
   * Tests whether this value is equal to 0.
   *
   * @return true if this value is equal to 0
   */
  bool isZero() {
    return hi == 0.0 && lo == 0.0;
  }

  /**
   * Tests whether this value is less than 0.
   *
   * @return true if this value is less than 0
   */
  bool isNegative() {
    return hi < 0.0 || (hi == 0.0 && lo < 0.0);
  }

  /**
   * Tests whether this value is greater than 0.
   *
   * @return true if this value is greater than 0
   */
  bool isPositive() {
    return hi > 0.0 || (hi == 0.0 && lo > 0.0);
  }

  /**
   * Tests whether this value is NaN.
   *
   * @return true if this value is NaN
   */
  bool isNaN() {
    return hi.isNaN;
  }

  /**
   * Tests whether this value is equal to another <tt>DoubleDouble</tt> value.
   *
   * @param y a DoubleDouble value
   * @return true if this value = y
   */
  bool equals(DD y) {
    return hi == y.hi && lo == y.lo;
  }

  /**
   * Tests whether this value is greater than another <tt>DoubleDouble</tt> value.
   * @param y a DoubleDouble value
   * @return true if this value &gt; y
   */
  bool gt(DD y) {
    return (hi > y.hi) || (hi == y.hi && lo > y.lo);
  }

  /**
   * Tests whether this value is greater than or equals to another <tt>DoubleDouble</tt> value.
   * @param y a DoubleDouble value
   * @return true if this value &gt;= y
   */
  bool ge(DD y) {
    return (hi > y.hi) || (hi == y.hi && lo >= y.lo);
  }

  /**
   * Tests whether this value is less than another <tt>DoubleDouble</tt> value.
   * @param y a DoubleDouble value
   * @return true if this value &lt; y
   */
  bool lt(DD y) {
    return (hi < y.hi) || (hi == y.hi && lo < y.lo);
  }

  /**
   * Tests whether this value is less than or equal to another <tt>DoubleDouble</tt> value.
   * @param y a DoubleDouble value
   * @return true if this value &lt;= y
   */
  bool le(DD y) {
    return (hi < y.hi) || (hi == y.hi && lo <= y.lo);
  }

  /**
   * Compares two DoubleDouble objects numerically.
   *
   * @return -1,0 or 1 depending on whether this value is less than, equal to
   * or greater than the value of <tt>o</tt>
   */
  int compareTo(dynamic o) {
    DD other = o;

    if (hi < other.hi) return -1;
    if (hi > other.hi) return 1;
    if (lo < other.lo) return -1;
    if (lo > other.lo) return 1;
    return 0;
  }

/*------------------------------------------------------------
   *   Output
   *------------------------------------------------------------
   */

  static final int MAX_PRINT_DIGITS = 32;
  static final DD TEN = DD.valueOf(10.0);
  static final DD ONE = DD.valueOf(1.0);
  static final String SCI_NOT_EXPONENT_CHAR = "E";
  static final String SCI_NOT_ZERO = "0.0E0";

  /**
   * Dumps the components of this number to a string.
   *
   * @return a string showing the components of the number
   */
  String dump() {
    return "DD<$hi, $lo>";
  }

  /**
   * Returns a string representation of this number, in either standard or scientific notation.
   * If the magnitude of the number is in the range [ 10<sup>-3</sup>, 10<sup>8</sup> ]
   * standard notation will be used.  Otherwise, scientific notation will be used.
   *
   * @return a string representation of this number
   */
  String toString() {
    int mag = magnitude(hi);
    if (mag >= -3 && mag <= 20) return toStandardNotation();
    return toSciNotation();
  }

  /**
   * Returns the string representation of this value in standard notation.
   *
   * @return the string representation in standard notation
   */
  String toStandardNotation() {
    String? specialStr = getSpecialNumberString();
    if (specialStr != null) return specialStr;

    List<int> magnitude = List.filled(1, 0);

    String sigDigits = extractSignificantDigits(true, magnitude);
    int decimalPointPos = magnitude[0] + 1;

    String num = sigDigits;
// add a leading 0 if the decimal point is the first char
    if (sigDigits[0] == '.') {
      num = "0" + sigDigits;
    } else if (decimalPointPos < 0) {
      num = "0." + stringOfChar('0', -decimalPointPos) + sigDigits;
    } else if (sigDigits.indexOf('.') == -1) {
// no point inserted - sig digits must be smaller than magnitude of number
// add zeroes to end to make number the correct size
      int numZeroes = decimalPointPos - sigDigits.length;
      String zeroes = stringOfChar('0', numZeroes);
      num = sigDigits + zeroes + ".0";
    }

    if (this.isNegative()) return "-" + num;
    return num;
  }

  /**
   * Creates a string of a given length containing the given character
   *
   * @param ch the character to be repeated
   * @param len the len of the desired string
   * @return the string
   */
  static String stringOfChar(String ch, int len) {
    StringBuffer buf = new StringBuffer();
    for (int i = 0; i < len; i++) {
      buf.write(ch);
    }
    return buf.toString();
  }

  /**
   * Returns the string representation of this value in scientific notation.
   *
   * @return the string representation in scientific notation
   */
  String toSciNotation() {
// special case zero, to allow as
    if (isZero()) return SCI_NOT_ZERO;

    String? specialStr = getSpecialNumberString();
    if (specialStr != null) return specialStr;

    List<int> magnitudeList = List.filled(1, 0);
    String digits = extractSignificantDigits(false, magnitudeList);
    String expStr = "$SCI_NOT_EXPONENT_CHAR${magnitudeList[0]}";

// should never have leading zeroes
// MD - is this correct?  Or should we simply strip them if they are present?
    if (digits[0] == '0') {
      throw ArgumentError("Found leading zero: $digits");
    }

// add decimal point
    String trailingDigits = "";
    if (digits.length > 1) trailingDigits = digits.substring(1);
    String digitsWithDecimal = digits[0] + "." + trailingDigits;

    if (this.isNegative()) return "-" + digitsWithDecimal + expStr;
    return digitsWithDecimal + expStr;
  }

  /**
   * Extracts the significant digits in the decimal representation of the argument.
   * A decimal point may be optionally inserted in the string of digits
   * (as long as its position lies within the extracted digits
   * - if not, the caller must prepend or append the appropriate zeroes and decimal point).
   *
   * @param y the number to extract ( >= 0)
   * @param decimalPointPos the position in which to insert a decimal point
   * @return the string containing the significant digits and possibly a decimal point
   */
  String extractSignificantDigits(bool insertDecimalPoint, List<int> magnitudeList) {
    DD y = this.abs();
// compute *correct* magnitude of y
    int mag = magnitude(y.hi);
    DD scale = TEN.pow(mag);
    y = y.divideDD(scale);

// fix magnitude if off by one
    if (y.gt(TEN)) {
      y = y.divideDD(TEN);
      mag += 1;
    } else if (y.lt(ONE)) {
      y = y.multiplyDD(TEN);
      mag -= 1;
    }

    int decimalPointPos = mag + 1;
    StringBuffer buf = new StringBuffer();
    int numDigits = MAX_PRINT_DIGITS - 1;
    for (int i = 0; i <= numDigits; i++) {
      if (insertDecimalPoint && i == decimalPointPos) {
        buf.write('.');
      }
      int digit = y.hi.round();
//      System.out.println("printDump: [" + i + "] digit: " + digit + "  y: " + y.dump() + "  buf: " + buf);

      /**
       * This should never happen, due to heuristic checks on remainder below
       */
      if (digit < 0 || digit > 9) {
//        System.out.println("digit > 10 : " + digit);
//        throw new IllegalStateException("Internal errror: found digit = " + digit);
      }
      /**
       * If a negative remainder is encountered, simply terminate the extraction.
       * This is robust, but maybe slightly inaccurate.
       * My current hypothesis is that negative remainders only occur for very small lo components,
       * so the inaccuracy is tolerable
       */
      if (digit < 0) {
        break;
// throw new IllegalStateException("Internal errror: found digit = " + digit);
      }
      bool rebiasBy10 = false;
      String digitChar = String.fromCharCode(0);
      if (digit > 9) {
// set flag to re-bias after next 10-shift
        rebiasBy10 = true;
// output digit will end up being '9'
        digitChar = '9';
      } else {
        digitChar = '0$digit';
      }
      buf.write(digitChar);
      y = (y.subtractDD(DD.valueOf(digit.toDouble())).multiplyDD(TEN));
      if (rebiasBy10) y.selfAddDD(TEN);

      bool continueExtractingDigits = true;
      /**
       * Heuristic check: if the remaining portion of
       * y is non-positive, assume that output is complete
       */
//      if (y.hi <= 0.0)
//        if (y.hi < 0.0)
//        continueExtractingDigits = false;
      /**
       * Check if remaining digits will be 0, and if so don't output them.
       * Do this by comparing the magnitude of the remainder with the expected precision.
       */
      int remMag = magnitude(y.hi);
      if (remMag < 0 && remMag.abs() >= (numDigits - i)) continueExtractingDigits = false;
      if (!continueExtractingDigits) break;
    }
    magnitudeList[0] = mag;
    return buf.toString();
  }

  /**
   * Returns the string for this value if it has a known representation.
   * (E.g. NaN or 0.0)
   *
   * @return the string for this special number
   * or null if the number is not a special number
   */
  String? getSpecialNumberString() {
    if (isZero()) return "0.0";
    if (isNaN()) return "NaN ";
    return null;
  }

  /**
   * Determines the decimal magnitude of a number.
   * The magnitude is the exponent of the greatest power of 10 which is less than
   * or equal to the number.
   *
   * @param x the number to find the magnitude of
   * @return the decimal magnitude of x
   */
  static int magnitude(double x) {
    double xAbs = x.abs();
    double xLog10 = math.log(xAbs) / math.log(10);
    int xMag = xLog10.floor();
    /**
     * Since log computation is inexact, there may be an off-by-one error
     * in the computed magnitude.
     * Following tests that magnitude is correct, and adjusts it if not
     */
    double xApprox = math.pow(10, xMag).toDouble();
    if (xApprox * 10 <= xAbs) xMag += 1;

    return xMag;
  }

/*------------------------------------------------------------
   *   Input
   *------------------------------------------------------------
   */

  /**
   * Converts a string representation of a real number into a DoubleDouble value.
   * The format accepted is similar to the standard Java real number syntax.
   * It is defined by the following regular expression:
   * <pre>
   * [<tt>+</tt>|<tt>-</tt>] {<i>digit</i>} [ <tt>.</tt> {<i>digit</i>} ] [ ( <tt>e</tt> | <tt>E</tt> ) [<tt>+</tt>|<tt>-</tt>] {<i>digit</i>}+
   * </pre>
   *
   * @param str the string to parse
   * @return the value of the parsed number
   * @throws NumberFormatException if <tt>str</tt> is not a valid representation of a number
   */
  static DD parse(String str) {
    int i = 0;
    int strlen = str.length;

// skip leading whitespace
    String strTrim = str.trimLeft();

//while ( Character.isWhitespace(str.charAt(i)))
//i++;

// check for sign
    bool isNegative = false;
    if (i < strlen) {
      String signCh = strTrim[i];
      if (signCh == '-' || signCh == '+') {
        i++;
        if (signCh == '-') isNegative = true;
      }
    }

// scan all digits and accumulate into an integral value
// Keep track of the location of the decimal point (if any) to allow scaling later
    DD val = new DD.empty();

    int numDigits = 0;
    int numBeforeDec = 0;
    int exp = 0;
    bool hasDecimalChar = false;
    while (true) {
      if (i >= strlen) break;
      String ch = str[i];
      i++;
      if (StringUtils.isDigit(ch, 0)) {
        double d = (ch.codeUnitAt(0) - '0'.codeUnitAt(0)).toDouble();
        val.selfMultiplyDD(TEN);
// MD: need to optimize this
        val.selfAdd(d);
        numDigits++;
        continue;
      }
      if (ch == '.') {
        numBeforeDec = numDigits;
        hasDecimalChar = true;
        continue;
      }
      if (ch == 'e' || ch == 'E') {
        String expStr = str.substring(i);
// this should catch any format problems with the exponent
        try {
          exp = int.parse(expStr);
        } catch (ex) {
          throw new ArgumentError("Invalid exponent $expStr in string $str");
        }
        break;
      }
      throw new ArgumentError("Unexpected character '$ch' at position $i in string $str");
    }
    DD val2 = val;

// correct number of digits before decimal sign if we don't have a decimal sign in the string
    if (!hasDecimalChar) numBeforeDec = numDigits;

// scale the number correctly
    int numDecPlaces = numDigits - numBeforeDec - exp;
    if (numDecPlaces == 0) {
      val2 = val;
    } else if (numDecPlaces > 0) {
      DD scale = TEN.pow(numDecPlaces);
      val2 = val.divideDD(scale);
    } else if (numDecPlaces < 0) {
      DD scale = TEN.pow(-numDecPlaces);
      val2 = val.multiplyDD(scale);
    }
// apply leading sign, if any
    if (isNegative) {
      return val2.negate();
    }
    return val2;
  }
}

/*
 * Copyright (c) 2016 Martin Davis.
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License 2.0
 * and Eclipse Distribution License v. 1.0 which accompanies this distribution.
 * The Eclipse Public License is available at http://www.eclipse.org/legal/epl-v20.html
 * and the Eclipse Distribution License is available at
 *
 * http://www.eclipse.org/org/documents/edl-v10.php.
 */
/**
 * Implements some 2D matrix operations 
 * (in particular, solving systems of linear equations).
 * 
 * @author Martin Davis
 *
 */
class Matrix {
  static void swapRows(List<List<double>> m, int i, int j) {
    if (i == j) return;
    for (int col = 0; col < m[0].length; col++) {
      double temp = m[i][col];
      m[i][col] = m[j][col];
      m[j][col] = temp;
    }
  }

  static void swapRowsList(List<double> m, int i, int j) {
    if (i == j) return;
    double temp = m[i];
    m[i] = m[j];
    m[j] = temp;
  }

  /**
   * Solves a system of equations using Gaussian Elimination.
   * In order to avoid overhead the algorithm runs in-place
   * on A - if A should not be modified the client must supply a copy.
   * 
   * @param a an nxn matrix in row/column order )modified by this method)
   * @param b a vector of length n
   * 
   * @return a vector containing the solution (if any)
   * or null if the system has no or no unique solution
   * 
   * @throws IllegalArgumentException if the matrix is the wrong size 
   */
  static List<double?>? solve(List<List<double>> a, List<double> b) {
    int n = b.length;
    if (a.length != n || a[0].length != n) throw ArgumentError("Matrix A is incorrectly sized");

    // Use Gaussian Elimination with partial pivoting.
    // Iterate over each row
    for (int i = 0; i < n; i++) {
      // Find the largest pivot in the rows below the current one.
      int maxElementRow = i;
      for (int j = i + 1; j < n; j++) if (a[j][i].abs() > a[maxElementRow][i].abs()) maxElementRow = j;

      if (a[maxElementRow][i] == 0.0) return null;

      // Exchange current row and maxElementRow in A and b.
      swapRows(a, i, maxElementRow);
      swapRowsList(b, i, maxElementRow);

      // Eliminate using row i
      for (int j = i + 1; j < n; j++) {
        double rowFactor = a[j][i] / a[i][i];
        for (int k = n - 1; k >= i; k--) a[j][k] -= a[i][k] * rowFactor;
        b[j] -= b[i] * rowFactor;
      }
    }

    /**
     * A is now (virtually) in upper-triangular form.
     * The solution vector is determined by back-substitution.
     */
    List<double> solution = List.filled(n, 0.0);
    for (int j = n - 1; j >= 0; j--) {
      double t = 0.0;
      for (int k = j + 1; k < n; k++) t += a[j][k] * solution[k];
      solution[j] = (b[j] - t) / a[j][j];
    }
    return solution;
  }
}
