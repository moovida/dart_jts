part of dart_sfs;

_require(cond, [msg]) {
  if (!cond) throw new ArgumentError(msg);
}

//abstract class ComparableMixin<T> implements Comparable {
//  bool operator ==(T other) => compareTo(other) == 0;
//  bool operator <(T other) => compareTo(other) == -1;
//  bool operator <=(T other) => compareTo(other) <= 0;
//  bool operator >(T other) => compareTo(other) == 1;
//  bool operator >=(T other) => compareTo(other) > 0;
//  int compareTo(T other);
//}

class FormattingUtils {}

/// Utilities for collections
class CollectionsUtils {
  /// Shift a list of items from a given [index] to the first position.
  static List<T> shiftToFirst<T>(List<T> list, int index) {
    if (list == null || list.isEmpty) return list;
    if (index > list.length - 1) throw ArgumentError("The shift index can't be bigger than the list size.");
    return list.sublist(index)..addAll(list.sublist(0, index));
  }

  /// Removes subsequent equal items from a [list].
  static List<T> removeRepeated<T>(List<T> list) {
    List<T> newList = [];
    for (int i = 1; i < list.length; i++) {
      if (list[i - 1] != list[i]) {
        newList.add(list[i - 1]);
      }
    }
    if (newList.last != list.last) {
      newList.add(list.last);
    }
    return newList;
  }

  /// CHecks if there are subsequent repeating items.
  static bool hasRepeated<T>(List<T> list) {
    for (int i = 1; i < list.length; i++) {
      if (list[i - 1] == list[i]) {
        return true;
      }
    }
    return false;
  }

  static bool addIfNotEqualToLast<T>(List<T> list, T item) {
    if (list.isEmpty || list.last == item) {
      list.add(item);
    }
  }

  static bool areEqual<T>(List<T> listA, List<T> listB) {
    Function eq = const ListEquality().equals;
    return eq(listA, listB);
  }
}

class StringUtils {

  static bool isDigit(String s, int idx) => (s.codeUnitAt(idx) ^ 0x30) <= 9;

  static bool equalsIgnoreCase(String string1, String string2) {
    return string1?.toLowerCase() == string2?.toLowerCase();
  }

  static String replaceCharAt(String oldString, int index, String newChar) {
    return oldString.substring(0, index) + newChar + oldString.substring(index + 1);
  }
}

class MatrixUtils {
  static List<List<int>> createIntMatrix(int rows, int cols) {
    final grid = List<List<int>>.generate(rows, (i) => List<int>.generate(cols, (j) => i * cols + j));
    return grid;
  }

  static List<List<T>> createMatrix<T>(int rows, int cols, T defaultValue) {
    final grid = List<List<T>>.generate(rows, (i) => List<T>.generate(cols, (j) => defaultValue));
    return grid;
  }
}

class MathUtils {
  /// Clamps a <tt>double</tt> value to a given range.
  /// @param x the value to clamp
  /// @param min the minimum value of the range
  /// @param max the maximum value of the range
  /// @return the clamped value
  static num clamp(num x, num min, num max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
  }

  static final double LOG_10 = math.log(10);

  /// Computes the base-10 logarithm of a <tt>double</tt> value.
  /// <ul>
  /// <li>If the argument is NaN or less than zero, then the result is NaN.
  /// <li>If the argument is positive infinity, then the result is positive infinity.
  /// <li>If the argument is positive zero or negative zero, then the result is negative infinity.
  /// </ul>
  ///
  /// @param x a positive number
  /// @return the value log a, the base-10 logarithm of the input value
  static double log10(double x) {
    double ln = math.log(x);
    if (ln.isInfinite) return ln;
    if (ln.isNaN) return ln;
    return ln / LOG_10;
  }

  /// Computes an index which wraps around a given maximum value.
  /// For values &gt;= 0, this is equals to <tt>val % max</tt>.
  /// For values &lt; 0, this is equal to <tt>max - (-val) % max</tt>
  ///
  /// @param index the value to wrap
  /// @param max the maximum value (or modulus)
  /// @return the wrapped index
  static int wrap(int index, int max) {
    if (index < 0) {
      return max - ((-index) % max);
    }
    return index % max;
  }

  /// Computes the average of two numbers.
  ///
  /// @param x1 a number
  /// @param x2 a number
  /// @return the average of the inputs
  static double average(double x1, double x2) {
    return (x1 + x2) / 2.0;
  }

  static double max(double v1, double v2, double v3) {
    double max = v1;
    if (v2 > max) max = v2;
    if (v3 > max) max = v3;
    return max;
  }

  static double max4(double v1, double v2, double v3, double v4) {
    double max = v1;
    if (v2 > max) max = v2;
    if (v3 > max) max = v3;
    if (v4 > max) max = v4;
    return max;
  }

  static double min(double v1, double v2, double v3, double v4) {
    double min = v1;
    if (v2 < min) min = v2;
    if (v3 < min) min = v3;
    if (v4 < min) min = v4;
    return min;
  }
}

class NumberUtils {
  static bool equalsWithTolerance(double x1, double x2, double tolerance) {
    return (x1 - x2).abs() <= tolerance;
  }
}
