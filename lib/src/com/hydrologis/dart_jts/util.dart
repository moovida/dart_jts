part of dart_jts;

class RuntimeException implements Exception {
  String msg;

  RuntimeException(this.msg);

  @override
  String toString() => 'RuntimeException: ' + msg;
}

class IOException implements Exception {
  String msg;

  IOException(this.msg);

  @override
  String toString() => 'IOException: ' + msg;
}

class ParseException implements Exception {
  String msg;

  ParseException(this.msg);

  @override
  String toString() => 'ParseException: ' + msg;
}

class FormattingUtils {}

/// Utilities for collections
class CollectionsUtils {
  /// Shift a list of items from a given [index] to the first position.
  static List<T>? shiftToFirst<T>(List<T>? list, int index) {
    if (list == null || list.isEmpty) return list;
    if (index > list.length - 1)
      throw ArgumentError(
          "The shift index can't be bigger than the list size.");
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
    for (var i = 1; i < list.length; i++) {
      if (list[i - 1] == list[i]) {
        return true;
      }
    }
    return false;
  }

  static bool addIfNotEqualToLast<T>(List<T> list, T item) {
    if (list.isEmpty || list.last == item) {
      list.add(item);
      return true;
    }
    return false;
  }

  static bool areEqual<T>(List<T> listA, List<T> listB) {
    Function eq = const ListEquality().equals;
    return eq(listA, listB);
  }
}

class StringUtils {
  static bool isDigit(String s, int idx) => (s.codeUnitAt(idx) ^ 0x30) <= 9;

  static bool equalsIgnoreCase(String string1, String string2) {
    return string1.toLowerCase() == string2.toLowerCase();
  }

  static String replaceCharAt(String oldString, int index, String newChar) {
    return oldString.substring(0, index) +
        newChar +
        oldString.substring(index + 1);
  }
}

class MatrixUtils {
  static List<List<int>> createIntMatrix(int rows, int cols) {
    final grid = List<List<int>>.generate(
        rows, (i) => List<int>.generate(cols, (j) => i * cols + j));
    return grid;
  }

  static List<List<int>> createIntMatrixWithDefault(
      int rows, int cols, int defaultValue) {
    final grid = List<List<int>>.generate(
        rows, (i) => List<int>.generate(cols, (j) => defaultValue));
    return grid;
  }

  static List<List<T?>> createMatrix<T>(int rows, int cols, T? defaultValue) {
    final grid = List<List<T?>>.generate(
        rows, (i) => List<T?>.generate(cols, (j) => defaultValue));
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
  static final int MAX_INT = 2E+53.toInt();
  static final int MIN_INT = -2E+53.toInt();

  static bool equalsWithTolerance(double x1, double x2, double tolerance) {
    return (x1 - x2).abs() <= tolerance;
  }
}

/// Lists used for data convertion (alias each other).
final Uint8List _unisgnedByteBuffer = Uint8List(8);
final Float64List _unsignedFloat64Buffer =
    Float64List.view(_unisgnedByteBuffer.buffer);
final Int8List _signedBytesBuffer = Int8List(8);
final Float64List _signedFloat64Buffer =
    Float64List.view(_signedBytesBuffer.buffer);

class Byteutils {
  static int doubleToLongBits(double value) {
    const pow52 = 4503599627370496.0; // 2^52
    const pow1022 = 4.49423283715579e+307; // 2^1022
    if (value.isNaN) {
      return 0x7FF8000000000000;
    }
    int signbit = 0;
    if (value.isNegative) {
      signbit = 0x8000000000000000;
      value = -value;
    }
    if (value.isInfinite) {
      return signbit | 0x7FF0000000000000;
    } else if (value < 2.2250738585072014e-308) {
      // Denormal or zero.
      // Multiply by 2^(1022+52) to get the bits into the correct position.
      int bits = (value * pow1022 * pow52).toInt();
      return signbit | bits;
    } else {
      // Slow linear search to move bits into correct position for mantissa.
      // Use binary search or something even smarter for speed.
      int exponent = 52;
      while (value < pow52) {
        value *= 2;
        exponent -= 1;
      }
      while (value >= pow52 * 2) {
        value /= 2;
        exponent += 1;
      }
      int mantissaBits = (value - pow52).toInt();
      int exponentBits = (exponent + 1023);
      return signbit | (exponentBits << 52) | mantissaBits;
    }
  }

  static double getFloat64(List<int> buf, Endian byteOrder, {position = 0}) {
    if (byteOrder == Endian.little) {
      _unisgnedByteBuffer[0] = buf[position + 7];
      _unisgnedByteBuffer[1] = buf[position + 6];
      _unisgnedByteBuffer[2] = buf[position + 5];
      _unisgnedByteBuffer[3] = buf[position + 4];
      _unisgnedByteBuffer[4] = buf[position + 3];
      _unisgnedByteBuffer[5] = buf[position + 2];
      _unisgnedByteBuffer[6] = buf[position + 1];
      _unisgnedByteBuffer[7] = buf[position + 0];
      return _signedFloat64Buffer[0];
    } else {
      _unisgnedByteBuffer[0] = buf[position + 0];
      _unisgnedByteBuffer[1] = buf[position + 1];
      _unisgnedByteBuffer[2] = buf[position + 2];
      _unisgnedByteBuffer[3] = buf[position + 3];
      _unisgnedByteBuffer[4] = buf[position + 4];
      _unisgnedByteBuffer[5] = buf[position + 5];
      _unisgnedByteBuffer[6] = buf[position + 6];
      _unisgnedByteBuffer[7] = buf[position + 7];
      return _signedFloat64Buffer[0];
    }
  }

  static void putFloat64(double value, List<int> buf, Endian byteOrder,
      {position = 0}) {
    ByteData bdata = ByteData.view(_signedBytesBuffer.buffer);
    bdata.setFloat64(0, value, byteOrder);
    buf.setAll(0, _signedBytesBuffer);
  }

  static int getInt32(List<int> buf, Endian byteOrder) {
    if (byteOrder == Endian.big) {
      return ((buf[0] & 0xff) << 24) |
          ((buf[1] & 0xff) << 16) |
          ((buf[2] & 0xff) << 8) |
          ((buf[3] & 0xff));
    } else {
      // LITTLE_ENDIAN
      return ((buf[3] & 0xff) << 24) |
          ((buf[2] & 0xff) << 16) |
          ((buf[1] & 0xff) << 8) |
          ((buf[0] & 0xff));
    }
  }

  static void putInt32(int intValue, List<int> buf, Endian byteOrder) {
    if (byteOrder == Endian.big) {
      buf[0] = (intValue >> 24) & 0xff;
      buf[1] = (intValue >> 16) & 0xff;
      buf[2] = (intValue >> 8) & 0xff;
      buf[3] = intValue & 0xff;
    } else {
      // LITTLE_ENDIAN
      buf[0] = intValue & 0xff;
      buf[1] = (intValue >> 8) & 0xff;
      buf[2] = (intValue >> 16) & 0xff;
      buf[3] = (intValue >> 24) & 0xff;
    }
  }
}
