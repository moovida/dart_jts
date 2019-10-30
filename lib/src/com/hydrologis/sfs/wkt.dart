part of dart_sfs;

/// Parses the WKT string [wkt] and replies the parsed [Geometry].
///
/// [wkt] must not be null.
///
/// Throws a [WKTError] if parsing fails.
Geometry parseWKT(String wkt) {
  _require(wkt != null);
  var parser = _WKTParser(wkt);
  parser.advanceMandatory();
  return parser.parseWKTObject();
}

class _WKTTokenType {
  static const QUOTED_NAME = _WKTTokenType._(0, "QUOTED_NAME");
  static const KEYWORD = _WKTTokenType._(1, "KEYWORD");
  static const NUMERIC_LITERAL = _WKTTokenType._(2, "NUMERIC_LITERAL");
  static const WHITESPACE = _WKTTokenType._(3, "WHITESPACE");
  static const LPAREN = _WKTTokenType._(4, "LPAREN");
  static const RPAREN = _WKTTokenType._(5, "RPAREN");
  static const LBRACKET = _WKTTokenType._(6, "LBRACKET");
  static const RBRACKET = _WKTTokenType._(7, "RBRACKET");
  static const ERROR = _WKTTokenType._(8, "ERROR");
  static const EOS = _WKTTokenType._(9, "EOS");
  static const COMMA = _WKTTokenType._(10, "COMMA");

  final _type;
  final _displayName;

  const _WKTTokenType._(this._type, this._displayName);

  String toString() => _displayName;
}

class _WKTToken {
  final type;
  final value;

  const _WKTToken(this.type, this.value);

  String toString() => "{${type}: $value}";

  bool get isKeyword => type == _WKTTokenType.KEYWORD;

  bool get isLParen => type == _WKTTokenType.LPAREN;

  bool get isRParen => type == _WKTTokenType.RPAREN;

  bool get isComma => type == _WKTTokenType.COMMA;

  bool get isNumber => type == _WKTTokenType.NUMERIC_LITERAL;

  bool get isEOS => type == _WKTTokenType.EOS;

  ensureKeyword([kw]) {
    if (!isKeyword) {
      throw WKTError("expected a keyword, got '$value'");
    }
    if (kw != null && value.toLowerCase() != kw.toLowerCase()) {
      throw WKTError("unexpected keyword: expected '$kw', got '$value'");
    }
  }

  ensureLParen() {
    if (!isLParen) throw WKTError("expected left paren '(', got '$value'");
  }

  ensureRParen() {
    if (!isRParen) throw WKTError("expected right paren ')', got '$value'");
  }

  ensureNumber() {
    if (!isNumber) throw WKTError("expected numeric literal, got '$value'");
  }

  ensureNotEOS() {
    if (isEOS) throw WKTError("unexpected end of input");
  }

  bool matchKeyword(kw) => isKeyword && value.toLowerCase() == kw.toLowerCase();
}

class _WKTTokenizer {
  static final CC_0 = "0".codeUnitAt(0);
  static final CC_9 = "9".codeUnitAt(0);
  static final CC_a = "a".codeUnitAt(0);
  static final CC_z = "z".codeUnitAt(0);
  static final CC_A = "A".codeUnitAt(0);
  static final CC_Z = "Z".codeUnitAt(0);
  static final CC_SPACE = " ".codeUnitAt(0);
  static final CC_CR = "\n".codeUnitAt(0);
  static final CC_LF = "\r".codeUnitAt(0);
  static final CC_LPAREN = "(".codeUnitAt(0);
  static final CC_RPAREN = ")".codeUnitAt(0);
  static final CC_DQUOTE = "\"".codeUnitAt(0);
  static final CC_QUOTE = "'".codeUnitAt(0);
  static final CC_MINUS = "-".codeUnitAt(0);
  static final CC_PLUS = "+".codeUnitAt(0);
  static final CC_UNDERSCORE = "_".codeUnitAt(0);
  static final CC_PERIOD = ".".codeUnitAt(0);
  static final CC_COMMA = ",".codeUnitAt(0);
  static final CC_DOT = ".".codeUnitAt(0);
  static final CC_E = "E".codeUnitAt(0);
  static final CC_e = "e".codeUnitAt(0);

  static bool isDigit(int c) => c >= CC_0 && c <= CC_9;

  static bool isLowerCaseLetter(int c) => c >= CC_a && c <= CC_z;

  static bool isUpperCaseLetter(int c) => c >= CC_A && c <= CC_Z;

  static bool isLetter(int c) => isLowerCaseLetter(c) || isUpperCaseLetter(c);

  static bool isWS(int c) => c == CC_SPACE || c == CC_CR || c == CC_LF;

  static bool isParen(int c) => c == CC_LPAREN || c == CC_RPAREN;

  static bool isDQuote(int c) => c == CC_DQUOTE;

  static bool isSpecial(int c) =>
      c == CC_RPAREN ||
      c == CC_LPAREN ||
      c == CC_MINUS ||
      c == CC_UNDERSCORE ||
      c == CC_PERIOD ||
      c == CC_QUOTE ||
      c == CC_SPACE;

  static bool isComma(int c) => c == CC_COMMA;

  static bool isSign(int c) => c == CC_MINUS || c == CC_PLUS;

  static bool isDot(int c) => c == CC_DOT;

  static bool isE(int c) => c == CC_E || c == CC_e;

  static bool isDelimiter(int c) =>
      isWS(c) || isParen(c) || isComma(c) || isDQuote(c) || isDot(c);

  static bool isLParen(int c) => c == CC_LPAREN;

  static bool isRParen(int c) => c == CC_RPAREN;

  static errorToken(value) => _WKTToken(_WKTTokenType.ERROR, value);

  static numericLiteralToken(value) =>
      _WKTToken(_WKTTokenType.NUMERIC_LITERAL, value);

  static eosToken() => _WKTToken(_WKTTokenType.EOS, null);

  static wsToken(value) => _WKTToken(_WKTTokenType.WHITESPACE, value);

  static keywordToken(value) => _WKTToken(_WKTTokenType.KEYWORD, value);

  static commaToken(value) => _WKTToken(_WKTTokenType.COMMA, value);

  static quotedNameToken(value) => _WKTToken(_WKTTokenType.QUOTED_NAME, value);

  Iterable<int> source;
  int head = 0;
  bool eos = false;
  bool skipWhitespace = true;

  _WKTTokenizer(source, {this.skipWhitespace = true}) {
    if (source is Iterable<int>) {
      this.source = source;
    } else if (source is String) {
      this.source = source.codeUnits;
    }
    try {
      this.source.elementAt(0);
    } on RangeError {
      eos = true;
    }
  }

  advance() {
    head++;
    try {
      source.elementAt(head);
      return true;
    } on RangeError {
      eos = true;
      return false;
    }
  }

  int get cur => source.elementAt(head);

  currentToken() => String.fromCharCodes(source.take(head));

  consumeToken() {
    var token = currentToken();
    source = source.skip(head);
    head = 0;
    return token;
  }

  consumeWhiteSpace() {
    assert(isWS(cur));
    while (advance()) {
      if (!isWS(cur)) break;
    }
    var token = consumeToken();
    return wsToken(token);
  }

  consumeComma() {
    assert(isComma(cur));
    advance();
    var token = consumeToken();
    return commaToken(token);
  }

  consumeKeyword() {
    assert(isLetter(cur));
    while (advance()) {
      var c = cur;
      if (isLetter(c)) continue;
      if (isWS(c)) break;
      if (isParen(c)) break;
      if (isComma(c)) break;
      return errorToken(
          "unexpected character in keyword at <${currentToken()}>");
    }
    var token = consumeToken();
    return keywordToken(token);
  }

  consumeQuotedName() {
    assert(isDQuote(cur));
    while (advance()) {
      var c = cur;
      if (isWS(c) || isLetter(c) || isDigit(c) || isSpecial(c)) continue;
      if (isDQuote(c)) {
        advance();
        var token = consumeToken();
        return quotedNameToken(token);
      }
      return errorToken(
          "unexpected character in quoted name at  <${currentToken()}>");
    }
    if (eos) {
      return errorToken("quoted name not terminated at <${currentToken()}>");
    }
  }

  scanUnsignedInteger() {
    assert(isDigit(cur));
    while (advance()) {
      if (!isDigit(cur)) break;
    }
  }

  scanSignedInteger() {
    assert(isDigit(cur) || isSign(cur));
    if (isSign(cur)) {
      if (!advance()) throw "unexpected end of stream";
      if (!isDigit(cur)) throw "unexpected character after sign";
    }
    scanUnsignedInteger();
  }

  scanExactNumericLiteral() {
    if (isDot(cur)) {
      if (!advance()) return;
      scanUnsignedInteger();
    } else if (isDigit(cur)) {
      scanUnsignedInteger();
      if (eos) return;
      if (isDot(cur)) {
        if (!advance()) return;
        scanUnsignedInteger();
      }
    }
  }

  scanUnsignedNumericLiteral() {
    scanExactNumericLiteral();
    if (eos) return;
    if (isE(cur)) {
      if (!advance()) throw "unexpected end of stream";
      if (!isDigit(cur) && !isSign(cur)) throw "unexpected character after 'E'";
      scanSignedInteger();
    }
  }

  scanSignedNumericLiteral() {
    if (isSign(cur)) {
      if (!advance()) throw "unexpected end of stream";
      if (!(isDot(cur) || isDigit(cur))) {
        throw "unexpected character after sign";
      }
    }
    scanUnsignedNumericLiteral();
  }

  consumeSignedNumericLiteral() {
    try {
      scanSignedNumericLiteral();
      var token = consumeToken();
      return numericLiteralToken(token);
    } catch (e, st) {
      print(st);
      return errorToken(e.toString());
    }
  }

  consumeLParen() {
    assert(isLParen(cur));
    advance();
    var token = consumeToken();
    return _WKTToken(_WKTTokenType.LPAREN, "(");
  }

  consumeRParen() {
    assert(isRParen(cur));
    advance();
    var token = consumeToken();
    return _WKTToken(_WKTTokenType.RPAREN, ")");
  }

  _WKTToken next() {
    if (eos) return eosToken();
    while (!eos) {
      var c = cur;
      if (isWS(c)) {
        var token = consumeWhiteSpace();
        if (!skipWhitespace) return token;
      } else if (isLetter(c)) {
        return consumeKeyword();
      } else if (isDQuote(c)) {
        return consumeQuotedName();
      } else if (isComma(c)) {
        return consumeComma();
      } else if (isLParen(c)) {
        return consumeLParen();
      } else if (isRParen(c)) {
        return consumeRParen();
      } else if (isDigit(c) || isDot(c) || isSign(c)) {
        return consumeSignedNumericLiteral();
      } else {
        errorToken("unexpected character at <${currentToken()}>");
      }
    }
  }
}

class WKTError implements Error {
  final String message;

  WKTError([this.message]);

  String toString() => "WKT format error" + (message != null ? message : "");

  @override
  StackTrace get stackTrace => StackTrace.fromString(toString());
}

class _CoordSpecification {
  bool withZ = false;
  bool withM = false;
  int ncoord = 2;

  _CoordSpecification({this.withZ: false, this.withM: false}) {
    ncoord = 2;
    if (withZ) ncoord++;
    if (withM) ncoord++;
  }

  _CoordSpecification.xy() : this(withZ: false, withM: false);

  _CoordSpecification.from(String s) {
    s = s.toLowerCase();
    if (s == "z") {
      withZ = true;
      ncoord = 3;
    } else if (s == "m") {
      withM = true;
      ncoord = 3;
    } else if (s == "zm" || s == "mz") {
      withM = true;
      withZ = true;
      ncoord = 4;
    } else {
      throw "unexpected specification for coordinate dimenisons, got '$s'";
    }
  }

  bool get withMZ => withZ && withM;
}

class _WKTParser {
  var tokenizer;
  var token;

  _WKTParser(source) {
    tokenizer = _WKTTokenizer(source);
  }

  parseNumbers(n) {
    decode(s) {
      try {
        return int.parse(s);
      } catch (e) {
        try {
          return double.parse(s);
        } catch (e) {
          throw WKTError("failed to parse numeric literal '${s}'");
        }
      }
    }

    var ret = [];
    for (int i = 0; i < n; i++) {
      token.ensureNotEOS();
      token.ensureNumber();
      ret.add(decode(token.value));
      advanceMandatory();
    }
    return ret;
  }

  parseCoordSpecificationOrEmpty() {
    var coordSpec;
    if (token.matchKeyword("empty")) {
      // return null
    } else if (token.isLParen) {
      coordSpec = _CoordSpecification.xy();
    } else {
      coordSpec = parseCoordSpecification();
      advanceMandatory();
      if (token.matchKeyword("empty")) {
        return null;
      } else {
        token.ensureLParen();
      }
    }
    return coordSpec;
  }

  /// pre: token is 'point'
  /// post: token is ')'
  parsePoint() {
    token.ensureKeyword("point");
    advanceMandatory();
    var coordSpec = parseCoordSpecificationOrEmpty();
    if (coordSpec == null) return Point.empty();
    advanceMandatory();
    var point = parsePointText(coordSpec);
    token.ensureRParen();
    return point;
  }

  advance() => token = tokenizer.next();

  ensureKeyword(kw) => token.ensureKeyword(kw);

  advanceMandatory() {
    advance();
    token.ensureNotEOS();
  }

  parseCoordSpecification() {
    try {
      return _CoordSpecification.from(token.value);
    } catch (e) {
      throw WKTError(e.toString());
    }
  }

  /// pre: token is a number
  /// post: token is the token following the coordSpec.ncoord parsed numbers
  parsePointText(coordSpec) {
    var numbers = parseNumbers(coordSpec.ncoord);
    var x = numbers[0];
    var y = numbers[1];
    var z = coordSpec.withZ ? numbers[2] : null;
    var m;
    if (coordSpec.withMZ) {
      m = numbers[3];
    } else if (coordSpec.withM) {
      m = numbers[2];
    } else {
      m = null;
    }
    return Point(x, y, z: z, m: m);
  }

  /// pre: token is 'multipoint'
  /// post: token is ')'
  parseMultiPoint() {
    token.ensureKeyword("multipoint");
    advanceMandatory();
    var coordSpec = parseCoordSpecificationOrEmpty();
    if (coordSpec == null) return new MultiPoint.empty();
    advanceMandatory();
    var points = [];
    while (true) {
      if (token.isRParen) {
        break;
      } else if (token.isLParen) {
        advanceMandatory();
        points.add(parsePointText(coordSpec));
        token.ensureRParen();
        advanceMandatory();
        if (token.isComma) {
          advanceMandatory();
        } else if (token.isRParen) {
          break;
        } else {
          throw WKTError("expected ',' or ')', got '${token.value}'");
        }
      } else {
        throw WKTError("expected '(' or ')', got '${token.value}'");
      }
    }
    return MultiPoint(points);
  }

  /// pre: token is 'linestring'
  /// post: token is ')'
  parseLineString() {
    token.ensureKeyword("linestring");
    advanceMandatory();
    var coordSpec = parseCoordSpecificationOrEmpty();
    if (coordSpec == null) return LineString.empty();
    List<Point> points = parseLineStringText(coordSpec);
    return LineString(points);
  }

  /// pre: token  is  '('
  /// post: token is  ')'
  parseLineStringText(coordSpec) {
    List<Point> points = [];
    token.ensureLParen();
    while (true) {
      advanceMandatory();
      points.add(parsePointText(coordSpec));
      if (token.isComma) {
        continue;
      } else if (token.isRParen) {
        break;
      } else {
        throw WKTError("expected ',' or ')', got '${token.value}'");
      }
    }
    return points;
  }

  parseCommaSeperatedList(parse()) {
    List<Geometry> elements = [];
    while (true) {
      elements.add(parse());
      advanceMandatory();
      if (token.isComma) {
        advanceMandatory();
        continue;
      } else if (token.isRParen) {
        break;
      } else {
        throw WKTError("expected ',' or ')', got '${token.value}'");
      }
    }
    return elements;
  }

  /// pre: token is 'multilinestring'
  /// post: token is ')'
  parseMultiLineString() {
    token.ensureKeyword("multilinestring");
    advanceMandatory();
    var coordSpec = parseCoordSpecificationOrEmpty();
    if (coordSpec == null) return MultiLineString.empty();
    advanceMandatory();
    var linestrings = parseCommaSeperatedList(() {
      if (token.matchKeyword("empty")) {
        return LineString.empty();
      } else {
        return LineString(parseLineStringText(coordSpec));
      }
    });
    return MultiLineString(linestrings);
  }

  ///  pre: token is '('  (...) (...)
  ///                 ^
  ///  post: token is '('  (...) (...) ')'
  ///                                   ^
  parsePolygonText(coordSpec) {
    token.ensureLParen();
    advanceMandatory();
    var rings = parseCommaSeperatedList(() {
      return LineString(parseLineStringText(coordSpec));
    });
    return rings;
  }

  parsePolygon() {
    token.ensureKeyword("polygon");
    advanceMandatory();
    var coordSpec = parseCoordSpecificationOrEmpty();
    if (coordSpec == null) return Polygon.empty();
    var rings = parsePolygonText(coordSpec);
    return Polygon(rings[0], rings.skip(1));
  }

  parsePolygonList(coordSpec) {
    var polygons = [];
    while (true) {
      var rings = parsePolygonText(coordSpec);
      polygons.add(Polygon(rings[0], rings.skip(1)));
      advanceMandatory();
      if (token.isComma) {
        advanceMandatory();
        continue;
      } else if (token.isRParen) {
        break;
      } else {
        throw WKTError("expected ',' or ')', got '${token.value}'");
      }
    }
    return polygons;
  }

  /// pre: token is 'multipolygon'
  parseMultiPolygon() {
    token.ensureKeyword("multipolygon");
    advanceMandatory();
    var coordSpec = parseCoordSpecificationOrEmpty();
    if (coordSpec == null) return MultiPolygon.empty();
    token.ensureLParen();
    advanceMandatory();
    var polygons = parsePolygonList(coordSpec);
    return MultiPolygon(polygons);
  }

  /// pre: token is 'polyhedralsurface'
  parsePolyhedralSurface() {
    token.ensureKeyword("polyhedralsurface");
    advanceMandatory();
    var coordSpec = parseCoordSpecificationOrEmpty();
    if (coordSpec == null) return PolyhedralSurface.empty();
    token.ensureLParen();
    advanceMandatory();
    var polygons = parsePolygonList(coordSpec);
    return PolyhedralSurface(polygons);
  }

  /// pre: token is 'tin'
  parseTin() {
    token.ensureKeyword("tin");
    advanceMandatory();
    var coordSpec = parseCoordSpecificationOrEmpty();
    if (coordSpec == null) return Tin.empty();
    token.ensureLParen();
    advanceMandatory();
    var polygons = parsePolygonList(coordSpec);
    return Tin(polygons);
  }

  /// pre: token is one of the keywords
  parseWKTObject() {
    token.ensureKeyword();
    switch (token.value.toLowerCase()) {
      case "point":
        return parsePoint();
      case "multipoint":
        return parseMultiPoint();
      case "linestring":
        return parseLineString();
      case "multilinestring":
        return parseMultiLineString();
      case "polygon":
        return parsePolygon();
      case "multipolygon":
        return parseMultiPolygon();
      case "geometrycollection":
        return parseGeometryCollection();
      case "polyhedralsurface":
        return parsePolyhedralSurface();
      case "tin":
        return parseTin();
      default:
        throw WKTError("unsupported keyword '${token.value}'"
            " in WKT string");
    }
  }

  /// pre: token is 'geometrycollection'
  parseGeometryCollection() {
    token.ensureKeyword("geometrycollection");
    advanceMandatory();
    if (token.matchKeyword("empty")) {
      return GeometryCollection.empty();
    } else {
      token.ensureLParen();
    }
    advanceMandatory();
    token.ensureKeyword();
    List<Geometry> geometries = parseCommaSeperatedList(() {
      return parseWKTObject();
    });
    return GeometryCollection(geometries);
  }
}

class _WKTWriter {
  final StringSink _sink;
  int _ident = 0;

  _WKTWriter(this._sink);

  lparen() => _sink.write("(");

  rparen() => _sink.write(")");

  comma() => _sink.write(",");

  blank() => _sink.write(" ");

  newline() => _sink.write("\n");

  position(num pos) {
    if (pos.toInt() == pos) {
      _sink.write(pos.toInt());
    } else {
      _sink.write(pos);
    }
  }

  incIdent() => _ident++;

  decIdent() {
    _ident--;
    if (_ident < 0) _ident = 0;
  }

  ident() {
    for (int i = 0; i < _ident; i++) _sink.write("  ");
  }

  coordinates(Point point) {
    position(point.x);
    blank();
    position(point.y);
    if (point.is3D) {
      blank();
      position(point.z);
    }
    if (point.isMeasured) {
      blank();
      position(point.m);
    }
  }

  ordinateSpecification({bool withZ: false, bool withM: false}) {
    var v = "";
    if (withZ) v += "Z";
    if (withM) v += "M";
    if (!v.isEmpty) {
      _sink.write(v);
      blank();
    }
  }

  empty() => _sink.write("EMPTY");

  write(String s) => _sink.write(s);
}
