import 'package:dart_jts/dart_jts.dart';
import 'package:flutter/material.dart';

void main() {
  runApp(MyApp());
}

class MyApp extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return MaterialApp(
      title: 'dart_jts example',
      home: ExamplePage(),
    );
  }
}

class ExamplePage extends StatefulWidget {
  @override
  _ExamplePageState createState() => _ExamplePageState();
}

class _ExamplePageState extends State<ExamplePage> {
  double buffer = 50;

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      appBar: AppBar(
        title: Text('Dart_JTS Example - Polygon Offset'),
      ),
      body: Center(
        child: ExampleShape(
          buffer: buffer,
        ),
      ),
      bottomNavigationBar: BottomAppBar(
          child: SizedBox(
        height: 50,
        child: Slider(
          min: 0,
          max: 100,
          onChanged: (val) {
            setState(() {
              buffer = val;
            });
          },
          value: buffer,
          divisions: 100,
          label: "${buffer.floor()}",
        ),
      )),
    );
  }
}

class ExampleShape extends StatelessWidget {
  final double buffer;
  ExampleShape({this.buffer});

  @override
  Widget build(BuildContext context) {
    return CustomPaint(
        painter: ExampleShapePainter(buffer: buffer, shape: [
      Offset(-90, -200),
      Offset(120, -25),
      Offset(90, 110),
      Offset(0, 0),
      Offset(-120, 25),
    ]));
  }
}

class ExampleShapePainter extends CustomPainter {
  final double buffer;
  final List<Offset> shape;

  ExampleShapePainter({this.buffer, this.shape})
      : assert(shape.length == 0 || shape.length >= 4);

  @override
  void paint(Canvas canvas, Size size) {
    Paint bodyPaint = Paint()
      ..color = Colors.blue
      ..style = PaintingStyle.fill;

    Paint inflatedPaint = Paint()
      ..color = Colors.orange
      ..style = PaintingStyle.fill;

    /// Draw inflated shape
    if (shape.length >= 4) {
      GeometryFactory geometryFactory = GeometryFactory.defaultPrecision();

      /// Convert [List<Offset>] to [List<Coordinate>] to work with dart_jts
      List<Coordinate> coords =
          shape.map((element) => Coordinate(element.dx, element.dy)).toList();
      coords.add(coords.first);

      /// First and last coord MUST be identical
      LinearRing linearRing = LinearRing.withFactory(coords, geometryFactory);

      /// Create the polygon with no holes (hence the empty list as the 2nd argument)
      Polygon polygon = geometryFactory.createPolygon(linearRing, []);

      polygon = polygon.buffer(buffer);
      List<Offset> inflated = polygon.shell.points
          .toCoordinateArray()
          .map((e) => Offset(e.x, e.y))
          .toList();

      Path inflatedPath = Path()..addPolygon(inflated, true);

      debugPrint(inflated.toString());
      canvas.drawPath(inflatedPath, inflatedPaint);
    }

    /// Draw normal shape
    Path path = Path();
    path.addPolygon(shape, true);
    canvas.drawPath(path, bodyPaint);
  }

  @override
  bool shouldRepaint(CustomPainter oldDelegate) => true;
}
