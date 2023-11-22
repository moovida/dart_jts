import 'package:dart_jts/dart_jts.dart';
import 'package:test/test.dart';
import "dart:math" as math;

void main() {
  final int NUM_ITEMS = 2000;
  final double MIN_EXTENT = -1000.0;
  final double MAX_EXTENT = 1000.0;

  EnvelopeList envList = EnvelopeList();
  Quadtree q = Quadtree();
  void createGrid(int nGridCells) {
    int gridSize = math.sqrt(nGridCells.toDouble()).toInt();
    gridSize += 1;
    double extent = MAX_EXTENT - MIN_EXTENT;
    double gridInc = extent / gridSize;
    double cellSize = 2 * gridInc;

    for (int i = 0; i < gridSize; i++) {
      for (int j = 0; j < gridSize; j++) {
        double x = MIN_EXTENT + gridInc * i;
        double y = MIN_EXTENT + gridInc * j;
        Envelope env = new Envelope(x, x + cellSize, y, y + cellSize);
        q.insert(env, env);
        envList.add(env);
      }
    }
  }

  void fill() {
    createGrid(NUM_ITEMS);
  }

  List getOverlapping(List items, Envelope searchEnv) {
    List result = List.empty(growable: true);
    for (int i = 0; i < items.length; i++) {
      Envelope env = items[i];
      if (env.intersectsEnvelope(searchEnv)) result.add(env);
    }
    return result;
  }

  void queryTest(Envelope env) {
    List<dynamic> candidateList = q.query(env);
    List finalList = getOverlapping(candidateList, env);

    List eList = envList.query(env);
//System.out.println(finalList.size());

    if (finalList.length != eList.length)
      throw new RuntimeException("queries do not match");
  }

  void queryGrid(int nGridCells, double cellSize) {
    Stopwatch sw = Stopwatch();
    sw.start();

    int gridSize = math.sqrt(nGridCells.toDouble()).toInt();
    gridSize += 1;
    double extent = MAX_EXTENT - MIN_EXTENT;
    double gridInc = extent / gridSize;

    for (int i = 0; i < gridSize; i++) {
      for (int j = 0; j < gridSize; j++) {
        double x = MIN_EXTENT + gridInc * i;
        double y = MIN_EXTENT + gridInc * j;
        Envelope env = new Envelope(x, x + cellSize, y, y + cellSize);
        queryTest(env);
        //queryTime(env);
      }
    }
    print('Time = ${sw.elapsed.inMicroseconds}');
  }

  void runQueries() {
    int nGridCells = 100;
    int cellSize = math.sqrt(NUM_ITEMS.toDouble()).toInt();
    double extent = MAX_EXTENT - MIN_EXTENT;
    double queryCellSize = 2.0 * extent / cellSize;

    queryGrid(nGridCells, queryCellSize);

    //queryGrid(200);
  }

  test('quadtreeCorrectTest', () {
    fill();
    print('depth = ${q.depth()} size = ${q.size()}');
    runQueries();
  });
}

class EnvelopeList {
  List<Envelope> envList = List.empty(growable: true);

  EnvelopeList();

  void add(Envelope env) {
    envList.add(env);
  }

  List<Envelope> query(Envelope searchEnv) {
    List<Envelope> result = List.empty(growable: true);
    for (Envelope env in envList) {
      if (env.intersectsEnvelope(searchEnv)) result.add(env);
    }
    return result;
  }
}
