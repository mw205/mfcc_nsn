import 'package:mcfcc_nsn/mcfcc_nsn.dart';
import 'package:test/test.dart';

void main() {
  group('A group of tests', () {
    List<List<double>> features = [];
    setUp(() {
      var sampleRate = 16000;
      var windowLength = 1024;
      var windowStride = 512;
      var fftSize = 512;
      var numFilter = 20;
      var numCoefs = 13;
      var mySignal = List<double>.generate(16000, (i) => i.toDouble());
      features = MFCC.mfccFeats(
        mySignal,
        sampleRate,
        windowLength,
        windowStride,
        fftSize,
        numFilter,
        numCoefs,
      );
    });

    test('MFCC test', () {
      expect(features.length, 30);
      expect(features[0].length, 13);
    });
  });
}
