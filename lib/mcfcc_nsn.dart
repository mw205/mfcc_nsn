/// Support for doing something awesome.
///
/// More dartdocs go here.
library;

import 'dart:async';
import 'dart:math';

import 'package:fft_nullsafety/fft.dart';

/// Converts hertz [freq] to mel
double hertzToMel(double freq) {
  return 1127 * log(1 + freq / 700);
}

/// Converts [mel] value to hertz
double melToHertz(double mel) {
  return 700 * (exp(mel / 1127) - 1);
}

/// Returns the log of a [value]... safely.
///
/// If [value] <=0 return the smallest possible value
double safeLog(double value) {
  if (value <= 0) {
    return log(double.minPositive);
  } else {
    return log(value);
  }
}

/// Class to extract MFCC features from a signal.
///
/// There are 3 ways to use this class:
/// - Use the static method mfccFeats() to extract features from a signal.
/// - Intantiate MFCC to process frames on the go with processFrame() or processFrames().
/// - Use the setStream method to process MFCC with Streams.
///
/// MFCC are generated on each windows by:
/// 1. (Optional) Aplying pre-Emphasis
/// 1. Computing Spectrum
/// 2. Applying triangular filter in the MEL domain.
/// 3. Computing log for each band
/// 4. Applying Discrete Cosinus Transform
/// 5. (Optional) Replace first value by the window log-Energy.
class MFCC {
  /// FFT size
  late int _fftSize;

  /// The number of mel filters
  late int _numFilters;

  /// Number of output values
  late int _numCoefs;

  /// The mel filters triangle windows indexes, generated once
  late List<List<double>> _fbanks;

  /// If set to true, replace mfcc first value with spectrum logenergy
  late bool _energy;

  /// PreEmphasis parameters
  late bool _useEmphasis;

  /// Emphasis factor
  late double? _emphasis;

  /// Previous frame last value for continue preEmphasis
  num _lastValue = 0.0;

  /// Stream input
  StreamSubscription<List<double>>? _audioInput;

  /// Stream output
  StreamController<List<double>>? _featureStream;

  /// Creates an MFCC processor.
  ///
  /// * [sampleRate] - The signal sampling rate in samples/s.
  /// * [fftSize] - The number of FFT points.
  /// * [numFilters] - The number of MEL filters to apply.
  /// * [numCoefs] - The number of cepstral coefficients to keep.
  /// * {[energy] = true} - If true, the first coefficient (C0) is replaced by the log-energy of the frame.
  /// * {[preEmphasis] = 0.97} - The factor for pre-emphasis. If null, no pre-emphasis is applied.
  MFCC({
    required int sampleRate,
    required int fftSize,
    required int numFilters,
    required int numCoefs,
    bool energy = true,
    double? preEmphasis = 0.97,
  }) {
    _fftSize = fftSize;
    _numFilters = numFilters;
    _numCoefs = numCoefs;

    if (!energy) {
      // When not using energy, the 0th DCT coefficient (DC offset) is discarded.
      // We compute one extra coefficient to compensate after sublisting.
      _numCoefs += 1;
    }

    _energy = energy;

    _useEmphasis = (preEmphasis != null);
    _emphasis = preEmphasis;

    _fbanks = MFCC.filterbanks(
      sampleRate,
      _numFilters,
      ((_fftSize / 2) + 1).toInt(),
    );
  }

  /// Apply preEmphasis filter on given signal.
  static List<double> preEmphasis(
    List<num> signal,
    double emphasisFactor, {
    num lastValue = 0.0,
  }) {
    var empSignal = List<double>.filled(signal.length, 0);
    var swSig = [lastValue] + signal.sublist(0, signal.length - 1);
    for (var i = 0; i < signal.length; i++) {
      empSignal[i] = signal[i] - swSig[i] * emphasisFactor;
    }
    return empSignal;
  }

  /// Returns the mel filters
  ///
  /// 1. Linearly splits the frequency interval using the mel scale.
  /// 2. Generates triangular overlapping windows
  /// 3. Generates filter coeficient for each window
  /// 4. Returns [numFilt] filters of length [nFFT]
  static List<List<double>> filterbanks(int samplerate, int numFilt, int nFFT) {
    var interval = hertzToMel(samplerate.toDouble()) / (numFilt + 1);
    var gridMels = List<double>.generate(numFilt + 2, (v) => v * interval);
    var gridHertz = gridMels.map((v) => melToHertz(v)).toList();
    var gridIndexes = gridHertz
        .map((v) => (v * nFFT / samplerate).floor())
        .toList();

    var filters = List<List<double>>.filled(numFilt, []);
    for (var i = 0; i < numFilt; i++) {
      var left = List<double>.generate(
        gridIndexes[i + 1] - gridIndexes[i],
        (v) => v / (gridIndexes[i + 1] - gridIndexes[i]),
      );
      var right = List<double>.generate(
        gridIndexes[i + 2] - gridIndexes[i + 1],
        (v) => v / (gridIndexes[i + 2] - gridIndexes[i + 1]),
      ).reversed.toList();
      var filter = [
        List<double>.filled(gridIndexes[i], 0.0),
        left,
        [1.0],
        right,
        List<double>.filled(nFFT - gridIndexes[i + 2], 0.0),
      ].expand((x) => x).toList();
      filters[i] = filter;
    }
    return filters;
  }

  /// Returns the power spectrum of a given [frame].
  /// The power spectrum is computed as `(real^2 + imag^2) / fftSize`.
  static List<double> powerSpectrum(List<double> frame, fftSize) {
    var fft = FFT().Transform(frame.sublist(0, fftSize));
    return fft
        .sublist(0, (fftSize / 2 + 1).round())
        .map((v) => (pow(v.real, 2) + pow(v.imaginary, 2)) / fftSize)
        .toList();
  }

  /// Maps the power spectrum over the mel [filters] to obtains a condensed spectrogram on the mel scale.
  static List<double> melCoefs(
    List<double> powerSpec,
    List<List<double>> filters,
  ) {
    var nFilt = filters.length;
    var result = List<double>.filled(nFilt, 0);
    for (var i = 0; i < nFilt; i++) {
      double sum = 0;
      for (var j = 0; j < powerSpec.length; j++) {
        sum += powerSpec[j] * filters[i][j];
      }
      result[i] = sum;
    }
    return result.map((v) => safeLog(v)).toList();
  }

  /// Returns the discrete cosinus transform.
  ///
  /// Uses the [scipy](https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.fftpack.dct.html) type-II dct implementation:
  ///
  /// `y(k) = 2 * sum{n ∈ [0,N-1]} x(n) * cos(pi * k * (2n+1)/(2 * N)), 0 <= k < N.`
  ///
  /// if norm is set to true apply a scaling factor f to y(k) as followed:
  ///
  /// - `f = sqrt(1/(4*N)) if k = 0`
  /// - `f = sqrt(1/(2*N)) otherwise.`
  ///
  static List<double> dct(List<double> x, bool norm) {
    var result = List<double>.filled(x.length, 0);
    var N = x.length;
    var sum = 0.0;
    var scalingFactor0 = sqrt(1 / (4 * N));
    var scalingFactor = sqrt(1 / (2 * N));
    for (var k = 0; k < N; k++) {
      sum = 0.0;
      for (var n = 0; n < N; n++) {
        sum += x[n] * cos(pi * k * (2 * n + 1) / (2 * N));
      }
      sum *= 2;
      if (norm) {
        if (k == 0) {
          sum = sum * scalingFactor0;
        } else {
          sum = sum * scalingFactor;
        }
      }
      result[k] = sum;
    }
    return result;
  }

  /// Returns the MFCC values for the given frame
  List<double> processFrame(List<double> frame) {
    if (_useEmphasis) {
      var v = frame.last;
      frame = preEmphasis(frame, _emphasis!, lastValue: _lastValue);
      _lastValue = v;
    }
    var powerSpectrum = MFCC.powerSpectrum(frame, _fftSize);
    var melCoefs = MFCC.melCoefs(powerSpectrum, _fbanks);
    var mfccs = MFCC.dct(melCoefs, true).sublist(0, _numCoefs);
    if (_energy) {
      mfccs[0] = safeLog(
        powerSpectrum.reduce((a, b) => a + b),
      ); // Replace first value with logenergy
    } else {
      return mfccs.sublist(1);
    }
    return mfccs;
  }

  /// Returns the MFCC values of a list of frames
  List<List<double>> processFrames(List<List<double>> frames) {
    List<List<double>> mfccs = [];
    for (List<double> frame in frames) {
      mfccs.add(processFrame(frame));
    }
    return mfccs;
  }

  /// Generates MFCC features from a [signal].
  ///
  /// * [signal] - The signal to extract the features from.
  /// * [sampleRate] - The signal sampling rate -sample/s- (> 0).
  ///
  /// MFCC are extracted using temporal sliding windows
  /// * [windowLength] - Window length in number of samples (> 0 && <= [signal].length).
  /// * [windowStride] - Window stride in number of sample (> 0)
  ///
  ///
  /// * [fftSize] - Number of fft generated by window (> 0)
  /// * [numFilters] - Number of MEL filters (> 0)
  /// * [numCoefs] - Number of cepstral coefficient to keep (> 0 && <= [numFilters])
  /// * {[energy] = true} - If True, replaces the first value by the window log-energy.
  /// * {[preEmphasis] = 0.97} - Apply signal preEmphasis. If the value is null, does nothing.
  ///
  ///
  /// Throws [ValueError] if any value is off limit.
  static List<List<double>> mfccFeats(
    List<double> signal,
    int sampleRate,
    int windowLength,
    int windowStride,
    int fftSize,
    int numFilters,
    int numCoefs, {
    bool energy = true,
    double preEmphasis = 0.97,
  }) {
    if (sampleRate <= 0) {
      throw ValueError('Sample rate must be > 0 (Got $sampleRate).');
    }
    if (windowLength <= 0) {
      throw ValueError('Window length must be > 0 (Got $windowLength).');
    }
    if (windowStride <= 0) {
      throw ValueError('Stride must be > 0 (Got $windowStride).');
    }
    if (windowLength > signal.length) {
      throw ValueError(
        'Window length cannot be greater than signal length.(Got $windowLength > ${signal.length})',
      );
    }
    if (numFilters <= 0) {
      throw ValueError('Number of filter must be positive (Got $numFilters).');
    }
    if (numCoefs <= 0) {
      throw ValueError(
        'Number of coefficient must be positive (Got $numCoefs).',
      );
    }
    if (numCoefs > numFilters) {
      throw ValueError(
        'Number of coefficient must be inferior to the number of filters (Got $numCoefs > $numFilters).',
      );
    }

    var frames = splitSignal(signal, windowLength, windowStride);
    var processor = MFCC(
      sampleRate: sampleRate,
      fftSize: fftSize,
      numFilters: numFilters,
      numCoefs: numCoefs,
      energy: energy,
      preEmphasis: preEmphasis,
    );
    return processor.processFrames(frames);
  }

  /// Splits a [signal] into overlapping frames.
  ///
  /// * [signal] - The input signal.
  /// * [windowLength] - The length of each frame in samples.
  /// * [windowStride] - The step size between the start of consecutive frames.
  static List<List<double>> splitSignal(
    List<double> signal,
    int windowLength,
    int windowStride,
  ) {
    var nFrames = ((signal.length - windowLength) / windowStride).floor() + 1;
    var frames = List<List<double>>.filled(nFrames, [], growable: false);
    for (var i = 0; i < nFrames; i++) {
      frames[i] = signal.sublist(
        i * windowStride,
        (i * windowStride) + windowLength,
      );
    }
    return frames;
  }

  /// Set input as `Stream<List<double>>`.
  /// Returns a `StreamController<List<double>>` on which features will be pushed.
  /// [audioInput] must provide frame of desired length.
  StreamController<List<double>>? setStream(Stream<List<double>> audioInput) {
    cancelStream();
    _featureStream = StreamController<List<double>>.broadcast();
    _audioInput = audioInput.listen((frame) {
      _featureStream?.add(processFrame(frame));
    });
    return _featureStream;
  }

  /// Cancel streamSubscription and closes egress feature stream.
  void cancelStream() {
    _audioInput?.cancel();
    _featureStream?.close();
    _audioInput = null;
    _featureStream = null;
  }
}

class ValueError implements Exception {
  String errMsg = '';
  ValueError(String msg) {
    errMsg = msg;
  }
}
