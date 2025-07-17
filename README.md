
# MFCC for Dart 🎵

A pure Dart, null-safe library for extracting **Mel-Frequency Cepstral Coefficients (MFCCs)** from an audio signal. MFCCs are the most widely used features for speech recognition and audio analysis.

This package provides a flexible and easy-to-use API for both one-off calculations and real-time stream processing, with no platform-specific dependencies.

## Features ✨

* **Static Method**: Quickly extract MFCCs from a complete signal.
* **Instance-based Processing**: Process audio frame-by-frame with a stateful processor.
* **Stream Processing**: Handle real-time audio streams seamlessly.
* **Highly Customizable**: Control all key parameters like sample rate, FFT size, number of filters, and coefficients.
* **Standard Options**: Includes pre-emphasis and log-energy replacement for the first coefficient.
* **Pure Dart**: Works on any platform where Dart runs (Mobile, Desktop, Web).
* **Null-safe and well-documented**.

---

## Installation 💻

Add this to your project's `pubspec.yaml` file:

```yaml
dependencies:
  mfcc: ^latest # Replace with the actual latest version
```

Then, run `dart pub get` or `flutter pub get`.

---

## How to Use 🚀

There are three main ways to use this library, depending on your needs.

### 1\. Static Method (Easiest)

This is the simplest way to get MFCCs from a complete audio signal.

```dart
import 'package:mfcc/mfcc.dart';
import 'dart:math';

void main() {
  // 1. Generate a sample audio signal (e.g., a 440 Hz sine wave)
  final int sampleRate = 16000;
  final signalDuration = 1; // seconds
  final signalLength = sampleRate * signalDuration;
  final signal = List<double>.generate(
      signalLength, (i) => sin(2 * pi * 440 * i / sampleRate));

  // 2. Define MFCC parameters
  final int windowLength = 400; // ~25 ms
  final int windowStride = 160; // ~10 ms
  final int fftSize = 512;
  final int numFilters = 40;
  final int numCoefs = 13;

  // 3. Extract features
  List<List<double>> features = MFCC.mfccFeats(
    signal,
    sampleRate,
    windowLength,
    windowStride,
    fftSize,
    numFilters,
    numCoefs,
    energy: true,
    preEmphasis: 0.97,
  );

  // 4. Print the results (num_frames x num_coeffs)
  print('Extracted ${features.length} frames.');
  print('Each frame has ${features[0].length} coefficients.');
  // Example: print the first frame's features
  print('First frame features: ${features[0]}');
}
```

### 2\. Instance-based Processing

This approach is ideal when you need to process audio frame by frame.

```dart
import 'package:mfcc/mfcc.dart';

void main() {
  // 1. Instantiate and configure the MFCC processor
  final processor = MFCC(
    sampleRate: 16000,
    fftSize: 512,
    numFilters: 40,
    numCoefs: 13,
  );

  // 2. Assume you have pre-framed audio data
  // (e.g., from a microphone buffer)
  final List<List<double>> audioFrames = [
    List<double>.filled(400, 0.5), // Frame 1
    List<double>.filled(400, 0.3), // Frame 2
    List<double>.filled(400, 0.8), // Frame 3
  ];

  // 3. Process all frames at once
  List<List<double>> features = processor.processFrames(audioFrames);
  print('Processed ${features.length} frames.');
  print('Features for frame 1: ${features[0]}');

  // Or process frames one by one
  List<double> singleFrameFeatures = processor.processFrame(audioFrames[0]);
  print('Features for frame 1 (processed individually): $singleFrameFeatures');
}
```

### 3\. Stream Processing

Use streams for real-time applications, like processing audio from a microphone.

```dart
import 'package:mfcc/mfcc.dart';
import 'dart:async';

void main() async {
  // 1. Instantiate the MFCC processor
  final processor = MFCC(
    sampleRate: 16000,
    fftSize: 512,
    numFilters: 40,
    numCoefs: 13,
  );

  // 2. Create a mock stream of audio frames (in a real app, this would
  // come from a source like `mic_stream`)
  final audioFrameStream = Stream.periodic(
    const Duration(milliseconds: 100),
    (i) => List<double>.filled(400, i * 0.1),
  ).take(5);

  // 3. Set the stream and get the output feature stream
  final featureStream = processor.setStream(audioFrameStream);

  // 4. Listen for new MFCC features as they are computed
  if (featureStream != null) {
    print('Listening for MFCC features...');
    await for (final features in featureStream) {
      print('New features computed: [${features.map((e) => e.toStringAsFixed(2)).join(', ')}]');
    }
  }
  
  // 5. Clean up when done
  processor.cancelStream();
  print('Stream processing finished.');
}
```

---

## API Overview 📖

* `MFCC` **class**: The main class for processing.
  * **Constructor**: `MFCC({ required int sampleRate, ... })` configures the processor.
  * `processFrame(List<double> frame)`: Processes a single frame.
  * `processFrames(List<List<double>> frames)`: Processes a list of frames.
  * `setStream(Stream<List<double>> audioInput)`: Sets up real-time stream processing.
  * `cancelStream()`: Stops and cleans up stream resources.
* `MFCC.mfccFeats(...)`: Static helper for a one-off calculation on a full signal.
* `ValueError`: Custom exception for invalid parameters.

---

## License

This project is licensed under the **GNU Affero General Public License v3.0**. See the [LICENSE](https://www.gnu.org/licenses/agpl-3.0.en.html) file for details.

---

## Contributing

Contributions are welcome\! If you find a bug or want to suggest a new feature, please open an issue.
