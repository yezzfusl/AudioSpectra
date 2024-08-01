# AudioSpectra
**Advanced real-time audio effects and processing system.**

## Overview
AudioSpectra is a high-performance audio processing system designed for real-time effects, mixing, and analysis. It offers a robust framework for implementing a variety of audio effects and provides a user-friendly interface for parameter adjustments. The system is built to ensure low-latency processing and high performance.
## Features

### Audio Processing
- **Real-Time Effects:** Includes real-time equalization, reverb, and compression effects.
- **Frequency-Domain Processing:** Utilizes FFTW for efficient frequency-domain transformations.
- **Dynamic Effects:** Supports time-stretching and pitch-shifting using SoundTouch.

### User Interface
- **Parameter Adjustment:** Allows real-time adjustment of audio parameters through an intuitive interface.
- **Effect Control:** Provides control over equalization, reverb, and compression settings.
### Data Storage
- [ ] **Persistent Storage:** Stores audio samples and presets in a SQLite database.
- [ ] **Serialization:** Utilizes Boost.Serialization for managing data.

## Technical Details

### Libraries and Frameworks

- [x] **PortAudio:** Cross-platform audio I/O library for handling audio input/output.
- [x] **JUCE:** Framework for audio processing and UI development.
- [x]  **FFTW:** Library for computing discrete Fourier transforms.
- [x] **SoundTouch:** Provides time-stretching and pitch-shifting capabilities.
- [ ] **SQLite:** Lightweight database engine for storing samples and presets.
- [ ] **Boost.Serialization:** For serialization and deserialization of data.
- [ ] **Intel IPP / ARM CMSIS-DSP:** Used for optimizing DSP operations to enhance performance and efficiency.

## Installation

### Prerequisites

- PortAudio
- JUCE
- FFTW
- SoundTouch
- SQLite
- Boost.Serialization
- Intel IPP / ARM CMSIS-DSP

### Build Instructions
1. **Clone the Repository:**
          ```bash
   git clone https://github.com/yezzfusl/AudioSpectra.git 
   cd AudioSpectra

2. Install Dependencies:
Follow the installation instructions for each library and framework.
    ```bash
    mkdir build
    cd build
    cmake ..
    make

## Usage:
1. Run the Application:
Execute the built application from the build directory.

1. Adjust Parameters:
Use the UI to adjust audio parameters and effects in real-time.


