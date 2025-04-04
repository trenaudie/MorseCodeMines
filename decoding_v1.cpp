#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <numeric>


// Windowed normalization function (translated from Python)
std::vector<float> windowed_normalization(const std::vector<float> &values,
                                          int window_size = 1024,
                                          float target_level = 0.5f,
                                          const std::string &method = "rms")
{
    // Create Hanning window
    std::vector<float> window(window_size);
    for (int i = 0; i < window_size; i++)
    {
        window[i] = 0.5f * (1.0f - std::cos(2.0f * M_PI * i / (window_size - 1)));
    }

    // Pad the input signal to handle edge cases
    int pad_size = window_size / 2;
    std::vector<float> padded_values(values.size() + 2 * pad_size);

    // Edge padding
    for (int i = 0; i < pad_size; i++)
    {
        padded_values[i] = values[0];
        padded_values[padded_values.size() - 1 - i] = values[values.size() - 1];
    }

    // Copy original values
    for (size_t i = 0; i < values.size(); i++)
    {
        padded_values[i + pad_size] = values[i];
    }

    // Initialize output arrays
    std::vector<float> output_accumulator(padded_values.size(), 0.0f);
    std::vector<float> window_accumulator(padded_values.size(), 0.0f);

    // Step size for 50% overlap
    int hop_size = window_size / 2;

    // Process the signal in overlapping windows
    for (size_t i = 0; i < values.size(); i += hop_size)
    {
        // Define the range for this window
        size_t start = i;
        size_t end = start + window_size;

        if (end > padded_values.size())
        {
            break;
        }

        // Extract the current frame
        std::vector<float> frame(window_size);
        for (size_t j = 0; j < window_size; j++)
        {
            frame[j] = padded_values[start + j];
        }

        // Apply window to the frame
        std::vector<float> windowed_frame(window_size);
        for (int j = 0; j < window_size; j++)
        {
            windowed_frame[j] = frame[j] * window[j];
        }

        // Compute normalization factor based on method
        float level = 0.0f;
        if (method == "rms")
        {
            // Root Mean Square level
            float sum_squared = 0.0f;
            for (float sample : windowed_frame)
            {
                sum_squared += sample * sample;
            }
            level = std::sqrt(sum_squared / window_size);
        }
        else
        { // peak
            // Peak level
            level = 0.0f;
            for (float sample : windowed_frame)
            {
                level = std::max(level, std::abs(sample));
            }
        }

        // Avoid division by zero
        float gain = (level < 1e-10f) ? 1.0f : target_level / level;

        // Apply gain and window
        for (size_t j = 0; j < window_size; j++)
        {
            float normalized_sample = frame[j] * gain;
            float windowed_normalized = normalized_sample * window[j];
            output_accumulator[start + j] += windowed_normalized;
            window_accumulator[start + j] += window[j];
        }
    }

    // Create final normalized output
    std::vector<float> normalized_output(values.size());
    for (size_t i = 0; i < values.size(); i++)
    {
        // Divide by window weights with small epsilon to avoid division by zero
        normalized_output[i] = output_accumulator[i + pad_size] /
                               (window_accumulator[i + pad_size] + 1e-10f);
    }

    return normalized_output;
}
std::vector<float> open_wav(const std::string &filename)
{
    std::cout << "Attempting to open WAV file: " << filename << std::endl;
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open " << filename << std::endl;
        return {};
    }

    std::string header(44, '\0');
    // istream& read (char* s, streamsize n); so must use the pointer to the first element of the string
    file.read(&header[0], 44);
    std::cout << "Header size: " << header.size() << std::endl;
    // Extract some useful information from the header

    // print header22
    std::cout << "Header[22]: " << static_cast<int>(header[22]) << std::endl;
    std::cout << "Header[23]: " << static_cast<int>(header[23]) << std::endl;

    // manually calculate the number of channels, using hexadecimal
    unsigned char byte1 = static_cast<unsigned char>(header[22]);
    unsigned char byte2 = static_cast<unsigned char>(header[23]);
    short channels = byte1 + (byte2 << 8);
    std::cout << "Channels: " << channels << std::endl;

    // automatically retrieving sample rate, bits per sample and file size, using pointer arithmetic
    int sample_rate = *(int *)&header[24];
    int bits_per_sample = *(short *)&header[34];
    int file_size = *(int *)&header[40];

    std::vector<float> values;
    // buffer to read samples, for 16bit audio
    char buffer[2];

    while (file.read(buffer, 2))
    {
        // convert to signed short (16-bit)
        short value = *(short *)buffer;
        float normalized = value / 32768.0f;

        values.push_back(normalized);
    }

    file.close();
    std::cout << "Number of samples: " << values.size() << std::endl;
    std::cout << "Sample rate: " << sample_rate << std::endl;
    std::cout << "Bits per sample: " << bits_per_sample << std::endl;
    return values;
}


#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>

/**
 * Compute the Discrete Fourier Transform (DFT) of a signal.
 * 
 * @param arr Input signal array
 * @return Vector of complex DFT values
 */
std::vector<std::complex<double>> dft(const std::vector<double>& arr) {
    const size_t N = arr.size();
    std::vector<std::complex<double>> dft_values(N, std::complex<double>(0, 0));
    
    const double TWO_PI = 2.0 * M_PI;
    
    // Display progress
    std::cout << "Computing DFT..." << std::endl;
    
    for (size_t k = 0; k < N; ++k) {
        // Show progress every 10%
        if (k % (N / 10) == 0) {
            std::cout << (k * 100) / N << "% complete" << std::endl;
        }
        
        for (size_t n = 0; n < N; ++n) {
            double angle = TWO_PI * k * n / N;
            dft_values[k] += arr[n] * std::complex<double>(cos(angle), -sin(angle));
        }
    }
    
    std::cout << "100% complete" << std::endl;
    return dft_values;
}

/**
 * Find non-zero segments in the signal that meet minimum length requirements.
 * 
 * @param values Input signal array
 * @param min_length Minimum required segment length
 * @return Pair of indices (start, end) of the first valid segment
 */
std::pair<size_t, size_t> findNonZeroSegment(const std::vector<double>& values, size_t min_length = 1000) {
    size_t non_null_index = 0;
    size_t range_end = 0;
    
    while (non_null_index < values.size() - 1) {
        std::cout << "Checking from index: " << non_null_index << std::endl;
        
        // Skip zeros
        while (non_null_index < values.size() && values[non_null_index] == 0) {
            non_null_index++;
        }
        
        // Find end of non-zero segment
        range_end = non_null_index + 1;
        while (range_end < values.size() - 1 && 
               (values[range_end] != 0 || values[range_end + 1] != 0)) {
            range_end++;
            
            if ((range_end - non_null_index) % 100 == 0) {
                std::cout << "Range length: " << range_end - non_null_index << std::endl;
            }
        }
        
        // Check if segment meets minimum length requirement
        if (range_end - non_null_index > min_length) {
            std::cout << "Found valid segment from " << non_null_index 
                      << " to " << range_end 
                      << " (length: " << range_end - non_null_index << ")" << std::endl;
            return {non_null_index, range_end};
        }
        
        // Move to next potential segment
        if (range_end < values.size() - 1 && values[range_end] == 0 && values[range_end + 1] == 0) {
            std::cout << "Found empty signal at index " << non_null_index << std::endl;
            non_null_index = range_end + 1;
        } else {
            // Safety check to avoid infinite loops
            non_null_index++;
        }
    }
    
    // Return invalid segment if nothing found
    return {0, 0};
}

/**
 * Calculate FFT frequencies for the given number of samples and sample rate.
 * 
 * @param n Number of samples
 * @param sample_rate Sampling rate in Hz
 * @return Vector of frequency values
 */
std::vector<double> calculateFrequencies(size_t n, double sample_rate) {
    std::vector<double> frequencies(n);
    
    for (size_t i = 0; i < n; ++i) {
        if (i <= n / 2) {
            frequencies[i] = i * sample_rate / n;
        } else {
            frequencies[i] = (i - n) * sample_rate / n;
        }
    }
    
    return frequencies;
}

/**
 * Find dominant frequencies in the DFT result.
 * 
 * @param dft_values DFT values
 * @param frequencies Frequency values
 * @param top_count Number of top frequencies to return
 * @return Vector of top frequencies
 */
std::vector<double> findDominantFrequencies(
    const std::vector<std::complex<double>>& dft_values,
    const std::vector<double>& frequencies,
    size_t top_count = 10
) {
    // Calculate magnitude for each DFT value
    std::vector<std::pair<double, size_t>> magnitudes;
    for (size_t i = 0; i < dft_values.size(); ++i) {
        double magnitude = std::abs(dft_values[i]);
        magnitudes.push_back({magnitude, i});
    }
    
    // Sort by magnitude (descending)
    std::sort(magnitudes.begin(), magnitudes.end(), 
        [](const auto& a, const auto& b) { return a.first > b.first; });
    
    // Extract top frequencies
    std::vector<double> top_frequencies;
    for (size_t i = 0; i < std::min(top_count, magnitudes.size()); ++i) {
        top_frequencies.push_back(frequencies[magnitudes[i].second]);
    }
    
    return top_frequencies;
}

// Example usage
void analyzeSignal(const std::vector<double>& values, double sample_rate) {
    // Find a suitable segment
    auto [start, end] = findNonZeroSegment(values);
    
    if (end - start < 1000) {
        std::cerr << "No suitable segment found!" << std::endl;
        return;
    }
    
    // Extract segment (limit to 3000 samples if larger)
    size_t length = std::min(static_cast<size_t>(3000), end - start - 2);
    std::vector<double> values_sinusoidal(values.begin() + start + 1, 
                                         values.begin() + start + 1 + length);
    
    // Compute DFT
    auto dft_values = dft(values_sinusoidal);
    
    // Calculate frequencies
    auto frequencies = calculateFrequencies(dft_values.size(), sample_rate);
    
    // Find dominant frequencies
    auto top_frequencies = findDominantFrequencies(dft_values, frequencies);
    
    // Print results
    std::cout << "Top frequencies:" << std::endl;
    for (auto freq : top_frequencies) {
        std::cout << freq << " Hz" << std::endl;
    }
    
    std::cout << "Major frequency: " << top_frequencies[0] << " Hz" << std::endl;
}

int main()
{
    std::string filename = "morse_code.wav";
    std::vector<float> audio_data = open_wav(filename);
    if (audio_data.empty())
    {
        std::cerr << "Error: No audio data read." << std::endl;
        return 1;
    }
    std::cout << "WAV file opened successfully." << std::endl;

    float min_val = *std::min_element(audio_data.begin(), audio_data.end());
    float max_val = *std::max_element(audio_data.begin(), audio_data.end());

    float sum = std::accumulate(audio_data.begin(), audio_data.end(), 0.0f,
                                [](float acc, float val)
                                { return acc + std::abs(val); });
    float avg_amplitude = sum / audio_data.size();

    std::cout << "Non-normalized audio stats:" << std::endl;
    std::cout << "  Min value: " << min_val << std::endl;
    std::cout << "  Max value: " << max_val << std::endl;
    std::cout << "  Average amplitude: " << avg_amplitude << std::endl;

    // viewing first 100 samples 
    std::cout << "First 100 samples:" << std::endl;
    for (int i = 0; i < 100 && i < audio_data.size(); i++)
    {
        std::cout << audio_data[i] << " ";
    }
    std::cout << std::endl;

    // Convert float audio data to double for DFT processing
    std::vector<double> audio_double(audio_data.begin(), audio_data.end());
    
    // For better performance, limit analysis to a reasonable segment
    // If audio is long, analyze a segment of 10000 samples
    const size_t MAX_ANALYSIS_SAMPLES = 10000;
    std::vector<double> analysis_segment;
    
    if (audio_double.size() > MAX_ANALYSIS_SAMPLES) {
        // Find a non-zero segment to analyze
        auto [start, end] = findNonZeroSegment(audio_double, 1000);
        
        if (end - start > 1000) {
            // Use found segment
            size_t length = std::min(MAX_ANALYSIS_SAMPLES, end - start - 2);
            analysis_segment = std::vector<double>(audio_double.begin() + start + 1, 
                                                 audio_double.begin() + start + 1 + length);
        } else {
            // If no suitable segment found, use first MAX_ANALYSIS_SAMPLES
            analysis_segment = std::vector<double>(audio_double.begin(), 
                                                 audio_double.begin() + MAX_ANALYSIS_SAMPLES);
        }
    } else {
        // If audio is short enough, use all of it
        analysis_segment = audio_double;
    }
    
    std::cout << "Analyzing segment of " << analysis_segment.size() << " samples" << std::endl;
    
    // Compute DFT
    std::vector<std::complex<double>> dft_values = dft(analysis_segment);
    
    // Get sampling rate (assuming 44100 Hz, adjust as needed)
    const double SAMPLE_RATE = 44100.0;
    
    // Calculate frequencies corresponding to DFT values
    std::vector<double> frequencies = calculateFrequencies(dft_values.size(), SAMPLE_RATE);
    
    // Find top frequencies
    std::vector<double> top_frequencies = findDominantFrequencies(dft_values, frequencies, 10);
    
    // Print top frequencies
    std::cout << "\nTop 10 frequencies detected:" << std::endl;
    for (size_t i = 0; i < top_frequencies.size(); ++i) {
        std::cout << i+1 << ". " << top_frequencies[i] << " Hz" << std::endl;
    }
    
    // The top frequency is the first one in the list
    double top_frequency = top_frequencies[0];
    std::cout << "\nMajor frequency: " << top_frequency << " Hz" << std::endl;
    
    // For Morse code analysis, this frequency might be the carrier tone
    std::cout << "This frequency is likely the carrier tone for the Morse code." << std::endl;
    
    


    return 0;
}