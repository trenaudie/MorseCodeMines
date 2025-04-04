#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <numeric>

#include "fft_copied/kiss_fft.h"
#include "fft_copied/kiss_fftr.h" // Include if you need real FFTs


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
    return 0;
}