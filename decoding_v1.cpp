#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <set>
#include <map>
#include "utils.h"


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
std::pair<std::vector<float>, float> open_wav(const std::string &filename)
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
    return {values, static_cast<float>(sample_rate)};
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
    return dft_values ;
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

float compute_top_frequency(const std::vector<float>& values, int sample_rate)
{

    // Convert float audio data to double for DFT processing
    std::vector<double> audio_double(values.begin(), values.end());
    
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
    // Filter out frequencies above 10,000 Hz
    std::vector<double> filtered_frequencies;
    for (const double& freq : top_frequencies) {
        if (freq <= 10000.0 && freq >= 0.01) {
            filtered_frequencies.push_back(freq);
        }
        std::cout << "Detected frequency: " << freq << " Hz" << std::endl;
    }
    // Sort frequencies (assuming you want them in descending order)
    std::sort(filtered_frequencies.begin(), filtered_frequencies.end(), std::greater<double>());

    
    // The top frequency is the first one in the list
    double top_frequency = top_frequencies[0];
    std::cout << "\nMajor frequency: " << top_frequency << " Hz" << std::endl;
    
    return top_frequency;
    
}



void remove_end_blanks(std::vector<float>& values)
{
    // Remove trailing zeros
    while (!values.empty() && values.back() == 0)
    {
        values.pop_back();
    }
}


// Define Morse code dictionary
std::map<std::string, char> morse_dict = {
    {".-", 'A'}, {"-...", 'B'}, {"-.-.", 'C'}, {"-..", 'D'}, {".", 'E'},
    {"..-.", 'F'}, {"--.", 'G'}, {"....", 'H'}, {"..", 'I'}, {".---", 'J'},
    {"-.-", 'K'}, {".-..", 'L'}, {"--", 'M'}, {"-.", 'N'}, {"---", 'O'},
    {".--.", 'P'}, {"--.-", 'Q'}, {".-.", 'R'}, {"...", 'S'}, {"-", 'T'},
    {"..-", 'U'}, {"...-", 'V'}, {".--", 'W'}, {"-..-", 'X'}, {"-.--", 'Y'},
    {"--..", 'Z'}, {"-----", '0'}, {".----", '1'}, {"..---", '2'}, {"...--", '3'},
    {"....-", '4'}, {".....", '5'}, {"-....", '6'}, {"--...", '7'}, {"---..", '8'},
    {"----.", '9'}, {".-.-.-", '.'}, {"--..--", ','}, {"..--..", '?'}, {".----.", '\''},
    {"-.-.--", '!'}, {"-..-.", '/'}, {"-.--.", '('}, {"-.--.-", ')'}, {".-...", '&'},
    {"---...", ':'}, {"-.-.-.", ';'}, {"-...-", '='}, {".-.-.", '+'}, {"-....-", '-'},
    {"..--.-", '_'}, {".-..-.", '"'}, {"...-..-", '$'}, {".--.-.", '@'}
};

std::string decoding(std::string symbols) {
    // example
    // std::string symbols = ".... . .-.. .-.. --- / .-- --- .-. .-.. -.."; //HELLO WORLD
    int ptr = 0; // position in the symbols string
    std::string decoded_text = "";
    std::string buffer = "";
    
    // Fixed issues:
    // 1. Cannot use find() on a vector with a char
    // 2. Need to increment ptr
    // 3. Need to handle end of string
    
    std::cout << "Morse code to decode: " << symbols << std::endl;
    
    while (ptr < symbols.size()) {
        char ptr_char = symbols[ptr];
        
        // Check if the current character is a space or slash
        if (ptr_char == ' ' || ptr_char == '/') {
            // If buffer is not empty, decode it
            if (!buffer.empty()) {
                if (morse_dict.find(buffer) != morse_dict.end()) {
                    decoded_text += morse_dict[buffer];
                } else {
                    decoded_text += "[" + buffer + "]";
                }
                buffer = "";
            }
            
            // If it's a slash, add a space to represent word break
            if (ptr_char == '/') {
                decoded_text += " ";
            }
        } else {
            // Add the character to the buffer
            buffer += ptr_char;
        }
        
        // Increment pointer to move to next character
        ptr++;
    }
    
    // Handle any remaining buffer at the end
    if (!buffer.empty()) {
        if (morse_dict.find(buffer) != morse_dict.end()) {
            decoded_text += morse_dict[buffer];
        } else {
            decoded_text += "[" + buffer + "]";
        }
    }    
    return decoded_text;
}


int main()
{
    std::string filename = "morse_code.wav";

    // std::vector<float>, audio_data = open_wav(filename);
    std::vector<float> audio_data;
    float sample_rate = 0.0f;
    std::tie(audio_data, sample_rate) = open_wav(filename);
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

    
    
    float top_frequency = compute_top_frequency(audio_data, 44100);
    float T_sinewave = 1.0f / top_frequency;
    float num_points_sinewave = T_sinewave * sample_rate;
    std::cout << "Number of points in sinewave: " << num_points_sinewave << std::endl;
    std::cout << "Sinewave period: " << T_sinewave << " seconds" << std::endl;
    
    // remove end blank
    remove_end_blanks(audio_data);
    
    // viewing first 100 samples 
    std::cout << "First 100 samples:" << std::endl;
    for (int i = 0; i < 100 && i < audio_data.size(); i++)
    {
        std::cout << audio_data[i] << " ";
    }
    std::cout << std::endl;


    float threshold = avg_amplitude * 0.5f; // Adjust this threshold as needed
    std::vector<int> above_threshold_noend(audio_data.size());
    for (size_t i = 0; i < audio_data.size(); i++) {
        above_threshold_noend[i] = (std::abs(audio_data[i]) > threshold) ? 1 : 0;
    }
    
    std::vector<int> signal_or_blank;
    int blank_counter = 0;
    int blank_threshold = num_points_sinewave * 20; // 20 wavelengths
    
    std::cout << "labelling a range of values as signal-OFF if the number of consecutive null samples exceeds " 
              << blank_threshold << std::endl;
    
              for (size_t i = 0; i < above_threshold_noend.size(); i++) {
                if (above_threshold_noend[i] == 1) {
                    if (blank_counter == 0) {
                        signal_or_blank.push_back(1);
                    } else {
                        if (blank_counter > blank_threshold) {
                            // Add zeros for long blank
                            for (int j = 0; j < blank_counter; j++) {
                                signal_or_blank.push_back(0);
                            }
                        } else {
                            // Add ones for short blank (treat as signal)
                            for (int j = 0; j < blank_counter; j++) {
                                signal_or_blank.push_back(1);
                            }
                        }
                    }
                    blank_counter = 0;
                } else if (above_threshold_noend[i] == 0) {
                    blank_counter++;
                }
            }
    if (blank_counter > 0) {
        if (blank_counter > blank_threshold) {
            // Add zeros for long blank
            for (int j = 0; j < blank_counter; j++) {
                signal_or_blank.push_back(0);
            }
        } else {
            // Add ones for short blank (treat as signal)
            for (int j = 0; j < blank_counter; j++) {
                signal_or_blank.push_back(1);
            }
        }
    }
    print_vector(signal_or_blank, "Signal or Blank", 100);
    // Print statistics about signal_or_blank
    int ones_count = std::count(signal_or_blank.begin(), signal_or_blank.end(), 1);
    int zeros_count = std::count(signal_or_blank.begin(), signal_or_blank.end(), 0);
    std::cout << "Signal/blank analysis:" << std::endl;
    std::cout << "  Total points: " << signal_or_blank.size() << std::endl;
    std::cout << "  Signal points (1): " << ones_count << std::endl;
    std::cout << "  Blank points (0): " << zeros_count << std::endl;
    // Add this after the signal_or_blank processing section in main()

    // Analyze signal blocks and durations
    std::vector<int> durations;
    std::vector<int> block_values;
    size_t ptr = 0;
    int current_value = signal_or_blank[0];
    int num_conseq = 0;

    while (ptr < signal_or_blank.size() - 1) {
        ptr++;
        if (signal_or_blank[ptr] == current_value) {
            num_conseq++;
        } else {
            block_values.push_back(current_value);
            durations.push_back(num_conseq);
            num_conseq = 0;
            current_value = signal_or_blank[ptr];
            num_conseq++;
        }
    }
    durations.push_back(num_conseq);
    block_values.push_back(current_value);

    // Print information about the blocks
    std::cout << "\nDetected " << durations.size() << " blocks of signals/blanks" << std::endl;
    std::cout << "First few blocks:" << std::endl;
    int display_limit = std::min(10, static_cast<int>(durations.size()));
    for (int i = 0; i < display_limit; i++) {
        std::cout << "  Block " << i << ": value=" << block_values[i] 
                << ", duration=" << durations[i] << std::endl;
    }

    // Find the minimum duration (unit length)
    int min_duration = *std::min_element(durations.begin(), durations.end());
    std::cout << "\nMinimum duration (unit length): " << min_duration << std::endl;

    // Calculate the number of unit lengths for each duration
    std::vector<int> num_unit_lens;
    for (int duration : durations) {
        float ratio = static_cast<float>(duration) / min_duration;
        int rounded = std::round(ratio);
        num_unit_lens.push_back(rounded);
    }
    print_vector(num_unit_lens, "Number of unit lengths", 100);

    // Find unique unit lengths (should be 3 or 4)
    std::set<int> unique_unit_lens(num_unit_lens.begin(), num_unit_lens.end());
    std::cout << "Number of unique unit lengths: " << unique_unit_lens.size() << std::endl;
    std::cout << "Unique unit lengths: ";
    for (int len : unique_unit_lens) {
        std::cout << len << " ";
    }
    std::cout << std::endl;
    // Sort the unique unit lengths
    std::vector<int> sorted_unique_unit_lens(unique_unit_lens.begin(), unique_unit_lens.end());
    std::sort(sorted_unique_unit_lens.begin(), sorted_unique_unit_lens.end());

    int len_unit = sorted_unique_unit_lens[0]; // The smallest unit length
    int len_letter_space_or_dash = sorted_unique_unit_lens[1]; // The second smallest unit length
    int len_word_space = sorted_unique_unit_lens[-1]; // The largest unit length 

    // Generate the Morse code symbols
    std::string symbols1 = "";
    for (size_t i = 0; i < num_unit_lens.size(); i++) {
        int block = num_unit_lens[i];
        int block_value = block_values[i];
        
        if (block == len_unit) {
            if (block_value == 1) {
                symbols1 += ".";
            }
            // If block_value is 0, add nothing (as in Python)
        }
        else if (block == len_letter_space_or_dash) {
            if (block_value == 1) {
                symbols1 += "-";
            }
            else {
                symbols1 += " ";
            }
        }
        else if (block == len_word_space && block_value == 0) {
            symbols1 += "/";
        }
    }
    std::cout << "Morse symbols found: " << symbols1 << std::endl;


    std::string decoded_text = decoding(symbols1);
    std::cout << "Decoded text: " << decoded_text << std::endl;
    return 0;
}