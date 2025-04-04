
#include <iostream>
#include <string>
#include <unordered_map>
#include <cctype>
#include <fstream> // For std::ofstream
#include <vector>  // For std::vector
#include <cmath>

using TAlphabet = std::unordered_map<std::string, char>;
const TAlphabet MorseCode = {
    {
        "01",
        'A',
    },
    {
        "1000",
        'B',
    },
    {
        "1010",
        'C',
    },
    {
        "100",
        'D',
    },
    {
        "0",
        'E',
    },
    {
        "0010",
        'F',
    },
    {
        "110",
        'G',
    },
    {
        "0000",
        'H',
    },
    {
        "00",
        'I',
    },
    {
        "0111",
        'J',
    },
    {
        "101",
        'K',
    },
    {
        "0100",
        'L',
    },
    {
        "11",
        'M',
    },
    {
        "10",
        'N',
    },
    {
        "111",
        'O',
    },
    {
        "0110",
        'P',
    },
    {
        "1101",
        'Q',
    },
    {
        "010",
        'R',
    },
    {
        "000",
        'S',
    },
    {
        "1",
        'T',
    },
    {
        "001",
        'U',
    },
    {
        "0001",
        'V',
    },
    {
        "011",
        'W',
    },
    {
        "1001",
        'X',
    },
    {
        "1011",
        'Y',
    },
    {
        "1100",
        'Z',
    },
    {
        "01111",
        '1',
    },
    {
        "00111",
        '2',
    },
    {
        "00011",
        '3',
    },
    {
        "00001",
        '4',
    },
    {
        "00000",
        '5',
    },
    {
        "10000",
        '6',
    },
    {
        "11000",
        '7',
    },
    {
        "11100",
        '8',
    },
    {
        "11110",
        '9',
    },
    {
        "11111",
        '0',
    },
    {
        "010101",
        '.',
    },
    {
        "110011",
        ',',
    },
    {
        "001100",
        '?',
    },
    {
        "011110",
        '\'',
    },
    {
        "101011",
        '!',
    },
    {
        "10010",
        '/',
    },
    {
        "10110",
        '(',
    },
    {
        "101101",
        ')',
    },
    {
        "01000",
        '&',
    },
    {
        "111000",
        ':',
    },
    {
        "101010",
        ';',
    },
    {
        "10001",
        '=',
    },
    {
        "01010",
        '+',
    },
    {
        "100001",
        '-',
    },
    {
        "001101",
        '_',
    },
    {
        "010010",
        '"',
    },
    {
        "0001001",
        '$',
    },
    {
        "011010",
        '@',
    },
};

// Data format of the audio samples
typedef enum
{
    MORSE_SAMPLE_FORMAT_UNDEFINED,
    MORSE_SAMPLE_FORMAT_U8,
    MORSE_SAMPLE_FORMAT_I8,
    MORSE_SAMPLE_FORMAT_U16,
    MORSE_SAMPLE_FORMAT_I16,
    MORSE_SAMPLE_FORMAT_F32,
} ggmorse_SampleFormat;

int bytesForSampleFormat(ggmorse_SampleFormat sample_format)
{
    switch (sample_format)
    {
    case MORSE_SAMPLE_FORMAT_UNDEFINED:
        return 0;
        break;
    case MORSE_SAMPLE_FORMAT_U8:
        return sizeof(uint8_t);
        break;
    case MORSE_SAMPLE_FORMAT_I8:
        return sizeof(int8_t);
        break;
    case MORSE_SAMPLE_FORMAT_U16:
        return sizeof(uint16_t);
        break;
    case MORSE_SAMPLE_FORMAT_I16:
        return sizeof(int16_t);
        break;
    case MORSE_SAMPLE_FORMAT_F32:
        return sizeof(float);
        break;
    };

    fprintf(stderr, "Invalid sample format: %d\n", (int)sample_format);

    return 0;
}

// lendot_ms returns the duration of a unit in ms
// duration of a unit (in samples)
float lendot_ms(float speed_wpm)
{
    // ms in a minute divided by speed in units per minute
    // ie ms per unit
    return 60000.0f / (50.0f * speed_wpm);
}

std::vector<float> encoding(std::string message, ggmorse_SampleFormat sample_format, int speed_wpm, int sample_rate_out)
{
    std::string symbols0; // 01234 encoding
    std::string symbols1; // dots and dashes encoding
    int num_samples;
    int nsamples_unit = sample_rate_out * 1e-3 * lendot_ms(speed_wpm);
    int nsamples_unit_space = nsamples_unit;
    int nsamples_letter_space = 3.0f * nsamples_unit_space;
    int nsamples_word_space = 7.0f * nsamples_unit_space;
    for (int i = 0; i < message.size(); ++i)
    {
        if (message[i] == ' ')
        {
            continue;
        }
        for (const auto &l : MorseCode)
        {
            if (std::tolower(l.second) == message[i])
            {
                std::cout << "found " << message[i] << std::endl;
                // letter found, lets add it to the symbols stack
                for (int k = 0; k < (int)l.first.size(); ++k)
                {
                    // l.first[k] is the 0 or 1
                    char letter_unit = l.first[k];
                    if (letter_unit == '0')
                    {
                        symbols0 += "0";
                        symbols1 += ".";
                        num_samples += 1 * nsamples_unit;
                    }
                    if (letter_unit == '1')
                    {
                        symbols0 += "1";
                        symbols1 += "_";
                        num_samples += 3 * nsamples_unit;
                    }
                    if (k < (int)l.first.size() - 1)
                    {
                        symbols0 += "2";
                        symbols1 += "";
                        num_samples += nsamples_unit_space;
                    }
                }
                break; // breaking to not try to match any other letters
            }
        }
        // one letter has been added
        // lets add a letter spacing, only if this was not the last letter
        // IF SAME WORD, but not last letter of word
        if (i < message.size() - 1)
        {
            // add 3 if inter word, 4 if between word
            if (
                message[i + 1] == ' ')
            {
                symbols0 += "4";
                symbols1 += "/";
                num_samples += nsamples_word_space;
            }
            else
            {
                symbols0 += "3";
                symbols1 += " ";
                num_samples += nsamples_letter_space;
            }
        }
    }
    std::cout << std::endl;
    std::cout << symbols0 << std::endl;
    std::cout << symbols1 << std::endl;
    std::cout << num_samples << std::endl;

    // Sine wave encoding
    const float frequency_hz = 550.0f; // Frequency of the signal
    std::vector<float> v(num_samples, 0.0f);
    int idx = 0;

    // We want dampFactor * num_samples to be 0.1 second
    float dampFactor = 1.0f / std::max(1.0f, 0.1f * nsamples_unit);
    float current_factor = 0.0f;

    for (int i = 0; i < symbols0.size(); ++i)
    {
        char s = symbols0[i];
        int snum_samples = 0;

        if (s == '0')
            snum_samples = nsamples_unit;
        else if (s == '1')
            snum_samples = 3 * nsamples_unit;
        else if (s == '2')
            snum_samples = nsamples_unit_space;
        else if (s == '3')
            snum_samples = nsamples_letter_space;
        else if (s == '4')
            snum_samples = nsamples_word_space;

        for (int j = 0; j < snum_samples && idx < num_samples; ++j)
        {
            // Generate sine wave
            float t = static_cast<float>(idx) / sample_rate_out;

            if (s == '0' || s == '1')
            {
                // For dots and dashes, generate tone with increasing amplitude
                v[idx] = current_factor * std::sin(2.0f * M_PI * frequency_hz * t);
                current_factor = std::min(1.0f, current_factor + dampFactor);
            }
            else
            {
                // For spaces, silence with decreasing amplitude
                v[idx] = 0.0f;
                current_factor = std::max(0.0f, current_factor - dampFactor);
            }

            idx++;
        }
    }
    return v;
}
void writeWavHeader(std::ofstream &file, const std::vector<float> &v, int sample_rate_out,
                    const int num_channels, int bits_per_sample)
{

    const int bytes_per_sample = bits_per_sample / 8;
    const int byte_rate = sample_rate_out * num_channels * bytes_per_sample;
    const int block_align = num_channels * bytes_per_sample;
    const int data_size = v.size() * bytes_per_sample;
    const int file_size = 36 + data_size;

    file.write("RIFF", 4);
    file.write(reinterpret_cast<const char *>(&file_size), 4);
    file.write("WAVE", 4);
    file.write("fmt ", 4);

    int subchunk1_size = 16;
    int audio_format = 1; // PCM

    file.write(reinterpret_cast<const char *>(&subchunk1_size), 4);
    file.write(reinterpret_cast<const char *>(&audio_format), 2);
    file.write(reinterpret_cast<const char *>(&num_channels), 2);
    file.write(reinterpret_cast<const char *>(&sample_rate_out), 4);
    file.write(reinterpret_cast<const char *>(&byte_rate), 4);
    file.write(reinterpret_cast<const char *>(&block_align), 2);
    file.write(reinterpret_cast<const char *>(&bits_per_sample), 2);

    file.write("data", 4);
    file.write(reinterpret_cast<const char *>(&data_size), 4);

    for (size_t i = 0; i < v.size(); ++i)
    {
        if (bits_per_sample == 8)
        {
            // 8-bit audio is unsigned [0, 255]
            uint8_t sample = static_cast<uint8_t>((v[i] + 1.0f) * 127.5f);
            file.write(reinterpret_cast<const char *>(&sample), 1);
        }
        else if (bits_per_sample == 16)
        {
            // 16-bit audio is signed [-32768, 32767]
            int16_t sample = static_cast<int16_t>(v[i] * 32767.0f);
            file.write(reinterpret_cast<const char *>(&sample), 2);
        }
        else if (bits_per_sample == 24)
        {
            // 24-bit audio is signed, 3 bytes
            int32_t sample_32 = static_cast<int32_t>(v[i] * 8388607.0f); // 2^23 - 1
            unsigned char bytes[3];
            bytes[0] = sample_32 & 0xFF;
            bytes[1] = (sample_32 >> 8) & 0xFF;
            bytes[2] = (sample_32 >> 16) & 0xFF;
            file.write(reinterpret_cast<const char *>(bytes), 3);
        }
        else if (bits_per_sample == 32)
        {
            // 32-bit audio is typically floating point
            int32_t sample = static_cast<int32_t>(v[i] * 2147483647.0f); // 2^31 - 1
            file.write(reinterpret_cast<const char *>(&sample), 4);
        }
    }
}

int main()
{
    // input pararms
    const std::string message = "hi in";
    const ggmorse_SampleFormat sample_format = MORSE_SAMPLE_FORMAT_F32;
    const int speed_wpm = 10;
    const int sample_rate_out = 44100;
    const int sample_rate_in = 44100;
    auto v = encoding(message, sample_format, speed_wpm, sample_rate_out);
    std::cout << "Encoded message! " << message << std::endl;

    // saving the vector itself
    std::cout << "Saving the vector to a file" << std::endl;
    std::ofstream outFile("morse_output.txt");
    for (const auto &sample : v)
    {
        outFile << sample << "\n";
    }
    outFile.close();
    std::cout << "Vector saved to morse_output.txt" << std::endl;

    std::ofstream wavFile("temp_out.wav", std::ios::binary);
    writeWavHeader(wavFile, v, sample_rate_out, 1, 16);
    return 0;
}
