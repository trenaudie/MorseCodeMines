#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <set>
#include "utils.h"
#include <map>
#include <sstream>
#include <tuple>

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
    
    std::cout << "Decoded text: " << decoded_text << std::endl;
    
    return decoded_text;
}

int main(){
    std::string symbols = ".... . .-.. .-.. --- / .-- --- .-. .-.. -.."; //HELLO WORLD
    std::string decoded_text = decoding(symbols);
}