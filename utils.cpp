#include <iostream>
#include <vector>
#include <string>
#include "utils.h"

// Function to print a vector of any printable type
template <typename T>
void print_vector(const std::vector<T>& vec, const std::string& label, size_t limit) {
    std::cout << label << " (size: " << vec.size() << "):" << std::endl;
    
    size_t print_count = (limit > 0 && limit < vec.size()) ? limit : vec.size();
    
    for (size_t i = 0; i < print_count; ++i) {
        std::cout << vec[i];
        
        // Add comma for all but the last element
        if (i < print_count - 1) {
            std::cout << ", ";
        }
        
        // Add line break every 10 elements for better readability
        if ((i + 1) % 10 == 0) {
            std::cout << std::endl;
        }
    }
    
    // Final newline if not already added
    if (print_count % 10 != 0) {
        std::cout << std::endl;
    }
    
    // Indicate if we truncated the output
    if (print_count < vec.size()) {
        std::cout << "... (and " << (vec.size() - print_count) << " more elements)" << std::endl;
    }
    
    std::cout << std::endl;
}

// Explicit template instantiations for common types
template void print_vector<int>(const std::vector<int>&, const std::string&, size_t);
template void print_vector<float>(const std::vector<float>&, const std::string&, size_t);
template void print_vector<double>(const std::vector<double>&, const std::string&, size_t);
template void print_vector<std::string>(const std::vector<std::string>&, const std::string&, size_t);
template void print_vector<char>(const std::vector<char>&, const std::string&, size_t);