#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>

// Function to print vector contents with an optional limit
template <typename T>
void print_vector(const std::vector<T>& vec, const std::string& label = "Vector", size_t limit = 0);

#endif // UTILS_H