#ifndef COLOR_H
#define COLOR_H

#include <array>
#include <string_view>

constexpr std::array<std::string_view, 12> colors{
    "LightCoral", "DarkMagenta", "DarkOrange",  "Yellow",
    "DarkKhaki",  "Lavender",    "LightPink",   "LawnGreen",
    "SteelBlue",  "Teal",        "SaddleBrown", "Cyan",
};

inline float opacity(int n) { return n % 2 ? 0.3f : 1.0f; }

#endif