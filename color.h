//
//  color.h
//  Ray_Tracking
//
//  Created by 張皓珂 on 2021/12/29.
//

#ifndef COLOR_H
#define COLOR_H

#include "vec3.h"

#include <iostream>
#include <fstream>

void write_color(std::ostream& out, color pixel_color, int samples_per_pixel, std::ofstream& file) {
    auto r = pixel_color.x();
    auto g = pixel_color.y();
    auto b = pixel_color.z();

    // Divide the color by the number of samples and using gamma-correction for gamma = 2.0.
    auto scale = 1.0 / samples_per_pixel;
    r = sqrt(scale * r);
    g = sqrt(scale * g);
    b = sqrt(scale * b);

    // Write the translated RGB[0,255] value of each color component.
    // out << static_cast<int>(256 * clamp(r, 0.0, 0.999)) << ' '
        // << static_cast<int>(256 * clamp(g, 0.0, 0.999)) << ' '
        // << static_cast<int>(256 * clamp(b, 0.0, 0.999)) << '\n';
    
    // Write the translated RGB[0,255] value of each color component to output image
    file << static_cast<int>(256 * clamp(r, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(g, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(b, 0.0, 0.999)) << '\n';
}

#endif /* color_h */
