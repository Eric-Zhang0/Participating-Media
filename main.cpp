//
//  main.cpp
//  Ray_Tracking
//
//  Created by 張皓珂 on 2021/12/29.
//

#include "ray_tracking.h"
#include "camera.h"
#include "color.h"
#include "hittable_list.h"
#include "sphere.h"
#include "material.h"

#include <iostream>
#include <fstream>
#include <chrono>
#include <future>

color ray_color(const ray& r, const hittable& world, const shared_ptr<sphere>& s ,int depth) {
    hit_record rec;
    
    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0,0,0);
    
    if (world.hit(r, 0.001, infinity, rec)) {
            ray scattered;
            color attenuation;
            if (rec.mat_ptr->scatter(r, rec, attenuation, scattered, s))
                return attenuation * ray_color(scattered, world, s, depth-1);
            return color(0,0,0);
    }
    
    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*color(1.0, 1.0, 1.0) + t*color(0.5, 0.7, 1.0);
}

hittable_list random_scene() {
    hittable_list world;

    auto ground_material = make_shared<lambertian>(color(0.5, 0.5, 0.5));
    world.add(make_shared<sphere>(point3(0,-1000,0), 1000, ground_material));

    for (int a = -5; a < 5; a++) {
        for (int b = -5; b < 5; b++) {
            auto choose_mat = random_double();
            point3 center(a + 0.9*random_double(), 0.2, b + 0.9*random_double());

            if ((center - point3(4, 0.2, 0)).length() > 0.9) {
                shared_ptr<material> sphere_material;

                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = color::random() * color::random();
                    sphere_material = make_shared<lambertian>(albedo);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                } else if (choose_mat < 0.95) {
                    // metal
                    // auto albedo = color::random(0.5, 1);
                    // auto fuzz = random_double(0, 0.5);
                    // sphere_material = make_shared<metal>(albedo, fuzz);
                    // world.add(make_shared<sphere>(center, 0.2, sphere_material));
                    auto albedo = color::random() * color::random();
                    sphere_material = make_shared<lambertian>(albedo);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                } else {
                    // glass
                    // sphere_material = make_shared<dielectric>(1.5);
		    auto albedo = color::random() * color::random();
                    sphere_material = make_shared<lambertian>(albedo);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }

    // auto material1 = make_shared<dielectric>(1.5);
    // world.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material1));

    auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
    world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));

    // auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
    // world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

    return world;
}

int main(int argc, const char * argv[]) {
    // insert code here...
    
    // Image

    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 400;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int max_depth = 50;
    
    int samplesPerPixel = argc == 2 ? atoi(argv[1]): 100;

    // World

    hittable_list world = random_scene();
    
    //auto material_ground = make_shared<lambertian>(color(0.8, 0.8, 0.0));
    auto material_center = make_shared<lambertian>(color(0.4, 0.2, 0.1));
    //auto material_left   = make_shared<metal>(color(0.8, 0.8, 0.8));
    //auto material_right  = make_shared<metal>(color(0.8, 0.6, 0.2));

    //world.add(make_shared<sphere>(point3( 0.0, -100.5, -1.0), 100.0, material_ground));
    //world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material_center));
    //world.add(make_shared<sphere>(point3(-1.0,    0.0, -1.0),   0.5, material_left));
    //world.add(make_shared<sphere>(point3( 1.0,    0.0, -1.0),   0.5, material_right));
    
    shared_ptr<sphere> center = make_shared<sphere>(point3(-4, 1, 0), 1.0, material_center);

    // Camera
    point3 lookfrom(13,2,3);
    point3 lookat(0,0,0);
    vec3 vup(0,1,0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.1;
    
    camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus);

    // Render
    std::cout << "P3\n" << image_width << " " << image_height << "\n255\n";
    std::ofstream file("image/reference_100.ppm");
    file << "P3\n"
         << image_width << ' ' << image_height << "\n255\n";

    auto start = std::chrono::system_clock::now();
    
    #pragma omp parallel
    {
        #pragma omp single
        {
    for (int j = 0; j < image_height; j++) {
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samplesPerPixel, 100. * j / (image_height - 1));
            for (int i = 0; i < image_width; i++) {
                #pragma omp task
                {
                color pixel_color(0, 0, 0);
                for (int s = 0; s < samplesPerPixel; ++s) {
                    auto u = (i + random_double()) / (image_width-1);
                    auto v = (image_height - 1 - j + random_double()) / (image_height-1);
                    ray r = cam.get_ray(u, v);
                    pixel_color += ray_color(r, world, center, max_depth);
                }
                write_color(std::cout, pixel_color, samplesPerPixel, file);
                }
                    }
            }
        }
    }
    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout <<  "The whole rendering costs "
         << double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den
         << " seconds." << std::endl;
    return 0;
}

