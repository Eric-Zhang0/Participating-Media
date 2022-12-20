//
//  material.h
//  Ray_Tracking
//
//  Created by 張皓珂 on 2022/01/01.
//

#ifndef material_h
#define material_h

#include "ray_tracking.h"
#include "hittable_list.h"
#include <algorithm>

struct hit_record;
struct sphere;
struct hittable_list;

class material {
    public:
        virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered, const shared_ptr<sphere>& center) const = 0;
};

inline double deltaTracking(){
    double sigma = 1.0;
    double random_number = random_double();
    
    double t = log( 1.0 - random_number ) / sigma;
    return t;
}

inline double getRadius(const shared_ptr<sphere>& s){
    return (*s).radius;
}

inline point3 getCenter(const shared_ptr<sphere>& s){
    return (*s).center;
}

inline double getTmax(const ray& r, const hit_record& rec, const shared_ptr<sphere>& s){
    double radius = getRadius(s);
    double direction_length = r.direction().length();
    double t_in = rec.t;
    
    return (2*radius/direction_length)+t_in;
}

inline double getTmin(const ray& r, const shared_ptr<sphere>& s){
    double radius = getRadius(s);
    point3 center = getCenter(s);
    double direction_length = r.direction().length();
    
    return (center.length() - radius) / direction_length;
}

inline double getTObj(const ray& r, const hit_record& rec, const shared_ptr<sphere>& s){
    double t_in = rec.t;
    point3 center = getCenter(s);
    vec3 pc_in = r.at(t_in) - center;
    double pc_in_length = pc_in.length();
    double direction_length = r.direction().length();
    
    double cos_angle = dot(r.direction(), pc_in) / (pc_in_length*direction_length);
    
    double t_obj = (2*getRadius(s)*cos_angle)/direction_length + t_in;
    return t_obj;
}

inline double getGrayValue(const color& c){
    double grayValue = c.x()*0.3 + c.y()*0.59 + c.z()*0.11;
    return grayValue;
}

inline double transmittance(double t){
    double sigma = 1.0;
    return exp(-sigma*t);
}

inline double luminosity(double t, const color& c){
    return getGrayValue(c)*transmittance(t)*(1/(4*pi));
}

inline double accept(double t1, double t2, const color& c){
    double grayValueT1 = luminosity(t1, c);
    double grayValueT2 = luminosity(t2, c);
    
    double a = std::min(1.0, grayValueT2/grayValueT1);
    return a;
}

inline double mutate(const ray& r, const hit_record& rec, const shared_ptr<sphere>& s){
    auto t_max = getTmax(r, rec, s);
    auto t_min = getTmin(r, s);
    double t = rec.t;
    double epsilon = 0.1;
    double length = t_max - t_min;
    
    t += epsilon * (2 * length * random_double() - length);
    
    if (t > t_max) {
        t -= length;
    }else if(t < t_min){
        t += length;
    }
    return t;
}
class metal : public material {
    public:
        metal(const color& a, double f) : albedo(a), fuzz(f < 1 ? f : 1) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered, const shared_ptr<sphere>& center
        ) const override {
            vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
            scattered = ray(rec.p, reflected + fuzz*random_in_unit_sphere());
            attenuation = albedo;
            return (dot(scattered.direction(), rec.normal) > 0);
        }

    public:
        color albedo;
        double fuzz;
};

class delta : public material {
    public:
        delta(const color& a) : albedo(a){}
        
    virtual bool scatter(
        const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered, const shared_ptr<sphere>& center
    ) const override {
        auto t_obj = getTObj(r_in, rec, center);
        auto t_dt = 0.0;
        auto t_max = getTmax(r_in, rec, center);
        
        auto r = random_double();
        while (true) {
            t_dt -= deltaTracking();
            
            if (t_dt >= t_max) {
                break;
            }
            if((t_dt/t_obj) > r){
                auto scatter_direction = rec.normal + random_unit_vector();
                
                // Catch degenerate scatter direction
                if (scatter_direction.near_zero())
                    scatter_direction = rec.normal;
                
                point3 p = r_in.at(t_dt);
                scattered = ray(p, scatter_direction);
                attenuation = albedo;
            }
            r = random_double();
        }
        return true;
    }
    public:
        color albedo;
};

class ratio : public material {
    public:
        ratio(const color& a) : albedo(a){}
        
    virtual bool scatter(
        const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered, const shared_ptr<sphere>& center
    ) const override {
        auto t_obj = getTObj(r_in, rec, center);
        auto t_dt = 0.0;
        auto t_max = getTmax(r_in, rec, center);
        auto T = 1.0;
        
        while (true) {
            t_dt -= deltaTracking();
            
            if (t_dt >= t_obj) break;
            T *= 1 - (t_dt/t_max);
            if(t_dt < t_obj){
                auto scatter_direction = rec.normal + random_unit_vector();
                
                // Catch degenerate scatter direction
                if (scatter_direction.near_zero())
                    scatter_direction = rec.normal;
                
                point3 p = r_in.at(t_dt);
                scattered = ray(p, scatter_direction);
                attenuation = albedo;
            }
        }
        return true;
    }
    public:
        color albedo;
};

class metropolis : public material {
    public:
        metropolis(const color& a) : albedo(a) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered, const shared_ptr<sphere>& center
        ) const override {
            auto scatter_direction = rec.normal + random_unit_vector();
            
            int n_mutations = 1000;
            point3 p;
            
            for (int i = 0; i != n_mutations; ++i) {
                double t_c = mutate(r_in, rec, center);
                double a = accept(rec.t, t_c, albedo);
                double r = random_double();
                
                if(r < a) p = r_in.at(t_c);
                
                scattered = ray(p, scatter_direction);
                attenuation = albedo;
            }
            return true;
        }

    public:
        color albedo;
};

class dtMixedmlt : public material {
    public:
        dtMixedmlt(const color& a) : albedo(a) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered, const shared_ptr<sphere>& center
        ) const override {
            auto scatter_direction = rec.normal + random_unit_vector();
            int n_mutations = 400;
            point3 p;
            auto t_obj = getTObj(r_in, rec, center);
            auto t_dt = 0.0;
            auto t_max = getTmax(r_in, rec, center);
            auto temp = t_dt;
            
            auto rand = random_double();
            while (true) {
                temp = t_dt;
                t_dt -= deltaTracking();
                if (t_dt >= t_max){
                    t_dt = temp;
                    break;
                }
                if((t_dt/t_obj) > rand){
                    auto scatter_direction = rec.normal + random_unit_vector();
                    
                    // Catch degenerate scatter direction
                    if (scatter_direction.near_zero())
                        scatter_direction = rec.normal;
                    
                    point3 p = r_in.at(t_dt);
                    scattered = ray(p, scatter_direction);
                    attenuation = albedo;
                }
                rand = random_double();
            }
            
            for (int i = 0; i != n_mutations; ++i) {
                double t_c = mutate(r_in, rec, center);
                double a = accept(t_dt, t_c, albedo);
                double r = random_double();
                
                if (r < a) p = r_in.at(t_c);
                
                scattered = ray(p, scatter_direction);
                attenuation = albedo;
            }
            return true;
        }

    public:
        color albedo;
};

class lambertian : public material {
    public:
        lambertian(const color& a) : albedo(a) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered, const shared_ptr<sphere>& center
        ) const override {
            auto scatter_direction = rec.normal + random_unit_vector();

            // Catch degenerate scatter direction
            if (scatter_direction.near_zero())
                scatter_direction = rec.normal;

            scattered = ray(rec.p, scatter_direction);
            attenuation = albedo;
            return true;
        }

    public:
        color albedo;
};

class dielectric : public material {
    public:
        dielectric(double index_of_refraction) : ir(index_of_refraction) {}

        virtual bool scatter(
            const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered, const shared_ptr<sphere>& center
        ) const override {
            attenuation = color(1.0, 1.0, 1.0);
            double refraction_ratio = rec.front_face ? (1.0/ir) : ir;

            vec3 unit_direction = unit_vector(r_in.direction());
            double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
            double sin_theta = sqrt(1.0 - cos_theta*cos_theta);

            bool cannot_refract = refraction_ratio * sin_theta > 1.0;
            vec3 direction;

            if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double())
                direction = reflect(unit_direction, rec.normal);
            else
                direction = refract(unit_direction, rec.normal, refraction_ratio);

            scattered = ray(rec.p, direction);
            return true;
        }

    public:
        double ir; // Index of Refraction

    private:
        static double reflectance(double cosine, double ref_idx) {
            // Use Schlick's approximation for reflectance.
            auto r0 = (1-ref_idx) / (1+ref_idx);
            r0 = r0*r0;
            return r0 + (1-r0)*pow((1 - cosine),5);
        }
};

#endif /* material_h */
