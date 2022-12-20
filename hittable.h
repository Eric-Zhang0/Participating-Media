//
//  hittable.h
//  Ray_Tracking
//
//  Created by 張皓珂 on 2021/12/31.
//

#ifndef hittable_h
#define hittable_h

#include "ray.h"
#include "ray_tracking.h"


class material;

struct hit_record {
    point3 p;
    vec3 normal;
    double t;
    bool front_face;
    shared_ptr<material>mat_ptr;

    inline void set_face_normal(const ray& r, const vec3& outward_normal) {
        front_face = dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal :-outward_normal;
    }
};

class hittable {
    public:
        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
};

#endif /* hittable_h */
