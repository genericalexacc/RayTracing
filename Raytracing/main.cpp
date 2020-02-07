//
//  main.cpp
//  Raytracing
//
//  Created by Alexander Shevchenko on 2020-02-06.
//  Copyright © 2020 Alexander Shevchenko. All rights reserved.
//

#include <iostream>
#include <math.h>
#include "float.h"

#include <stdlib.h>
#include <fstream>
#include <cstdlib>



using namespace std;

class vec3 {
public:
    vec3(){}
    
    vec3(float e1, float e2, float e3){ e[0]=e1; e[1]=e2; e[2]=e3; };
    
    inline float x() const { return e[0]; }
    inline float y() const { return e[1]; }
    inline float z() const { return e[2]; }
    
    inline float r() const { return e[0]; }
    inline float g() const { return e[1]; }
    inline float b() const { return e[2]; }
    
    inline const vec3& operator+() const { return *this; }
    inline vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
    inline float operator[](int i) const { return e[i]; }
    inline float& operator[](int i) { return e[i]; }
    
    inline vec3& operator+=(const vec3 &v2);
    inline vec3& operator-=(const vec3 &v2);
    inline vec3& operator*=(const vec3 &v2);
    inline vec3& operator/=(const vec3 &v2);
    inline vec3& operator*=(const float t);
    inline vec3& operator/=(const float t);
    
    inline float length() const{
        return sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
    }
    inline float squaredLength() const{
        return e[0]*e[0]+e[1]*e[1]+e[2]*e[2];
    }
    inline void make_unit_vector();
    
    float e[3];
};

inline std::istream& operator>>(std::istream &is, vec3 &t){
    is >> t.e[0] >> t.e[1] >> t.e[2];
    return is;
}

inline std::ostream& operator>>(std::ostream &os, vec3 &t){
    os << t.e[0] << " " << t.e[1] << " " << t.e[2];
    return os;
}

inline void vec3::make_unit_vector(){
    float k = 1.0 / sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
    e[0] *= k; e[1] *= k; e[2] *= k;
}

inline vec3 operator+(const vec3 &v1, const vec3 &v2){
    return vec3(v1.e[0]+v2.e[0], v1.e[1]+v2.e[1], v1.e[2]+v2.e[2]);
}

inline vec3 operator-(const vec3 &v1, const vec3 &v2){
    return vec3(v1.e[0]-v2.e[0], v1.e[1]-v2.e[1], v1.e[2]-v2.e[2]);
}

inline vec3 operator*(const vec3 &v1, const vec3 &v2){
    return vec3(v1.e[0]*v2.e[0], v1.e[1]*v2.e[1], v1.e[2]*v2.e[2]);
}

inline vec3 operator/(const vec3 &v1, const vec3 &v2){
    return vec3(v1.e[0]/v2.e[0], v1.e[1]/v2.e[1], v1.e[2]/v2.e[2]);
}

inline vec3 operator*(float t, const vec3 &v){
    return vec3(t*v.e[0], t*v.e[1], t*v.e[2]);
}

inline vec3 operator/(const vec3 &v, float t){
    return vec3(v.e[0]/t, v.e[1]/t, v.e[2]/t);
}

inline vec3 operator*(const vec3 &v, float t){
    return vec3(v.e[0]*t, v.e[1]*t, v.e[2]*t);
}

inline float dot(const vec3 &v1, const vec3 &v2){
    return v1.e[0]*v2.e[0] + v1.e[1]*v2.e[1] + v1.e[2]*v2.e[2];
}

inline vec3 cross(const vec3 &v1, const vec3 &v2){
    return vec3( (v1.e[1] * v2.e[2] - v1.e[2] * v2.e[1]),
                (v1.e[0] * v2.e[2] - v1.e[2] * v2.e[0]),
                (v1.e[0] * v2.e[1] - v1.e[1] * v2.e[0]));
}

inline vec3& vec3::operator+=(const vec3 &v){
    e[0] += v.e[0];
    e[1] += v.e[1];
    e[2] += v.e[2];
    return *this;
}

inline vec3& vec3::operator-=(const vec3 &v){
    e[0] -= v.e[0];
    e[1] -= v.e[1];
    e[2] -= v.e[2];
    return *this;
}

inline vec3& vec3::operator*=(const vec3 &v){
    e[0] *= v.e[0];
    e[1] *= v.e[1];
    e[2] *= v.e[2];
    return *this;
}


inline vec3& vec3::operator*=(const float t){
    e[0] *= t;
    e[1] *= t;
    e[2] *= t;
    return *this;
}

inline vec3& vec3::operator/=(const vec3 &v){
    e[0] /= v.e[0];
    e[1] /= v.e[1];
    e[2] /= v.e[2];
    return *this;
}

inline vec3& vec3::operator/=(const float t){
    float k = 1.0 / t;
    e[0] *= k;
    e[1] *= k;
    e[2] *= k;
    return *this;
}

inline vec3 unit_vector(vec3 v){
    return v / v.length();
}



class ray{
public:
    ray(){}
    ray(const vec3& a, const vec3& b){ A = a; B = b; }
    vec3 origin() const { return A; }
    vec3 direction() const { return B; }
    vec3 point_at_parameter(float t) const {return A + t*B; }
    
    vec3 A;
    vec3 B;
};

//bool hit_sphere(const vec3& center, float radius, const ray& r){
//    vec3 oc = r.origin() - center;
//    float a = dot(r.direction(), r.direction());
//    float b = 2.0 * dot(oc, r.direction());
//    float c = dot(oc, oc) - radius*radius;
//    float discriminant = b*b - 4*a*c;
//    if (discriminant < 0) {
//        return -1.0;
//    }else{
//        return (-b - sqrt(discriminant)) / (2.0 * a);
//    }
//}

class material;

struct hit_record{
    float t;
    vec3 p;
    vec3 normal;
    material *mat_ptr;
};

class hitable{
public:
    virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const = 0;
};

inline double random_double() {
    return rand() / (RAND_MAX + 1.0);
}

vec3 random_in_unit_sphere(){
    vec3 p;
    do {
        p = 2.0 * vec3(random_double(), random_double(), random_double()) - vec3(1,1,1);
    }while (p.squaredLength() >= 1.0);
    return p;
}

class material{
public:
    virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const = 0;
};

class lambertian : public material{
public:
    lambertian(const vec3& a) : albedo(a) {}
    virtual bool scatter(const ray &r_in, const hit_record &rec, vec3 &attenuation, ray &scattered) const {
        vec3 target = rec.p + rec.normal + random_in_unit_sphere();
        scattered = ray(rec.p, target-rec.p);
        attenuation = albedo;
        return true;
    }
    vec3 albedo;
};

vec3 reflect(const vec3& v, const vec3& n){
    return v - 2*dot(v, n) * n;
}

class metal : public material {
public:
    metal(const vec3& a) : albedo(a) {}
    virtual bool scatter(const ray &r_in, const hit_record &rec, vec3 &attenuation, ray &scattered) const {
        vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
        scattered = ray(rec.p, reflected);
        attenuation = albedo;
        return (dot(scattered.direction(), rec.normal) > 0);
    }
    vec3 albedo;
};

class sphere : public hitable{
public:
    sphere(){}
    sphere(vec3 cen, float r, material *m) : center(cen), radius(r), mat_ptr(m) {};
    virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const;
    vec3 center;
    float radius;
    material *mat_ptr;
};

bool sphere::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
    vec3 oc = r.origin() - center;
    float a = dot(r.direction(), r.direction());
    float b = dot(oc, r.direction());
    float c = dot(oc, oc) - radius*radius;
    float discriminant = b*b - a*c;
    if (discriminant > 0) {
        float temp = (-b - sqrt(discriminant))/a;
        if(temp < t_max && temp > t_min){
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            rec.normal = (rec.p - center) / radius;
            rec.mat_ptr = mat_ptr;
            return true;
        }
        temp = (-b + sqrt(discriminant))/a;
        if(temp < t_max && temp > t_min){
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            rec.normal = (rec.p - center) / radius;
            rec.mat_ptr = mat_ptr;
            return true;
        }
    }
    return false;
};

class hitable_list: public hitable{
public:
    hitable_list(){};
    hitable_list(hitable **l, int n){list = l; list_size = n;}
    virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const;
    hitable **list;
    int list_size;
};

bool hitable_list::hit(const ray& r, float t_min, float t_max, hit_record& rec) const{
    hit_record temp_rec;
    bool hit_anything = false;
    double closest_so_far = t_max;
    for(int i = 0; i < list_size; i++){
        if(list[i] -> hit(r, t_min, closest_so_far, temp_rec)){
            hit_anything=true;
            closest_so_far=temp_rec.t;
            rec=temp_rec;
        }
    }
    return hit_anything;
}

class camera{
public:
    camera(){
        lower_left_corner = vec3 (-2.0, -1.0, -1.0);
        horizontal = vec3 (4.0, 0.0, 0.0);
        vertical = vec3(0.0, 2.0, 0.0);
        origin = vec3(0.0, 0.0, 0.0);
    }
    
    ray get_ray(float u, float v){
        auto r = ray(origin,
                     lower_left_corner + u*horizontal + v*vertical - origin);
        return r;
    }
    
    vec3 origin;
    vec3 lower_left_corner;
    vec3 horizontal;
    vec3 vertical;
};


vec3 color(const ray& r, hitable *world, int depth){
    hit_record rec;
    if(world -> hit(r, 0.001, MAXFLOAT, rec)){
        ray scattered;
        vec3 attenuation;
        if(depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)){
            return attenuation*color(scattered, world, depth+1);
        }else{
            return vec3(0,0,0);
        }
    }else{
        vec3 unit_direction = unit_vector(r.direction());
        float t = 0.5 * (unit_direction.y() + 1.0);
        return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
    }
}



int main(int argc, const char * argv[]) {
    int nx = 400;
    int ny = 200;
    int ns = 100;
    
    ofstream myfile;
    myfile.open ("/Users/alexandershevchenko/Desktop/out.ppm");
    myfile << "P3\n" << nx << " " << ny << "\n255\n";
    
    hitable *list[4];
    list[0] = new sphere(vec3(0, 0, -1), 0.5, new lambertian(vec3(0.8, 0.3, 0.3)));
    list[1] = new sphere(vec3(0, -100.5, -1), 100, new lambertian(vec3(0.8, 0.8, 0.3)));
    list[2] = new sphere(vec3(1,0,-1), 0.5, new metal(vec3(0.8, 0.6, 0.2)));
    list[3] = new sphere(vec3(-1,0,-1), 0.5, new lambertian(vec3(0.8, 0.8, 0.8)));
    
    hitable *world = new hitable_list(list, 4);
    camera cam;
    
    for (int j = ny - 1; j >= 0; j--) {
        for (int i = 0; i < nx; i++) {
            vec3 col(0, 0, 0);
            for (int s = 0; s < ns; s++) {
                float u = float(i + random_double()) / float(nx);
                float v = float(j + random_double()) / float(ny);
                ray r = cam.get_ray(u, v);
                col += color(r, world, 0);
            }
            col /= float(ns);
            col = vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
            
            int ir = int(255.99*col[0]);
            int ig = int(255.99*col[1]);
            int ib = int(255.99*col[2]);
            
            myfile << ir << " " << ig << " " << ib << "\n";
        }
    }
}
