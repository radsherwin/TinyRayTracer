#include <limits>
#include <iostream>
#include <fstream>
#include <vector>

#include "MathLibrary.h"

struct Material
{
    Material(const Vec2f &a, const Vec3f &color, const float &spec) 
        : albedo(a), diffuseColor(color), specularExponent(spec)
    {
    }

    Material() 
        : albedo(1,0), diffuseColor(), specularExponent()
    {
    }

    Vec2f albedo;
    Vec3f diffuseColor;
    float specularExponent;
};


struct Light
{
    Light(const Vec3f &p, const float &i) : position(p), intensity(i)
    {
    }

    Vec3f position;
    float intensity;
};

struct Sphere
{
    Vec3f center;
    float radius;
    Material material;

    Sphere(const Vec3f &c, const float &r, const Material &m) : center(c), radius(r), material(m)
    {
    }

    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &t0) const
    {
        Vec3f L = center - orig;
        float tca = L.dot(dir);
        float d2 = L.dot(L) - tca * tca;
        if (d2 > radius * radius) return false;
        float thc = sqrtf(radius * radius - d2);
        t0 = tca - thc;
        float t1 = tca + thc;
        if (t0 > t1) std::swap(t0, t1);
        if (t0 < 0)
        {
            t0 = t1;
            return false;
        }
        return true;
    }
};

Vec3f reflect(const Vec3f &I, const Vec3f &N)
{
    return I-N*2.f*(I*N);
}

bool scene_intersect(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, Vec3f &hit, Vec3f &N, Material &material)
{
    float spheres_dist = std::numeric_limits<float>::max();
    for (Sphere sphere : spheres)
    {
        float dist_i;
        if (sphere.ray_intersect(orig, dir, dist_i) && dist_i < spheres_dist)
        {
            spheres_dist = dist_i;
            hit = orig + dir * dist_i;
            N = (hit - sphere.center).getNorm();
            material = sphere.material;
        }
    }

    return spheres_dist < 1000;
}
// the dot product of Vectors of unit length give the intensity of surface illumination
Vec3f backgroundColor{ 0.2f, 0.7f, 0.8f };

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, const std::vector<Light> &lights)
{
    Vec3f point, N;
    Material material;

    if (!scene_intersect(orig, dir, spheres, point, N, material))
    {
        return backgroundColor;
    }

    float diffuseLightIntensity = 0.0f;
    float specularLightIntensity = 0.0f;
    for (Light light : lights)
    {
        Vec3f lightDir = (light.position - point).getNorm();
        diffuseLightIntensity += light.intensity * std::max(0.f, lightDir.dot(N));
        specularLightIntensity += powf(std::max(0.f, -reflect(-lightDir, N).dot(dir)), material.specularExponent) * (light.intensity);
    }

    return material.diffuseColor * diffuseLightIntensity * material.albedo[0] + Vec3f(1.f, 1.f, 1.f) * specularLightIntensity * material.albedo[1];
}

// Draw a "ray" from center of camera, through exch pixel, and check if that ray intersects with the sphere
void render(const std::vector<Sphere> &spheres, const std::vector<Light> &lights)
{
    const int width = 1024;
    const int height = 768;
    const int fov = 3.14159 / 2;
    std::vector<Vec3f> framebuffer(width * height);

    #pragma omp parallel for
    for (size_t j = 0; j < height; j++)
    {
        for (size_t i = 0; i < width; i++)
        {
            float x = (2 * (i + 0.5f) / (float)width - 1.0f) * tanf(fov / 2.0f) * width / (float)height;
            float y = -(2 * (j + 0.5f) / (float)width - 1.0f) * tanf(fov / 2.0f);
            Vec3f dir = Vec3f(x, y, -1.0f).getNorm();
            framebuffer[i + j * width] = cast_ray(Vec3f(0, 0, 0), dir, spheres, lights);
        }
    }

    std::ofstream ofs;
    ofs.open("./out.ppm", std::ofstream::out | std::ofstream::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i< height * width; ++i)
    {
        Vec3f &c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if(max > 1) c = c*(1.f/max);
        for (size_t j = 0; j < 3; j++)
        {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
}

int main()
{
    Material ivory(Vec2f(0.6f, 0.3f), Vec3f(0.4f, 0.4f, 0.3f), 50.f);
    Material redRubber(Vec2f(0.9f, 0.1f), Vec3f(0.3f, 0.1f, 0.1f), 10.f);

    std::vector<Sphere> spheres
    {
        Sphere(Vec3f(-3,0,-16),          2.0f, ivory),
        Sphere(Vec3f(-1,-1.5,-12),       2.0f, redRubber),
        Sphere(Vec3f(1.5f, -0.5f,-18.f), 3.f,  redRubber),
        Sphere(Vec3f(7,5,-18),           4.f,  ivory),
    };

    std::vector<Light> lights
    {
        Light(Vec3f(-20.f, 20.f, 20.f), 1.5f),
        Light(Vec3f(30.f, 50.f, -25.f), 1.8f),
        Light(Vec3f(30.f, 20.f, 30.f), 1.7f),
    };

    render(spheres, lights);
    return 0;
}