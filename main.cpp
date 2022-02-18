#include <limits>
#include <iostream>
#include <fstream>
#include <vector>

#include "MathLibrary.h"

struct Material
{
    Material(const float &r, const Vec4f &a, const Vec3f &color, const float &spec) 
        : refractiveIndex(r), albedo(a), diffuseColor(color), specularExponent(spec)
    {
    }

    Material() 
        : refractiveIndex(1), albedo(1,0,0,0), diffuseColor(), specularExponent()
    {
    }

    float refractiveIndex;
    Vec4f albedo;
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

    // checking if ray coming out of camera hits this sphere
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
    return I-N*2.f*(I.dot(N));
}

Vec3f refract(const Vec3f &I, const Vec3f &N, const float &refractiveIndex)
{
    float cosi = -std::max(-1.0f, std::min(1.0f, I.dot(N)));
    float etai = 1;
    float etat = refractiveIndex;
    Vec3f n = N;
    // If the ray is inside the object, swap the indices and invert the normal to get the correct result
    if(cosi < 0)
    {
        cosi = -cosi;
        std::swap(etai, etat); 
        n = -N;
    }
    const float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? Vec3f(0,0,0) : I*eta + n*(eta * cosi - sqrtf(k));
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
constexpr int recursionDepth = 8;
Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, const std::vector<Light> &lights, size_t depth = 0)
{
    Vec3f point, N;
    Material material;
    // POint is the "t" in line formula. The exact point at which intersects the sphere. 
    if (depth > recursionDepth || !scene_intersect(orig, dir, spheres, point, N, material))
    {
        return backgroundColor;
    }

    Vec3f reflectDir = reflect(dir, N).getNorm();
    Vec3f refractDir = refract(dir, N, material.refractiveIndex).getNorm();

    Vec3f reflectOrig = reflectDir.dot(N) < 0 ? point - (N*1e-3) : point + N*1e-3; // slight offset because the point is directly on sphere
    Vec3f refractOrig = refractDir.dot(N) < 0 ? point - N*1e-3 : point + N*1e-3;

    Vec3f reflectColor = cast_ray(reflectOrig, reflectDir, spheres, lights, depth+1);
    Vec3f refractColor = cast_ray(refractOrig, refractDir, spheres, lights, depth+1);


    float diffuseLightIntensity = 0.0f;
    float specularLightIntensity = 0.0f;
    for (Light light : lights)
    {
        // vector from lightOrigin to the point the ray intersected with the sphere
        const Vec3f lightDirOrig = light.position - point;
        Vec3f lightDir = lightDirOrig.getNorm();
        float lightNormalDot = lightDir.dot(N);
        // Shadow detection
        float lightDistance = (light.position - point).mag();
        Vec3f shadowOrig = lightNormalDot < 0 ? point - N*1e-3 : point + N*1e-3; // checking if the point lies in the shadow of the light
        Vec3f shadowPt, shadowN;
        Material tmpMaterial;
        if(scene_intersect(shadowOrig, lightDir, spheres, shadowPt, shadowN, tmpMaterial) && (shadowPt-shadowOrig).mag() < lightDistance)
            continue;

        diffuseLightIntensity += light.intensity * std::max(0.f, lightNormalDot);
        specularLightIntensity += powf(std::max(0.f, -reflect(-lightDir, N).dot(dir)), material.specularExponent) * (light.intensity);
    }

    return material.diffuseColor * diffuseLightIntensity * material.albedo[0] +
        Vec3f(1.f, 1.f, 1.f) * specularLightIntensity * material.albedo[1] +
        reflectColor * material.albedo[2] +
        refractColor * material.albedo[3];
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
            Vec3f dir = Vec3f(x, y, -1.0f).getNorm(); // cast ray from (0,0,0) to the pixel (x,y,-1.0f), -1.0f as is out on the Zaxis
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
    Material      ivory(1.0f, Vec4f(0.6f, 0.3f, 0.1f, 0.0f), Vec3f(0.4f, 0.4f, 0.3f), 50.f);
    Material      glass(1.5f, Vec4f(0.0f, 0.5f, 0.1f, 0.8f), Vec3f(0.6f, 0.7f, 0.8f), 125.f);
    Material redRubber(1.0f, Vec4f(0.9f, 0.1f, 0.0f, 0.0f), Vec3f(0.3f, 0.1f, 0.1f), 10.f);
    Material     mirror(1.0f, Vec4f(0.0f, 10.0f, 0.8f, 0.0f), Vec3f(1.0f, 1.0f, 1.0f), 1425.f);

    std::vector<Sphere> spheres
    {
        Sphere(Vec3f(-3,0,-16),          2.0f, ivory),
        Sphere(Vec3f(-1,-1.5,-12),       2.0f, glass),
        Sphere(Vec3f(1.5f, -0.5f,-18.f), 3.f,  redRubber),
        Sphere(Vec3f(7,5,-18),           4.f,  mirror),
    };

    std::vector<Light> lights
    {
        Light(Vec3f(-20.f, 20.f, 20.f), 1.5f),
        //Light(Vec3f(30.f, 50.f, -25.f), 1.8f),
        //Light(Vec3f(30.f, 20.f, 30.f), 1.7f),
    };

    render(spheres, lights);
    return 0;
}