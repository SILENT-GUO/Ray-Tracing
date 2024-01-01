#include <iostream>
#include <vector>
#include "Common.h"
#include "Scene.h"
#include "Camera.h"
#include "Material.h"
#include "Hittable.h"
#include "Utils/lodepng.h"

const int kMaxTraceDepth = 5;

Color TraceRay(const Ray& ray,
               const std::vector<LightSource>& light_sources,
               const Hittable& scene,
               int trace_depth);

Color Shade(const std::vector<LightSource>& light_sources,
            const Hittable& hittable_collection,
            const HitRecord& hit_record,
            int trace_depth) {
    // TODO: Add your code here.
    Color color(0.f, 0.f, 0.f);
    HitRecord record = hit_record;
    HitRecord dummy;

    //ambient
    color = hit_record.material.ambient * hit_record.material.k_a;

    //recurse through all light sources
    for (int i = 0; i < light_sources.size(); i++){
        //get light source
        LightSource light = light_sources[i];
        //get light direction
        Vec shadow_ray = glm::normalize(light.position - hit_record.position);
        if(glm::dot(hit_record.normal, shadow_ray) > 0){
            Ray temp_ray(hit_record.position, shadow_ray);
            if(!hittable_collection.Hit(temp_ray, &dummy)){
                Vec shadow_ray_reflection = glm::normalize(glm::reflect(-shadow_ray, hit_record.normal));
                color += 1.0f * light.intensity * (hit_record.material.diffuse * hit_record.material.k_d * glm::dot(hit_record.normal, shadow_ray)
                        + hit_record.material.specular * hit_record.material.k_s * glm::pow(glm::dot(shadow_ray_reflection, -glm::normalize(record.in_direction)), hit_record.material.sh));
            }
        }
    }
    //ray tracing for kMaxTraceDepth
    if(trace_depth < kMaxTraceDepth){
        if(hit_record.material.k_s > 0){
            Ray reflection_ray(hit_record.position, record.reflection);
            color += hit_record.material.k_s * TraceRay(reflection_ray, light_sources, hittable_collection, trace_depth + 1);
        }
    }
    //clamping color
    for (int i = 0; i < 3; i++){
        if(color[i] > 1.0f){
            color[i] = 1.0f;
        }
    }
    return color;
}

Color TraceRay(const Ray& ray,
               const std::vector<LightSource>& light_sources,
               const Hittable& hittable_collection,
               int trace_depth) {
    // TODO: Add your code here.
    HitRecord record;
    Color color(0.0f, 0.0f, 0.0f);
    if(hittable_collection.Hit(ray, &record)){
        color = Shade(light_sources, hittable_collection, record, trace_depth);
    }
    return color;
}

int main() {
    // TODO: Set your workdir (absolute path) here.
    const std::string work_dir("/Users/silent/Documents/Year4/Semester 1/COMP3271/PAs/PA2_Release/");

    // Construct scene
    Scene scene(work_dir, "scene/teapot.toml");
    const Camera& camera = scene.camera_;
    int width = camera.width_;
    int height = camera.height_;

    std::vector<unsigned char> image(width * height * 4, 0);

    float progress = 0.f;

    // Traverse all pixels
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            Color color(0.f, 0.f, 0.f);
            int count = 0;
            for (float bias_x = 0.25f; bias_x < 1.f; bias_x += .5f) {
                for (float bias_y = 0.25f; bias_y < 1.f; bias_y += .5f) {
                    Ray ray = camera.RayAt(float(x) + bias_x, float(y) + bias_y);
                    color += TraceRay(ray, scene.light_sources_, scene.hittable_collection_, 1);
                    count++;
                }
            }
            color /= float(count);
            int idx = 4 * ((height - y - 1) * width + x);
            for (int i = 0; i < 3; i++) {
                image[idx + i] = (uint8_t) (glm::min(color[i], 1.f - 1e-5f) * 256.f);
            }
            image[idx + 3] = 255;

            float curr_progress = float(x * height + y) / float(height * width);
            if (curr_progress > progress + 0.05f) {
                progress += 0.05f;
                std::cout << "Progress: " << progress << std::endl;
            }
        }
    }

    // Save result as png file
    std::vector<unsigned char> png;
    unsigned error = lodepng::encode(png, image, width, height);
    lodepng::save_file(png, work_dir + "outputs/output.png");
}
