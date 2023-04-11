#include <chrono>
#include <iostream>

#include <FastNoiseLite.h>

#include "mve/detail/defs.hpp"
#include "mve/math/math.hpp"
#include "mve/renderer.hpp"
#include "mve/shader.hpp"
#include "mve/vertex_data.hpp"
#include "mve/window.hpp"

#include "camera.hpp"
#include "fluid.hpp"
#include "logger.hpp"
#include "text_buffer.hpp"
#include "text_pipeline.hpp"
#include "util/fixed_loop.hpp"

int main()
{
    initLogger();

    mve::Window window("Fluid Sim 3D", mve::Vector2i(600, 600), true);
    mve::Renderer renderer(window, "Fluid Sim 3D", 1, 0, 0);

    mve::VertexLayout vertex_layout;
    vertex_layout.push_back(mve::VertexAttributeType::vec3); // pos

    mve::VertexData data(vertex_layout);
    data.push_back(mve::Vector3(-0.5f, -0.5f, 0.5f));
    data.push_back(mve::Vector3(0.5f, -0.5f, 0.5f));
    data.push_back(mve::Vector3(0.5f, -0.5f, -0.5f));
    data.push_back(mve::Vector3(-0.5f, -0.5f, -0.5f));

    data.push_back(mve::Vector3(0.5f, 0.5f, 0.5f));
    data.push_back(mve::Vector3(-0.5f, 0.5f, 0.5f));
    data.push_back(mve::Vector3(-0.5f, 0.5f, -0.5f));
    data.push_back(mve::Vector3(0.5f, 0.5f, -0.5f));

    data.push_back(mve::Vector3(-0.5f, 0.5f, 0.5f));
    data.push_back(mve::Vector3(-0.5f, -0.5f, 0.5f));
    data.push_back(mve::Vector3(-0.5f, -0.5f, -0.5f));
    data.push_back(mve::Vector3(-0.5f, 0.5f, -0.5f));

    data.push_back(mve::Vector3(0.5f, -0.5f, 0.5f));
    data.push_back(mve::Vector3(0.5f, 0.5f, 0.5f));
    data.push_back(mve::Vector3(0.5f, 0.5f, -0.5f));
    data.push_back(mve::Vector3(0.5f, -0.5f, -0.5f));

    data.push_back(mve::Vector3(-0.5f, 0.5f, 0.5f));
    data.push_back(mve::Vector3(0.5f, 0.5f, 0.5f));
    data.push_back(mve::Vector3(0.5f, -0.5f, 0.5f));
    data.push_back(mve::Vector3(-0.5f, -0.5f, 0.5f));

    data.push_back(mve::Vector3(0.5f, 0.5f, -0.5f));
    data.push_back(mve::Vector3(-0.5f, 0.5f, -0.5f));
    data.push_back(mve::Vector3(-0.5f, -0.5f, -0.5f));
    data.push_back(mve::Vector3(0.5f, -0.5f, -0.5f));

    std::array<uint32_t, 6> base_indices = { 0, 3, 2, 0, 2, 1 };
    std::vector<uint32_t> indices;
    indices.reserve(6 * 6);
    for (int f = 0; f < 6; f++) {
        for (int i = 0; i < 6; i++) {
            indices.push_back(base_indices[i] + (f * 4));
        }
    }

    mve::VertexBuffer vertex_buffer = renderer.create_vertex_buffer(data);
    mve::IndexBuffer index_buffer = renderer.create_index_buffer(indices);

    mve::Shader vert_shader("../res/bin/shader/cloud.vert.spv");
    mve::Shader frag_shader("../res/bin/shader/cloud.frag.spv");
    mve::GraphicsPipeline pipeline = renderer.create_graphics_pipeline(vert_shader, frag_shader, vertex_layout, true);

    mve::DescriptorSet global_descriptor = pipeline.create_descriptor_set(vert_shader.descriptor_set(0));
    mve::UniformBuffer global_ubo = renderer.create_uniform_buffer(vert_shader.descriptor_set(0).binding(0));
    global_descriptor.write_binding(vert_shader.descriptor_set(0).binding(0), global_ubo);
    mve::Matrix4 proj = mve::perspective(90.0f, 1.0f, 0.001f, 100.0f);
    global_ubo.update(vert_shader.descriptor_set(0).binding(0).member("proj").location(), proj);

    Camera camera;

    // Clouds Volume
    const int width = 64;
    const int height = 64;
    const int depth = 64;
    const float scale = 4.0f;

    const int voxel_per_slice = width * height;
    const int voxel_count = voxel_per_slice * depth;

    std::vector<std::byte> buffer(voxel_count);

    FastNoiseLite noise;
    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2S);

    for (int i = 0; i < voxel_count; ++i) {
        const int x = i % width;
        const auto y = (i % voxel_per_slice) / width;
        const auto z = i / voxel_per_slice;
        const float p = noise.GetNoise((float)x * scale, (float)y * scale, (float)z * scale);
        const float rand = (p + 1.0f) * 0.5f;
        buffer[i] = (std::byte)mve::clamp((int)mve::round(rand * 16), 0, 16);
    }

    mve::Texture texture = renderer.create_texture(mve::TextureFormat::r, width, height, depth, buffer.data());

    global_descriptor.write_binding(frag_shader.descriptor_set(0).binding(1), texture);

    util::FixedLoop fixed_loop(60.0f);

    TextPipeline text_pipeline(renderer, 36);
    TextBuffer fps_text = text_pipeline.create_text_buffer("FPS:", { 0, 0 }, 1.0f, { 1.0f, 1.0f, 1.0f });

    window.set_resize_callback([&](mve::Vector2i new_size) {
        renderer.resize(window);
        text_pipeline.resize();
        mve::Matrix4 proj = mve::perspective(90.0f, (float)new_size.x / (float)new_size.y, 0.001f, 100.0f);
        global_ubo.update(vert_shader.descriptor_set(0).binding(0).member("proj").location(), proj);
    });
    text_pipeline.resize();

    window.disable_cursor();
    bool cursor_captured = true;
    int current_frame_count = 0;
    int frame_count = 0;
    auto begin_time = std::chrono::steady_clock::time_point();
    while (!window.should_close()) {
        window.poll_events();

        if (window.is_key_pressed(mve::Key::c)) {
            if (cursor_captured) {
                window.enable_cursor();
                cursor_captured = false;
            }
            else {
                window.disable_cursor();
                cursor_captured = true;
            }
        }

        if (window.is_key_pressed(mve::Key::f)) {
            if (window.is_fullscreen()) {
                window.windowed();
            }
            else {
                window.fullscreen(true);
            }
        }

        camera.update(window);
        fixed_loop.update(20, [&] { camera.fixed_update(window); });

        global_ubo.update(
            vert_shader.descriptor_set(0).binding(0).member("view").location(), camera.view_matrix(fixed_loop.blend()));

        renderer.begin_frame(window);
        renderer.begin_render_pass_present();

        renderer.bind_graphics_pipeline(pipeline);
        renderer.bind_descriptor_set(global_descriptor);
        renderer.bind_vertex_buffer(vertex_buffer);
        renderer.draw_index_buffer(index_buffer);

        text_pipeline.draw(fps_text);

        renderer.end_render_pass_present();
        renderer.end_frame(window);

        std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();

        if (std::chrono::duration_cast<std::chrono::microseconds>(end_time - begin_time).count() >= 1000000) {
            begin_time = std::chrono::high_resolution_clock::now();
            frame_count = current_frame_count;
            fps_text.update("FPS:" + std::to_string(frame_count));
            current_frame_count = 0;
        }
        current_frame_count++;
    }
}