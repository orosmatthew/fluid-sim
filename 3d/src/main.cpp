#include <iostream>

#include "mve/detail/defs.hpp"
#include "mve/renderer.hpp"
#include "mve/shader.hpp"
#include "mve/vertex_data.hpp"
#include "mve/window.hpp"

#include "logger.hpp"

int main()
{
    initLogger();

    mve::Window window("Fluid Sim 3D", mve::Vector2i(800, 600), false);
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
    global_ubo.update(vert_shader.descriptor_set(0).binding(0).member("proj").location(), mve::Matrix4::identity());
    global_ubo.update(
        vert_shader.descriptor_set(0).binding(0).member("view").location(),
        mve::Matrix4::identity().translate({ 0, 0, 1 }));

    while (!window.should_close()) {
        window.poll_events();

        renderer.begin_frame(window);
        renderer.begin_render_pass_present();

        renderer.bind_graphics_pipeline(pipeline);
        renderer.bind_descriptor_set(global_descriptor);
        renderer.bind_vertex_buffer(vertex_buffer);
        renderer.draw_index_buffer(index_buffer);

        renderer.end_render_pass_present();
        renderer.end_frame(window);
    }
}