#pragma once

#include "mve/detail/defs.hpp"
#include "mve/renderer.hpp"

class SimplePipeline {
public:
    SimplePipeline(mve::Renderer& renderer);

    void set_view(const mve::Matrix4 mat);

    void resize(mve::Vector2i extent);

    inline const mve::DescriptorSet& global_descriptor_set() const
    {
        return m_global_descriptor_set;
    }

    inline static mve::VertexLayout vertex_layout()
    {
        return {
            mve::VertexAttributeType::vec3, // Position
            mve::VertexAttributeType::vec3, // Color
        };
    }

    inline mve::GraphicsPipeline& pipeline()
    {
        return m_pipeline;
    }

    inline const mve::GraphicsPipeline& pipeline() const
    {
        return m_pipeline;
    }

    mve::ShaderDescriptorSet model_descriptor_set() const;

    mve::ShaderDescriptorBinding model_uniform_binding() const;

private:
    mve::Renderer* m_renderer;
    mve::Shader m_vert_shader;
    mve::Shader m_frag_shader;
    mve::GraphicsPipeline m_pipeline;
    mve::DescriptorSet m_global_descriptor_set;
    mve::UniformBuffer m_global_ubo;
    mve::UniformLocation m_view_loc;
    mve::UniformLocation m_proj_loc;
};