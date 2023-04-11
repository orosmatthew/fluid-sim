#pragma once

#include <raylib-cpp.hpp>

#include "fluid.hpp"
#include "util/fixed_loop.hpp"

namespace rl = raylib;

class App {
public:
    inline App()
        : m_window(1200, 1200, "Fluid Sim")
        , m_fixed_loop(60.0f)
        , m_fluid(512, 0.0f, 0.0f, 4)
        , m_prev_pos()
        , m_image(512, 512, rl::Color::Black())
        , m_texture(m_image)
        , m_scale(1200.0f / 512.0f)
    {
    }

    inline void update()
    {
        m_fixed_loop.update(5, [&]() {
            if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
                rl::Vector2 pos = GetMousePosition();
                pos /= m_scale;
                pos = pos.Clamp(
                    { 0.0f, 0.0f },
                    { static_cast<float>(m_window.GetWidth()) / m_scale - 1,
                      static_cast<float>(m_window.GetHeight()) / m_scale - 1 });
                m_prev_pos /= m_scale;
                m_prev_pos = m_prev_pos.Clamp(
                    { 0.0f, 0.0f },
                    { static_cast<float>(m_window.GetWidth()) / m_scale - 1,
                      static_cast<float>(m_window.GetHeight()) / m_scale - 1 });
                float amount_x = pos.x - m_prev_pos.x;
                float amount_y = pos.y - m_prev_pos.y;
                m_fluid.add_density(static_cast<int>(pos.x), static_cast<int>(pos.y), 10000.0f);
                m_fluid.add_velocity(static_cast<int>(pos.x), static_cast<int>(pos.y), amount_x, amount_y);
            }

            if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
                rl::Vector2 pos = GetMousePosition();
                pos /= m_scale;
                pos = pos.Clamp(
                    { 0.0f, 0.0f },
                    { static_cast<float>(m_window.GetWidth()) / m_scale - 1,
                      static_cast<float>(m_window.GetHeight()) / m_scale - 1 });
                m_fluid.add_density(static_cast<int>(pos.x), static_cast<int>(pos.y), -10000.0f);
            }

            if (IsMouseButtonDown(MOUSE_BUTTON_MIDDLE)) {
                rl::Vector2 pos = GetMousePosition();
                pos /= m_scale;
                pos = pos.Clamp(
                    { 0.0f, 0.0f },
                    { static_cast<float>(m_window.GetWidth()) / m_scale - 1,
                      static_cast<float>(m_window.GetHeight()) / m_scale - 1 });
                m_prev_pos /= m_scale;
                m_prev_pos = m_prev_pos.Clamp(
                    { 0.0f, 0.0f },
                    { static_cast<float>(m_window.GetWidth()) / m_scale - 1,
                      static_cast<float>(m_window.GetHeight()) / m_scale - 1 });
                float amount_x = pos.x - m_prev_pos.x;
                float amount_y = pos.y - m_prev_pos.y;
                m_fluid.add_velocity(static_cast<int>(pos.x), static_cast<int>(pos.y), amount_x, amount_y);
            }
            m_prev_pos = GetMousePosition();

            m_fluid.step(0.1f);
        });

        BeginDrawing();

        ClearBackground(BLACK);
        draw_fluid();
        DrawFPS(10, 10);

        EndDrawing();
    }

    [[nodiscard]] inline bool should_close() const
    {
        return m_window.ShouldClose();
    }

private:
    struct Vector2i {
        int x;
        int y;
    };

    rl::Window m_window;
    BS::thread_pool m_thread_pool {};
    util::FixedLoop m_fixed_loop;
    Fluid m_fluid;
    rl::Vector2 m_prev_pos;
    rl::Image m_image;
    rl::Texture m_texture;
    float m_scale;

    inline static Vector2i index_to_pos(size_t index, int size)
    {
        return { .x = static_cast<int>(index % size), .y = static_cast<int>(index / size) };
    }

    inline void draw_fluid()
    {
        m_thread_pool
            .parallelize_loop(
                m_fluid.size() * m_fluid.size(),
                [&](int begin, int end) {
                    for (int i = begin; i < end; i++) {
                        Vector2i pos = index_to_pos(i, m_fluid.size());
                        float d = m_fluid.density_at(pos.x, pos.y);
                        float c = map(d, 0.0f, 10000.0f, 0.0f, 255.0f);
                        c = std::clamp(c, 0.0f, 255.0f);
                        m_image.DrawPixel(pos.x, pos.y, rl::Color(0, static_cast<char>(c), static_cast<char>(c), 255));
                    }
                })
            .wait();
        m_texture.Update(m_image.data);
        m_texture.Draw(rl::Vector2(0, 0), 0.0f, m_scale);
    }

    inline static float map(float value, float in_min, float in_max, float out_min, float out_max)
    {
        return (value - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
    }
};