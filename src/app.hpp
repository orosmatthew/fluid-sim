#pragma once

#include <raylib-cpp.hpp>

#ifdef ENABLE_MULTITHREADING
#define FLUID_MULTITHREADING
#endif

#include "fluid.hpp"
#include "util/fixed_loop.hpp"

namespace rl = raylib;

constexpr int g_sim_size = 256;

class App {
public:
    App()
        : m_window(1200, 1200, "Fluid Sim")
        , m_fixed_loop(60.0f)
        , m_fluid(g_sim_size, 0.0f, 0.0f, 4)
        , m_image(g_sim_size, g_sim_size, rl::Color::Black())
        , m_texture(m_image)
        , m_scale(1200.0f / static_cast<float>(g_sim_size))
    {
    }

    void update()
    {
        m_fixed_loop.update(5, [&] {
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
                const float amount_x = pos.x - m_prev_pos.x;
                const float amount_y = pos.y - m_prev_pos.y;
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
                const float amount_x = pos.x - m_prev_pos.x;
                const float amount_y = pos.y - m_prev_pos.y;
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

    [[nodiscard]] bool should_close() const
    {
        return m_window.ShouldClose();
    }

private:
    struct Vector2i {
        int x;
        int y;
    };

    rl::Window m_window;
#ifdef ENABLE_MULTITHREADING
    BS::thread_pool m_thread_pool {};
#endif
    util::FixedLoop m_fixed_loop;
    Fluid m_fluid;
    rl::Vector2 m_prev_pos;
    rl::Image m_image;
    rl::Texture m_texture;
    float m_scale;

    static Vector2i index_to_pos(const size_t index, const int size)
    {
        return { .x = static_cast<int>(index % size), .y = static_cast<int>(index / size) };
    }

    void draw_fluid()
    {
#ifdef ENABLE_MULTITHREADING
        m_thread_pool
            .parallelize_loop(
                m_fluid.size() * m_fluid.size(),
                [&](const int begin, const int end) {
                    for (int i = begin; i < end; i++) {
                        auto [x, y] = index_to_pos(i, m_fluid.size());
                        const float d = m_fluid.density_at(x, y);
                        float c = map(d, 0.0f, 10000.0f, 0.0f, 255.0f);
                        c = std::clamp(c, 0.0f, 255.0f);
                        m_image.DrawPixel(x, y, rl::Color(0, static_cast<char>(c), static_cast<char>(c), 255));
                    }
                })
            .wait();
#else
        for (int i = 0; i < m_fluid.size() * m_fluid.size(); i++) {
            Vector2i pos = index_to_pos(i, m_fluid.size());
            float d = m_fluid.density_at(pos.x, pos.y);
            float c = map(d, 0.0f, 10000.0f, 0.0f, 255.0f);
            c = std::clamp(c, 0.0f, 255.0f);
            m_image.DrawPixel(pos.x, pos.y, rl::Color(0, static_cast<char>(c), static_cast<char>(c), 255));
        }
#endif
        m_texture.Update(m_image.data);
        m_texture.Draw(rl::Vector2(0, 0), 0.0f, m_scale);
    }

    static float map(
        const float value, const float in_min, const float in_max, const float out_min, const float out_max)
    {
        return (value - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
    }
};