#pragma once

#include <raylib-cpp.hpp>

#include "fluid.hpp"
#include "util/fixed_loop.hpp"

namespace rl = raylib;

class App {
public:
    inline App()
        : m_window(640, 640, "Fluid Sim")
        , m_fixed_loop(60.0f)
        , m_fluid(128, 0.0f, 0.0f, 4)
        , m_prev_pos()
        , m_image(128, 128, rl::Color::Black())
        , m_texture(m_image)
    {
    }

    inline void update()
    {
        m_fixed_loop.update(5, [&]() {
            if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
                rl::Vector2 pos = GetMousePosition();
                pos /= 5.0f;
                pos = pos.Clamp(
                    { 0.0f, 0.0f },
                    { static_cast<float>(m_window.GetWidth()) / 5.0f - 1,
                      static_cast<float>(m_window.GetHeight()) / 5.0f - 1 });
                m_prev_pos /= 5.0f;
                m_prev_pos = m_prev_pos.Clamp(
                    { 0.0f, 0.0f },
                    { static_cast<float>(m_window.GetWidth()) / 5.0f - 1,
                      static_cast<float>(m_window.GetHeight()) / 5.0f - 1 });
                float amount_x = pos.x - m_prev_pos.x;
                float amount_y = pos.y - m_prev_pos.y;
                m_fluid.add_density(static_cast<int>(pos.x), static_cast<int>(pos.y), 10000.0f);
                m_fluid.add_velocity(static_cast<int>(pos.x), static_cast<int>(pos.y), amount_x, amount_y);
            }

            if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
                rl::Vector2 pos = GetMousePosition();
                pos /= 5.0f;
                pos = pos.Clamp(
                    { 0.0f, 0.0f },
                    { static_cast<float>(m_window.GetWidth()) / 5.0f - 1,
                      static_cast<float>(m_window.GetHeight()) / 5.0f - 1 });
                m_fluid.add_density(static_cast<int>(pos.x), static_cast<int>(pos.y), -10000.0f);
            }

            if (IsMouseButtonDown(MOUSE_BUTTON_MIDDLE)) {
                rl::Vector2 pos = GetMousePosition();
                pos /= 5.0f;
                pos = pos.Clamp(
                    { 0.0f, 0.0f },
                    { static_cast<float>(m_window.GetWidth()) / 5.0f - 1,
                      static_cast<float>(m_window.GetHeight()) / 5.0f - 1 });
                m_prev_pos /= 5.0f;
                m_prev_pos = m_prev_pos.Clamp(
                    { 0.0f, 0.0f },
                    { static_cast<float>(m_window.GetWidth()) / 5.0f - 1,
                      static_cast<float>(m_window.GetHeight()) / 5.0f - 1 });
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
    rl::Window m_window;
    util::FixedLoop m_fixed_loop;
    Fluid m_fluid;
    rl::Vector2 m_prev_pos;
    rl::Image m_image;
    rl::Texture m_texture;

    inline void draw_fluid()
    {
        for (int i = 0; i < m_fluid.size(); i++) {
            for (int j = 0; j < m_fluid.size(); j++) {
                float d = m_fluid.density_at(i, j);
                float c = map(d, 0.0f, 10000.0f, 0.0f, 255.0f);
                c = std::clamp(c, 0.0f, 255.0f);
                m_image.DrawPixel(
                    i, j, rl::Color(static_cast<char>(c), static_cast<char>(c), static_cast<char>(c), 255));
            }
        }
        m_texture.Update(m_image.data);
        m_texture.Draw(rl::Vector2(0, 0), 0.0f, 5.0f);
    }

    inline static float map(float value, float in_min, float in_max, float out_min, float out_max)
    {
        return (value - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
    }
};