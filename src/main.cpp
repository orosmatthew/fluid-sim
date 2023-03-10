#include <cassert>
#include <iostream>
#include <vector>

#include <raylib-cpp.hpp>
#if defined(PLATFORM_WEB)
#include <emscripten.h>
#endif

#include "util/fixed_loop.hpp"

namespace rl = raylib;

const int N = 64;
const int iter = 16;

struct FluidCube {
    int size {};
    float dt {};
    float diff {};
    float visc {};

    std::vector<float> s;
    std::vector<float> density;

    std::vector<float> vx;
    std::vector<float> vy;

    std::vector<float> vx0;
    std::vector<float> vy0;
};

static size_t IX(int x, int y)
{
    return static_cast<size_t>(x) + static_cast<size_t>(y) * N;
}

static FluidCube fluid_cube_create(int size, int diffusion, int viscosity, float dt)
{
    assert(size > 0);
    FluidCube cube;
    cube.size = size;
    cube.dt = dt;
    cube.diff = static_cast<float>(diffusion);
    cube.visc = static_cast<float>(viscosity);

    std::vector<float> vec(size * size, 0.0f);
    cube.s = vec;
    cube.density = vec;

    cube.vx = vec;
    cube.vy = vec;

    cube.vx0 = vec;
    cube.vy0 = vec;

    return cube;
}

static void fluid_cube_add_density(FluidCube& cube, int x, int y, float amount)
{
    cube.density[IX(x, y)] += amount;
}

static void fluid_cube_add_velocity(FluidCube& cube, int x, int y, float amount_x, float amount_y)
{
    size_t index = IX(x, y);
    cube.vx[index] += amount_x;
    cube.vy[index] += amount_y;
}

static void set_bnd(int b, std::vector<float>& x)
{
    for (int i = 1; i < N - 1; i++) {
        x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N - 1)] = b == 2 ? -x[IX(i, N - 2)] : x[IX(i, N - 2)];
    }
    for (int j = 1; j < N - 1; j++) {
        x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
        x[IX(N - 1, j)] = b == 1 ? -x[IX(N - 2, j)] : x[IX(N - 2, j)];
    }

    x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N - 1)] = 0.5f * (x[IX(1, N - 1)] + x[IX(0, N - 2)]);
    x[IX(N - 1, 0)] = 0.5f * (x[IX(N - 2, 0)] + x[IX(N - 1, 1)]);
    x[IX(N - 1, N - 1)] = 0.5f * (x[IX(N - 2, N - 1)] + x[IX(N - 1, N - 2)]);
}

static void lin_solve(int b, std::vector<float>& x, const std::vector<float>& x0, float a, float c)
{
    const float c_inv = 1.0f / c;
    for (int t = 0; t < iter; t++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                x[IX(i, j)]
                    = (x0[IX(i, j)] + a * (x[IX(i + 1, j)] + x[IX(i - 1, j)] + x[IX(i, j + 1)] + x[IX(i, j - 1)]))
                    * c_inv;
            }
        }
        set_bnd(b, x);
    }
}

static void project(
    std::vector<float>& vel_x, std::vector<float>& vel_y, std::vector<float>& p, std::vector<float>& div)
{
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            div[IX(i, j)]
                = (-0.5f * (vel_x[IX(i + 1, j)] - vel_x[IX(i - 1, j)] + vel_y[IX(i, j + 1)] - vel_y[IX(i, j - 1)])) / N;
            p[IX(i, j)] = 0;
        }
    }
    set_bnd(0, div);
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 6);

    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            vel_x[IX(i, j)] -= 0.5f * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) * N;
            vel_y[IX(i, j)] -= 0.5f * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) * N;
        }
    }
    set_bnd(1, vel_x);
    set_bnd(2, vel_y);
}

static void advect(
    int b,
    std::vector<float>& d,
    const std::vector<float>& d0,
    const std::vector<float>& vel_x,
    const std::vector<float>& vel_y,
    float dt)
{
    float i0, i1, j0, j1;

    float dtx = dt * (N - 2);
    float dty = dtx;

    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;

    float Nfloat = N - 2;
    float ifloat, jfloat;
    int i, j;

    for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
        for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
            tmp1 = dtx * vel_x[IX(i, j)];
            tmp2 = dty * vel_y[IX(i, j)];
            x = ifloat - tmp1;
            y = jfloat - tmp2;

            if (x < 0.5f)
                x = 0.5f;
            if (x > Nfloat + 0.5f)
                x = Nfloat + 0.5f;
            i0 = floorf(x);
            i1 = i0 + 1.0f;
            if (y < 0.5f)
                y = 0.5f;
            if (y > Nfloat + 0.5f)
                y = Nfloat + 0.5f;
            j0 = floorf(y);
            j1 = j0 + 1.0f;

            s1 = x - i0;
            s0 = 1.0f - s1;
            t1 = y - j0;
            t0 = 1.0f - t1;

            int i0i = static_cast<int>(i0);
            int i1i = static_cast<int>(i1);
            int j0i = static_cast<int>(j0);
            int j1i = static_cast<int>(j1);

            d[IX(i, j)] = s0 * (t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)])
                + s1 * (t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)]);
        }
    }

    set_bnd(b, d);
}

static void diffuse(int b, std::vector<float>& x, const std::vector<float>& x0, float diff, float dt)
{
    float a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a);
}

static void fluid_cube_step(FluidCube& cube)
{
    //    auto visc = cube.visc;
    //    auto diff = cube.diff;
    //    auto dt = cube.dt;
    //    auto Vx = cube.vx;
    //    auto Vy = cube.vy;
    //    auto Vx0 = cube.vx0;
    //    auto Vy0 = cube.vy0;
    //    auto s = cube.s;
    //    auto density = cube.density;

    for (size_t i = 0; i < N; i++) {
        cube.density[i] = std::clamp(cube.density[i], 0.0f, 10.0f);
    }

    diffuse(1, cube.vx0, cube.vx, cube.visc, cube.dt);
    diffuse(2, cube.vy0, cube.vy, cube.visc, cube.dt);

    project(cube.vx0, cube.vy0, cube.vx, cube.vy);

    advect(1, cube.vx, cube.vx0, cube.vx0, cube.vy0, cube.dt);
    advect(2, cube.vy, cube.vy0, cube.vx0, cube.vy0, cube.dt);

    project(cube.vx, cube.vy, cube.vx0, cube.vy0);
    diffuse(0, cube.s, cube.density, cube.diff, cube.dt);
    advect(0, cube.density, cube.s, cube.vx, cube.vy, cube.dt);

    for (size_t i = 0; i < N; i++) {
        cube.density[i] = std::clamp(cube.density[i], 0.0f, 10.0f);
    }
}

void draw_fluid(const FluidCube& cube)
{
    for (int i = 0; i < cube.size; i++) {
        for (int j = 0; j < cube.size; j++) {
            float x = static_cast<float>(i) * 10.0f;
            float y = static_cast<float>(j) * 10.0f;
            float d = cube.density[IX(i, j)];
            d = std::clamp(d, 0.0f, 255.0f);
            DrawRectangle(
                static_cast<int>(x),
                static_cast<int>(y),
                10.0f,
                10.0f,
                rl::Color(static_cast<char>(d), static_cast<char>(d), static_cast<char>(d), 255));
        }
    }
}

struct State {
    util::FixedLoop fixed_loop;
    FluidCube cube;
    rl::Vector2 prev_pos;
};

const int screen_width = 640;
const int screen_height = 640;

void main_loop(void* state_ptr)
{
    State& state = *static_cast<State*>(state_ptr);
    if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
        rl::Vector2 pos = GetMousePosition();
        pos /= 10.0f;
        pos = pos.Clamp(
            { 0.0f, 0.0f },
            { static_cast<float>(screen_width) / 10.0f - 1, static_cast<float>(screen_height) / 10.0f - 1 });
        state.prev_pos /= 10.0f;
        state.prev_pos = state.prev_pos.Clamp(
            { 0.0f, 0.0f },
            { static_cast<float>(screen_width) / 10.0f - 1, static_cast<float>(screen_height) / 10.0f - 1 });
        float amount_x = pos.x - state.prev_pos.x;
        float amount_y = pos.y - state.prev_pos.y;
        fluid_cube_add_density(state.cube, static_cast<int>(pos.x), static_cast<int>(pos.y), 100.0f);
        fluid_cube_add_velocity(state.cube, static_cast<int>(pos.x), static_cast<int>(pos.y), amount_x, amount_y);
    }

    if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
        rl::Vector2 pos = GetMousePosition();
        pos /= 10.0f;
        pos = pos.Clamp(
            { 0.0f, 0.0f },
            { static_cast<float>(screen_width) / 10.0f - 1, static_cast<float>(screen_height) / 10.0f - 1 });
        fluid_cube_add_density(state.cube, static_cast<int>(pos.x), static_cast<int>(pos.y), -100.0f);
    }
    state.prev_pos = GetMousePosition();

    state.fixed_loop.update(5, [&]() { fluid_cube_step(state.cube); });

    BeginDrawing();

    ClearBackground(BLACK);
    draw_fluid(state.cube);
    DrawFPS(10, 10);

    EndDrawing();
}

int main()
{
    SetConfigFlags(FLAG_VSYNC_HINT | FLAG_MSAA_4X_HINT);

    rl::Window window(screen_width, screen_height, "Fluid Sim");

    State state { .fixed_loop = util::FixedLoop(60.0f), .cube = fluid_cube_create(64, 0, 0, 0.1f), .prev_pos {} };

#if defined(PLATFORM_WEB)
    emscripten_set_main_loop_arg(main_loop, &state, 0, 1);
#else
    while (!window.ShouldClose()) {
        main_loop(&state);
    }
#endif

    return EXIT_SUCCESS;
}
