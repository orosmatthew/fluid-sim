#include <raylib-cpp.hpp>

#define LOGGER_RAYLIB
#include "util/logger.hpp"

namespace rl = raylib;

int main()
{
    util::init_logger();

#if NDEBUG
    LOG->set_level(spdlog::level::err);
#else
    LOG->set_level(spdlog::level::debug);
#endif

    SetTraceLogCallback(util::logger_callback_raylib);

    SetConfigFlags(FLAG_VSYNC_HINT | FLAG_MSAA_4X_HINT);

    const int screen_width = 800;
    const int screen_height = 600;

    rl::Window window(screen_width, screen_height, "Fluid Sim");

    while (!window.ShouldClose()) {
        BeginDrawing();

        ClearBackground(BLACK);
        rl::DrawText("Hello World!", 200, 200, 24, rl::Color::LightGray());
        DrawFPS(10, 10);

        EndDrawing();
    }

    return EXIT_SUCCESS;
}
