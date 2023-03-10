#if defined(PLATFORM_WEB)
#include <emscripten.h>
#endif

#include "app.hpp"

#if defined(PLATFORM_WEB)
void web_update(void* ptr)
{
    App& app = *static_cast<App*>(ptr);
    app.update();
}
#endif

int main()
{
    SetConfigFlags(FLAG_VSYNC_HINT | FLAG_MSAA_4X_HINT);
    App app;
#if defined(PLATFORM_WEB)
    emscripten_set_main_loop_arg(web_update, &app, 0, 1);
#else
    while (!app.should_close()) {
        app.update();
    }
#endif

    return EXIT_SUCCESS;
}
