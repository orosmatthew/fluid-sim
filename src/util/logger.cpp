#include "logger.hpp"


#include "spdlog/sinks/daily_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"

namespace util {

void init_logger()
{
    std::vector<spdlog::sink_ptr> sinks;
    sinks.push_back(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
    auto logger = std::make_shared<spdlog::logger>("main", std::begin(sinks), std::end(sinks));
    spdlog::register_logger(logger);
}

}
