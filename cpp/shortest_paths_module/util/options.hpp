#pragma once

#include <optional>
#include <string_view>

#include <mgp.hpp>

#include "util/map.hpp"

namespace shortest_paths {
namespace util {

class Options {
public:
    Options(): _map()
    {}
    Options(mgp::Map&& map): _map(std::move(map))
    {}
    Options(const mgp::Map& map): _map(map)
    {}

    std::optional<std::string_view> StringView(std::string_view key) const {
        return MapGetStringViewOpt(_map, key);
    }

    std::optional<std::string> String(std::string_view key) const {
        return MapGetStringOpt(_map, key);
    }

    std::optional<double> Numeric(std::string_view key) const {
        return MapGetNumericOpt(_map, key);
    }

    std::optional<mgp::List> List(std::string_view key) const {
        return MapGetListOpt(_map, key);
    }
private:
    mgp::Map _map;
};

} // namespace util
} // namespace shortest_path