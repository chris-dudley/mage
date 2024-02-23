#pragma once

#include <mgp.hpp>
#include <mg_utils.hpp>

#include <fmt/core.h>

#include <optional>
#include <string_view>

namespace shortest_paths {
namespace util {

std::optional<double> MapGetNumericOpt(const mgp::Map& map, std::string_view key) {
    if (!map.KeyExists(key)) {
        return std::nullopt;
    }

    auto value = map.At(key);
    if (value.IsNumeric()) {
        return value.ValueNumeric();
    }

    throw mgp::ValueException(fmt::format("Config parameter \"{}\" not a numeric value", key));
}

std::optional<int64_t> MapGetIntOpt(const mgp::Map& map, std::string_view key) {
    if (!map.KeyExists(key)) {
        return std::nullopt;
    }

    auto value = map.At(key);
    if (value.IsInt()) {
        return value.ValueInt();
    }

    throw mgp::ValueException(fmt::format("Config parameter \"{}\" not an integer value", key));
}

std::optional<std::string_view> MapGetStringViewOpt(const mgp::Map& map, std::string_view key) {
    if (!map.KeyExists(key)) {
        return std::nullopt;
    }

    auto value = map.At(key);
    if (value.IsString()) {
        return value.ValueString();
    }

    throw mgp::ValueException(fmt::format("Config parameter \"{}\" not a string", key));
}

std::optional<std::string> MapGetStringOpt(const mgp::Map& map, std::string_view key) {
    // Sadly, optional::transform and optional::and_then not available in c++20
    auto opt = MapGetStringViewOpt(map, key);
    if (opt.has_value()) {
        return std::string(opt.value());
    }
    return {};
}

std::optional<bool> MapGetBoolOpt(const mgp::Map& map, std::string_view key) {
    if (!map.KeyExists(key)) {
        return std::nullopt;
    }

    auto value = map.At(key);
    if (value.IsBool()) {
        return value.ValueBool();
    }

    throw mgp::ValueException(fmt::format("Config parameter \"{}\" not a bool", key));
}

std::optional<mgp::List> MapGetListOpt(const mgp::Map& map, std::string_view key) {
    if (!map.KeyExists(key)) {
        return std::nullopt;
    }

    auto value = map.At(key);
    if (value.IsList()) {
        return value.ValueList();
    }

    throw mgp::ValueException(fmt::format("Config parameter \"{}\" not a list", key));
}

} // namespace utils
} // namespace shortest_paths