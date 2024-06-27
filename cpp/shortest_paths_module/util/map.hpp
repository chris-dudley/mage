#pragma once

#include <mgp.hpp>

#include <optional>
#include <string_view>

namespace shortest_paths {
namespace util {

std::optional<double> MapGetNumericOpt(const mgp::Map &map, std::string_view key);

std::optional<int64_t> MapGetIntOpt(const mgp::Map &map, std::string_view key);

std::optional<std::string_view> MapGetStringViewOpt(const mgp::Map &map, std::string_view key);

std::optional<std::string> MapGetStringOpt(const mgp::Map &map, std::string_view key);

std::optional<bool> MapGetBoolOpt(const mgp::Map &map, std::string_view key);

std::optional<mgp::List> MapGetListOpt(const mgp::Map &map, std::string_view key);

}  // namespace util
}  // namespace shortest_paths
