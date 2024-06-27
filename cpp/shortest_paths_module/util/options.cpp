#include "options.hpp"
#include "map.hpp"

#include <optional>
#include <stdexcept>
#include <string_view>
#include <vector>

#include <fmt/core.h>

namespace shortest_paths {
namespace util {

Options::Options() : _map() {}
Options::Options(mgp::Map &&map) : _map(std::move(map)) {}
Options::Options(const mgp::Map &map) : _map(map) {}

std::optional<std::string_view> Options::StringView(std::string_view key) const {
  return MapGetStringViewOpt(_map, key);
}

std::optional<std::string> Options::String(std::string_view key) const { return MapGetStringOpt(_map, key); }

std::optional<double> Options::Numeric(std::string_view key) const { return MapGetNumericOpt(_map, key); }

std::optional<int64_t> Options::Integer(std::string_view key) const { return MapGetIntOpt(_map, key); }

std::optional<mgp::List> Options::List(std::string_view key) const { return MapGetListOpt(_map, key); }

std::optional<bool> Options::Bool(std::string_view key) const { return MapGetBoolOpt(_map, key); }

std::optional<std::vector<double>> Options::NumericList(std::string_view key) const {
  auto list = List(key);
  if (not list.has_value()) {
    return std::nullopt;
  }

  std::vector<double> result;
  for (const auto &value : *list) {
    if (not value.IsNumeric()) {
      throw std::invalid_argument(fmt::format("option {} index {} is not numeric", key, result.size()));
    }

    result.push_back(value.ValueNumeric());
  }

  return result;
}

}  // namespace util
}  // namespace shortest_paths
