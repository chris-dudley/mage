#pragma once

#include <optional>
#include <string_view>
#include <vector>

#include <mgp.hpp>

namespace shortest_paths {
namespace util {

class Options {
 public:
  Options();
  Options(mgp::Map &&map);
  Options(const mgp::Map &map);

  std::optional<std::string_view> StringView(std::string_view key) const;
  std::optional<std::string> String(std::string_view key) const;
  std::optional<double> Numeric(std::string_view key) const;
  std::optional<int64_t> Integer(std::string_view key) const;
  std::optional<mgp::List> List(std::string_view key) const;
  std::optional<bool> Bool(std::string_view key) const;
  std::optional<std::vector<double>> NumericList(std::string_view key) const;

 private:
  mgp::Map _map;
};

}  // namespace util
}  // namespace shortest_paths
