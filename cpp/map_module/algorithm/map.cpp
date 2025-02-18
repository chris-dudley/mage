#include <fmt/format.h>
#include <list>
#include <sstream>

#include "map.hpp"

const auto number_of_elements_in_pair = 2;

void Map::FromNodes(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory) {
  mgp::memory = memory;
  const auto arguments = mgp::List(args);
  const auto record_factory = mgp::RecordFactory(result);

  try {
    const auto label{arguments[0].ValueString()};
    const auto property{arguments[1].ValueString()};
    mgp::Map map_result{};

    const auto all_nodes = mgp::Graph(memgraph_graph).Nodes();
    for (const auto node : all_nodes) {
      if (!node.HasLabel(label) || !node.Properties().contains(std::string(property))) continue;

      std::ostringstream oss;
      oss << node.GetProperty(std::string(property));
      const auto key = oss.str();

      mgp::Map map{};
      map.Update("identity", mgp::Value(node.Id().AsInt()));

      mgp::List labels{};
      for (const auto &label : node.Labels()) {
        labels.AppendExtend(mgp::Value(label));
      }
      map.Update("labels", mgp::Value(std::move(labels)));

      const auto property_map = node.Properties();
      mgp::Map properties{};
      for (const auto &[key, value] : property_map) {
        properties.Insert(key, value);
      }
      map.Update("properties", mgp::Value(std::move(properties)));

      map_result.Update(key, mgp::Value(std::move(map)));
    }

    auto record = record_factory.NewRecord();
    record.Insert(std::string(kResultFromNodes).c_str(), map_result);

  } catch (const std::exception &e) {
    record_factory.SetErrorMessage(e.what());
    return;
  }
}

void Map::FromValues(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory) {
  mgp::memory = memory;
  const auto arguments = mgp::List(args);
  const auto record_factory = mgp::RecordFactory(result);

  try {
    const auto values{arguments[0].ValueList()};
    mgp::Map map{};

    if (values.Size() % 2) {
      throw mgp::ValueException("List needs to have an even number of elements");
    }

    auto iterator = values.begin();
    while (iterator != values.end()) {
      std::ostringstream oss;
      oss << *iterator;
      const auto key = oss.str();

      ++iterator;
      map.Update(key, *iterator);
      ++iterator;
    }

    auto record = record_factory.NewRecord();
    record.Insert(std::string(kResultFromValues).c_str(), map);

  } catch (const std::exception &e) {
    record_factory.SetErrorMessage(e.what());
    return;
  }
}

void Map::SetKey(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory) {
  mgp::memory = memory;
  const auto arguments = mgp::List(args);
  const auto record_factory = mgp::RecordFactory(result);

  try {
    auto map = arguments[0].ValueMap();
    const auto key{arguments[1].ValueString()};
    const auto value{arguments[2]};
    map.Update(key, std::move(value));

    auto record = record_factory.NewRecord();
    record.Insert(std::string(kResultSetKey).c_str(), std::move(map));

  } catch (const std::exception &e) {
    record_factory.SetErrorMessage(e.what());
    return;
  }
}

void Map::RemoveRecursion(mgp::Map &result, bool recursive, std::string_view key) {
  for (auto element : result) {
    if (element.key == key) {
      result.Erase(element.key);
      continue;
    }
    if (element.value.IsMap() && recursive) {
      // TO-DO no need for non_const_value_map in new version of memgraph
      mgp::Map non_const_value_map = mgp::Map(std::move(element.value.ValueMap()));
      RemoveRecursion(non_const_value_map, recursive, key);
      if (non_const_value_map.Empty()) {
        result.Erase(element.key);
        continue;
      }
      result.Update(element.key, mgp::Value(std::move(non_const_value_map)));
    }
  }
}

void Map::RemoveKey(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory) {
  mgp::memory = memory;
  const auto arguments = mgp::List(args);
  const auto record_factory = mgp::RecordFactory(result);
  try {
    const auto map = arguments[0].ValueMap();
    const auto key = arguments[1].ValueString();
    const auto recursive = arguments[2].ValueBool();

    mgp::Map map_removed = mgp::Map(std::move(map));

    RemoveRecursion(map_removed, recursive, key);

    auto record = record_factory.NewRecord();
    record.Insert(std::string(kResultRemoveKey).c_str(), std::move(map_removed));

  } catch (const std::exception &e) {
    record_factory.SetErrorMessage(e.what());
    return;
  }
}

void Map::FromPairs(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory) {
  mgp::memory = memory;
  const auto arguments = mgp::List(args);
  const auto record_factory = mgp::RecordFactory(result);
  try {
    const auto list = arguments[0].ValueList();

    mgp::Map pairs_map;

    for (const auto inside_list : list) {
      if (inside_list.ValueList().Size() != number_of_elements_in_pair) {
        throw mgp::IndexException(
            fmt::format("Pairs must consist of {} elements exactly.", number_of_elements_in_pair));
      }
      if (!inside_list.ValueList()[0].IsString()) {
        throw mgp::ValueException("All keys have to be type string.");
      }
      pairs_map.Update(inside_list.ValueList()[0].ValueString(), std::move(inside_list.ValueList()[1]));
    }

    auto record = record_factory.NewRecord();
    record.Insert(std::string(kResultFromPairs).c_str(), std::move(pairs_map));

  } catch (const std::exception &e) {
    record_factory.SetErrorMessage(e.what());
    return;
  }
}

void Map::Merge(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory) {
  mgp::memory = memory;
  const auto arguments = mgp::List(args);
  const auto record_factory = mgp::RecordFactory(result);
  try {
    const auto map1 = arguments[0].ValueMap();
    const auto map2 = arguments[1].ValueMap();

    mgp::Map merged_map = mgp::Map(std::move(map2));

    for (const auto element : map1) {
      if (merged_map.At(element.key).IsNull()) {
        merged_map.Insert(element.key, element.value);
      }
    }

    auto record = record_factory.NewRecord();
    record.Insert(std::string(kResultMerge).c_str(), std::move(merged_map));

  } catch (const std::exception &e) {
    record_factory.SetErrorMessage(e.what());
    return;
  }
}

void Map::FlattenRecursion(mgp::Map &result, const mgp::Map &input, const std::string &key,
                           const std::string &delimiter) {
  for (auto element : input) {
    std::string el_key(element.key);
    if (element.value.IsMap()) {
      FlattenRecursion(result, element.value.ValueMap(), key + el_key + delimiter, delimiter);
    } else {
      result.Insert(key + el_key, element.value);
    }
  }
}

void Map::Flatten(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory) {
  mgp::memory = memory;
  const auto arguments = mgp::List(args);
  const auto record_factory = mgp::RecordFactory(result);
  try {
    const mgp::Map map = arguments[0].ValueMap();
    const std::string delimiter(arguments[1].ValueString());
    mgp::Map result_map = mgp::Map();
    FlattenRecursion(result_map, map, "", delimiter);
    auto record = record_factory.NewRecord();
    record.Insert(std::string(Map::kReturnValueFlatten).c_str(), std::move(result_map));

  } catch (const std::exception &e) {
    record_factory.SetErrorMessage(e.what());
    return;
  }
}

void Map::FromLists(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory) {
  mgp::memory = memory;
  auto arguments = mgp::List(args);
  const auto record_factory = mgp::RecordFactory(result);
  try {
    mgp::List list1 = arguments[0].ValueList();
    mgp::List list2 = arguments[1].ValueList();

    const auto expected_list_size = list1.Size();
    if (expected_list_size != list2.Size() || expected_list_size == 0) {
      throw mgp::ValueException("Lists must be of same size and not empty");
    }
    mgp::Map result = mgp::Map();
    for (size_t i = 0; i < expected_list_size; i++) {
      result.Update(std::move(list1[i].ValueString()), std::move(list2[i]));
    }
    auto record = record_factory.NewRecord();
    record.Insert(std::string(Map::kReturnListFromLists).c_str(), std::move(result));

  } catch (const std::exception &e) {
    record_factory.SetErrorMessage(e.what());
    return;
  }
}

void Map::RemoveRecursionSet(mgp::Map &result, bool recursive, std::unordered_set<std::string_view> &set) {
  for (auto element : result) {
    bool inSet = false;
    if (set.find(element.key) != set.end()) {
      inSet = true;
    }
    if (inSet) {
      result.Erase(element.key);
      continue;
    }
    if (element.value.IsMap() && recursive) {
      mgp::Map non_const_value_map = mgp::Map(std::move(element.value.ValueMap()));
      RemoveRecursionSet(non_const_value_map, recursive, set);
      if (non_const_value_map.Empty()) {
        result.Erase(element.key);
        continue;
      }
      result.Update(element.key, mgp::Value(std::move(non_const_value_map)));
    }
  }
}

void Map::RemoveKeys(mgp_list *args, mgp_graph *memgraph_graph, mgp_result *result, mgp_memory *memory) {
  mgp::memory = memory;
  auto arguments = mgp::List(args);
  const auto record_factory = mgp::RecordFactory(result);
  try {
    mgp::Map map = arguments[0].ValueMap();
    const mgp::List list = arguments[1].ValueList();
    bool recursive = arguments[2].ValueBool();
    std::unordered_set<std::string_view> set;
    for (auto elem : list) {
      set.insert(std::move(elem.ValueString()));
    }
    RemoveRecursionSet(map, recursive, set);
    auto record = record_factory.NewRecord();
    record.Insert(std::string(Map::kReturnRemoveKeys).c_str(), std::move(map));

  } catch (const std::exception &e) {
    record_factory.SetErrorMessage(e.what());
    return;
  }
}
