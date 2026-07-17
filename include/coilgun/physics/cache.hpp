/**
 * @file cache.hpp
 * @brief Fixed-capacity LRU cache for memoising expensive physics computations.
 * @author Winston Meursault
 */

#pragma once

#include <cstddef>
#include <list>
#include <unordered_map>
#include <utility>

namespace coilgun::physics {

/**
 * @brief Fixed-capacity LRU (Least Recently Used) cache.
 * @tparam Key Cache key type.
 * @tparam Value Cache value type.
 * @tparam MaxSize Maximum number of entries (default 4096).
 * @tparam KeyHash Hash function for the key type.
 */
template<typename Key, typename Value, std::size_t MaxSize = 4096,
         typename KeyHash = std::hash<Key>>
class LRUCache {
public:
    /**
     * @brief Retrieve a value from the cache.
     * @param key Lookup key.
     * @param[out] value Set to the cached value on hit.
     * @return true if the key was found, false otherwise.
     */
    bool get(const Key& key, Value& value) {
        auto it = map_.find(key);
        if (it == map_.end()) {
            return false;
        }
        list_.splice(list_.begin(), list_, it->second);
        value = it->second->second;
        return true;
    }

    /**
     * @brief Insert or update a key-value pair.
     * @param key Lookup key.
     * @param value Value to store.
     *
     * Evicts the least recently used entry if at capacity.
     */
    void put(const Key& key, const Value& value) {
        auto it = map_.find(key);
        if (it != map_.end()) {
            list_.splice(list_.begin(), list_, it->second);
            it->second->second = value;
            return;
        }
        if (list_.size() >= MaxSize) {
            auto last = list_.end();
            --last;
            map_.erase(last->first);
            list_.pop_back();
        }
        list_.emplace_front(key, value);
        map_[key] = list_.begin();
    }

    /**
     * @brief Current number of entries in the cache.
     * @return Entry count.
     */
    std::size_t size() const {
        return list_.size();
    }

    /**
     * @brief Remove all entries from the cache.
     */
    void clear() {
        list_.clear();
        map_.clear();
    }

private:
    using ListIterator =
        typename std::list<std::pair<Key, Value>>::iterator;

    std::list<std::pair<Key, Value>> list_;
    std::unordered_map<Key, ListIterator, KeyHash> map_;
};

} // namespace coilgun::physics
