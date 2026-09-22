#ifndef QTPLANET_FILTERCHAIN_H
#define QTPLANET_FILTERCHAIN_H

#include "chrono_planet/ChApiPlanet.h"

#include <functional>
#include <type_traits>
#include <utility>

namespace qtplanet {

// Eager, synchronous value transformations. Each filter consumes the current
// value and returns an owned result of the next type. No shared context, stored
// callbacks or type erasure. Exceptions propagate and stop the chain.
template <typename Value>
class FilterChain {
    static_assert(std::is_object_v<Value> && !std::is_reference_v<Value>,
                  "A filter chain must own its value");

public:
    explicit FilterChain(Value value) : value_(std::move(value)) {}
    FilterChain(const FilterChain&) = delete;
    FilterChain& operator=(const FilterChain&) = delete;
    FilterChain(FilterChain&&) = default;
    FilterChain& operator=(FilterChain&&) = default;

    // Only temporary or explicitly moved chains can advance, avoiding implicit
    // copies of large buffers. A filter may accept its input by value or T&&.
    template <typename Filter>
    [[nodiscard]] auto then(Filter&& filter) && {
        using Result = std::invoke_result_t<Filter, Value&&>;
        static_assert(std::is_object_v<Result> && !std::is_reference_v<Result>,
                      "A filter must return an owned value, not void or a reference");
        return FilterChain<Result>(std::invoke(std::forward<Filter>(filter), std::move(value_)));
    }

    [[nodiscard]] Value take() && { return std::move(value_); }

private:
    Value value_;
};

// Copies lvalue input; use std::move to transfer an existing buffer or resource.
template <typename Value>
[[nodiscard]] FilterChain<Value> filterChain(Value value) {
    return FilterChain<Value>(std::move(value));
}

}   // namespace qtplanet

#endif   // QTPLANET_FILTERCHAIN_H
