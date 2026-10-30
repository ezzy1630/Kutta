#pragma once

#include <QString>

#include <optional>
#include <utility>
#include <variant>

namespace kapp {

// Value or a plain-language error. Parsing at the process/file boundary
// returns this; everything past the boundary works with parsed types only.
template <class T>
class Result {
public:
    Result(T value) : m_v(std::move(value)) {}
    static Result failure(QString message) { return Result(Error{std::move(message)}); }

    bool ok() const { return std::holds_alternative<T>(m_v); }
    explicit operator bool() const { return ok(); }
    const T& value() const& { return std::get<T>(m_v); }
    T& value() & { return std::get<T>(m_v); }
    T&& value() && { return std::get<T>(std::move(m_v)); }
    const T* operator->() const { return &std::get<T>(m_v); }
    T* operator->() { return &std::get<T>(m_v); }
    QString error() const { return ok() ? QString() : std::get<Error>(m_v).message; }

private:
    struct Error { QString message; };
    explicit Result(Error e) : m_v(std::move(e)) {}
    std::variant<T, Error> m_v;
};

} // namespace kapp
