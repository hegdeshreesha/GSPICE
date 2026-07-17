#ifndef GSPICE_EXPRESSION_HPP
#define GSPICE_EXPRESSION_HPP

#include "matrix.hpp"
#include "utils.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace gspice {

class BehavioralExpression {
public:
    struct VoltageRef {
        std::string key;
        int node_pos = -1;
        int node_neg = -1;
    };

    struct CurrentRef {
        std::string name;
        int branch_index = -1;
    };

    BehavioralExpression() = default;

    BehavioralExpression(
        std::string expression,
        const std::function<int(const std::string&)>& nodeResolver)
        : expression_(stripOuterBraces(trim(expression))) {
        collectReferences(nodeResolver);
    }

    double evaluate(const VectorReal& x, double time) const {
        Evaluator evaluator(*this, x, time);
        return evaluator.parse();
    }

    const std::vector<VoltageRef>& voltageRefs() const { return voltage_refs_; }
    const std::vector<CurrentRef>& currentRefs() const { return current_refs_; }

    void bindBranch(const std::string& name, int branchIndex) {
        for (auto& ref : current_refs_) {
            if (canonical(ref.name) == canonical(name)) {
                ref.branch_index = branchIndex;
            }
        }
    }

    bool allBranchesBound(std::string& missing) const {
        for (const auto& ref : current_refs_) {
            if (ref.branch_index < 0) {
                missing = ref.name;
                return false;
            }
        }
        return true;
    }

    std::vector<int> dependencyIndices() const {
        std::vector<int> deps;
        auto add = [&deps](int idx) {
            if (idx < 0) return;
            if (std::find(deps.begin(), deps.end(), idx) == deps.end()) deps.push_back(idx);
        };
        for (const auto& ref : voltage_refs_) {
            add(ref.node_pos);
            add(ref.node_neg);
        }
        for (const auto& ref : current_refs_) {
            add(ref.branch_index);
        }
        return deps;
    }

    const std::string& text() const { return expression_; }

private:
    std::string expression_;
    std::vector<VoltageRef> voltage_refs_;
    std::vector<CurrentRef> current_refs_;

    static std::string trim(const std::string& text) {
        const auto first = text.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) return "";
        const auto last = text.find_last_not_of(" \t\r\n");
        return text.substr(first, last - first + 1);
    }

    static std::string stripOuterBraces(std::string text) {
        text = trim(text);
        if (text.size() >= 2 && ((text.front() == '{' && text.back() == '}') ||
            (text.front() == '\'' && text.back() == '\'') ||
            (text.front() == '"' && text.back() == '"'))) {
            return trim(text.substr(1, text.size() - 2));
        }
        return text;
    }

    static std::string upper(std::string text) {
        std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c) {
            return static_cast<char>(std::toupper(c));
        });
        return text;
    }

    static std::string canonical(std::string text) {
        text = trim(text);
        std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        });
        return text;
    }

    static std::vector<std::string> splitProbePayload(const std::string& payload) {
        std::vector<std::string> parts;
        size_t start = 0;
        int depth = 0;
        for (size_t i = 0; i < payload.size(); ++i) {
            const char c = payload[i];
            if (c == '(') ++depth;
            if (c == ')') --depth;
            if (c == ',' && depth == 0) {
                parts.push_back(trim(payload.substr(start, i - start)));
                start = i + 1;
            }
        }
        parts.push_back(trim(payload.substr(start)));
        return parts;
    }

    static bool identifierStart(char c) {
        return std::isalpha(static_cast<unsigned char>(c)) || c == '_' || c == '$';
    }

    static bool identifierChar(char c) {
        return std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '$' || c == ':' || c == '.';
    }

    static size_t matchingParen(const std::string& text, size_t open) {
        int depth = 0;
        for (size_t i = open; i < text.size(); ++i) {
            if (text[i] == '(') ++depth;
            if (text[i] == ')') {
                --depth;
                if (depth == 0) return i;
            }
        }
        return std::string::npos;
    }

    void addVoltageRef(const std::string& payload, const std::function<int(const std::string&)>& nodeResolver) {
        auto parts = splitProbePayload(payload);
        if (parts.empty() || parts[0].empty()) return;
        const std::string pos = trim(parts[0]);
        const std::string neg = parts.size() > 1 && !parts[1].empty() ? trim(parts[1]) : "0";
        const std::string key = canonical(pos) + "," + canonical(neg);
        for (const auto& ref : voltage_refs_) {
            if (ref.key == key) return;
        }
        voltage_refs_.push_back({key, nodeResolver(pos), nodeResolver(neg)});
    }

    void addCurrentRef(const std::string& payload) {
        const std::string name = trim(payload);
        if (name.empty()) return;
        for (const auto& ref : current_refs_) {
            if (canonical(ref.name) == canonical(name)) return;
        }
        current_refs_.push_back({name, -1});
    }

    void collectReferences(const std::function<int(const std::string&)>& nodeResolver) {
        for (size_t i = 0; i < expression_.size();) {
            if (!identifierStart(expression_[i])) {
                ++i;
                continue;
            }
            const size_t idStart = i;
            while (i < expression_.size() && identifierChar(expression_[i])) ++i;
            std::string id = upper(expression_.substr(idStart, i - idStart));
            size_t j = i;
            while (j < expression_.size() && std::isspace(static_cast<unsigned char>(expression_[j]))) ++j;
            if (j >= expression_.size() || expression_[j] != '(') continue;
            const size_t close = matchingParen(expression_, j);
            if (close == std::string::npos) continue;
            const std::string payload = expression_.substr(j + 1, close - j - 1);
            if (id == "V") addVoltageRef(payload, nodeResolver);
            if (id == "I") addCurrentRef(payload);
            i = close + 1;
        }
    }

    double voltageValue(const std::string& payload, const VectorReal& x) const {
        auto parts = splitProbePayload(payload);
        if (parts.empty() || parts[0].empty()) return 0.0;
        const std::string pos = trim(parts[0]);
        const std::string neg = parts.size() > 1 && !parts[1].empty() ? trim(parts[1]) : "0";
        const std::string key = canonical(pos) + "," + canonical(neg);
        for (const auto& ref : voltage_refs_) {
            if (ref.key == key) {
                const double vp = ref.node_pos >= 0 ? x[ref.node_pos] : 0.0;
                const double vn = ref.node_neg >= 0 ? x[ref.node_neg] : 0.0;
                return vp - vn;
            }
        }
        return 0.0;
    }

    double currentValue(const std::string& payload, const VectorReal& x) const {
        const std::string name = canonical(payload);
        for (const auto& ref : current_refs_) {
            if (canonical(ref.name) == name) {
                return ref.branch_index >= 0 ? x[ref.branch_index] : 0.0;
            }
        }
        return 0.0;
    }

    class Evaluator {
    public:
        Evaluator(const BehavioralExpression& owner, const VectorReal& x, double time)
            : owner_(owner), x_(x), time_(time) {}

        double parse() {
            pos_ = 0;
            const double value = parseTernary();
            skipSpace();
            if (pos_ != owner_.expression_.size()) {
                throw std::runtime_error("Unexpected expression text near '" + owner_.expression_.substr(pos_) + "'");
            }
            return value;
        }

    private:
        const BehavioralExpression& owner_;
        const VectorReal& x_;
        double time_ = 0.0;
        size_t pos_ = 0;

        void skipSpace() {
            while (pos_ < owner_.expression_.size() && std::isspace(static_cast<unsigned char>(owner_.expression_[pos_]))) ++pos_;
        }

        bool consume(char c) {
            skipSpace();
            if (pos_ < owner_.expression_.size() && owner_.expression_[pos_] == c) {
                ++pos_;
                return true;
            }
            return false;
        }

        bool consumeText(const std::string& text) {
            skipSpace();
            if (owner_.expression_.substr(pos_, text.size()) == text) {
                pos_ += text.size();
                return true;
            }
            return false;
        }

        double parseTernary() {
            double condition = parseComparison();
            if (consume('?')) {
                const double whenTrue = parseTernary();
                if (!consume(':')) throw std::runtime_error("Expected ':' in ternary expression");
                const double whenFalse = parseTernary();
                return condition != 0.0 ? whenTrue : whenFalse;
            }
            return condition;
        }

        double parseComparison() {
            double lhs = parseAdditive();
            while (true) {
                if (consumeText(">=")) lhs = lhs >= parseAdditive() ? 1.0 : 0.0;
                else if (consumeText("<=")) lhs = lhs <= parseAdditive() ? 1.0 : 0.0;
                else if (consumeText("==")) lhs = lhs == parseAdditive() ? 1.0 : 0.0;
                else if (consumeText("!=")) lhs = lhs != parseAdditive() ? 1.0 : 0.0;
                else if (consume('>')) lhs = lhs > parseAdditive() ? 1.0 : 0.0;
                else if (consume('<')) lhs = lhs < parseAdditive() ? 1.0 : 0.0;
                else break;
            }
            return lhs;
        }

        double parseAdditive() {
            double value = parseMultiplicative();
            while (true) {
                if (consume('+')) value += parseMultiplicative();
                else if (consume('-')) value -= parseMultiplicative();
                else break;
            }
            return value;
        }

        double parseMultiplicative() {
            double value = parsePower();
            while (true) {
                if (consume('*')) value *= parsePower();
                else if (consume('/')) {
                    const double denom = parsePower();
                    value /= denom;
                } else {
                    break;
                }
            }
            return value;
        }

        double parsePower() {
            double value = parseUnary();
            if (consume('^')) {
                value = std::pow(value, parsePower());
            }
            return value;
        }

        double parseUnary() {
            if (consume('+')) return parseUnary();
            if (consume('-')) return -parseUnary();
            return parsePrimary();
        }

        double parsePrimary() {
            skipSpace();
            if (consume('(')) {
                const double value = parseTernary();
                if (!consume(')')) throw std::runtime_error("Expected ')' in expression");
                return value;
            }
            if (pos_ >= owner_.expression_.size()) throw std::runtime_error("Unexpected end of expression");
            if (std::isdigit(static_cast<unsigned char>(owner_.expression_[pos_])) || owner_.expression_[pos_] == '.') {
                return parseNumber();
            }
            if (identifierStart(owner_.expression_[pos_])) {
                return parseIdentifier();
            }
            throw std::runtime_error("Unexpected character in expression");
        }

        double parseNumber() {
            const size_t start = pos_;
            bool sawExp = false;
            while (pos_ < owner_.expression_.size()) {
                char c = owner_.expression_[pos_];
                if (std::isdigit(static_cast<unsigned char>(c)) || c == '.') {
                    ++pos_;
                } else if ((c == 'e' || c == 'E') && !sawExp) {
                    sawExp = true;
                    ++pos_;
                    if (pos_ < owner_.expression_.size() && (owner_.expression_[pos_] == '+' || owner_.expression_[pos_] == '-')) ++pos_;
                } else {
                    break;
                }
            }
            while (pos_ < owner_.expression_.size() && std::isalpha(static_cast<unsigned char>(owner_.expression_[pos_]))) ++pos_;
            return Utils::parseValue(owner_.expression_.substr(start, pos_ - start));
        }

        double parseIdentifier() {
            const size_t start = pos_;
            while (pos_ < owner_.expression_.size() && identifierChar(owner_.expression_[pos_])) ++pos_;
            const std::string id = owner_.expression_.substr(start, pos_ - start);
            const std::string idUpper = upper(id);
            skipSpace();
            if (consume('(')) {
                if (idUpper == "V" || idUpper == "I") {
                    const size_t payloadStart = pos_;
                    int depth = 1;
                    while (pos_ < owner_.expression_.size() && depth > 0) {
                        if (owner_.expression_[pos_] == '(') ++depth;
                        else if (owner_.expression_[pos_] == ')') --depth;
                        if (depth > 0) ++pos_;
                    }
                    if (depth != 0) throw std::runtime_error("Unmatched probe parentheses in expression");
                    const std::string payload = owner_.expression_.substr(payloadStart, pos_ - payloadStart);
                    ++pos_;
                    return idUpper == "V" ? owner_.voltageValue(payload, x_) : owner_.currentValue(payload, x_);
                }
                std::vector<double> args;
                if (!consume(')')) {
                    do {
                        args.push_back(parseTernary());
                    } while (consume(','));
                    if (!consume(')')) throw std::runtime_error("Expected ')' after function arguments");
                }
                return applyFunction(idUpper, args);
            }
            if (idUpper == "TIME" || idUpper == "T") return time_;
            if (idUpper == "PI") return 3.14159265358979323846;
            if (idUpper == "E") return 2.71828182845904523536;
            throw std::runtime_error("Unknown expression identifier '" + id + "'");
        }

        static double arg(const std::vector<double>& args, size_t i, const std::string& name) {
            if (i >= args.size()) throw std::runtime_error("Function " + name + " has too few arguments");
            return args[i];
        }

        static double applyFunction(const std::string& name, const std::vector<double>& args) {
            if (name == "SIN") return std::sin(arg(args, 0, name));
            if (name == "COS") return std::cos(arg(args, 0, name));
            if (name == "TAN") return std::tan(arg(args, 0, name));
            if (name == "EXP") return std::exp(arg(args, 0, name));
            if (name == "LOG" || name == "LN") return std::log(arg(args, 0, name));
            if (name == "LOG10") return std::log10(arg(args, 0, name));
            if (name == "SQRT") return std::sqrt(arg(args, 0, name));
            if (name == "ABS") return std::abs(arg(args, 0, name));
            if (name == "POW") return std::pow(arg(args, 0, name), arg(args, 1, name));
            if (name == "MIN") return std::min(arg(args, 0, name), arg(args, 1, name));
            if (name == "MAX") return std::max(arg(args, 0, name), arg(args, 1, name));
            if (name == "LIMIT" || name == "CLAMP") {
                return std::clamp(arg(args, 0, name), arg(args, 1, name), arg(args, 2, name));
            }
            if (name == "IF") return arg(args, 0, name) != 0.0 ? arg(args, 1, name) : arg(args, 2, name);
            if (name == "SGN" || name == "SIGN") return (arg(args, 0, name) > 0.0) - (arg(args, 0, name) < 0.0);
            if (name == "U" || name == "STEP") return arg(args, 0, name) >= 0.0 ? 1.0 : 0.0;
            if (name == "URAMP") return std::max(0.0, arg(args, 0, name));
            if (name == "FLOOR") return std::floor(arg(args, 0, name));
            if (name == "CEIL" || name == "CEILING") return std::ceil(arg(args, 0, name));
            if (name == "ROUND") return std::round(arg(args, 0, name));
            throw std::runtime_error("Unsupported expression function '" + name + "'");
        }
    };
};

} // namespace gspice

#endif // GSPICE_EXPRESSION_HPP
