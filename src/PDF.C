#include "PDF.h"

nontefonte::Function::Function(const std::string &expr,
                           const std::vector<std::string> &variable_names)
        : _expression(expr), _pos(0), _variables(variable_names) {}

double nontefonte::Function::parseExpression() {
    auto value = parseTerm();
    while (true) {
        skipWhitespace();
        if (match('+'))
            value += parseTerm();
        else if (match('-'))
            value -= parseTerm();
        else
            break;
    }
    return value;
}

double nontefonte::Function::parseTerm() {
    double value = parseFactor();
    while (true) {
        skipWhitespace();
        if (match('*')) {
            value *= parseFactor();
        } else if (match('/')) {
            value /= parseFactor();
        } else {
            break;
        }
    }
    return value;
}

double nontefonte::Function::parseFactor() {
    double base = parseUnary();
    skipWhitespace();
    if (match('^')) {
        return std::pow(base, parseFactor()); // right associative
    }
    return base;
}

double nontefonte::Function::parseUnary() {
    skipWhitespace();
    if (match('+'))
        return parseUnary();
    if (match('-'))
        return -parseUnary();
    return parsePrimary();
}

double nontefonte::Function::parsePrimary() {
    skipWhitespace();

    if (match('(')) {
        double value = parseExpression();
        if (!match(')'))
            throw std::runtime_error("Missing ')'");
        return value;
    }

    if (std::isalpha(peek())) {
        std::string name = parseIdentifier();

        skipWhitespace();
        if (match('(')) {
            double arg = parseExpression();
            if (!match(')'))
                throw std::runtime_error("Missing ')'");
            return callFunction(name, arg);
        }

        return getVariable(name);
    }

    return parseNumber();
}

double nontefonte::Function::parseNumber() {
    skipWhitespace();
    size_t start = _pos;

    while (_pos < _expression.size() &&
           (std::isdigit(_expression[_pos]) || _expression[_pos] == '.'))
        _pos++;

    if (start == _pos)
        throw std::runtime_error("Expected number");

    return std::stod(_expression.substr(start, _pos - start));
}

std::string nontefonte::Function::parseIdentifier() {
    size_t start = _pos;
    while (_pos < _expression.size() &&
           (std::isalnum(_expression[_pos]) || _expression[_pos] == '_'))
        _pos++;
    return _expression.substr(start, _pos - start);
}

double nontefonte::Function::getVariable(const std::string &name) {
    for (size_t i = 0; i < _variables.size(); ++i)
        if (_variables[i] == name)
            return _values[i];

    throw std::runtime_error("Unknown variable: " + name);
}

double nontefonte::Function::callFunction(const std::string &name, double arg) {
    return _functions.at(name)(arg);
}

bool nontefonte::Function::match(char c) {
    if (_pos < _expression.size() && _expression[_pos] == c) {
        ++_pos;
        return true;
    }
    return false;
}

char nontefonte::Function::peek() const {
    return _pos < _expression.size() ? _expression[_pos] : '\0';
}

void nontefonte::Function::skipWhitespace() {
    while (_pos < _expression.size() && std::isspace(_expression[_pos]))
        ++_pos;
}