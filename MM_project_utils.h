#pragma once

#include <string>
#include <tuple>
#include <regex>
#include <sstream>
#include <iomanip>

#include "MM_project_const.h"

/**
 * @brief Генерация стандартного имени ячейки вида "CELL i-j-k-l"
 */
inline std::string makeCellName(int i, int j, int k, int l) {
    return "CELL " + std::to_string(i + 1) + "-" + std::to_string(j + 1) + "-" +
           std::to_string(k + 1) + "-" + std::to_string(l + 1);
}

/**
 * @brief Извлекает индексы (i, j, k, l) из строки формата "CELL i-j-k-l".
 * Возвращает true при успешном извлечении, false в противном случае.
 */
inline bool parseCellName(const std::string& name, int& i, int& j, int& k, int& l) {
    std::smatch match;
    std::regex pattern(R"(.*?(\d+)-(\d+)-(\d+)-(\d+))");

    if (std::regex_search(name, match, pattern) && match.size() == 5) {
        i = std::stoi(match[1]) - 1;
        j = std::stoi(match[2]) - 1;
        k = std::stoi(match[3]) - 1;
        l = std::stoi(match[4]) - 1;
        return true;
    }
    return false;
}


/**
 * @brief Преобразует число с плавающей точкой в строку с заданной точностью.
 * 
 * @param value Значение для форматирования
 * @param precision Количество знаков после запятой (по умолчанию 2)
 * @return std::string Отформатированное число
 */
inline std::string formatFloat(double value, int precision = 2) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << value;
    return out.str();
}

auto cellGenerator = [i = 0, j = 0, k = 0, l = 0]() mutable -> std::function<bool(std::string& nextName)> {
    return [=](std::string& name) mutable -> bool {
        if (i >= NUMBER_Q2) return false;
        name = makeCellName(i, j, k, l);
        // обновляем индексы
        if (++l >= NUMBER_PHI) {
            l = 0;
            if (++k >= NUMBER_COS_THETA) {
                k = 0;
                if (++j >= NUMBER_W) {
                    j = 0;
                    ++i;
                }
            }
        }
        return true;
    };
}();
