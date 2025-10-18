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

std::string generateTitleForYieldPlot(int iq2, int iw, int ict) {
    // Границы по Q²
    double q2Min = STEPS_Q2[2 * iq2];
    double q2Max = STEPS_Q2[2 * iq2 + 1];

    // Границы по W
    double wMin = W_MIN + STEP_W * iw;
    double wMax = W_MIN + STEP_W * (iw + 1);

    // Границы по cos(θ)
    double cosThetaMin = -1.0 + STEP_COS_THETA * ict;
    double cosThetaMax = -1.0 + STEP_COS_THETA * (ict + 1);

    /*
    return "Yield vs #phi: Q^{2} = " + formatFloat(q2Min, 1) + "-" + formatFloat(q2Max, 1) +
        " GeV^{2}, W = " + formatFloat(wMin, 2) + "-" + formatFloat(wMax, 2) + " GeV, " +
        "cos(#theta) = " + formatFloat(cosThetaMin, 1) + " to " + formatFloat(cosThetaMax, 1) + "; #phi, rad; Yield";
    */
    return "Yield vs #phi: cos(#theta) = " + formatFloat(cosThetaMin, 1) + " to " + formatFloat(cosThetaMax, 1) + "; #phi, rad; Yield";
}

std::string generateTitleForEfficiencyPlot(int iq2, int iw, int ict) {
    // Границы по cos(θ)
    double cosThetaMin = -1.0 + STEP_COS_THETA * ict;
    double cosThetaMax = -1.0 + STEP_COS_THETA * (ict + 1);

    return "Efficiency vs #phi: cos(#theta) = " + formatFloat(cosThetaMin, 1) + " to " + formatFloat(cosThetaMax, 1) + "; #phi, rad; Yield";
}