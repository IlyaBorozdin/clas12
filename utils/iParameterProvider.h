#pragma once

#include <string>
#include <limits>
#include <cmath>

#include "TF1.h"

#include "iHistogramAnalyzer.h"

/**
 * @brief Структура, описывающая параметр фита.
 * Содержит значение, флаг фиксации и границы.
 */
struct FitParameter {
    std::string name; // Имя параметра (для логов и идентификации)
    double val;       // Значение
    bool fix;         // Флаг фиксации
    double lower;     // Нижняя граница
    double upper;     // Верхняя граница

    // Статическая константа для "пустой" границы
    static constexpr double NO_LIMIT = std::numeric_limits<double>::infinity();

    FitParameter(std::string n = "unnamed", double v = 0.0, bool f = false, double low = -NO_LIMIT, double high = NO_LIMIT) 
        : name(n), val(v), fix(f), lower(low), upper(high) {}

    // --- Хелперы для проверки границ (пункт 1) ---

    bool hasLowerLimit() const { return lower != -NO_LIMIT && !std::isnan(lower); }
    bool hasUpperLimit() const { return upper !=  NO_LIMIT && !std::isnan(upper); }
    
    // Проверка для ROOT: заданы ли обе границы
    bool hasBothLimits() const { return hasLowerLimit() && hasUpperLimit(); }

    /**
     * @brief Применяет настройки данного параметра к функции ROOT TF1.
     * @param func Указатель на функцию ROOT.
     * @param parIdx Индекс параметра в функции.
     */
    void apply(TF1* func, int parIdx) const {
        if (!func) return;

        // 1. Устанавливаем имя (помогает при отладке в ROOT Browser)
        func->SetParName(parIdx, name.c_str());

        // 2. Устанавливаем значение
        func->SetParameter(parIdx, val);

        // 3. Обработка фиксации или лимитов
        if (fix) {
            func->FixParameter(parIdx, val);
        } else {
            // Если есть хотя бы одна граница, вызываем SetParLimits
            if (hasLowerLimit() || hasUpperLimit()) {
                // ROOT требует оба числа. Используем очень большие значения 
                // вместо бесконечности, так как MINUIT иногда плохо переваривает inf
                double l = hasLowerLimit() ? lower : -1e10; 
                double u = hasUpperLimit() ? upper : 1e10;
                func->SetParLimits(parIdx, l, u);
            } else {
                // Если границ нет совсем, снимаем лимиты (на случай, если они были там раньше)
                func->ReleaseParameter(parIdx);
            }
        }
    }
};

/**
 * @brief Интерфейс для получения параметров фитирования.
 * Позволяет абстрагироваться от источника данных (JSON, БД, алгоритм).
 */
class IParameterProvider {
protected:
    IHistogramAnalyzer* m_analyzer; // Инструментарий для анализа
public:
    IParameterProvider(IHistogramAnalyzer* analyzer = nullptr) 
        : m_analyzer(analyzer) {}

    virtual ~IParameterProvider() = default;

    void setAnalyzer(IHistogramAnalyzer* analyzer) { m_analyzer = analyzer; }

    // Специфические параметры (например, начальные приближения для формул)
    virtual FitParameter getParameter(const std::string& name, int i, int j, int k, int l, TH1F* hist = nullptr) const = 0;
    
    // Проверка наличия параметра
    virtual bool hasParameter(const std::string& name, int i, int j, int k, int l) const = 0;

    virtual void getBounds(int i, int j, int k, int l, double& low, double& high) const = 0;
};

// class DefaultParameterProvider : public IParameterProvider {
// public:
//     DefaultParameterProvider() = default;
//     virtual ~DefaultParameterProvider() = default;

//     FitParameter getParameter(const std::string& name, int i, int j, int k, int l) const override {
//         return FitParameter(); 
//     }

//     bool hasParameter(const std::string& name, int i, int j, int k, int l) const override {
//         return false;
//     }
// };
