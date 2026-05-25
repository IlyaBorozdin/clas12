#pragma once

#include <TH1F.h>

/**
 * @brief Интерфейс для извлечения характеристик из гистограмм.
 */
class IHistogramAnalyzer {
public:
    virtual ~IHistogramAnalyzer() = default;

    // Счётчики и интегралы
    virtual double getIntegral(TH1F* h, double xmin, double xmax) const = 0;
    virtual int countEntries(TH1F* h, double xmin, double xmax) const = 0;

    // Работа с пиками
    // Мы возвращаем void, но передаем ссылки для записи позиции и высоты
    virtual void findPeak(TH1F* h, double xmin, double xmax, double& pos, double& val) const = 0;
    virtual void getGlobalMax(TH1F* hist, double xmin, double xmax, double& pos, double& val) const = 0;

    virtual void constrainToAxis(TH1F* hist, double& low, double& high) const = 0;
};