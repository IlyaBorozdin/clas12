#pragma once

#include <fstream>
#include <nlohmann/json.hpp>
#include "TH1F.h"
#include "TF1.h"

#include "fittingStrategy.h"

using json = nlohmann::json;

class HistogramFitterUniversal : public FittingStrategy {
private:
    json config;

    // Рекурсивное слияние объектов для гранулярного обновления параметров
    void deepMerge(json& target, const json& source) const {
        for (auto it = source.begin(); it != source.end(); ++it) {
            if (it.value().is_object() && target.contains(it.key()) && target[it.key()].is_object()) {
                deepMerge(target[it.key()], it.value());
            } else {
                target[it.key()] = it.value();
            }
        }
    }

    json getEffectiveParams(int i, int j, int k, int l) const {
        double q2, w, cos_theta, phi;
        convertNumbersToVariable(i, j, k, l, q2, w, cos_theta, phi);

        // Карта для сбора параметров: priority -> список объектов параметров
        std::map<int, std::vector<json>> priorityGroups;

        // 1. Проходим по верхнеуровневым блокам (обычно это Priority 0)
        for (const auto& block : config) {
            if (isInsideBoundaries(block["boundaries"], q2, w, cos_theta, phi)) {
                
                // Добавляем параметры самого блока (P0)
                if (block.contains("parameters")) {
                    priorityGroups[block.value("priority", 0)].push_back(block["parameters"]);
                }

                // 2. Если внутри есть вложенные правила (rules) — проверяем их
                if (block.contains("rules") && block["rules"].is_array()) {
                    for (const auto& rule : block["rules"]) {
                        if (isInsideBoundaries(rule["boundaries"], q2, w, cos_theta, phi)) {
                            priorityGroups[rule.value("priority", 1)].push_back(rule["parameters"]);
                        }
                    }
                }
            }
        }

        // Сливаем всё в один объект (Deep Merge)
        json finalParams = json::object();
        for (auto const& [priority, rulesStack] : priorityGroups) {
            for (const auto& params : rulesStack) {
                deepMerge(finalParams, params);
            }
        }
        return finalParams;
    }

    // Вспомогательный метод для настройки параметров TF1
    void configureParameter(TF1* func, int parIdx, const std::string& name, const json& pData) const {
        if (pData.is_null()) return;

        // Если это просто число (как в downEdge), а не объект
        // if (!pData.is_object()) {
        //     func->SetParameter(parIdx, pData.get<double>());
        //     return;
        // }

        // Установка значения
        if (pData.contains("val")) {
            func->SetParameter(parIdx, pData["val"].get<double>());
        }

        // Установка лимитов
        // if (pData.contains("limits")) {
        //     double low = pData["limits"][0].is_null() ? -1e9 : pData["limits"][0].get<double>();
        //     double high = pData["limits"][1].is_null() ? 1e9 : pData["limits"][1].get<double>();
        //     func->SetParLimits(parIdx, low, high);
        // }

        // Фиксация
        if (pData.value("fix", false)) {
            func->FixParameter(parIdx, func->GetParameter(parIdx));
        }
    }

    bool isInsideBoundaries(const json& bounds, double q2, double w, double cos, double phi) const {
        auto checkSingle = [&](const json& b) {
            auto rangeMatch = [](const json& range_list, double val) {
                if (range_list.is_null() || !range_list.is_array()) return true;
                for (const auto& range : range_list) {
                    if (val >= range[0].get<double>() && val <= range[1].get<double>()) return true;
                }
                return false;
            };
            return rangeMatch(b["Q2"], q2) && rangeMatch(b["W"], w) && 
                rangeMatch(b["cos_theta"], cos) && rangeMatch(b["phi"], phi);
        };

        if (bounds.is_array()) {
            for (const auto& b : bounds) {
                if (checkSingle(b)) return true;
            }
            return false;
        }
        return checkSingle(bounds);
    }

public:
    HistogramFitterUniversal(const std::string& configPath) {
        std::ifstream f(configPath);
        f >> config;
    }

    TF1* fit(TH1F* hist, int i, int j, int k, int l) const override {
        json p = getEffectiveParams(i, j, k, l);
        
        double xMin = getDownEdge(hist, i, j, k, l);
        double xMax = getUpEdge(hist, i, j, k, l);

        double maxPosition, maxValue;
        getNeutronPeak(hist, i, j, k, l, maxPosition, maxValue);

        // Создаем универсальную функцию (например, Gauss + Gauss + Parabolic Pol2)
        // Название и формула зависят от вашей модели
        TF1* f = new TF1("fit_func", "[0]*exp(-0.5*((x-[1])/[2]*(x < [1] ? [3] : 1))^2) + TMath::Max(0., [4]*(x - [5])*(x - [6])) + [7]*exp(-0.5*((x-[8])/[9])^2)", xMin, xMax);

        // Пример маппинга параметров из JSON на индексы TF1
        f->SetParameter(0, maxValue);
        f->SetParameter(1, maxPosition);

        configureParameter(f, 2, "stdDevGauss", p["stdDevGauss"]);
        configureParameter(f, 3, "asymGauss", p["asymGauss"]);

        configureParameter(f, 4, "ampBack", p["ampBack"]);
        configureParameter(f, 5, "x1", p["x1"]);
        configureParameter(f, 6, "x2", p["x2"]);

        configureParameter(f, 7, "ampGaussDelta", p["ampGaussDelta"]);
        configureParameter(f, 8, "meanGaussDelta", p["meanGaussDelta"]);
        configureParameter(f, 9, "stdDevGaussDelta", p["stdDevGaussDelta"]);

        f->SetParLimits(0, 0.05 * maxValue, 1.1 * maxValue);
        f->SetParLimits(1, xMin, 1.05);
        f->SetParLimits(2, 0.0085, 0.085);
        f->SetParLimits(3, 0.33, 3.00);
        f->SetParLimits(7, 0.00, hist->GetBinContent(hist->FindBin(getDeltaPeak(hist, i, j, k, l))));
        f->SetParLimits(9, 0.01, 0.1);

        hist->Fit(f, "RQ");
        return f;
    }

    // Методы Get... теперь просто дергают общую функцию слияния
    double getDownEdge(TH1F* hist, int i, int j, int k, int l) const override {
        double downEdge =  getEffectiveParams(i, j, k, l).value("downEdge", 0.8);
        double x_min = hist->GetXaxis()->GetXmin();

        if (downEdge < x_min) downEdge = x_min;
        return downEdge;
    }

    double getUpEdge(TH1F* hist, int i, int j, int k, int l) const override {
        double upEdge = getEffectiveParams(i, j, k, l).value("upEdge", 1.1);
        double x_max = hist->GetXaxis()->GetXmax();

        if (upEdge > x_max) upEdge = x_max;
        return upEdge;
    }

    double getDeltaPeak(TH1F* hist, int i, int j, int k, int l) const override {
        json p = getEffectiveParams(i, j, k, l);
        if (p.contains("meanGaussDelta")) return p["meanGaussDelta"].value("val", 0.0);
        return 0.0;
    }
};