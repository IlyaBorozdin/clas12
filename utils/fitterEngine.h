#pragma once

#include <fstream>
#include <nlohmann/json.hpp>
#include "TH1F.h"
#include "TF1.h"

#include "fittingStrategy.h"

using json = nlohmann::json;

class FitterEngine : public FittingStrategy {
protected:
    json config;

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
    FitterEngine(const std::string& configPath) {
        std::ifstream f(configPath);
        f >> config;
    }

    // Эти методы общие для всех — границы бина зависят от кинематики, а не от функции
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
        return 1.232;
    }
    
    // Чисто виртуальный метод — конкретика будет в наследниках
    virtual TF1* fit(TH1F* hist, int i, int j, int k, int l) const override = 0;
};