#pragma once

#include "iParameterProvider.h"
#include "../MM_project_const.h"

#include <nlohmann/json.hpp>
#include <fstream>

using json = nlohmann::json;

class JsonParameterProvider : public IParameterProvider {
private:
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
    JsonParameterProvider(const std::string& path)
        : IParameterProvider() {
        std::ifstream f(path);
        f >> config;
    }

    FitParameter getParameter(const std::string& name, int i, int j, int k, int l, TH1F* hist = nullptr) const override {
        // Получаем объединенный JSON объект для данных индексов (Q2, W, etc.)
        // Эта логика у тебя уже реализована в FitterEngine
        json p = getEffectiveParams(i, j, k, l);

        if (p.contains(name)) {
            auto& item = p[name];

            // Случай 1: В JSON просто число (например, "sigma": 0.05)
            if (item.is_number()) {
                return FitParameter(name, item.get<double>());
            }

            // Случай 2: В JSON объект (например, "sigma": {"val": 0.05, "fix": true, "limits": [0.01, 0.1]})
            if (item.is_object()) {
                double val = item.value("val", 0.0);
                bool fix  = item.value("fix", false);
                double low = -FitParameter::NO_LIMIT;
                double high = FitParameter::NO_LIMIT;

                if (item.contains("limits") && item["limits"].is_array() && item["limits"].size() == 2) {
                    auto& l_json = item["limits"];
                    
                    low  = l_json[0].is_null() ? -FitParameter::NO_LIMIT : l_json[0].get<double>();
                    high = l_json[1].is_null() ?  FitParameter::NO_LIMIT : l_json[1].get<double>();
                }

                return FitParameter(name, val, fix, low, high);
            }
        }

        // Если параметра нет, возвращаем дефолтный (val=0, fix=false, no limits)
        return FitParameter(name); 
    }

    bool hasParameter(const std::string& name, int i, int j, int k, int l) const override {
        return getEffectiveParams(i, j, k, l).contains(name);
    }

    void getBounds(int i, int j, int k, int l, double& low, double& high) const override {
        double defLow = 0.8, defHigh = 1.1;
        
        low = hasParameter("downEdge", i, j, k, l) 
            ? getParameter("downEdge", i, j, k, l).val 
            : defLow;
            
        high = hasParameter("upEdge", i, j, k, l) 
            ? getParameter("upEdge", i, j, k, l).val 
            : defHigh;
    }
};