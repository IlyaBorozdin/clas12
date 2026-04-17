#pragma once

#include "fitterEngine.h"
#include "TMath.h"

class LegacyNeutronFitter : public FitterEngine {
public:
    using FitterEngine::FitterEngine;

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

    double calculateYield(TF1* f, double binWidth) const override {
        double amp = f->GetParameter(0);
        double sigma = std::abs(f->GetParameter(2));
        double asym = std::abs(f->GetParameter(3));
        
        // Твоя текущая формула
        return std::sqrt(TMath::Pi() / 2.0) * amp * sigma * (1.0 + asym) / binWidth;
    }

    double calculateYieldError(TF1* f, double binWidth) const override {
        double yieldVal = calculateYield(f, binWidth);
        if (yieldVal <= 0) return 1.15 / binWidth; // Наша страховка для пустых бинов

        double yieldErr;
        double amp = f->GetParameter(0);
        double sigma = std::abs(f->GetParameter(2));
        double asym = std::abs(f->GetParameter(3));
        double eAmpGauss = f->GetParError(0);
        double eStdDevGauss = f->GetParError(2);
        double eAsymGauss = f->GetParError(3);
        
        if (amp != 0 && sigma != 0) {
                yieldErr = yieldVal * std::sqrt(
                    std::pow(eAmpGauss / amp, 2) +
                    std::pow(eStdDevGauss / sigma, 2) +
                    std::pow(eAsymGauss / (1.0 + std::abs(asym)), 2)
                );
            } else {
                yieldErr = 1.15 / binWidth; 
            }
        
        return yieldErr;
    }
};

class EMGNeutronFitter : public FitterEngine {
public:
    using FitterEngine::FitterEngine;

    static double evaluateFullModel(double* x, double* p) {
        double val = x[0];
        double A     = p[0]; // Интеграл
        double mu    = p[1];
        double sigma = p[2];
        double tau   = p[3]; // Может быть и положительным, и отрицательным

        double emg = 0.0;
        if (std::abs(tau) < 1e-7) {
            // Предел EMG при tau -> 0. 
            // Убираем std::abs, чтобы сохранить знак хвоста и непрерывность производной.
            double arg = (val - mu) / sigma;
            double gauss = (A / (sigma * std::sqrt(2.0 * TMath::Pi()))) * std::exp(-0.5 * arg * arg);
            
            // Асимптотическая поправка: 1 / (1 + z)
            // double correction = 1.0 + arg * tau / sigma;
            
            // Защита от деления на 0 (хотя при tau < 1e-7 это почти невозможно)
            // emg = (std::abs(correction) > 1e-9) ? gauss / correction : gauss;
            emg = gauss;
        } else {
            // Используем стабильную форму с erfcx
            double absTau = std::abs(tau);
            double argErfcx = (1.0 / std::sqrt(2.0)) * (sigma / absTau + (val - mu) * (tau > 0 ? 1.0 : -1.0) / sigma);
            
            double gaussPart = std::exp(-0.5 * std::pow((val - mu) / sigma, 2));
            emg = (A / (2.0 * absTau)) * gaussPart * my_erfcx(argErfcx);
        }

        // Фон и Дельта остаются как были
        double back = TMath::Max(0., p[4] * (val - p[5]) * (val - p[6]));
        double argDelta = (p[9] != 0) ? (val - p[8]) / p[9] : 0;
        double delta = p[7] * TMath::Exp(-0.5 * argDelta * argDelta);

        double res = emg + back + delta;
        return (std::isnan(res) || std::isinf(res)) ? 1e30 : res; // Возвращаем очень большое число для минимизатора
    }

    TF1* fit(TH1F* hist, int i, int j, int k, int l) const override {
        json p = getEffectiveParams(i, j, k, l);
        
        double xMin = getDownEdge(hist, i, j, k, l);
        double xMax = getUpEdge(hist, i, j, k, l);

        double maxPosition, maxValue, integral;
        getNeutronPeak(hist, i, j, k, l, maxPosition, maxValue);
        getNeutronIntegral(hist, i, j, k, l, integral);

        // Создаем TF1, передавая указатель на метод текущего класса
        // 10 - количество параметров
        TF1* f = new TF1("emg_fit", &EMGNeutronFitter::evaluateFullModel, xMin, xMax, 10);

        // Присваиваем имена параметрам (очень помогает при отладке в TBrowser)
        f->SetParNames("h_EMG", "mu_EMG", "sigma_EMG", "tau_EMG", "amp_BG", "x1_BG", "x2_BG", "h_Delta", "mu_Delta", "sigma_Delta");

        // --- Инициализация параметров из JSON ---
        f->SetParameter(0, 0.8 * integral * hist->GetXaxis()->GetBinWidth(1));      // h_EMG
        f->SetParameter(1, maxPosition);   // mu_EMG
        
        configureParameter(f, 2, "stdDevGauss", p["stdDevGauss"]);
        f->SetParameter(3, 0.01);
        
        configureParameter(f, 4, "ampBack", p["ampBack"]);
        configureParameter(f, 5, "x1", p["x1"]);
        configureParameter(f, 6, "x2", p["x2"]);

        configureParameter(f, 7, "ampGaussDelta", p["ampGaussDelta"]);
        configureParameter(f, 8, "meanGaussDelta", p["meanGaussDelta"]);
        configureParameter(f, 9, "stdDevGaussDelta", p["stdDevGaussDelta"]);

        // --- Установка лимитов (аналогично Legacy) ---
        f->SetParLimits(0, 0.0, 1.8 * integral * hist->GetXaxis()->GetBinWidth(1));
        f->SetParLimits(1, xMin, 1.05);
        f->SetParLimits(2, 0.005, 0.1); 
        f->SetParLimits(3, -0.5, 0.5); // Лимит на тау (хвост)

        // Для дельты используем твой метод getDeltaPeak
        double deltaPos = getDeltaPeak(hist, i, j, k, l);
        f->SetParLimits(7, 0.00, hist->GetBinContent(hist->FindBin(getDeltaPeak(hist, i, j, k, l))));
        f->SetParLimits(9, 0.01, 0.1);

        hist->Fit(f, "RQ");
        return f;
    }

    double calculateYield(TF1* f, double binWidth) const override {
        return f->GetParameter(0) / binWidth; 
    }

    double calculateYieldError(TF1* f, double binWidth) const override {
        return (f->GetParameter(0) <= 0.0 ? 1.15 : f->GetParError(0))  / binWidth;
    }

private:
    static double my_erfcx (double x) {
        double a, d, e, m, p, q, r, s, t;

        a = fmax (x, 0.0 - x); // NaN preserving absolute value computation

        /* Compute q = (a-4)/(a+4) accurately. [0,INF) -> [-1,1] */
        m = a - 4.0;
        p = a + 4.0;
        r = 1.0 / p;
        q = m * r;
        t = fma (q + 1.0, -4.0, a); 
        e = fma (q, -a, t); 
        q = fma (r, e, q); 

        /* Approximate (1+2*a)*exp(a*a)*erfc(a) as p(q)+1 for q in [-1,1] */
        s = q * q;
        p =             0x1.edcad78fc8044p-31;  //  8.9820305531190140e-10
        t =             0x1.b1548f14735d1p-30;  //  1.5764464777959401e-09
        p = fma (p, s, -0x1.a1ad2e6c4a7a8p-27); // -1.2155985739342269e-08
        t = fma (t, s, -0x1.1985b48f08574p-26); // -1.6386753783877791e-08
        p = fma (p, s,  0x1.c6a8093ac4f83p-24); //  1.0585794011876720e-07
        t = fma (t, s,  0x1.31c2b2b44b731p-24); //  7.1190423171700940e-08
        p = fma (p, s, -0x1.b87373facb29fp-21); // -8.2040389712752056e-07
        t = fma (t, s,  0x1.3fef1358803b7p-22); //  2.9796165315625938e-07
        p = fma (p, s,  0x1.7eec072bb0be3p-18); //  5.7059822144459833e-06
        t = fma (t, s, -0x1.78a680a741c4ap-17); // -1.1225056665965572e-05
        p = fma (p, s, -0x1.9951f39295cf4p-16); // -2.4397380523258482e-05
        t = fma (t, s,  0x1.3be1255ce180bp-13); //  1.5062307184282616e-04
        p = fma (p, s, -0x1.a1df71176b791p-13); // -1.9925728768782324e-04
        t = fma (t, s, -0x1.8d4aaa0099bc8p-11); // -7.5777369791018515e-04
        p = fma (p, s,  0x1.49c673066c831p-8);  //  5.0319701025945277e-03
        t = fma (t, s, -0x1.0962386ea02b7p-6);  // -1.6197733983519948e-02
        p = fma (p, s,  0x1.3079edf465cc3p-5);  //  3.7167515521269866e-02
        t = fma (t, s, -0x1.0fb06dfedc4ccp-4);  // -6.6330365820039094e-02
        p = fma (p, s,  0x1.7fee004e266dfp-4);  //  9.3732834999538536e-02
        t = fma (t, s, -0x1.9ddb23c3e14d2p-4);  // -1.0103906603588378e-01
        p = fma (p, s,  0x1.16ecefcfa4865p-4);  //  6.8097054254651804e-02
        t = fma (t, s,  0x1.f7f5df66fc349p-7);  //  1.5379652102610957e-02
        p = fma (p, q, t);
        p = fma (p, q, -0x1.1df1ad154a27fp-3);  // -1.3962111684056208e-01
        p = fma (p, q,  0x1.dd2c8b74febf6p-3);  //  2.3299511862555250e-01

        /* Divide (1+p) by (1+2*a) ==> exp(a*a)*erfc(a) */
        d = a + 0.5;
        r = 1.0 / d;
        r = r * 0.5;
        q = fma (p, r, r); // q = (p+1)/(1+2*a)
        t = q + q;
        e = (p - q) + fma (t, -a, 1.0); // residual: (p+1)-q*(1+2*a)
        r = fma (e, r, q);

        /* Handle argument of infinity */
        if (a > 0x1.fffffffffffffp1023) r = 0.0;

        /* Handle negative arguments: erfcx(x) = 2*exp(x*x) - erfcx(|x|) */
        if (x < 0.0) {
            s = x * x;
            d = fma (x, x, -s);
            e = exp (s);
            r = e - r;
            r = fma (e, d + d, r); 
            r = r + e;
            if (e > 0x1.fffffffffffffp1023) r = e; // avoid creating NaN
        }
        return r;
    }
};