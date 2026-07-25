#ifndef GSPICE_RESULT_DATABASE_HPP
#define GSPICE_RESULT_DATABASE_HPP

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <iomanip>

namespace gspice {

struct AnalysisResult {
    std::string analysisType;
    std::vector<std::string> variableNames;
    std::vector<double> timeOrFreq;
    std::vector<std::vector<double>> values; // values[sample_index][var_index]
};

class ResultDatabase {
public:
    void clear() {
        results_.clear();
    }

    void setVariableNames(const std::vector<std::string>& names) {
        currentResult_.variableNames = names;
    }

    void setAnalysisType(const std::string& type) {
        currentResult_.analysisType = type;
    }

    void addSample(double timeOrFreq, const std::vector<double>& nodeValues) {
        currentResult_.timeOrFreq.push_back(timeOrFreq);
        currentResult_.values.push_back(nodeValues);
    }

    void finalizeResult() {
        results_.push_back(currentResult_);
        currentResult_ = AnalysisResult();
    }

    const std::vector<AnalysisResult>& results() const { return results_; }

    bool exportCsv(const std::string& filename) const {
        if (results_.empty()) return false;
        std::ofstream out(filename);
        if (!out.is_open()) return false;

        const auto& res = results_.back();
        out << "time";
        for (const auto& var : res.variableNames) {
            out << "," << var;
        }
        out << "\n";

        for (size_t i = 0; i < res.timeOrFreq.size(); ++i) {
            out << std::scientific << std::setprecision(9) << res.timeOrFreq[i];
            for (size_t j = 0; j < res.values[i].size(); ++j) {
                out << "," << res.values[i][j];
            }
            out << "\n";
        }
        return true;
    }

private:
    AnalysisResult currentResult_;
    std::vector<AnalysisResult> results_;
};

} // namespace gspice

#endif // GSPICE_RESULT_DATABASE_HPP
