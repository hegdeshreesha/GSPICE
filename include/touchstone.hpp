#ifndef GSPICE_TOUCHSTONE_HPP
#define GSPICE_TOUCHSTONE_HPP

#include <string>
#include <vector>
#include <complex>
#include <fstream>
#include <sstream>
#include <iostream>
#include "matrix.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace gspice {

class Touchstone {
public:
    struct FreqPoint {
        double freq;
        MatrixComplex S;
    };

    Touchstone(const std::string& filename, int num_ports) : num_ports_(num_ports) {
        load(filename);
    }

    MatrixComplex getS(double freq) const {
        if (data_.empty()) return MatrixComplex(num_ports_);
        
        // Find interpolation bounds
        if (freq <= data_.front().freq) return data_.front().S;
        if (freq >= data_.back().freq) return data_.back().S;

        for (size_t i = 0; i < data_.size() - 1; ++i) {
            if (freq >= data_[i].freq && freq <= data_[i+1].freq) {
                double f1 = data_[i].freq;
                double f2 = data_[i+1].freq;
                double t = (freq - f1) / (f2 - f1);
                
                MatrixComplex S_interp(num_ports_);
                for (int r = 0; r < num_ports_; ++r) {
                    for (int c = 0; c < num_ports_; ++c) {
                        // Linear interpolation of complex values
                        S_interp(r, c) = data_[i].S(r, c) + t * (data_[i+1].S(r, c) - data_[i].S(r, c));
                    }
                }
                return S_interp;
            }
        }
        return data_.front().S;
    }

    double getZ0() const { return z0_; }

private:
    int num_ports_;
    double z0_ = 50.0;
    std::vector<FreqPoint> data_;

    void load(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Warning: Could not open Touchstone file " << filename << std::endl;
            return;
        }

        std::string line;
        bool format_RI = false; 
        bool format_MA = true;
        double freq_mult = 1.0;

        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '!') continue;
            if (line[0] == '#') {
                std::stringstream ss(line.substr(1));
                std::string token;
                while (ss >> token) {
                    for (auto & c: token) c = toupper(c);
                    if (token == "HZ") freq_mult = 1.0;
                    else if (token == "KHZ") freq_mult = 1e3;
                    else if (token == "MHZ") freq_mult = 1e6;
                    else if (token == "GHZ") freq_mult = 1e9;
                    else if (token == "MA") { format_MA = true; format_RI = false; }
                    else if (token == "RI") { format_RI = true; format_MA = false; }
                    else if (token == "DB") { format_MA = false; format_RI = false; }
                }
                continue;
            }

            std::stringstream ss(line);
            double f;
            if (!(ss >> f)) continue;

            FreqPoint pt;
            pt.freq = f * freq_mult;
            pt.S = MatrixComplex(num_ports_);

            std::vector<double> vals;
            double v;
            while (ss >> v) vals.push_back(v);

            int expected_vals = num_ports_ * num_ports_ * 2;
            while (vals.size() < expected_vals && std::getline(file, line)) {
                if (line.empty() || line[0] == '!') continue;
                std::stringstream ss2(line);
                while (ss2 >> v) vals.push_back(v);
            }

            if (vals.size() >= expected_vals) {
                int idx = 0;
                if (num_ports_ == 2) {
                    // 2-port standard: S11, S21, S12, S22
                    double v1 = vals[0], v2 = vals[1], v3 = vals[2], v4 = vals[3];
                    double v5 = vals[4], v6 = vals[5], v7 = vals[6], v8 = vals[7];
                    
                    if (format_RI) {
                        pt.S(0,0) = {v1, v2}; pt.S(1,0) = {v3, v4};
                        pt.S(0,1) = {v5, v6}; pt.S(1,1) = {v7, v8};
                    } else if (format_MA) {
                        pt.S(0,0) = std::polar(v1, v2*M_PI/180.0); pt.S(1,0) = std::polar(v3, v4*M_PI/180.0);
                        pt.S(0,1) = std::polar(v5, v6*M_PI/180.0); pt.S(1,1) = std::polar(v7, v8*M_PI/180.0);
                    } else { // DB
                        pt.S(0,0) = std::polar(std::pow(10.0, v1/20.0), v2*M_PI/180.0);
                        pt.S(1,0) = std::polar(std::pow(10.0, v3/20.0), v4*M_PI/180.0);
                        pt.S(0,1) = std::polar(std::pow(10.0, v5/20.0), v6*M_PI/180.0);
                        pt.S(1,1) = std::polar(std::pow(10.0, v7/20.0), v8*M_PI/180.0);
                    }
                } else {
                    for (int r = 0; r < num_ports_; ++r) {
                        for (int c = 0; c < num_ports_; ++c) {
                            double v1 = vals[idx++];
                            double v2 = vals[idx++];
                            if (format_RI) pt.S(r,c) = {v1, v2};
                            else if (format_MA) pt.S(r,c) = std::polar(v1, v2*M_PI/180.0);
                            else pt.S(r,c) = std::polar(std::pow(10.0, v1/20.0), v2*M_PI/180.0); // DB
                        }
                    }
                }
                data_.push_back(pt);
            }
        }
        std::cout << "Loaded Touchstone file " << filename << " with " << data_.size() << " frequency points." << std::endl;
    }
};

} // namespace gspice
#endif // GSPICE_TOUCHSTONE_HPP
