#include "devices/bsim4_dc_core.hpp"

#include <iomanip>
#include <initializer_list>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

struct Point {
    double vg;
    double vd;
    double vs;
    double vb;
    double temp_c;
};

struct Deck {
    std::unordered_map<std::string, double> params;
    std::vector<Point> points;
    double width;
    double length;
};

std::unordered_map<std::string, double> makeParams(
    std::initializer_list<std::pair<const char*, double>> values) {
    std::unordered_map<std::string, double> out;
    for (const auto& [name, value] : values) out[name] = value;
    return out;
}

Deck deckA() {
    std::vector<Point> points;
    for (int i = 0; i <= 10; ++i) points.push_back({0.1 * i, 0.8, 0.0, 0.0, 27.0});
    return {makeParams({
        {"LEVEL", 54.0}, {"VTH0", 0.40}, {"U0", 500.0},
        {"TOXE", 1.0e-8}, {"XJ", 0.15e-6}, {"K1", 0.0}, {"K2", 0.0},
        {"DVT0", 0.0}, {"DVT1", 0.53}, {"DVT2", -0.032}, {"UA", 2.0e-9},
        {"UB", 5.0e-19}, {"VSAT", 1.0e5}
    }), points, 1.0e-6, 1.0e-6};
}

Deck deckB() {
    std::vector<Point> points;
    const double vgs[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2};
    const double vbs[] = {0.0, -0.3, 0.3};
    for (double vgsV : vgs)
        for (double vbsV : vbs)
            points.push_back({vgsV, 0.6, 0.0, vbsV, 27.0});
    return {makeParams({
        {"LEVEL", 54.0}, {"VTH0", 0.50}, {"TOXE", 5.0e-9}, {"TOXM", 5.0e-9},
        {"NDEP", 2.0e17}, {"NSD", 1.0e20}, {"K1", 0.6}, {"K2", -0.05},
        {"K3", 120.0}, {"K3B", 0.5}, {"W0", 0.6e-6}, {"LPE0", 1.2e-7},
        {"LPEB", 0.3}, {"DVT0", 1.8}, {"DVT1", 0.42}, {"DVT2", -0.02},
        {"DVT0W", 0.4}, {"DVT1W", 1.2e6}, {"DVT2W", -0.03}, {"DSUB", 0.35},
        {"DROUT", 0.45}, {"NFACTOR", 1.2}, {"CDSC", 1.2e-4}, {"CDSCB", 0.02},
        {"CDSCD", 0.005}, {"CIT", 1.0e-5}, {"ETA0", 0.15}, {"ETAB", -0.08},
        {"XJ", 1.2e-7}, {"VSAT", 1.1e5}, {"U0", 280.0}, {"A0", 0.9},
        {"AGS", 0.2}, {"B0", 5.0e-8}, {"B1", 1.0e-8}, {"KETA", 0.03},
        {"A1", 0.0}, {"A2", 1.0}, {"PCLM", 1.2}, {"PDIBL1", 0.2},
        {"PDIBL2", 0.003}, {"PVAG", 0.1}, {"PDIBLB", 0.02}, {"FPROUT", 0.05},
        {"RDSMOD", 0.0}, {"RDSW", 150.0}, {"WR", 1.0}, {"PRWG", 1.0},
        {"PRWB", 0.1}, {"DELTA", 0.01}, {"VOFF", -0.08}, {"MINV", 0.2},
        {"UA", 2.0e-9}, {"UB", 6.0e-19}
    }), points, 1.0e-6, 1.0e-6};
}

Deck deckC() {
    std::vector<Point> points;
    const double vds[] = {0.02, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5};
    for (double vdsV : vds) points.push_back({0.8, vdsV, 0.0, 0.0, 27.0});
    return {makeParams({
        {"LEVEL", 54.0}, {"VTH0", 0.40}, {"U0", 500.0},
        {"TOXE", 1.0e-8}, {"XJ", 0.15e-6}, {"K1", 0.0}, {"K2", 0.0},
        {"DVT0", 0.0}, {"DVT1", 0.53}, {"DVT2", -0.032}, {"UA", 2.0e-9},
        {"UB", 5.0e-19}, {"VSAT", 1.0e5}
    }), points, 5.0e-6, 0.3e-6};
}

Deck deckD() {
    std::vector<Point> points;
    const double vgs[] = {0.2, 0.5, 0.8, 1.1};
    const double temps[] = {0.0, 100.0};
    for (double vgV : vgs)
        for (double tempC : temps)
            points.push_back({vgV, 0.8, 0.0, 0.0, tempC});
    return {makeParams({
        {"LEVEL", 54.0}, {"VTH0", 0.40}, {"U0", 500.0},
        {"TOXE", 1.0e-8}, {"XJ", 0.15e-6}, {"K1", 0.0}, {"K2", 0.0},
        {"DVT0", 0.0}, {"DVT1", 0.53}, {"DVT2", -0.032}, {"UA", 2.0e-9},
        {"UB", 5.0e-19}, {"VSAT", 1.0e5}, {"AT", 3.3e-4}
    }), points, 1.0e-6, 1.0e-6};
}

} // namespace

int main(int argc, char** argv) {
    const std::string name = argc > 1 ? argv[1] : "a";
    Deck deck;
    if (name == "a") deck = deckA();
    else if (name == "b") deck = deckB();
    else if (name == "c") deck = deckC();
    else if (name == "d") deck = deckD();
    else { std::cerr << "unknown deck " << name << "\n"; return 2; }
    std::cout << std::scientific << std::setprecision(10);
    const bool internals = argc > 2 && std::string(argv[2]) == "--internals";
    for (const auto& point : deck.points) {
        const auto model = gspice::Bsim4ParameterSet::from(deck.params)
            .prepare(deck.width, deck.length, point.temp_c);
        const auto result =
            gspice::bsim4EvaluateDc(model, {point.vd, point.vg, point.vs, point.vb});
        if (!result.valid) { std::cerr << result.reason << "\n"; return 2; }
        if (internals) {
            const auto& s = result.state;
            std::cout << point.vg << " " << point.vb << " id=" << result.current[0]
                      << " vth=" << s.vth << " vgt=" << s.vgtEff
                      << " vdsat=" << s.vdsat << " abulk=" << s.abulk
                      << " mu=" << s.mobility << " beta=" << s.effectiveBeta
                      << " idl=" << s.idl << " clm=" << s.clmBoost
                      << " vadibl=" << s.vadibl << "\n";
            continue;
        }
        std::cout << result.current[0] << "\n";
    }
}