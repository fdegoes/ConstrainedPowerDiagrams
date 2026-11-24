// Runs a 2D example of the original Power Particles code (Before2015/pp,
// data/2D) with one of the integrators of core/integrator.h, fluid only
// (the fluid fills the domain: no air or solid particles), and writes
// the state and its energy and momentum per frame.
//
//   runPP2D <domain.line> <particles.pp|.dat> <param.txt> <out_dir>
//           [--integrator powerparticles|constrained|softened|gm]
//           [--frames N] [--every K] [--alpha A] [--stiffness k]
//           [--maxiter N] [--volumetol T] [--lloydtol T] [--lloyditers N] [--log]
//           [--no-csv] [--no-energy]
//
// Inputs (the formats of the 2015 code, pp2/SimIO.cpp, ParticlesIO.cpp,
// Polygon.cpp):
//   .line   "v x y" vertices and "e i j" edges (1-based) of one closed polygon
//   .pp     <Fluid> ... </End> blocks of <Particle> records with <pos x y />,
//           <vel vx vy />, <mass m />, <pres p />, <weight w />,
//           <targetVol V />; <Air> and <Solid> blocks are rejected
//   .dat    a header value, then "x y vx vy w p density" per particle; the
//           target volumes are the cell volumes of the loaded weights
//   param   "key: value" lines; used: dt, cfl, gravity (pointing down -y),
//           cap_tol, cap_max_iters. viscosity (see notes/viscosity.md),
//           buoyancy, stiffness, reseed and solids are not supported
//           (warning if nonzero).
//
// The frames follow the 2015 code's Sim::run: frame f is the state after f
// advances of dt, each split into CFL substeps (cfl * r_min / max |v_k|,
// ComputeCFLStep). Outputs in out_dir:
//   frame_NNNN.csv  one row per particle: x,y,vx,vy,ux,uy,w,p,volume,cx,cy,
//                   cellVolume (u and w of our state: generators q + u,
//                   power weights w + |u|^2); skipped with --no-csv
//   energy.csv      one FlowMeasures row per frame (time, energies, momenta)
//                   plus the substeps, Newton iterations and errors; skipped
//                   with --no-energy
//   run.txt         the parameters of the run (always written)

#include "core/domain.h"
#include "core/integrator.h"
#include "core/shiftedDiagram.h"
#include "core/triangulate2D.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace
{

using PP::Integrator;
using PP::ParticleState;

struct Params
{
    double _dt = 0.01;
    double _cfl = 1.0;
    double _gravity = 0.0;
    double _capTolerance = 0.0;
    int    _capIterations = 2;
    std::map<std::string, std::string> _ignored;  // unsupported keys with nonzero values
};

bool
_ReadLine(const std::string& path, std::vector<double>& loop, std::string& err)
{
    std::ifstream in(path);
    if (not in) {
        err = "cannot open " + path;
        return false;
    }
    std::vector<double> pts;
    std::vector<std::pair<int, int>> edges;
    std::string tag;
    while (in >> tag) {
        if (tag == "v") {
            double x, y;
            in >> x >> y;
            pts.push_back(x);
            pts.push_back(y);
        } else if (tag == "e") {
            int a, b;
            in >> a >> b;
            edges.emplace_back(a - 1, b - 1);
        } else {
            std::string rest;
            std::getline(in, rest);
        }
    }
    if (edges.size() < 3) {
        err = path + ": fewer than 3 edges";
        return false;
    }
    // Chain the edges into one loop, from the first edge.
    std::map<int, int> next;
    for (const auto& e : edges) {
        next[e.first] = e.second;
    }
    int v = edges[0].first;
    for (size_t k = 0; k < edges.size(); ++k) {
        loop.push_back(pts[2 * size_t(v)]);
        loop.push_back(pts[2 * size_t(v) + 1]);
        const auto it = next.find(v);
        if (it == next.end()) {
            err = path + ": the edges do not form one closed loop";
            return false;
        }
        v = it->second;
    }
    if (v != edges[0].first) {
        err = path + ": the edges do not form one closed loop";
        return false;
    }
    return true;
}

// Values of "<tag a b />".
bool
_ReadTag(const std::string& line, const char* tag, std::vector<double>& values)
{
    const std::string open = std::string("<") + tag + " ";
    if (line.compare(0, open.size(), open) != 0) {
        return false;
    }
    std::istringstream s(line.substr(open.size()));
    values.clear();
    double x;
    while (s >> x) {
        values.push_back(x);
    }
    return true;
}

bool
_ReadPP(const std::string& path, ParticleState& state, std::vector<double>& volumes, std::string& err)
{
    std::ifstream in(path);
    if (not in) {
        err = "cannot open " + path;
        return false;
    }
    std::string line;
    bool fluid = false;
    std::vector<double> v, masses;
    while (std::getline(in, line)) {
        if (line == "<Fluid>") {
            fluid = true;
        } else if (line == "<Air>" or line == "<Solid>") {
            err = path + ": " + line + " particles are not supported (the fluid must fill the domain)";
            return false;
        } else if (line == "</End>") {
            fluid = false;
        } else if (not fluid) {
            continue;
        } else if (_ReadTag(line, "pos", v)) {
            state._positions.insert(state._positions.end(), v.begin(), v.end());
        } else if (_ReadTag(line, "vel", v)) {
            state._velocities.insert(state._velocities.end(), v.begin(), v.end());
        } else if (_ReadTag(line, "pres", v)) {
            state._pressures.push_back(v.at(0));
        } else if (_ReadTag(line, "weight", v)) {
            state._weights.push_back(v.at(0));
        } else if (_ReadTag(line, "targetVol", v)) {
            volumes.push_back(v.at(0));
        } else if (_ReadTag(line, "mass", v)) {
            masses.push_back(v.at(0));
        }
    }
    const size_t n = volumes.size();
    if (n == 0 or state._positions.size() != 2 * n or state._velocities.size() != 2 * n or
        state._weights.size() != n) {
        err = path + ": inconsistent particle records";
        return false;
    }
    // Our integrators have unit density (mass = volume), as all the
    // fluid-only 2D examples.
    for (size_t i = 0; i < n and i < masses.size(); ++i) {
        if (std::fabs(masses[i] - volumes[i]) > 1e-6 * volumes[i]) {
            err = path + ": particle " + std::to_string(i) + " has density " + std::to_string(masses[i] / volumes[i]) +
                  " (only unit density is supported)";
            return false;
        }
    }
    state._shifts.assign(2 * n, 0.0);
    return true;
}

bool
_ReadDat(const std::string& path, ParticleState& state, std::string& err)
{
    std::ifstream in(path);
    if (not in) {
        err = "cannot open " + path;
        return false;
    }
    double header;
    in >> header;
    double x, y, vx, vy, w, p, d;
    while (in >> x >> y >> vx >> vy >> w >> p >> d) {
        state._positions.insert(state._positions.end(), {x, y});
        state._velocities.insert(state._velocities.end(), {vx, vy});
        state._weights.push_back(w);
        state._pressures.push_back(p);
    }
    if (state._weights.empty()) {
        err = path + ": no particles";
        return false;
    }
    state._shifts.assign(state._positions.size(), 0.0);
    return true;
}

bool
_ReadParams(const std::string& path, Params& p, std::string& err)
{
    std::ifstream in(path);
    if (not in) {
        err = "cannot open " + path;
        return false;
    }
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream s(line);
        std::string key, value;
        if (not (s >> key >> value)) {
            continue;
        }
        if (key == "dt:") {
            p._dt = std::stod(value);
        } else if (key == "cfl:") {
            p._cfl = std::stod(value);
        } else if (key == "gravity:") {
            p._gravity = std::stod(value);
        } else if (key == "cap_tol:") {
            p._capTolerance = std::stod(value);
        } else if (key == "cap_max_iters:") {
            p._capIterations = std::stoi(value);
        } else if (key == "viscosity:" or key == "buoyancy:" or key == "stiffness:" or key == "diffusion:" or
                   key == "source_rate:") {
            if (std::stod(value) != 0.0) {
                p._ignored[key] = value;
            }
        } else if (key == "reseed:") {
            if (std::stoi(value) >= 0) {
                p._ignored[key] = value;
            }
        }
    }
    return true;
}

void
_WriteFrame(const std::string& path, const ParticleState& s, const std::vector<double>& volumes,
            const Integrator::Result& r)
{
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (not f) {
        return;
    }
    std::fprintf(f, "x,y,vx,vy,ux,uy,w,p,volume,cx,cy,cellVolume\n");
    for (size_t i = 0; i < volumes.size(); ++i) {
        const double cx = r._centroids.empty() ? NAN : r._centroids[2 * i];
        const double cy = r._centroids.empty() ? NAN : r._centroids[2 * i + 1];
        const double cv = r._cellVolumes.empty() ? NAN : r._cellVolumes[i];
        std::fprintf(f, "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                     s._positions[2 * i], s._positions[2 * i + 1], s._velocities[2 * i], s._velocities[2 * i + 1],
                     s._shifts[2 * i], s._shifts[2 * i + 1], s._weights[i],
                     s._pressures.empty() ? 0.0 : s._pressures[i], volumes[i], cx, cy, cv);
    }
    std::fclose(f);
}

int
_Usage()
{
    std::cerr << "usage: runPP2D <domain.line> <particles.pp|.dat> <param.txt> <out_dir> "
                 "[--integrator powerparticles|constrained|softened|gm] [--frames N] [--every K] "
                 "[--alpha A] [--stiffness k] [--maxiter N] [--volumetol T] [--lloydtol T] [--lloyditers N] "
                 "[--velocity-projection none|initial|every] [--no-csv] [--no-energy] [--log]"
              << std::endl;
    return 2;
}

} // namespace

int
main(int argc, char** argv)
{
    if (argc < 5) {
        return _Usage();
    }
    const std::string linePath = argv[1], particlesPath = argv[2], paramPath = argv[3], out = argv[4];
    std::string kindName = "powerparticles";
    int frames = 100, every = 1;
    bool writeFrames = true, writeEnergy = true;
    PP::IntegratorOptions options;
    Integrator::Parameters params;
    for (int k = 5; k < argc; ++k) {
        const std::string a = argv[k];
        if (a == "--no-csv") {
            writeFrames = false;
            continue;
        }
        if (a == "--no-energy") {
            writeEnergy = false;
            continue;
        }
        if (a == "--log") {
            params._solver._log = &std::cerr;  // Newton iterations of every solve
            continue;
        }
        if (k + 1 >= argc) {
            return _Usage();
        }
        const std::string v = argv[++k];
        if (a == "--integrator") {
            kindName = v;
        } else if (a == "--frames") {
            frames = std::stoi(v);
        } else if (a == "--every") {
            every = std::max(1, std::stoi(v));
        } else if (a == "--alpha") {
            options._alpha = std::stod(v);
        } else if (a == "--stiffness") {
            options._stiffness = std::stod(v);
        } else if (a == "--maxiter") {
            params._solver._maxIterations = std::stoi(v);
        } else if (a == "--lloydtol") {
            params._lloydTolerance = std::stod(v);
        } else if (a == "--lloyditers") {
            params._lloydMaxIterations = std::stoi(v);
        } else if (a == "--velocity-projection") {
            // Constrained: none (SHAKE), initial (the initial velocities only), every (RATTLE).
            const std::map<std::string, PP::IntegratorOptions::VelocityProjection> modes = {
                {"none", PP::IntegratorOptions::VelocityProjection::None},
                {"initial", PP::IntegratorOptions::VelocityProjection::Initial},
                {"every", PP::IntegratorOptions::VelocityProjection::EveryStep},
            };
            if (not modes.count(v)) {
                return _Usage();
            }
            options._velocityProjection = modes.at(v);
        } else if (a == "--volumetol") {
            params._solver._volumeTolerance = params._solver._shiftTolerance = std::stod(v);
        } else {
            return _Usage();
        }
    }
    const std::map<std::string, PP::IntegratorKind> kinds = {
        {"powerparticles", PP::IntegratorKind::PowerParticles},
        {"constrained", PP::IntegratorKind::Constrained},
        {"softened", PP::IntegratorKind::Softened},
        {"gm", PP::IntegratorKind::GM},
    };
    if (not kinds.count(kindName)) {
        return _Usage();
    }

    std::string err;
    Params pp;
    std::vector<double> loop;
    if (not _ReadParams(paramPath, pp, err) or not _ReadLine(linePath, loop, err)) {
        std::cerr << err << std::endl;
        return 2;
    }
    for (const auto& kv : pp._ignored) {
        std::cerr << "warning: " << kv.first << " " << kv.second << " is not supported and ignored" << std::endl;
    }
    std::vector<double> triPts;
    std::vector<int> tris;
    PP::Domain domain;
    if (not PP::Triangulate2D(loop, {int(loop.size() / 2)}, triPts, tris, err) or
        not PP::Domain::Make2D(triPts, tris, domain, err)) {
        std::cerr << err << std::endl;
        return 2;
    }

    ParticleState state;
    std::vector<double> volumes;
    const bool dat = particlesPath.size() > 4 and particlesPath.substr(particlesPath.size() - 4) == ".dat";
    if (dat ? not _ReadDat(particlesPath, state, err) : not _ReadPP(particlesPath, state, volumes, err)) {
        std::cerr << err << std::endl;
        return 2;
    }
    if (dat) {
        // As the 2015 code: target volumes are the cells of the loaded weights.
        PP::ShiftedDiagram d;
        if (not d.Build(domain, state._positions, state._positions, state._shifts, state._weights, err)) {
            std::cerr << err << std::endl;
            return 2;
        }
        volumes = d.GetCells()._volumes;
    }
    // The targets must sum to the domain area (the fluid fills it): the
    // stored volumes are rounded (e.g. sum 1 + 1.2e-8), which the volume
    // solves cannot reach to 1e-9, so they are rescaled.
    double sum = 0.0;
    for (double V : volumes) {
        sum += V;
    }
    if (std::fabs(sum - domain.GetVolume()) > 1e-3 * domain.GetVolume()) {
        std::cerr << "warning: target volumes sum to " << sum << ", domain area " << domain.GetVolume() << std::endl;
    }
    for (double& V : volumes) {
        V *= domain.GetVolume() / sum;
    }

    params._gravity = {0.0, -pp._gravity};
    options._weightIterations = pp._capIterations;
    options._weightTolerance = pp._capTolerance;
    const auto integrator = PP::MakeIntegrator(kinds.at(kindName), domain, volumes, params, options);

    if (std::system(("mkdir -p '" + out + "'").c_str()) != 0) {
        std::cerr << "cannot create " << out << std::endl;
        return 2;
    }
    {
        std::ofstream run(out + "/run.txt");
        run << "domain " << linePath << "\nparticles " << particlesPath << "\nparams " << paramPath
            << "\nintegrator " << integrator->GetName() << "\nparticles_count " << volumes.size() << "\ndt " << pp._dt
            << "\ncfl " << pp._cfl << "\ngravity " << pp._gravity << "\ncap_tol " << pp._capTolerance
            << "\ncap_max_iters " << pp._capIterations << "\nalpha " << options._alpha << "\nstiffness "
            << options._stiffness << "\nframes " << frames << "\nevery " << every << "\nframe_csv "
            << (writeFrames ? "yes" : "no") << "\nenergy_csv " << (writeEnergy ? "yes" : "no") << "\n";
    }

    // Our integrators other than Power Particles start from their own
    // consistent state (the loaded positions are kept as the initial guess).
    Integrator::Result r;
    if (kinds.at(kindName) != PP::IntegratorKind::PowerParticles) {
        r = integrator->Initialize(state);
        if (not r.HasConverged()) {
            std::cerr << "Initialize failed: " << r._error << std::endl;
            return 1;
        }
    } else if (not integrator->Measure(state, r._measures, err)) {
        std::cerr << err << std::endl;
        return 1;
    }
    {
        // The cells of the initial state, for frame 0.
        PP::ShiftedDiagram d;
        if (d.Build(domain, state._positions, state._positions, state._shifts, state._weights, err)) {
            r._centroids = d.GetCells()._centroids;
            r._cellVolumes = d.GetCells()._volumes;
        }
    }

    std::FILE* energy = writeEnergy ? std::fopen((out + "/energy.csv").c_str(), "w") : nullptr;
    if (energy) {
        std::fprintf(energy, "%s,frame,substeps,iterations,maxVolumeError,maxCentroidGap,projected\n",
                     PP::FlowMeasures::GetCSVHeader(2).c_str());
    }
    auto log = [&](int frame, double t, int substeps, int iterations, const Integrator::Result& res,
                   const PP::FlowMeasures& measures) {
        if (not energy) {
            return;
        }
        std::fprintf(energy, "%s,%d,%d,%d,%.6g,%.6g,%d\n", measures.GetCSVRow(t).c_str(), frame, substeps,
                     iterations, res._maxVolumeError, res._maxCentroidGap, res._numProjected);
        std::fflush(energy);
    };
    // Our integrators measure the energy of a state at the start of the
    // step from it (Integrator::Result::_startMeasures), so that its kinetic
    // and potential energies are of the same time level: the row of a frame
    // is written at the first step of the next frame, and the last frame is
    // measured with one more solve (Integrator::Measure()).
    const bool lagged = kinds.at(kindName) != PP::IntegratorKind::PowerParticles;
    struct Pending
    {
        bool               _active = false;
        int                _frame = 0, _substeps = 0, _iterations = 0;
        double             _t = 0.0;
        Integrator::Result _res;
    } pending;
    char name[64];
    auto writeFrame = [&](int frame) {
        if (writeFrames) {
            std::snprintf(name, sizeof(name), "/frame_%04d.csv", frame);
            _WriteFrame(out + name, state, volumes, r);
        }
    };
    writeFrame(0);
    log(0, 0.0, 0, 0, r, r._measures);

    double t = 0.0;
    for (int frame = 1; frame <= frames; ++frame) {
        int substeps = 0, iterations = 0, projected = 0;
        for (int e = 0; e < every; ++e) {
            // One advance of dt in CFL substeps, as Sim::advance.
            double done = 0.0;
            while (done < pp._dt) {
                double h = PP::ComputeCFLStep(state._velocities, volumes, 2, pp._cfl);
                if (h == 0.0) {
                    h = pp._dt;
                }
                h = std::min(h, pp._dt - done);
                r = integrator->Step(state, h);
                if (r.HasConverged() and pending._active) {
                    log(pending._frame, pending._t, pending._substeps, pending._iterations, pending._res,
                        r._startMeasures);
                    pending._active = false;
                }
                if (not r.HasConverged()) {
                    std::cerr << "frame " << frame << ": step failed (" << r._status << ") " << r._error << std::endl;
                    if (energy) {
                        std::fclose(energy);
                    }
                    return 1;
                }
                done += h;
                ++substeps;
                iterations += r._iterations;
                projected += r._numProjected;
            }
            t += pp._dt;
        }
        r._numProjected = projected;
        writeFrame(frame);
        if (lagged) {
            pending = {true, frame, substeps, iterations, t, r};
        } else {
            log(frame, t, substeps, iterations, r, r._measures);
        }
        std::cout << "frame " << frame << "/" << frames << "  t " << t << "  substeps " << substeps
                  << "  max vol err " << r._maxVolumeError << std::endl;
    }
    if (pending._active) {
        PP::FlowMeasures last;
        if (integrator->Measure(state, last, err)) {
            log(pending._frame, pending._t, pending._substeps, pending._iterations, pending._res, last);
        } else {
            std::cerr << "measure of the last frame failed: " << err << std::endl;
        }
    }
    if (energy) {
        std::fclose(energy);
    }
    return 0;
}
