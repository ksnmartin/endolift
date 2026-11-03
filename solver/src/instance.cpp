#include "instance.h"
#include "nprng.h"
#include <cstdio>
#include "recourse.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <random>
#include <sstream>
#include <stdexcept>

void Instance::finalize() {
  bpos.assign(S, 0.0);
  double bsum = 0.0;
  for (int s = 0; s < S; ++s) {
    bpos[s] = b[s];
    for (int j = 0; j < n; ++j) bpos[s] += std::min(a(s, j), 0.0);
    if (bpos[s] < -1e-12) throw std::runtime_error("p_s(y) negative at some binary y");
    bsum += b[s];
  }
  if (std::fabs(bsum - 1.0) > 1e-9) throw std::runtime_error("sum_s b_s != 1");
  for (int j = 0; j < n; ++j) {
    double c = 0.0;
    for (int s = 0; s < S; ++s) c += a(s, j);
    if (std::fabs(c) > 1e-9) throw std::runtime_error("columns of A must sum to 0");
  }
}

double Instance::p(int s, const double* y) const {
  double v = b[s];
  for (int j = 0; j < n; ++j) v += a(s, j) * y[j];
  return v;
}

std::string Instance::describe() const {
  int nz = 0;
  for (double v : A) nz += v != 0.0;
  double mb = bpos.empty() ? 0.0 : *std::min_element(bpos.begin(), bpos.end());
  std::ostringstream os;
  os << (name.empty() ? "instance" : name) << ": n=" << n << " m=" << m << " S=" << S << " nonzero a_sj=" << nz
     << " min bpos=" << std::fixed << std::setprecision(4) << mb;
  return os.str();
}

// ---------------------------------------------------------------------------------------------------- file format
// Plain text: "endog-instance 1", name, endo, "n m S pen", then cap, f, C (n rows of m), d (S rows of m), b, A (S rows of n).
static void read_vec(std::istream& in, std::vector<double>& v, size_t count) {
  v.resize(count);
  for (size_t i = 0; i < count; ++i)
    if (!(in >> v[i])) throw std::runtime_error("instance file truncated");
}

Instance load_instance(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open instance file " + path);
  std::string magic; int ver;
  in >> magic >> ver;
  if (magic != "endog-instance" || ver != 1) throw std::runtime_error("bad instance header in " + path);
  Instance I;
  in >> I.name >> I.endo >> I.n >> I.m >> I.S >> I.pen;
  read_vec(in, I.cap, I.n); read_vec(in, I.f, I.n);
  read_vec(in, I.C, (size_t)I.n * I.m); read_vec(in, I.d, (size_t)I.S * I.m);
  read_vec(in, I.b, I.S); read_vec(in, I.A, (size_t)I.S * I.n);
  I.finalize();
  return I;
}

static void write_vec(std::ostream& out, const std::vector<double>& v, size_t per_line) {
  for (size_t i = 0; i < v.size(); ++i) out << v[i] << ((i + 1) % per_line == 0 || i + 1 == v.size() ? '\n' : ' ');
}

void save_instance(const Instance& I, const std::string& path) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write instance file " + path);
  out << std::setprecision(17);
  out << "endog-instance 1\n" << I.name << '\n' << I.endo << '\n' << I.n << ' ' << I.m << ' ' << I.S << ' ' << I.pen << '\n';
  write_vec(out, I.cap, I.n); write_vec(out, I.f, I.n); write_vec(out, I.C, I.m); write_vec(out, I.d, I.m);
  write_vec(out, I.b, I.S); write_vec(out, I.A, I.n);
}

// ---------------------------------------------------------------------------------------------------- generator
static double round2(double v) { return std::nearbyint(v * 100.0) / 100.0; }   // numpy.round: half to even
static std::string pyfloat(double v) {                                               // 1 -> "1.0": the form used in the instance names
  char b[32]; std::snprintf(b, sizeof b, "%.15g", v); std::string s(b);
  if (s.find_first_of(".e") == std::string::npos) s += ".0";
  return s;
}

Instance make_instance(int n, int m, int S, int seed, double strength, double density, double fmult, const std::string& endo, bool bpos0,
                       const std::string& rng_kind) {
  if (rng_kind != "pcg64" && rng_kind != "mt") throw std::runtime_error("unknown rng kind " + rng_kind + " (pcg64|mt)");
  const bool pcg = rng_kind == "pcg64";
  nprng::PCG64 np(seed);                 // numpy.random.default_rng(seed): the stream of the paper's table instances
  std::mt19937_64 mt(seed);              // the stream of the cgen_* instances (sensitivity, sparse, bpos0, customers)
  auto U = [&](double lo, double hi) { return pcg ? np.uniform(lo, hi) : std::uniform_real_distribution<double>(lo, hi)(mt); };
  Instance I; I.n = n; I.m = m; I.S = S; I.pen = 10.0; I.endo = endo;
  I.C.resize((size_t)n * m);
  for (double& c : I.C) c = round2(U(1.0, 5.0));
  std::vector<double> scale(S);
  for (int s = 0; s < S; ++s) scale[s] = S == 1 ? 0.5 : 0.5 + 1.1 * s / (S - 1);
  I.d.resize((size_t)S * m);
  double mean_total = 0.0;
  for (int s = 0; s < S; ++s)
    for (int k = 0; k < m; ++k) { I.d[(size_t)s * m + k] = round2(U(1.5, 4.5) * scale[s]); mean_total += I.d[(size_t)s * m + k]; }
  mean_total /= S;
  I.cap.resize(n);
  for (double& c : I.cap) c = round2(U(0.7, 1.3) * 1.8 * mean_total / n);
  // Calibrate fixed costs: f_i = (recourse saved by site i when all others are open, at the mean demand) x U(0.6, 1.6).
  Instance tmp = I; tmp.f.assign(n, 0.0); tmp.b.assign(S, 1.0 / S); tmp.A.assign((size_t)S * n, 0.0); tmp.finalize();
  RecourseDual dual(tmp);
  std::vector<double> dbar(m, 0.0);
  for (int s = 0; s < S; ++s) for (int k = 0; k < m; ++k) dbar[k] += I.d[(size_t)s * m + k] / S;
  std::vector<double> al, be, c(I.cap);
  double q_all = dual.solve(dbar.data(), c.data(), al, be);
  std::vector<double> saving(n);
  for (int i = 0; i < n; ++i) {
    c = I.cap; c[i] = 0.0;
    saving[i] = std::max(dual.solve(dbar.data(), c.data(), al, be) - q_all, 0.05 * q_all / n);
  }
  I.f.resize(n);
  for (int i = 0; i < n; ++i) I.f[i] = round2(saving[i] * U(0.6, 1.6) * fmult);
  std::vector<char> active(n);
  bool any = false;
  for (int j = 0; j < n; ++j) { active[j] = U(0.0, 1.0) < density; any |= active[j]; }
  if (!any) active[pcg ? (int)np.integers(n) : std::uniform_int_distribution<int>(0, n - 1)(mt)] = 1;
  std::vector<double> u(n);
  for (int j = 0; j < n; ++j) u[j] = U(0.3, 1.0) * active[j];
  I.b.assign(S, 0.0); I.A.assign((size_t)S * n, 0.0);
  if (endo == "weak") {
    for (int s = 0; s < S; ++s) I.b[s] = 1.0 / S;
    std::vector<double> shift(S);
    double mean = 0.0, mx = 0.0;
    for (int s = 0; s < S; ++s) { shift[s] = S == 1 ? -1.0 : -1.0 + 2.0 * s / (S - 1); mean += shift[s] / S; }
    for (int s = 0; s < S; ++s) { shift[s] -= mean; mx = std::max(mx, std::fabs(shift[s])); }
    for (int s = 0; s < S; ++s) for (int j = 0; j < n; ++j) I.A[(size_t)s * n + j] = (shift[s] / mx) * u[j] * strength / (n * S);
  } else if (endo == "strong") {
    std::vector<double> gain(S), lost(S);
    double bsum = 0.0, gsum = 0.0, lsum = 0.0, usum = 0.0;
    for (int s = 0; s < S; ++s) { bool low = s < S / 2; I.b[s] = low ? 1.0 : (bpos0 ? 0.0 : 0.1); bsum += I.b[s]; gain[s] = low ? 0.0 : scale[s] * scale[s]; gsum += gain[s]; }
    for (int s = 0; s < S; ++s) { I.b[s] /= bsum; gain[s] /= gsum; lost[s] = s < S / 2 ? I.b[s] : 0.0; lsum += lost[s]; }
    for (int j = 0; j < n; ++j) usum += u[j];
    for (int s = 0; s < S; ++s)
      for (int j = 0; j < n; ++j) I.A[(size_t)s * n + j] = strength * (-lost[s] + gain[s] * lsum) * (u[j] / usum);
  } else {
    throw std::runtime_error("unknown endo mode " + endo);
  }
  std::ostringstream tag;
  if (pcg) tag << "gen_n" << n << "_m" << m << "_S" << S << "_seed" << seed << "_str" << pyfloat(strength) << "_den" << pyfloat(density);
  else tag << "cgen_n" << n << "_m" << m << "_S" << S << "_seed" << seed << "_str" << strength << "_den" << density;
  if (fmult != 1.0) tag << "_fm" << fmult;
  if (bpos0) tag << "_bpos0";
  if (endo != "weak") tag << "_" << endo;
  I.name = tag.str();
  I.finalize();
  return I;
}
