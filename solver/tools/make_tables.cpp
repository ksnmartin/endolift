// Writes Tables 1, 2, 4 and 5 of the paper into paper/main.tex (between the % TABk-BEGIN / % TABk-END markers) from
// the result CSVs in experiments/solver/.   usage: solver/build/make_tables [repo_root]   (default: the current directory)
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <regex>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using Row = std::map<std::string, std::string>;
using Key = std::tuple<int, int, std::string>;
static std::string E, ROOT;

static std::vector<Row> read_csv(const std::string& path) {
  std::vector<Row> rows; std::ifstream in(path); std::string line;
  if (!in || !std::getline(in, line)) return rows;
  std::vector<std::string> head; { std::stringstream ss(line); std::string c; while (std::getline(ss, c, ',')) head.push_back(c); }
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::vector<std::string> f; { std::stringstream ss(line); std::string c; while (std::getline(ss, c, ',')) f.push_back(c); }
    Row r; for (size_t i = 0; i < head.size() && i < f.size(); ++i) r[head[i]] = f[i];
    rows.push_back(r);
  }
  return rows;
}
static double num(const Row& r, const char* k) { auto it = r.find(k); return it == r.end() ? NAN : std::atof(it->second.c_str()); }
static std::string commas(long v) {        // 12345 -> "12,345"
  std::string s = std::to_string(v), o; int c = 0;
  for (int i = (int)s.size() - 1; i >= 0; --i) { o.insert(o.begin(), s[i]); if (++c % 3 == 0 && i > 0 && s[i - 1] != '-') o.insert(o.begin(), ','); }
  return o;
}
static std::string f0(double v) { return commas((long)std::llround(std::nearbyint(v))); }   // "{v:,.0f}"
static std::string pct(double v) { char b[32]; std::snprintf(b, sizeof b, "%.2f\\%%", v); return b; }
static std::string gap_of(const Row& r) { return pct((num(r, "OPT") - num(r, "lb")) / std::fabs(num(r, "OPT")) * 100.0); }
static bool proven(const Row& r) { auto it = r.find("proven"); return it != r.end() && it->second == "1"; }

// (n|m, seed, master) -> row; later files / later rows win. pattern must have two capture groups.
static void load(std::map<Key, Row>& rows, const std::string& file, const std::regex& pat, const std::vector<std::string>& masters, bool split_root = false) {
  for (const Row& r : read_csv(E + "/" + file)) {
    std::smatch m; const std::string& name = r.at("name");
    if (!std::regex_search(name, m, pat)) continue;
    std::string master = r.at("master");
    if (!masters.empty()) { bool ok = false; for (auto& x : masters) ok |= x == master; if (!ok) continue; }
    if (split_root && master == "mcbdyn-cpp" && r.at("nodes") == "0") master = "mcbdyn-root";   // --root-only run: root value only
    rows[Key(std::stoi(m[1]), std::stoi(m[2]), master)] = r;
  }
}
static const Row* get(const std::map<Key, Row>& rows, int a, int b, const std::string& master) {
  auto it = rows.find(Key(a, b, master)); return it == rows.end() ? nullptr : &it->second;
}
static std::pair<std::string, std::string> tree(const Row* r) {   // nodes, time | unproven (gap)
  if (!r) return {"---", "---"};
  std::string nodes = commas((long)num(*r, "nodes"));
  if (proven(*r)) return {nodes, f0(num(*r, "total_time"))};
  return {nodes, "unproven (" + gap_of(*r) + ")"};
}
static std::string splice(const std::string& tex, const std::string& tag, const std::vector<std::string>& lines) {
  size_t i = tex.find("% " + tag + "-BEGIN"); if (i == std::string::npos) return tex;
  i = tex.find('\n', i) + 1; size_t j = tex.find("% " + tag + "-END");
  std::string body; for (auto& l : lines) body += l + "\n";
  return tex.substr(0, i) + body + tex.substr(j);
}

int main(int argc, char** argv) {
  ROOT = argc > 1 ? argv[1] : "."; E = ROOT + "/experiments/solver";
  const std::regex main_pat("gen_n(\\d+)_m15_S20_seed(\\d+)_"), bp_pat("cgen_n(\\d+)_m15_S20_seed(\\d+)_"), ms_pat("cgen_n30_m(\\d+)_S20_seed(\\d+)_");
  std::map<Key, Row> rows;
  load(rows, "validate.csv", main_pat, {"dec-cpp"});                       // loaded machine, overwritten below where idle runs exist
  load(rows, "validate50.csv", main_pat, {"dec-cpp"});
  load(rows, "mpi_4050b_single.csv", main_pat, {"dec-cpp"});               // idle, 40/50 seeds 1-3
  load(rows, "extra_seeds.csv", main_pat, {"dec-cpp"});                    // idle, seeds 4-5
  load(rows, "mcb_cpp.csv", main_pat, {"dec-cpp", "mcb-cpp", "mcbdyn-cpp"});   // idle, 25-35 dec seeds 1-3 + all McCormick
  load(rows, "unlinked_cpp.csv", main_pat, {"c2r-cpp"});                   // Table 1 unlinked column: dec --unlinked --root-only
  // Table 2: nodes and time of the three masters, root time of the decomposed solver
  std::vector<std::string> out2;
  for (int n : {25, 30, 35, 40, 50})
    for (int s : {1, 2, 3, 4, 5}) {
      auto a = tree(get(rows, n, s, "mcb-cpp")), b = tree(get(rows, n, s, "mcbdyn-cpp")), c = tree(get(rows, n, s, "dec-cpp"));
      if (s >= 4 && a.first == "---") a = {"", ""};                        // static not run on seeds 4-5
      const Row* d = get(rows, n, s, "dec-cpp");
      out2.push_back(std::to_string(n) + " & " + std::to_string(s) + " & " + a.first + " & " + a.second + " & " + b.first + " & " + b.second +
                     " & " + c.first + " & " + (d ? f0(num(*d, "root_time")) : "---") + " & " + c.second + "\\\\");
    }
  // Table 1: root gaps (F* - v)/F* of MC-B, unlinked and linked at the converged root LP, seeds 1-3, two column blocks
  auto gap = [&](int n, int s, const std::string& master) -> std::string {
    const Row* r = get(rows, n, s, master); const Row* o = get(rows, n, s, "dec-cpp");
    if (!r || !o || !proven(*o)) return "---";
    double F = num(*o, "OPT"); return pct((F - num(*r, "root_value")) / std::fabs(F) * 100.0);
  };
  std::vector<std::pair<int, int>> cells; for (int n : {25, 30, 35, 40, 50}) for (int s : {1, 2, 3}) cells.push_back({n, s});
  size_t half = (cells.size() + 1) / 2;
  std::vector<std::string> out1;
  for (size_t k = 0; k < half; ++k) {
    auto L = cells[k];
    std::string row = std::to_string(L.first) + " & " + std::to_string(L.second) + " & " + gap(L.first, L.second, "mcb-cpp") + " & " +
                      gap(L.first, L.second, "c2r-cpp") + " & " + gap(L.first, L.second, "dec-cpp") + " & ";
    if (k + half < cells.size()) {
      auto R = cells[k + half];
      row += std::to_string(R.first) + " & " + std::to_string(R.second) + " & " + gap(R.first, R.second, "mcb-cpp") + " & " +
             gap(R.first, R.second, "c2r-cpp") + " & " + gap(R.first, R.second, "dec-cpp");
    } else row += "& & & & ";
    out1.push_back(row + "\\\\");
  }
  // Table 4: b_s = 0 for the high-demand scenarios; roots of the three masters, dec trees, mcb-dyn trees at 25-30
  std::map<Key, Row> bp; load(bp, "bpos0.csv", bp_pat, {}, true);
  auto bgap = [&](int n, int s, const std::string& master) -> std::string {
    const Row* r = master == "mcbdyn-cpp" ? get(bp, n, s, "mcbdyn-root") : nullptr;
    if (!r) r = get(bp, n, s, master);
    const Row* o = get(bp, n, s, "dec-cpp");
    if (!r || !o) return "---";
    std::string dag = proven(*o) ? "" : "$^\\dagger$";        // unproven: gap against the incumbent (an upper bound on the true gap)
    double F = num(*o, "OPT"); return pct((F - num(*r, "root_value")) / std::fabs(F) * 100.0) + dag;
  };
  std::vector<std::string> out4;
  for (int n : {25, 30, 35, 40})
    for (int s : {1, 2, 3}) {
      auto d = tree(get(bp, n, s, "dec-cpp"));
      const Row* mr = get(bp, n, s, "mcbdyn-cpp");
      auto m = mr ? tree(mr) : (n > 30 ? std::make_pair(std::string(""), std::string("")) : std::make_pair(std::string("---"), std::string("---")));
      out4.push_back(std::to_string(n) + " & " + std::to_string(s) + " & " + bgap(n, s, "mcbdyn-cpp") + " & " + bgap(n, s, "c2r-cpp") + " & " +
                     bgap(n, s, "dec-cpp") + " & " + d.first + " & " + d.second + " & " + m.first + " & " + m.second + "\\\\");
    }
  // Table 5: n = 30, m in {15, 50, 100, 200}, decomposed solver
  std::map<Key, Row> ms; load(ms, "msens.csv", ms_pat, {}, true);
  std::vector<std::string> out5;
  for (int m : {15, 50, 100, 200})
    for (int s : {1, 2, 3}) {
      const Row* r = get(ms, m, s, "dec-cpp");
      if (!r) { out5.push_back(std::to_string(m) + " & " + std::to_string(s) + " & --- & --- & --- & --- & --- & ---\\\\"); continue; }
      std::string g = pct(num(*r, "root_gap")) + (proven(*r) ? "" : "$^\\dagger$");
      std::string tm = proven(*r) ? f0(num(*r, "total_time")) : "unproven (" + gap_of(*r) + ")";
      out5.push_back(std::to_string(m) + " & " + std::to_string(s) + " & " + g + " & " + f0(num(*r, "root_time")) + " & " +
                     commas((long)num(*r, "root_cuts")) + " & " + commas((long)(num(*r, "cuts") - num(*r, "root_cuts"))) + " & " +
                     commas((long)num(*r, "nodes")) + " & " + tm + "\\\\");
    }
  std::string mp = ROOT + "/paper/main.tex";
  std::ifstream in(mp); if (!in) { std::fprintf(stderr, "cannot read %s\n", mp.c_str()); return 1; }
  std::stringstream buf; buf << in.rdbuf(); std::string tex = buf.str(); in.close();
  tex = splice(tex, "TAB1", out1); tex = splice(tex, "TAB2", out2); tex = splice(tex, "TAB4", out4); tex = splice(tex, "TAB5", out5);
  std::ofstream o(mp); o << tex; o.close();
  for (auto* v : {&out1, &out2, &out4, &out5}) { for (auto& l : *v) std::puts(l.c_str()); std::puts(""); }
  return 0;
}
