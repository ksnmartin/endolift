// endog: branch-and-cut for two-stage facility location with endogenous probabilities (decomposed linked lifted-cut bound).
// usage: endog [--master dec|mcb|mcb-dyn|mcb-mip] [--unlinked] [--mip-time sec] [--seed-c2l] [--bpos0] [--gen-only] [--inst file | --n N --m M --S S --seed K [--endo strong|weak --strength x --density x --fmult x --rng pcg64|mt]]
//              [--workers T] [--wall sec] [--tol x] [--node-tol x] [--root-loose x] [--max-age k] [--lam x]
//              [--no-heur] [--root-only] [--out file.csv] [--save-inst file] [-v]
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include "bnc.h"
#include "instance.h"
#include "mcb.h"

int main(int argc, char** argv) {
  setvbuf(stdout, nullptr, _IOLBF, 0);   // progress lines must reach a redirected log as they happen
  std::string inst_path, out, save_inst, endo = "strong", master = "dec";   // master: dec | mcb | mcb-dyn | mcb-mip
  std::string rng = "pcg64";                                                // generator stream: pcg64 (table instances) | mt (cgen_*)
  bool bpos0 = false, seed_c2l = false, gen_only = false;
  int n = 15, m = 15, S = 20, seed = 1; double strength = 0.9, density = 1.0, fmult = 1.0;
  BncOptions opt;
  auto need = [&](int& i) { if (++i >= argc) { std::fprintf(stderr, "missing value for %s\n", argv[i - 1]); std::exit(2); } return argv[i]; };
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--inst") inst_path = need(i);
    else if (a == "--master") master = need(i);
    else if (a == "--n") n = std::atoi(need(i));
    else if (a == "--m") m = std::atoi(need(i));
    else if (a == "--S") S = std::atoi(need(i));
    else if (a == "--seed") seed = std::atoi(need(i));
    else if (a == "--endo") endo = need(i);
    else if (a == "--strength") strength = std::atof(need(i));
    else if (a == "--density") density = std::atof(need(i));
    else if (a == "--fmult") fmult = std::atof(need(i));
    else if (a == "--rng") rng = need(i);
    else if (a == "--workers") opt.workers = std::atoi(need(i));
    else if (a == "--wall") opt.wall = std::atof(need(i));
    else if (a == "--tol") opt.tol = std::atof(need(i));
    else if (a == "--node-tol") opt.node_tol = std::atof(need(i));
    else if (a == "--root-loose") opt.root_loose = std::atof(need(i));
    else if (a == "--max-age") opt.max_age = std::atoi(need(i));
    else if (a == "--lam") opt.lam = std::atof(need(i));
    else if (a == "--no-heur") opt.heur = false;
    else if (a == "--sparse") opt.sparse = true;
    else if (a == "--lazy-mc") opt.lazy_mc = true;
    else if (a == "--root-only") opt.root_only = true;
    else if (a == "--unlinked") opt.unlinked = true;
    else if (a == "--bpos0") bpos0 = true;
    else if (a == "--seed-c2l") seed_c2l = true;
    else if (a == "--mip-time") opt.mip_time = std::atof(need(i));
    else if (a == "--out") out = need(i);
    else if (a == "--save-inst") save_inst = need(i);
    else if (a == "--gen-only") gen_only = true;                       // generate + --save-inst, no solve
    else if (a == "-v") opt.verbose = true;
    else { std::fprintf(stderr, "unknown argument %s\n", argv[i]); return 2; }
  }
  Instance inst = inst_path.empty() ? make_instance(n, m, S, seed, strength, density, fmult, endo, bpos0, rng) : load_instance(inst_path);
  if (!save_inst.empty()) save_instance(inst, save_inst);
  if (gen_only) { std::printf("%s saved to %s\n", inst.name.c_str(), save_inst.c_str()); return 0; }
  std::printf("%s\n", inst.describe().c_str()); std::fflush(stdout);
  const char* tag = master == "dec" ? (opt.unlinked ? "c2r-cpp" : "dec-cpp") : master == "mcb" ? "mcb-cpp" : master == "mcb-dyn" ? "mcbdyn-cpp" : master == "mcb-mip" ? "mcbmip-cpp" : nullptr;
  if (!tag) { std::fprintf(stderr, "unknown master %s (dec | mcb | mcb-dyn | mcb-mip)\n", master.c_str()); return 2; }
  BncStats st;
  if (master == "dec") st = branch_and_cut_dec(inst, opt);
  else if (master == "mcb-mip") st = mcb_mip_loop(inst, opt);
  else if (seed_c2l) {   // McCormick master started from the converged root pool of the linked master
    BncOptions ro = opt; ro.root_only = true; BncStats rs; rs.name = inst.name;
    NodeWorker w(inst, ro); w.run_root(rs, 0.0);
    std::printf("  [seed-c2l] linked root pool: %zu vertices\n", w.pool.total());
    st = branch_and_cut_mcb(inst, opt, master == "mcb-dyn", &w.pool);
    tag = master == "mcb" ? "mcbseed-cpp" : "mcbdynseed-cpp";
  } else st = branch_and_cut_mcb(inst, opt, master == "mcb-dyn");
  std::printf("%s\n", st.summary(tag).c_str());
  if (!out.empty()) {
    std::filesystem::path op(out);
    if (op.has_parent_path()) std::filesystem::create_directories(op.parent_path());
    bool fresh = !std::ifstream(out).good();
    std::ofstream f(out, std::ios::app);
    if (fresh) f << st.csv_header() << '\n';
    f << st.csv_row(tag) << '\n';
  }
  return 0;
}
