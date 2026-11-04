#include "radar/output.hpp"
#include "test_helpers.hpp"
#include "test_framework.hpp"

#include <fstream>
#include <sstream>

using namespace radar;

namespace {

std::string slurp(const std::filesystem::path& p) {
    std::ifstream in(p, std::ios::binary);
    std::ostringstream s;
    s << in.rdbuf();
    return s.str();
}

}  // namespace

TEST(base64_reference_vectors) {
    auto enc = [](const std::string& s) { return base64_encode(std::vector<std::uint8_t>(s.begin(), s.end())); };
    CHECK(enc("").empty());
    CHECK(enc("f") == "Zg==");
    CHECK(enc("fo") == "Zm8=");
    CHECK(enc("foo") == "Zm9v");
    CHECK(enc("foobar") == "Zm9vYmFy");
}

TEST(writes_all_artifacts) {
    Scenario s = testing::small_scenario();
    Target t;
    t.range_m = 1200.0;
    t.velocity_mps = 12.0;
    t.rcs_m2 = 1.0;
    s.targets = {t};
    s.cpis = 2;
    StageTimer timer;
    const SimulationResult r = run_simulation(s, timer);
    const auto dir = std::filesystem::temp_directory_path() / "radar_sim_test_output";
    std::filesystem::create_directories(dir);

    write_ppm(dir / "map.ppm", r.cpis.back().rd_map, r.cpis.back().detections, {});
    const std::string ppm = slurp(dir / "map.ppm");
    CHECK(ppm.rfind("P6\n512 ", 0) == 0);

    write_rd_csv(dir / "map.csv", r.cpis.back().rd_map, r.axes);
    write_detections_csv(dir / "det.csv", r.detections);
    write_tracks_csv(dir / "trk.csv", r.tracks);
    write_timing_csv(dir / "time.csv", timer, s.cpis);
    write_html_viewer(dir / "viewer.html", s, r);
    CHECK(slurp(dir / "det.csv").rfind("cpi,range_m,velocity_mps", 0) == 0);
    CHECK(slurp(dir / "trk.csv").rfind("id,status", 0) == 0);
    CHECK(slurp(dir / "time.csv").find("pulse compression") != std::string::npos);
    const std::string csv = slurp(dir / "map.csv");
    CHECK(static_cast<std::size_t>(std::count(csv.begin(), csv.end(), '\n')) == r.cpis.back().rd_map.rows() + 1);
    const std::string html = slurp(dir / "viewer.html");
    CHECK(html.find("__RADAR_DATA__") == std::string::npos);
    CHECK(html.find("\"frames\":[\"") != std::string::npos);
    CHECK(html.find("<script src") == std::string::npos);  // self-contained

    std::ostringstream os;
    print_timing(os, timer, s.cpis, r.raw_samples);
    print_detections(os, r.cpis.back().detections, 1);
    print_tracks(os, r.tracks);
    CHECK(os.str().find("Msamples/s") != std::string::npos);
    std::filesystem::remove_all(dir);
    CHECK_THROWS(write_ppm(dir / "missing" / "x.ppm", r.cpis.back().rd_map, {}, {}));
}

TEST(monte_carlo_csv) {
    MonteCarloResult r;
    PdPoint p;
    p.snr_db = 10.0;
    p.trials = 10;
    p.detections = 5;
    p.pd = 0.5;
    p.pd_cfar = std::nan("");
    r.points.push_back(p);
    const auto path = std::filesystem::temp_directory_path() / "radar_pd_test.csv";
    write_pd_csv(path, r);
    CHECK(slurp(path).find("10,10,5,0.5,0,\n") != std::string::npos);
    std::ostringstream os;
    print_monte_carlo(os, r, Scenario{});
    CHECK(os.str().find("Pfa design") != std::string::npos);
    std::filesystem::remove(path);
}

TEST_MAIN()
