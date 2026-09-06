#include <algorithm>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr const char* kSchemaMarker = "Benchmark schema: gpu-benchmark/v2";
constexpr const char* kTimingRelation = "wall-outer;report-timers-nested-not-additive";

std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::vector<std::string> split_markdown_row(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, '|')) fields.push_back(trim(field));
    if (!fields.empty() && fields.front().empty()) fields.erase(fields.begin());
    if (!fields.empty() && fields.back().empty()) fields.pop_back();
    return fields;
}

bool is_separator(const std::vector<std::string>& fields) {
    return !fields.empty() && std::all_of(fields.begin(), fields.end(), [](const auto& field) {
        return field.find_first_not_of("-:") == std::string::npos;
    });
}

bool parse_finite_number(const std::string& value, double& parsed) {
    if (value.empty()) return false;
    char* end = nullptr;
    errno = 0;
    parsed = std::strtod(value.c_str(), &end);
    return errno == 0 && end == value.c_str() + value.size() && std::isfinite(parsed);
}

bool is_finite_number(const std::string& value) {
    double parsed = 0.0;
    return parse_finite_number(value, parsed);
}

bool parse_integer(const std::string& value, int& parsed) {
    if (value.empty()) return false;
    const auto result = std::from_chars(value.data(), value.data() + value.size(), parsed);
    return result.ec == std::errc{} && result.ptr == value.data() + value.size();
}

bool is_one_of(const std::string& value, std::initializer_list<const char*> allowed) {
    return std::any_of(allowed.begin(), allowed.end(), [&](const char* candidate) {
        return value == candidate;
    });
}

std::string join_key(std::initializer_list<std::string> fields) {
    std::string key;
    for (const auto& field : fields) {
        if (!key.empty()) key += '\x1f';
        key += field;
    }
    return key;
}

struct ExpectedRequest {
    const char* workload;
    const char* backend;
    const char* solver;
    const char* precision;
    const char* thermal;
    int batch;
    const char* active_ratio;
    const char* runtime_mask_change;
};

const std::vector<ExpectedRequest>& expected_requests() {
    static const std::vector<ExpectedRequest> requests = {
        {"small-single", "direct", "eigen", "full", "disabled", 1, "1.000", "no"},
        {"small-single", "graph", "eigen", "full", "disabled", 1, "1.000", "no"},
        {"medium-multi", "direct", "eigen", "full", "disabled", 1, "1.000", "no"},
        {"medium-multi", "graph", "cusolver", "full", "disabled", 1, "1.000", "no"},
        {"medium-multi", "persistent", "eigen", "full", "disabled", 1, "1.000", "no"},
        {"medium-multi", "fallback", "eigen", "full", "disabled", 1, "1.000", "no"},
        {"medium-multi", "graph", "cusolver", "full", "disabled", 1, "1.000", "yes"},
        {"large-single", "direct", "cusolver", "full", "disabled", 1, "1.000", "no"},
        {"large-single", "graph", "cusolver", "full", "disabled", 1, "1.000", "no"},
        {"medium-multi", "direct", "eigen", "full", "disabled", 128, "1.000", "no"},
        {"medium-multi-thermal", "graph", "cusolver", "full", "gpu", 1, "1.000", "no"},
    };
    return requests;
}

const std::vector<std::string>& expected_windows() {
    static const std::vector<std::string> windows = {
        "setup", "cold-first-step", "cold-replay", "warm-up", "steady-state"
    };
    return windows;
}

std::string request_key(const ExpectedRequest& request) {
    return join_key({request.workload, request.backend, request.solver, request.precision,
                     request.thermal, std::to_string(request.batch), request.active_ratio,
                     request.runtime_mask_change});
}

std::string composite_key(const std::string& request, int repeat,
                          const std::string& window) {
    return join_key({request, std::to_string(repeat), window});
}

struct Check {
    int failures = 0;

    void require(bool condition, const std::string& message) {
        if (condition) return;
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
};

int validate(const std::string& path) {
    Check check;
    std::ifstream input(path);
    check.require(input.good(), "cannot open benchmark artifact: " + path);
    if (!input) return check.failures;

    std::vector<std::string> lines;
    std::string line;
    bool marker_found = false;
    while (std::getline(input, line)) {
        marker_found = marker_found || line == kSchemaMarker;
        lines.push_back(line);
    }
    check.require(marker_found, std::string("missing schema marker: ") + kSchemaMarker);
    const auto find_metadata = [&](const std::string& prefix) {
        return std::find_if(lines.begin(), lines.end(), [&](const auto& candidate) {
            return candidate.rfind(prefix, 0) == 0 && candidate.size() > prefix.size();
        });
    };
    const auto commit = find_metadata("Commit: ");
    const auto compiler = find_metadata("Compiler: ");
    check.require(commit != lines.end(), "missing commit metadata");
    check.require(commit == lines.end() || *commit != "Commit: unrecorded",
                  "commit metadata was not supplied by the reproduction command");
    check.require(compiler != lines.end(), "missing compiler metadata");

    const std::vector<std::string> expected_cpu_header = {
        "Workload", "Requested", "Thermal", "Phase", "Iterations",
        "Wall ms", "Per-step ms", "Finite"
    };
    std::vector<std::string> cpu_header;
    std::size_t cpu_header_line = lines.size();
    for (std::size_t index = 0; index < lines.size(); ++index) {
        if (lines[index].rfind("| Workload | Requested | Thermal |", 0) == 0) {
            cpu_header = split_markdown_row(lines[index]);
            cpu_header_line = index;
            break;
        }
    }
    check.require(!cpu_header.empty(), "missing CPU reference table header");
    check.require(cpu_header == expected_cpu_header, "CPU columns do not match the v2 order");
    std::set<std::string> cpu_phase_keys;
    if (cpu_header == expected_cpu_header) {
        int cpu_rows = 0;
        for (std::size_t index = cpu_header_line + 1; index < lines.size(); ++index) {
            if (lines[index] == "GPU execution phases:") break;
            if (lines[index].empty() || lines[index].front() != '|') continue;
            const auto fields = split_markdown_row(lines[index]);
            if (is_separator(fields)) continue;
            check.require(fields.size() == cpu_header.size(), "CPU row field count mismatch");
            if (fields.size() != cpu_header.size()) continue;
            ++cpu_rows;
            int iterations = -1;
            check.require(is_one_of(fields[0], {"small-single", "medium-multi", "large-single",
                                                "medium-multi-thermal"}),
                          "CPU row has invalid Workload");
            check.require(fields[1] == "cpu-reference", "CPU row has invalid Requested value");
            check.require(is_one_of(fields[2], {"cold", "thermal"}),
                          "CPU row has invalid Thermal value");
            check.require(is_one_of(fields[3], {"setup", "cold-first-step", "cold-replay",
                                                "steady-state"}),
                          "CPU row has invalid Phase value");
            check.require(parse_integer(fields[4], iterations) && iterations >= 0,
                          "CPU Iterations is not a non-negative integer");
            double wall_ms = -1.0;
            double per_step_ms = -1.0;
            check.require(parse_finite_number(fields[5], wall_ms) && wall_ms >= 0.0,
                          "CPU Wall ms is not a non-negative finite number");
            check.require(parse_finite_number(fields[6], per_step_ms) && per_step_ms >= 0.0,
                          "CPU Per-step ms is not a non-negative finite number");
            check.require(fields[7] == "yes", "CPU reference row is not finite");
            check.require(cpu_phase_keys.insert(join_key({fields[0], fields[3]})).second,
                          "duplicate CPU workload/phase row");
        }
        check.require(cpu_rows == 16, "expected exactly 16 CPU reference rows");
    }

    const std::vector<std::string> required_columns = {
        "Workload", "Requested backend", "Requested solver", "Requested precision",
        "Requested thermal", "Batch", "Active ratio", "Runtime mask change", "Mask updates",
        "Repeat", "Phase", "Measurement window", "Capture/replay", "Execution kind",
        "Iterations", "Resolved backend", "Resolved solver", "Precision", "Resolved thermal",
        "Timing relation", "Setup wall ms", "Cold first-step wall ms", "Warm-up wall ms",
        "Steady-state wall ms", "Capture-inclusive wall ms", "Replay wall ms",
        "Fallback wall ms", "Wall ms", "Per-step ms", "Steps/s", "Simulations/s",
        "GPU ms", "Transfer ms", "Mutual ms", "Assembly ms", "Solver ms",
        "State update ms", "Force ms", "Thermal ms", "Control/status ms", "Sync ms",
        "Solver status", "Residual", "Graph rebuild delta", "Graph rebuild total",
        "Fallback delta", "Fallback total", "CPU/GPU speedup", "GPU executed", "Finite",
        "Fallback reason"
    };

    std::vector<std::string> header;
    std::size_t header_line = lines.size();
    for (std::size_t index = 0; index < lines.size(); ++index) {
        if (lines[index].rfind("| Workload | Requested backend |", 0) == 0) {
            header = split_markdown_row(lines[index]);
            header_line = index;
            break;
        }
    }
    check.require(!header.empty(), "missing GPU result table header");
    if (header.empty()) return check.failures;

    std::map<std::string, std::size_t> column;
    for (std::size_t index = 0; index < header.size(); ++index) column[header[index]] = index;
    for (const auto& name : required_columns)
        check.require(column.count(name) == 1, "missing required column: " + name);
    check.require(header == required_columns, "GPU columns do not match the v2 order");
    if (check.failures != 0) return check.failures;

    const std::vector<std::string> optional_numeric_columns = {
        "Setup wall ms", "Cold first-step wall ms", "Warm-up wall ms", "Steady-state wall ms",
        "Capture-inclusive wall ms", "Replay wall ms", "Fallback wall ms", "Mutual ms",
        "Assembly ms", "State update ms", "Force ms", "Control/status ms", "Sync ms", "Residual"
    };
    const std::vector<std::string> required_numeric_columns = {
        "Wall ms", "Per-step ms", "Steps/s", "Simulations/s", "GPU ms",
        "Transfer ms", "Solver ms", "Thermal ms"
    };

    std::set<std::string> expected_composite_keys;
    for (const auto& request : expected_requests())
        for (int repeat = 0; repeat < 3; ++repeat)
            for (const auto& window : expected_windows())
                expected_composite_keys.insert(
                    composite_key(request_key(request), repeat, window));

    std::set<std::string> actual_composite_keys;
    bool graph_capture = false;
    bool graph_replay = false;
    bool fallback_row = false;
    int rows = 0;

    for (std::size_t index = header_line + 1; index < lines.size(); ++index) {
        if (lines[index].empty() || lines[index].front() != '|') continue;
        const auto fields = split_markdown_row(lines[index]);
        if (is_separator(fields)) continue;
        check.require(fields.size() == header.size(),
                      "row " + std::to_string(index + 1) + " has " +
                          std::to_string(fields.size()) + " fields; expected " +
                          std::to_string(header.size()));
        if (fields.size() != header.size()) continue;
        ++rows;
        const auto value = [&](const std::string& name) -> const std::string& {
            return fields[column.at(name)];
        };
        const auto row_label = "row " + std::to_string(index + 1);

        check.require(is_one_of(value("Workload"),
                                {"small-single", "medium-multi", "large-single",
                                 "medium-multi-thermal"}),
                      row_label + " has invalid Workload");
        check.require(is_one_of(value("Requested backend"),
                                {"direct", "graph", "persistent", "fallback"}),
                      row_label + " has invalid Requested backend");
        check.require(is_one_of(value("Requested solver"), {"eigen", "cusolver"}),
                      row_label + " has invalid Requested solver");
        check.require(is_one_of(value("Requested precision"),
                                {"full", "standard", "aggressive"}),
                      row_label + " has invalid Requested precision");
        check.require(is_one_of(value("Requested thermal"), {"disabled", "cpu", "gpu"}),
                      row_label + " has invalid Requested thermal");
        check.require(is_one_of(value("Runtime mask change"), {"yes", "no"}),
                      row_label + " has invalid Runtime mask change");
        check.require(is_one_of(value("Phase"),
                                {"setup", "first-step/capture-inclusive", "replay-only",
                                 "warm-up", "steady-state"}),
                      row_label + " has invalid Phase");
        check.require(is_one_of(value("Measurement window"),
                                {"setup", "cold-first-step", "cold-replay", "warm-up",
                                 "steady-state"}),
                      row_label + " has invalid Measurement window");
        check.require(is_one_of(value("Capture/replay"),
                                {"capture-inclusive", "replay-only", "not-applicable"}),
                      row_label + " has invalid Capture/replay");
        check.require(is_one_of(value("Execution kind"),
                                {"setup", "gpu", "fallback", "not-executed"}),
                      row_label + " has invalid Execution kind");
        check.require(is_one_of(value("Resolved backend"),
                                {"direct", "graph", "persistent", "fallback"}),
                      row_label + " has invalid Resolved backend");
        check.require(is_one_of(value("Resolved solver"), {"eigen", "cusolver"}),
                      row_label + " has invalid Resolved solver");
        check.require(is_one_of(value("Precision"), {"full", "standard", "aggressive"}),
                      row_label + " has invalid Precision");
        check.require(is_one_of(value("Resolved thermal"), {"disabled", "cpu", "gpu"}),
                      row_label + " has invalid Resolved thermal");
        check.require(is_one_of(value("Solver status"), {"not-run", "success"}),
                      row_label + " has invalid Solver status");
        check.require(is_one_of(value("GPU executed"), {"yes", "no"}),
                      row_label + " has invalid GPU executed");
        check.require(is_one_of(value("Finite"), {"yes", "no"}),
                      row_label + " has invalid Finite");
        check.require(value("Timing relation") == kTimingRelation,
                      row_label + " has invalid timing relation");
        check.require(value("Finite") == "yes",
                      row_label + " is not finite");

        int batch = -1;
        int mask_updates = -1;
        int repeat = -1;
        int iterations = -1;
        int graph_rebuild_delta = -1;
        int graph_rebuild_total = -1;
        int fallback_delta = -1;
        int fallback_total = -1;
        check.require(parse_integer(value("Batch"), batch) && batch > 0,
                      row_label + " Batch is not a positive integer");
        check.require(parse_integer(value("Mask updates"), mask_updates) && mask_updates >= 0,
                      row_label + " Mask updates is not a non-negative integer");
        check.require(parse_integer(value("Repeat"), repeat) && repeat >= 0 && repeat <= 2,
                      row_label + " Repeat is not an integer from 0 to 2");
        check.require(parse_integer(value("Iterations"), iterations) && iterations >= 0,
                      row_label + " Iterations is not a non-negative integer");
        check.require(parse_integer(value("Graph rebuild delta"), graph_rebuild_delta) &&
                          graph_rebuild_delta >= 0,
                      row_label + " Graph rebuild delta is not a non-negative integer");
        check.require(parse_integer(value("Graph rebuild total"), graph_rebuild_total) &&
                          graph_rebuild_total >= 0,
                      row_label + " Graph rebuild total is not a non-negative integer");
        check.require(parse_integer(value("Fallback delta"), fallback_delta) &&
                          fallback_delta >= 0,
                      row_label + " Fallback delta is not a non-negative integer");
        check.require(parse_integer(value("Fallback total"), fallback_total) &&
                          fallback_total >= 0,
                      row_label + " Fallback total is not a non-negative integer");
        check.require(graph_rebuild_delta <= graph_rebuild_total,
                      row_label + " Graph rebuild delta exceeds total");
        check.require(fallback_delta <= fallback_total,
                      row_label + " Fallback delta exceeds total");

        double active_ratio = -1.0;
        check.require(parse_finite_number(value("Active ratio"), active_ratio) &&
                          active_ratio >= 0.0 && active_ratio <= 1.0,
                      row_label + " Active ratio is not a number from 0 to 1");
        for (const auto& name : required_numeric_columns) {
            double parsed = -1.0;
            check.require(parse_finite_number(value(name), parsed) && parsed >= 0.0,
                          row_label + " has invalid " + name);
        }
        for (const auto& name : optional_numeric_columns) {
            double parsed = -1.0;
            check.require(value(name) == "n/a" ||
                              (parse_finite_number(value(name), parsed) && parsed >= 0.0),
                          row_label + " has invalid " + name);
        }

        const bool gpu_executed = value("GPU executed") == "yes";
        const bool fallback = value("Execution kind") == "fallback";
        const bool setup = value("Execution kind") == "setup";
        const bool gpu = value("Execution kind") == "gpu";
        check.require(gpu_executed == gpu,
                      row_label + " has inconsistent Execution kind and GPU executed");
        check.require(!setup || value("Measurement window") == "setup",
                      "setup execution kind appears outside setup window");
        if (fallback) {
            fallback_row = true;
            check.require(!gpu_executed, "fallback row claims GPU execution");
            check.require(value("CPU/GPU speedup") == "n/a", "fallback row has numeric speedup");
            check.require(value("Fallback reason") != "none", "fallback row lacks a reason");
            check.require(is_finite_number(value("Fallback wall ms")),
                          "fallback row lacks fallback latency");
        } else {
            check.require(value("Fallback wall ms") == "n/a",
                          "non-fallback row exposes fallback latency");
        }
        if (gpu) check.require(value("Fallback reason") == "none", "GPU row has fallback reason");

        const bool cpu_phase_exists =
            cpu_phase_keys.count(join_key({value("Workload"), value("Measurement window")})) == 1;
        const bool identical_cpu_work = value("Runtime mask change") == "no";
        const bool speedup_allowed = gpu && value("Finite") == "yes" && cpu_phase_exists &&
            identical_cpu_work && value("Measurement window") != "setup" &&
            value("Measurement window") != "warm-up";
        double speedup = -1.0;
        check.require(speedup_allowed
                          ? parse_finite_number(value("CPU/GPU speedup"), speedup) && speedup > 0.0
                          : value("CPU/GPU speedup") == "n/a",
                      row_label + " has invalid CPU/GPU speedup eligibility or value");

        const auto check_wall_alias = [&](const std::string& name, bool present) {
            check.require(present ? value(name) == value("Wall ms") : value(name) == "n/a",
                          row_label + " has invalid " + name + " alias");
        };
        check_wall_alias("Setup wall ms", value("Measurement window") == "setup");
        check_wall_alias("Cold first-step wall ms",
                         value("Measurement window") == "cold-first-step");
        check_wall_alias("Warm-up wall ms", value("Measurement window") == "warm-up");
        check_wall_alias("Steady-state wall ms", value("Measurement window") == "steady-state");
        check_wall_alias("Capture-inclusive wall ms",
                         value("Capture/replay") == "capture-inclusive");
        check_wall_alias("Replay wall ms", value("Capture/replay") == "replay-only");
        const bool resolved_graph_execution = gpu && value("Resolved backend") == "graph";
        const std::string expected_capture_replay = !resolved_graph_execution
            ? "not-applicable"
            : (graph_rebuild_delta > 0 ? "capture-inclusive" : "replay-only");
        check.require(value("Capture/replay") == expected_capture_replay,
                      row_label + " Capture/replay does not match resolved backend");
        graph_capture = graph_capture ||
            (resolved_graph_execution && value("Capture/replay") == "capture-inclusive");
        graph_replay = graph_replay ||
            (resolved_graph_execution && value("Capture/replay") == "replay-only");
        for (const auto* name : {"Mutual ms", "Assembly ms", "State update ms", "Force ms",
                                "Control/status ms", "Sync ms", "Residual"})
            check.require(value(name) == "n/a", std::string(name) + " must be n/a until instrumented");

        const int expected_iterations = value("Measurement window") == "setup" ? 0 :
            (value("Measurement window") == "warm-up" ? 5 :
             (value("Measurement window") == "steady-state" ? 10 : 1));
        const std::string expected_phase = value("Measurement window") == "cold-first-step"
            ? "first-step/capture-inclusive" :
            (value("Measurement window") == "cold-replay" ? "replay-only" :
             value("Measurement window"));
        const int expected_mask_updates = value("Runtime mask change") == "yes" &&
                value("Measurement window") != "setup" &&
                value("Measurement window") != "cold-first-step"
            ? expected_iterations : 0;
        check.require(iterations == expected_iterations,
                      row_label + " Iterations does not match Measurement window");
        check.require(value("Phase") == expected_phase,
                      row_label + " Phase does not match Measurement window");
        check.require(mask_updates == expected_mask_updates,
                      row_label + " Mask updates does not match request/window");

        const std::string actual_request_key = join_key({
            value("Workload"), value("Requested backend"), value("Requested solver"),
            value("Requested precision"), value("Requested thermal"), value("Batch"),
            value("Active ratio"), value("Runtime mask change")
        });
        const auto key = composite_key(actual_request_key, repeat, value("Measurement window"));
        check.require(actual_composite_keys.insert(key).second,
                      row_label + " duplicates a request/repeat/window key");
    }

    check.require(rows == static_cast<int>(expected_composite_keys.size()),
                  "expected exactly 165 GPU rows");
    check.require(actual_composite_keys.size() == expected_composite_keys.size(),
                  "expected exactly 165 unique request/repeat/window keys");
    for (const auto& key : expected_composite_keys)
        check.require(actual_composite_keys.count(key) == 1,
                      "missing expected request/repeat/window key");
    for (const auto& key : actual_composite_keys)
        check.require(expected_composite_keys.count(key) == 1,
                      "unexpected request/repeat/window key");
    check.require(graph_capture, "missing parseable Graph capture-inclusive row");
    check.require(graph_replay, "missing parseable Graph replay-only row");
    check.require(fallback_row, "missing fallback execution row");

    if (check.failures == 0)
        std::cout << "PASS: " << path << " (" << rows << " GPU rows)\n";
    return check.failures;
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: check_gpu_benchmark_schema <raw-output> [raw-output ...]\n";
        return 2;
    }
    int failures = 0;
    for (int index = 1; index < argc; ++index) failures += validate(argv[index]);
    return failures == 0 ? 0 : 1;
}
