#include "fortran_io.h"

#include <cctype>
#include <map>
#include <sstream>

namespace {
std::map<int, FILE*> g_units;
int g_next_unit = 10;

const char* mode_from_status(const std::string& status) {
    if (status == "replace" || status == "new" || status == "w") return "w+";
    if (status == "old" || status == "r") return "r";
    if (status == "readwrite" || status == "rw") return "r+";
    return "a+";
}
}

int fortran_open_newunit(const std::string& filename,
                         const std::string& status,
                         const std::string&) {
    FILE* fp = std::fopen(filename.c_str(), mode_from_status(status));
    if (!fp && status == "old") {
        fp = std::fopen(filename.c_str(), "r");
    }
    if (!fp) return -1;
    int unit = g_next_unit++;
    g_units[unit] = fp;
    return unit;
}

int fortran_open_unit(int unit,
                      const std::string& filename,
                      const std::string& status,
                      const std::string&) {
    FILE* fp = std::fopen(filename.c_str(), mode_from_status(status));
    if (!fp && status == "old") {
        fp = std::fopen(filename.c_str(), "r");
    }
    if (!fp) return -1;
    auto existing = g_units.find(unit);
    if (existing != g_units.end()) {
        std::fclose(existing->second);
    }
    g_units[unit] = fp;
    return 0;  // Fortran IOSTAT=0 means success
}

int fortran_open_newunit(const std::string& filename,
                         const std::string& status,
                         const std::string& form,
                         int* iostat) {
    int unit = fortran_open_newunit(filename, status, form);
    if (iostat) *iostat = unit < 0 ? 1 : 0;
    return unit;
}

int fortran_open_unit(int unit,
                      const std::string& filename,
                      const std::string& status,
                      const std::string& form,
                      int* iostat) {
    int opened = fortran_open_unit(unit, filename, status, form);
    if (iostat) *iostat = opened < 0 ? 1 : 0;
    return opened;
}

FILE* fortran_get_unit(int unit) {
    auto it = g_units.find(unit);
    if (it != g_units.end()) return it->second;
    // Fortran default units: 5=stdin, 6=stdout
    if (unit == 5) return stdin;
    if (unit == 6) return stdout;
    return nullptr;
}

int fortran_close_unit(int unit) {
    auto it = g_units.find(unit);
    if (it == g_units.end()) return 0;
    int rc = std::fclose(it->second);
    g_units.erase(it);
    return rc;
}

int fortran_close_unit(int unit, int* iostat) {
    int rc = fortran_close_unit(unit);
    if (iostat) *iostat = rc;
    return rc;
}

int fortran_write_unit(int unit, const std::string& text) {
    FILE* fp = fortran_get_unit(unit);
    if (!fp) return -1;
    if (text.empty()) return std::fprintf(fp, "\n") < 0 ? -1 : 0;

    const char last = text.back();
    const bool raw_fragment =
        text.find('\n') != std::string::npos ||
        last == ' ' ||
        last == '\t' ||
        last == '\r';
    const int written = raw_fragment
        ? std::fprintf(fp, "%s", text.c_str())
        : std::fprintf(fp, "%s\n", text.c_str());
    return written < 0 ? -1 : 0;
}

int fortran_flush_unit(int unit) {
    FILE* fp = fortran_get_unit(unit);
    return fp ? std::fflush(fp) : 0;
}

void rewind(int unit) {
    FILE* fp = fortran_get_unit(unit);
    if (fp) std::rewind(fp);
}

bool fortran_file_exists(const std::string& filename) {
    FILE* fp = std::fopen(filename.c_str(), "r");
    if (fp) { std::fclose(fp); return true; }
    return false;
}

// =====================================================================
// Fortran list-directed READ implementation
// =====================================================================

std::vector<std::string> fortran_tokenize_line(FILE* fp) {
    std::vector<std::string> tokens;
    char buf[4096];
    if (!fgets(buf, sizeof(buf), fp)) return tokens;

    // Fortran list-directed READ can span multiple lines (continuation).
    // A line ending with a comma or that starts an implied-DO may continue.
    // We read at least one line; the caller can request more tokens later.
    std::string line(buf);
    // Remove trailing newline
    while (!line.empty() && (line.back() == '\n' || line.back() == '\r'))
        line.pop_back();

    size_t i = 0;
    while (i < line.size()) {
        // Skip whitespace and commas (Fortran list-directed separators)
        while (i < line.size() && (line[i] == ' ' || line[i] == '\t' || line[i] == ','))
            i++;
        if (i >= line.size()) break;

        // Quoted string
        if (line[i] == '\'' || line[i] == '"') {
            char q = line[i++];
            std::string tok;
            while (i < line.size() && line[i] != q) tok += line[i++];
            if (i < line.size()) i++; // skip closing quote
            tokens.push_back(tok);
            continue;
        }

        // Regular token
        std::string tok;
        while (i < line.size() && line[i] != ' ' && line[i] != '\t' &&
               line[i] != ',' && line[i] != '\'' && line[i] != '"')
            tok += line[i++];

        // Check for Fortran repeat notation: "3*0.02"
        auto star = tok.find('*');
        if (star != std::string::npos && star > 0 && star < tok.size() - 1) {
            // Verify left side is a positive integer
            std::string left = tok.substr(0, star);
            bool is_repeat = true;
            for (char c : left) { if (!std::isdigit(c)) { is_repeat = false; break; } }
            if (is_repeat) {
                int rep = std::stoi(left);
                std::string val = tok.substr(star + 1);
                for (int r = 0; r < rep; r++) tokens.push_back(val);
                continue;
            }
        }

        tokens.push_back(tok);
    }
    return tokens;
}

void fortran_assign_token(const std::string& tok, int* out) {
    if (tok.empty()) { *out = 0; return; }
    // FortranArray1D<bool> stores as int; handle T/F/.TRUE./.FALSE.
    char c = std::toupper(tok[0]);
    if (c == 'T' || c == 'F' || c == '.') {
        if (c == '.') c = (tok.size() > 1) ? std::toupper(tok[1]) : 'F';
        *out = (c == 'T' || c == '1') ? 1 : 0;
        return;
    }
    try { *out = std::stoi(tok); } catch (...) { *out = 0; }
}

void fortran_assign_token(const std::string& tok, float* out) {
    // Handle Fortran D exponent: 1.23D+04 → 1.23E+04
    std::string s = tok;
    for (auto& c : s) { if (c == 'd' || c == 'D') c = 'E'; }
    try { *out = std::stof(s); } catch (...) { *out = 0.0f; }
}

void fortran_assign_token(const std::string& tok, double* out) {
    std::string s = tok;
    for (auto& c : s) { if (c == 'd' || c == 'D') c = 'E'; }
    try { *out = std::stod(s); } catch (...) { *out = 0.0; }
}

void fortran_assign_token(const std::string& tok, bool* out) {
    if (tok.empty()) { *out = false; return; }
    char c = std::toupper(tok[0]);
    // Handle: T, .TRUE., 1, F, .FALSE., 0
    if (c == '.') {
        // .TRUE. or .FALSE.
        c = (tok.size() > 1) ? std::toupper(tok[1]) : 'F';
    }
    *out = (c == 'T' || c == '1');
}

void fortran_assign_token(const std::string& tok, std::string* out) {
    *out = tok;
}
