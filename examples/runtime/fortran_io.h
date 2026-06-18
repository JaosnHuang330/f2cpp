#ifndef FORTRAN_IO_H
#define FORTRAN_IO_H

#include "fortran_array.h"

#include <algorithm>
#include <cstdio>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#if defined(__has_include)
#  if __has_include(<mpi.h>)
#    include <mpi.h>
#    define F2CPP_HAS_MPI 1
#  endif
#endif

int fortran_open_newunit(const std::string& filename,
                         const std::string& status = "unknown",
                         const std::string& form = "formatted");
int fortran_open_newunit(const std::string& filename,
                         const std::string& status,
                         const std::string& form,
                         int* iostat);
int fortran_open_unit(int unit,
                      const std::string& filename,
                      const std::string& status = "unknown",
                      const std::string& form = "formatted");
int fortran_open_unit(int unit,
                      const std::string& filename,
                      const std::string& status,
                      const std::string& form,
                      int* iostat);
FILE* fortran_get_unit(int unit);
int fortran_close_unit(int unit);
int fortran_close_unit(int unit, int* iostat);
int fortran_write_unit(int unit, const std::string& text = "");
template<typename... Args>
int fortran_write_unit(int unit, const std::string& fmt, Args&&... args);
template<typename... Args>
int fprintf(int unit, const char* fmt, Args... args);
int fortran_flush_unit(int unit);
void rewind(int unit);
bool fortran_file_exists(const std::string& filename);

// =====================================================================
// Fortran internal WRITE:  WRITE(str, fmt) value
//
// Safe C++ equivalent of Fortran's "WRITE to character variable".
// These use snprintf internally — no buffer overflow risk.
// =====================================================================

// Format an integer with a Fortran-style field width (right-justified).
// Equivalent to: WRITE(str, '(I<width>)') ival
inline std::string fortran_format_int(int val, int width = 0) {
    if (width <= 0) return std::to_string(val);
    char buf[64];
    snprintf(buf, sizeof(buf), "%*d", width, val);
    return std::string(buf);
}

// Format a double with Fortran-style F/E format.
// fmt examples: "f8.3", "e12.4", "g15.6"
inline std::string fortran_format_real(double val, const std::string& fmt = "") {
    char buf[128];
    if (fmt.empty()) {
        snprintf(buf, sizeof(buf), "%g", val);
        return std::string(buf);
    }
    // Parse "fW.D" or "eW.D" or "gW.D"
    char type = 'f';
    int w = 0, d = 6;
    if (!fmt.empty()) {
        size_t pos = 0;
        // Skip optional '(' and find format letter
        while (pos < fmt.size() && (fmt[pos] == '(' || fmt[pos] == ' ')) pos++;
        if (pos < fmt.size() && (fmt[pos] == 'f' || fmt[pos] == 'F' ||
                                  fmt[pos] == 'e' || fmt[pos] == 'E' ||
                                  fmt[pos] == 'g' || fmt[pos] == 'G')) {
            type = std::tolower(fmt[pos]);
            pos++;
        }
        // Parse width
        while (pos < fmt.size() && fmt[pos] >= '0' && fmt[pos] <= '9') {
            w = w * 10 + (fmt[pos] - '0');
            pos++;
        }
        // Parse decimal places
        if (pos < fmt.size() && fmt[pos] == '.') {
            pos++;
            d = 0;
            while (pos < fmt.size() && fmt[pos] >= '0' && fmt[pos] <= '9') {
                d = d * 10 + (fmt[pos] - '0');
                pos++;
            }
        }
    }
    char spec[32];
    if (w > 0) {
        snprintf(spec, sizeof(spec), "%%%d.%d%c", w, d, type);
    } else {
        snprintf(spec, sizeof(spec), "%%.%d%c", d, type);
    }
    snprintf(buf, sizeof(buf), spec, val);
    return std::string(buf);
}

// ADJUSTL equivalent: remove leading whitespace
inline std::string fortran_adjustl(const std::string& s) {
    size_t pos = s.find_first_not_of(" \t");
    return (pos == std::string::npos) ? "" : s.substr(pos);
}

// TRIM equivalent: remove trailing whitespace
inline std::string fortran_trim(const std::string& s) {
    size_t pos = s.find_last_not_of(" \t");
    return (pos == std::string::npos) ? "" : s.substr(0, pos + 1);
}

// =====================================================================
// Fortran list-directed READ:  READ(unit,*) var1, var2, ...
//
// Semantics (Fortran standard):
//   - Record-based: each READ consumes at least one full line.
//   - Tokens separated by commas or whitespace.
//   - Repeat notation: "3*0.02" means 0.02 0.02 0.02.
//   - Logical values: T/.TRUE./1 → true,  F/.FALSE./0 → false.
//   - Strings: 'quoted text' or unquoted token.
//   - Unit 5 = stdin (default), Unit 6 = stdout (default).
//
// Usage:
//   fortran_read_list(5, &LPLOT, &KG1, &NGR);
//   fortran_read_list(5, &LREAD, &LWRITE);   // bool*
//   fortran_read_list(5, &X(1), &X(2), ...);  // or use array overload
// =====================================================================

// Tokenizer: reads one line, expands repeats, handles T/F
std::vector<std::string> fortran_tokenize_line(FILE* fp);

// Per-type assignment from a string token
void fortran_assign_token(const std::string& tok, int* out);
void fortran_assign_token(const std::string& tok, float* out);
void fortran_assign_token(const std::string& tok, double* out);
void fortran_assign_token(const std::string& tok, bool* out);
void fortran_assign_token(const std::string& tok, std::string* out);

// --- Template implementation ---
#include <sstream>

namespace detail {

// Base case: no more pointers to assign
inline void assign_tokens(const std::vector<std::string>&, size_t) {}

// Recursive variadic: assign token[idx] to *head, then recurse
template<typename T, typename... Rest>
void assign_tokens(const std::vector<std::string>& toks, size_t idx,
                   T* head, Rest*... rest) {
    if (idx < toks.size()) fortran_assign_token(toks[idx], head);
    assign_tokens(toks, idx + 1, rest...);
}

} // namespace detail

template<typename... Args>
void fortran_read_list(int unit, Args*... args) {
    FILE* fp = fortran_get_unit(unit);
    if (!fp) return;
    auto toks = fortran_tokenize_line(fp);
    // If we need more tokens than one line provides, read continuation lines
    constexpr size_t needed = sizeof...(args);
    while (toks.size() < needed && !feof(fp)) {
        auto more = fortran_tokenize_line(fp);
        if (more.empty()) break;
        toks.insert(toks.end(), more.begin(), more.end());
    }
    detail::assign_tokens(toks, 0, args...);
}

// Array overload: READ(unit,*) (arr(I), I=lb, ub)
// Reads values into a contiguous range of a FortranArray.
// One call reads as many lines as needed to fill lb..ub.
template<typename T>
void fortran_read_list_array(int unit, FortranArray1D<T>& arr, int lb, int ub) {
    FILE* fp = fortran_get_unit(unit);
    if (!fp) return;
    std::vector<std::string> toks;
    int count = ub - lb + 1;
    while (static_cast<int>(toks.size()) < count && !feof(fp)) {
        auto more = fortran_tokenize_line(fp);
        if (more.empty()) break;
        toks.insert(toks.end(), more.begin(), more.end());
    }
    for (int i = 0; i < count && i < static_cast<int>(toks.size()); i++) {
        fortran_assign_token(toks[i], &arr(lb + i));
    }
}

// Overload for raw pointer arrays
template<typename T>
void fortran_read_list_array(int unit, T* arr, int count) {
    FILE* fp = fortran_get_unit(unit);
    if (!fp) return;
    std::vector<std::string> toks;
    while (static_cast<int>(toks.size()) < count && !feof(fp)) {
        auto more = fortran_tokenize_line(fp);
        if (more.empty()) break;
        toks.insert(toks.end(), more.begin(), more.end());
    }
    for (int i = 0; i < count && i < static_cast<int>(toks.size()); i++) {
        fortran_assign_token(toks[i], &arr[i]);
    }
}

namespace detail {
inline void write_args(std::ostream&) {}
template<typename T, typename... Rest>
void write_args(std::ostream& os, const T& first, const Rest&... rest) {
    os << first;
    if constexpr (sizeof...(rest) > 0) { os << " "; write_args(os, rest...); }
}
} // namespace detail

template<typename... Args>
int fortran_write_unit(int unit, const std::string& fmt, Args&&... args) {
    if (fmt.find('%') != std::string::npos) {
        FILE* fp = fortran_get_unit(unit);
        return std::fprintf(fp ? fp : stdout, fmt.c_str(), std::forward<Args>(args)...);
    }

    std::ostringstream oss;
    detail::write_args(oss, std::forward<Args>(args)...);
    return fortran_write_unit(unit, oss.str());
}

template<typename... Args>
int fprintf(int unit, const char* fmt, Args... args) {
    FILE* fp = fortran_get_unit(unit);
    return std::fprintf(fp ? fp : stdout, fmt, args...);
}

#ifdef F2CPP_HAS_MPI
template<typename View>
int fortran_mpi_file_write_at(MPI_File fh, MPI_Offset offset, View&& view,
                              int count, MPI_Datatype dtype, MPI_Status* status) {
    auto packed = fortran_pack(view);
    const int n = count >= 0
        ? std::min(count, static_cast<int>(packed.size()))
        : static_cast<int>(packed.size());
    return MPI_File_write_at(fh, offset, packed.data(), n, dtype,
                             status ? status : MPI_STATUS_IGNORE);
}

template<typename View>
int fortran_mpi_file_write_at_all(MPI_File fh, MPI_Offset offset, View&& view,
                                  int count, MPI_Datatype dtype, MPI_Status* status) {
    auto packed = fortran_pack(view);
    const int n = count >= 0
        ? std::min(count, static_cast<int>(packed.size()))
        : static_cast<int>(packed.size());
    return MPI_File_write_at_all(fh, offset, packed.data(), n, dtype,
                                 status ? status : MPI_STATUS_IGNORE);
}

template<typename View>
int fortran_mpi_file_read_at(MPI_File fh, MPI_Offset offset, View&& view,
                             int count, MPI_Datatype dtype, MPI_Status* status) {
    const size_t target = view.total_size();
    const size_t requested = count >= 0 ? static_cast<size_t>(count) : target;
    using Value = std::remove_cv_t<std::remove_reference_t<decltype(*view.data())>>;
    std::vector<Value> packed(std::max(target, requested));
    const int n = static_cast<int>(requested);
    const int rc = MPI_File_read_at(fh, offset, packed.data(), n, dtype,
                                    status ? status : MPI_STATUS_IGNORE);
    fortran_unpack(view, packed.data(), std::min(target, requested));
    return rc;
}

template<typename View>
int fortran_mpi_file_read_at_all(MPI_File fh, MPI_Offset offset, View&& view,
                                 int count, MPI_Datatype dtype, MPI_Status* status) {
    const size_t target = view.total_size();
    const size_t requested = count >= 0 ? static_cast<size_t>(count) : target;
    using Value = std::remove_cv_t<std::remove_reference_t<decltype(*view.data())>>;
    std::vector<Value> packed(std::max(target, requested));
    const int n = static_cast<int>(requested);
    const int rc = MPI_File_read_at_all(fh, offset, packed.data(), n, dtype,
                                        status ? status : MPI_STATUS_IGNORE);
    fortran_unpack(view, packed.data(), std::min(target, requested));
    return rc;
}
#endif

#endif
