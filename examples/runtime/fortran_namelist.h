#ifndef FORTRAN_NAMELIST_H
#define FORTRAN_NAMELIST_H

#include <string>
#include <map>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <functional>
#include <cstdio>
#include <initializer_list>
#include <iomanip>
#include <utility>
#include "fortran_io.h"

/**
 * Fortran NAMELIST parser for C++.
 *
 * Supports two usage patterns:
 *
 * Pattern A (parse/get): explicit fetch after parse
 *   NamelistGroup ng;
 *   ng.parse("file.prm", "GroupName");
 *   ng.get_double("key", val);
 *
 * Pattern B (add/read): registration-based, auto-populates variables
 *   NamelistGroup ng("GroupName");
 *   ng.add("key", val);
 *   ng.read(file);   // or ng.read("file.prm")
 */
class NamelistGroup {
public:
    NamelistGroup() {}
    explicit NamelistGroup(const std::string& group_name) : groupName_(group_name) {}

    // ---- Pattern B: Registration API ----

    void add(const std::string& key, double& v) {
        register_writer(key, [&v]() { return num_to_string(v); });
        setters_[to_lower_s(key)] = [&v](const std::vector<std::string>& vals) {
            if (!vals.empty()) try { v = std::stod(fix_exp(vals[0])); } catch(...) {}
        };
    }
    void add(const std::string& key, float& v) {
        register_writer(key, [&v]() { return num_to_string(v); });
        setters_[to_lower_s(key)] = [&v](const std::vector<std::string>& vals) {
            if (!vals.empty()) try { v = (float)std::stod(fix_exp(vals[0])); } catch(...) {}
        };
    }
    void add(const std::string& key, int& v) {
        register_writer(key, [&v]() { return std::to_string(v); });
        setters_[to_lower_s(key)] = [&v](const std::vector<std::string>& vals) {
            if (!vals.empty()) try { v = (int)std::stoll(vals[0]); } catch(...) {}
        };
    }
    void add(const std::string& key, long& v) {
        register_writer(key, [&v]() { return std::to_string(v); });
        setters_[to_lower_s(key)] = [&v](const std::vector<std::string>& vals) {
            if (!vals.empty()) try { v = std::stol(vals[0]); } catch(...) {}
        };
    }
    void add(const std::string& key, bool& v) {
        register_writer(key, [&v]() { return v ? ".true." : ".false."; });
        setters_[to_lower_s(key)] = [&v](const std::vector<std::string>& vals) {
            if (!vals.empty()) {
                std::string s = to_lower_s(vals[0]);
                v = (s == "t" || s == ".true." || s == "true" || s == "1");
            }
        };
    }
    void add(const std::string& key, std::string& v) {
        register_writer(key, [&v]() { return "'" + v + "'"; });
        setters_[to_lower_s(key)] = [&v](const std::vector<std::string>& vals) {
            if (!vals.empty()) {
                std::string s = vals[0];
                // Strip surrounding quotes
                if (s.size() >= 2 && (s.front() == '"' || s.front() == '\'') && s.back() == s.front())
                    s = s.substr(1, s.size() - 2);
                v = s;
            }
        };
    }

    // Compatibility helpers used by generated translators.  They forward to
    // the registration API above and keep the runtime contract stable.
    void add_int(const std::string& key, int* v) { if (v) add(key, *v); }
    void add_int(const std::string& key, int& v) { add(key, v); }
    void add_double(const std::string& key, double* v) { if (v) add(key, *v); }
    void add_double(const std::string& key, double& v) { add(key, v); }
    void add_bool(const std::string& key, bool* v) { if (v) add(key, *v); }
    void add_bool(const std::string& key, bool& v) { add(key, v); }
    void add_string(const std::string& key, std::string* v, int = 0) { if (v) add(key, *v); }
    void add_string(const std::string& key, std::string& v, int = 0) { add(key, v); }
    // Raw pointer arrays
    void add(const std::string& key, double* arr, int n) {
        register_writer(key, [arr, n]() { return join_array(arr, n); });
        setters_[to_lower_s(key)] = [arr, n](const std::vector<std::string>& vals) {
            for (int i = 0; i < n && i < (int)vals.size(); ++i)
                try { arr[i] = std::stod(fix_exp(vals[i])); } catch(...) {}
        };
    }
    void add(const std::string& key, int* arr, int n) {
        register_writer(key, [arr, n]() { return join_array(arr, n); });
        setters_[to_lower_s(key)] = [arr, n](const std::vector<std::string>& vals) {
            for (int i = 0; i < n && i < (int)vals.size(); ++i)
                try { arr[i] = (int)std::stoll(vals[i]); } catch(...) {}
        };
    }
    void add_double_array(const std::string& key, double* arr, int n) { add(key, arr, n); }
    void add_int_array(const std::string& key, int* arr, int n) { add(key, arr, n); }

    // std::vector overloads
    void add(const std::string& key, std::vector<double>& v) {
        register_writer(key, [&v]() {
            std::ostringstream os;
            for (size_t i = 0; i < v.size(); ++i) { if (i) os << ", "; os << std::setprecision(16) << v[i]; }
            return os.str();
        });
        setters_[to_lower_s(key)] = [&v](const std::vector<std::string>& vals) {
            v.resize(vals.size());
            for (size_t i = 0; i < vals.size(); ++i)
                try { v[i] = std::stod(fix_exp(vals[i])); } catch(...) {}
        };
    }
    void add(const std::string& key, std::vector<int>& v) {
        register_writer(key, [&v]() {
            std::ostringstream os;
            for (size_t i = 0; i < v.size(); ++i) { if (i) os << ", "; os << v[i]; }
            return os.str();
        });
        setters_[to_lower_s(key)] = [&v](const std::vector<std::string>& vals) {
            v.resize(vals.size());
            for (size_t i = 0; i < vals.size(); ++i)
                try { v[i] = (int)std::stoll(vals[i]); } catch(...) {}
        };
    }
    void add(const std::string& key, std::vector<float>& v) {
        register_writer(key, [&v]() {
            std::ostringstream os;
            for (size_t i = 0; i < v.size(); ++i) { if (i) os << ", "; os << v[i]; }
            return os.str();
        });
        setters_[to_lower_s(key)] = [&v](const std::vector<std::string>& vals) {
            v.resize(vals.size());
            for (size_t i = 0; i < vals.size(); ++i)
                try { v[i] = (float)std::stod(fix_exp(vals[i])); } catch(...) {}
        };
    }
    void add(const std::string& key, std::vector<std::string>& v) {
        register_writer(key, [&v]() {
            std::ostringstream os;
            for (size_t i = 0; i < v.size(); ++i) { if (i) os << ", "; os << "'" << v[i] << "'"; }
            return os.str();
        });
        setters_[to_lower_s(key)] = [&v](const std::vector<std::string>& vals) {
            v.resize(vals.size());
            for (size_t i = 0; i < vals.size(); ++i) {
                std::string s = vals[i];
                if (s.size() >= 2 && (s.front() == '"' || s.front() == '\'') && s.back() == s.front())
                    s = s.substr(1, s.size() - 2);
                v[i] = s;
            }
        };
    }
    template<size_t N>
    void add(const std::string& key, std::array<double, N>& v) {
        register_writer(key, [&v]() {
            std::ostringstream os;
            for (size_t i = 0; i < N; ++i) { if (i) os << ", "; os << num_to_string(v[i]); }
            return os.str();
        });
        setters_[to_lower_s(key)] = [&v](const std::vector<std::string>& vals) {
            for (size_t i = 0; i < N && i < vals.size(); ++i)
                try { v[i] = std::stod(fix_exp(vals[i])); } catch(...) {}
        };
    }
    template<size_t N>
    void add(const std::string& key, std::array<float, N>& v) {
        register_writer(key, [&v]() {
            std::ostringstream os;
            for (size_t i = 0; i < N; ++i) { if (i) os << ", "; os << v[i]; }
            return os.str();
        });
        setters_[to_lower_s(key)] = [&v](const std::vector<std::string>& vals) {
            for (size_t i = 0; i < N && i < vals.size(); ++i)
                try { v[i] = static_cast<float>(std::stod(fix_exp(vals[i]))); } catch(...) {}
        };
    }
    template<size_t N>
    void add(const std::string& key, std::array<int, N>& v) {
        register_writer(key, [&v]() {
            std::ostringstream os;
            for (size_t i = 0; i < N; ++i) { if (i) os << ", "; os << v[i]; }
            return os.str();
        });
        setters_[to_lower_s(key)] = [&v](const std::vector<std::string>& vals) {
            for (size_t i = 0; i < N && i < vals.size(); ++i)
                try { v[i] = static_cast<int>(std::stoll(vals[i])); } catch(...) {}
        };
    }
    template<size_t N>
    void add(const std::string& key, std::array<bool, N>& v) {
        register_writer(key, [&v]() {
            std::ostringstream os;
            for (size_t i = 0; i < N; ++i) { if (i) os << ", "; os << (v[i] ? ".true." : ".false."); }
            return os.str();
        });
        setters_[to_lower_s(key)] = [&v](const std::vector<std::string>& vals) {
            for (size_t i = 0; i < N && i < vals.size(); ++i) {
                std::string s = to_lower_s(vals[i]);
                v[i] = (s == "t" || s == ".true." || s == "true" || s == "1");
            }
        };
    }
    template<size_t N>
    void add(const std::string& key, std::array<std::string, N>& v) {
        register_writer(key, [&v]() {
            std::ostringstream os;
            for (size_t i = 0; i < N; ++i) { if (i) os << ", "; os << "'" << v[i] << "'"; }
            return os.str();
        });
        setters_[to_lower_s(key)] = [&v](const std::vector<std::string>& vals) {
            for (size_t i = 0; i < N && i < vals.size(); ++i) {
                std::string s = vals[i];
                if (s.size() >= 2 && (s.front() == '"' || s.front() == '\'') && s.back() == s.front())
                    s = s.substr(1, s.size() - 2);
                v[i] = s;
            }
        };
    }
    template<typename T>
    void add_double_array(const std::string& key, T& arr, int n) {
        auto* ptr = arr.data();
        int count = std::min(n, arr.size(1));
        add(key, ptr, count);
    }
    template<typename T>
    void add_int_array(const std::string& key, T& arr, int n) {
        auto* ptr = arr.data();
        int count = std::min(n, arr.size(1));
        add(key, ptr, count);
    }
    template<typename T>
    void add_variable(const std::string& key, T& v) { add(key, v); }
    // Generic template fallback for FortranArray and other types - tries .data()/.size()
    template<typename T>
    void add(const std::string& key, T& v) {
        // For array types with .data() and .size(int) (like FortranArray1D<double>)
        _add_array_if_possible(key, v, 0);
    }

    // Canonical entry points for reading namelist input.
    //
    // CANONICAL (use these in translated code):
    //   read(int unit)      -- Fortran-style: open(unit,file=...); read(unit, nml=...)
    //   read(int* unit)     -- same, pointer variant
    //   read(string path)   -- direct path; opens its own ifstream internally
    //
    // INTERNAL (kept for the legacy load_namelist_file() helpers below; new
    // translations should NOT call this directly — use read(int unit) instead):
    //   read(FILE*)
    //
    // REMOVED (do NOT add back without a deprecation cycle):
    //   read(std::ifstream&)  -- caller-managed stream was a footgun: order-
    //                            dependent semantics differed from FILE* path,
    //                            triggered rst_21 `nycSection Wrong 1` abort.
    //
    bool read(FILE* file) {
        if (!file) return false;
        std::rewind(file);
        std::stringstream ss;
        char buffer[4096];
        while (std::fgets(buffer, sizeof(buffer), file)) {
            ss << buffer;
        }
        return _parse_stream_into_entries(ss) && _apply_setters();
    }
    bool read(int unit) {
        return read(fortran_get_unit(unit));
    }
    bool read(int* unit) {
        return unit ? read(*unit) : false;
    }
    bool read(const std::string& filename) {
        std::ifstream f(filename);
        if (!f.is_open()) return false;
        return _parse_stream_into_entries(f) && _apply_setters();
    }

    void write_namelist(int unit) const {
        std::stringstream ss(to_string());
        std::string line;
        while (std::getline(ss, line)) {
            fortran_write_unit(unit, line);
        }
    }
    void write_to_unit(int unit) const { write_namelist(unit); }
    std::string to_string() const {
        std::ostringstream os;
        os << "&" << groupName_ << "\n";
        for (const auto& key : writer_order_) {
            auto it = writers_.find(to_lower_s(key));
            if (it != writers_.end()) {
                os << "  " << key << " = " << it->second() << "\n";
            }
        }
        os << "/";
        return os.str();
    }

    // ---- Pattern A: parse/get API ----
    bool parse(const std::string& filename, const std::string& group_name);

    bool get_int(const std::string& key, int& val) const;
    bool get_double(const std::string& key, double& val) const;
    bool get_bool(const std::string& key, bool& val) const;
    bool get_string(const std::string& key, std::string& val) const;
    bool get_int_array(const std::string& key, int* arr, int n) const;
    bool get_double_array(const std::string& key, double* arr, int n) const;
    bool has_key(const std::string& key) const;

private:
    std::string groupName_;
    std::map<std::string, std::string> entries_;
    std::map<std::string, std::function<void(const std::vector<std::string>&)>> setters_;
    std::map<std::string, std::function<std::string()>> writers_;
    std::vector<std::string> writer_order_;

    static std::string to_lower(const std::string& s);
    static std::string to_lower_s(const std::string& s) {
        std::string r = s;
        std::transform(r.begin(), r.end(), r.begin(),
                       [](unsigned char c){ return std::tolower(c); });
        return r;
    }
    static std::string fix_exp(const std::string& s) {
        // Fortran D/d exponent notation -> E
        std::string r = s;
        for (char& c : r) if (c == 'd' || c == 'D') { c = 'e'; break; }
        return r;
    }
    static std::string trim_s(const std::string& s) {
        size_t a = s.find_first_not_of(" \t\r\n");
        if (a == std::string::npos) return "";
        size_t b = s.find_last_not_of(" \t\r\n");
        return s.substr(a, b - a + 1);
    }
    static std::string num_to_string(double v) {
        std::ostringstream os;
        os << std::setprecision(16) << v;
        return os.str();
    }
    template<typename T>
    static std::string join_array(T* arr, int n) {
        std::ostringstream os;
        for (int i = 0; i < n; ++i) {
            if (i) os << ", ";
            os << arr[i];
        }
        return os.str();
    }
    void register_writer(const std::string& key, std::function<std::string()> fn) {
        std::string lower = to_lower_s(key);
        if (writers_.find(lower) == writers_.end()) {
            writer_order_.push_back(key);
        }
        writers_[lower] = std::move(fn);
    }
    static std::vector<std::string> split_vals(const std::string& s);

    bool _parse_stream_into_entries(std::istream& stream);
    bool _apply_setters();

    static std::string trim(const std::string& s);
    static std::vector<std::string> split_values(const std::string& s);
    static bool parse_bool_value(const std::string& s, bool& val);

    // SFINAE helper: if T has .data() returning double* and .size() returning int
    template<typename T>
    auto _add_array_if_possible(const std::string& key, T& v, int) ->
        decltype(v.data(), v.size(1), void()) {
        auto* ptr = v.data();
        int n = v.size(1);
        register_writer(key, [ptr, n]() { return join_array(ptr, n); });
        setters_[to_lower_s(key)] = [ptr, n](const std::vector<std::string>& vals) {
            for (int i = 0; i < n && i < (int)vals.size(); ++i)
                try { ptr[i] = std::stod(fix_exp(vals[i])); } catch(...) {}
        };
    }
    template<typename T>
    void _add_array_if_possible(const std::string& key, T& v, long) {
        // Fallback: do nothing for unrecognized types
        (void)key; (void)v;
    }
};

// Compatibility facade for generated code that treats a namelist file as an
// object with a current group.  It delegates to NamelistGroup so both API
// styles remain supported.
class FortranNamelist {
public:
    bool open(const std::string& filename) {
        filename_ = filename;
        return true;
    }

    void close() {}

    bool read_group(const std::string& group_name) {
        current_ = NamelistGroup();
        return current_.parse(filename_, group_name);
    }

    int get_int(const std::string& key, int default_value = 0) const {
        int value = default_value;
        current_.get_int(key, value);
        return value;
    }

    double get_real(const std::string& key, double default_value = 0.0) const {
        double value = default_value;
        current_.get_double(key, value);
        return value;
    }

    double get_double(const std::string& key, double default_value = 0.0) const {
        return get_real(key, default_value);
    }

    bool get_logical(const std::string& key, bool default_value = false) const {
        bool value = default_value;
        current_.get_bool(key, value);
        return value;
    }

    bool get_bool(const std::string& key, bool default_value = false) const {
        return get_logical(key, default_value);
    }

    std::string get_string(const std::string& key,
                           const std::string& default_value = "") const {
        std::string value = default_value;
        current_.get_string(key, value);
        return value;
    }

    void get_real_array(const std::string& key, double* arr, int n) const {
        current_.get_double_array(key, arr, n);
    }

    void get_double_array(const std::string& key, double* arr, int n) const {
        get_real_array(key, arr, n);
    }

    void get_int_array(const std::string& key, int* arr, int n) const {
        current_.get_int_array(key, arr, n);
    }

private:
    std::string filename_;
    NamelistGroup current_;
};

inline bool load_namelist_file(FILE* file, NamelistGroup& group) {
    return group.read(file);
}

inline bool load_namelist_file(int unit, NamelistGroup& group) {
    return group.read(fortran_get_unit(unit));
}

inline bool load_namelist_file(int* unit, NamelistGroup& group) {
    return unit ? load_namelist_file(*unit, group) : false;
}

inline bool load_namelist_file(FILE* file, std::initializer_list<NamelistGroup*> groups) {
    bool ok = true;
    for (auto* group : groups) {
        if (!group) continue;
        if (file) std::rewind(file);
        ok = group->read(file) && ok;
    }
    return ok;
}

inline bool load_namelist_file(int unit, std::initializer_list<NamelistGroup*> groups) {
    return load_namelist_file(fortran_get_unit(unit), groups);
}

inline bool load_namelist_file(int* unit, std::initializer_list<NamelistGroup*> groups) {
    return unit ? load_namelist_file(*unit, groups) : false;
}

#endif // FORTRAN_NAMELIST_H
