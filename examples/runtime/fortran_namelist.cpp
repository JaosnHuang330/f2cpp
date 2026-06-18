#include "fortran_namelist.h"

namespace {

bool is_key_start(char c) {
    return std::isalpha(static_cast<unsigned char>(c)) || c == '_';
}

bool is_key_char(char c) {
    return std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '%';
}

size_t skip_balanced_parens(const std::string& s, size_t pos) {
    if (pos >= s.size() || s[pos] != '(') return pos;
    int depth = 0;
    for (size_t i = pos; i < s.size(); ++i) {
        if (s[i] == '(') depth++;
        else if (s[i] == ')') {
            depth--;
            if (depth == 0) return i + 1;
        }
    }
    return s.size();
}

bool looks_like_key_assignment(const std::string& s, size_t start, size_t* key_pos = nullptr) {
    size_t p = start;
    while (p < s.size() && (std::isspace(static_cast<unsigned char>(s[p])) || s[p] == ',')) p++;
    if (p >= s.size() || !is_key_start(s[p])) return false;
    size_t candidate = p;
    while (p < s.size()) {
        if (is_key_char(s[p])) {
            p++;
            continue;
        }
        if (s[p] == '(') {
            p = skip_balanced_parens(s, p);
            continue;
        }
        break;
    }
    while (p < s.size() && std::isspace(static_cast<unsigned char>(s[p]))) p++;
    if (p < s.size() && s[p] == '=') {
        if (key_pos) *key_pos = candidate;
        return true;
    }
    return false;
}

std::string strip_group_terminator(const std::string& line, bool& ended) {
    bool in_str = false;
    char quote_char = 0;
    for (size_t i = 0; i < line.size(); ++i) {
        char c = line[i];
        if (!in_str && (c == '\'' || c == '"')) {
            in_str = true;
            quote_char = c;
        } else if (in_str && c == quote_char) {
            in_str = false;
        } else if (!in_str && c == '/') {
            // Fortran standard: '/' terminates a namelist group; anything after
            // it on the same line is a comment (e.g. `/End of NAMELIST "..."`).
            // Only require that '/' is not embedded inside a token.
            bool before_ok = (i == 0) || std::isspace(static_cast<unsigned char>(line[i - 1])) ||
                             line[i - 1] == ',';
            if (before_ok) {
                ended = true;
                return line.substr(0, i);
            }
        } else if (!in_str && c == '&') {
            std::string rest = line.substr(i);
            std::transform(rest.begin(), rest.end(), rest.begin(),
                           [](unsigned char ch){ return std::tolower(ch); });
            if (rest.rfind("&end", 0) == 0) {
                ended = true;
                return line.substr(0, i);
            }
        }
    }
    return line;
}

size_t find_next_key_assignment(const std::string& s, size_t start) {
    bool in_q = false;
    char qc = 0;
    int paren_depth = 0;
    for (size_t i = start; i < s.size(); ++i) {
        char c = s[i];
        if (!in_q && (c == '\'' || c == '"')) {
            in_q = true;
            qc = c;
        } else if (in_q && c == qc) {
            in_q = false;
        } else if (!in_q && c == '(') {
            paren_depth++;
        } else if (!in_q && c == ')' && paren_depth > 0) {
            paren_depth--;
        } else if (!in_q && paren_depth == 0 &&
                   (c == ',' || std::isspace(static_cast<unsigned char>(c)))) {
            if (looks_like_key_assignment(s, i + 1)) return i;
        }
    }
    return std::string::npos;
}

}

std::string NamelistGroup::to_lower(const std::string& s) {
    std::string r = s;
    std::transform(r.begin(), r.end(), r.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return r;
}

std::string NamelistGroup::trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

// Internal: parse stream into entries_ for a given group name
bool NamelistGroup::_parse_stream_into_entries(std::istream& file) {
    entries_.clear();
    std::string target = to_lower("&" + groupName_);
    std::string line;
    bool in_group = false;
    std::string accumulated;

    while (std::getline(file, line)) {
        std::string stripped;
        bool in_str = false;
        char quote_char = 0;
        for (size_t i = 0; i < line.size(); ++i) {
            char c = line[i];
            if (!in_str && (c == '\'' || c == '"')) {
                in_str = true; quote_char = c; stripped += c;
            } else if (in_str && c == quote_char) {
                in_str = false; stripped += c;
            } else if (!in_str && c == '!') {
                break;
            } else { stripped += c; }
        }
        std::string trimmed = trim(stripped);
        if (!in_group) {
            std::string lower_trimmed = to_lower(trimmed);
            bool target_match = lower_trimmed.rfind(target, 0) == 0 &&
                (lower_trimmed.size() == target.size() ||
                 std::isspace(static_cast<unsigned char>(lower_trimmed[target.size()])) ||
                 lower_trimmed[target.size()] == ',');
            if (target_match) {
                in_group = true;
                trimmed = trim(trimmed.substr(target.size()));
                if (trimmed.empty()) continue;
            } else {
                continue;
            }
        }
        bool ended = false;
        std::string usable = trim(strip_group_terminator(trimmed, ended));
        if (!usable.empty()) accumulated += " " + usable;
        if (ended) break;
    }
    if (!in_group) return false;

    // Parse key=value pairs
    std::string buf = accumulated;
    size_t pos = 0;
    while (pos < buf.size()) {
        size_t key_pos = std::string::npos;
        if (!looks_like_key_assignment(buf, pos, &key_pos)) {
            pos++;
            continue;
        }
        size_t eq_pos = buf.find('=', key_pos);
        if (eq_pos == std::string::npos) break;
        std::string key = trim(buf.substr(key_pos, eq_pos - key_pos));
        key = to_lower(key);
        size_t value_start = eq_pos + 1;
        size_t next = find_next_key_assignment(buf, value_start);
        size_t value_end = (next == std::string::npos) ? buf.size() : next;
        std::string value = trim(buf.substr(value_start, value_end - value_start));
        if (!value.empty() && value.back() == ',')
            value = trim(value.substr(0, value.size() - 1));
        entries_[key] = value;
        pos = value_end;
    }
    return true;
}

bool NamelistGroup::_apply_setters() {
    for (auto& kv : setters_) {
        auto it = entries_.find(kv.first);
        if (it != entries_.end()) {
            auto vals = split_values(it->second);
            kv.second(vals);
        }
    }
    return true;
}

bool NamelistGroup::parse(const std::string& filename, const std::string& group_name) {
    groupName_ = group_name;
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    return _parse_stream_into_entries(file);
}

std::vector<std::string> NamelistGroup::split_values(const std::string& s) {
    std::vector<std::string> raw;
    std::string current;
    bool in_str = false;
    char qc = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        char c = s[i];
        if (!in_str && (c == '\'' || c == '"')) { in_str = true; qc = c; current += c; }
        else if (in_str && c == qc) { in_str = false; current += c; }
        else if (!in_str && (c == ',' || c == ' ' || c == '\t')) {
            std::string t = trim(current);
            if (!t.empty()) raw.push_back(t);
            current.clear();
        } else { current += c; }
    }
    std::string t = trim(current);
    if (!t.empty()) raw.push_back(t);

    // Expand Fortran repeat notation: "3*1.0" → "1.0", "1.0", "1.0"
    std::vector<std::string> result;
    for (const auto& token : raw) {
        size_t star = token.find('*');
        if (star != std::string::npos && star > 0 && star < token.size() - 1) {
            // Check that everything before * is a positive integer
            std::string count_str = token.substr(0, star);
            bool all_digit = true;
            for (char ch : count_str) {
                if (!std::isdigit(static_cast<unsigned char>(ch))) { all_digit = false; break; }
            }
            if (all_digit) {
                int count = std::atoi(count_str.c_str());
                std::string val = token.substr(star + 1);
                for (int r = 0; r < count; ++r) result.push_back(val);
                continue;
            }
        }
        result.push_back(token);
    }
    return result;
}

bool NamelistGroup::has_key(const std::string& key) const {
    return entries_.find(to_lower(key)) != entries_.end();
}

bool NamelistGroup::get_int(const std::string& key, int& val) const {
    auto it = entries_.find(to_lower(key));
    if (it == entries_.end()) return false;
    val = std::atoi(it->second.c_str());
    return true;
}

bool NamelistGroup::get_double(const std::string& key, double& val) const {
    auto it = entries_.find(to_lower(key));
    if (it == entries_.end()) return false;
    std::string s = it->second;
    for (auto& c : s) { if (c == 'd' || c == 'D') c = 'e'; }
    val = std::atof(s.c_str());
    return true;
}

bool NamelistGroup::parse_bool_value(const std::string& s, bool& val) {
    std::string lower = to_lower(trim(s));
    if (lower == "t" || lower == ".true." || lower == "true") { val = true; return true; }
    if (lower == "f" || lower == ".false." || lower == "false") { val = false; return true; }
    return false;
}

bool NamelistGroup::get_bool(const std::string& key, bool& val) const {
    auto it = entries_.find(to_lower(key));
    if (it == entries_.end()) return false;
    return parse_bool_value(it->second, val);
}

bool NamelistGroup::get_string(const std::string& key, std::string& val) const {
    auto it = entries_.find(to_lower(key));
    if (it == entries_.end()) return false;
    std::string s = it->second;
    if (s.size() >= 2 && ((s.front() == '\'' && s.back() == '\'') ||
                          (s.front() == '"' && s.back() == '"'))) {
        s = s.substr(1, s.size() - 2);
    }
    val = s;
    return true;
}

bool NamelistGroup::get_int_array(const std::string& key, int* arr, int n) const {
    auto it = entries_.find(to_lower(key));
    if (it == entries_.end()) return false;
    auto vals = split_values(it->second);
    for (int i = 0; i < n && i < (int)vals.size(); ++i) {
        arr[i] = std::atoi(vals[i].c_str());
    }
    return true;
}

bool NamelistGroup::get_double_array(const std::string& key, double* arr, int n) const {
    auto it = entries_.find(to_lower(key));
    if (it == entries_.end()) return false;
    auto vals = split_values(it->second);
    for (int i = 0; i < n && i < (int)vals.size(); ++i) {
        std::string s = vals[i];
        for (auto& c : s) { if (c == 'd' || c == 'D') c = 'e'; }
        arr[i] = std::atof(s.c_str());
    }
    return true;
}
