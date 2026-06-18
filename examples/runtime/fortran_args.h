#ifndef FORTRAN_ARGS_H
#define FORTRAN_ARGS_H
/*
 * fortran_args.h - Fortran command_argument_count / get_command_argument runtime.
 *
 * In Fortran, command_argument_count() and get_command_argument(i, arg)
 * are global intrinsics callable from any scope. In C++, argc/argv are
 * only available in main(). This header bridges the gap by storing
 * argc/argv globally at program start.
 *
 * Usage in main():
 *   #include "fortran_args.h"
 *   int main(int argc, char** argv) {
 *       fortran_args_init(argc, argv);  // call ONCE before anything else
 *       ...
 *   }
 *
 * Usage anywhere:
 *   int n = command_argument_count();
 *   std::string arg = get_command_argument(i);  // 0 = program name, 1..n = args
 */

#include <string>

namespace fortran_args_detail {
    inline int    g_argc = 0;
    inline char** g_argv = nullptr;
}

inline void fortran_args_init(int argc, char** argv) {
    fortran_args_detail::g_argc = argc;
    fortran_args_detail::g_argv = argv;
}

inline int command_argument_count() {
    return fortran_args_detail::g_argc - 1;  // Fortran excludes program name
}

inline std::string get_command_argument(int i) {
    if (i < 0 || i >= fortran_args_detail::g_argc) return "";
    return std::string(fortran_args_detail::g_argv[i]);
}

#endif // FORTRAN_ARGS_H
