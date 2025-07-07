// The MIT License (MIT)
//
// Copyright (c) 2019 Luigi Pertoldi
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
//
// ============================================================================
//  ___   ___   ___   __    ___   ____  __   __   ___    __    ___
// | |_) | |_) / / \ / /`_ | |_) | |_  ( (` ( (` | |_)  / /\  | |_)
// |_|   |_| \ \_\_/ \_\_/ |_| \ |_|__ _)_) _)_) |_|_) /_/--\ |_| \_
//
// Very simple progress bar for c++ loops with internal running variable
//
// Author: Luigi Pertoldi
// Created: 3 dic 2016
//
// Notes: The bar must be used when there's no other possible source of output
//        inside the for loop
//

#ifndef __PROGRESSBAR_HPP
#define __PROGRESSBAR_HPP

#include <string>
#include <stdexcept>
#include <Rcpp.h>    // for Rcerr
using Rcpp::Rcerr;

class progressbar {
public:
    ~progressbar() = default;
    progressbar(const progressbar&) = delete;
    progressbar& operator=(const progressbar&) = delete;
    progressbar(progressbar&&) = delete;
    progressbar& operator=(progressbar&&) = delete;

    // constructor: write to Rcerr
    progressbar();
    progressbar(int n, bool showbar=true);

    void reset();
    void set_niter(int iter);
    inline void set_done_char(const std::string& sym) { done_char = sym; }
    inline void set_todo_char(const std::string& sym) { todo_char = sym; }
    inline void show_bar(bool flag=true) { do_show_bar = flag; }
    inline void set_opening_bracket_char(const std::string& sym) { opening_bracket_char = sym; }
    inline void set_closing_bracket_char(const std::string& sym) { closing_bracket_char = sym; }
    void update();

private:
    int progress;
    int n_cycles;
    int last_perc;
    bool do_show_bar;
    bool update_is_called;
    std::string done_char, todo_char;
    std::string opening_bracket_char, closing_bracket_char;
};

inline progressbar::progressbar()
    : progress(0), n_cycles(0), last_perc(0),
      do_show_bar(true), update_is_called(false),
      done_char("#"), todo_char(" "),
      opening_bracket_char("["), closing_bracket_char("]") {}

inline progressbar::progressbar(int n, bool showbar)
    : progress(0), n_cycles(n), last_perc(0),
      do_show_bar(showbar), update_is_called(false),
      done_char("#"), todo_char(" "),
      opening_bracket_char("["), closing_bracket_char("]") {}

inline void progressbar::reset() {
    progress = 0;
    update_is_called = false;
    last_perc = 0;
}

inline void progressbar::set_niter(int iter) {
    if (iter <= 0) throw std::invalid_argument("progressbar::set_niter: iter<=0");
    n_cycles = iter;
}

inline void progressbar::update() {
    if (n_cycles == 0) throw std::runtime_error("progressbar::update: no cycles set");

    if (!update_is_called) {
        if (do_show_bar) {
            Rcerr << opening_bracket_char;
            for (int i = 0; i < 50; ++i) Rcerr << todo_char;
            Rcerr << closing_bracket_char << " 0%";
        } else {
            Rcerr << "0%";
        }
        Rcerr << "\n";
    }
    update_is_called = true;

    int perc = progress * 100 / (n_cycles-1);
    if (perc < last_perc) return;

    if (perc > last_perc) {
        Rcerr << perc << "%\n";
        last_perc = perc;
    }
    ++progress;
}

#endif // __PROGRESSBAR_HPP
