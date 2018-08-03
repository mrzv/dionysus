#pragma once

#include <iostream>
#include <pybind11/iostream.h>
#include <dionysus/dlog/progress.h>
namespace py = pybind11;

struct Progress
{
    virtual void operator()() const = 0;
    virtual ~Progress() {}
};

struct ShowProgress: public Progress
{
        ShowProgress(size_t total):
            progress(total)                     {}

    void    operator()() const override         { ++progress; }

    py::scoped_ostream_redirect stream = py::scoped_ostream_redirect(
        std::cout,                               // std::ostream&
        py::module::import("sys").attr("stdout") // Python output
    );
    mutable dlog::progress  progress;
};

struct NoProgress: public Progress
{
    void    operator()() const override         {}
};
