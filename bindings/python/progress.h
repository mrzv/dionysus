#pragma once

#include <dionysus/dlog/progress.h>

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

    mutable dlog::progress  progress;
};

struct NoProgress: public Progress
{
    void    operator()() const override         {}
};
