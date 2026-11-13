#include "core/geogramContext.h"

#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/common.h>
#include <geogram/basic/logger.h>

namespace PP
{

namespace
{

std::mutex&
_GetGeogramMutex()
{
    static std::mutex m;
    return m;
}

void
_InitGeogram()
{
    static bool initialized = false;
    if (initialized) {
        return;
    }
    
    // No signal / FPE handlers: we may be running inside a host application.
    GEO::initialize(GEO::GEOGRAM_INSTALL_NONE);
    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("algo");
    GEO::Logger::instance()->set_quiet(true);
    initialized = true;
}

} // namespace

GeogramLock::GeogramLock()
    : _lock(_GetGeogramMutex())
{
    _InitGeogram();
}

} // namespace PP
