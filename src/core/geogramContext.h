#ifndef PP_CORE_GEOGRAMCONTEXT_H
#define PP_CORE_GEOGRAMCONTEXT_H

#include <mutex>

namespace PP
{

// Scoped access to geogram. geogram keeps global state (command line,
// factories, logger), so all calls into it are serialized. The first lock
// also initializes geogram.
class GeogramLock
{
public:
    GeogramLock();
    ~GeogramLock() = default;

    GeogramLock(const GeogramLock&) = delete;
    GeogramLock& operator=(const GeogramLock&) = delete;

private:
    std::lock_guard<std::mutex> _lock;
};

} // namespace PP

#endif // PP_CORE_GEOGRAMCONTEXT_H
